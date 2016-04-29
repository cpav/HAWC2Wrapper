import copy
import numpy as np
import os
import shutil

from wetb.prepost import dlcdefs
from openmdao.api import Component, Group, ParallelGroup

from hawc2_inputreader import HAWC2InputReader
from hawc2_inputwriter import HAWC2InputWriter
from hawc2_wrapper import HAWC2Wrapper
from hawc2_output import HAWC2OutputCompact
from hawc2_geometry import HAWC2GeometryBuilder


class HAWC2Workflow(Component):
    """OpenMDAO component to run the HAWC2S workflow.

    Parameters
    ----------
    config: dict
        Configuration dictionary. It has to contain the following entries:

        * 'with_structure': bool to add structural properties to the parameters
        * 'with_geom': bool to add the blade geometry to the parameters
        * 'master_file': str with the name of the master file.
        * 'aerodynamic_sections': int of the number of aerodynamic sections.
        * 'HAWC2Outputs': dict of the outputs required.
        * 'HAWC2InputWriter': dict for initialization of HAWC2SInputWriter
          parameters.
        * 'HAWC2Wrapper': dict for initialization of HAWC2Wrapper parameters.
        * 'HAWC2GeometryBuilder': dict for initialization of
          HAWC2GeometryBuilder parameters

    case_id: str
        Name of the HAWC2 case to create and run.

    case: list/dict



    Returns
    -------

    """
    def __init__(self, config, case, cssize, pfsize):
        super(HAWC2Workflow, self).__init__()

        case_id = case['[case_id]']
        self.case = case
        self.basedir = os.getcwd()
        self.keep_work_dirs = False

        self.with_structure = config['with_structure']
        self.with_geom = config['with_geom']
        self.with_tsr = config['with_tsr']

        self.reader = HAWC2InputReader(config['master_file'])
        self.reader.execute()

        self.writer = HAWC2InputWriter(**config['HAWC2InputWriter'])

        self.writer.turb_directory = config['HAWC2InputWriter']['turb_directory']
        self.writer.data_directory = config['HAWC2InputWriter']['data_directory']

        self.writer.vartrees = copy.copy(self.reader.vartrees)
        self.writer.case_id = case_id
        self.writer.vartrees.aero.ae_filename = \
            os.path.join(self.writer.data_directory, case_id+'_ae.dat')
        self.writer.vartrees.aero.pc_filename = \
            os.path.join(self.writer.data_directory, case_id+'_pc.dat')

        self.writer.vartrees.aero.aerosections = config['aerodynamic_sections']
        # Update vartrees with variables included in "case"
        self.writer.vartrees.tags2var(case)

        self.wrapper = HAWC2Wrapper(**config['HAWC2Wrapper'])
        self.wrapper.case_id = case_id
        self.wrapper.log_directory = case['[log_dir]']

        self.output = HAWC2OutputCompact(config['HAWC2Outputs'])

        self.add_output('outputs_statistics', shape=[self.output.Nstat,
                                                     self.output.Nch])
        self.add_output('outputs_fatigue', shape=[len(self.output.m),
                                                  self.output.Nch])

        if self.with_tsr:
            self.add_param('tsr', 0.)

        if self.with_structure:
            self.add_param('blade_beam_structure', shape=cssize)

        self.geom = HAWC2GeometryBuilder(config['structural_sections'], **config['BladeGeometryBuilder'])
        self.geom.c12axis_init = self.reader.vartrees.main_bodies.blade1.c12axis.copy()
        self.geom.c12axis_init[:, :3] /= self.geom.c12axis_init[-1, 2]

        if self.with_geom:
            self.geom_var = ['s', 'x', 'y', 'z', 'rot_x', 'rot_y', 'rot_z',
                             'chord', 'rthick', 'p_le']
            for v in self.geom_var:
                self.add_param(v, shape=[pfsize])
            self.add_param('blade_length', 0.)
            self.geom.interp_from_htc = False
        else:
            self.geom.interp_from_htc = True

    def solve_nonlinear(self, params, unknowns, resids):

        workdir = 'hawc2_model_%s_%i' % (self.writer.case_id, self.__hash__())

        try:
            os.mkdir(workdir)
        except:
            pass

        shutil.copytree('control', workdir+'/control')
        os.chdir(workdir)

        vt = self.writer.vartrees
        if self.with_tsr:
            vt.dlls.risoe_controller.dll_init.designTSR = \
                params['tsr']

        if self.with_geom:
            blade_length = params['blade_length']
            for v in self.geom_var:
                setattr(self.geom.bladegeom, v, params[v])
        else:
            blade_length = vt.main_bodies.blade1.c12axis[-1, 2]
            # conversion from HAWC2BladeGeometry to BladeGeometryVT
            self.geom.bladegeom.s = vt.blade_ae.s
            self.geom.bladegeom.rthick = vt.blade_ae.rthick/100.
            self.geom.bladegeom.chord = vt.blade_ae.chord/blade_length

        self.geom.blade_length = blade_length
        self.geom.execute()

        vt.blade_ae = self.geom.blade_ae
        vt.main_bodies.blade1.c12axis = self.geom.c12axis

        if self.with_structure:
            body = vt.main_bodies.blade1
            body.array2hawc2beamstructure(blade_length, params['blade_beam_structure'])

        self.writer.execute()
        self.wrapper.compute()

        self.output.execute(self.case)

        unknowns['outputs_statistics'] = self.output.outputs_statistics
        unknowns['outputs_fatigue'] = self.output.outputs_fatigue

        os.chdir(self.basedir)


class OutputsAggregator(Component):
    """
    OpenMDAO component to gather the outputs from the ParallelGroup.

    parameters
    ----------

    config: dict
        Configuration dictionary. Requires:
        * 'aerodynamic_sections'. Number of aerodynamic sections.
        * 'HAWC2SOutputs'. Dictionary of the outputs.

    n_cases: int
        Number of cases evaluated.

    returns
    -------
    Output Arrays

    """
    def __init__(self, config, n_cases):
        super(OutputsAggregator, self).__init__()

        self.n_cases = n_cases
        self.naes = config['aerodynamic_sections']-2
        self.channel = config['HAWC2Outputs']['channels']

        # Initialize to get all the required constants
        Nch = len(config['HAWC2Outputs']['channels'])
        m = config['HAWC2Outputs']['m']
        stat_list = config['HAWC2Outputs']['stat_list']
        Nstat = len(stat_list)
        # Add parameters coming from ParallelGroup
        for i in range(n_cases):
            p_name = 'outputs_statistics_%d' % i
            self.add_param(p_name, shape=[Nstat, Nch])
        for i in range(n_cases):
            p_name = 'outputs_fatigue_%d' % i
            self.add_param(p_name, shape=[len(m), Nch])

        # Add outputs
        self.stat_var = []
        for stat in stat_list:
            name = 'stat_'+stat
            self.stat_var.append(name)
            self.add_output(name, shape=[n_cases, Nch])

        self.fatigue_var = []
        for im in m:
            name = 'fatigue_m%i' % im
            self.fatigue_var.append(name)
            self.add_output(name, shape=[n_cases, Nch])

    def solve_nonlinear(self, params, unknowns, resids):

        for j, var in enumerate(self.stat_var):
            for i in range(self.n_cases):
                p_name = 'outputs_statistics_%d' % i
                out = params[p_name]
                unknowns[var][i, :] = out[j, :]

        for j, var in enumerate(self.fatigue_var):
            for i in range(self.n_cases):
                p_name = 'outputs_fatigue_%d' % i
                out = params[p_name]
                unknowns[var][i, :] = out[j, :]


class HAWC2AeroElasticSolver(Group):
    """
    OpenMDAO group to execute HAWC2 in parallel.

    parameters
    -----------
    config: dict
        Configuration dictionary.

    """
    def __init__(self, config, dlcs_folder, cssize, pfsize):
        super(HAWC2AeroElasticSolver, self).__init__()

        # check that the config is ok
        self._check_config(config)

        # load cases and their tags
        dlcs = dlcdefs.excel_stabcon(dlcs_folder, fext='xls', silent=True)

        promote = []
        if config['with_tsr']:
            promote.append('tsr')

        if config['with_structure']:
            promote.append('blade_beam_structure')

        if config['with_geom']:
            var = ['s', 'x', 'y', 'z', 'rot_x', 'rot_y', 'rot_z',
                             'chord', 'rthick', 'p_le']
            for v in var:
                promote.append(v)
            promote.append('blade_length')

        # Add outputs
        agg_promo = []
        for stat in config['HAWC2Outputs']['stat_list']:
            name = 'stat_'+stat
            agg_promo.append(name)
        for m in config['HAWC2Outputs']['m']:
            name = 'fatigue_m%i' % m
            agg_promo.append(name)

        self.add('aggregate', OutputsAggregator(config, len(dlcs)),
                 promotes=agg_promo)
        pg = self.add('pg', ParallelGroup(), promotes=promote)

        for icase, case in enumerate(dlcs):
            case_id = case['[case_id]'].replace('-', '_')
            pg.add(case_id, HAWC2Workflow(config, case, cssize, pfsize),
                   promotes=promote)

            self.connect('pg.%s.outputs_statistics' % case_id,
                         'aggregate.outputs_statistics_%d' % icase)
            self.connect('pg.%s.outputs_fatigue' % case_id,
                         'aggregate.outputs_fatigue_%d' % icase)

    def _check_config(self, config):

        if 'master_file' not in config.keys():
            raise RuntimeError('You need to supply the name of the master' +
                               'file in the configuration dictionary.')

        if 'BladeGeometryBuilder' not in config.keys():
            config['BladeGeometryBuilder'] = {}

        if 'HAWC2InputWriter' not in config.keys():
            raise RuntimeError('You need to supply a config dict' +
                               'for HAWC2InputWriter.')

        if 'HAWC2Wrapper' not in config.keys():
            raise RuntimeError('You need to supply a config dict' +
                               'for HAWC2Wrapper.')

        if 'HAWC2Outputs' not in config.keys():
            raise RuntimeError('You need to supply a config dict' +
                               'for HAWC2Outputs.')

        if 'aerodynamic_sections' not in config.keys():
            config['aerodynamic_sections'] = 40
            print 'The number of aerodynamic section has not been specified' +\
                  ' in the configuration dictionary. Default value of 40 selected.'

        if 'with_tsr' not in config.keys():
            config['with_tsr'] = False
            print 'Tip-speed-ratio is not set as parameters because' +\
                  ' no option "with_tsr" was given in the configuration.'

        if 'with_structure' not in config.keys():
            config['with_structure'] = False
            print 'Structural properties are not set as parameters because' +\
                  ' no option "with_structure" was given in the configuration.'

        if 'with_geom' not in config.keys():
            config['with_geom'] = False
            print 'Blade planform properties are not set as parameters because' +\
                  ' no option "with_geom" was given in the configuration.'

        if 'HAWC2SInputWriter' not in config.keys():
            config['HAWC2SInputWriter'] = {}
            print 'No configuration dictionary given for HAWC2SInputWriter' +\
                  'proceeding with default values.'

        if 'channels' not in config['HAWC2Outputs']:
            raise RuntimeError('You need to supply a list of output channels ' +
                               'in HAWC2Outputs["channels"].')

        if 'no_bins' not in config['HAWC2Outputs'].keys():
            config['HAWC2Outputs']['no_bins'] = 128
        if 'neq' not in config['HAWC2Outputs'].keys():
            config['HAWC2Outputs']['neq'] = 600
        if 'm' not in config['HAWC2Outputs'].keys():
            config['HAWC2Outputs']['m'] = [3, 4, 6, 8, 10, 12]
        if 'stat_list' not in config['HAWC2Outputs'].keys():
            config['HAWC2Outputs']['stat_list'] = ['std', 'rms', 'min', 'int',
                                                   'max', 'range', 'absmax',
                                                   'mean']
        if 'data_directory' not in config['HAWC2InputWriter'].keys():
            config['HAWC2InputWriter']['data_directory'] = 'data'
        if 'turb_directory' not in config['HAWC2InputWriter'].keys():
            config['HAWC2InputWriter']['turb_directory'] = './../turb'


if __name__ == '__main__':

    pass
