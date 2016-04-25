import copy
import numpy as np
import os
import shutil

from openmdao.api import Component, Group, ParallelGroup

from hawc2_inputreader import HAWC2InputReader
from hawc2_inputwriter import HAWC2SInputWriter
from hawc2_wrapper import HAWC2Wrapper
from hawc2_output import HAWC2SOutputCompact, FreqDampTarget
from hawc2_geometry import HAWC2GeometryBuilder


class HAWC2SWorkflow(Component):
    """OpenMDAO component to run the HAWC2S workflow.

    Parameters
    ----------
    config: dict
        Configuration dictionary. It has to contain the following entries:

        * 'with_structure': bool to add structural properties to the parameters

        * 'with_geom': bool to add the blade geometry to the parameters

        * 'master_file': str with the name of the master file.

        * 'aerodynamic_sections': int of the number of aerodynamic sections.

        * 'HAWC2SOutputs': dict of the outputs required. It has to include two\
            dictionaries 'rotor' and 'blade' with the list of rotor sensor and\
            blade sensors.

        * 'HAWC2SInputWriter': dict for initialization of HAWC2SInputWriter\
            parameters.

        * 'HAWC2Wrapper': dict for initialization of HAWC2Wrapper parameters.

        * 'HAWC2GeometryBuilder': dict for initialization of \
          HAWC2GeometryBuilder parameters

    case_id: str
        Name of the HAWC2s case to create and run. If case_id contains the word
        "user" then it is assumed that the user is providing the operational
        points. Therefore "case" has to be a dictionary.

    case: list/dict
        This variable can be a list or a dictionary. If the operational points
        need to be computed "case" is a list with the wind speeds to evaluate.
        If more than one wind speed is given the wind speeds evaluated are
        actually the first, the last, and those obtained with a linspace of the
        same length. If the operational points do not need to be computed
        "case" is a dictionary with keys "wsp", "rpm", and "pitch". Each
        parameter is a list containing wind speed, rotorspeed and pitch angle
        respectively.

    Returns
    -------

    """
    def __init__(self, config, case_id, case):
        super(HAWC2SWorkflow, self).__init__()

        self.basedir = os.getcwd()
        self.keep_work_dirs = False

        self.with_structure = config['with_structure']
        self.with_geom = config['with_geom']
        self.with_tsr = config['with_tsr']
        self.with_ctr_tuning = 0
        if 'FreqDampTarget' in config.keys():
            self.with_freq_placement = True
        else:
            self.with_freq_placement = False

        self.reader = HAWC2InputReader(config['master_file'])
        self.reader.execute()

        self.writer = HAWC2SInputWriter(**config['HAWC2SInputWriter'])
        self.writer.data_directory = config['HAWC2SInputWriter']['data_directory']
        self.writer.vartrees = copy.copy(self.reader.vartrees)
        self.writer.case_id = case_id
        self.writer.vartrees.aero.ae_filename = \
            os.path.join(self.writer.data_directory, case_id+'_ae.dat')
        self.writer.vartrees.aero.pc_filename = \
            os.path.join(self.writer.data_directory, case_id+'_pc.dat')
        self.writer.vartrees.aero.aerosections = config['aerodynamic_sections']

        nws = self._check_cases(self.writer.vartrees, case_id, case)

        self.wrapper = HAWC2Wrapper(**config['HAWC2Wrapper'])
        self.wrapper.case_id = case_id
        naes = config['aerodynamic_sections'] - 2
        nsts = config['structural_sections']

        self.output = HAWC2SOutputCompact(config['HAWC2SOutputs'])
        self.output.case_id = case_id
        self.output.commands = self.reader.vartrees.h2s.commands

        n = len(config['HAWC2SOutputs']['rotor'])
        self.add_output('outputs_rotor', shape=[nws, n])
        n = len(config['HAWC2SOutputs']['blade'])
        self.add_output('outputs_blade', shape=[nws, n*naes])
        self.add_output('outputs_blade_fext', shape=[nws, 6*nsts])

        if self.with_tsr:
            self.add_param('tsr', 0.)

        if self.with_structure:
            if config['hawc2_FPM']:
                cssize = (config['structural_sections'], 30)
            else:
                cssize = (config['structural_sections'], 19)
            self.add_param('blade_beam_structure', shape=cssize)

        self.geom = HAWC2GeometryBuilder(config['structural_sections'], **config['BladeGeometryBuilder'])
        self.geom.c12axis_init = self.reader.vartrees.main_bodies.blade1.c12axis.copy()
        self.geom.c12axis_init[:, :3] /= self.geom.c12axis_init[-1, 2]
        if self.with_geom:
            self.geom_var = ['s', 'x', 'y', 'z', 'rot_x', 'rot_y', 'rot_z',
                             'chord', 'rthick', 'p_le']
            for v in self.geom_var:
                self.add_param(v, shape=(config['aerodynamic_sections']))
            self.add_param('blade_length', 0.)
            self.geom.interp_from_htc = False
        else:
            self.geom.interp_from_htc = True

        if self.with_freq_placement:
            freq_nf = config['FreqDampTarget']['mode_freq'].shape[1]-1

            self.freq = FreqDampTarget(**config['FreqDampTarget'])
            self.add_output('freq_factor', shape=[nws, 2*freq_nf])

    def solve_nonlinear(self, params, unknowns, resids):

        workdir = 'hawc2s_model_%s_%i' % (self.writer.case_id, self.__hash__())

        try:
            os.mkdir(workdir)
        except:
            pass
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

        if self.with_ctr_tuning:
            pass


        self.writer.execute()
        self.wrapper.compute()
        if self.wrapper.success:
            self.output.execute()

        try:
            unknowns['outputs_rotor'] = self.output.outputs_rotor
        except:
            pass
        try:
            unknowns['outputs_blade'] = self.output.outputs_blade
        except:
            pass
        try:
            unknowns['outputs_blade_fext'] = self.output.outputs_blade_fext
        except:
            pass

        if self.with_freq_placement:
            self.freq.freqdamp = self.output.aeroelasticfreqdamp
            self.freq.execute()
            unknowns['freq_factor'] = self.freq.freq_factor

        os.chdir(self.basedir)

    def _check_cases(self, vt, case_id, case):

        if 'user' not in case_id:
            if isinstance(case, int) or isinstance(case, float):
                vt.dlls.risoe_controller.dll_init.Vin = case
                vt.dlls.risoe_controller.dll_init.Vout = case
                vt.dlls.risoe_controller.dll_init.nV = 1
                vt.h2s.wsp_cases = np.asarray([case])
            else:
                vt.dlls.risoe_controller.dll_init.Vin = case[0]
                vt.dlls.risoe_controller.dll_init.Vout = case[-1]
                vt.dlls.risoe_controller.dll_init.nV = len(case)
                vt.h2s.wsp_cases = np.asarray(case)
            nws = vt.dlls.risoe_controller.dll_init.nV

        else:

            vt.h2s.wsp_curve = case['wsp']
            vt.h2s.pitch_curve = case['pitch']
            vt.h2s.rpm_curve = case['rpm']

            if isinstance(vt.h2s.wsp_curve, (int, float)):
                nws = 1
            else:
                nws = len(vt.h2s.wsp_curve)

            try:
                vt.h2s.commands.remove('compute_optimal_pitch_angle')
            except:
                pass

            if np.min(vt.h2s.rpm_curve) < 0.1:
                try:
                    vt.h2s.commands.remove('compute_stability_analysis')
                except:
                    pass
                try:
                    vt.h2s.commands.remove('compute_structural_modal_analysis')
                except:
                    pass
                vt.h2s.options.induction = 'noinduction'
                vt.h2s.options.tipcorrect = 'notipcorrect'
        return nws


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
        self.nsts = config['structural_sections']

        # Add parameters coming from ParallelGroup
        n = len(config['HAWC2SOutputs']['rotor'])
        for i in range(n_cases):
            p_name = 'outputs_rotor_%d' % i
            self.add_param(p_name, shape=[1, n])
        n = len(config['HAWC2SOutputs']['blade'])
        for i in range(n_cases):
            p_name = 'outputs_blade_%d' % i
            self.add_param(p_name, shape=[1, n*self.naes])
        for i in range(n_cases):
            p_name = 'outputs_blade_fext_%d' % i
            self.add_param(p_name, shape=[1, 6*self.nsts])

        # Add outputs
        self.sensor_rotor = []
        for sensor in config['HAWC2SOutputs']['rotor']:
            self.sensor_rotor.append(sensor)
            self.add_output(sensor, shape=[n_cases])

        self.sensor_blade = []
        for sensor in config['HAWC2SOutputs']['blade']:
            self.sensor_blade.append(sensor)
            self.add_output(sensor, shape=[n_cases, self.naes])

        self.sensor_fext_blade = []
        if config['with_sectional_forces']:
            for sensor in ['Fx_e', 'Fy_e', 'Fz_e', 'Mx_e', 'My_e', 'Mz_e']:
                self.sensor_fext_blade.append(sensor)
                self.add_output(sensor, shape=[n_cases, self.nsts])

        if 'FreqDampTarget' in config.keys():
            self.with_freq_placement = True
        else:
            self.with_freq_placement = False

        if self.with_freq_placement:

            freq_nf = config['FreqDampTarget']['mode_freq'].shape[1]-1

            for i in range(n_cases):
                p_name = 'freq_factor_%d' % i
                self.add_param(p_name, shape=[1, 2*freq_nf])

            self.add_output('freq_factor', shape=[n_cases, 2*freq_nf])


    def solve_nonlinear(self, params, unknowns, resids):

        for j, sensor in enumerate(self.sensor_rotor):
            for i in range(self.n_cases):
                p_name = 'outputs_rotor_%d' % i
                out = params[p_name]
                unknowns[sensor][i] = out[0, j]
        for j, sensor in enumerate(self.sensor_blade):
            for i in range(self.n_cases):
                p_name = 'outputs_blade_%d' % i
                out = params[p_name]
                unknowns[sensor][i, :] = out[0, j*self.naes:(j+1)*self.naes]
        for j, sensor in enumerate(self.sensor_fext_blade):
            for i in range(self.n_cases):
                p_name = 'outputs_blade_fext_%d' % i
                out = params[p_name]
                unknowns[sensor][i, :] = out[0, j*self.nsts:(j+1)*self.nsts]

        if self.with_freq_placement:
            for i in range(self.n_cases):
                p_name = 'freq_factor_%d' % i
                unknowns['freq_factor'][i, :] = params[p_name][0, :]

        fid = open('Rotor_loads.dat', 'w')
        fmt = '#'+'%23s '+(len(self.sensor_rotor)-1)*'%24s '+'\n'
        fid.write(fmt % tuple(self.sensor_rotor))
        for i in range(self.n_cases):
            for j, sensor in enumerate(self.sensor_rotor):
                fid.write('%24.15e ' % unknowns[sensor][i])
            fid.write('\n')
        fid.close()


class HAWC2SAeroElasticSolver(Group):
    """
    OpenMDAO group to execute HAWC2s in parallel.

    parameters
    -----------
    config: dict
        Configuration dictionary.

    """
    def __init__(self, config, cssize=None, pfsize=None):
        super(HAWC2SAeroElasticSolver, self).__init__()

        if cssize != None:
             print 'Warning: cssize should be set in config["structural_sections"]'
             config['structural_sections'] = cssize

        if pfsize != None:
             print 'Warning: pfsize should be set in config["aerodynamic_sections"]'
             config['aerodynamic_sections'] = pfsize

        # check that the config is ok
        self._check_config(config)

        cases = {}
        cases_list = []
        for ws in config['cases']['wsp']:
            name = 'wsp_%2.2f' % ws
            cases[name.replace('.', '_')] = ws
            cases_list.append(name.replace('.', '_'))

        if not isinstance(config['cases']['user'], list):
            print "Error! config['cases']['user'] has to be a list of dictionaries"
        for icase, case in enumerate(config['cases']['user']):
            name = 'user_%i' % icase
            cases[name] = case
            cases_list.append(name)

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
        # these should be converted to FUSED-Wind variables
        # but we'll just promote them for now
        promotions = config['HAWC2SOutputs']['blade'] + config['HAWC2SOutputs']['rotor']
        if config['with_sectional_forces']:
            promotions += ['Fx_e', 'Fy_e', 'Fz_e', 'Mx_e', 'My_e', 'Mz_e']
        if 'FreqDampTarget' in config.keys():
            promotions += ['freq_factor']
        self.add('aggregate', OutputsAggregator(config, len(cases_list)), promotes=promotions)
        pg = self.add('pg', ParallelGroup(), promotes=promote)

        for i, case_id in enumerate(cases_list):
            pg.add(case_id, HAWC2SWorkflow(config, case_id, cases[case_id]),
                                           promotes=promote)

            self.connect('pg.%s.outputs_rotor' % case_id,
                         'aggregate.outputs_rotor_%d' % i)
            self.connect('pg.%s.outputs_blade' % case_id,
                         'aggregate.outputs_blade_%d' % i)
            self.connect('pg.%s.outputs_blade_fext' % case_id,
                         'aggregate.outputs_blade_fext_%d' % i)
            if 'FreqDampTarget' in config.keys():
                self.connect('pg.%s.freq_factor' % case_id,
                             'aggregate.freq_factor_%d' % i)

    def _check_config(self, config):

        if 'master_file' not in config.keys():
            raise RuntimeError('You need to supply the name of the master' +
                               'file in the configuration dictionary.')

        if 'BladeGeometryBuilder' not in config.keys():
            config['BladeGeometryBuilder'] = {}

        if 'HAWC2Wrapper' not in config.keys():
            raise RuntimeError('You need to supply a config dict' +
                               'for HAWC2Wrapper.')

        if 'aerodynamic_sections' not in config.keys():
            raise RuntimeError('You need to supply the aerodynamic_sections' +
                               ' parameter for HAWC2Wrapper.')

        if 'structural_sections' not in config.keys():
            raise RuntimeError('You need to supply the structural_sections' +
                               ' parameter for HAWC2Wrapper.')

        if 'with_tsr' not in config.keys():
            config['with_tsr'] = False
            print 'Tip-speed-ratio is not added as parameter because' +\
                  ' the option "with_tsr" was not provided in the configuration.'

        if 'with_structure' not in config.keys():
            config['with_structure'] = False
            print 'Structural properties are not set as parameters because' +\
                  ' no option "with_structure" was given in the configuration.'

        if 'with_sectional_forces' not in config.keys():
            config['with_sectional_forces'] = False
            print 'config parameter "with_sectional_forces" not set' +\
                  'set to True to output sectional forces (fext files)'

        if 'hawc2_FPM' not in config.keys():
            config['hawc2_FPM'] = False
            print 'hawc2_FPM not supplied: structural properties assumed to be ' +\
                  'in standard HAWC2 format.'

        if 'with_geom' not in config.keys():
            config['with_geom'] = False
            print 'Blade planform properties are not set as parameters because' +\
                  ' no option "with_geom" was given in the configuration.'

        if 'HAWC2SInputWriter' not in config.keys():
            config['HAWC2SInputWriter'] = {}
            print 'No configuration dictionary given for HAWC2SInputWriter' +\
                  'proceeding with default values.'

        if 'HAWC2SOutputs' not in config.keys():
            print 'No HAWC2SOutputs output dictionary given proceeding' +\
                  ' with default values.'
            config['HAWC2SOutputs'] = {}
            config['HAWC2SOutputs']['rotor'] = ['wsp', 'pitch', 'P', 'T']
            config['HAWC2SOutputs']['blade'] = ['aoa', 'cl', 'Fn']

        if 'data_directory' not in config['HAWC2SInputWriter'].keys():
            config['HAWC2SInputWriter']['data_directory'] = 'data'

if __name__ == '__main__':

    pass
