import copy
import numpy as np

from openmdao.api import Component, Group, ParallelGroup, IndepVarComp

from hawc2_inputreader import HAWC2InputReader
from hawc2_inputwriter import HAWC2SInputWriter
from hawc2_wrapper import HAWC2Wrapper
from hawc2_output import HAWC2SOutputCompact
from hawc2_geomIDO import HAWC2GeometryBuilder


class HAWC2SWorkflow(Component):
    """
    OpenMDAO component to run the HAWC2S workflow.

    parameters
    ----------
    config: dict
        Configuration dictionary. It has to contain the following entries:
        * 'with_structure': bool to add structural properties to the parameters
        * 'with_geom': bool to add the blade geometry to the parameters
        * 'master_file': str with the name of the master file.
        * 'aerodynamic_sections': int of the number of aerodynamic sections.
        * 'HAWC2SOutputs': dict of the outputs required. It has to include two
            dictionaries 'rotor' and 'blade' with the list of rotor sensor and
            blade sensors.

        * 'HAWC2SInputWriter': dict for initialization of HAWC2SInputWriter
            parameters.
        * 'HAWC2Wrapper': dict for initialization of HAWC2Wrapper parameters.
        * 'HAWC2GeometryBuilder': dict for initialization of
            HAWC2GeometryBuilder parameters

    """
    def __init__(self, config, case_id, case, cs_size, pfsize):
        super(HAWC2SWorkflow, self).__init__()
        self.with_structure = config['with_structure']
        self.with_geom = config['with_geom']
        self.with_ctr_tuning = 0
        self.with_freq_placement = 0

        self.reader = HAWC2InputReader(config['master_file'])
        self.reader.execute()

        self.writer = HAWC2SInputWriter(**config['HAWC2SInputWriter'])
        self.writer.vartrees = copy.copy(self.reader.vartrees)
        self.writer.case_id = case_id
        self.writer.vartrees.aero.aerosections = config['aerodynamic_sections']

        nws = self._check_cases(self.writer.vartrees, case_id, case)

        self.wrapper = HAWC2Wrapper(**config['HAWC2Wrapper'])
        self.wrapper.case_id = case_id
        naes = config['aerodynamic_sections']-2

        self.output = HAWC2SOutputCompact(config['HAWC2SOutputs'])
        self.output.case_id = case_id
        self.output.commands = self.reader.vartrees.h2s.commands

        n = len(config['HAWC2SOutputs']['rotor'])
        self.add_output('outputs_rotor', shape=[nws, n])
        n = len(config['HAWC2SOutputs']['blade'])
        self.add_output('outputs_blade', shape=[nws, n*naes])

        if self.with_structure:
            self.add_param('blade_beam_structure', shape=cs_size)

        if self.with_geom:
            self.geom = HAWC2GeometryBuilder(**config['HAWC2GeometryBuilder'])
            self.geom.c12axis_init = self.reader.vartrees.main_bodies.blade1.c12axis.copy()
            self.geom_var = ['s', 'x', 'y', 'z', 'rot_x', 'rot_y', 'rot_z',
                             'chord', 'rthick', 'p_le']
            for v in self.geom_var:
                self.add_param(v, shape=[pfsize])
            self.add_param('blade_length', 0.)

    def solve_nonlinear(self, params, unknowns, resids):

        if self.with_structure:
            body = self.writer.vartrees.main_bodies.blade1
            self._array2hawc2beamstructure(body, params['blade_beam_structure'])

        if self.with_geom:
            for v in self.geom_var:
                setattr(self.geom.bladegeom, v, params[v])
            self.geom.blade_length = params['blade_length']
            self.geom.execute()

            self.writer.vartrees.blade_ae = self.geom.blade_ae
            self.writer.vartrees.main_bodies.blade1.c12axis = self.geom.c12axis

        if self.with_ctr_tuning:
            pass
        if self.with_freq_placement:
            pass

        self.writer.execute()
        self.wrapper.compute()
        self.output.execute()

        unknowns['outputs_rotor'] = self.output.outputs_rotor
        unknowns['outputs_blade'] = self.output.outputs_blade

    def _array2hawc2beamstructure(self, body, body_st):

        bset = body.body_set[1] - 1
        if body.st_input_type is 0:
            body.beam_structure[bset].s = body_st[:, 0]
            body.beam_structure[bset].dm = body_st[:, 1]
            body.beam_structure[bset].x_cg = body_st[:, 2]
            body.beam_structure[bset].y_cg = body_st[:, 3]
            body.beam_structure[bset].ri_x = body_st[:, 4]
            body.beam_structure[bset].ri_y = body_st[:, 5]
            body.beam_structure[bset].x_sh = body_st[:, 6]
            body.beam_structure[bset].y_sh = body_st[:, 7]
            body.beam_structure[bset].E = body_st[:, 8]
            body.beam_structure[bset].G = body_st[:, 9]
            body.beam_structure[bset].I_x = body_st[:, 10]
            body.beam_structure[bset].I_y = body_st[:, 11]
            body.beam_structure[bset].K = body_st[:, 12]
            body.beam_structure[bset].k_x = body_st[:, 13]
            body.beam_structure[bset].k_y = body_st[:, 14]
            body.beam_structure[bset].A = body_st[:, 15]
            body.beam_structure[bset].pitch = body_st[:, 16]
            body.beam_structure[bset].x_e = body_st[:, 17]
            body.beam_structure[bset].y_e = body_st[:, 18]
        else:
            body.beam_structure[bset].s = body_st[:, 0]
            body.beam_structure[bset].dm = body_st[:, 1]
            body.beam_structure[bset].x_cg = body_st[:, 2]
            body.beam_structure[bset].y_cg = body_st[:, 3]
            body.beam_structure[bset].ri_x = body_st[:, 4]
            body.beam_structure[bset].ri_y = body_st[:, 5]
            body.beam_structure[bset].pitch = body_st[:, 6]
            body.beam_structure[bset].x_e = body_st[:, 7]
            body.beam_structure[bset].y_e = body_st[:, 8]
            body.beam_structure[bset].K_11 = body_st[:, 9]
            body.beam_structure[bset].K_12 = body_st[:, 10]
            body.beam_structure[bset].K_13 = body_st[:, 11]
            body.beam_structure[bset].K_14 = body_st[:, 12]
            body.beam_structure[bset].K_15 = body_st[:, 13]
            body.beam_structure[bset].K_16 = body_st[:, 14]
            body.beam_structure[bset].K_22 = body_st[:, 15]
            body.beam_structure[bset].K_23 = body_st[:, 16]
            body.beam_structure[bset].K_24 = body_st[:, 17]
            body.beam_structure[bset].K_25 = body_st[:, 18]
            body.beam_structure[bset].K_26 = body_st[:, 19]
            body.beam_structure[bset].K_33 = body_st[:, 20]
            body.beam_structure[bset].K_34 = body_st[:, 21]
            body.beam_structure[bset].K_35 = body_st[:, 22]
            body.beam_structure[bset].K_36 = body_st[:, 23]
            body.beam_structure[bset].K_44 = body_st[:, 24]
            body.beam_structure[bset].K_45 = body_st[:, 25]
            body.beam_structure[bset].K_46 = body_st[:, 26]
            body.beam_structure[bset].K_55 = body_st[:, 27]
            body.beam_structure[bset].K_56 = body_st[:, 28]
            body.beam_structure[bset].K_66 = body_st[:, 29]

    def _check_cases(self, vt, case_id, case):

        if 'user' not in case_id:
            if isinstance(case, int):
                vt.dlls.risoe_controller.dll_init.Vin = case
                vt.dlls.risoe_controller.dll_init.Vout = case
                vt.dlls.risoe_controller.dll_init.nV = 1
            else:
                vt.dlls.risoe_controller.dll_init.Vin = case[0]
                vt.dlls.risoe_controller.dll_init.Vout = case[-1]
                vt.dlls.risoe_controller.dll_init.nV = len(case)
            nws = vt.dlls.risoe_controller.dll_init.nV

        else:

            vt.h2s.wsp_curve = case['wsp']
            vt.h2s.pitch_curve = case['pitch']
            vt.h2s.rpm_curve = case['rpm']

            if isinstance(vt.h2s.wsp_curve, int):
                nws = 1
            else:
                nws = len(vt.h2s.wsp_curve)

            vt.h2s.commands.remove('compute_optimal_pitch_angle')

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
    """
    def __init__(self, config, nws):
        super(OutputsAggregator, self).__init__()

        self.nws = nws
        self.naes = config['aerodynamic_sections']-2

        # Add parameters coming from ParallelGroup
        n = len(config['HAWC2SOutputs']['rotor'])
        for i in range(nws):
            p_name = 'outputs_rotor_%d' % i
            self.add_param(p_name, shape=[1, n])
        n = len(config['HAWC2SOutputs']['blade'])
        for i in range(nws):
            p_name = 'outputs_blade_%d' % i
            self.add_param(p_name, shape=[1, n*self.naes])

        # Add outputs
        self.sensor_rotor = []
        for sensor in config['HAWC2SOutputs']['rotor']:
            self.sensor_rotor.append(sensor)
            self.add_output(sensor, shape=[nws])

        self.sensor_blade = []
        for sensor in config['HAWC2SOutputs']['blade']:
            self.sensor_blade.append(sensor)
            self.add_output(sensor, shape=[nws, self.naes])

    def solve_nonlinear(self, params, unknowns, resids):

        for j, sensor in enumerate(self.sensor_rotor):
            for i in range(self.nws):
                p_name = 'outputs_rotor_%d' % i
                out = params[p_name]
                unknowns[sensor][i] = out[0, j]
        for j, sensor in enumerate(self.sensor_blade):
            for i in range(self.nws):
                p_name = 'outputs_blade_%d' % i
                out = params[p_name]
                unknowns[sensor][i, :] = out[0, j*self.naes:(j+1)*self.naes]


class HAWC2SAeroElasticSolver(Group):
    """
    Assembly for running HAWC2S in parallel using a CaseIteratorDriver

    parameters
    -----------
    wsp: array-like
        array of wind speeds. pitch and RPM will either be interpolated from
        the opt file or computed on the fly
    htc_master_file: string
        name of HAWC2S master file.
    hawc2bin: string
        absolute path to HAWC2S executable

    optional parameters:
    --------------------
    bladegeom: BladeGeometryVT
        IDOtools blade geometry VarTree
    beamprops: BeamStructureVT
        IDOtools blade beam structural properties
    radius: float
        blade length
    """
    def __init__(self, config, cs_size, pfsize):
        super(HAWC2SAeroElasticSolver, self).__init__()

        # check that the config is ok
        self._check_config(config)

        cases = {}
        cases_list = []
        for ws in config['cases']['wsp']:
            name = 'wsp_%2.2f' % ws
            cases[name.replace('.', '_')] = ws
            cases_list.append(name.replace('.', '_'))
        for icase, case in enumerate(config['cases']['user']):
            name = 'user_%i' % icase
            cases[name] = case
            cases_list.append(name)

        promote = []
        if config['with_structure']:
            promote.append('blade_beam_structure')

        if config['with_geom']:
            var = ['s', 'x', 'y', 'z', 'rot_x', 'rot_y', 'rot_z',
                             'chord', 'rthick', 'p_le']
            for v in var:
                promote.append(v)
            promote.append('blade_length')

        self.add('aggregate', OutputsAggregator(config, len(cases_list)))
        pg = self.add('pg', ParallelGroup(), promotes=promote)

        for i, case_id in enumerate(cases_list):
            pg.add(case_id, HAWC2SWorkflow(config, case_id, cases[case_id],
                                           cs_size, pfsize), promotes=promote)

            self.connect('pg.%s.outputs_rotor' % case_id,
                         'aggregate.outputs_rotor_%d' % i)
            self.connect('pg.%s.outputs_blade' % case_id,
                         'aggregate.outputs_blade_%d' % i)

    def _check_config(self, config):

        if 'master_file' not in config.keys():
            raise RuntimeError('You need to supply the name of the master' +
                               'file in the configuration dictionary.')

        if 'HAWC2GeometryBuilder' not in config.keys():
            raise RuntimeError('You need to supply a config dict' +
                               'for HAWC2GeometryBuilder.')

        if 'HAWC2Wrapper' not in config.keys():
            raise RuntimeError('You need to supply a config dict' +
                               'for HAWC2Wrapper.')

        if 'aerodynamic_sections' not in config.keys():
            config['aerodynamic_sections'] = 40
            print 'The number of aerodynamic section has not been specified' +\
                  ' in the configuration dictionary. Default value of 40 selected.'

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

        if 'HAWC2SOutputs' not in config.keys():
            print 'No HAWC2SOutputs output dictionary given proceeding' +\
                  ' with default values.'
            config['HAWC2SOutputs'] = {}
            config['HAWC2SOutputs']['rotor'] = ['wsp', 'pitch', 'P', 'T']
            config['HAWC2SOutputs']['blade'] = ['aoa', 'cl', 'Fn']
#
#        self.connect('designTSR',
#                     'casegen.vartrees.dlls.risoe_controller.dll_init.designTSR')
#        self.connect('designTSR',
#                     'vartrees_out.dlls.risoe_controller.dll_init.designTSR')
#
#    def configure_freq_placement(self, freq_type='ae'):
#
#        self.h2.configure_freq_placement_cid(freq_type=freq_type)
#        self.connect('h2.freq_factor', 'h2post.freq_factor_cid')
#        self.create_passthrough('h2post.freq_factor')
#        self.create_passthrough('h2.mode_freq')
#        self.create_passthrough('h2.mode_damp')
#        self.create_passthrough('h2.mode_target_freq')
#        self.create_passthrough('h2.mode_target_damp')
#
#    def configure_controller_tuning(self):
#
#        self.controller = self.reader.vartrees.dlls.risoe_controller.dll_init.copy()
#        for att in self.controller.list_vars():
#            if att == 'designTSR':
#                continue
#            self.connect('controller.'+att,
#                         'casegen.vartrees.dlls.risoe_controller.dll_init.'+att)
#            self.connect('controller.'+att,
#                         'vartrees_out.dlls.risoe_controller.dll_init.'+att)
#
#    """

if __name__ == '__main__':

    pass
