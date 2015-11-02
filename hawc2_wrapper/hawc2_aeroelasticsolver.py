import copy
import numpy as np

from openmdao.api import Component, Group, ParallelGroup, IndepVarComp

from hawc2_inputreader import HAWC2InputReader
from hawc2_inputwriter import HAWC2SInputWriter
from hawc2_wrapper import HAWC2Wrapper
from hawc2_output import HAWC2SOutputCompact
from hawc2_geomIDO import HAWC2GeometryBuilder


class HAWC2SWorkflow(Component):

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

        nws = self._check_cases(self.writer.vartrees, case_id, case)

        self.wrapper = HAWC2Wrapper(**config['HAWC2Wrapper'])
        self.wrapper.case_id = case_id
        naes = self.writer.vartrees.aero.aerosections-2

        self.output = HAWC2SOutputCompact(config['HAWC2SOutputs'])
        self.output.case_id = case_id
        self.output.commands = self.reader.vartrees.h2s.commands

        n = len(config['HAWC2SOutputs']['rotor'])
        self.add_output('outputs_rotor', np.zeros([nws, n]))
        n = len(config['HAWC2SOutputs']['blade'])
        self.add_output('outputs_blade', np.zeros([nws, n*naes]))

        if self.with_structure:
            self.add_param('cs_props', np.zeros(cs_size))

        if self.with_geom:
            self.geom = HAWC2GeometryBuilder(**config['HAWC2GeometryBuilder'])
            self.geom.c12axis_init = self.reader.vartrees.main_bodies.blade1.c12axis.copy()

            self.add_param('pfIn', np.zeros([pfsize, 11]))
            self.add_param('blade_length', 0.)
            self.add_param('blade_ni_span', 0)


    def solve_nonlinear(self, params, unknowns, resids):

        if self.with_structure:
            body = self.writer.vartrees.main_bodies.blade1
            self._array2hawc2beamstructure(body, params['cs_props'])

        if self.with_geom:
            self.geom.blade_length = params['blade_length']
            self.geom.blade_ni_span = params['blade_ni_span']
            bladegeom = self.geom.bladegeom
            self._array2bladegeometry(bladegeom, params['pfIn'])
            self.geom.execute()
            self.writer.vartrees.blade_ae = self.geom.blade_ae
            self.writer.vartrees.mainbodies.blade1.c12axis = self.geom.c12axis

        if self.with_ctr_tuning:
            pass
        if self.with_freq_placement:
            pass

        self.writer.execute()

        self.wrapper.compute()

        self.output.execute()

        unknowns['outputs_rotor'] = self.output.outputs_rotor
        unknowns['outputs_blade'] = self.output.outputs_blade

    def _array2bladegeometry(bladegeom, pfIn):

        bladegeom.s = pfIn[:, 0]

        bladegeom.chord = pfIn[:, 1]
        bladegeom.rthick = pfIn[:, 2]
        bladegeom.athick = pfIn[:, 3]
        bladegeom.p_le = pfIn[:, 4]

        bladegeom.main_axis = pfIn[:, 5:8]
        bladegeom.rot_x = pfIn[:, 8]
        bladegeom.rot_y = pfIn[:, 9]
        bladegeom.rot_z = pfIn[:, 10]

    def _array2hawc2beamstructure(body, body_st):

        bset = body.body_set[1]
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
    def __init__(self, config, nws, naes):
        super(OutputsAggregator, self).__init__()

        self.nws = nws
        self.naes = naes
        self.sensor_rotor = []
        for sensor in config['HAWC2SOutputs']['rotor']:
            self.sensor_rotor.append(sensor)
            self.add_output(sensor, np.array([nws, 1]))

        self.sensor_blade = []
        for sensor in config['HAWC2SOutputs']['blade']:
            self.sensor_blade.append(sensor)
            self.add_output(sensor, np.array([nws, naes]))

        n = len(self.sensor_rotor)
        for i in range(nws):
            self.add_param('outputs_rotor_%i' % i, np.array([1, n]))
        n = len(self.sensor_blade)
        for i in range(nws):
            self.add_param('outputs_blade_%i' % i, np.array([1, n*naes]))

    def solve_nonlinear(self, params, unknowns, resids):

        for j, sensor in enumerate(self.sensor_rotor):
            for i in range(self.nws):
                out = params['%s_%i' % (sensor, i)]
                unknowns[sensor][i, 0] = out[:, j]
        for j, sensor in enumerate(self.sensor_blade):
            for i in range(self.nws):
                out = params['%s_%i' % (sensor, i)]
                unknowns[sensor][i, :] = out[:, j*self.naes, (j+1)*self.naes]


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
    def __init__(self, config, cs_size, pfsize, naes):
        super(HAWC2SAeroElasticSolver, self).__init__()

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

        self.add('aggregate', OutputsAggregator(config, len(cases_list), naes))
        pg = self.add('pg', ParallelGroup())

        for i, case_id in enumerate(cases_list):
            pg.add(case_id, HAWC2SWorkflow(config, case_id, cases[case_id], cs_size, pfsize))

            self.connect('pg.%s.outputs_rotor' % case_id,
                         'aggregate.outputs_rotor_%i' % i)
            self.connect('pg.%s.outputs_blade' % case_id,
                         'aggregate.outputs_blade_%i' % i)

#
#    """
#    def configure_hawc2s(self, htc_master_file=''):
#
#        # Generate simple CID cases
#        self.add('casegen', MakeCases())
#        self.driver.workflow.add('casegen')
#        self.create_passthrough('casegen.user_cases')
#
#        self.htc_master_file = htc_master_file
#
#        if not self.htc_master_file == '':
#            # Read Input file for data initialization
#            self.reader = HAWC2InputReader()
#            self.reader.htc_master_file = self.htc_master_file
#            self.reader.execute()
#
#            self.casegen.vartrees = self.reader.vartrees.copy()
#            self.vartrees_out = self.reader.vartrees.copy()
#
#        # connect FUSED-Wind inflow variables used by HAWC2S
#        self.connect('inflow.vhub', 'casegen.wsp')
#        self.connect('inflow.density[0]', 'casegen.vartrees.wind.density')
#
#        self.connect('designTSR',
#                     'casegen.vartrees.dlls.risoe_controller.dll_init.designTSR')
#        self.connect('designTSR',
#                     'vartrees_out.dlls.risoe_controller.dll_init.designTSR')
#
#        # add case iterator with HAWC2 wrapper
#        self.add('h2', HAWC2SCaseIter())
#        self.driver.workflow.add('h2')
#        self.connect('model_name', 'h2.model_name')
#        self.connect('hawc2bin', 'h2.hawc2bin')
#        self.connect('casegen.cases', 'h2.vartrees')
#        self.connect('casegen.case_ids', 'h2.case_ids')
#
#        self.create_passthrough('h2.sequential')
#        self.create_passthrough('h2.set_tsr_flag')
#        self.h2.output.commands = self.reader.vartrees.h2s.commands
#
#        # postprocess CID cases
#        self.add('h2post', H2SCIDPostProcess())
#        self.driver.workflow.add('h2post')
#        self.connect('h2.rotor_loads', 'h2post.rotor_loads_cid')
#        self.connect('h2.blade_loads', 'h2post.blade_loads_cid')
#        self.connect('h2.hub_loads', 'h2post.hub_loads_cid')
#        self.connect('h2.blade_disps', 'h2post.blade_disps_cid')
#        self.connect('h2.oper', 'h2post.oper_cid')
#
#        self.connect('h2post.rotor_loads', 'rotor_loads')
#        self.connect('h2post.blade_loads', 'blade_loads')
#        self.connect('h2post.blade_disps', 'blade_disps')
#        self.connect('h2post.oper', 'oper')
#
#        self.create_passthrough('h2post.hub_loads')
#
#        self.log_level = logging.DEBUG
#
#
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

    config = {}
    config['master_file'] = 'main_hs2.htc'
    cf = {}
    config['HAWC2SInputWriter'] = cf
    cf = {}
    cf['dry_run'] = False
    cf['copyback_results'] = False
    cf['hawc2bin'] = 'hawc2s.exe'
    config['HAWC2Wrapper'] = cf




    a = HAWC2SWorkflow(config, 0, 0, 0, 0)
    a.solve_nonlinear()