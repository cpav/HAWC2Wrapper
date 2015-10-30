import copy
import numpy as np

from openmdao.api import Component, Group, ParallelGroup, IndepVarComp

from hawc2_inputreader import HAWC2InputReader
from hawc2_inputwriter import HAWC2SInputWriter
from hawc2_wrapper import HAWC2Wrapper
from hawc2_output import HAWC2SOutput
from hawc2_geomIDO import HAWC2GeometryBuilder
#class MakeCases(Component):
#    """
#    build case vartrees
#
#    inputs
#    -------
#
#    wsp: array-like
#        array of wind speeds for which to compute loads under normal operational conditions
#    user_cases: list of dict
#        list of off design cases e.g. stand still extreme wind conditions. the syntax is:
#        {"wsp": 70, "pitch": 0., "rpm": 0.001}
#
#    outputs
#    --------
#    cases: list of vartrees
#    case_ids: list of case_ids
#    """
#
#    vartrees = VarTree(HAWC2VarTrees(), iotype='in')
#    wsp = Array(iotype='in', desc='array of wind speeds for which to compute the power curve')
#    user_cases = List(iotype='in', desc='List of user defined off-design cases containing'
#                                        'the format is a list of dictionaries with e.g.:'
#                                        '{"wsp": 70, "pitch": 0., "rpm": 0.001}')
#
#    cases = List(iotype='out')
#    case_ids = List(iotype='out')
#
#    def execute(self):
#
#        self.cases = []
#        self.case_ids = []
#        for w in self.wsp:
#            vt = self.vartrees.copy()
#            vt.h2s.wsp_cases = [w]
#            self.cases.append(vt)
#            self.case_ids.append('wsp_%2.2f' % w)
#
#        for i, case in enumerate(self.user_cases):
#            vt = self.vartrees.copy()
#            vt.h2s.cases = [case]
#            vt.h2s.wsp_cases = []
#            if 'pitch' in case:
#                try:
#                    vt.h2s.commands.remove('compute_optimal_pitch_angle')
#                except:
#                    pass
#
#            if case['rpm'] < 0.1:
#                try:
#                    vt.h2s.commands.remove('compute_optimal_pitch_angle')
#                except:
#                    pass
#                try:
#                    vt.h2s.commands.remove('compute_stability_analysis')
#                except:
#                    pass
#                try:
#                    vt.h2s.commands.remove('compute_structural_modal_analysis')
#                except:
#                    pass
#                vt.h2s.options.induction = 'noinduction'
#                vt.h2s.options.tipcorrect = 'notipcorrect'
#            self.cases.append(vt)
#            self.case_ids.append('user_%i' % i)
#

class HAWC2SWorkflow(Component):

    def __init__(self, config, case_id, cs_size, pfsize):
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
        
        self.wrapper = HAWC2Wrapper(**config['HAWC2Wrapper'])
        self.wrapper.case_id = case_id
        ns = self.writer.vartrees.aero.aerosections-2
        nws = self.writer.vartrees.dlls.risoe_controller.dll_init.nV
        self.output = HAWC2SOutput()
        self.output.case_id = case_id
        self.output.commands = self.reader.vartrees.h2s.commands
        for name in self.output.outlist1:
            self.add_output(name, np.zeros(nws))
        for name in self.output.outlist2:
            self.add_output(name, np.zeros([nws, ns]))
            
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

        #self.wrapper.case_id = case_id
        self.wrapper.compute()

        self.output.execute()

        for name in self.output.outlist1:
            print name
            unknowns[name] = getattr(self.output, name)
        for name in self.output.outlist2:
            print name
            unknowns[name] = getattr(self.output, name)

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

    def array2hawc2bladegeometry(blade_ae, blade_geom):

            blade_ae.s = blade_geom[:, 0]
            blade_ae.chord = blade_geom[:, 1]
            blade_ae.rthick = blade_geom[:, 2]

    def check_options(vt):
            vt.h2s.cases = [case]
            vt.h2s.wsp_cases = []
            if 'pitch' in case:
                try:
                    vt.h2s.commands.remove('compute_optimal_pitch_angle')
                except:
                    pass

            if case['rpm'] < 0.1:
                try:
                    vt.h2s.commands.remove('compute_optimal_pitch_angle')
                except:
                    pass
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
    def __init__(self, config, vartrees):
        super(HAWC2SAeroElasticSolver, self).__init__()



        cid = self.add('cid', ParallelGroup())

        for w in self.cases:
            case_id = 'wsp_%2.2f' % w
            cid.add(case_id, HAWC2SWorkflow(config, case_id, cs_size, pfsize))
            if config['with_structure']:
                self.add('cs_props_c', IndepVarComp('cs_props', cs_props), promotes=['*'])
            if config['with_structure']:
                self.add('blade_length_c', IndepVarComp('blade_length', blade_length), promotes=['*'])
                self.add('blade_span_ni_c', IndepVarComp('blade_span_ni', blade_span_ni), promotes=['*'])
                self.add('pfIn_c', IndepVarComp('pfIn', pfIn), promotes=['*'])

            
            # create connections



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