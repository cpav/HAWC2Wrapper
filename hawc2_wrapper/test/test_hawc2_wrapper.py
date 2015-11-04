# -*- coding: utf-8 -*-
import unittest
import numpy as np
from openmdao.api import Problem, Group, IndepVarComp
from hawc2_wrapper.hawc2_inputreader import HAWC2InputReader
from hawc2_wrapper.hawc2_inputwriter import HAWC2SInputWriter, HAWC2InputWriter
from hawc2_wrapper.hawc2_wrapper import HAWC2Wrapper
from hawc2_wrapper.hawc2_output import HAWC2SOutput
from hawc2_wrapper.hawc2_aeroelasticsolver import HAWC2SWorkflow, HAWC2SAeroElasticSolver
import filecmp


class TestIO():

    def run(self):
        a = HAWC2InputReader()
        a.htc_master_file = 'main_hs2.htc'
        a.execute()
        self.vartrees = a.vartrees
        b = HAWC2SInputWriter()
        b.vartrees = a.vartrees
        b.case_id = 'h2s'
        b.execute()
        print filecmp.cmp('main_hs2.htc', 'h2s.htc')

        a = HAWC2InputReader()
        a.htc_master_file = 'main_h2.htc'
        a.execute()

        b = HAWC2InputWriter()
        b.vartrees = a.vartrees
        b.turb_directory = a.vartrees.wind.turb_directory
        b.case_id = 'h2'
        b.execute()
        print filecmp.cmp('main_h2.htc', 'h2.htc')


class TestRun():

    def run_hs2(self):

        t2 = HAWC2Wrapper()
        t2.hawc2bin = 'HAWC2s.exe'
        t2.case_id = 'h2s'
        t2.compute()

    def run_h2(self):

        t2 = HAWC2Wrapper()
        t2.hawc2bin = 'HAWC2mb.exe'
        t2.case_id = 'h2'
        t2.compute()

class TestROutput():

    def run(self):

        a = HAWC2InputReader()
        a.htc_master_file = 'main_hs2.htc'
        a.execute()

        b = HAWC2SInputWriter()
        b.vartrees = a.vartrees
        b.case_id = 'h2s'
        b.vartrees.h2s.options.include_torsiondeform = 0
        b.vartrees.h2s.commands = ['compute_optimal_pitch_angle',
                                   'compute_steady_states',
                                   'save_power',
                                   'save_induction']
        b.execute()

        c = HAWC2Wrapper()
        c.hawc2bin = 'HAWC2s.exe'
        c.case_id = b.case_id
        c.compute()

        t3 = HAWC2SOutput()
        t3.case_id = b.case_id
        t3.commands = b.vartrees.h2s.commands
        t3.execute()
        print 'Pitch angle:'
        print t3.pitch
        print 'Radial positions:'
        print t3.s

class TestHAWC2SWorkflow(object):
    
    def run(self):
        
        config = {}
        config['master_file'] = 'main_hs2.htc'
        config['with_structure'] = 0
        config['with_geom'] = 0
        cf = {}
        config['HAWC2SInputWriter'] = cf
        cf = {}
        cf['dry_run'] = False
        cf['copyback_results'] = False
        cf['hawc2bin'] = 'hawc2s.exe'
        config['HAWC2Wrapper'] = cf
        cf = {}
        cf['blade'] = ['aoa', 'Ft', 'Fn', 'cl', 'cd']
        cf['rotor'] = ['wsp', 'pitch', 'P', 'rpm', 'T']
        config['HAWC2SOutputs'] = cf
      
        root = Group()
        #root.add('indep_var', IndepVarComp('x', 7.0))
        #case = {'wsp': [5, 6], 'pitch': [8, 54], 'rpm': [0.001, 1]}
        case = [3, 4 ,5 ,67 ,8]
        root.add('my_comp', HAWC2SWorkflow(config, 'ws_id', case, 0, 0))
        #root.connect('indep_var.x', 'my_comp.x_input')
        root.my_comp.writer.vartrees.h2s.options.include_torsiondeform = 0
        root.my_comp.writer.vartrees.h2s.options.bladedeform = 'nobladedeform'
        prob = Problem(root)
        prob.setup()
        prob.run()
        
        print prob['my_comp.outputs_rotor']
        print prob['my_comp.outputs_blade']
        #print prob['my_comp.pitch']
        #print prob['my_comp.rpm']


class TestHAWC2SAeroElasticSolver(object):
    
    def run(self):
        
        config = {}
        config['master_file'] = 'main_hs2.htc'
        config['with_structure'] = 1
        config['with_geom'] = 1
        #config['aerodynamic_sections'] = 30
        cf = {}
        #config['HAWC2SInputWriter'] = cf
        cf = {}
        cf['blade_ni_span'] = 10
        cf['interp_from_htc'] = False
        cf['hub_radius'] = 3.
        config['HAWC2GeometryBuilder'] = cf
        cf = {}
        cf['dry_run'] = False
        cf['copyback_results'] = False
        cf['hawc2bin'] = 'hawc2s.exe'
        config['HAWC2Wrapper'] = cf

        cf = {}
        cf['wsp'] = [4, 5, 7]
        cf['user'] = [{'wsp': 70, 'pitch': 80, 'rpm': 0.01}]
        config['cases'] = cf

        cf = {}
        cf['blade'] = ['aoa', 'Ft', 'Fn', 'cl', 'cd']
        cf['rotor'] = ['wsp', 'pitch', 'P', 'rpm', 'T']
        config['HAWC2SOutputs'] = cf
        n = 5
        root = HAWC2SAeroElasticSolver(config, [10, 19], n)
        prob = Problem(root)
        
        if config['with_structure']:
            prob.root.add('blade_beam_structure', IndepVarComp('blade_beam_structure', np.zeros([10, 19])), promotes=['*'])

        if config['with_geom']:
            geom_var = ['s', 'x', 'y', 'z', 'rot_x', 'rot_y', 'rot_z',
                        'chord', 'rthick', 'p_le']

            for v in geom_var:
                prob.root.add(v, IndepVarComp(v, np.zeros(n)), promotes=['*'])

            prob.root.add('blade_length', IndepVarComp('blade_length', 0.), promotes=['*'])    

        prob.setup()

        if config['with_structure']:
            mat = np.zeros([10, 19])
            mat[:, 0] = np.linspace(0, 86, 10)
            mat[:, 1] = 1.
            mat[:, 4] = 0.1
            mat[:, 5] = 0.1
            mat[:, 8] = 1e9
            mat[:, 9] = 1e9
            for i in [10, 11, 12, 13, 14, 15]:
                mat[:, i] = 0.1
            prob['blade_beam_structure'] = mat

        if config['with_geom']:
            for v in geom_var[:]:
                prob[v] = np.zeros(n)
            prob['z'] = np.array([0., 0.25, 0.5, 0.75, 1.])
            prob['s'] = np.array([0., 0.25, 0.5, 0.75, 1.])
            prob['chord'] = np.array([1., 2., 2.5, 0.75, 0.5])/86.
            prob['rthick'] = np.array([1., .3, .25, 0.25, 0.2])
            prob['blade_length'] = 86.
        prob.run()
        print prob['aggregate.wsp']
        print prob['aggregate.P']
        print prob['aggregate.Ft']
  
if __name__ == '__main__':

    print 'Test 1'
    t1 = TestIO()
    t1.run()

    t2 = TestRun()
    # t2.run_hs2()

    print 'Test 3'
    t3 = TestROutput()
    #t3.run()
    
    print 'Test 4'
    t4 = TestHAWC2SWorkflow()
    #t4.run()
    
    print 'Test 5'
    t5 = TestHAWC2SAeroElasticSolver()
    t5.run()