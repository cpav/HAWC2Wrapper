# -*- coding: utf-8 -*-
import unittest

from hawc2_wrapper.hawc2_inputreader import HAWC2InputReader
from hawc2_wrapper.hawc2_inputwriter import HAWC2SInputWriter, HAWC2InputWriter

import numpy as np
import copy


def assert_differet_types(self, val, var, vt):

    if type(val) == str:
        self.assertEqual(val, getattr(vt, var))
    elif type(val) == np.ndarray:
        for i in range(len(val)):
            self.assertAlmostEqual(val[i], getattr(vt, var)[i], places=14)
    elif type(val) == list:
        if type(val[0]) == np.ndarray:
            for i in range(len(val)):
                for j in range(len(val[i])):
                    self.assertAlmostEqual(val[i][j], getattr(vt, var)[i][j], places=14)
                    
        else:
            self.assertAlmostEqual(val, getattr(vt, var), places=14)
    else:
        
        self.assertAlmostEqual(val, getattr(vt, var), places=14)


def read_htc(filename):

    read = HAWC2InputReader()
    read.htc_master_file = filename
    read.execute()

    return read


class Test(unittest.TestCase):

    def _test_IO_Aero(self):

        init = read_htc('main_hs2.htc')

        writer = HAWC2SInputWriter()

        for var in init.vartrees.aero.var:
            writer.vartrees = copy.deepcopy(init.vartrees)
            val = round(np.random.rand()*10)
            if var in ['ae_filename']:
                val = './data/DTU_10MW_RWT_ae.dat'
            if var in ['pc_filename']:
                val = './data/DTU_10MW_RWT_pc.dat'
            setattr(writer.vartrees.aero, var, val)
            writer.case_id = 'h2s'
            writer.execute()

            reader = read_htc('h2s.htc')

            assert_differet_types(self, val, var, reader.vartrees.aero)

    def _test_IO_Wind(self):

        init = read_htc('main_hs2.htc')

        writer = HAWC2SInputWriter()

        for var in init.vartrees.wind.var:

            writer.vartrees = copy.deepcopy(init.vartrees)
            val = np.ceil(np.random.rand()*10)

            if var == 'iec_gust':
                val = True
            if var == 'iec_gust_type':
                val = 'abc'
                setattr(writer.vartrees.wind, 'iec_gust', True)

            setattr(writer.vartrees.wind, var, val)

            if var in ['wind_ramp_factor0', 'wind_ramp_factor1']:
                setattr(writer.vartrees.wind, 'wind_ramp_t1', val)
            if var in ['G_A', 'G_phi0', 'G_t0', 'G_T']:
                setattr(writer.vartrees.wind, 'iec_gust', True)
            writer.case_id = 'h2s'
            writer.execute()

            reader = read_htc('h2s.htc')

            assert_differet_types(self, val, var, reader.vartrees.wind)

    def _test_IO_Sim(self):

        init = read_htc('main_h2.htc')

        writer = HAWC2InputWriter()

        for var in init.vartrees.sim.var:

            writer.vartrees = copy.deepcopy(init.vartrees)
            val = np.ceil(np.random.rand()*10)
            if var == 'convergence_limits':
                val = np.ceil(np.random.rand(3)*100)/100
            if var in ['on_no_convergence', 'logfile']:
                val = 'zaq'

            setattr(writer.vartrees.sim, var, val)

            writer.case_id = 'h2s'
            writer.execute()

            reader = read_htc('h2s.htc')

            assert_differet_types(self, val, var, reader.vartrees.sim)

    def _test_IO_Aerodrag(self):

        init = read_htc('main_h2.htc')

        writer = HAWC2InputWriter()

        for var in init.vartrees.aerodrag.elements[0].var:
            writer.vartrees = copy.deepcopy(init.vartrees)
            val = np.ceil(np.random.rand()*10)

            if var == 'sections':
                val = [[0,0, 8],[0,1, 5]]
            setattr(writer.vartrees.aerodrag.elements[0], var, val)

            writer.case_id = 'h2s'
            writer.execute()

            reader = read_htc('h2s.htc')
            assert_differet_types(self, val, var, reader.vartrees.aerodrag.elements[0])

    def _test_IO_Bodies(self):

        init = read_htc('main_h2.htc')

        writer = HAWC2InputWriter()

        for var in init.vartrees.main_bodies.blade1.var:
            print var
            writer.vartrees = copy.deepcopy(init.vartrees)

            val = np.ceil(np.random.rand()*10)
            if var in ['body_name', 'body_type', 'node_distribution']:
                val = 'blade1'
            if var == 'st_input_type':
                val = 0
            if var == 'body_set':
                val = [1, 2]
            if var == 'damping_aniso':
                val = np.ceil(np.random.rand(6)*10)
                writer.vartrees.main_bodies.blade1.damping_type = 'ani'
            if var == 'damping_posdef':
                val = np.ceil(np.random.rand(6)*10)
            if var == 'concentrated_mass':
                val = [np.ceil(np.random.rand(8)*10)]
            setattr(writer.vartrees.main_bodies.blade1, var, val)

            writer.case_id = 'h2s'
            writer.execute()

            reader = read_htc('h2s.htc')
            assert_differet_types(self, val, var, reader.vartrees.main_bodies.blade1)

    def _test_IO_BeamStructure(self):

        init = read_htc('main_h2.htc')

        writer = HAWC2InputWriter()
        writer.vartrees = copy.deepcopy(init.vartrees)

        for iset in range(2):
            for var in writer.vartrees.main_bodies.blade1.beam_structure[iset].var:
                n = len(writer.vartrees.main_bodies.blade1.beam_structure[iset].s)
                writer.vartrees = copy.deepcopy(init.vartrees)

                val = np.random.rand(n)

                setattr(writer.vartrees.main_bodies.blade1.beam_structure[iset], var, val)

                writer.case_id = 'h2s'
                writer.execute()

                reader = read_htc('h2s.htc')
                assert_differet_types(self, val, var, reader.vartrees.main_bodies.blade1.beam_structure[iset])

    def _test_IO_BeamStructure_FPM(self):

        init = read_htc('main_h2_FPM.htc')

        writer = HAWC2InputWriter()
        writer.vartrees = copy.deepcopy(init.vartrees)

        for iset in range(2):
            for var in writer.vartrees.main_bodies.blade1.beam_structure[iset].var:
                n = len(writer.vartrees.main_bodies.blade1.beam_structure[iset].s)
                writer.vartrees = copy.deepcopy(init.vartrees)

                val = np.random.rand(n)

                setattr(writer.vartrees.main_bodies.blade1.beam_structure[iset], var, val)

                writer.case_id = 'h2s'
                writer.execute()

                reader = read_htc('h2s.htc')
                assert_differet_types(self, val, var, reader.vartrees.main_bodies.blade1.beam_structure[iset])

    def test_IO_Controller(self):

        init = read_htc('main_h2.htc')

        writer = HAWC2InputWriter()
        for var in init.vartrees.dlls.risoe_controller.dll_init.init_dic.keys():
            tag = init.vartrees.dlls.risoe_controller.dll_init.init_dic[var][0]
            
            writer.vartrees = copy.deepcopy(init.vartrees)
            val = np.ceil(np.random.rand()*10)

            setattr(writer.vartrees.dlls.risoe_controller.dll_init, tag, val)

            writer.case_id = 'h2s'
            writer.execute()

            reader = read_htc('h2s.htc')
            
            print dir(reader.vartrees.dlls.risoe_controller)
            assert_differet_types(self, val, tag, reader.vartrees.dlls.risoe_controller.dll_init)
            

if __name__ == '__main__':

    unittest.main()
