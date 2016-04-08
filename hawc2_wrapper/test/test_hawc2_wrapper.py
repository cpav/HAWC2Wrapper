# -*- coding: utf-8 -*-
import unittest

from hawc2_wrapper.hawc2_inputreader import HAWC2InputReader
from hawc2_wrapper.hawc2_inputwriter import HAWC2SInputWriter, HAWC2InputWriter
from hawc2_wrapper.hawc2_wrapper import HAWC2Wrapper
from hawc2_wrapper.hawc2_output import HAWC2OutputBase, HAWC2SOutputBase, FreqDampTarget

from wetb.prepost import dlcdefs

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
                    self.assertAlmostEqual(val[i][j],
                                           getattr(vt, var)[i][j], places=14)

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

    def compare_lists(self, a, b):

        for val, ref in zip(a, b):

            if isinstance(val, list):
                val = val[0]
            if abs(val) < 1e-6:
                val = 0.
                ref = 0.

            if ref == 0.0:
                self.assertAlmostEqual(val, ref)
            elif ref is np.nan:
                self.assertIs(val, np.nan)
            else:
                self.assertLess(abs(val/ref - 1.), 1e-4)

    def test_IO_Aero(self):

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

    def test_IO_Wind(self):

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

    def test_IO_Sim(self):

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

    def test_IO_Aerodrag(self):

        init = read_htc('main_h2.htc')

        writer = HAWC2InputWriter()

        for var in init.vartrees.aerodrag.elements[0].var:
            writer.vartrees = copy.deepcopy(init.vartrees)
            val = np.ceil(np.random.rand()*10)

            if var == 'sections':
                val = [[0, 0, 8], [0, 1, 5]]
            setattr(writer.vartrees.aerodrag.elements[0], var, val)

            writer.case_id = 'h2s'
            writer.execute()

            reader = read_htc('h2s.htc')
            assert_differet_types(self, val, var,
                                  reader.vartrees.aerodrag.elements[0])

    def test_IO_Bodies(self):

        init = read_htc('main_h2.htc')

        writer = HAWC2InputWriter()

        for var in init.vartrees.main_bodies.blade1.var:
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
            assert_differet_types(self, val, var,
                                  reader.vartrees.main_bodies.blade1)

    def test_IO_BeamStructure(self):

        init = read_htc('main_h2.htc')
        ivt = init.vartrees
        writer = HAWC2InputWriter()
        writer.vartrees = copy.deepcopy(ivt)

        for iset in range(2):
            for var in ivt.main_bodies.blade1.beam_structure[iset].var:
                n = len(ivt.main_bodies.blade1.beam_structure[iset].s)
                writer.vartrees = copy.deepcopy(ivt)

                val = np.random.rand(n)

                setattr(writer.vartrees.main_bodies.blade1.beam_structure[iset],
                        var, val)

                writer.case_id = 'h2s'
                writer.execute()

                reader = read_htc('h2s.htc')
                assert_differet_types(self, val, var,
                                      reader.vartrees.main_bodies.blade1.beam_structure[iset])

    def test_IO_BeamStructure_FPM(self):

        init = read_htc('main_h2_FPM.htc')
        ivt = init.vartrees
        writer = HAWC2InputWriter()
        writer.vartrees = copy.deepcopy(ivt)

        for iset in range(2):
            for var in ivt.main_bodies.blade1.beam_structure[iset].var:
                n = len(ivt.main_bodies.blade1.beam_structure[iset].s)
                writer.vartrees = copy.deepcopy(ivt)

                val = np.random.rand(n)

                setattr(writer.vartrees.main_bodies.blade1.beam_structure[iset],
                        var, val)

                writer.case_id = 'h2s'
                writer.execute()

                reader = read_htc('h2s.htc')
                assert_differet_types(self, val, var,
                                      reader.vartrees.main_bodies.blade1.beam_structure[iset])

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

            assert_differet_types(self, val, tag,
                                  reader.vartrees.dlls.risoe_controller.dll_init)

    def test_IO_Servo(self):

        init = read_htc('main_h2.htc')
        ivt = init.vartrees
        writer = HAWC2InputWriter()
        for var in ivt.dlls.servo_with_limits.dll_init.init_dic.keys():
            tag = ivt.dlls.servo_with_limits.dll_init.init_dic[var][0]

            writer.vartrees = copy.deepcopy(ivt)
            val = np.ceil(np.random.rand()*10)

            setattr(writer.vartrees.dlls.servo_with_limits.dll_init, tag, val)

            writer.case_id = 'h2s'
            writer.execute()

            reader = read_htc('h2s.htc')

            assert_differet_types(self, val, tag,
                                  reader.vartrees.dlls.servo_with_limits.dll_init)

    def test_IO_HS2(self):

        init = read_htc('main_hs2.htc')
        ivt = init.vartrees
        writer = HAWC2SInputWriter()
        for var in ['pgTorque', 'igTorque', 'Qg', 'pgPitch', 'igPitch', 'KK1',
                    'KK2', 'generatorFreq', 'generatorDamping', 'ffFreq',
                    'generatorSwitch', 'Kp2', 'Ko1', 'Ko2', 'ratedAeroPower',
                    'designTSR']:
            writer.vartrees = copy.deepcopy(ivt)
            val = np.ceil(np.random.rand()*10)

            setattr(writer.vartrees.dlls.risoe_controller.dll_init, var, val)

            writer.case_id = 'h2s'
            writer.execute()

            reader = read_htc('h2s.htc')

            assert_differet_types(self, val, var,
                                  reader.vartrees.dlls.risoe_controller.dll_init)

    def test_DLC(self):

        opt_tags = dlcdefs.excel_stabcon('./DLCs/', fext='xls', silent=True)

        reader = HAWC2InputReader()
        reader.htc_master_file = 'main_h2.htc'
        reader.execute()
        nan = np.nan
        std_ref = [[1.15468610453,20.623820042,0.0740770493071,0.0,1.63883455995e-13,0.0,2.92527680742e-13,0.0,1.05186187846e-13,0.00448084986512,0.0,0.0,0.0,0.0667195480442,0.0399570220043,0.0472457197536,8274.98592144,2262.75202353,1757.4399563,1312.15413729,160.871998807,1728.97639082,2342.04778773,2037.83218151,642.926932085,528.301580539,1863.59205836,28.5365363394,623.07691339,3799.86728189,97.9153652396,1184.56277166,4311.55945802,91.3815558524,337.08531561,716.108650974,9.8802695592,393.601651419,2052.93701803,35.1914829007,566.294446379,2164.10365741,31.2861603289,0.00519990875458,0.0317126196777,0.152721894284,0.112067971908,0.158007178759,0.181269098788,0.010676652669,0.433596301372,0.175455345668,0.013016978806,0.434643362013,0.359175518897,0.0214229824982,18.6265089573,2.06803508288,25.0941539507,0.212244149063,1.31535369178,1.02342877799,2.57083866195,0.131488392814,0.0976950118386,0.233573131863,0.00221016870252,0.00374349100002,0.00730561627894,0.0,0.0,0.0,0.0,0.0,0.0158365766493,0.000246378358914,0.000246373116992,0.0,17355.3928633,17298.5867446,0.0,0.0,0.0,0.000246376550094,1825.73496434,0.00514181322312,0.00511408917708,0.0,2.38418579102e-07,0.0,0.0,0.0,6.76987012183e-19,0.0,0.0,0.00607777508278,62.6460368858,19.4451050235,0.0,62.6460368858,0.0206862416792,642931.690171,1.78813934326e-07,0.0,0.0,0.0,0.0,0.0,12.84568112],
                   [1.15468610453,20.623820042,0.0740770493071,0.0,1.63883455995e-13,0.0,2.92527680742e-13,0.0,1.05186187846e-13,0.00448084986512,0.0,0.0,0.0,0.0618272310876,0.136897588111,0.0416334075112,8274.98592144,2262.75202353,1757.4399563,1312.15413729,160.871998807,1728.97639082,2342.04778773,2037.83218151,642.926932085,528.301580539,1863.59205836,28.5365363394,623.07691339,3799.86728189,97.9153652396,1184.56277166,4311.55945802,91.3815558524,337.08531561,716.108650974,9.8802695592,393.601651419,2052.93701803,35.1914829007,566.294446379,2164.10365741,31.2861603289,0.00519990875458,0.0317126196777,0.152721894284,0.112067971908,0.158007178759,0.181269098788,0.010676652669,0.433596301372,0.175455345668,0.013016978806,0.434643362013,0.359175518897,0.0214229824982,18.6265089573,2.06803508288,25.0941539507,0.0772392568926,0.922723143347,1.26687366066,2.22751007797,0.0900142147177,0.12071187503,0.202693610237,0.00358832154051,0.00394294579261,0.00862773061716,0.0,0.0,0.0,0.0,0.0,0.0237388360613,0.000246378358914,0.000246373116992,0.0,17355.3928633,17298.5867446,0.0,0.0,0.0,0.000246376550094,1825.73496434,0.00514181322312,0.00511408917708,0.0,2.38418579102e-07,0.0,0.0,0.0,6.76987012183e-19,0.0,0.0,0.00607777508278,62.6460368858,19.4451050235,0.0,62.6460368858,0.0206862416792,642931.690171,1.78813934326e-07,0.0,0.0,0.0,0.0,0.0,12.84568112],
                   [1.15468610453,20.623820042,0.0740770493071,0.0,1.63883455995e-13,0.0,2.92527680742e-13,0.0,1.05186187846e-13,0.00448084986512,0.0,0.0,0.0,0.0304850999051,0.0493848314637,0.018868019291,8274.98592144,2262.75202353,1757.4399563,1312.15413729,160.871998807,1728.97639082,2342.04778773,2037.83218151,642.926932085,528.301580539,1863.59205836,28.5365363394,623.07691339,3799.86728189,97.9153652396,1184.56277166,4311.55945802,91.3815558524,337.08531561,716.108650974,9.8802695592,393.601651419,2052.93701803,35.1914829007,566.294446379,2164.10365741,31.2861603289,0.00519990875458,0.0317126196777,0.152721894284,0.112067971908,0.158007178759,0.181269098788,0.010676652669,0.433596301372,0.175455345668,0.013016978806,0.434643362013,0.359175518897,0.0214229824982,18.6265089573,2.06803508288,25.0941539507,0.0290628277852,0.995506567128,1.02838432033,2.19556951898,0.0950211339582,0.093436891305,0.203007308904,0.00434149608136,0.00352270188921,0.00959668564148,0.0,0.0,0.0,0.0,0.0,0.0171455782373,0.000246378358914,0.000246373116992,0.0,17355.3928633,17298.5867446,0.0,0.0,0.0,0.000246376550094,1825.73496434,0.00514181322312,0.00511408917708,0.0,2.38418579102e-07,0.0,0.0,0.0,6.76987012183e-19,0.0,0.0,0.00607777508278,62.6460368858,19.4451050235,0.0,62.6460368858,0.0206862416792,642931.690171,1.78813934326e-07,0.0,0.0,0.0,0.0,0.0,12.84568112]]

        min_ref = [[1.02,18.2528,2.71844,0.0,-4.37326e-13,0.0,-8.74653e-13,0.0,-3.57812e-13,0.295914,0.0,0.0,0.0,-0.163887,0.10136,-0.0705254,-20387.8,-4770.09,-3219.08,-8557.29,-331.309,-3094.91,-8273.69,-10273.5,-1914.08,-1927.65,1908.95,14.4236,-2336.55,-15279.6,-307.368,-3830.56,-7509.66,-154.559,-1508.59,559.788,-7.70657,-505.606,-6453.2,-113.679,-1734.36,-3971.15,-52.0381,-0.0103338,-0.093089,-0.280412,-0.170996,0.0414103,-3.35387,86.3581,-1.38289,-3.61419,86.3357,-0.771033,-3.6373,86.3418,-88.9843,-13.8366,-203.745,-0.80437,-2.58734,-0.0269182,-4.64975,0.0696768,0.342503,-0.0503314,0.00624713,0.00382888,-0.00087769,0.0,0.0,0.0,0.0,0.0,0.00716154,0.311282,-0.316718,0.0,-2.24362e+07,2.22475e+07,0.0,5.13464e+06,0.0,-0.693718,-1.00068e+07,-5.94435,5.88228,0.0,1.43117,0.0,0.0,0.0,-1.96412e-18,0.0,0.0,0.13194,-397.661,-179.295,0.0,-397.661,-0.190739,-2.7327e+06,0.94,0.0,0.0,0.0,0.0,0.0,19.3039],
                   [1.02,18.2528,2.71844,0.0,-4.37326e-13,0.0,-8.74653e-13,0.0,-3.57812e-13,0.295914,0.0,0.0,0.0,-0.0786841,0.052737,-0.0416812,-20387.8,-4770.09,-3219.08,-8557.29,-331.309,-3094.91,-8273.69,-10273.5,-1914.08,-1927.65,1908.95,14.4236,-2336.55,-15279.6,-307.368,-3830.56,-7509.66,-154.559,-1508.59,559.788,-7.70657,-505.606,-6453.2,-113.679,-1734.36,-3971.15,-52.0381,-0.0103338,-0.093089,-0.280412,-0.170996,0.0414103,-3.35387,86.3581,-1.38289,-3.61419,86.3357,-0.771033,-3.6373,86.3418,-88.9843,-13.8366,-203.745,-0.0582635,0.618724,-0.90401,-2.7911,0.436152,0.290172,0.147358,0.00412495,0.00295085,-0.00239838,0.0,0.0,0.0,0.0,0.0,0.00323438,0.311282,-0.316718,0.0,-2.24362e+07,2.22475e+07,0.0,5.13464e+06,0.0,-0.693718,-1.00068e+07,-5.94435,5.88228,0.0,1.43117,0.0,0.0,0.0,-1.96412e-18,0.0,0.0,0.13194,-397.661,-179.295,0.0,-397.661,-0.190739,-2.7327e+06,0.94,0.0,0.0,0.0,0.0,0.0,19.3039],
                   [1.02,18.2528,2.71844,0.0,-4.37326e-13,0.0,-8.74653e-13,0.0,-3.57812e-13,0.295914,0.0,0.0,0.0,-0.020951,0.0971278,-0.0560079,-20387.8,-4770.09,-3219.08,-8557.29,-331.309,-3094.91,-8273.69,-10273.5,-1914.08,-1927.65,1908.95,14.4236,-2336.55,-15279.6,-307.368,-3830.56,-7509.66,-154.559,-1508.59,559.788,-7.70657,-505.606,-6453.2,-113.679,-1734.36,-3971.15,-52.0381,-0.0103338,-0.093089,-0.280412,-0.170996,0.0414103,-3.35387,86.3581,-1.38289,-3.61419,86.3357,-0.771033,-3.6373,86.3418,-88.9843,-13.8366,-203.745,-0.0530781,0.461831,-0.272344,-2.18937,0.428871,0.344658,0.211512,0.0052319,0.0040271,-0.00574167,0.0,0.0,0.0,0.0,0.0,0.00633216,0.311282,-0.316718,0.0,-2.24362e+07,2.22475e+07,0.0,5.13464e+06,0.0,-0.693718,-1.00068e+07,-5.94435,5.88228,0.0,1.43117,0.0,0.0,0.0,-1.96412e-18,0.0,0.0,0.13194,-397.661,-179.295,0.0,-397.661,-0.190739,-2.7327e+06,0.94,0.0,0.0,0.0,0.0,0.0,19.3039]]

        m12_ref = [[3.34676773271,59.7548542476,0.468059943101,nan,7.76083600193e-13,nan,1.49054185128e-12,nan,5.62240977778e-13,0.0248237883755,nan,nan,nan,0.214395394717,0.143972797391,0.164388856447,23175.7579878,8263.14369097,5629.17003676,4144.20611741,598.314219059,5506.35573909,7120.33780651,7408.99315106,3919.978416,2002.69311878,6421.48847582,116.400749013,2322.37864321,12700.9957402,283.781644151,4139.58768945,14634.4856958,287.905841024,1263.10491073,2599.44022781,42.4278111763,1310.31010198,5982.42125634,118.97185467,1783.03402345,7017.21799324,108.920729664,0.017825860875,0.0764152462997,0.515002924131,0.411861132403,0.606506071365,0.624728509324,0.0371483988949,1.20104703452,0.576028926024,0.0458456186893,1.27743259762,1.02915696571,0.0638579250184,51.3860015957,5.78773035326,70.463754011,0.772368825253,4.46989677658,3.09977639807,8.32214929559,0.419068094375,0.296971698047,0.722744001694,0.00882688351975,0.0127246329558,0.0263877197342,nan,nan,nan,nan,nan,0.0448986227494,0.0010587898629,0.00105881494617,nan,74317.6676684,74065.1727613,nan,nan,nan,0.00105881494617,9995.68097918,0.0223128089364,0.0222036427514,nan,nan,nan,nan,nan,1.87259744285e-18,nan,nan,0.0181154282518,810.407157775,260.90408433,nan,810.407157775,0.277556892867,3920433.62154,nan,nan,nan,nan,nan,nan,38.0249079495],
                  [3.34676773271,59.7548542476,0.468059943101,nan,7.76083600193e-13,nan,1.49054185128e-12,nan,5.62240977778e-13,0.0248237883755,nan,nan,nan,0.206887549082,0.386371713228,0.13328230517,23175.7579878,8263.14369097,5629.17003676,4144.20611741,598.314219059,5506.35573909,7120.33780651,7408.99315106,3919.978416,2002.69311878,6421.48847582,116.400749013,2322.37864321,12700.9957402,283.781644151,4139.58768945,14634.4856958,287.905841024,1263.10491073,2599.44022781,42.4278111763,1310.31010198,5982.42125634,118.97185467,1783.03402345,7017.21799324,108.920729664,0.017825860875,0.0764152462997,0.515002924131,0.411861132403,0.606506071365,0.624728509324,0.0371483988949,1.20104703452,0.576028926024,0.0458456186893,1.27743259762,1.02915696571,0.0638579250184,51.3860015957,5.78773035326,70.463754011,0.313304240049,3.29159645251,3.98032516206,7.2300707241,0.317039631799,0.361070772872,0.629128964369,0.0134177222875,0.0141452199232,0.0295225026225,nan,nan,nan,nan,nan,0.0659000592328,0.0010587898629,0.00105881494617,nan,74317.6676684,74065.1727613,nan,nan,nan,0.00105881494617,9995.68097918,0.0223128089364,0.0222036427514,nan,nan,nan,nan,nan,1.87259744285e-18,nan,nan,0.0181154282518,810.407157775,260.90408433,nan,810.407157775,0.277556892867,3920433.62154,nan,nan,nan,nan,nan,nan,38.0249079495],
                  [3.34676773271,59.7548542476,0.468059943101,nan,7.76083600193e-13,nan,1.49054185128e-12,nan,5.62240977778e-13,0.0248237883755,nan,nan,nan,0.10360011273,0.152019126971,0.0618277615999,23175.7579878,8263.14369097,5629.17003676,4144.20611741,598.314219059,5506.35573909,7120.33780651,7408.99315106,3919.978416,2002.69311878,6421.48847582,116.400749013,2322.37864321,12700.9957402,283.781644151,4139.58768945,14634.4856958,287.905841024,1263.10491073,2599.44022781,42.4278111763,1310.31010198,5982.42125634,118.97185467,1783.03402345,7017.21799324,108.920729664,0.017825860875,0.0764152462997,0.515002924131,0.411861132403,0.606506071365,0.624728509324,0.0371483988949,1.20104703452,0.576028926024,0.0458456186893,1.27743259762,1.02915696571,0.0638579250184,51.3860015957,5.78773035326,70.463754011,0.127498417605,3.86710462969,3.42210684433,7.35967487136,0.362446030724,0.292646731593,0.660302374728,0.0162588393672,0.0133652432195,0.0340428561448,nan,nan,nan,nan,nan,0.047025988233,0.0010587898629,0.00105881494617,nan,74317.6676684,74065.1727613,nan,nan,nan,0.00105881494617,9995.68097918,0.0223128089364,0.0222036427514,nan,nan,nan,nan,nan,1.87259744285e-18,nan,nan,0.0181154282518,810.407157775,260.90408433,nan,810.407157775,0.277556892867,3920433.62154,nan,nan,nan,nan,nan,nan,38.0249079495]]

        for icase, case in enumerate(opt_tags):
            case['[run_dir]'] = ''
            writer = HAWC2InputWriter()
            writer.vartrees = reader.vartrees
            writer.vartrees.aero.ae_filename = 'data/'+case['[Case id.]']+'_ae.dat'
            writer.vartrees.aero.pc_filename = 'data/'+case['[Case id.]']+'_pc.dat'
            writer.case_id = case['[case_id]']
            writer.vartrees.tags2var(case)
            writer.execute()

            wrapper = HAWC2Wrapper()
            wrapper.copyback_results = False
            wrapper.hawc2bin = 'HAWC2mb.exe'
            wrapper.log_directory = case['[log_dir]']
            wrapper.case_id = writer.case_id
            wrapper.verbose = False
            wrapper.compute()

            config = {}
            config['neq'] = 4
            config['no_bins'] = 46
            config['m'] = [12]
            output = HAWC2OutputBase(config)
            output.execute(case)

            self.compare_lists(output.stats['min'], min_ref[icase])
            self.compare_lists(output.stats['std'], std_ref[icase])
            self.compare_lists([eq[:1] for eq in output.eq], m12_ref[icase])

    def test_DLC_envelope(self):

        opt_tags = dlcdefs.excel_stabcon('./DLCs/', fext='xls', silent=True)

        reader = HAWC2InputReader()
        reader.htc_master_file = 'main_h2_envelope.htc'
        reader.execute()

        for icase, case in enumerate(opt_tags):
            case['[run_dir]'] = ''
            writer = HAWC2InputWriter()
            writer.vartrees = reader.vartrees
            writer.vartrees.aero.ae_filename = 'data/'+case['[case_id]']+'_ae.dat'
            writer.vartrees.aero.pc_filename = 'data/'+case['[case_id]']+'_pc.dat'
            writer.case_id = case['[case_id]']
            writer.vartrees.tags2var(case)
            writer.execute()

            wrapper = HAWC2Wrapper()
            wrapper.copyback_results = False
            wrapper.hawc2bin = 'HAWC2mb.exe'
            wrapper.log_directory = case['[log_dir]']
            wrapper.case_id = writer.case_id
            wrapper.verbose = False
            wrapper.compute()

            config = {}
            config['neq'] = 4
            config['no_bins'] = 46
            config['m'] = [12]

            ch_list = []
            string = 'blade%i-blade%i-node-%3.3i-'
            for iblade in range(1, 4):
                for i in range(1, 6):
                    ch_list.append([string % (iblade, iblade, i)+'momentvec-x',
                                    string % (iblade, iblade, i)+'momentvec-y',
                                    string % (iblade, iblade, i)+'momentvec-z',
                                    string % (iblade, iblade, i)+'forcevec-x',
                                    string % (iblade, iblade, i)+'forcevec-y',
                                    string % (iblade, iblade, i)+'forcevec-z'])

            config['ch_envelope'] = ch_list
            output = HAWC2OutputBase(config)
            output.execute(case)

    def test_hs2_workflow(self):

        reader = HAWC2InputReader()
        reader.htc_master_file = 'main_hs2.htc'
        reader.execute()
        writer = HAWC2SInputWriter()
        writer.vartrees = copy.copy(reader.vartrees)
        writer.case_id = 'h2s'
        writer.vartrees.dlls.risoe_controller.dll_init.nV = 21
        writer.execute()
        wrapper = HAWC2Wrapper()
        wrapper.hawc2bin = 'HAWC2s.exe'
        wrapper.case_id = writer.case_id
        wrapper.verbose = False
        wrapper.copyback_results = False
        wrapper.compute()
        output = HAWC2SOutputBase()
        output.case_id = wrapper.case_id
        output.commands = writer.vartrees.h2s.commands
        output.execute()

        freq = FreqDampTarget()

        freq.freqdamp = output.aeroelasticfreqdamp

        freq.mode_freq = np.array([[15., 0.25, 0.64],
                                   [20., 0.25, 0.64]])
        freq.mode_damp = np.array([[15., 0.8, 80.],
                                   [20., 0.8, 80.]])

        freq.mode_target_freq = np.array([[15., 0.3, -1],
                                          [20., 0.25, 0.6]])
        freq.mode_target_damp = np.array([[15., -1,  85.],
                                          [20., 0.4, 80.]])
        freq.execute()

        opt_ref = [[5.0, 0.1856139007801657E+01, 5.997000000000000],
                   [6.0, 0.7679248486373157E+00, 5.997000000000000],
                   [7.0, 0.2756222542925286E-02, 6.388794823610514],
                   [8.0, 0.2674377573694913E-02, 7.301487463662347],
                   [9.0, 0.2674368819285188E-02, 8.214173397307714],
                   [10.0, 0.2674368818350105E-02, 9.126859330341988],
                   [11.0, 0.2674368818350105E-02, 9.591314000000001],
                   [12.0, 0.4128043801574695E+01, 9.591314000000001],
                   [13.0, 0.6998870452044216E+01, 9.591314000000001],
                   [14.0, 0.9047112663285953E+01, 9.591314000000001],
                   [15.0, 0.1075513280195886E+02, 9.591314000000001],
                   [16.0, 0.1228078643556466E+02, 9.591314000000001],
                   [17.0, 0.1368494154839761E+02, 9.591314000000001],
                   [18.0, 0.1499512279716439E+02, 9.591314000000001],
                   [19.0, 0.1623682089352936E+02, 9.591314000000001],
                   [20.0, 0.1742409870354544E+02, 9.591314000000001],
                   [21.0, 0.1856907953949775E+02, 9.591314000000001],
                   [22.0, 0.1967436315021541E+02, 9.591314000000001],
                   [23.0, 0.2074408302584453E+02, 9.591314000000001],
                   [24.0, 0.2178409809683554E+02, 9.591314000000001],
                   [25.0, 0.2279885895128390E+02, 9.591314000000001]]

        for a, b in zip(opt_ref, output.operational_data):
            self.compare_lists(a, b)

if __name__ == '__main__':

    unittest.main()
