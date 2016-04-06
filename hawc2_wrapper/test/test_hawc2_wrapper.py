# -*- coding: utf-8 -*-
import unittest

from hawc2_wrapper.hawc2_inputreader import HAWC2InputReader
from hawc2_wrapper.hawc2_inputwriter import HAWC2SInputWriter, HAWC2InputWriter
from hawc2_wrapper.hawc2_wrapper import HAWC2Wrapper
from hawc2_wrapper.hawc2_output import HAWC2Output

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

        opt_tags = dlcdefs.excel_stabcon('./DLCs/', fext='xls')

        reader = HAWC2InputReader()
        reader.htc_master_file = 'main_h2.htc'
        reader.execute()
        nan = np.nan
        std_ref = [[25.980761472,102.87939865,1.19775424624,0.0,3.03156092694e-13,0.0,3.68954008942e-13,0.0,3.56182331568e-13,0.125455502805,963.241608551,628.797692956,160.55166562,1.3623849977,1.50991069111,0.769961429747,24301.5717971,6716.46410255,2258.73118208,1790.97062498,1116.03998335,2251.13797277,6180.04979484,6146.48601557,1024.05955382,3192.04443968,6017.77328002,78.0817064932,3256.18515444,6077.98340358,79.0705479654,3175.89455856,5991.7837718,77.6575559749,1910.61016525,2218.52061058,13.285367891,1943.38634302,2267.08410901,13.4352407662,1851.09598863,2260.94285937,13.7239389389,0.027165820368,0.0948445669269,0.0823223793887,0.182379294415,0.492759808358,1.1733210985,0.0265623075192,0.499066511412,1.19190315799,0.0262389055654,0.486917108345,1.1392949818,0.0270738751451,63.0888333137,5.41651634407,62.5097467956,1.89578405385,2.34541357627,2.49101667247,2.37351086356,0.265584045045,0.281613889787,0.267209505064,0.00682162964308,0.00709379540838,0.00693394070565,968999.065271,0.0,0.0,0.0,614339.048497,1.31459843824,0.124974267642,0.124974267642,0.0,8852592.91795,8120249.56273,181909.709002,0.0,0.0,0.124974267642,612756.391754,1.82722888768,1.81635413978,0.0,1.19209289551e-07,0.0,0.0,0.0,3.73018475478e-18,0.0,0.0,0.105710712194,1029842.62893,613735.275282,0.0,1029842.62893,652.90988182,1023171.61029,0.0,0.0,0.0,0.0,0.0,0.0,22.0531290297],
                   [25.980761472,103.895331945,1.3946664121,0.0,3.11003905197e-13,0.0,3.74381184814e-13,0.0,3.85232000459e-13,0.146059934858,2250.44118934,1519.68685429,258.912872518,0.812955901742,2.31816545118,0.602928891009,36490.1250167,6455.87837418,3767.03006908,2837.74697737,2452.37101126,3755.48960678,7609.8379703,7655.36386351,2243.1512564,5126.60019029,6352.02040187,81.4479956782,4952.79621386,6308.24264227,79.4887105891,5137.03158349,6291.02601771,79.5255598325,3022.09719794,2268.41718131,15.0274018232,2878.51228735,2331.53004298,15.4053958801,2972.93345472,2334.13107582,17.0208540392,0.02754486047,0.142761749146,0.0853524414762,0.282302068727,0.569309048733,1.77235192239,0.0328822900435,0.548094151981,1.7098182063,0.0301981292679,0.551003612029,1.76678421877,0.0322470835074,62.9678507035,5.67041787604,62.8154298101,0.805183403021,2.76391853348,2.70321407207,2.80000418526,0.306358377047,0.301456695601,0.311302545173,0.00933585968056,0.00850948676356,0.00892545017303,2115264.68704,0.0,0.0,0.0,1376482.14997,1.88710479214,0.145761323135,0.145761323135,0.0,10325348.2465,8519379.31329,1393134.68752,70826.0465065,0.0,0.145761317319,1374907.03972,2.04939481243,2.03727003218,0.0,1.19209289551e-07,0.0,0.0,0.0,3.8688377115e-18,0.0,0.0,0.137767763765,2249448.88274,1375956.07077,0.0,2249448.88274,1463.78297845,2242368.50684,0.0,0.0,0.0,0.0,0.0,0.0,21.9873195779],
                   [25.980761472,102.921778281,1.4114312728,0.0,3.24067351781e-13,0.0,3.81976361011e-13,0.0,3.95277676401e-13,0.147827161474,2469.80193106,1770.66461044,279.532460449,0.763694739374,1.89240601374,0.530599785914,34839.3868642,7997.42934698,2360.08681422,2127.14580015,2340.63068576,2353.89074721,6283.83217712,6447.27923184,2378.85684143,5128.43749217,6073.57264951,76.170195396,5037.6697936,6149.21082734,74.5352265711,5144.39291536,6012.89088454,73.9953451989,2993.22054249,2217.86200793,24.8235148375,2949.2448739,2278.80475805,22.4700294219,2965.13626204,2273.58972865,25.1429623465,0.0347200145443,0.142072471882,0.0905932842312,0.117899504402,0.530715008398,1.81381950294,0.0303311344649,0.533861178816,1.79753007383,0.0302639293985,0.51667144445,1.81438455666,0.0312613555621,63.0543645507,5.71964177928,62.6491925992,0.989449662025,2.36115197064,2.34549083497,2.3633816697,0.272629734464,0.270246685348,0.271565329407,0.00562004053215,0.00567565433354,0.005871163272,2238640.75684,0.0,0.0,0.0,1510407.42092,1.94154573121,0.147313101368,0.147313095613,0.0,10435069.8997,8531062.4412,1754359.06011,248398.181881,0.0,0.147313107122,1507926.72816,2.00722274358,1.99542184345,0.0,1.19209289551e-07,0.0,0.0,0.0,3.3365445746e-18,0.0,0.0,0.0559485659318,2380416.41888,1509660.48001,0.0,2380416.41888,1606.02180586,2377812.90473,0.0,0.0,0.0,0.0,0.0,0.0,22.1821625035]]

        min_ref = [[10.02,0.007374,2.93084,0.0,-1.23246e-12,0.0,-1.56046e-12,0.0,-1.31198e-12,0.307094,-111.43,-34.5224,35.2633,-5.63649,1.10979,-3.44384,-12636.3,-22547.5,-8584.39,-13673.7,-604.216,-8488.27,-10531.7,-14209.3,-3533.36,-15645.2,-9326.74,-161.635,-16093.0,-10031.9,-179.831,-14888.4,-10102.1,-196.386,-8716.64,-4312.66,-57.3807,-9209.12,-4234.28,-52.3428,-8174.98,-4531.72,-56.5243,-0.103772,-0.065552,-0.283331,-0.666991,-0.783977,-2.94692,86.3862,-0.878581,-2.9593,86.3822,-0.945276,-2.94507,86.3845,-89.2771,-21.267,-208.229,-9.07854,-1.39373,-3.00256,-1.4435,0.213931,0.0314179,0.212227,-0.0105143,-0.00891211,-0.010392,0.0,0.0,0.0,0.0,0.0,0.313935,0.307369,-0.320631,0.0,-2.27135e+07,1.51122e+06,0.0,5.13464e+06,0.0,-0.697631,-1e+07,-6.03438,1.27378,0.0,1.43117,0.0,0.0,0.0,-9.07979e-18,0.0,0.0,0.0371245,-3.49597e+06,-6.35176e-32,0.0,-3.49597e+06,-6.75719e-35,-55701.0,0.94,0.0,0.0,0.0,0.0,0.0,13.6325],
                   [10.02,0.0136868,2.93111,0.0,-1.67973e-12,0.0,-2.07606e-12,0.0,-1.75924e-12,0.306902,-107.258,-33.2245,34.4584,-1.41093,0.911508,-1.72546,-12604.0,-18354.0,-16162.6,-17526.6,-455.889,-16071.4,-18635.7,-20801.9,-6193.43,-23198.1,-9214.01,-163.594,-21822.8,-10028.0,-178.707,-21660.5,-10124.6,-202.34,-12942.8,-4249.18,-74.1218,-12712.0,-4230.05,-63.088,-12394.6,-4670.25,-60.6649,-0.105407,-0.0652236,-0.246421,-0.98226,-0.78978,-3.0446,86.3097,-0.867526,-3.04786,86.3378,-0.957038,-3.02007,86.3311,-89.3235,-21.148,-208.224,-2.68495,-0.562902,-1.54243,-1.08345,0.295658,0.299077,0.269726,-0.0283761,-0.0134863,-0.0131455,0.0,0.0,0.0,0.0,0.0,0.292177,0.307434,-0.320566,0.0,-2.27089e+07,789065.0,0.0,5.13464e+06,0.0,-0.697566,-1e+07,-6.03284,0.859205,0.0,1.43117,0.0,0.0,0.0,-1.06389e-17,0.0,0.0,0.0442893,-6.08975e+06,-6.3545e-32,0.0,-6.08975e+06,-6.76011e-35,-56033.5,0.94,0.0,0.0,0.0,0.0,0.0,14.3767],
                   [10.02,0.0614566,2.94181,0.0,-1.98785e-12,0.0,-1.86858e-12,0.0,-2.17172e-12,0.308258,-94.3728,-29.2383,40.802,-2.71504,0.546443,-1.99343,-11707.2,-19681.9,-5656.57,-10488.1,-147.878,-5585.83,-14092.1,-10920.1,-7250.94,-25005.7,-9213.52,-164.467,-21217.3,-10078.7,-175.765,-24073.4,-10133.9,-200.089,-13247.2,-4532.57,-70.7574,-12460.8,-4280.26,-73.1965,-13140.8,-4519.88,-61.8976,-0.140256,-0.0618888,-0.316896,-0.54958,-0.787699,-2.98123,86.2923,-0.868469,-2.94488,86.3158,-0.954463,-2.97615,86.271,-89.3288,-21.1504,-208.226,-3.71445,0.210976,0.479871,0.131402,0.391109,0.420219,0.39089,-0.00782037,-0.00739773,-0.00767305,0.0,0.0,0.0,0.0,0.0,0.186196,0.308662,-0.319338,0.0,-2.26219e+07,-2.19485e+06,0.0,5.13464e+06,0.0,-0.696338,-1e+07,-6.00476,0.584721,0.0,1.43117,0.0,0.0,0.0,-8.91359e-18,0.0,0.0,0.037966,-7.10505e+06,-6.35488e-32,0.0,-7.10505e+06,-6.76051e-35,-56745.5,0.94,0.0,0.0,0.0,0.0,0.0,12.9002]]

        max_ref = [[100.0,359.974,6.08952,0.0,1.19271e-12,0.0,1.55052e-12,0.0,1.74931e-12,0.639168,3879.29,2456.82,626.611,3.23988,7.86517,1.57835,95624.1,24279.0,5092.76,499.942,4229.33,5002.48,16926.1,15062.8,55.701,-1036.12,9391.36,120.201,-830.451,11785.8,179.588,-818.857,10818.4,162.186,-278.16,3268.48,17.2821,-322.218,4844.78,30.8788,-198.589,4303.98,27.6126,0.0789184,0.357567,0.277738,0.726607,1.0242,2.15298,86.4755,1.05813,2.40106,86.4753,0.996645,1.68774,86.4753,89.1424,-2.63827,-30.2652,8.08816,10.9799,10.7082,10.5377,1.51757,1.52892,1.49865,0.0434512,0.0472945,0.0416291,3.2862e+06,0.0,0.0,0.0,2.08441e+06,6.12466,0.636206,0.00820559,0.0,581258.0,2.26115e+07,733813.0,5.13464e+06,0.0,-0.368794,-7.91679e+06,-1.28185,5.99817,0.0,1.43117,0.0,0.0,0.0,5.72093e-18,0.0,0.0,0.370709,2.18295e-31,2.08383e+06,0.0,2.18295e-31,2216.84,3.53336e+06,0.94,0.0,0.0,0.0,0.0,0.0,87.7595],
                   [100.0,359.974,6.68242,0.0,1.86858e-12,0.0,1.97791e-12,0.0,1.68967e-12,0.697762,9066.0,6134.74,926.235,2.8196,11.0493,1.83653,142665.0,22543.1,7944.08,2466.78,6940.04,7894.55,21376.1,21071.6,56.0335,-876.063,12601.9,192.94,-940.674,11874.0,179.065,-639.951,11821.7,183.88,-196.483,3300.53,35.5697,-90.9386,4843.37,61.3506,26.9575,4269.65,74.9323,0.0532082,0.519867,0.290918,0.759182,1.54488,4.06629,86.4758,1.37078,3.73832,86.4755,1.3948,3.83104,86.4756,89.1645,-1.14921,-30.3694,3.49083,14.5857,12.9525,12.5531,1.88972,1.74555,1.70931,0.0717285,0.046438,0.0551287,5.72437e+06,0.0,0.0,0.0,3.98867e+06,7.23954,0.697983,0.0699828,0.0,4.95753e+06,2.2607e+07,5.72437e+06,5.72437e+06,0.0,-0.307017,-6.01394e+06,-0.864577,5.99664,0.0,1.43117,0.0,0.0,0.0,6.86847e-18,0.0,0.0,0.435756,2.18295e-31,3.98765e+06,0.0,2.18295e-31,4242.19,6.18441e+06,0.94,0.0,0.0,0.0,0.0,0.0,87.9016],
                   [100.0,359.92,7.2904,0.0,2.00773e-12,0.0,1.67973e-12,0.0,2.13197e-12,0.761335,10043.2,7585.74,1130.18,2.54286,7.57967,1.01652,149031.0,29158.2,10742.8,4481.79,7830.23,10716.5,14218.3,13712.1,56.7455,-1184.44,10589.0,190.231,-1242.57,11801.8,178.324,-783.275,10878.7,176.899,-198.903,3218.34,136.677,-331.354,4823.57,91.1343,-28.3344,4244.52,113.661,0.058821,0.575637,0.30969,0.415156,1.34501,4.48025,86.4755,1.36173,3.98288,86.4757,1.26762,4.67083,86.4757,89.2161,-0.900187,-30.2011,4.37199,11.9315,10.5781,10.8485,1.64101,1.53694,1.55442,0.0415635,0.0409956,0.0406153,6.68346e+06,0.0,0.0,0.0,5.09818e+06,6.3419,0.754192,0.126192,0.0,8.93097e+06,2.25204e+07,6.68346e+06,6.68346e+06,0.0,-0.250808,-4.91971e+06,-0.58733,5.96872,0.0,1.43117,0.0,0.0,0.0,5.80763e-18,0.0,0.0,0.229449,2.18295e-31,5.09458e+06,0.0,2.18295e-31,5419.77,7.25094e+06,0.94,0.0,0.0,0.0,0.0,0.0,87.8494]]

        m3_ref = [[58.372274137,285.633464484,2.0491146106,nan,1.80267811094e-12,nan,2.19078936538e-12,nan,2.14860115464e-12,0.215424711215,2588.90590519,1616.21215778,383.622044865,5.81762063534,4.38256217222,3.32708754168,75392.575835,33554.178112,9293.94705981,9849.24512705,3135.64212276,9172.18175503,18688.7633343,19653.1098558,2328.31361268,9493.50161612,14139.3832288,193.291735198,9902.91170964,15050.7464796,236.061717822,9147.68635324,15032.1949976,238.66728356,5493.78354972,5453.28331684,51.6403677393,5771.83441443,6081.46003438,54.5173192177,5191.34336446,6110.26408592,57.136155084,0.129928178856,0.287604581374,0.401262051848,0.974227180895,1.2419482972,3.31484151251,0.0579368639397,1.33172596343,3.48139628399,0.0604103451072,1.31448781487,3.01163013346,0.0589219053584,142.375563138,12.7397363351,142.172097826,11.444058934,8.6725402527,9.06133443981,8.65881336039,0.917381205689,0.982575480158,0.893882851027,0.0372919758987,0.0381700435128,0.0378665699118,2131840.04522,nan,nan,nan,1352208.84567,3.76956248986,0.213324763236,0.213324510087,nan,15111891.5307,13688278.8234,487111.644259,nan,nan,0.213324782569,1351430.37569,3.08308493847,3.06482977325,nan,nan,nan,nan,nan,1.02021416759e-17,nan,nan,0.216404907968,2267923.08529,1351832.58518,nan,2267923.08529,1438.11955244,2328313.54286,nan,nan,nan,nan,nan,nan,62.3132286236],
                  [58.372274137,289.164660982,2.43356843968,nan,2.38507489214e-12,nan,2.76604434158e-12,nan,2.52550744371e-12,0.253560661277,5951.05069319,4001.35781766,578.517763534,2.79963122198,6.57927963846,2.36332606925,103377.417374,28189.6772879,15856.2004987,14366.0022296,4797.92400192,15748.4366872,27550.9408849,27753.286004,4054.18322884,14532.6072198,15473.2327859,233.79525088,13593.5646161,15975.1798398,244.242180221,13681.2551189,16025.1462025,256.310158608,8299.8040736,5644.22677297,73.2788529221,8205.73413501,6150.78803989,81.724406254,8077.07190658,6290.88812397,88.0225740646,0.109649215165,0.39255277588,0.382129494977,1.23910521806,1.5488682569,4.62107882923,0.115662054237,1.50959756553,4.40730507074,0.0968597108665,1.55660039924,4.45044367946,0.103844625199,143.383645243,13.8518476495,143.759447149,4.27185422031,10.2811728748,9.72203670644,9.32134290515,1.04754860308,0.973229061302,0.963352343111,0.066438498524,0.0443145675256,0.0485288133116,3713541.84154,nan,nan,nan,2587549.88533,4.50692802725,0.253358916088,0.25335876142,nan,17947904.3827,14153839.5524,3713543.31877,382572.585317,nan,0.253358877421,2585856.71312,3.35278119348,3.33278257089,nan,nan,nan,nan,nan,1.16276429575e-17,nan,nan,0.25396055961,3950572.9765,2586888.18584,nan,3950572.9765,2752.01460548,4048331.61501,nan,nan,nan,nan,nan,nan,62.7069227979],
                  [58.372274137,289.253170275,2.82103906983,nan,2.75865538922e-12,nan,2.55057097321e-12,nan,2.97327576223e-12,0.293922381808,6576.49693,4940.02684525,706.706762599,3.45485352196,4.5626987689,1.9595643513,104274.97232,34562.7771144,10854.5985056,10077.4852287,5175.59799846,10789.612343,19234.8466768,18041.1419654,4740.67812641,15453.4899557,14517.4681092,230.888836913,12958.1415934,15274.2856173,234.838220055,15108.9304904,15241.3993174,250.483188455,8464.80004205,5556.68490942,134.758970194,7868.68888229,6064.92855283,108.028808378,8506.39580188,6071.16122712,116.380934408,0.137514265987,0.413578992045,0.446946072342,0.652449763789,1.39795326607,4.84045623343,0.118882826493,1.47236787384,4.49421439499,0.104197635366,1.46022246216,4.96078936269,0.133775171846,144.317624702,13.55954418,143.908949966,5.38808690558,7.61546589869,6.60284372376,6.99498221165,0.811192604394,0.726402983366,0.756088591224,0.0350661620214,0.0325383090314,0.03406107431,4335727.48726,nan,nan,nan,3307316.74327,3.9933589014,0.289026440781,0.289026440781,nan,20469135.1128,16033400.1819,4335727.48747,1004758.23104,nan,0.289026440781,3295711.05329,3.51442206806,3.4927346246,nan,nan,nan,nan,nan,1.0043841911e-17,nan,nan,0.124221356067,4609223.45362,3304981.33333,nan,4609223.45362,3515.94022382,4740678.16529,nan,nan,nan,nan,nan,nan,63.6229509798]]

        for icase, case in enumerate(opt_tags):
            case['[run_dir]'] = ''
            writer = HAWC2InputWriter()
            writer.vartrees = reader.vartrees
            writer.vartrees.aero.ae_filename = 'data/'+case['[Case id.]']+'_ae.dat'
            writer.vartrees.aero.pc_filename = 'data/'+case['[Case id.]']+'_pc.dat'
            writer.case_id = case['[Case id.]']
            writer.vartrees.tags2var(case)
            writer.execute()

            wrapper = HAWC2Wrapper()
            wrapper.copyback_results = False
            wrapper.hawc2bin = 'HAWC2mb.exe'
            wrapper.log_directory = case['[log_dir]']
            wrapper.case_id = writer.case_id
            #wrapper.compute()

            config = {}
            config['neq'] = 90
            config['no_bins'] = 46
            config['m'] = [12]
            output = HAWC2Output(config)
            output.execute(case)

            self.compare_lists(output.stats['max'], max_ref[icase])
            self.compare_lists(output.stats['min'], min_ref[icase])
            self.compare_lists(output.stats['std'], std_ref[icase])
            self.compare_lists([eq[:1] for eq in output.eq], m3_ref[icase])

    def test_DLC_envelope(self):

        opt_tags = dlcdefs.excel_stabcon('./DLCs/', fext='xls')

        reader = HAWC2InputReader()
        reader.htc_master_file = 'main_h2_envelope.htc'
        reader.execute()

        for icase, case in enumerate(opt_tags):
            case['[run_dir]'] = ''
            writer = HAWC2InputWriter()
            writer.vartrees = reader.vartrees
            writer.vartrees.aero.ae_filename = 'data/'+case['[Case id.]']+'_ae.dat'
            writer.vartrees.aero.pc_filename = 'data/'+case['[Case id.]']+'_pc.dat'
            writer.case_id = case['[Case id.]']
            writer.vartrees.tags2var(case)
            writer.execute()

            wrapper = HAWC2Wrapper()
            wrapper.copyback_results = False
            wrapper.hawc2bin = 'HAWC2mb.exe'
            wrapper.log_directory = case['[log_dir]']
            wrapper.case_id = writer.case_id
            #wrapper.compute()

            config = {}
            config['neq'] = 90
            config['no_bins'] = 46
            config['m'] = [12]
            
            ch_list = []
            for iblade in range(1, 4):
                for i in range(1, 6):
                    ch_list.append(['blade%i-blade%i-node-%3.3i-momentvec-x'%(iblade, iblade, i),
                                    'blade%i-blade%i-node-%3.3i-momentvec-y'%(iblade, iblade, i),
                                    'blade%i-blade%i-node-%3.3i-momentvec-z'%(iblade, iblade, i), 
                                    'blade%i-blade%i-node-%3.3i-forcevec-x'%(iblade, iblade, i),
                                    'blade%i-blade%i-node-%3.3i-forcevec-y'%(iblade, iblade, i),
                                    'blade%i-blade%i-node-%3.3i-forcevec-z'%(iblade, iblade, i)])
                                
            config['ch_envelope'] = ch_list
            output = HAWC2Output(config)
            output.execute(case)
            print output.envelope
            
            


if __name__ == '__main__':

    unittest.main()
