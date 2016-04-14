"""
Module containing classes to read HAWC2 and HAWCStab2 results and perform basic
postprocessing.
"""
import numpy as np
import os
import glob
import re
from wetb.prepost import Simulations

from hawc2_vartrees import DTUBasicControllerVT


class HAWC2OutputBase(object):
    """
    Class that reads an HAWC2 output file and computes statistics and damage
    equivalent load of all the channels in the file.

    Parameters
    ----------

    case_tags: dict
        Dictionary with tags as keys. Keys required [run_dir], [res_dir], and
        [case_id].
    config: dict
        * m: list. Values of the slope of the SN curve.
        * no_bins: int. Number of bins for the binning of the amplitudes.
        * neq: int. Number of equivalent cycles

    Returns
    -------
    stats: dict
        Dictionary containing lists of the statistics of each channel. The keys
        are 'std', 'rms', 'min', 'int', 'max', 'range', 'absmax', and 'mean'.
    eq: list
        List containing a list of each channel containing the damage equivalent
        load for each value of m.

    Example
    -------
    >>> case = {}
    >>> case['[case_id]'] = 'dlc12_wsp04_wdir000_s1001'
    >>> case['[res_dir]'] = 'res/dlc12_iec61400-1ed3'
    >>> config = {}
    >>> config['neq'] = 600
    >>> config['no_bins'] = 2**7
    >>> config['m'] = [12]
    >>> output = HAWC2OutputBase(config)
    >>> output.execute(case)

    """
    def __init__(self, config):

        self.eq = []
        self.stats = []

        self.no_bins = config['no_bins']
        self.neq = config['neq']
        self.m = config['m']
        if 'ch_envelope' not in config.keys():
            self.ch_envelope = []
        else:
            self.ch_envelope = config['ch_envelope']

    def execute(self, case_tags):

        case_tags['[run_dir]'] = ''
        case = Simulations.Cases(case_tags)
        case.load_result_file(case_tags)
        self.ch_dict = case.res.ch_dict
        self.stats = case.res.calc_stats(case.sig, i0=0, i1=None)
        self.eq = []
        for s in case.sig.T:
            self.eq.append(case.res.calc_fatigue(s, no_bins=self.no_bins,
                                                 neq=self.neq,
                                                 m=self.m))
        if self.ch_envelope != []:
            self.envelope = case.compute_envelope(case.sig, self.ch_envelope)
        else:
            self.envelope = {}


class HAWC2OutputCompact(HAWC2OutputBase):
    """
    HAWC2SOutput: HAWC2Output class that organize results in
    compact arrays to minimize the number of outputs. The outputs are selected
    with a dictionary passed in the initialization.

    Parameters
    ----------

    Returns
    -------
    outputs_statistics: array
        Statistics outputs.

    outputs_fatigue: array
        Fatigue outputs.
    """
    def __init__(self, config):
        super(HAWC2OutputCompact, self).__init__(config)

        self.channel = config['channels']

        if 'stat_list' in config.keys():
            self.stat_list = config['stat_list']
        else:
            self.stat_list = ['std', 'rms', 'min', 'int', 'max', 'range',
                              'absmax', 'mean']
        self.Nch = len(self.channel)
        self.Nstat = len(self.stat_list)

    def execute(self, case_tags):
        super(HAWC2OutputCompact, self).execute(case_tags)

        self.outputs_statistics = np.zeros([self.Nstat, self.Nch])
        self.outputs_fatigue = np.zeros([len(self.m), self.Nch])

        for ich, ch in enumerate(self.channel):
            self.outputs_statistics[:, ich] = \
                np.array([self.stats[s][self.ch_dict[ch]['chi']] for s in self.stat_list])
            self.outputs_fatigue[:, ich] = \
                np.array(self.eq[self.ch_dict[ch]['chi']])

        for ch in self.ch_envelope:
            setattr(self, 'env_'+ch[0].replace('-', '_'), np.array(self.envelope[ch[0]]))


class HAWC2SOutputBase(object):
    """
    HAWC2SOutputBase: class that reads HAWC2s output files.

    Parameters
    ----------
    case_id: str
        Name of solution to open.
    commands: list
        List containing the strings of the HAWC2s commands that have been
        executed. Only files associated with these commands are read.

    Returns
    -------
    operational_data: list
        List containing results included in the .opt file.

    rotor_loads_data: list
        List containing the results included in the .pwr file.

    blade_loads_data: list
        List containing the results included in all the .ind files.

    structuralfreqdamp: list
        List containing the structural frequencies and damping ratios.

    aeroelasticfreqdamp: list
        List containing the aeroelastic frequencies and damping ratios.

    aeroservoelasticfreqdamp: list
        List containing the aeroservoelastic frequencies and damping ratios.

    controller_data: DTUBasicControllerVT
        Variable tree containing the controller tuning inputs.

    Exalmple
    --------

    >>> output = HAWC2SOutputBase()
    >>> output.case_id = wrapper.case_id
    >>> output.commands = writer.vartrees.h2s.commands
    >>> output.execute()

    """
    def __init__(self):
        self.case_id = ''
        self.commands = []
        self.structuralfreqdamp = np.array([0.])
        self.aeroelasticfreqdamp = np.array([0.])
        self.aeroservoelasticfreqdamp = np.array([0.])
        self.controller_data = DTUBasicControllerVT()
        self.cl_matrices = []

    def execute(self):

        for name in self.commands:
            if name == 'compute_optimal_pitch_angle':

                data = np.loadtxt(self.case_id + '.opt', skiprows=1)
                if len(data.shape) == 1:
                    data = data.reshape(1, data.shape[0])
                self.operational_data = data[:, :3]

            elif name == 'compute_steady_states':
                # To read the opt file even when they are not computed
                # We do it in 'compute_steadystate' because this command
                # requires the opt file, so it should be there!
                if 'compute_optimal_pitch_angle' not in self.commands:
                    data = np.loadtxt(self.case_id + '.opt', skiprows=1)
                    # read first line
                    if len(data.shape) == 1:
                        data = data.reshape(1, data.shape[0])
                    self.operational_data = data[:, :3]

            elif name == 'compute_stability_analysis':

                data = np.loadtxt(self.case_id + '.cmb', skiprows=1)
                if len(data.shape) == 1:
                    data = data.reshape(1, data.shape[0])
                self.aeroelasticfreqdamp = data

            elif name == 'compute_aeroservoelastic':

                data = np.loadtxt(self.case_id + '_Servo.cmb', skiprows=1)

                if len(data.shape) == 1:
                    data = data.reshape(1, data.shape[0])
                self.aeroservoelasticfreqdamp = data

            elif name == 'compute_controller_input':
                fid = open('controller_input.txt', 'r')
                line = fid.readline()
                line = fid.readline()
                temp = line[line.find('K =')+4:line.rfind('[')]
                self.controller_data.Qg = float(temp.strip())
                line = fid.readline()
                line = fid.readline()
                line = fid.readline()
                temp = line[line.find('Kp =')+5:line.rfind('[')]
                self.controller_data.pgTorque = float(temp.strip())
                line = fid.readline()
                temp = line[line.find('Ki =')+5:line.rfind('[')]
                self.controller_data.igTorque = float(temp.strip())
                line = fid.readline()
                line = fid.readline()
                temp = line[line.find('Kp =')+5:line.rfind('[')]
                self.controller_data.pgPitch = float(temp.strip())
                line = fid.readline()
                temp = line[line.find('Ki =')+5:line.rfind('[')]
                self.controller_data.igPitch = float(temp.strip())
                line = fid.readline()
                temp = line[line.find('K1 =')+5:line.rfind('[deg]')]
                self.controller_data.KK1 = float(temp.strip())
                temp = line[line.find('K2 =')+5:line.rfind('[deg^2]')]
                self.controller_data.KK2 = float(temp.strip())
                fid.close()

            elif name == 'save_beam_data':
                print 'not implemented yet'

            elif name == 'save_blade_geometry':
                print 'not implemented yet'

            elif name == 'save_aero_point_data':
                print 'not implemented yet'

            elif name == 'save_profile_coeffs':
                print 'not implemented yet'

            elif name == 'compute_structural_modal_analysis':
                data = np.loadtxt(self.case_id + '_struc.cmb', skiprows=1)
                if len(data.shape) == 1:
                    data = data.reshape(1, data.shape[0])
                self.structuralfreqdamp = data
            elif name == 'save_power':
                # read pwr file
                data = np.loadtxt(self.case_id+'.pwr', skiprows=1)

                if len(data.shape) == 1:
                    data = data.reshape(1, data.shape[0])
                self.rotor_loads_data = data.copy()

            elif name == 'save_induction':
                self.blade_loads_data = []
                wsp_array = []
                wsp_files = glob.glob(self.case_id+'_u*.ind')
                if len(wsp_files) > 0:
                    for f in wsp_files:
                        w = float(re.sub('%s' % self.case_id + '_u',
                                         '', f).strip('.ind'))/1000.
                        wsp_array.append(w)
                    wsp_array.sort()

                    for wsp in wsp_array:
                        filename = self.case_id + '_u' + \
                                   str(int(wsp*1000)) + '.ind'
                        data = np.loadtxt(filename)
                        self.blade_loads_data.append(data)
                else:
                    self.blade_loads_data.append(np.zeros((30, 34)))

                # read fext files
                self.blade_fext_loads_data = []
                wsp_array = []
                wsp_files = glob.glob(self.case_id+'_fext_u*.ind')
                if len(wsp_files) > 0:
                    for f in wsp_files:
                        w = float(re.sub('%s' % self.case_id + '_fext_u',
                                         '', f).strip('.ind'))/1000.
                        wsp_array.append(w)
                    wsp_array.sort()

                    for wsp in wsp_array:
                        filename = self.case_id + '_fext_u' + \
                                   str(int(wsp*1000)) + '.ind'
                        data = np.loadtxt(filename)
                        self.blade_fext_loads_data.append(data)
                else:
                    # fixme
                    self.blade_fext_loads_data.append(np.zeros((30, 34)))

            elif name == 'save_cl_matrices_all':
                wsp_files = glob.glob(self.case_id+'_amat_ase_ops_*.dat')
                nwsp = len(wsp_files)

                for iws in range(1,nwsp+1):
                    name_file = self.case_id+'_amat_ase_ops_%i.dat' % iws
                    Aase = np.matrix(np.loadtxt(name_file))

                    #name_file = self.case_id+'_bvmat_ase_ops_%i.dat' % iws
                    #Bvase = np.loadtxt(name_file)

                    name_file = self.case_id+'_bvmat_loc_v_ase_ops_%i.dat' % iws
                    Bvlocase = np.matrix(np.loadtxt(name_file))

                    name_file = self.case_id+'_emat_ase_ops_%i.dat' % iws
                    Case = np.matrix(np.loadtxt(name_file))

                    #name_file = self.case_id+'_fmat_ase_ops_%i.dat' % iws
                    Dase = np.matrix(np.zeros([Case.shape[0], Bvlocase.shape[1]]))

                    self.cl_matrices.append([Aase, Bvlocase, Case, Dase])

class HAWC2SOutput(HAWC2SOutputBase):
    """
    HAWC2SOutput: HAWC2SOutputBase class that organize results in more general
    arrays.

    Parameters
    ----------

    Returns
    -------
    wsp : array [nW]
        Wind speed [m/s].
    pitch : array [nW]
        Pitch angle [deg].
    rpm : array [nW]
        Rotational speed [rpm].
    P : array [nW]
        Aerodynamic power [W].
    Q : array [nW]
        Aerodynamic torque [Nm].
    T : array [nW]
        Thrust [N].
    CP : array [nW]
        Power coefficient [-].
    CT : array [nW]
        Thrust coefficient [-].
    Mx : array [nW]
        Blade root in-plane bending moment [Nm].
    My : array [nW]
        Blade root out-of-plane bending moment [Nm].
    Mz : array [nW]
        Blade root torsional moment [Nm].
    tip_pos : array [nW]
        Blade tip position [m].
    tip_rot : array [nW]
        Blade tip rotation [deg].
    disp_x : array [nW]
        In plane blade deflection [m].
    disp_y : array [nW]
        Out of plane blade deflection [m].
    disp_z : array [nW]
        Blade deflection along the blade axis [m].
    disp_rot_z : array [nW]
        Blade sections rotation [deg].
    s: array [nS]
        Position of radial sections [-].
    aoa : array [nW, nS]
        Sections angle of attack  [deg].
    Ft : array [nW, nS]
        Sections tangential force [N].
    Fn : array [nW, nS]
        Sections normal force  [N].
    cl : array [nW, nS]
        Lift coefficients [-].
    cd : array [nW, nS]
        Drag coefficient [-].
    cm : array [nW, nS]
        Moment coefficient [-].
    ct : array [nW, nS]
        Thrust coefficient [-].
    cp : array [nW, nS]
        Power coefficient [-].
    v_a : array [nW, nS]
        Axial induced velocity [m/s].
    v_t : array [nW, nS]
        Tangential induced velocity [m/s].
    Fx : array [nW, nS]
        Integrated lateral force [N].
    Fy : array [nW, nS]
        Intagrated longitudinal force [N].
    Fx_e: array [nW, nS]
        Intagrated sectional in plane force [N]
        in the element coordinate system.
    Fy_e: array [nW, nS]
        Intagrated sectional out of plane force [N]
        in the element coordinate system.
    Fz_e: array [nW, nS]
        Intagrated sectional radial force [N]
        in the element coordinate system.
    Mx_e: array [nW, nS]
        Intagrated sectional in plane moment [Nm]
        in the element coordinate system.
    My_e: array [nW, nS]
        Intagrated sectional out of plane moment [Nm]
        in the element coordinate system.
    Mz_e: array [nW, nS]
        Intagrated sectional radial moment [Nm]
        in the element coordinate system.
    Fx_r: array [nW, nS]
        Intagrated sectional in plane force [N]
        in the rotor coordinate system.
    Fy_r: array [nW, nS]
        Intagrated sectional out of plane force [N]
        in the rotor coordinate system.
    Fz_r: array [nW, nS]
        Intagrated sectional radial force [N]
        in the rotor coordinate system.
    Mx_r: array [nW, nS]
        Intagrated sectional in plane moment [Nm]
        in the rotor coordinate system.
    My_r: array [nW, nS]
        Intagrated sectional out of plane moment [Nm]
        in the rotor coordinate system.
    Mz_r: array [nW, nS]
        Intagrated sectional radial moment [Nm]
        in the rotor coordinate system.
    outlist : list
        List with the names of all the outputs
    """
    def __init__(self):
        super(HAWC2SOutput, self).__init__()
        self.outlist1 = ['wsp', 'pitch', 'rpm', 'P', 'Q', 'T', 'CP', 'CT', 'Mx',
                         'My', 'Mz', 'Fx', 'Fy']
        self.outlist2 = ['aoa', 'Ft', 'Fn', 'cl', 'cd', 'cm', 'ct', 'cp',
                         'v_a', 'v_t', 'disp_x', 'disp_y',
                         'disp_z', 'disp_rot_z']
        self.outlist3 = ['Fx_e', 'Fy_e', 'Fz_e', 'Mx_e', 'My_e', 'Mz_e',
                         'Fx_r', 'Fy_r', 'Fz_r', 'Mx_r', 'My_r', 'Mz_r']

    def execute(self):
        if ('save_power' not in self.commands) or \
           ('save_induction' not in self.commands):
            raise RuntimeError('HAWC2SOutput can only run if pwr and ind ' +
                               ' files have been computed.')
        super(HAWC2SOutput, self).execute()

        data = np.array(self.operational_data)

        self.wsp = data[:, 0]
        self.pitch = data[:, 1]
        self.rpm = data[:, 2]

        data = np.array(self.rotor_loads_data)
        self.P = data[:, 1] * 1000.
        self.Q = data[:, 1] * 1000. / (data[:, 9] * 2. * np.pi / 60.)
        self.T = data[:, 2] * 1000.
        self.CP = data[:, 3]
        self.CT = data[:, 4]
        self.Mx = data[:, 6] * 1000.
        self.My = data[:, 7] * 1000.
        self.Mz = data[:, 5]

        nW = len(self.wsp)
        nS = self.blade_loads_data[0].shape[0]

        self.tip_rot = np.zeros(nW)
        self.disp_x = np.zeros((nW, nS))
        self.disp_y = np.zeros((nW, nS))
        self.disp_z = np.zeros((nW, nS))
        self.disp_rot_z = np.zeros((nW, nS))
        self.aoa = np.zeros((nW, nS))
        self.Ft = np.zeros((nW, nS))
        self.Fn = np.zeros((nW, nS))
        self.cl = np.zeros((nW, nS))
        self.cd = np.zeros((nW, nS))
        self.cm = np.zeros((nW, nS))
        self.ct = np.zeros((nW, nS))
        self.cp = np.zeros((nW, nS))
        self.v_a = np.zeros((nW, nS))
        self.v_t = np.zeros((nW, nS))

        self.Fx = np.zeros(nW)
        self.Fy = np.zeros(nW)

        for iw, wsp in enumerate(self.wsp):
            try:
                data = self.blade_loads_data[iw]
                if len(data.shape) == 1:
                    data = data.reshape(1, data.shape[0])
                self.s = data[:, 0] / data[-1, 0]
                self.aoa[iw, :] = data[:, 4] * 180. / np.pi
                self.Ft[iw, :] = data[:, 6]
                self.Fn[iw, :] = data[:, 7]
                self.cl[iw, :] = data[:, 16]
                self.cd[iw, :] = data[:, 17]
                self.cm[iw, :] = data[:, 18]
                self.ct[iw, :] = data[:, 32]
                self.cp[iw, :] = data[:, 33]
                self.v_a[iw, :] = data[:, 26]
                self.v_t[iw, :] = data[:, 27]
                self.Fx[iw] = np.trapz(data[:, 6], x=data[:, 0])
                self.Fy[iw] = np.trapz(data[:, 7], x=data[:, 0])
                main_axis = data[:, 13:16]
                main_axis[:, 2] *= -1.
                self.disp_x[iw, :] = main_axis[:, 0]
                self.disp_y[iw, :] = main_axis[:, 1]
                self.disp_z[iw, :] = main_axis[:, 2]
                self.disp_rot_z[iw, :] = data[:, 28] * 180. / np.pi
            except:
                pass

        nfS = self.blade_fext_loads_data[0].shape[0]

        self.Fx_e = np.zeros((nW, nfS))
        self.Fy_e = np.zeros((nW, nfS))
        self.Fz_e = np.zeros((nW, nfS))
        self.Mx_e = np.zeros((nW, nfS))
        self.My_e = np.zeros((nW, nfS))
        self.Mz_e = np.zeros((nW, nfS))
        self.Fx_r = np.zeros((nW, nfS))
        self.Fy_r = np.zeros((nW, nfS))
        self.Fz_r = np.zeros((nW, nfS))
        self.Mx_r = np.zeros((nW, nfS))
        self.My_r = np.zeros((nW, nfS))
        self.Mz_r = np.zeros((nW, nfS))

        for iw, wsp in enumerate(self.wsp):
            try:
                data = self.blade_fext_loads_data[iw]
                if len(data.shape) == 1:
                    data = data.reshape(1, data.shape[0])
                self.s_e = data[:, 0] / data[-1, 0]
                self.Fx_e[iw, :] = data[:, 2]
                self.Fy_e[iw, :] = data[:, 3]
                self.Fz_e[iw, :] = data[:, 4]
                self.Mx_e[iw, :] = data[:, 5]
                self.My_e[iw, :] = data[:, 6]
                self.Mz_e[iw, :] = data[:, 7]
                self.Fx_r[iw, :] = data[:, 8]
                self.Fy_r[iw, :] = data[:, 9]
                self.Fz_r[iw, :] = data[:, 10]
                self.Mx_r[iw, :] = data[:, 11]
                self.My_r[iw, :] = data[:, 12]
                self.Mz_r[iw, :] = data[:, 13]
            except:
                pass


class HAWC2SOutputCompact(HAWC2SOutput):
    """
    HAWC2SOutput: HAWC2SOutputBase class that organize results in
    compact arrays to minimize the number of outputs. The outputs are selected
    with a dictionary passed in the initialization.

    Parameters
    ----------

    Returns
    -------
    outputs_rotor: array
        Rotor outputs.

    outputs_blade: array
        Blade outputs.
    outputs_blade: array
        Blade outputs.
    outputs_blade_fext: array
        Blade outputs.
    """
    def __init__(self, config):
        super(HAWC2SOutputCompact, self).__init__()

        self.sensor_rotor = []
        for name in config['rotor']:
            self.sensor_rotor.append(name)
            if name not in self.outlist1:
                raise RuntimeError('Rotor output required not supported.' +
                                   ' Wrong sensor: %s.' % name)
        self.sensor_blade = []
        for name in config['blade']:
            self.sensor_blade.append(name)
            if name not in self.outlist2:
                raise RuntimeError('Blade output required not supported.' +
                                   ' Wrong sensor: %s.' % name)

    def execute(self):
        super(HAWC2SOutputCompact, self).execute()

        nW = len(self.wsp)
        nS = self.blade_loads_data[0].shape[0]
        nfS = self.blade_fext_loads_data[0].shape[0] + 1

        self.outputs_rotor = np.zeros([nW, len(self.sensor_rotor)])
        for i, sensor in enumerate(self.sensor_rotor):
            self.outputs_rotor[:, i] = getattr(self, sensor)

        self.outputs_blade = np.zeros([nW, len(self.sensor_blade)*nS])
        for i, sensor in enumerate(self.sensor_blade):
            self.outputs_blade[:, i*nS:(i+1)*nS] = getattr(self, sensor)

        # fext loads are extrapolated to the inner-most section using polyfit
        self.outputs_blade_fext = np.zeros([nW, 6*nfS])
        for i, sensor in enumerate(['Fx_e', 'Fy_e', 'Fz_e', 'Mx_e', 'My_e', 'Mz_e']):
            d = getattr(self, sensor)
            f0 = np.polyval(np.polyfit(self.s_e[:2], d[0, :2], 1), 0.)
            self.outputs_blade_fext[:, i*nfS] = f0
            self.outputs_blade_fext[:, (i*nfS+1):(i+1)*nfS] = d[0, :]


class FreqDampTarget(object):
    """
    Component to compute th cost function for freqeuncies and dampings
    placement given the indexed of the modes.

    Parameters
    ----------
    freqdamp: array
        Two dimensional array containing the freqeuncies and dampings at
        different operational points. Same structure as in aeroelasticfreqdamp.

    mode_freq: array
        Two dimenstional array containing the target values of the frequencies
        at operational points. With the wind speed in the first column.

    mode_damp: array
        Two dimenstional array containing the target values of the dampings at
        operational points.  With the wind speed in the first column. It has
        to be of the same size as mode_freq.

    mode_target_freq: array
        Two deimensional array containing the target values of the modes
        freqeuncies.

    mode_target_damp: array
        Two deimensional array containing the target values of the modes
        dampings.

    Returns
    -------
    freq_factor: array
        RMS of the errors. Array with dimensions equal to freqdamp but without
        the column with the wind speeds.

    Example
    -------

    >>> freq = FreqDampTarget()
    >>> freq.freqdamp = output.aeroelasticfreqdamp
    >>> freq.mode_freq = np.array([[15., 0.25, 0.64],
    >>>                            [20., 0.25, 0.64]])
    >>> freq.mode_damp = np.array([[15., 0.8, 80.],
    >>>                            [20., 0.8, 80.]])
    >>> freq.mode_target_freq = np.array([[15., 0.3, -1],
    >>>                                   [20., 0.25, 0.6]])
    >>> freq.mode_target_damp = np.array([[15., -1,  85.],
    >>>                                   [20., 0.4, 80.]])
    >>> freq.execute()

    """
    def __init__(self, **kwargs):
        self.freqdamp = np.array([0.])
        self.mode_freq = np.array([0.])
        self.mode_damp = np.array([0.])
        self.mode_target_freq = np.array([0.])
        self.mode_target_damp = np.array([0.])

        for k, w in kwargs.iteritems():
            try:
                setattr(self, k, w)
            except:
                pass

    def mode_tracking_freq_damp(self):
        """
        Component to compute the indexes of modes for given frequencies and
        dampings
        """
        ws = self.mode_freq[:, 0]
        Nmodes = (self.freqdamp.shape[1]-1)/2  # remove the wind speed
        Nop = self.mode_freq.shape[0]
        Nm = self.mode_freq.shape[1]
        self.mode_index = np.zeros([Nop, Nm])
        iop = -1
        # Loop the operational points
        for freqdamp in self.freqdamp:
            # select right wind speed of reference values
            if freqdamp[0] in ws:
                iop += 1
                for w, iw in zip(ws, range(len(ws))):
                    if w == freqdamp[0]:
                        break
                mode_freq = self.mode_freq[iw, 1:]
                mode_damp = self.mode_damp[iw, 1:]

                allfreq = freqdamp[1:Nmodes+1]
                alldamp = freqdamp[Nmodes+1:]
                self.mode_index[iop, 0] = freqdamp[0]
                # Loop the modes to be tracked
                for im, (freq, damp) in enumerate(zip(mode_freq, mode_damp)):

                    err = np.sqrt((allfreq/freq - 1.)**2 +
                                  max(allfreq)/max(alldamp)*(alldamp/damp - 1.)**2)

                    if err[err.argmin()] > 1:
                        print 'Distance between computed mode freqeuncies ' +\
                              'and dampings and reference values for ' +\
                              'tracking is high! %f' % err[err.argmin()]

                    self.mode_index[iop, im+1] = int(err.argmin() + 1)

    def freq_damp_target_index(self, verbose=False):
        """
        Execute the computation of the objective function for frequency and
        dampings placement.
        """

        freq_factor = np.zeros([self.freqdamp.shape[0],
                                (self.mode_index.shape[1]-1)*2])

        Nmodes = int((self.freqdamp.shape[1]-1)/2)
        # for loop along the different operational points
        for i, freqdamp in enumerate(self.freqdamp):
            for iop, (mode_index, mode_target_freq, mode_target_damp) in \
                    enumerate(zip(self.mode_index, self.mode_target_freq,
                                  self.mode_target_damp)):
                if freqdamp[0] == mode_index[0]:
                    # for loop along the different target modes
                    for im, (index, target_freq, target_damp) in \
                            enumerate(zip(mode_index[1:], mode_target_freq[1:],
                                          mode_target_damp[1:])):

                        freq_factor[i, im] = abs(freqdamp[index] /
                                                 target_freq - 1)

                        freq_factor[i, im+len(mode_index[1:])] = \
                            abs(freqdamp[index+Nmodes] / target_damp - 1)
                        if verbose:
                            print 'Freq: ', freqdamp[index],\
                                  'Target Freq: ', target_freq,\
                                  'Diff.:', freq_factor[i, im]
                            print 'Damp: ', freqdamp[index+Nmodes],\
                                  'Target Damp:', target_damp,\
                                  'Diff.:', freq_factor[i, im+len(mode_index[1:])]
        self.freq_factor = freq_factor

    def execute(self):

        self.mode_tracking_freq_damp()
        self.freq_damp_target_index()


if __name__ == '__main__':

    pass
