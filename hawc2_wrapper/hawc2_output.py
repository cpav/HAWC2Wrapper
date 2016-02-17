"""
"""
import numpy as np
import os
import glob
import re

from hawc2_vartrees import DTUBasicControllerVT


class HAWC2Dataset(object):
    """

    """
    # TODO: link with pdap!!!!
    def __init__(self, casename=None, readonly=False):

        self.casename = casename
        self.readonly = readonly
        self.iknown = []
        self.data = None

        if os.path.isfile(self.casename + ".sel"):
            self._read_hawc2_sel_file()

        else:
            print "file not found: ", casename + '.sel'

    def _read_hawc2_sel_file(self):
        """
        Some title
        ==========

        Parameters
        ----------
        signal : ndarray
            some description

        Returns
        -------
        output : int
            describe variable
        """

        # read *.sel hawc2 output file for result info
        fid = open(self.casename + '.sel', 'r')
        lines = fid.readlines()
        fid.close()

        # find general result info (number of scans, number of channels,
        # simulation time and file format)
        temp = lines[8].split()
        self.scans = int(temp[0])
        self.channels = int(temp[1])
        self.total_time = float(temp[2])
        self.frequency = self.scans / self.total_time
        self.timestep = 1. / self.frequency
        self.time = np.linspace(0, self.total_time, self.scans)
        self.format = temp[3]
        # reads channel info (name, unit and description)
        name = []
        unit = []
        description = []
        for i in range(0, self.channels):
            temp = str(lines[i + 12][12:43])
            name.append(temp.strip())
            temp = str(lines[i + 12][43:50])
            unit.append(temp.strip())
            temp = str(lines[i + 12][54:100])
            description.append(temp.strip())
        self.chinfo = {'name': name, 'unit': unit, 'desc': description}
        # if binary file format, scaling factors are read
        if self.format.lower() == 'binary':
            self.scale_factor = np.zeros(self.channels)
            self.fileformat = 'HAWC2_BINARY'
            for i in range(0, self.channels):
                self.scale_factor[i] = float(lines[i + 12 + self.channels + 2])
        else:
            self.fileformat = 'HAWC2_ASCII'

    def _read_binary(self, channel_ids=[]):
        """Read results in binary format"""

        if not channel_ids:
            channel_ids = range(0, self.channels)
        fid = open(self.casename + '.dat', 'rb')
        data = np.zeros((self.scans, len(channel_ids)))
        j = 0
        for i in channel_ids:
            fid.seek(i * self.scans * 2, 0)
            data[:, j] = np.fromfile(fid, 'int16', self.scans) * self.scale_factor[i]
            j += 1
        fid.close()
        return data

    def _read_ascii(self, channel_ids=[]):
        """Read results in ASCII format"""

        if not channel_ids:
            channel_ids = range(0, self.channels)
        temp = np.loadtxt(self.casename + '.dat', usecols=channel_ids)
        return temp.reshape((self.scans, len(channel_ids)))

    def _read(self, channel_ids=[]):

        if not channel_ids:
            channel_ids = range(0, self.channels)
        if self.fileformat == 'HAWC2_BINARY':
            return self._read_binary(channel_ids)
        elif self.fileformat == 'HAWC2_ASCII':
            return self._read_ascii(channel_ids)

    def read_all(self):

        self.data = self._read()
        self.iknown = range(0, self.channels)

    def get_channels(self, channel_ids=[]):
        """  """

        if not channel_ids:
            channel_ids = range(0, self.channels)
        elif max(channel_ids) >= self.channels:
            print "channel number out of range"
            return

        if self.readonly:
            return self._read(channel_ids)

        else:
            # sort into known channels and channels to be read
            I1 = []
            I2 = []   # I1=Channel mapping, I2=Channels to be read
            for i in channel_ids:
                try:
                    I1.append(self.iknown.index(i))
                except:
                    self.iknown.append(i)
                    I2.append(i)
                    I1.append(len(I1))
            # read new channels
            if I2:
                temp = self._read(I2)
                # add new channels to Data
                if isinstance(self.data, np.ndarray):
                    self.data = np.append(self.data, temp, axis=1)
                # if first call, data is empty
                else:
                    self.data = temp
            return self.data[:, tuple(I1)]


class HAWC2OutputBase(object):

    def __init__(self):

        self.datasets = HAWC2Dataset()
        self.case_id = 'hawc2_case'
        self.res_directory = ''

    def execute(self):

        self._logger.info('reading outputs for case %s ...' % (self.case_id))
        # this is fragile since it will grab all files
        # which may be old ...

        self.datasets = []
        sets = glob.glob(self.res_directory+'\*.sel')

        for s in sets:
            name = re.sub("\.sel$", '', s)
            case = HAWC2Dataset(casename=name)
            case.read_all()
            self.datasets.append(case)


class HAWC2SOutputBase(object):
    """
    HAWC2SOutputBase: class that reads HAWC2s output files.

    parameters
    ----------
    case_id: str
        Name of solution to open.
    commands: list
        List containing the strings of the HAWC2s commands that have been
        executed. Only files associated with these commands are read.

    returns
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

    """
    def __init__(self):
        self.case_id = ''
        self.commands = []
        self.structuralfreqdamp = np.array([0.])
        self.aeroelasticfreqdamp = np.array([0.])
        self.aeroservoelasticfreqdamp = np.array([0.])
        self.controller_data = DTUBasicControllerVT()

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
            else:
                print 'Command "%s" not known.' % name


class HAWC2SOutput(HAWC2SOutputBase):
    """
    HAWC2SOutput: HAWC2SOutputBase class that organize results in more general
    arrays.

    parameters
    ----------

    returns
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


class HAWC2SOutputCompact(HAWC2SOutput):
    """
    HAWC2SOutput: HAWC2SOutputBase class that organize results in
    compact arrays to minimize the number of outputs. The outputs are selected
    with a dictionary passed in the initialization.

    parameters
    ----------

    returns
    -------
    outputs_rotor: array
        Rotor outputs.

    outputs_blade: array
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
        nS = len(self.blade_loads_data[0][:, 0])

        self.outputs_rotor = np.zeros([nW, len(self.sensor_rotor)])
        for i, sensor in enumerate(self.sensor_rotor):
            self.outputs_rotor[:, i] = getattr(self, sensor)

        self.outputs_blade = np.zeros([nW, len(self.sensor_blade)*nS])
        for i, sensor in enumerate(self.sensor_blade):
            self.outputs_blade[:, i*nS:(i+1)*nS] = getattr(self, sensor)


class FreqDampTargetByIndex(object):
    """
    Component to compute th cost function for freqeuncies and dampings
    placement given the indexed of the modes

    parameters:
    -----------
    freqdamp: array
        Two dimensional array containing the freqeuncies and dampings at
        different operational points.

    mode_index: array
        Two dimensional array containing the indexed of the freqeuncies and
        dampings to be placed at different operational points.

    mode_target_freq: array
        Two dimenstional array containing the target values of the freqeuncies
        at operational points. Has to be of the same size as mode_index.

    mode_target_damp: array
        Two dimenstional array containing the target values of the dampings at
        operational points. Has to be of the same size as mode_index.

    freq_factor: list
        RMS of the errors

    example
    --------

    """
    def __init__(self):
        self.freqdamp = np.array([0.])
        self.mode_index = np.array([0.])
        self.mode_target_freq = np.array([0.])
        self.mode_target_damp = np.array([0.])

        self.freq_factor = []

    def execute(self):
        """
        Execute the computation of the objective function for frequency and
        dampings placement.
        """
        freq_factor = []
        Nmodes = (self.freqdamp.shape[1]-1)/2
        # for loop along the different operational points
        for freqdamp in self.freqdamp:
            for mode_index, mode_target_freq, mode_target_damp in \
                    zip(self.mode_index, self.mode_target_freq,
                        self.mode_target_damp):
                if freqdamp[0] == mode_index[0]:
                    # for loop along the different target modes
                    for index, target_freq, target_damp in \
                            zip(mode_index[1:], mode_target_freq[1:],
                                mode_target_damp[1:]):

                        if target_freq != -1:
                            freq_factor.append(abs(freqdamp[index] /
                                                   target_freq - 1))
                            print 'Freq: ', freqdamp[index],\
                                  'Target Freq: ', target_freq,\
                                  'Diff.:', freq_factor[-1]
                        if target_damp != -1:
                            freq_factor.append(abs(freqdamp[index+Nmodes] /
                                                   target_damp - 1))
                            print 'Damp: ', freqdamp[index+Nmodes],\
                                  'Target Damp:', target_damp,\
                                  'Diff.:', freq_factor[-1]
        self.freq_factor.freq_factor = freq_factor


class ModeTrackingByFreqDamp(object):
    """
    Component to compute the indexes of modes for given frequencies and
    dampings

    parameters:
    -----------
    freqdamp: array
        One dimensional array containing the freqeuncies and dampings.

    mode_freq: array
        One dimenstional array containing the reference values of the mode
        freqeuncies.

    mode_damp: array
        One dimenstional array containing the reference values of the mode
        dampings.The dampings has to be of the same mode indicated by the
        freqeuncies in mode_freq.

    mode_index: array
        One dimensional array containing the indexed of the tracked modes.

    example
    --------

    """
    def __init__(self):
        self.freqdamp = np.array([0.])
        self.mode_freq = np.array([0.])
        self.mode_damp = np.array([0.])

        self.mode_index = np.array([0.])

    def execute(self):
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
                for freq, damp, im in zip(mode_freq, mode_damp, range(1, Nm)):

                    err = np.sqrt(((allfreq - freq) / freq)**2 +
                                  ((alldamp - damp) / damp)**2)
                    if err[err.argmin()] > 1:
                        print 'Distance between computed mode freqeuncies ' +\
                              'and dampings and reference values for ' +\
                              'tracking is high! %f' % err[err.argmin()]

                    self.mode_index[iop, im] = err.argmin() + 1


class FreqDampTarget(object):
    """
    Component to compute th cost function for freqeuncies and dampings
    placement given the indexed of the modes

    parameters:
    -----------
    freqdamp: array
        Two dimensional array containing the freqeuncies and dampings at
        different operational points.

    mode_freq: array
        Two dimenstional array containing the target values of the frequencies
        at operational points.

    mode_damp: array
        Two dimenstional array containing the target values of the dampings at
        operational points. Has to be of the same size as mode_freq.

    results
    -------
    freq_factor: array
        RMS of the errors

    example
    --------

    """
    def __init__(self):
        self.freqdamp = np.array([0.])
        self.mode_freq = np.array([0.])
        self.mode_damp = np.array([0.])
        self.mode_target_freq = np.array([0.])
        self.mode_target_damp = np.array([0.])
        self.freq_factor = []

    def execute(self):

        modetrack = ModeTrackingByFreqDamp()
        modetrack.mode_freq = self.mode_freq
        modetrack.mode_damp = self.mode_damp
        modetrack.freqdamp = self.freqdamp
        modetrack.execute()

        freqtarget = FreqDampTargetByIndex()
        freqtarget.mode_target_freq = self.mode_target_freq
        freqtarget.mode_target_damp = self.mode_target_damp
        freqtarget.freqdamp = self.freqdamp
        freqtarget.mode_index = modetrack.mode_index
        freqtarget.execute()

        self.freq_factor = freqtarget.freq_factor
