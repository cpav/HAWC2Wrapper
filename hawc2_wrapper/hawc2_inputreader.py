""""""
from hawc2_inputdict import HAWC2InputDict, read_hawc2_st_file,\
                            read_hawc2_pc_file, read_hawc2_ae_file
from hawc2_vartrees import HAWC2VarTrees, HAWC2Simulation, HAWC2Wind,\
                           HAWC2Aero, HAWC2AirfoilDataset, HAWC2AirfoilPolar,\
                           HAWC2AeroDrag, HAWC2AeroDragElement, HAWC2MainBody,\
                           HAWC2BeamStructure, HAWC2Type2DLL, HAWC2OutputVT,\
                           HAWC2SBody

from numpy import array, loadtxt


class HAWC2InputReader(object):
    """
    Class to read HAWC2 files and store the data in variables trees.

    Parameters
    ----------
    htc_master_file: str
        Name of the htc file to read.

    Creates
    -------
    vartrees: HAWC2VarTrees
        Attribute corresponding to the variable tree with all the information
        of the model.

    Example
    -------
    >>> from hawc2_inputreader import HAWC2InputReader
    >>> reader = HAWC2InputReader('hawc2_master.htc')
    >>> reader.execute()

    """
    def __init__(self, htc_master_file='hawc_master.htc'):
        self.htc_master_file = htc_master_file
        self.vartrees = HAWC2VarTrees()

    def execute(self):
        self.dict = HAWC2InputDict()
        self.dict.read(self.htc_master_file)
        self.htc = self.dict.htc
        self.vartrees.body_order = self.dict.body_order

        for section in self.htc:
            if section.name == 'simulation':
                self._add_simulation(section)
            elif section.name == 'wind':
                self._add_wind(section)
            elif section.name == 'aero':
                self._add_aero(section)
            elif section.name == 'aerodrag':
                self._add_aerodrag(section)
            elif section.name == 'new_htc_structure':
                self._add_structure(section)
            elif section.name == 'output':
                self._add_output(section)
            elif section.name == 'dll':
                self._add_dlls(section)
            elif section.name == 'hawcstab2':
                self._add_hawcstab2(section)

    def set_entry(self, vt, section, name, h2name=None, required=False):

        if h2name is None:
            var = section.get_entry(name)

        else:
            var = section.get_entry(h2name)

        if var is None and required:
            raise RuntimeError('Missing input variable %s in section %s' %
                               (name, section.name))

        elif var is not None:
            if isinstance(getattr(vt, name), list) and isinstance(var, str):
                setattr(vt, name, [var])

            else:
                setattr(vt, name, var)

        return vt

    def read_operational_data_file(self):
        data = loadtxt(self.vartrees.h2s.operational_data_filename, skiprows=1)
        if len(data.shape) == 1:
            data = data.reshape([1, data.shape[0]])
        self.vartrees.h2s.wsp_curve = data[:, 0]
        self.vartrees.h2s.pitch_curve = data[:, 1]
        self.vartrees.h2s.rpm_curve = data[:, 2]

    def add_pc_data(self):

        pcdata = read_hawc2_pc_file(self.vartrees.aero.pc_filename)

        desc = pcdata[0]
        data = pcdata[1]

        self.vartrees.airfoildata.nset = len(data)
        self.vartrees.airfoildata.desc = desc

        for dataset in data:
            pcset = HAWC2AirfoilDataset()
            pcset.np = len(dataset['polars'])
            rthick = []
            for p in dataset['polars']:
                polar = HAWC2AirfoilPolar()
                polar.rthick = p['rthick']
                polar.desc = p['desc']
                polar.aoa = p['aoa']
                polar.cl = p['cl']
                polar.cd = p['cd']
                polar.cm = p['cm']
                rthick.append(polar.rthick)
                pcset.polars.append(polar)
            pcset.rthick = rthick
            self.vartrees.airfoildata.pc_sets.append(pcset)

    def add_ae_data(self):

        blade_ae = read_hawc2_ae_file(self.vartrees.aero.ae_filename)

        for var in ['s', 'chord', 'rthick', 'aeset']:
            setattr(self.vartrees.blade_ae, var, blade_ae[var])

    def _add_simulation(self, section):

        vt = HAWC2Simulation()
        for var in vt.var:
            vt = self.set_entry(vt, section, var)

        newmark = section.get_entry('newmark')
        vt.newmark_deltat = newmark.get_entry('deltat')
        vt = self.set_entry(vt, section, 'eig_out')

        self.vartrees.sim = vt

    def _add_wind(self, section):

        vt = HAWC2Wind()
        vt = self.set_entry(vt, section, 'density')
        vt = self.set_entry(vt, section, 'wsp', required=True)
        vt = self.set_entry(vt, section, 'tint')
        vt = self.set_entry(vt, section, 'horizontal_input')
        vt = self.set_entry(vt, section, 'center_pos0')
        vt = self.set_entry(vt, section, 'windfield_rotations')
        vt = self.set_entry(vt, section, 'tint')
        vt = self.set_entry(vt, section, 'scale_time_start')
        shear = section.get_entry('shear_format')
        vt.shear_type = shear[0]
        vt.shear_factor = shear[1]
        vt = self.set_entry(vt, section, 'tower_shadow_method')
        vt = self.set_entry(vt, section, 'turb_format')
        pot = section.get_entry('tower_shadow_potential_2')
        if pot is not None:
            vt.add_shadow('tower_potential')
            vt.tower_potential.tower_mbdy_link =\
                pot.get_entry('tower_mbdy_link')
            vt.tower_potential.nsec = pot.get_entry('nsec')
            vt.tower_potential.sections = pot.get_entry('radius')

        ramp = section.get_entry('wind_ramp_factor')
        if ramp is not None:
            vt.wind_ramp_t0 = ramp[0]
            vt.wind_ramp_t1 = ramp[1]
            vt.wind_ramp_factor0 = ramp[2]
            vt.wind_ramp_factor1 = ramp[3]
        ramp = section.get_entry('wind_ramp_abs')
        if ramp is not None:
            if type(ramp[0]) is not list:
                vt.wind_ramp_abs = [ramp]
            else:
                for rm in ramp:
                    vt.wind_ramp_abs.append(rm)

        gust = section.get_entry('iec_gust')
        if gust is not None:
            vt.iec_gust = True
            vt.iec_gust_type = gust[0]
            vt.G_A = gust[1]
            vt.G_phi0 = gust[2]
            vt.G_t0 = gust[3]
            vt.G_T = gust[4]
        if vt.turb_format == 1:
            mann = section.get_entry('mann')
            vt.add_turbulence('mann')

            temp = mann.get_entry('create_turb_parameters')
            if temp is not None:
                vt.mann.create_turb = True
                vt.mann.L = temp[0]
                vt.mann.alfaeps = temp[1]
                vt.mann.gamma = temp[2]
                vt.mann.seed = temp[3]
                vt.mann.highfrq_compensation = temp[4]
            turb = mann.get_entry('filename_u').split('/')
            if len(turb) == 1:
                turb = mann.get_entry('filename_u').split('\\')
                if len(turb) == 1:
                    pass
                else:
                    vt.turb_directory = '\\'.join(turb[:-1])
            else:
                vt.turb_directory = '/'.join(turb[:-1])

            vt.mann.turb_base_name = turb[-1].strip('u.bin')
            vt.mann.box_nu = mann.get_entry('box_dim_u')[0]
            vt.mann.box_nv = mann.get_entry('box_dim_v')[0]
            vt.mann.box_nw = mann.get_entry('box_dim_w')[0]
            vt.mann.box_du = mann.get_entry('box_dim_u')[1]
            vt.mann.box_dv = mann.get_entry('box_dim_v')[1]
            vt.mann.box_dw = mann.get_entry('box_dim_w')[1]
            vt.mann.std_scaling = mann.get_entry('std_scaling')

        self.vartrees.wind = vt

    def _add_aero(self, section):

        vt = HAWC2Aero()
        for var in vt.var:
            vt = self.set_entry(vt, section, var)

        aero_dist = section.get_entry('aero_distribution')
        if aero_dist is not None:
            vt.aero_distribution_file = aero_dist[0]
            vt.aero_distribution_set = int(aero_dist[1])

        hub_vec = section.get_entry('hub_vec')
        if hub_vec is not None:
            vt.hub_vec_mbdy_name = hub_vec[0]
            vt.hub_vec_coo = hub_vec[1]
        for b in range(vt.nblades):
            link = [b + 1, 'mbdy_c2_def', 'blade%i' % (b + 1)]
            vt.links.append(link)

        self.vartrees.aero = vt

        # read ae and pc data
        self.add_pc_data()
        self.add_ae_data()

    def _add_aerodrag(self, section):

        vt = HAWC2AeroDrag()
        for entry in section.entries:
            e = HAWC2AeroDragElement()
            e = self.set_entry(e, entry, 'mbdy_name')
            dist = entry.get_entry('aerodrag_sections')
            e.dist = dist[0]
            e.calculation_points = int(dist[1])
            e = self.set_entry(e, entry, 'sections', h2name='sec')
            vt.elements.append(e)

        self.vartrees.aerodrag = vt

    def _add_structure(self, section):

        for sec in section.entries:
            if sec.name == 'main_body':
                b = self._add_main_body(sec)
                self.vartrees.main_bodies.add_main_body(b.body_name, b)

            elif sec.name == 'orientation':
                self._add_orientations(sec)

            elif sec.name == 'constraint':
                self._add_constraints(sec)

    def _add_main_body(self, section):

        b = HAWC2MainBody()
        b = self.set_entry(b, section, 'body_name', h2name='name',
                           required=True)

        copy = section.get_entry('copy_main_body')
        if copy is not None:
            b.copy_main_body = copy
            return b
        b = self.set_entry(b, section, 'body_type', h2name='type')
        b = self.set_entry(b, section, 'nbodies')
        b = self.set_entry(b, section, 'node_distribution')

        d_ani = section.get_entry('damping_aniso')
        if d_ani is not None:
            b.damping_aniso = d_ani
            b.damping_type = 'ani'
        else:
            b = self.set_entry(b, section, 'damping_posdef')

        cm = section.get_entry('concentrated_mass')
        if cm is not None:
            if isinstance(cm[0], (int, float)):
                cm = [cm]
            b.concentrated_mass = cm
        timo = section.get_entry('timoschenko_input')
        st_type = timo.get_entry('FPM')
        b.body_set = timo.get_entry('set')

        b.st_input_type = 0
        if st_type is not None:
            b.st_input_type = st_type

        st = HAWC2BeamStructure(b.st_input_type)
        stdic = read_hawc2_st_file(timo.get_entry('filename'),
                                   st.var, b.body_set[0])
        b.body_set[0] = 1

        for stset in stdic:
            st = HAWC2BeamStructure(b.st_input_type)
            for k, w in stset.iteritems():
                setattr(st, k, w)
            b.beam_structure.append(st)
        c2def = section.get_entry('c2_def')
        b.c12axis = array(c2def.get_entry('sec'))[:, 1:5]
        return b

    def _add_orientations(self, section):

        for sec in section.entries:
            if sec.name == 'base':
                body_name = sec.get_entry('body')
                if body_name is None:
                    # try new input name
                    body_name = sec.get_entry('mbdy')
                b = self.vartrees.main_bodies.get_main_body(body_name)
                b.orientations = []
                o = b.add_orientation(sec.name)
                o = self.set_entry(o, sec, 'inipos')
                orien = sec.get_entry('body_eulerang')
                if orien is not None:
                    if isinstance(orien[0], float):
                        orien = [orien]
                    o.body_eulerang = orien

            elif sec.name == 'relative':
                body1 = sec.get_entry('body1')
                body2 = sec.get_entry('body2')
                if body1[0] is None:
                    # try new input name
                    body1 = sec.get_entry('mbdy1')
                    body2 = sec.get_entry('mbdy2')
                b = self.vartrees.main_bodies.get_main_body(body2[0])
                o = b.add_orientation(sec.name)
                o = self.set_entry(o, sec, 'body1')
                o = self.set_entry(o, sec, 'body2')
                orien = sec.get_entry('body2_eulerang')
                if orien is not None:
                    if isinstance(orien[0], float):
                        orien = [orien]
                    o.body2_eulerang = orien
                o = self.set_entry(o, sec, 'mbdy2_ini_rotvec_d1',
                                   h2name='body2_ini_rotvec_d1')
                o = self.set_entry(o, sec, 'mbdy2_ini_rotvec_d1')
                # these inputs aren't used in HAWC2 as far as I know...
                o = self.set_entry(o, sec, 'initial_speed')
                o = self.set_entry(o, sec, 'rotation_dof')

    def _add_constraints(self, section):

        for sec in section.entries:
            if sec.name in ['fix0', 'fix2', 'fix3']:
                body_name = sec.get_entry('body')
                if body_name is None:
                    # try new input name
                    body_name = sec.get_entry('mbdy')
                b = self.vartrees.main_bodies.get_main_body(body_name)
                c = b.add_constraint_new(sec.name)
                c.con_type = sec.name
                c.mbdy = body_name
                if sec.name in ['fix2', 'fix3']:
                    c.dof = self.set_entry(c, sec, 'dof')

            elif sec.name in ['fix1', 'fix4']:
                # body1 and body2 need to be lists to deal with either
                # ["string" "string"] or ["string" "int"] used in
                # HAWC2 constraints

                body1 = sec.get_entry('body1')
                body2 = sec.get_entry('body2')
                if body1[0] is None:
                    # try new input name
                    body1 = sec.get_entry('mbdy1')
                    body2 = sec.get_entry('mbdy2')
                b = self.vartrees.main_bodies.get_main_body(body2[0])
                c = b.add_constraint_new(sec.name)
                c.con_type = sec.name
                c = self.set_entry(c, sec, 'mbdy1')
                c = self.set_entry(c, sec, 'mbdy2')
                c.mbdy1 = body1
                c.mbdy2 = body2
                if sec.name == 'fix1':
                    c = self.set_entry(c, sec, 'disable_at')
                elif sec.name == 'fix4':
                    c = self.set_entry(c, sec, 'time')

            elif 'bearing' in sec.name:
                # body1 and body2 need to be lists to deal with either
                # ["string" "string"] or ["string" "int"] used in
                # HAWC2 constraints

                body1 = sec.get_entry('body1')
                body2 = sec.get_entry('body2')
                if body1[0] is None:
                    # try new input name
                    body1 = sec.get_entry('mbdy1')
                    body2 = sec.get_entry('mbdy2')
                b = self.vartrees.main_bodies.get_main_body(body2[0])
                c = b.add_constraint_new(sec.name)
                c.mbdy1 = body1
                c.mbdy2 = body2
                c = self.set_entry(c, sec, 'name')
                c = self.set_entry(c, sec, 'bearing_vector')
                if sec.name == 'bearing3':
                    c = self.set_entry(c, sec, 'omegas')
                else:
                    c = self.set_entry(c, sec, 'disable_at')

    def _add_dlls(self, section):

        for sec in section.entries:
            if sec.name == 'type2_dll':
                dll = self._add_type2_dll(sec)
                self.vartrees.dlls.add_dll(dll.name, dll)
                self.vartrees.dlls_order.append(dll.name)
            elif sec.name == 'hawc_dll':
                raise NotImplemented('%s: hawc_dll type not implemented, \
                                     use type2 dlls' % sec.get_entry('name'))

    def _add_type2_dll(self, sec):

        dll = HAWC2Type2DLL()
        dll = self.set_entry(dll, sec, 'name', required=True)
        dll = self.set_entry(dll, sec, 'filename', required=True)
        dll = self.set_entry(dll, sec, 'dll_subroutine_init', required=True)
        dll = self.set_entry(dll, sec, 'dll_subroutine_update', required=True)
        dll = self.set_entry(dll, sec, 'arraysizes_init', required=True)
        dll = self.set_entry(dll, sec, 'arraysizes_update', required=True)
        dll = self.set_entry(dll, sec, 'deltat')
        dll.set_init(dll.name)
        constants = sec.get_entry('init').entries

        dll.dll_init.set_constants(constants)

        io = sec.get_entry('output').entries
        dll.output.set_outputs(io)

        if sec.get_entry('actions') is not None:
            io = sec.get_entry('actions').entries
            dll.actions.set_actions(io)
        return dll

    def _add_output(self, section):

        o = HAWC2OutputVT()
        o = self.set_entry(o, section, 'filename', required=True)
        o = self.set_entry(o, section, 'out_buffer', h2name='buffer')
        o = self.set_entry(o, section, 'out_format', h2name='data_format')
        time = section.get_entry('time')
        o.time_start = time[0]
        o.time_stop = time[1]
        o.set_outputs(section.entries)
        self.vartrees.output = o

    def _add_hawcstab2(self, section):

        dll = HAWC2Type2DLL()
        dll.set_init('dtu_we_controller')
        self.vartrees.dlls.add_dll('dtu_we_controller', dll)
        for sec in section.entries:
            if sec.name == 'ground_fixed_substructure':
                self.vartrees.h2s.ground_fixed = self._add_hawc2s_body(sec)

            elif sec.name == 'rotating_axissym_substructure':
                self.vartrees.h2s.rotating_axissym =\
                    self._add_hawc2s_body(sec)

            elif sec.name == 'rotating_threebladed_substructure':
                self.vartrees.h2s.rotating_threebladed =\
                    self._add_hawc2s_body(sec)
                try:
                    b = sec.get_entry('second_order_actuator')
                    self.vartrees.h2s.second_order_actuator.name = b[0]
                    self.vartrees.h2s.second_order_actuator.frequency = b[1]
                    self.vartrees.h2s.second_order_actuator.damping = b[2]
                except:
                    pass

            elif sec.name == 'operational_data':
                self._add_operational_data(sec)

            elif sec.name == 'controller_tuning':
                self._add_controller_tuning(sec)

            elif sec.name == 'controller':
                self._add_controller(sec)

            elif sec.name == 'print_full_precision':
                self.vartrees.h2s.commands.append('print_full_precision')

            elif sec.name == 'compute_optimal_pitch_angle':
                self.vartrees.h2s.commands.append('compute_optimal_pitch_angle')

            elif sec.name == 'degrees_of_freedom':
                self.vartrees.h2s.commands.append(
                    'degrees_of_freedom' + 5*' %s' % tuple(sec.val[:5]))

            elif sec.name == 'steady_state_convergence_limits':
                self.vartrees.h2s.commands.append(
                    'steady_state_convergence_limits' + 9*' %g'
                    % tuple(sec.val[:9]))

            elif sec.name == 'compute_steady_states':
                self.vartrees.h2s.commands.append('compute_steady_states')
                self.vartrees.h2s.options.bladedeform = sec.val[0]
                self.vartrees.h2s.options.tipcorrect = sec.val[1]
                self.vartrees.h2s.options.induction = sec.val[2]
                self.vartrees.h2s.options.gradients = sec.val[3]

            elif sec.name == 'compute_steadystate':
                self.vartrees.h2s.commands.append('compute_steadystate')
                self.vartrees.h2s.options.bladedeform = sec.val[0]
                self.vartrees.h2s.options.tipcorrect = sec.val[1]
                self.vartrees.h2s.options.induction = sec.val[2]
                self.vartrees.h2s.options.gradients = sec.val[3]

            elif sec.name == 'compute_stability_analysis':
                self.vartrees.h2s.commands.append('compute_stability_analysis')
                self.vartrees.h2s.options.matrixwriteout = sec.val[0]
                self.vartrees.h2s.options.eigenvaluewriteout = sec.val[1]
                self.vartrees.h2s.options.number_of_modes = sec.val[2]
                self.vartrees.h2s.options.maximum_damping = sec.val[3]
                self.vartrees.h2s.options.minimum_frequency = sec.val[4]
                self.vartrees.h2s.options.zero_pole_threshold = sec.val[5]
                self.vartrees.h2s.options.aero_deflect_ratio = sec.val[6]
                self.vartrees.h2s.options.frequencysorting = sec.val[7]

            elif sec.name == 'compute_aeroservoelastic':
                self.vartrees.h2s.commands.append('compute_aeroservoelastic')
                self.vartrees.h2s.options.matrixwriteout = sec.val[0]
                self.vartrees.h2s.options.eigenvaluewriteout = sec.val[1]
                self.vartrees.h2s.options.number_of_modes = sec.val[2]
                self.vartrees.h2s.options.maximum_damping = sec.val[3]
                self.vartrees.h2s.options.minimum_frequency = sec.val[4]
                self.vartrees.h2s.options.zero_pole_threshold = sec.val[5]
                self.vartrees.h2s.options.aero_deflect_ratio = sec.val[6]
                self.vartrees.h2s.options.frequencysorting = sec.val[7]

            elif sec.name == 'basic_dtu_we_controller':
                self.vartrees.h2s.commands.append('basic_dtu_we_controller')
                dll_init = self.vartrees.dlls.dtu_we_controller.dll_init
                dll_init.pgTorque = sec.val[0]
                dll_init.igTorque = sec.val[1]
                dll_init.Qg = sec.val[2]
                dll_init.pgPitch = sec.val[3]
                dll_init.igPitch = sec.val[4]
                dll_init.KK1 = sec.val[5]
                dll_init.KK2 = sec.val[6]
                dll_init.generatorFreq = sec.val[7]
                dll_init.generatorDamping = sec.val[8]
                dll_init.ffFreq = sec.val[9]

                if sec.val[10] == 1:
                    dll_init.generatorSwitch = 1

                elif sec.val[10] == 0:
                    dll_init.generatorSwitch = 2

                if len(sec.val) > 11:
                    dll_init.Kp2 = sec.val[11]
                    dll_init.Ko1 = sec.val[12]
                    dll_init.Ko2 = sec.val[13]

            elif sec.name == 'compute_controller_input':
                self.vartrees.h2s.commands.append('compute_controller_input')

            elif sec.name == 'save_beam_data':
                self.vartrees.h2s.commands.append('save_beam_data')

            elif sec.name == 'save_blade_geometry':
                self.vartrees.h2s.commands.append('save_blade_geometry')

            elif sec.name == 'save_aero_point_data':
                self.vartrees.h2s.commands.append('save_aero_point_data')

            elif sec.name == 'save_profile_coeffs':
                self.vartrees.h2s.commands.append('save_profile_coeffs')

            elif sec.name == 'save_power':
                self.vartrees.h2s.commands.append('save_power')

            elif sec.name == 'save_induction':
                self.vartrees.h2s.commands.append('save_induction')

            elif sec.name == 'save_ol_matrices':
                self.vartrees.h2s.commands.append('save_ol_matrices')

            elif sec.name == 'save_cl_matrices_all':
                self.vartrees.h2s.commands.append('save_cl_matrices_all')
                self.vartrees.h2s.options.vloc_out = sec.val == 'vloc_out'

            elif sec.name == 'compute_structural_modal_analysis':
                self.vartrees.h2s.commands.append(
                    'compute_structural_modal_analysis')
                self.vartrees.h2s.options.blade_only =\
                    sec.val[0] == 'bladeonly'
                self.vartrees.h2s.options.number_of_modes = sec.val[1]

        self.set_entry(self.vartrees.h2s, section, 'operational_data_filename')

        try:
            self.read_operational_data_file()
        except:
            print 'Warning: couldn''t read operational file %s' % \
                self.vartrees.h2s.operational_data_filename

    def _add_hawc2s_body(self, section):

        b = HAWC2SBody()
        b = self.set_entry(b, section, 'main_body')
        b = self.set_entry(b, section, 'log_decrements')
        return b

    def _add_operational_data(self, section):

        c = section.get_entry('windspeed')
        self.vartrees.dlls.dtu_we_controller.dll_init.Vin = c[0]
        self.vartrees.dlls.dtu_we_controller.dll_init.Vout = c[1]
        self.vartrees.dlls.dtu_we_controller.dll_init.nV = c[2]
        c = section.get_entry('genspeed')
        self.vartrees.dlls.dtu_we_controller.dll_init.minRPM = c[0]
        self.vartrees.dlls.dtu_we_controller.dll_init.maxRPM = c[1]

        self.vartrees.dlls.dtu_we_controller.dll_init = \
            self.set_entry(self.vartrees.dlls.dtu_we_controller.dll_init,
                           section, 'ratedAeroPower', h2name='maxpow')

        self.vartrees.dlls.dtu_we_controller.dll_init = \
            self.set_entry(self.vartrees.dlls.dtu_we_controller.dll_init,
                           section, 'designTSR', h2name='opt_lambda')

        self.vartrees.dlls.dtu_we_controller.dll_init = \
            self.set_entry(self.vartrees.dlls.dtu_we_controller.dll_init,
                           section, 'minPitch', h2name='minpitch')

        self.vartrees.dlls.dtu_we_controller.dll_init = \
            self.set_entry(self.vartrees.dlls.dtu_we_controller.dll_init,
                           section, 'gearRatio', h2name='gearratio')

        self.vartrees.dlls.dtu_we_controller.dll_init = \
            self.set_entry(self.vartrees.dlls.dtu_we_controller.dll_init,
                           section, 'prvs_turbine')

        self.vartrees.h2s.options = \
            self.set_entry(self.vartrees.h2s.options,
                           section, 'include_torsiondeform')

        self.vartrees.h2s.options = \
            self.set_entry(self.vartrees.h2s.options,
                           section, 'set_torque_limit')

    def _add_controller_tuning(self, section):

        c = section.get_entry('partial_load')
        self.vartrees.dlls.dtu_we_controller.dll_init.poleFreqTorque = c[0]
        self.vartrees.dlls.dtu_we_controller.dll_init.poleDampTorque = c[1]
        c = section.get_entry('full_load')
        self.vartrees.dlls.dtu_we_controller.dll_init.poleFreqPitch = c[0]
        self.vartrees.dlls.dtu_we_controller.dll_init.poleDampPitch = c[1]

        self.vartrees.dlls.dtu_we_controller.dll_init = \
            self.set_entry(self.vartrees.dlls.dtu_we_controller.dll_init,
                           section, 'gainScheduling', h2name='gain_scheduling')

        self.vartrees.dlls.dtu_we_controller.dll_init = \
            self.set_entry(self.vartrees.dlls.dtu_we_controller.dll_init,
                           section, 'generatorSwitch', h2name='constant_power')

        self.vartrees.dlls.dtu_we_controller.dll_init = \
            self.set_entry(self.vartrees.dlls.dtu_we_controller.dll_init,
                           section, 'rotorspeed_gs', h2name='rotorspeed_gs')

        self.vartrees.h2s.options = \
            self.set_entry(self.vartrees.h2s.options,
                           section, 'regions', h2name='regions')

    def _add_controller(self, section):

        o = HAWC2OutputVT()
        i = section.get_entry('input')
        o.set_outputs(i.entries)
        self.vartrees.h2s.ch_list_in = o

        o = HAWC2OutputVT()
        i = section.get_entry('output')
        o.set_outputs(i.entries)
        self.vartrees.h2s.ch_list_out = o

if __name__ == '__main__':

    pass
