import os
import numpy as np
from hawc2_vartrees import HAWC2VarTrees


def _get_fmt(val):

    def _get_fmt1(val):
        if isinstance(val, str):
            return '%s'
        elif isinstance(val, int):
            return '%i'
        elif isinstance(val, float):
            return '%24.15f'

    if isinstance(val, list):
        return ' '.join([_get_fmt1(v) for v in val])
    else:
        return _get_fmt1(val)


def write_pcfile(path, pc):
    """write the blade airfoils"""

    fid = open(path, 'w')
    fid.write('%i %s\n' % (pc.nset, pc.desc))
    for i in range(pc.nset):
        pcset = pc.pc_sets[i]
        fid.write('%i\n' % pcset.np)
        for j in range(pcset.np):
            polar = pcset.polars[j]
            fid.write('%i %i %f %s\n' %
                      (j + 1, polar.aoa.shape[0], polar.rthick, polar.desc))
            for k in range(polar.aoa.shape[0]):
                fid.write((4*'%23.15e '+'\n') % (polar.aoa[k], polar.cl[k],
                                               polar.cd[k], polar.cm[k]))
    fid.close()


def write_aefile(path, b):
    """write the blade shape to the ae_filename"""

    fid = open(path, 'w')

    fid.write('1 HAWC_AE data\n')
    fid.write("1 %i r,c,t/c prin.set\n" % (b.s.shape[0]))
    data = np.array([b.s,
                     b.chord,
                     np.minimum(100., b.rthick),
                     b.aeset]).T
    np.savetxt(fid, data, fmt="%23.15e %23.15e %23.15e %i")
    fid.close()


def write_stfile(body, case_id):

    """write the beam structural data to an st_filename"""
    exp_prec = 15             # exponential precesion
    col_width = exp_prec + 9  # column width required for exp precision
    if body.st_input_type is 0:
        header = ['r', 'm', 'x_cg', 'y_cg', 'ri_x', 'ri_y', 'x_sh', 'y_sh',
                  'E', 'G', 'I_x', 'I_y', 'K', 'k_x', 'k_y', 'A', 'pitch',
                  'x_e', 'y_e']
        # for readable files with headers above the actual data column

        header_full = '='*20*col_width + '\n'
        header_full += ''.join([(hh + ' [%i]').center(col_width + 1) % i
                               for i, hh in enumerate(header)])+'\n'
        header_full += '='*20*col_width + '\n'

    else:
        header = ['r', 'm', 'x_cg', 'y_cg', 'ri_x', 'ri_y', 'pitch', 'x_e',
                  'y_e', 'K_11', 'K_12', 'K_13', 'K_14', 'K_15', 'K_16',
                  'K_22', 'K_23', 'K_24', 'K_25', 'K_26', 'K_33', 'K_34',
                  'K_35', 'K_36', 'K_44', 'K_45', 'K_46', 'K_55', 'K_56',
                  'K_66']
        # for readable files with headers above the actual data column
        header_full = '='*32*col_width + '\n'
        header_full += ''.join([(hh + ' [%i]').center(col_width + 1) % i
                                for i, hh in enumerate(header)])+'\n'
        header_full += '='*32*col_width + '\n'

    fid = open(body.st_filename, 'w')
    fid.write('%i  number of sets, Nset\n' % body.body_set[1])
    fid.write('-----------------\n')
    fid.write('#1 written using the HAWC2 OpenMDAO wrapper\n')
    fid.write('Case ID: %s\n' % case_id)

    for i, st in enumerate(body.beam_structure):
        fid.write(header_full)
        fid.write('$%i %i\n' % (i + 1, st.s.shape[0]))
        data = np.array([getattr(st, var) for var in st.var]).T
        np.savetxt(fid, data, fmt='%'+' %i.%ie' % (col_width, exp_prec))

    fid.close()


class HAWC2InputWriter(object):
    """
    Class to write HAWC2 input files.

    parameters
    ----------
    case_id: str
        Name of the file to write.
    vartrees: HAWC2VarTrees
        Variable tree containing all the model infomation.
    data_directory: str
        Name of data directory.
    res_directory: str
        Name of results directory.
    turb_directory: str
        Name of turbulence files directory.
    log_directory: str
        Name of log files directory.
    control_directory: str
        Name of controller files directory.
    returns
    -------
        nothing
    """
    def __init__(self, **kwargs):

        self.case_id = 'hawc2_case'
        self.vartrees = HAWC2VarTrees()

        self.data_directory = 'data'
        self.res_directory = 'res'
        self.turb_directory = 'turb'
        self.log_directory = 'logfile'
        self.control_directory = 'control'

        for k, w in kwargs.iteritems():
            try:
                setattr(self, k, w)
            except:
                pass

        # if not os.path.exists(self.data_directory):
        #     os.mkdir(self.data_directory)

    def execute(self):

        if not os.path.exists(self.data_directory):
            os.mkdir(self.data_directory)
        self.case_idout = self.case_id
        self._write_all()
        self._write_master()
        self._write_pcfile()
        self._write_aefile()

    def _write_all(self):

        self.htc_master = []

        self._write_simulation()
        self._write_wind()
        self._write_aero()
        self._write_aerodrag()
        self._write_structure()
        self._write_dlls()
        self._write_output()

    def _write_master(self):
        """
        write htc_master list to file
        """
        fid = open(self.case_id + '.htc', 'w')

        for i, line in enumerate(self.htc_master):
                line = line.rstrip()+'\n'
                self.htc_master[i] = line
                fid.write(line)
        fid.close()

    def _write_simulation(self):
        """ write simulation block """

        sim = []
        sim.append('begin simulation')
        sim.append('time_stop    %s' % self.vartrees.sim.time_stop)
        sim.append('solvertype   %i' %
                   self.vartrees.sim.solvertype)
        sim.append('convergence_limits %9.2e %9.2e %9.2e' %
                   tuple(self.vartrees.sim.convergence_limits))
        sim.append('on_no_convergence %s'% self.vartrees.sim.on_no_convergence)
        sim.append('max_iterations %i' % self.vartrees.sim.max_iterations)
        sim.append('logfile %s' % self.vartrees.sim.logfile)
        sim.append('begin newmark')
        sim.append('  deltat    %1.3f' % self.vartrees.sim.newmark_deltat)
        sim.append('end newmark')
        sim.append('end simulation')
        # Adding indent and semicolons
        sim[0] += ';'
        sim[-1] += ';'
        for iw, w in enumerate(sim[1:-1]):
            sim[iw+1] = '  ' + w + ';'
        self.htc_master.extend(sim)

    def _write_wind(self):
        """ write wind definition block """

        wind_vt = self.vartrees.wind
        fmt = ' %12.6e'
        wind = []
        wind.append('begin wind')
        wind.append('density' + fmt % wind_vt.density)
        wind.append('wsp' + fmt % wind_vt.wsp)
        wind.append('tint' + fmt % wind_vt.tint)
        wind.append('horizontal_input  %i' % wind_vt.horizontal_input)
        wind.append('windfield_rotations' + 3*fmt %
                    tuple(wind_vt.windfield_rotations))
        wind.append('center_pos0' + 3*fmt % tuple(wind_vt.center_pos0))
        wind.append('shear_format         %i' % wind_vt.shear_type +
                    fmt % wind_vt.shear_factor)
        wind.append('turb_format          %i' % wind_vt.turb_format)
        wind.append('tower_shadow_method  %i' % wind_vt.tower_shadow_method)

        for wind_ramp in wind_vt.wind_ramp_abs:
            wind.append('wind_ramp_abs' + 4*fmt % tuple(wind_ramp))

        if wind_vt.scale_time_start > 0:
            wind.append('scale_time_start' + fmt % wind_vt.scale_time_start)

        if wind_vt.wind_ramp_t1 > 0 or wind_vt.wind_ramp_t0 > 0 :
            wind.append('wind_ramp_factor' + 4*fmt %
                        (wind_vt.wind_ramp_t0, wind_vt.wind_ramp_t1,
                         wind_vt.wind_ramp_factor0, wind_vt.wind_ramp_factor1))

        if wind_vt.iec_gust:
            wind.append('iec_gust %s' % wind_vt.iec_gust_type + 4*fmt %
                        (wind_vt.G_A, wind_vt.G_phi0, wind_vt.G_t0,
                         wind_vt.G_T))

        if wind_vt.turb_format == 1:
            wind.extend(self._write_mann_turbulence())

        if wind_vt.tower_shadow_method > 0:
            wind.extend(self._write_tower_potential())

        wind.append('end wind')
        # Adding indent and semicolons
        wind[0] += ';'
        wind[-1] += ';'
        for iw, w in enumerate(wind[1:-1]):
            wind[iw+1] = '  ' + w + ';'
        self.htc_master.extend(wind)

    def _write_mann_turbulence(self):
        """ write mann turbulence model"""

        fmt = ' %12.6e'
        mann_vt = self.vartrees.wind.mann
        mann = []
        mann.append('begin mann')
        if mann_vt.create_turb:
            mann.append(('create_turb_parameters' + 3*fmt + ' %i' + fmt) %
                        (mann_vt.L, mann_vt.alfaeps, mann_vt.gamma,
                         mann_vt.seed, mann_vt.highfrq_compensation))
        mann.append('filename_u %s' % (os.path.join(self.turb_directory,
                                       mann_vt.turb_base_name + 'u.bin')))
        mann.append('filename_v %s' % (os.path.join(self.turb_directory,
                                       mann_vt.turb_base_name + 'v.bin')))
        mann.append('filename_w %s' % (os.path.join(self.turb_directory,
                                       mann_vt.turb_base_name + 'w.bin')))
        mann.append(('box_dim_u %5i'+fmt) % (mann_vt.box_nu, mann_vt.box_du))
        mann.append(('box_dim_v %5i'+fmt) % (mann_vt.box_nv, mann_vt.box_dv))
        mann.append(('box_dim_w %5i'+fmt) % (mann_vt.box_nw, mann_vt.box_dw))
        mann.append('std_scaling'+3*fmt % tuple(mann_vt.std_scaling))
        mann.append('end mann')

        for im, m in enumerate(mann[1:-1]):
            mann[im+1] = '  ' + m
        return mann

    def _write_tower_potential(self):
        """ write tower shadow with potential method"""

        fmt = ' %12.6e'
        tower_pot = []
        if hasattr(self.vartrees.wind, 'tower_potential'):
            tp = self.vartrees.wind.tower_potential
            tower_pot.append('begin tower_shadow_potential_2')
            tower_pot.append('tower_mbdy_link %s' % tp.tower_mbdy_link)
            tower_pot.append('nsec %d' % tp.nsec)
            for sec in tp.sections:
                tower_pot.append('radius' + 2*fmt % tuple(sec))
            tower_pot.append('end tower_shadow_potential_2')

        for im, m in enumerate(tower_pot[1:-1]):
            tower_pot[im+1] = '  ' + m
        return tower_pot

    def _write_aero(self):
        """ write aerodynamic block """

        aerovt = self.vartrees.aero
        aero = []
        aero.append('begin aero')
        aero.append('nblades  %d' % aerovt.nblades)
        aero.append('hub_vec shaft %i' % aerovt.hub_vec_coo)
        for link in aerovt.links:
            aero.append('link %i mbdy_c2_def %s' % (link[0], link[2]))
        aero.append('ae_filename %s' % aerovt.ae_filename)
        aero.append('pc_filename %s' % aerovt.pc_filename)
        aero.append('induction_method %i' % aerovt.induction_method)
        aero.append('aerocalc_method  %i' % aerovt.aerocalc_method)
        if aerovt.aero_distribution_file != '':
            aero.append('aero_distribution %s %i' %
                        (aerovt.aero_distribution_file,
                         aerovt.aero_distribution_set))
        else:
            aero.append('aerosections %i' % aerovt.aerosections)
        aero.append('ae_sets         %s' % ' '.join(map(str, aerovt.ae_sets)))
        aero.append('tiploss_method  %i' % aerovt.tiploss_method)
        aero.append('dynstall_method %i' % aerovt.dynstall_method)
        aero.append('end aero')

        aero = [a + ';' for a in aero]
        for ia, a in enumerate(aero[1:-1]):
            aero[ia+1] = '  ' + a

        self.htc_master.extend(aero)

    def _write_aerodrag(self):
        """ write aerodrag block """

        aerodrag = []
        fmt = ' %12.6e'
        if len(self.vartrees.aerodrag.elements) > 0:
            aerodrag.append('begin aerodrag')
            for i, e in enumerate(self.vartrees.aerodrag.elements):
                aerodrag.append('begin aerodrag_element')
                aerodrag.append('  mbdy_name %s' % e.mbdy_name)
                aerodrag.append('  aerodrag_sections %s %i' % (
                                     e.dist, e.calculation_points))
                aerodrag.append('  nsec %d' % len(e.sections))
                for j in range(len(e.sections)):
                    aerodrag.append('  sec' + 3*fmt % tuple(e.sections[j][:]))
                aerodrag.append('end aerodrag_element')
            aerodrag.append('end aerodrag')
            aerodrag = [a + ';' for a in aerodrag]
            for ia, a in enumerate(aerodrag[1:-1]):
                aerodrag[ia+1] = '  ' + a
            self.htc_master.extend(aerodrag)

    def _write_structure_res(self):
        """ write result files for eigen analysis """

        structure_out = []
        if self.vartrees.sim.eig_out:
            structure_out.append('beam_output_file_name ./info/%s_beam.dat;' %
                                 self.case_id)
            structure_out.append('body_output_file_name ./info/%s_body.dat;' %
                                 self.case_id)
            structure_out.append('struct_inertia_output_file_name \
                                 ./info/%s_inertia.dat;' % self.case_id)
            structure_out.append('body_eigenanalysis_file_name \
                                 ./info/%s_body_eigs.dat;' % self.case_id)
            structure_out.append('constraint_output_file_name \
                                 ./info/%s_constraints.dat;' % self.case_id)
            structure_out.append('structure_eigenanalysis_file_name \
                                 ./info/%s_struct_eigs.dat 1;' % self.case_id)

        self.structure.extend(structure_out)

    def _write_structure(self):
        """ write model structure"""

        self.structure = []
        self.structure.append('begin new_htc_structure;')
        self._write_structure_res()
        self._write_main_bodies()
        self._write_orientations()
        self._write_constraints()
        self.structure.append('end new_htc_structure;')
        self.htc_master.extend(self.structure)

    def _write_main_bodies(self):
        """ write all main bodies """

        for name in self.vartrees.body_order:
            self.structure.extend(self._write_main_body(name))

    def _write_main_body(self, body_name):
        """ write one main body """

        fmt = ' %23.15e'
        body = self.vartrees.main_bodies.get_main_body(body_name)
        main_body = []
        main_body.append('begin main_body')
        if body.copy_main_body is not '':
            main_body.append('name           %s' % body.body_name)
            main_body.append('copy_main_body %s' % body.copy_main_body)
        else:
            main_body.append('name        %s' % body.body_name)
            main_body.append('type        %s' % body.body_type)
            if body.nbodies < body.c12axis.shape[0]:
                main_body.append('nbodies     %d' % body.nbodies)
            else:
                main_body.append('nbodies     %d' % (body.c12axis.shape[0]-1))
            main_body.append('node_distribution %s'% body.node_distribution)

            if body.damping_type is 'ani':
                main_body.append('damping_aniso' + 6*fmt %
                                 tuple(body.damping_aniso))
            else:
                main_body.append('damping_posdef' + 6*fmt %
                                 tuple(body.damping_posdef))

            for i in range(len(body.concentrated_mass)):
                main_body.append('concentrated_mass' + 8*fmt %
                                 tuple(body.concentrated_mass[i]))
            main_body.append('begin timoschenko_input')
            tmpname = ''.join([i for i in body.body_name if not i.isdigit()])
            body.st_filename = (os.path.join(self.data_directory, self.case_id +
                                           '_' + tmpname + '_st.dat'))
            main_body.append('  filename %s' % body.st_filename)
            if body.st_input_type is not 0:
                main_body.append('  FPM %d' % body.st_input_type)
            main_body.append('  set %d %d' % tuple(body.body_set))
            main_body.append('end timoschenko_input')
            main_body.append('begin c2_def')
            main_body.append('  nsec %i' % body.c12axis.shape[0])
            for i in range(body.c12axis.shape[0]):
                main_body.append('  sec %2i' % (i+1) + 4*' %22.15e' %
                                 tuple(body.c12axis[i, :]))
            main_body.append('end c2_def')
            if len(body.beam_structure) > 0:
                self._write_stfile(body)

        main_body.append('end main_body')

        # add indent and semicolon
        main_body[0] = '  ' + main_body[0] + ';'
        main_body[-1] = '  ' + main_body[-1] + ';'
        for ib, b in enumerate(main_body[1:-1]):
            main_body[ib+1] = '    ' + b + ';'

        return main_body

    def _write_orientations(self):

        orientations = []
        orientations.append('begin orientation')
        fmt = ' %22.15e'
        for name in self.vartrees.body_order:
            body = self.vartrees.main_bodies.get_main_body(name)
            for orien in body.orientations:
                if orien.type == 'base':
                    orientations.append('begin base')
                    orientations.append('  body %s' % body.body_name)
                    orientations.append('  inipos' + 3*fmt %
                                        tuple(orien.inipos))
                    for eulerang in orien.body_eulerang:
                        orientations.append('  body_eulerang'+3*fmt %
                                            tuple(eulerang))
                    orientations.append('end base')

                elif orien.type == 'relative':
                    orientations.append('begin relative')
                    fmt2 = '  body1 '+_get_fmt(orien.body1)
                    orientations.append(fmt2 % tuple(orien.body1))
                    fmt2 = '  body2 '+_get_fmt(orien.body2)
                    orientations.append(fmt2 % tuple(orien.body2))
                    for eulerang in orien.body2_eulerang:
                        orientations.append('  body2_eulerang'+3*fmt %
                                            tuple(eulerang))
                    if orien.mbdy2_ini_rotvec_d1[3] != 0.:
                        orientations.append('  body2_ini_rotvec_d1'+4*fmt %
                                            tuple(orien.mbdy2_ini_rotvec_d1))
                    orientations.append('end relative')
        orientations.append('end orientation')
        # add indent and semicolon
        orientations = ['  ' + o + ';' for o in orientations]
        for io, o in enumerate(orientations[1:-1]):
            orientations[io+1] = '  ' + o
        self.structure.extend(orientations)

    def _write_constraints(self):

        constraints = []
        fmt = ' %22.15e'
        constraints.append('begin constraint')
        for name in self.vartrees.body_order:
            body = self.vartrees.main_bodies.get_main_body(name)
            for con in body.constraints:
                if con.con_type in ['fix0', 'fix2', 'fix3']:
                    constraints.append('begin %s' % con.con_type)
                    constraints.append('  body %s' % con.mbdy)
                    if con.disable_at > 0.:
                        constraints.append('  disable_at %s' %
                                           con.disable_at)
                    if con.con_type in ['fix2', 'fix3']:
                        constraints.append('  dof %i %i %i' % tuple(con.dof))
                    constraints.append('end %s' % con.con_type)

                elif con.con_type in ['fix1', 'fix4']:
                    constraints.append('begin %s' % con.con_type)
                    fmt2 = '  body1 ' + _get_fmt(con.mbdy1)
                    constraints.append(fmt2 % tuple(con.mbdy1))
                    fmt2 = '  body2 ' + _get_fmt(con.mbdy2)
                    constraints.append(fmt2 % tuple(con.mbdy2))
                    constraints.append('end %s' % con.con_type)

                elif 'bearing' in con.con_type:
                    constraints.append('begin %s' % con.con_type)
                    constraints.append('  name %s' % con.name)
                    fmt2 = '  body1 '+_get_fmt(con.mbdy1)
                    constraints.append(fmt2 % tuple(con.mbdy1))
                    fmt2 = '  body2 '+_get_fmt(con.mbdy2)
                    constraints.append(fmt2 % tuple(con.mbdy2))
                    constraints.append(('  bearing_vector %i'+3*fmt) %
                                       tuple(con.bearing_vector))
                    if con.con_type == 'bearing3':
                        constraints.append('  omegas'+3*fmt % con.omegas)
                    else:
                        if con.disable_at > 0:
                            constraints.append('  disable_at %s' %
                                               con.disable_at)
                    constraints.append('end %s' % con.con_type)

        constraints.append('end constraint')
        # add indent and semicolon
        constraints = ['  ' + c + ';' for c in constraints]
        for ic, c in enumerate(constraints[1:-1]):
            constraints[ic+1] = '  ' + c
        self.structure.extend(constraints)

    def _write_dlls(self):
        """
        write all the dlls
        """
        dlls = []
        dlls.append('begin dll;')
        for name in self.vartrees.dlls_order:
            dlls.extend(self._write_dll(name))
        dlls.append('end dll;')
        # add indent
        for i, c in enumerate(dlls[1:-1]):
            dlls[i+1] = '  ' + c
        self.htc_master.extend(dlls)

    def _write_dll(self, dll_name):
        """ write general type2 dll"""

        fmt = '%2i %22.15e'
        dll = []
        dll_vt = getattr(self.vartrees.dlls, dll_name)
        dll_init = dll_vt. dll_init
        dll.append('begin type2_dll')
        dll.append('name %s' % dll_vt.name)
        dll.append('filename %s' % dll_vt.filename)
        dll.append('dll_subroutine_init %s' % dll_vt.dll_subroutine_init)
        dll.append('dll_subroutine_update %s' % dll_vt.dll_subroutine_update)
        dll.append('arraysizes_init %i %i' % tuple(dll_vt.arraysizes_init))
        dll.append('arraysizes_update %i %i' % tuple(dll_vt.arraysizes_update))
        dll.append('begin init')

        for i in range(len(dll_init.init_dic.keys())):
            val = getattr(dll_init,
                          dll_init.init_dic[i+1][0])/dll_init.init_dic[i+1][1]
            dll.append('  constant ' + fmt % (i+1, val))

        dll.append('end init')

        dll.append('begin output')
        for i in range(len(dll_vt.output.out_dic.keys())):
            val = getattr(dll_vt.output, 'out_%i' % (i + 1))
            dll.append('  %s' % val._print())
        dll.append('end output')

        if len(dll_vt.actions.action_dic.keys()) > 0:
            dll.append('begin actions')
            for i in range(len(dll_vt.actions.action_dic.keys())):
                val = getattr(dll_vt.actions, 'action_%i' % (i + 1))
                dll.append('  %s' % val._print())
            dll.append('end actions')

        dll.append('end type2_dll')

        # add indent and semicolon
        dll[0] += ';'
        dll[-1] += ';'
        for i, d in enumerate(dll[1:-1]):
            dll[i+1] = '  ' + d + ';'

        return dll

    def _write_output(self):

        sns = []
        sns.append('begin output')
        sns.append('  filename %s' % (os.path.join(self.vartrees.output.res_dir,
                                                   self.case_id)))
        sns.append('  time %3.6f %3.6f' % (self.vartrees.output.time_start,
                                           self.vartrees.sim.time_stop))
        sns.append('  data_format %s' % self.vartrees.output.out_format)
        sns.append('  buffer 1')

        for i in range(len(self.vartrees.output.sensor_list)):
            sns.append('  ' + self.vartrees.output.sensor_list[i])
        sns.append('end output')

        sns = [s+';' for s in sns]

        self.htc_master.extend(sns)

    def _write_aefile(self):

        write_aefile(self.vartrees.aero.ae_filename, self.vartrees.blade_ae)

    def _write_stfile(self, body):

        write_stfile(body, self.case_id)

    def _write_pcfile(self):

        write_pcfile(self.vartrees.aero.pc_filename, self.vartrees.airfoildata)


class HAWC2SInputWriter(HAWC2InputWriter):
    """
    HAWC2InputWriter-type class to write HAWC2s files.

    parameters
    ----------
        same as for HAWC2InputWriter

    returns
    -------
        nothing
    """
    def __init__(self, **kwargs):
        super(HAWC2SInputWriter, self).__init__()

        for k, w in kwargs.iteritems():
            try:
                setattr(self, k, w)
            except:
                pass

    def execute(self):

        if not os.path.exists(self.data_directory):
            os.mkdir('data')
        self.case_idout = self.case_id

        self._write_all()
        self._write_master()
        self._write_pcfile()
        self._write_aefile()

    def _write_all(self):

        self.htc_master = []

        self._write_aero()
        self._write_wind()
        self._write_structure()
        self._write_hawcstab2()

    def _write_hawcstab2(self):

        self.h2s = []
        self.h2s.append('begin hawcstab2;')
        self._write_hawcstab2_structure()
        self._write_h2s_operational_data()
        self._write_h2s_control()
        self._write_h2s_commands()
        self.h2s.append('end hawcstab2;')
        self.htc_master.extend(self.h2s)

        self._write_operational_data_file()

    def _write_h2s_commands(self):

        opt_vt = self.vartrees.h2s.options
        cmd = []
        for name in self.vartrees.h2s.commands:
            if name == 'compute_optimal_pitch_angle':
                cmd.append('compute_optimal_pitch_angle use_operational_data')

            elif name == 'compute_steady_states':
                cmd.append('compute_steady_states %s %s %s %s' %
                           (opt_vt.bladedeform, opt_vt.tipcorrect,
                            opt_vt.induction, opt_vt.gradients))

            elif name == 'compute_steadystate':
                cmd.append('compute_steadystate %s %s %s %s' %
                           (opt_vt.bladedeform, opt_vt.tipcorrect,
                            opt_vt.induction, opt_vt.gradients))

            elif name == 'compute_stability_analysis':
                cmd.append(('compute_stability_analysis %s %s %i'+4*' %12.6e' +
                           ' %s') % (opt_vt.matrixwriteout,
                                     opt_vt.eigenvaluewriteout,
                                     opt_vt.number_of_modes,
                                     opt_vt.maximum_damping,
                                     opt_vt.minimum_frequency,
                                     opt_vt.zero_pole_threshold,
                                     opt_vt.aero_deflect_ratio,
                                     opt_vt.frequencysorting))

            elif name == 'compute_aeroservoelastic':
                cmd.append(('compute_aeroservoelastic %s %s %i'+4*'%12.6e'+'%s'
                            ) % (opt_vt.matrixwriteout,
                                 opt_vt.eigenvaluewriteout,
                                 opt_vt.number_of_modes,
                                 opt_vt.maximum_damping,
                                 opt_vt.minimum_frequency,
                                 opt_vt.zero_pole_threshold,
                                 opt_vt.aero_deflect_ratio,
                                 opt_vt.frequencysorting))

            elif name == 'save_cl_matrices_all':
                if opt_vt.vloc_out:
                    cmd.append('save_cl_matrices_all vloc_out')
                else:
                    cmd.append('save_cl_matrices_all')
            elif name == 'compute_structural_modal_analysis':
                if opt_vt.blade_only:
                    cmd.append('compute_structural_modal_analysis '
                               'bladeonly %i' % opt_vt.number_of_modes)
                else:
                    cmd.append('compute_structural_modal_analysis '
                               'nobladeonly %i' % opt_vt.number_of_modes)

            elif name == 'basic_dtu_we_controller':
                init = self.vartrees.dlls.risoe_controller.dll_init

                cmd.append(('basic_dtu_we_controller' + 10*' %20.15e' + ' %i' +
                            3*' %20.15e') % (init.pgTorque, init.igTorque,
                                             init.Qg, init.pgPitch,
                                             init.igPitch, init.KK1, init.KK2,
                                             init.generatorFreq,
                                             init.generatorDamping,
                                             init.ffFreq, init.generatorSwitch,
                                             init.Kp2, init.Ko1, init.Ko2))
            else:
                cmd.append(name)
        cmd = ['  ' + c + ';' for c in cmd]
        self.h2s.extend(cmd)

    def _write_operational_data_file(self):

        h2s = self.vartrees.h2s

        # Write opt file only if it is not computed by HS2
        wsp = []
        rpm = []
        pitch = []
        if 'compute_optimal_pitch_angle' not in h2s.commands:
            if len(h2s.wsp_cases) > 0:
                for w in h2s.wsp_cases:
                    rpm.append(np.interp(w, h2s.wsp_curve, h2s.rpm_curve))
                    pitch.append(np.interp(w, h2s.wsp_curve, h2s.pitch_curve))
                    wsp.append(w)
        for case in h2s.cases:
            try:
                wsp.append(case['wsp'])
                pitch.append(case['pitch'])
                rpm.append(case['rpm'])
            except:
                pass
        if len(wsp) > 0:
            data = np.array([wsp, pitch, rpm]).T
        else:
            data = np.array([h2s.wsp_curve, h2s.pitch_curve, h2s.rpm_curve]).T
        if len(data.shape) == 1:
            data = np.array([data])
        fid = open(self.case_id + '.opt', 'w')
        fid.write(('%i Wind speed [m/s]          Pitch [deg]     ' +
                  'Rot. speed [rpm]\n') % data.shape[0])
        np.savetxt(fid, data)
        fid.close()

    def _write_hawcstab2_structure(self):

        h2s_vt = self.vartrees.h2s
        st = []
        st.append('begin ground_fixed_substructure')
        for name in h2s_vt.ground_fixed.main_body:
            st.append('  main_body %s' % name)
        if h2s_vt.ground_fixed.log_decrements[0] != 0:
            st.append('  log_decrements %12.6f %12.6f' %
                      tuple(h2s_vt.ground_fixed.log_decrements))
        st.append('end ground_fixed_substructure')

        st.append('begin rotating_axissym_substructure')
        for name in h2s_vt.rotating_axissym.main_body:
            st.append('  main_body %s' % name)
        if h2s_vt.rotating_axissym.log_decrements[0] != 0:
            st.append('  log_decrements %12.6f %12.6f' %
                      tuple(h2s_vt.rotating_axissym.log_decrements))
        st.append('end rotating_axissym_substructure')

        st.append('begin rotating_threebladed_substructure')
        for name in h2s_vt.rotating_threebladed.main_body:
            st.append('  main_body %s' % name)
        if h2s_vt.rotating_threebladed.log_decrements[0] != 0:
            st.append('  log_decrements' + 6*' %12.6f' %
                      tuple(h2s_vt.rotating_threebladed.log_decrements))
        st.append('  second_order_actuator %s %12.6e %12.6e' %
                  (h2s_vt.second_order_actuator.name,
                   h2s_vt.second_order_actuator.frequency,
                   h2s_vt.second_order_actuator.damping))
        st.append('end rotating_threebladed_substructure')

        st = ['  ' + s + ';' for s in st]
        self.h2s.extend(st)

    def _write_h2s_control(self):

        ctr = []
        dll_init = self.vartrees.dlls.risoe_controller.dll_init
        ctr.append('begin controller_tuning')
        ctr.append('  partial_load %3.6f %3.6f' % (dll_init.poleFreqTorque,
                                                   dll_init.poleDampTorque))
        ctr.append('  full_load %3.6f %3.6f' % (dll_init.poleFreqPitch,
                                                dll_init.poleDampPitch))
        ctr.append('  gain_scheduling %d' % dll_init.gainScheduling)
        ctr.append('  constant_power %d' % dll_init.generatorSwitch)

        if len(self.vartrees.h2s.options.regions) > 1:
            ctr.append('  regions %i %i %i %i' %
                       tuple(self.vartrees.h2s.options.regions))
        ctr.append('end controller_tuning')

        ctr.append('begin controller')
        ctr.append('  begin input')
        for i in range(len(self.vartrees.h2s.ch_list_in.sensor_list)):
            ctr.append('    '+self.vartrees.h2s.ch_list_in.sensor_list[i])
        ctr.append('  end input')
        ctr.append('  begin output')
        for i in range(len(self.vartrees.h2s.ch_list_out.sensor_list)):
            ctr.append('    '+self.vartrees.h2s.ch_list_out.sensor_list[i])
        ctr.append('  end output')

        ctr.append('end controller')

        ctr = ['  ' + c + ';' for c in ctr]
        self.h2s.extend(ctr)

    def _write_h2s_operational_data(self):

        dll_init = self.vartrees.dlls.risoe_controller.dll_init
        opt = []
        opt.append('operational_data_filename %s' % self.case_id + '.opt')

        opt.append('begin operational_data')
        opt.append('  windspeed %12.6e %12.6e %d' % (dll_init.Vin, dll_init.Vout,
                                                     dll_init.nV))
        opt.append('  genspeed %12.6e %12.6e' % (dll_init.minRPM,
                                                 dll_init.maxRPM))
        opt.append('  gearratio %12.6e' % dll_init.gearRatio)
        opt.append('  minpitch %12.6e' % dll_init.minPitch)
        opt.append('  opt_lambda %22.15e' % dll_init.designTSR)
        opt.append('  maxpow %12.6e' % dll_init.ratedAeroPower)
        opt.append('  prvs_turbine %d' % dll_init.prvs_turbine)
        opt.append('  include_torsiondeform %d' %
                   self.vartrees.h2s.options.include_torsiondeform)
        opt.append('  set_torque_limit %d' %
                   self.vartrees.h2s.options.set_torque_limit)
        opt.append('end operational_data')

        opt = ['  ' + o + ';' for o in opt]
        self.h2s.extend(opt)
