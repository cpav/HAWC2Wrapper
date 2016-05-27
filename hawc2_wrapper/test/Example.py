# -*- coding: utf-8 -*-
"""
Created on Thu Apr 07 13:28:26 2016

@author: tlbl
"""
import numpy as np
import time

from openmdao.core.mpi_wrap import MPI
from openmdao.api import Problem, Group, Component, ScipyOptimizer, IndepVarComp, ParallelFDGroup, ExecComp
from openmdao.solvers.ln_gauss_seidel import LinearGaussSeidel
from hawc2_wrapper.hawc2_inputdict import read_hawc2_st_file
from hawc2_wrapper.hawc2s_aeroelasticsolver import HAWC2SAeroElasticSolver

if MPI:
    from openmdao.core.petsc_impl import PetscImpl as impl
else:
    from openmdao.core.basic_impl import BasicImpl as impl


class NormalizeDesVar(Component):

    def __init__(self, init_blade_beam_structure):
        super(NormalizeDesVar, self).__init__()

        nsec_st = init_blade_beam_structure.shape[0]

        self.init_blade_beam_structure = init_blade_beam_structure.copy()
        self.add_param('Ixnorm', val=np.zeros(nsec_st))
        self.add_param('Iynorm', val=np.zeros(nsec_st))
        self.add_output('blade_beam_structure', shape=(nsec_st, 19))

    def solve_nonlinear(self, params, unknowns, resids):

        blade_beam_structure = self.init_blade_beam_structure.copy()
        print 'Ixnorm', params['Ixnorm']
        print 'Iynorm', params['Iynorm']
        blade_beam_structure[:, 10] *= (1.+params['Ixnorm'])
        blade_beam_structure[:, 11] *= (1.+params['Iynorm'])

        unknowns['blade_beam_structure'] = blade_beam_structure


class CostFunction(Component):

    def __init__(self, nsec_st):
        super(CostFunction, self).__init__()
        self.add_param('disp_y', np.zeros((1, nsec_st)))
        self.add_output('cf', shape=(1))

    def solve_nonlinear(self, params, unknowns, resids):

        ref = -4.7
        print 'Current blade tip position: %1.4f m and target: %1.4fm ' % (params['disp_y'][-1][-1], ref)
        unknowns['cf'] = np.sqrt((params['disp_y'][-1][-1]/ref - 1.)**2)


class CostFunction2(Component):

    def __init__(self, nmodes):
        super(CostFunction2, self).__init__()
        self.add_param('freq_factor', np.zeros((1, 2*nmodes)))
        self.add_output('cf', shape=(1))

    def solve_nonlinear(self, params, unknowns, resids):

        print 'Current freq_factor: %1.4f ' % (params['freq_factor'][0,0])
        unknowns['cf'] = params['freq_factor'][0,0]


def setup_config_dict(n):

    config = {}
    config['master_file'] = './main_hs2_example.htc'
    config['with_structure'] = True
    config['with_geom'] = False
    config['structural_sections'] = n
    config['aerodynamic_sections'] = 20
    cf = {}
    config['HAWC2SInputWriter'] = cf
    cf = {}
    cf['blade_ni_span'] = 10
    cf['interp_from_htc'] = False
    cf['hub_radius'] = 3.
    config['HAWC2GeometryBuilder'] = cf
    cf = {}
    cf['dry_run'] = False
    cf['copyback_results'] = True
    cf['hawc2bin'] = 'HAWC2S.exe'
    cf['verbose'] = False
    config['HAWC2Wrapper'] = cf

    cf2 = {}
    cf2['wsp'] = [15.]
    cf2['pitch'] = [10.]
    cf2['rpm'] = [13.]

    cf = {}
    cf['wsp'] = []
    cf['pitch'] = []
    cf['rpm'] = []
    cf['user'] = [cf2]
    config['cases'] = cf
    cf = {}
    cf['blade'] = ['aoa', 'Ft', 'Fn', 'cl', 'cd',  'cp', 'ct',
                   'disp_x', 'disp_y', 'disp_z', 'disp_rot_z']
    cf['rotor'] = ['wsp', 'pitch', 'rpm', 'P', 'Q', 'T', 'CP',  'CT', 'Mx']
    config['HAWC2SOutputs'] = cf

    return config


def setup_init_structure(n):

    blade_beam_structure = []
    var = ['s', 'dm', 'x_cg', 'y_cg', 'ri_x', 'ri_y', 'x_sh', 'y_sh', 'E', 'G',
           'I_x', 'I_y', 'K', 'k_x', 'k_y', 'A', 'pitch', 'x_e', 'y_e']
    tmp = read_hawc2_st_file('data/DTU_10MW_RWT_Blade_st.dat', var)

    for v in var:
        blade_beam_structure.append(tmp[0][v])
    blade_beam_structure = np.array(blade_beam_structure).T

    init_blade_beam_structure = np.zeros((n, blade_beam_structure.shape[1]))
    s = np.linspace(0, blade_beam_structure[-1, 0], n)
    for i in range(blade_beam_structure.shape[1]):
        init_blade_beam_structure[:, i] = np.interp(s, blade_beam_structure[:, 0], blade_beam_structure[:, i])
    return init_blade_beam_structure


def set_top(optimize_flag):
    par_fd = 5

    if optimize_flag:
        top = Problem(impl=impl, root=ParallelFDGroup(par_fd))
    else:
        top = Problem(impl=impl, root=Group())

    # Set the Driver
    top.driver = ScipyOptimizer()
    top.driver.options['optimizer'] = 'SLSQP'
    top.driver.options['tol'] = 1.0e-5
    # from openmdao.drivers.pyoptsparse_driver import pyOptSparseDriver
    # top.driver = pyOptSparseDriver()
    # top.driver.options['optimizer'] = 'IPOPT'
    # top.driver.opt_settings['mu_strategy'] = 'adaptive'
    # top.driver.opt_settings['linear_solver'] = 'ma27'
    # top.driver.opt_settings['max_iter'] = 40
    # top.driver.opt_settings['tol'] = 1.e-5
    # Some options and parameters
    top.root.fd_options['force_fd'] = True
    top.root.fd_options['step_size'] = 1.e-3
    return top


def example_blade_tip_target():
    optimize_flag = True

    top = set_top(optimize_flag)

    n = 5  # number of section in the st file
    init_blade_beam_structure = setup_init_structure(n)
    config = setup_config_dict(n)

    # Add input variables
    top.root.add('Ixnorm_c', IndepVarComp('Ixnorm', np.zeros(n)), promotes=['*'])
    top.root.add('Iynorm_c', IndepVarComp('Iynorm', np.zeros(n)), promotes=['*'])

    # Add the Components and Groups forming the workflow
    top.root.add('normalize',
                 NormalizeDesVar(init_blade_beam_structure), promotes=['*'])
    top.root.add('loads',
                 HAWC2SAeroElasticSolver(config), promotes=['*'])
    top.root.add('cfuntion',
                 CostFunction(config['aerodynamic_sections']-2), promotes=['*'])

    # Add cost function
    top.driver.add_objective('cf')
    # Add design variable
    top.driver.add_desvar('Ixnorm', lower=-0.2, upper=0.2)
    return top


def example_frequency_placement():
    optimize_flag = True

    top = set_top(optimize_flag)

    n = 5  # number of section in the st file
    init_blade_beam_structure = setup_init_structure(n)
    config = setup_config_dict(n)
    cf = {}
    cf['mode_freq'] = np.array([[15., 0.71]])
    cf['mode_damp'] = np.array([[15., 2.7]])

    cf['mode_target_freq'] = np.array([[15., 0.75]])
    cf['mode_target_damp'] = np.array([[15., -1]])
    config['FreqDampTarget'] = cf

    # Add input variables
    top.root.add('Ixnorm_c', IndepVarComp('Ixnorm', np.zeros(n)), promotes=['*'])
    top.root.add('Iynorm_c', IndepVarComp('Iynorm', np.zeros(n)), promotes=['*'])

    # Add the Components and Groups forming the workflow
    top.root.add('normalize',
                 NormalizeDesVar(init_blade_beam_structure), promotes=['*'])
    top.root.add('loads',
                 HAWC2SAeroElasticSolver(config), promotes=['*'])
    top.root.add('cfuntion', CostFunction2(cf['mode_freq'].shape[1]-1), promotes=['*'])

    # Add cost function
    top.driver.add_objective('cf')
    # Add design variable
    top.driver.add_desvar('Iynorm', lower=-0.4, upper=0.4)
    return top

# SELECT HERE THE EXAMPLE TO RUN
# top = example_blade_tip_target()
top = example_frequency_placement()

# add the recorder
from openmdao.api import SqliteRecorder
recorder = SqliteRecorder('optimization.sqlite')

recorder.options['record_metadata'] = True
recorder.options['includes'] = ['Ixnorm', 'Iynorm', 'cf']

top.driver.add_recorder(recorder)

top.setup()

t0 = time.time()
top.run()
print('Total time', time.time()-t0)
