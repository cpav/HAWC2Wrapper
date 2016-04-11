# -*- coding: utf-8 -*-
import unittest

from openmdao.core.mpi_wrap import MPI
from openmdao.api import Problem, Group
from openmdao.solvers.ln_gauss_seidel import LinearGaussSeidel
import numpy as np
if MPI:
    from openmdao.core.petsc_impl import PetscImpl as impl
else:
    from openmdao.core.basic_impl import BasicImpl as impl


from hawc2_wrapper.hawc2_aeroelasticsolver import HAWC2AeroElasticSolver
from hawc2_wrapper.hawc2s_aeroelasticsolver import HAWC2SAeroElasticSolver

H2 = False

top = Problem(impl=impl, root=Group())
root = top.root
root.ln_solver = LinearGaussSeidel()

config = {}

# HAWC2GeometryBuilder
cf = {}
cf['blade_ni_span'] = 27
cf['interp_from_htc'] = False
cf['hub_radius'] = 3.
config['HAWC2GeometryBuilder'] = cf

# HAWC2Wrapper
cf = {}
cf['dry_run'] = False
cf['copyback_results'] = False
if H2:
    cf['hawc2bin'] = 'HAWC2mb.exe'
else:
    cf['hawc2bin'] = 'HAWC2s.exe'

cf['verbose'] = False
config['HAWC2Wrapper'] = cf

# HAWC2Outputs
cf = {}
if H2:
    cf['channels'] = ['local-blade1-node-013-momentvec-z']
    config['HAWC2Outputs'] = cf
else:
    cf['blade'] = ['aoa', 'Ft', 'Fn', 'cl', 'cd',  'cp', 'ct',
                   'disp_x', 'disp_y', 'disp_z', 'disp_rot_z']
    cf['rotor'] = ['wsp', 'pitch', 'rpm', 'P', 'Q', 'T', 'CP',  'CT', 'Mx']
    config['HAWC2SOutputs'] = cf

# HAWC2InputWriter
cf = {}
config['HAWC2InputWriter'] = cf

# General options
if H2:
    config['master_file'] = 'main_h2.htc'
else:
    config['master_file'] = 'main_hs2.htc'

config['with_tsr'] = False
config['with_structure'] = False
config['with_geom'] = False
config['aerodynamic_sections'] = 30
config['structural_sections'] = 20

if not H2:
    # Cases
    cf2 = {}
    cf2['wsp']   = [8.]
    cf2['pitch'] = [5.]
    cf2['rpm']   = [13.]

    cf = {}
    cf['wsp'] = [10., 12., 14., 16., 18., 20.]
    cf['pitch'] = [5., 7., 9., 11., 13., 15.]
    cf['rpm'] = [13., 13., 13., 13., 13., 13.]
    cf['user'] = [cf2]
    config['cases'] = cf

    # Frequency placement
    cf = {}
    cf['mode_freq'] = np.array([[15., 0.25, 0.64],
                                [20., 0.25, 0.64]])
    cf['mode_damp'] = np.array([[15., 0.8, 80.],
                                [20., 0.8, 80.]])
    cf['mode_target_freq'] = np.array([[15., 0.3, -1],
                                       [20., 0.25, 0.6]])
    cf['mode_target_damp'] = np.array([[15., -1,  85.],
                                       [20., 0.4, 80.]])
    config['FreqDampTarget'] = cf


if H2:
    root.add('loads', HAWC2AeroElasticSolver(config, './DLCs_longer/', None, None), promotes=['*'])
else:
    root.add('loads', HAWC2SAeroElasticSolver(config), promotes=['*'])



top.setup()

top.run()

