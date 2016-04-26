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

from wetb.prepost import dlcdefs
from fatfreq import fatfreq

from hawc2_wrapper.hawc2_aeroelasticsolver import HAWC2AeroElasticSolver
from hawc2_wrapper.hawc2s_aeroelasticsolver import HAWC2SAeroElasticSolver


def compute_wind_input():

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
    cf['hawc2bin'] = 'HAWC2mb.exe'
    cf['verbose'] = False
    config['HAWC2Wrapper'] = cf

    # HAWC2Outputs
    cf = {}
    cf['channels'] = ['bearing-shaft_rot-angle-deg']
    config['HAWC2Outputs'] = cf

    # HAWC2InputWriter
    cf = {}
    config['HAWC2InputWriter'] = cf

    # General options
    config['master_file'] = 'main_h2_wind.htc'

    config['with_tsr'] = False
    config['with_structure'] = False
    config['with_geom'] = False
    config['aerodynamic_sections'] = 30
    config['structural_sections'] = 20

    root.add('loads', HAWC2AeroElasticSolver(config, './DLCs_fatigue/', None, None), promotes=['*'])

    top.setup()
    top.run()


###############################################################################        
# Not required if wind_structure already exists in the directory
###############################################################################
#compute_wind_input()

cases_list = [['hawc2_model_dlc12_wsp10_wdir000_s1001_24036903/res/dlc12_iec61400-1ed3/dlc12_wsp10_wdir000_s1001'],
              ['hawc2_model_dlc12_wsp15_wdir000_s1002_29122925/res/dlc12_iec61400-1ed3/dlc12_wsp15_wdir000_s1002'],
              ['hawc2_model_dlc12_wsp20_wdir000_s1003_42782675/res/dlc12_iec61400-1ed3/dlc12_wsp20_wdir000_s1003']]

#fatfreq.ComputeWindResp(cases_list, 'wind_structure', Fmax=2.5)  
###############################################################################

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
cf['hawc2bin'] = 'HAWC2s.exe'

cf['verbose'] = False
config['HAWC2Wrapper'] = cf

# HAWC2Outputs
cf = {}
cf['blade'] = ['aoa', 'Ft', 'Fn', 'cl', 'cd',  'cp', 'ct',
               'disp_x', 'disp_y', 'disp_z', 'disp_rot_z']
cf['rotor'] = ['wsp', 'pitch', 'rpm', 'P', 'Q', 'T', 'CP',  'CT', 'Mx']
config['HAWC2SOutputs'] = cf

# HAWC2InputWriter
cf = {}
config['HAWC2InputWriter'] = cf

# General options
config['master_file'] = 'main_hs2.htc'

config['with_tsr'] = False
config['with_structure'] = False
config['with_geom'] = False
config['aerodynamic_sections'] = 30
config['structural_sections'] = 20

# Cases
cf2 = {}

cf = {}
cf['wsp'] = [10., 15., 20.]
cf['pitch'] = [5., 11., 15.]
cf['rpm'] = [13., 13., 13.]
cf['user'] = []
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

cf = {}
cf['Sensors'] = [[0], [1]]
cf['m'] = [3, 3]

config['Fatigue'] = cf

root.add('loads', HAWC2SAeroElasticSolver(config), promotes=['*'])

top.setup()

top.run()