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

top = Problem(impl=impl, root=Group())
root = top.root
root.ln_solver = LinearGaussSeidel()

config = {}

# HAWC2GeometryBuilder
cf = {}
cf['interp_from_htc'] = False
config['HAWC2GeometryBuilder'] = cf

# HAWC2Wrapper
cf = {}
cf['dry_run'] = False
cf['copyback_results'] = False
cf['hawc2bin'] = 'HAWC2mb.exe'

cf['verbose'] = False
config['HAWC2Wrapper'] = cf

# HAWC2InputWriter
cf = {}
config['HAWC2InputWriter'] = cf

# General options
config['master_file'] = 'main_h2_envelope.htc'

config['with_tsr'] = False
config['with_structure'] = False
config['with_geom'] = False
config['aerodynamic_sections'] = 30
config['structural_sections'] = 20

# Cases
cf2 = {}
cf2['wsp']   = [8.]
cf2['pitch'] = [5.]
cf2['rpm']   = [13.]

cf = {}
cf['wsp'] = [10., 15., 20.]
cf['pitch'] = [5., 11., 15.]
cf['rpm'] = [13., 13., 13.]
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

# HAWC2Outputs
cf = {}

cf['channels'] = ['local-blade1-node-013-momentvec-z']
cf['ch_envelope'] = []
cf['nr_blades_out'] = 3

for nsec in range(config['structural_sections']):
    for nrb in range(cf['nr_blades_out']):
        cf['ch_envelope'].append(\
    ['local-blade'+'{:01d}'.format(nrb+1)+'-node-'+'{:03d}'.format(nsec+1)+'-momentvec-x',
     'local-blade'+'{:01d}'.format(nrb+1)+'-node-'+'{:03d}'.format(nsec+1)+'-momentvec-y',
     'local-blade'+'{:01d}'.format(nrb+1)+'-node-'+'{:03d}'.format(nsec+1)+'-momentvec-z',
     'local-blade'+'{:01d}'.format(nrb+1)+'-node-'+'{:03d}'.format(nsec+1)+'-forcevec-x',
     'local-blade'+'{:01d}'.format(nrb+1)+'-node-'+'{:03d}'.format(nsec+1)+'-forcevec-y',
     'local-blade'+'{:01d}'.format(nrb+1)+'-node-'+'{:03d}'.format(nsec+1)+'-forcevec-z'])
                    
                                                            
config['HAWC2Outputs'] = cf

root.add('loads', HAWC2AeroElasticSolver(config, './DLCs_longer/', None, None), promotes=['*'])


top.setup()

top.run()

#print(top.root.loads.unknowns['aggregate.envelope_sec1'])
