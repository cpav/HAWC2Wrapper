# -*- coding: utf-8 -*-
import unittest

from openmdao.core.mpi_wrap import MPI
from openmdao.api import Problem, Group
from openmdao.solvers.ln_gauss_seidel import LinearGaussSeidel

if MPI:
    from openmdao.core.petsc_impl import PetscImpl as impl
else:
    from openmdao.core.basic_impl import BasicImpl as impl


from hawc2_wrapper.hawc2_aeroelasticsolver import HAWC2AeroElasticSolver


top = Problem(impl=impl, root=Group())
root = top.root
root.ln_solver = LinearGaussSeidel()

config = {}

cf = {}
cf['blade_ni_span'] = 27
cf['interp_from_htc'] = False
cf['hub_radius'] = 3.
config['HAWC2GeometryBuilder'] = cf

cf = {}
cf['dry_run'] = False
cf['copyback_results'] = False
cf['hawc2bin'] = 'HAWC2mb.exe'
cf['verbose'] = False
config['HAWC2Wrapper'] = cf

cf = {}
cf['channels'] = ['local-blade1-node-013-momentvec-z']
config['HAWC2Outputs'] = cf

cf = {}
config['HAWC2InputWriter'] = cf

config['master_file'] = 'main_h2.htc'
config['with_tsr'] = False
config['with_structure'] = False
config['with_geom'] = False
config['aerodynamic_sections'] = 50
root.add('loads', HAWC2AeroElasticSolver(config, './DLCs_longer/', None, None), promotes=['*'])


top.setup()

top.run()

