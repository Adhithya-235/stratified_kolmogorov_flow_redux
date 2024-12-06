#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Usage: 
  galerkin_kolmogorov_flow.py [--Rb=<buoyancy_reynolds> --Pr=<prandtl> --Fr=<froude> \
  --Nz=<Nz> --Tend=<stop_time>] 
  
Options:
  --Rb=<buoyancy_reynolds>          Buoyancy Reynolds Number [default: 50.0]
  --Pr=<prandtl>                    Prandtl Number [default: 1.0]
  --Fr=<froude>                     Froude Number [default: 0.02]
  --Nz=<Nz>                         Number of vertical modes [default: 128]
  --Tend=<stop_time>                Simulation stop time [default: 10.0]  
"""

"""
Dedalus script for the numerical simulation of a galerkin-truncated version of the 
2D Boussinesq Equations to explore the dynamics of strongly stratified Kolmogorov flow.

This script uses a Fourier basis in z for the out-of-plane vorticity 'xi',
and buoyancy 'b'. 

This script can be run serially or in parallel, and uses the built-in analysis 
framework to save data in HDF5 files. The 'merge_procs' command can be used to 
merge distributed analysis sets from parallel runs. 

To run and merge using 4 processes, for instance, you could use: 
    $ mpiexec -n 4 python3 galerkin_kolmogorov_flow.py
    $ mpiexec -n 4 python3 -m dedalus merge_procs rb1
    
This script can restart the simulation from the last save of the original
output to extend the integration. This requires that the output files from the
original simulation are merged, and the last is symlinked or copied to 
'restart.h5'

"""
# SET UP ENVIRONMENT

import pathlib
import time
import h5py
import numpy as np
from mpi4py import MPI
from docopt import docopt
import os
from dedalus import public as de
from dedalus.extras import flow_tools
import logging
logger = logging.getLogger(__name__)

# DEFINE PARAMETERS

args = docopt(__doc__)
ReynoldsB    = float(args['--Rb'])                                                             # Buoyancy Reynolds Number
Prandtl      = float(args['--Pr'])                                                             # Prandtl Number
Froude       = float(args['--Fr'])                                                             # Froude Number
Lx, Lz       = (Froude*6.0*np.pi/0.34, 4.0*np.pi/3.0)                                          # Box Size
Nz           = int(args['--Nz'])                                                               # No. of Gridpoints
stop_time    = float(args['--Tend'])                                                           # Sim. stop time

# LOGGER: RECORD INPUT PARAMETERS

if MPI.COMM_WORLD.rank == 0:
    logger.info("Running 2DBSQ galerkin truncated simulation for Rb={:.3e}, Pr={:.3e}, Fr={:.3e}".format(ReynoldsB, Prandtl, Froude))

# CREATE RESULTS FOLDER

path = 'results_galerkin/'
if MPI.COMM_WORLD.rank == 0:
    if not os.path.exists(path):
        os.mkdir(path)    

# CREATE BASES AND DOMAIN

z_basis = de.Fourier("z", Nz, interval=(0, Lz), dealias=3/2)
domain  = de.Domain([z_basis], grid_dtype=np.complex128)
z       = domain.grid(0)

# FORCING TERM

def perforce(*args):
    z = args[0].data
    m = 3.0
    R = ReynoldsB
    return 1.0*((m*m*m)/R)*np.cos(m*z)

def Forcing(*args, domain=domain, F=perforce):
    return de.operators.GeneralFunction(domain, layout='g', func=F, args=args)

de.operators.parseables['Kolmogorov'] = Forcing

# PROBLEM SETUP

problem = de.IVP(domain, variables=['xi0','psi0','b0','xi19','psi19','b19'])

# EQUATION ENTRY SUBSTITUTIONS

problem.parameters['Rb']                    = ReynoldsB
problem.parameters['Pr']                    = Prandtl
problem.parameters['Fr']                    = Froude
problem.parameters['Lx']                    = Lx
problem.parameters['Lz']                    = Lz
problem.parameters['k']                     = 2.0*np.pi/Lx
problem.substitutions['Lap(A, n)']          = 'dz(dz(A)) - ((n*k*Fr)**2.0)*A'
problem.substitutions['dx(A, n)']           = '1j*n*k*A'
problem.substitutions['Jp(A, B, C, D)']     = '1j*k*(A*dz(B) + C*dz(D))'
problem.substitutions['Jm(A, B, C, D)']     = '1j*k*(A*dz(B) - C*dz(D))'
problem.substitutions['J0(A, B, n)']        = 'n*Jp(conj(A), B, B, conj(A)) - n*Jp(A, conj(B), conj(B), A)' 

# VELOCITY SUBSTITUTIONS

problem.substitutions['u0']  = 'dz(psi0)'
problem.substitutions['w0']  = '-1.0*dx(psi0, 0)'
problem.substitutions['u19'] = 'dz(psi19)'
problem.substitutions['w19'] = '-1.0*dx(psi19, 19)'

# EVOLUTION EQUATIONS 
    
problem.add_equation("dt(xi0) - (1/Rb)*Lap(xi0, 0) = Kolmogorov(z) + J0(xi19, psi19, 19)")
problem.add_equation("dt(b0)  - (1/(Pr*Rb))*Lap(b0, 0) = J0(b19, psi19, 19)")
problem.add_equation("dt(xi19) + dx(b19, 19) - (1/Rb)*Lap(xi19, 19) = -19*Jm(xi19, psi0, psi19, xi0)")
problem.add_equation("dt(b19) - dx(psi19, 19) - (1/(Pr*Rb))*Lap(b19, 19) = -19*Jm(b19, psi0, psi19, b0)") 

# CONSTRAINT EQUATIONS
    
problem.add_equation("xi0 - Lap(psi0, 0) = 0", condition="(nz != 0)")
problem.add_equation("psi0 = 0", condition="(nz == 0)")
problem.add_equation("xi19 - Lap(psi19, 19) = 0")

# BUILD SOLVER

ts = de.timesteppers.SBDF4
solver = problem.build_solver(ts)
logger.info('Solver built')

# INITIAL CONDITIONS 

xi0   = solver.state['xi0']
psi0  = solver.state['psi0']
b0    = solver.state['b0']

xi19  = solver.state['xi19']
psi19 = solver.state['psi19']
b19   = solver.state['b19']

# INITIALIZE NOISE IN PARALLEL SAFE MANNER

gshape = domain.dist.grid_layout.global_shape(scales=1)
slices = domain.dist.grid_layout.slices(scales=1)
rand   = np.random.RandomState(seed=42)
noise  = rand.standard_normal(gshape)[slices]
    
# INITIALIZE MEAN FIELDS WITH BASIC STATE

xi0['g']  = 3.0*np.cos(3.0*z)
psi0['g'] = -(1/3)*np.cos(3.0*z)

# GET EIGENFUNCTIONS FROM FILE

eigenfuncs_real = h5py.File('initialize_galerkin_real.h5','r')
eigenfuncs_imag = h5py.File('initialize_galerkin_imag.h5','r')

# INITIALIZE MODE 19 WITH EIGENFUNCTIONS

xi19['g']  = eigenfuncs_real.get('xi19')[slices] + 1j*eigenfuncs_imag.get('xi19')[slices] 
psi19['g'] = eigenfuncs_real.get('psi19')[slices] + 1j*eigenfuncs_imag.get('psi19')[slices]
b19['g']   = eigenfuncs_real.get('b19')[slices] + 1j*eigenfuncs_imag.get('b19')[slices]
    
# INTEGRATION PARAMETERS
    
dt = 0.00001
solver.stop_sim_time = stop_time
fh_mode = 'append'
    
# ANALYSIS

snapshot = solver.evaluator.add_file_handler(path+"/field_snapshots", sim_dt=0.01, max_writes=100, mode=fh_mode)
snapshot.add_task("u0", name  = 'u0')
snapshot.add_task("w0", name  = 'w0')
snapshot.add_task("b0", name  = 'b0')
snapshot.add_task("xi0", name = 'xi0')
snapshot.add_task("u19", name  = 'u19')
snapshot.add_task("w19", name  = 'w19')
snapshot.add_task("b19", name  = 'b19')
snapshot.add_task("xi19", name = 'xi19')

# CFL 

CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=10, safety=0.1, max_change=2, min_change=0.5, max_dt=0.05, threshold=0.05)
CFL.add_velocity('w19',0)

# FLOW TOOLS

flow = flow_tools.GlobalFlowProperty(solver, cadence = 100)
flow.add_property("xi0*conj(xi0)", name="vort")

# MAIN LOOP

try:  
    logger.info('Starting loop')
    start_time = time.time()
    while solver.proceed:
        
        timestep = CFL.compute_dt()
        solver.step(timestep)
        if (solver.iteration-1) % 5 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e, max(xi0): %e' %(solver.iteration, solver.sim_time, timestep, flow.max('vort'))) 
             
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
    logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))


