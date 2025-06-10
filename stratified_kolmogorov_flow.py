#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Usage: 
  stratified_kolmogorov_flow.py [--Rb=<buoyancy_reynolds> --Pr=<prandtl> --Fr=<froude> \
  --Nx=<Nx> --Nz=<Nz> --Tend=<stop_time>] 
  
Options:
  --Rb=<buoyancy_reynolds>          Buoyancy Reynolds Number [default: 50.0]
  --Pr=<prandtl>                    Prandtl Number [default: 1.0]
  --Fr=<froude>                     Froude Number [default: 0.02]
  --Nx=<Nx>                         Number of downwind modes [default: 128]
  --Nz=<Nz>                         Number of vertical modes [default: 128]
  --Tend=<stop_time>                Simulation stop time [default: 100.0]  
"""

"""
Dedalus script for the direct numerical simulation of the 2D Boussinesq Equations
to explore the dynamics of strongly stratified Kolmogorov flow.

This script uses a Fourier basis in x (PBCs) and a Chebyshev basis in z 
with nonhomogeneous Dirichlet BCs for the horizontal velocity 'u', vertical
velocity 'w', and buoyancy 'b'. 

This script can be run serially or in parallel, and uses the built-in analysis 
framework to save data in HDF5 files. The 'merge_procs' command can be used to 
merge distributed analysis sets from parallel runs. 

To run and merge using 4 processes, for instance, you could use: 
    $ mpiexec -n 4 python3 stratified_kolmogorov_flow.py
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
Lx, Lz       = (0.03, 2.0*np.pi/3.0)                                          # Box Size
Nx, Nz       = (int(args['--Nx']), int(args['--Nz']))                                          # No. of Gridpoints
stop_time    = float(args['--Tend'])                                                           # Sim. stop time

# LOGGER: RECORD INPUT PARAMETERS

if MPI.COMM_WORLD.rank == 0:
    logger.info("Running 2DBSQ simulation for Rb={:.3e}, Pr={:.3e}, Fr={:.3e}".format(ReynoldsB, Prandtl, Froude))

# CREATE RESULTS FOLDER

path = 'results_ecs/'
if MPI.COMM_WORLD.rank == 0:
    if not os.path.exists(path):
        os.mkdir(path)    

# CREATE BASES AND DOMAIN

x_basis = de.Fourier("x", Nx, interval=(0, Lx), dealias=3/2)
z_basis = de.Fourier("z", Nz, interval=(0, Lz), dealias=3/2)
domain  = de.Domain([x_basis, z_basis], grid_dtype=np.float64, mesh=[4,1])
x, z    = domain.grids(scales=1)

# FORCING TERM

def perforce(*args):
    z = args[0].data
    m = 3.0
    R = ReynoldsB
    return -1.0*((m*m*m)/R)*np.sin(m*z)

def Forcing(*args, domain=domain, F=perforce):
    return de.operators.GeneralFunction(domain, layout='g', func=F, args=args)

de.operators.parseables['Kolmogorov'] = Forcing

# PROBLEM SETUP

problem = de.IVP(domain, variables=['zeta','psi','b'])

# EQUATION ENTRY SUBSTITUTIONS

problem.parameters['Rb']          = ReynoldsB
problem.parameters['Pr']          = Prandtl
problem.parameters['Fr']          = Froude
problem.parameters['Lx']          = Lx
problem.parameters['Lz']          = Lz

problem.substitutions['Lap(A)']   = 'Fr*Fr*dx(dx(A)) + dz(dz(A))'
problem.substitutions['J(A, B)']  = 'dz(A)*dx(B) - dx(A)*dz(B)'

# ANALYSIS SUBSTITUTIONS

#---- OPERATORS

problem.substitutions['Vavg(A)'] = 'integ(integ(A, "x"), "z")/(Lx*Lz)'

#---- VELOCITIES

problem.substitutions['u'] = 'dz(psi)'
problem.substitutions['w'] = '-1.0*dx(psi)'

#---- TOTAL ENERGIES

problem.substitutions['KE']   = 'Vavg((u*u +  Fr*Fr*w*w)/2.0)'
problem.substitutions['PE']   = 'Vavg((b*b)/2.0)'
problem.substitutions['TE']   = 'KE + PE'

# EVOLUTION EQUATIONS 
    
problem.add_equation("dt(zeta) - (1/Rb)*Lap(zeta) + dx(b) = -J(psi, zeta) + Kolmogorov(z)")
problem.add_equation("dt(b) - (1/(Pr*Rb))*Lap(b) - dx(psi) = -J(psi, b)")

# CONSTRAINT EQUATIONS
    
problem.add_equation("zeta - Lap(psi) = 0", condition="(nx != 0) or (nz != 0)")
problem.add_equation("psi = 0", condition="(nx == 0) and (nz == 0)")

# BUILD SOLVER

ts = de.timesteppers.SBDF4
solver = problem.build_solver(ts)
logger.info('Solver built')

# INITIAL CONDITIONS OR RESTART CODE

if not pathlib.Path('restart.h5').exists():
    
    zeta = solver.state['zeta']
    psi  = solver.state['psi']
    b    = solver.state['b']
    
    # INITIALIZE NOISE IN PARALLEL SAFE MANNER

    gshape = domain.dist.grid_layout.global_shape(scales=1)
    slices = domain.dist.grid_layout.slices(scales=1)
    rand   = np.random.RandomState(seed=42)
    noise  = rand.standard_normal(gshape)[slices]

    # GET EIGENFUNCTIONS FROM FILE

    eigenfuncs = h5py.File('initialize_ecs.h5','r')

    # BACKGROUND CONFIGURATIONS + EIGENFUNCTIONS
    
    zeta['g'] = eigenfuncs.get('zeta')[slices]
    psi['g']  = eigenfuncs.get('psi')[slices]
    b['g']    = eigenfuncs.get('b')[slices]

    # INTEGRATION PARAMETERS
    
    dt = 1e-5
    solver.stop_sim_time = stop_time
    fh_mode = 'overwrite'
    
else:
    
    # RESTART

    write, last_dt = solver.load_state('restart.h5', -1)
    
    # INTEGRATION PARAMETERS

    dt = last_dt
    solver.stop_sim_time = 100
    fh_mode = 'append'
    
# ANALYSIS

snapshot = solver.evaluator.add_file_handler(path+"/field_snapshots", sim_dt=0.01, max_writes=100, mode=fh_mode)
snapshot.add_task("u", name = 'u')
snapshot.add_task("w", name = 'w')
snapshot.add_task("b", name = 'b')
snapshot.add_task("zeta", name = 'zeta')

globalp = solver.evaluator.add_file_handler(path+"/energy_timeseries", sim_dt=0.01, max_writes=10000000, mode=fh_mode)
globalp.add_task("KE", name = 'KE')
globalp.add_task("PE", name = 'PE')
globalp.add_task("TE", name = 'TE')

photo = solver.evaluator.add_file_handler(path+"/checkpointing_data", sim_dt=25, max_writes=100, mode=fh_mode)
photo.add_system(solver.state)

# CFL 

CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=10, safety=0.1, max_change=2, min_change=0.5, max_dt=0.05, threshold=0.05)
CFL.add_velocities(('u', 'w'))

# FLOW TOOLS

flow = flow_tools.GlobalFlowProperty(solver, cadence = 100)
flow.add_property("w", name="w2")
flow.add_property("TE", name="VTE")

# MAIN LOOP

try:  
    logger.info('Starting loop')
    start_time = time.time()
    while solver.proceed:
        
        timestep = CFL.compute_dt()
        # timestep = dt
        solver.step(timestep)
        if (solver.iteration-1) % 100 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e, max(w): %e, VTE: %e' %(solver.iteration, solver.sim_time, timestep, flow.max('w2'), flow.max('VTE')))     
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
    logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))


