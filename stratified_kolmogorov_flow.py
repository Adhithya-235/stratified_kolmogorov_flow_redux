#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Usage: 
  stratified_shear.py [--Re=<reynolds> --Pr=<prandtl> --Ri=<bulk_richardson> \
  --R=<sld_dld_ratio> --Sp=<layer_separation> --Nx=<Nx> --Nz=<Nz> --Tend=<stop_time>] 
  
Options:
  --Re=<reynolds>                   Reynolds Number [default: 1000.0]
  --Pr=<prandtl>                    Prandtl Number [default: 1.0]
  --Ri=<bulk_richardson>            Bulk Richardson Number [default: 1.0]
  --R=<sld_dld_ratio>               SLD DLD Ratio [default: 3]
  --Sp=<layer_separation>           Shear Layer Separation [default: 1]
  --Nx=<Nx>                         Number of downwind modes [default: 128]
  --Nz=<Nz>                         Number of vertical modes [default: 128]
  --Tend=<stop_time>                Simulation stop time [default: 100.0]  
"""

"""
Dedalus script for the direct numerical simulation of the 2D Boussinesq Equations
to explore the dynamics of stacked stratified shear layers.

This script uses a Fourier basis in x (PBCs) and a Chebyshev basis in z 
with nonhomogeneous Dirichlet BCs for the horizontal velocity 'u', vertical
velocity 'w', and buoyancy 'b'. 

This script can be run serially or in parallel, and uses the built-in analysis 
framework to save data in HDF5 files. The 'merge_procs' command can be used to 
merge distributed analysis sets from parallel runs. 

To run and merge using 4 processes, for instance, you could use: 
    $ mpiexec -n 4 python3 stratified_shear.py
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
Reynolds     = float(args['--Re'])                                                             # Reynolds Number
Prandtl      = float(args['--Pr'])                                                             # Prandtl Number
Richardson   = float(args['--Ri'])                                                             # Bulk Richardson Number
R            = float(args['--R'])                                                              # Layer Thickness Ratio
Sp           = float(args['--Sp'])                                                             # Shear Layer Separation
Lx, Lz       = (4.0*np.pi, 20.0)                                                               # Box Size
Nx, Nz       = (int(args['--Nx']), int(args['--Nz']))                                          # No. of Gridpoints
stop_time    = float(args['--Tend'])                                                           # Sim. stop time

# LOGGER: RECORD INPUT PARAMETERS

if MPI.COMM_WORLD.rank == 0:
    logger.info("Running 2DBSQ simulation for Re={:.3e}, Pr={:.3e}, Ri={:.3e}, R={:.3e}, Sp={:.3e}".format(Reynolds, Prandtl, Richardson, R, Sp))

# CREATE RESULTS FOLDER

path = 'results_mode2/'
if MPI.COMM_WORLD.rank == 0:
    if not os.path.exists(path):
        os.mkdir(path)    

# CREATE BASES AND DOMAIN

x_basis = de.Fourier("x", Nx, interval=(0, Lx), dealias=3/2)
z_basis = de.Chebyshev("z", Nz, interval=(-0.5*Lz, 0.5*Lz), dealias=3/2)
domain  = de.Domain([x_basis, z_basis], grid_dtype=np.float64, mesh=[24,1])
x, z    = domain.grids(scales=1)
xd, zd  = domain.grids(scales=3/2)

# PROBLEM SETUP

problem = de.IVP(domain, variables=['u','w','b','p','uz','wz','bz'])
problem.meta['u', 'w', 'b', 'p']['z']['dirichlet'] = True

# EQUATION ENTRY SUBSTITUTIONS

problem.parameters['Re']   = Reynolds
problem.parameters['Pr']   = Prandtl
problem.parameters['Ri']   = Richardson
problem.parameters['Rd']   = R
problem.parameters['Hh']   = Sp
problem.parameters['Lx']   = Lx
problem.parameters['Lz']   = Lz

problem.substitutions['Lap(A, Az)']       = 'dx(dx(A)) + dz(Az)'
problem.substitutions['Adv(A, B, C, Cz)'] = 'A*dx(C) + B*Cz'

# ANALYSIS SUBSTITUTIONS

#---- OPERATORS

problem.substitutions['savg(A)'] = 'integ(A, "x")/Lx'
problem.substitutions['Vavg(A)'] = 'integ(integ(A, "x"), "z")/(Lx*Lz)'
problem.substitutions['sper(A)'] = 'A - savg(A)'
problem.substitutions['rms(A)']  = 'sqrt(savg((sper(A))**2))'

#---- AVERAGES AND PERTURBATIONS

problem.substitutions['ub'] = 'savg(u)'
problem.substitutions['wb'] = 'savg(w)'
problem.substitutions['bb'] = 'savg(b)'
problem.substitutions['up'] = 'sper(u)'
problem.substitutions['wp'] = 'sper(w)'
problem.substitutions['bp'] = 'sper(b)'

#---- VORTICITY

problem.substitutions['omega_y'] = 'uz - dx(w)'

#---- TOTAL ENERGIES

problem.substitutions['KE']   = 'Vavg((u*u +  w*w)/2.0)'
problem.substitutions['PE']   = 'Vavg((b*b)/2.0)'
problem.substitutions['TE']   = 'KE + Ri*PE'

#---- PERTURBATION ENERGIES

problem.substitutions['PKE']  = 'Vavg((up*up +  wp*wp)/2.0)'
problem.substitutions['PPE']  = 'Vavg((bp*bp)/2.0)'
problem.substitutions['PTE']  = 'PKE + Ri*PPE'

# EVOLUTION EQUATIONS 
    
problem.add_equation("dt(u) + dx(p) - (Re**(-1.0))*Lap(u, uz) = -Adv(u, w, u, uz)") 
problem.add_equation("dt(w) + dz(p) - (Re**(-1.0))*Lap(w, wz) -Ri*b = - Adv(u, w, w, wz)")
problem.add_equation("dt(b) - ((Pr*Re)**(-1.0))*Lap(b, bz) = - Adv(u, w, b, bz)")

# CONSTRAINT EQUATIONS
    
problem.add_equation("dx(u) + wz = 0")

# AUXILIARY EQUATIONS (DEFINE z-DERIVATIVES)

problem.add_equation("uz - dz(u) = 0")
problem.add_equation("wz - dz(w) = 0")
problem.add_equation("bz - dz(b) = 0")

# BOUNDARY CONDITIONS

problem.add_bc("left(u)  = 0")
problem.add_bc("right(u) = 2")
problem.add_bc("left(b)  = 0")
problem.add_bc("right(b) = 2")
problem.add_bc("left(w)  = 0")
problem.add_bc("right(w) = 0", condition="(nx != 0)")
problem.add_bc("left(p)  = 0", condition="(nx == 0)")

# BUILD SOLVER

ts = de.timesteppers.SBDF3
solver = problem.build_solver(ts)
logger.info('Solver built')

# INITIAL CONDITIONS OR RESTART CODE

if not pathlib.Path('restart.h5').exists():
    
    u   = solver.state['u']
    uz  = solver.state['uz']
    w   = solver.state['w']
    wz  = solver.state['wz']
    b   = solver.state['b']
    bz  = solver.state['bz']

    # INITIALIZE NOISE IN PARALLEL SAFE MANNER

    gshape = domain.dist.grid_layout.global_shape(scales=1)
    slices = domain.dist.grid_layout.slices(scales=1)
    rand   = np.random.RandomState(seed=42)
    noise  = rand.standard_normal(gshape)[slices]

    # GET EIGENFUNCTIONS FROM FILE

    eigenfuncs = h5py.File('initialize_mode2.h5','r');

    # BACKGROUND CONFIGURATIONS + SMALL VERTICAL VELOCITY PERTURBATIONS LOCALIZED AROUND SHEAR
    
    # u['g'] = np.tanh(z-Sp-1.0) + np.tanh(z+Sp+1.0)
    # b['g'] = (np.tanh(R*(z-Sp-1.0)) + np.tanh(R*(z+Sp+1.0)))
    # w['g'] = 0.01*noise*(np.exp(-((z-Sp-1.0)**2)) + np.exp(-((z+Sp+1.0)**2)))

    # BACKGROUND CONFIGURATIONS + EIGENFUNCTIONS
    
    u['g'] = eigenfuncs.get('u')[slices]
    b['g'] = eigenfuncs.get('b')[slices]
    w['g'] = eigenfuncs.get('w')[slices]

    # DIFFERENTIATE IN z

    u.differentiate('z',out=uz)
    w.differentiate('z',out=wz)
    b.differentiate('z',out=bz)

    # INTEGRATION PARAMETERS
    
    dt = 1e-2
    solver.stop_sim_time = stop_time
    fh_mode = 'overwrite'
    
else:
    
    # RESTART

    write, last_dt = solver.load_state('restart.h5', -1)
    
    # INTEGRATION PARAMETERS

    dt = last_dt
    solver.stop_sim_time = 700
    fh_mode = 'append'
    
# ANALYSIS

snapshot = solver.evaluator.add_file_handler(path+"/field_snapshots", sim_dt=0.1, max_writes=100, mode=fh_mode)
snapshot.add_task("u", name = 'u')
snapshot.add_task("w", name = 'w')
snapshot.add_task("b", name = 'b')
snapshot.add_task("p", name = 'p')
snapshot.add_task("omega_y", name = 'o')

globalp = solver.evaluator.add_file_handler(path+"/energy_timeseries", sim_dt=0.02, max_writes=10000000, mode=fh_mode)
globalp.add_task("KE", name = 'KE')
globalp.add_task("PE", name = 'PE')
globalp.add_task("TE", name = 'TE')

photo = solver.evaluator.add_file_handler(path+"/checkpointing_data", sim_dt=25, max_writes=100, mode=fh_mode)
photo.add_system(solver.state)

# CFL 

CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=10, safety=0.2, max_change=1.5, min_change=0.5, max_dt=1e-2, threshold=0.1)
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


