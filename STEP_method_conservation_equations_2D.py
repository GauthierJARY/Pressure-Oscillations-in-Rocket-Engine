# -*- coding: utf-8 -*-
"""
Created on Sat May 20 20:02:38 2023

@author: gauth
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import sys


"""
Fonction definition 
"""

def time_integration(u_prev, u_next, v_prev, v_next, H_prev, H_next, rho_prev, rho_next, p_prev, p_next):
    for i in range(1,nx-1):
        for j in range (1,ny-1):
            if True: #bordures[i,j]==0: 
                rho_next[i,j] = rho_prev[i,j]- dt * ( (u_prev[i+1,j]*rho_prev[i+1,j]-u_prev[i-1,j]*rho_prev[i-1,j])/(2*dx) + (v_prev[i,j+1]*rho_prev[i,j+1]-v_prev[i,j-1]*rho_prev[i,j-1])/(2*dy) ) 
                p_next[i,j] = s0[i,j] * ( rho_next[i,j]**(gamma) )
                u_next[i,j] = u_prev[i,j] - dt * ( u_prev[i,j] * (u_prev[i+1,j]-u_prev[i-1,j])/(2*dx) + v_prev[i,j] * (u_prev[i,j+1]-u_prev[i,j-1])/(2*dy) ) - 1/rho_prev[i,j] * dt *( ( p_prev[i+1,j]-p_prev[i-1,j] )/(2*dx) )
                v_next[i,j] = v_prev[i,j] - dt * ( u_prev[i,j] * (v_prev[i+1,j]-v_prev[i-1,j])/(2*dx) + v_prev[i,j] * (v_prev[i,j+1]-v_prev[i,j-1])/(2*dy) ) - 1/rho_prev[i,j] * dt *( ( p_prev[i,j+1]-p_prev[i,j-1] )/(2*dy) )           
                # H_next[i,j] = 1/rho_prev[i,j] * ( p_prev[i,j]-p_prev[i,j]) + 1

def boundary(u_next,v_next,rho_next, p_next,U0,R,gamma,s0,Dm):
    # Dirichlet conditions
    u_next[-1,:]=0
    u_next[0,:]=0
    v_next[:,-1]=0
    v_next[:,0]=0
    u_next[round(nx/2)-R:round(nx/2)+R,0]=0
    v_next[round(nx/2)-R:round(nx/2)+R,0]=-U0
    # Neumann conditions 
    rho_next[:,0]=rho_next[:,1]
    rho_next[0,:]=rho_next[1,:]
    rho_next[-1,:]=rho_next[-2,:]
    rho_next[:,-1]=rho_next[:,-2]
    # rho_next[round(nx/2)-R:round(nx/2)+R,0]=1
    # actualisation of pressure with new rho 
    p_next[:,:] = s0 * rho_next[:,:]**(gamma) 
   
"""
Loop main 
"""

nx=40 # mesh points on x axe
ny=40 # mesh points on y axe
dt = 0.1 # time step
dx = 1 # spatial step
dy = 1 # spatial step
s0=1 # initial entropy 
gamma=9/5 # reel gas model
U0=10 # inflow velocity
R=round(nx/8) # half diameter of inflow hole
print('nombre de maille = ' + str(nx*nx)+'\n')
cavity = np.zeros( (nx,ny) ) # space of work 
time_vector=0 # time stepper 
Dm=1 # mass inflow (kg/s)

# velocities on x and y 
u_prev = np.zeros( (nx,ny) )
u_next = u_prev.copy()
v_prev = np.zeros( (nx,ny) )
v_next = v_prev.copy()
# pressure field
p_prev = np.zeros( (nx,ny) )
p_next = p_prev.copy()
# rho field (surfacique)
rho_prev = 0.0708 * np.ones( (nx,ny) ) # with H2 volumical mass
rho_next = rho_prev.copy()
# energy field
H_prev = np.zeros( (nx,ny) )
H_next = H_prev.copy()

# mask for boundaries 
bordures=np.zeros((nx,ny))
bordures[0,:]=2 # left side of the cavity 
bordures[-1,:]=3 # right
bordures[:round(nx/2)-R,0]=1 #♥ top left side of the inflow hole 
bordures[round(nx/2)+R:,0]=1 # top  right side of the inflow hole 
bordures[:,-1]=4 # bottom
bordures[round(nx/2)-R:round(nx/2)+R,0]=5 # inflow hole of the cavity 

# initialisation 
boundary(u_prev,v_prev,rho_prev, p_prev,U0,R,gamma,s0,Dm)
s0=p_prev * rho_prev**(-gamma)
record=[[],[]]
time_max=100*dt

# loop 
while time_vector < time_max:
    
    sys.stdout.write("\rSimulation time advancement: {0:.2f} %".format(100*(time_vector)/time_max))
    sys.stdout.flush()
    time_vector+=dt # time incrementation along the time step 
    # call for fonction 
    time_integration(u_prev, u_next, v_prev, v_next, H_prev, H_next, rho_prev, rho_next, p_prev, p_next)
    # call for fonction
    boundary(u_next,v_next,rho_next, p_next,U0,R,gamma,s0,Dm) 
    # time looping and discreted time integration 
    u_prev = u_next # advancing in time 
    v_prev = u_next
    H_prev = H_next
    rho_prev = rho_next
    p_prev = p_next
    # data stocks for analysis along time advancement 
    record[0].append([p_next[10,5],p_next[30,-1],p_next[0,30],p_next[-1,30]])
    record[1].append(rho_next.sum())
print("\n end Simulation !")


# analysis and data plotting 
plt.figure()
plt.pcolormesh(v_next)
plt.colorbar()
plt.title('v_next repartition in cavity')
plt.axis('equal')

plt.figure()
plt.pcolormesh(u_next)
plt.colorbar()
plt.title('u_next repartition in cavity')
plt.axis('equal')

plt.figure()
plt.pcolormesh(rho_next/0.0708)
plt.colorbar()
plt.title('rho_next repartition in cavity')
plt.axis('equal')

plt.figure()
plt.pcolormesh(p_next)
plt.colorbar()
plt.title('pressure repartition in cavity')
plt.axis('equal')

plt.figure()
plt.semilogy(record[0])
plt.legend(['point 1','point 2','point 3','point 4'])
plt.title('pressure_evolution on specific point in the cavity')

plt.figure()
plt.plot(record[1])
plt.title('kind of mass evolution in the cavity')


fig=plt.figure() 
x=np.linspace(0,nx,nx)
y=np.linspace(0,ny,ny)
plt.streamplot(x,y,u_next.T,v_next.T)
plt.xlabel(r'$X$')
plt.ylabel(r'$Y$')
plt.title('flow of the fluid incoming in the cavity')
plt.axis('equal')

plt.figure() 
mesh=[np.arange(1,nx,1)]*nx
plt.plot(mesh, marker='o',markersize=0.2, color='k', linestyle='none')
plt.title('mesh of the cavity')
plt.axis('equal')

