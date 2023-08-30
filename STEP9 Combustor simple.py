# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 15:10:47 2023

@author: gauth
"""

"""
we call it simple because : 
    no delay in combustion 
    we consider temperature of function of what is put in, what is bruned and that temperature is instenous given everywhere
"""

# =============================================================================
# Modules
# =============================================================================


import scipy as scipy
import numpy as np
import math
import matplotlib.pyplot as plt
import sys
from scipy.optimize import fsolve
import time 

a = 1355**2
rho = 1197

# =============================================================================
# Geometry
# =============================================================================


l_pipe = 2
d_pipe = 0.01
A_pipe = (d_pipe/2)**2 * 3.14
V_pipe = l_pipe * A_pipe

l_combustor = 0.2
d_combustor = 0.15 
A_combustor = (d_combustor/2)**2 * 3.1
V_combustor = (0.15/2)**2 * 3.14 *0.2

A_throat = A_combustor/8

# =============================================================================
# Constants
# =============================================================================

R = 8.734
R_cc = R/29


T3 = 2300 #K

t = 0
dt = 1e-4

type_simulation = 1
gamma = 1.3
mu = 1e-5
eta = 1e-3
tau = 0.03

P1i = 1e6
P2i = 1e6
P3i = 0.1e6

m1i = 0.
m2i = 0.
m3i = 0.

m1 = m1i
m2 = m2i
m3 = m3i

P1 = P1i
P2 = P2i
P3 = P3i
# dynamic viscosity
eta_oxygene = 42e-6
eta_propane = 117e-6 # Pa/s

density_oxygene = 1141
density_fuel = 500

bulk_modulus_oxygene = 0.2e9
bulk_modulus_fuel = 0.36e9

material_pipe = 200e9 #Pa , Young Modulus steel
material_dome = 200e9 #Pa , steel
material_injector = 200e9 #Pa , steel
material_combustor = 117e9 #Pa , copper

# =============================================================================
# Functions
# =============================================================================

def a_sound(density, bulk_modulus,length, material):
    # if liquid, option to return something cste
    return ( bulk_modulus / density ) / (1+(length/0.1)*(bulk_modulus/material)) # corrected speed of sound
    # return  bulk_modulus / density

def calculate_rho(Pressure_current,Pressure_init,rho_init):
    # if liquid, option to return something cste
    if type_simulation == 0 : 
        return rho_init * (Pressure_current/Pressure_init)**(1/gamma)
    else :
        return rho_init

def Reynolds(m_flow, length, section):
    m=m_flow
    if m_flow==0:
        m=1e-12
    # return m / (length*mu)
    return (m * length ) / ( section * eta )

def function(S):
    
    return np.array([
        S[0] - P1 ,
        S[1] - P2 , 
        S[3] - m1 - dt * A_pipe/l_pipe * (S[0] - S[2])  ,
        S[4] - m2 - dt * A_pipe/l_pipe * (S[1] - S[2])  ,
        # (S[2] * V_combustor ) / (  R_cc * S[6] ) * (S[2] - P3)/dt + S[5] * S[2] - np.sqrt(R_cc * S[6])/(np.sqrt(gamma * (2/(gamma+1) )**((gamma+1)/(gamma-1)) ) * A_throat) * (S[3] + S[4]) ,
        V_combustor / (R_cc * S[6]) * (S[2] - P3)/dt - 2 * S[3] + S[5]  ,
        S[5] - np.sqrt(gamma * (2/(gamma+1) )**((gamma+1)/(gamma-1)) ) * A_throat * S[2] / np.sqrt(R_cc * S[6])  ,
        S[6] - T3
        ])


list_pressure_1 = []
list_pressure_2 = []
list_pressure_3 = []

list_mass_flow_1 = []
list_mass_flow_2 = []
list_mass_flow_3 = []

time_storage = []

start = time.time() 
print("Start computing") 

# =============================================================================
# Loop
# =============================================================================

while t<2:
    
    ############################################################  
    #### LOOP AND CALCULUS 
    ############################################################
    
    t = t + dt

    root = fsolve(function, [P1, P2, P3, m1, m2, m3, T3])
    P1 = root[0]
    P2 = root[1]
    P3 = root[2]
    m1 = root[3]
    m2 = root[4]
    m3 = root[5]
    T3 = root[6]
    
    ############################################################  
    #### SAVE 
    ############################################################
    
    list_pressure_1.append(P1)
    list_pressure_2.append(P2)
    list_pressure_3.append(P3)

    list_mass_flow_1.append(m1)
    list_mass_flow_2.append(m2)
    list_mass_flow_3.append(m3)

    time_storage.append(t)
    
end = time.time() 
print(f'time of simulation : {end - start}')  
  
############################################################  
#### PLOTS 
############################################################


plt.figure()
plt.plot(time_storage[::2],list_pressure_1[::2])
plt.plot(time_storage[::2],list_pressure_2[::2])
plt.plot(time_storage[::2],list_pressure_3[::2])

plt.title(f'Pressure evolution in the filling of a cavity Combustor \n fed at constant pressure with water \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['P1 tank supply ox','P2 tank supply fuel','P3 combustion chamber pressure'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time_storage[::2],list_mass_flow_1[::2])
plt.plot(time_storage[::2],list_mass_flow_2[::2])
plt.plot(time_storage[::2],list_mass_flow_3[::2])

plt.title(f'Massflow evolution in the filling of a cavity Combustor \n fed at constant pressure with water \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['Tank n°1','Tank n°2','Outflow by nozzle'])
plt.ylabel(r'$mass flow$ (kg/s)')
# plt.ylim(-0, 50) 