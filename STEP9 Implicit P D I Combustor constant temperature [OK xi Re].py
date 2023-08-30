# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 15:17:00 2023

@author: gauth


Modele with 
"""

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
# Geometric parametrisation 
# =============================================================================

l_pipe = 2
l_dome = 0.1
l_injector = 0.05
l_combustor = 0.2

Nombre_injector = 5

d_injector = 0.003
d_pipe = 0.08
d_dome = 0.15
d_combustor = 0.15 

A_pipe = (d_pipe/2)**2 * 3.14
A_dome = (d_dome/2)**2 * 3.14
A_injector = (d_injector/2)**2 * 3.14

V_injector = l_injector * A_injector
V_pipe = l_pipe * A_pipe
V_dome = l_dome * A_dome

V_combustor = (0.15/2)**2 * 3.14 *0.2
R_combustor = 8.734
A_throat = 0.1

# =============================================================================
# Constant values 
# =============================================================================

T3 = 2500 #K

t = 0
dt = 0.00001

type_simulation = 1
gamma = 1.3
mu = 1e-5
eta = 1e-3

zeta1 = (1-d_pipe/d_dome)**2
zeta2 = 30
zeta3 =  2.85

# =============================================================================
# Initial values 
# =============================================================================

P0i = 11e6
P1i = 1e6
P2i = 1e6
P3i = 0.1e6

m0i = 0
m1i = 0
m2i = 0
m3i = 0
m4i = 0

m0 = m0i
m1 = m1i
m2 = m2i
m3 = m3i
m4 = m4i

P0 = P0i
P1 = P1i
P2 = P2i
P3 = P3i

def a_sound(Pressure_current,Pressure_init,rho_init):
    # if liquid, option to return something cste
    if type_simulation == 0 :
        return gamma * Pressure_current  / (rho_init * (Pressure_current/Pressure_init)**(1/gamma) )
    else : 
        # return ( 2.2e9 / rho_init ) / (1+(l_pipe/0.1)*(2.2e9/50e9))
        return 2.2e9 / rho_init

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
    
    lambda1 = 64/Reynolds(m1,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi1 = 0.5 * A_pipe**(-2) * ( lambda1 + zeta1 )
    lambda2 = 64/Reynolds(m2,l_dome,A_dome)*1/d_dome *l_dome
    xi2 = 0.5 * A_dome**(-2) * ( lambda2 + zeta2 )
    lambda3 = 64/Reynolds(m3/ Nombre_injector ,l_injector,A_injector)*1/d_injector *l_injector
    xi3 = 0.5 * A_injector**(-2) * ( lambda3 + zeta3 )  
    
    return np.array([
        S[0] - P0,
        S[4] - S[5],
        S[1] - P1 - dt * a / V_dome *(S[5] - S[6]) ,
        S[2] - P2 - dt * a / V_injector *(S[6] - S[7]) ,
        (S[3] * V_combustor) / (S[8] * R_combustor * T3) * (S[3]-P3)/dt + S[3] - S[3]*S[7]/S[8] ,
        S[5] - m1 - dt * A_pipe / l_pipe * ( S[0] - S[1] - xi1 / rho * S[5]*np.abs(S[5]) ) ,
        S[6] - m2 - dt * A_dome / l_dome * ( S[1] - S[2] - xi2 / rho * S[6]*np.abs(S[6]) )  ,
        ( S[7] - m3 ) / Nombre_injector - dt * A_injector/l_injector * ( S[2] - S[3] - xi3 / rho * (S[7]*np.abs(S[7]) ) / (Nombre_injector**2) ) ,
        S[8] - np.sqrt(gamma * (2/(gamma+1) )**((gamma+1)/(gamma-1)) ) * S[3]*A_throat / np.sqrt( R_combustor * T3 )
        ])


list_pressure_0 = []
list_pressure_1 = []
list_pressure_2 = []
list_pressure_3 = []

list_mass_flow_0 = []
list_mass_flow_1 = []
list_mass_flow_2 = []
list_mass_flow_3 = []
list_mass_flow_4 = []

time_storage = []

start = time.time() 
print("Start computing") 

while t<0.15:
    
    ############################################################  
    #### LOOP AND CALCULUS 
    ############################################################
    
    t = t + dt
    root = fsolve(function, [P0, P1, P2, P3, m0, m1, m2, m3, m4+1])
    P0 = root[0]
    P1 = root[1]
    P2 = root[2]
    P3 = root[3]
    m0 = root[4]
    m1 = root[5]
    m2 = root[6]
    m3 = root[7]
    m4 = root[8]
    
    ############################################################  
    #### SAVE 
    ############################################################
    
    list_pressure_0.append(P0)
    list_pressure_1.append(P1)
    list_pressure_2.append(P2)
    list_pressure_3.append(P3)

    list_mass_flow_0.append(m0)
    list_mass_flow_1.append(m1)
    list_mass_flow_2.append(m2)
    list_mass_flow_3.append(m3)
    list_mass_flow_4.append(m4)

    
    time_storage.append(t)
    
end = time.time() 
print(f'time of simulation : {end - start}')    
############################################################  
#### PLOTS 
############################################################


plt.figure()
plt.plot(time_storage[::10],list_pressure_0[::10])
plt.plot(time_storage[::10],list_pressure_1[::10])
plt.plot(time_storage[::10],list_pressure_2[::10])
plt.plot(time_storage[::10],list_pressure_3[::10])

plt.title(f'Liquid Oxygen model : Pipe Dome Injector Combustor \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['P0','P1','P2','P3'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time_storage[::10],list_mass_flow_0[::10])
plt.plot(time_storage[::10],list_mass_flow_1[::10])
plt.plot(time_storage[::10],list_mass_flow_2[::10])
plt.plot(time_storage[::10],list_mass_flow_3[::10])
plt.plot(time_storage[1::10],list_mass_flow_4[1::10])

plt.title(f'mass evolution in Pipe Dome Injector Combustor System \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['m0','m1','m2','m3','m4'])
plt.ylabel(r'$mass flow$ (kg/s)')
plt.ylim(0, 5) 