# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 11:26:32 2023

@author: gauth

"""
"""
Comment : 
    modele with equation for the injector 0D => 1 equation without capacitance and inductance hypothesis
    full implicit 
    comparable to STEP7 but this time all is implicit and using fsolve !
    it is basically STEP7 improvement
    BOUNDARY CONDITION IS PRESSURE
"""

import scipy as scipy
import numpy as np
import math
import matplotlib.pyplot as plt
import sys
from scipy.optimize import fsolve

a = 1355**2
rho = 1197

Nombre_injector = 5

l_pipe = 2
l_dome = 0.1
l_injector = 0.05

d_injector = 0.003
d_pipe = 0.1
d_dome = 0.15

A_pipe = (d_pipe/2)**2 * 3.14
A_dome = (d_dome/2)**2 * 3.14
A_injector = (d_injector/2)**2 * 3.14

V_injector = l_injector * A_injector
V_pipe = l_pipe * A_pipe
V_dome = l_dome * A_dome

t = 0
dt = 0.00001

type_simulation = 1
gamma = 1
mu = 1e-5
eta = 1e-4
zeta0 = 1
zeta1 = (1-d_pipe/d_dome)**2
zeta2 = 2.85

P0i = 0.6e6
P1i = 0.1e6
P2i = 0.1e6
P3i = 0.1e6

m0i = 0
m1i = 0
m2i = 0
m3i = 0

m0 = m0i
m1 = m1i
m2 = m2i
m3 = m3i

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
    # xi0 = 64/Reynolds(m0,l_pipe,A_pipe)*1/d_pipe * 0.5 * A_pipe**(-2) + zeta0 * 0.5 * A_pipe**(-2)
    lambda0 = 64/Reynolds(m0,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi0 = 0.5 * A_pipe**(-2) * ( lambda0 + zeta0 )
    lambda1 = 64/Reynolds(m1,l_dome,A_dome)*1/d_dome *l_dome
    xi1 = 0.5 * A_dome**(-2) * ( lambda1 + zeta1 )
    # xi1 = 64/Reynolds(m1,l_dome,A_dome)*1/d_dome * 0.5 * A_dome**(-2) + zeta1 * 0.5 * A_dome**(-2)
    # QINJ = A_injector * zeta2 * np.sign(P2-P3) * np.sqrt(2*rho*np.abs(P2-P3))
    # m2 = Nombre_injector * QINJ
    
    return np.array([
        S[0] - P0,
        S[1] - P1 - dt * a / V_pipe *(S[4] - S[5]) ,
        S[2] - P2 - dt * a / V_dome *(S[5] - S[6]) ,
        S[3] - P3,
        S[4] - m0 - dt * A_pipe / l_pipe * ( S[0] - S[1] - xi0 / rho * S[4]*np.abs(S[4]) ) ,
        S[5] - m1 - dt * A_dome / l_dome * ( S[1] - S[2] - xi1 / rho * S[5]*np.abs(S[5]) )  ,
        S[6] - Nombre_injector * A_injector * zeta2 * np.sign(S[2]-S[3]) * np.sqrt(2*rho*np.abs(S[2]-S[3])) ,
        S[7] - S[6]
        ])


list_pressure_0 = []
list_pressure_1 = []
list_pressure_2 = []
list_pressure_3 = []

list_mass_flow_0 = []
list_mass_flow_1 = []
list_mass_flow_2 = []
list_mass_flow_3 = []

time = []


while t<0.03:
    
    ############################################################  
    #### LOOP AND CALCULUS 
    ############################################################
    
    t = t + dt
    root = fsolve(function, [P0, P1, P2, P3, m0, m1, m2, m3])
    P0 = root[0]
    P1 = root[1]
    P2 = root[2]
    P3 = root[3]
    m0 = root[4]
    m1 = root[5]
    m2 = root[6]
    m3 = root[7]
    
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
    
    time.append(t)
    
    
############################################################  
#### PLOTS 
############################################################


plt.figure()
plt.plot(time[::10],list_pressure_0[::10])
plt.plot(time[::10],list_pressure_1[::10])
plt.plot(time[::10],list_pressure_2[::10])
plt.plot(time[::10],list_pressure_3[::10])

plt.title(f'Pressure evolution Tank-Pipe-Dome_Injector 0D modelling\n Filled with water and constant output presssure in the chamber \n')
plt.title(f'Liquid Oxygen model : Pipe Dome Injector \n with constant pressure out {P3} Pa \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['Tank $P_0$','Pipe $P_1$','Dome $P_2$','Chamber $P_3$'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time[::10],list_mass_flow_0[::10])
plt.plot(time[::10],list_mass_flow_1[::10])
plt.plot(time[::10],list_mass_flow_2[::10])

plt.title(f'Massflow evolution Tank-Pipe-Dome_Injector 0D modelling\n Filled with water and constant output presssure in the chamber \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['Pipe $\dot{m}_0$','Dome $\dot{m}_1$','Injector & Chamber $\dot{m}_2$'])
plt.ylabel(r'$mass flow$ (kg/s)')