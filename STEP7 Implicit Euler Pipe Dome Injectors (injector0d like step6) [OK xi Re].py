# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 11:25:55 2023

@author: gauth


"""

"""
Comment : 
    improvement of step6 where it was explicit and could not run because of overflow or lack of convergence
    Modele with equation for the injector 0D => 1 equation without capacitance and inductance hypothesis
    partly implicit 
    partly implicit because m_injecotr calculated out of the fsolve solver because m_injector is approximated by mn and not mn+1 like the others
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
d_pipe = 0.2
d_dome = 0.4

A_pipe = (d_pipe/2)**2 * 3.14
A_dome = (d_dome/2)**2 * 3.14
A_injector = (d_injector/2)**2 * 3.14

V_injector = l_injector * A_injector
V_pipe = l_pipe * A_pipe
V_dome = l_dome * A_dome

t = 0
dt = 0.0001

type_simulation = 1
gamma = 1
mu = 1e-5
eta = 1e-3

zeta0 = 1
zeta1 = 0.36
zeta2 = 1 + 0.36

P0i = 0.6e6
P1i = 0.4e6
P2i = 0.4e6
P3i = 1e5

m0i = 0
m1i = 0
m2i = 0
m3i = 0

m0 = m0i
m1 = m1i
m2 = m2i

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
    return (m * length ) / ( section * eta )

def function(S):
    lambda0 = 64/Reynolds(m0,l_pipe,A_pipe)*1/d_pipe *l_pipe
    xi0 = 0.5 * A_pipe**(-2) * ( lambda0 + zeta0 )
    lambda1 = 64/Reynolds(m1,l_dome,A_dome)*1/d_dome *l_dome
    xi1 = 0.5 * A_dome**(-2) * ( lambda1 + zeta1 )
    QINJ = A_injector * zeta2 * np.sign(P2-P3) * np.sqrt(2*rho*np.abs(P2-P3))
    m2 = Nombre_injector * QINJ
    
    return np.array([
        S[0] - P1 - dt * a / V_pipe *(S[2] - S[3]) ,
        S[1] - P2 - dt * a / V_dome *(S[3] - m2) ,
        S[2] - m0 - dt * A_pipe / l_pipe * ( P0 - S[0] - xi0 / rho * S[2]*np.abs(S[2]) ) ,
        S[3] - m1 - dt * A_dome / l_dome * ( S[0] - S[1] - xi1 / rho * S[3]*np.abs(S[3]) )  ])


list_pressure_0 = []
list_pressure_1 = []
list_pressure_2 = []
list_pressure_3 = []

list_mass_flow_0 = []
list_mass_flow_1 = []
list_mass_flow_2 = []
list_mass_flow_3 = []

time = []


while t<0.2:
    
    ############################################################  
    #### LOOP AND CALCULUS 
    ############################################################
    
    t = t + dt
    root = fsolve(function, [1, 1,1,1])
    P1 = root[0]
    P2 = root[1]
    m0 = root[2]
    m1 = root[3]
    QINJ = A_injector * zeta2 * np.sign(P2-P3) * np.sqrt(2*rho*np.abs(P2-P3))
    m2 = Nombre_injector * QINJ
    
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

plt.title(f'Liquid Oxygen model : Pipe Dome \n with constant ouput pressure P3 ={P3} Pa \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['P0','P1','P2','P3'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time[::10],list_mass_flow_0[::10])
plt.plot(time[::10],list_mass_flow_1[::10])
plt.plot(time[::10],list_mass_flow_2[::10])

plt.title(f'mass evolution in Pipe Dome Injector System \n with constant ouput pressure P3 ={P3} Pa \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['m0','m1','m2','m3'])
plt.ylabel(r'$mass flow$ (kg/s)')