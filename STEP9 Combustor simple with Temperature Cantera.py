# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 13:22:53 2023

@author: gauth
"""



import scipy as scipy
import numpy as np
import math
import matplotlib.pyplot as plt
import sys
from scipy.optimize import fsolve
import time 
import cantera as ct 
import numpy as np


a = 1355**2
rho = 1197


l_pipe = 2
l_combustor = 0.2

d_pipe = 0.01
d_combustor = 0.15 

A_pipe = (d_pipe/2)**2 * 3.14

V_pipe = l_pipe * A_pipe

V_combustor = (0.15/2)**2 * 3.14 *0.2
R = 8.734
R_cc = R/16
A_throat = 0.01

T3 = 300 #K

t = 0
dt = 1e-5

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

list_temperature_3 = []

time_storage = []

start = time.time() 
print("Start computing") 

gas1 = ct.Solution('gri30.yaml')

while t<0.6:
    
    ############################################################  
    #### LOOP AND CALCULUS 
    ############################################################
    
    t = t + dt
    if m2 == 0 :
        Mixture_ratio = 1
    else:
        Mixture_ratio = np.abs(m1/m2)
    gas1.Y = f'H2:1 , O2:{Mixture_ratio}'
    gas1.TP = 300, P3
    gas1.equilibrate('HP')
    T3n = gas1.T
    P1n = P1
    P2n = P2
    m1n = m1 + dt * A_pipe/l_pipe * (P1 - P3)
    m2n = m2 + dt * A_pipe/l_pipe * (P2 - P3)
    # P3n = P3 + dt * R_cc * T3 / V_combustor * ( m1 + m2 - m3)
    P3n = P3 + dt * R_cc * T3 / V_combustor * ( m1 + m2 - m3 + V_combustor * P3 / R_cc * (T3n - T3)/dt * (+1 * 1/T3**2) )
    m3n = A_throat * np.sqrt(gamma * (2/(gamma+1) )**((gamma+1)/(gamma-1)) ) / np.sqrt(R_cc * T3n) * P3n
    P1 = P1n
    P2 = P2n
    P3 = P3n
    m1 = m1n
    m2 = m2n
    m3 = m3n
    T3 = T3n
    
    ############################################################  
    #### SAVE 
    ############################################################
    
    list_pressure_1.append(P1)
    list_pressure_2.append(P2)
    list_pressure_3.append(P3)

    list_mass_flow_1.append(m1)
    list_mass_flow_2.append(m2)
    list_mass_flow_3.append(m3)
    
    list_temperature_3.append(T3)
    
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

plt.title(f'Pressure evolution Combustor \n fed at constant P1 = {P1}, P2 = {P2}, and changing T3 \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['P1 tank supply ox','P2 tank supply fuel','P3 combustion chamber pressure'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time_storage[::2],list_mass_flow_1[::2])
plt.plot(time_storage[::2],list_mass_flow_2[::2])
plt.plot(time_storage[::2],list_mass_flow_3[::2])

plt.title(f'massflow evolution in Combustor \n fed at constant P1 = {P1}, P2 = {P2}, and changing T3 \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['m1 ox','m2 fuel','m3 out'])
plt.ylabel(r'$mass flow$ (kg/s)')
# plt.ylim(-0, 50)

plt.figure()
plt.plot(time_storage[::2],list_temperature_3[::2]) 
plt.title(f'temperature evolution in Combustor \n fed at constant P1 = {P1}, P2 = {P2}, and changing T3 \n')
plt.xlabel(r'$time$ (s)')
plt.ylabel(r'$T$ (K)')
