# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 17:40:56 2023

@author: gauth
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


l_pipe = 2
l_combustor = 0.2

d_pipe = 0.01
d_combustor = 0.15 

A_pipe = (d_pipe/2)**2 * 3.14

V_pipe = l_pipe * A_pipe

V_combustor = (0.15/2)**2 * 3.14 *0.2
R = 8.734
R_cc = R/16
A_throat = 0.1

T3 = 3500 #K

t = 0
dt = 0.00001

type_simulation = 1
gamma = 1.3
mu = 1e-5
eta = 1e-3
tau = 0.03

P1i = 11e6
P2i = 11e6
P3i = 0.1e6

m1i = 0.
m2i = 0.
m3i = 0.
mbi = 0.

m1 = m1i
m2 = m2i
m3 = m3i
mb = mbi

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
        (S[2] * V_combustor ) / (  R_cc * S[7] ) * (S[2] - P3)/dt + S[5] * S[2] - S[2]*S[6] ,
        S[6] - ( (S[3] - m1)/dt * (t - tau) + S[3] ) - ( (S[4] - m2)/dt * (t - tau) + S[4] ) ,
        S[5] - np.sqrt(gamma * (2/(gamma+1) )**((gamma+1)/(gamma-1)) ) * A_throat * S[2] / np.sqrt(R_cc * S[7])  ,
        S[7] - T3
        ])


list_pressure_1 = []
list_pressure_2 = []
list_pressure_3 = []

list_mass_flow_1 = []
list_mass_flow_2 = []
list_mass_flow_3 = []
list_mass_flow_b = []

time_storage = []

start = time.time() 
print("Start computing") 

while t<0.6:
    
    ############################################################  
    #### LOOP AND CALCULUS 
    ############################################################
    
    t = t + dt

    root = fsolve(function, [P1, P2, P3, m1, m2, m3, mb , T3])
    P1 = root[0]
    P2 = root[1]
    P3 = root[2]
    m1 = root[3]
    m2 = root[4]
    m3 = root[5]
    mb = root[6]
    T3 = root[7]
    
    ############################################################  
    #### SAVE 
    ############################################################
    
    list_pressure_1.append(P1)
    list_pressure_2.append(P2)
    list_pressure_3.append(P3)

    list_mass_flow_1.append(m1)
    list_mass_flow_2.append(m2)
    list_mass_flow_3.append(m3)
    list_mass_flow_b.append(mb)

    
    time_storage.append(t)
    
end = time.time() 
print(f'time of simulation : {end - start}')    
############################################################  
#### PLOTS 
############################################################


plt.figure()
plt.plot(time_storage[::10],list_pressure_1[::10])
plt.plot(time_storage[::10],list_pressure_2[::10])
plt.plot(time_storage[::10],list_pressure_3[::10])

plt.title(f'Pressure evolution Combustor fed at constant pressure \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['P1 tank supply fuel','P2 tank supply ox','P3 combustion chamber pressure'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time_storage[::10],list_mass_flow_1[::10])
plt.plot(time_storage[::10],list_mass_flow_2[::10])
plt.plot(time_storage[::10],list_mass_flow_3[::10])
plt.plot(time_storage[::10],list_mass_flow_b[::10])

plt.title(f'mass evolution in Combustor fed at constant pressure \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['m1 fuel','m2 Ox','m3 out','mb burnt gases'])
plt.ylabel(r'$mass flow$ (kg/s)')
# plt.ylim(-0, 50) 