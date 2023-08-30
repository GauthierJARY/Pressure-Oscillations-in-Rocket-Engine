# -*- coding: utf-8 -*-
"""
Created on Fri May 26 12:14:34 2023

@author: gauth
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import sys

P0 = 1e5
V0 = 0.1**2 * 3.14 * 0.2
rho_in = 1141
T_in = 90
P_in = 100e5
A1 = 0.025**2 * 3.14
gamma=1.3
t=0
dt=0.000001
P1=[]
P_prev = P0
C=1
M=32
R=8.138
A2 = 0.1**2 * 3.14
n_iteration = int(18e3)
time_simulation = n_iteration * dt
stockQ1=[]
stockQ2=[]
P2=[]
P_prev1 = 100e5
P_prev2 = 1e5

for i in range(200):
    t = t + dt
    beta1 = round(gamma * P_prev1)
    beta2 = round(gamma * P_prev2)
    P_pipe_entry = P_prev1
    P_pipe_end = P_prev2


    Q_in1 = C/np.sqrt(  (A2/A1)**2-1  ) * A2 * np.sign(-P_prev1+P_prev2) * np.sqrt( 2 * rho_in * np.abs(-P_prev1+P_prev2) )
    Q_in2 = C/np.sqrt(  (A2/A1)**2-1  ) * A2 * np.sign(-P_prev2+P_prev1) * np.sqrt( 2 * rho_in * np.abs(-P_prev2+P_prev1) )
    P_loss = Q_in2 * 64/300 * 1/0.1 * A1
    P_next1 = P_prev1 + dt * beta1 / V0 * Q_in1
    P_next2 = P_prev2 + dt * beta2 / V0 * Q_in2
        #P_next = P_prev + dt * beta/V0 * C * np.sign(-P_prev+P_in) * A2 / rho_in / (np.sqrt(-1+(A2/A1)**2)) * np.sqrt(2*rho_in*np.abs(P_prev-P_in))
        #P_next2 = P_prev2 + dt * beta/V0 * C * np.sign(-P_prev2+P_in) * A2 / rho_in / (np.sqrt(-1+(A2/A12)**2)) * np.sqrt(2*rho_in*np.abs(P_prev2-P_in))
 
    P_prev1 = P_next1
    P_prev2 = P_next2
    P1.append(P_next1)
    P2.append(P_next2)
    stockQ1.append(Q_in1)
    stockQ2.append(Q_in2)
    


    
plt.figure()   
plt.plot(P1)
plt.plot(P2)
plt.title(f'presssure evolution in filling of LOx of a tank of volume ={V0} m^3 with P_in={P_in} Pa \n')
plt.xlabel(r'$time$ (micro_s)')
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()   
plt.plot(stockQ1)
plt.plot(stockQ2)
plt.title(f'inflow evolution in filling of LOx of a tank of volume ={V0} m^3 with P_in={P_in} Pa \n')
plt.xlabel(r'$time$ (micro_s)')
plt.ylabel(r'$volumic flow$ ($m^3$/s)')