# -*- coding: utf-8 -*-
"""
Created on Thu May 25 12:04:40 2023

@author: gauth
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import sys

P0 = 1e5

V0 = 0.1**2 * 3.14 * 0.2
rho_in = 1.14e-3
T_in = 273
P_in = 200e5
A1 = 0.025**2 * 3.14
gamma=1.3
t=0
dt=1
P=[]
P_prev = P0
M=32
T_prev=T_in 
R=8.138
def beta():
    return 1/rho_in * M / R / T_prev

C=1
A2 = 0.1**2 * 3.14
n_iteration = int(4e2)
time_simulation = n_iteration * dt
for i in range(n_iteration):
    t = t + dt
    T_prev = T_in * (P0/P_prev)**((1-gamma)/gamma)
    P_next = P_prev + dt * beta() /V0 * C * np.sign(-P_prev+P_in) * A2 / rho_in / (np.sqrt(-1+(A2/A1)**2)) * np.sqrt(2*rho_in*np.abs(P_prev-P_in))
    P_prev = P_next
    P.append(P_next)

    
plt.figure()   
A1=round(A1*1e5)
plt.plot(P, label=f'{A1} 10^-5 $m^2$')
plt.title(f'presssure evolution in filling of gas of a tank of volume ={V0} m^3 with P_in={P_in} Pa \n')
plt.xlabel(r'$time$ (s)')
plt.legend(title='surface hole entry')
plt.ylabel(r'$Pressure$ (Pa)')





