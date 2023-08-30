# -*- coding: utf-8 -*-
"""
Created on Fri May 26 10:44:44 2023

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
dt=0.0000001
P=[]
P_prev = P0
C=1
M=32
R=8.138
A2 = 0.1**2 * 3.14
n_iteration = int(18e3)
time_simulation = n_iteration * dt
stockQ=[]
while t<0.00010:
    t = t + dt
    beta = round(gamma * P_prev)
    Q_in = C/np.sqrt(  (A2/A1)**2-1  ) * A2 * np.sign(-P_prev+P_in) * np.sqrt( 2 * (gamma/(-1+gamma))* rho_in * np.abs(-P_prev+P_in) )
    rho2 = (P_prev/P_in)**(1/gamma) * rho_in
    Q_in = 100
    # Q_in = np.sqrt( ( (gamma/(1-gamma)) *2 * (P_prev/rho2 - P_in/rho_in) ) / ( -( (rho2*A2) / (rho_in * A1 ) )**2 + 1 ) )   
    stockQ.append(Q_in)
    P_next = P_prev + dt * beta / V0 * Q_in 
    #P_next = P_prev + dt * beta/V0 * C * np.sign(-P_prev+P_in) * A2 / rho_in / (np.sqrt(-1+(A2/A1)**2)) * np.sqrt(2*rho_in*np.abs(P_prev-P_in))
        #P_next2 = P_prev2 + dt * beta/V0 * C * np.sign(-P_prev2+P_in) * A2 / rho_in / (np.sqrt(-1+(A2/A12)**2)) * np.sqrt(2*rho_in*np.abs(P_prev2-P_in))
 
    P_prev = P_next
    P.append(P_next)


    
plt.figure()   
beta=round(beta/1e6)
plt.plot(P, label=f'{beta} MPa')
plt.title(f'presssure evolution in filling of gas Ox of a tank of volume ={V0} m^3 with P_in={P_in} Pa \n')
plt.xlabel(r'$time$ (micro_s)')
plt.legend(title='final bulk modulus')
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()   
beta=round(beta/1e6)
plt.plot(stockQ, label=f'volumic flow')
plt.title(f'inflow evolution in filling of gas Ox of a tank of volume ={V0} m^3 with P_in={P_in} Pa \n')
plt.xlabel(r'$time$ (micro_s)')
plt.ylabel(r'$volumic flow$ ($m^3$/s)')
