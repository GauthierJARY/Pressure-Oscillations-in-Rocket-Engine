# -*- coding: utf-8 -*-
"""
Created on Wed May 24 09:10:37 2023

@author: gauth
"""


import numpy as np
import math
import matplotlib.pyplot as plt
import sys

P0 = 1e5
T0 = 273
V0 = 0.1**2 * 3.14 * 0.2
V0 = 0.00628

P_next = P0
P_prev =  P0
Q_in = 0.00001  

kb=1.38e-23
t=0
tmax=1
dt=0.001
R=8.138
rho0 = P0 / (R * T0)
rho_next = rho0
rho_prev = rho0
gamma=1.3

const = P0/(rho0**gamma)

T_prev=T0
T_next=T0

P=[]
T=[]
for i in range(465):
    t = t + dt
    const_A = Q_in * rho0**2
    const_A = const_A / V0
    const_B = P0**(1/gamma)
    const_B = 1/const_B
    P_next = P_prev + dt * const_A * const_B * P_prev**(1.76)
    #P_next = P_prev + dt * Q_in / V0 * P0**(-1/gamma) * P_prev**((1-gamma)/gamma)
    P_prev = P_next 
    T.append( ( (P0**(gamma-1)*T0**(-gamma) ) / (P_next**(gamma-1)) )**(1/(-gamma)) )
    P.append(P_next)
    
plt.figure()   
plt.semilogy(P)
plt.title(f'presssure evolution in filling of gas of a tank of volume ={V0} m^3 with Qin={Q_in} m^3/s \n')
plt.xlabel(r'$time$ (ms)')
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()   
plt.semilogy(T)
plt.title(f'temperature evolution in filling of gas of a tank of volume ={V0} m^3 with Qin={Q_in} m^3/s \n')
plt.xlabel(r'$time$ (ms)')
plt.ylabel(r'$Temperature$ (K)')

# remark : almost incompressible fluid
# the pressure rise fast to infinity as it is exponential like 
# maybe the inflow is not realistic compared to the volume
# i do have first V0= 0.0063 m^3 for Qin=1^-5 m^3/s
