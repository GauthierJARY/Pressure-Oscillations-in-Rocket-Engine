# -*- coding: utf-8 -*-
"""
Created on Tue May 23 11:39:34 2023

@author: gauth
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import sys

########################
## Example 1 : liquid 
########################


P=[]
bulk_modulus_beta = 2e9 
Volume = 10
t=0
dt=1
tmax=60
rho_0=860
A=0.1 # section - m^-2
U0= 10 # m^3/s - inflow velocity
Qin = A * U0
Pressure_next=0
Pressure_prev = 0 
while  t<tmax:
    t+=dt
    Pressure_next = Pressure_prev + dt * Qin * bulk_modulus_beta / Volume
    P.append(Pressure_next)
    Pressure_prev=Pressure_next
plt.plot(P)
plt.title(f'presssure evolution in filling of oil of a tank with Qin={Qin}\n')
plt.xlabel(r'$time$ (s)')
plt.ylabel(r'$Pressure$ (Pa)')


########################
## Example 2 : gaz 
########################


V = 10
n = 1.3
R = 8.31
T = 2000 # Kelvin
C = V / (n*R*T)
Cd = 0.8
A0 = 0.1
Ps = 2e6
gamma = 1.3
P=[]
Pressure_prev = 1e5
Pressure_next = 1e5 # 1 bar 
t=0
tmax=60
dt=1e-3


while  t<tmax:
    t+=dt
    KK = (2*gamma)/((gamma-1)*R*T)
    if ((Pressure_prev/Ps)**(2/gamma) - (Pressure_prev/Ps)**((gamma+1)/gamma) )<0:
        break
    Pressure_next = Pressure_prev + dt * Cd * A0 * Ps * (1/C) * ( np.sqrt( KK *( (Pressure_prev/Ps)**(2/gamma) - (Pressure_prev/Ps)**((gamma+1)/gamma) ) ) )
    P.append(Pressure_next)
    Pressure_prev=Pressure_next

    
plt.figure()
plt.plot(P)
plt.title(f'presssure evolution in filling of gaz of a tank with Ps={Ps/1e5} bar\n')
plt.xlabel(r'$time$ (ms)')
plt.ylabel(r'$Pressure$ (Pa)')


########################
## Example 3 : delay gas 
########################


V = 10
n = 1.3
R = 8.31
T = 2000 # Kelvin
C = V / (n*R*T)
Cd = 0.8
A0 = 0.1
Ps = 2e6
gamma = 1.3
P=[]
Pressure_prev = 1e5
Pressure_next = 1e5 # 1 bar 
t=0
tmax=60
dt=1e-3
P_in = 2e6

list_time=np.arange(1,2500,1)
a=(P_in-1e5)/1700
def Ps(t):
    if t<200:
        return  1e5
    else : 
        return P_in

for t in list_time:
    KK = (2*gamma)/((gamma-1)*R*T)
    if ((Pressure_prev/Ps(t))**(2/gamma) - (Pressure_prev/Ps(t))**((gamma+1)/gamma) )<0:
        print('error')
        time_error = t
        break
    Pressure_next = Pressure_prev + dt * Cd * A0 * Ps(t) * (1/C) * ( np.sqrt( KK *( (Pressure_prev/Ps(t))**(2/gamma) - (Pressure_prev/Ps(t))**((gamma+1)/gamma) ) ) )
    P.append(Pressure_next)
    Pressure_prev=Pressure_next

    
plt.figure()
plt.plot(P)
plt.title(f'presssure evolution in filling of gaz of a tank with retarded Ps={P_in/1e5} bar\n')
x=[]
for i in list_time:
    x.append(Ps(i))  
plt.plot(x)
plt.vlines(time_error,0,P_in*(1+0.2),colors='r', linestyles='-',label='time error')
plt.xlabel(r'$time$ (ms)')
plt.ylabel(r'$Pressure$ (Pa)')
plt.legend(['pressure in tank (answer)', 'pressure supply (order)'])


