# -*- coding: utf-8 -*-
"""
Created on Tue May 30 11:19:52 2023

@author: gauth
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import sys


t = 0 
dt = 0.000001
P0_prev = 0.6e6
P1_prev = 0.4e6
P2_prev = 0.5e6
beta0 = 50000

A0 = 0.5**2 * 3.14
rho1 = 1100
m_flow1_prev = 0
m_flow2_prev = 0
A1 = 0.25**2 * 3.14
l = 1
xi = 0.32 / (2 * A1**2)
xi=0
V_tank0 = A0 * 10
P111 = []
P222 = []
P2222 = []
N = 500000
time_simu = N * dt
while t<15:
    t = t + dt 
    # sys.stdout.write("\rSimulation time advancement: {0:.2f} %".format(100* t/time_simu))
    # sys.stdout.flush()
    P0_next = P0_prev + dt * ( -beta0 / (V_tank0 * rho1) * m_flow1_prev )
    P2_next = P2_prev + dt * ( beta0 / (V_tank0 * rho1) * m_flow2_prev )
    # P0_next = P0_prev + dt * (-10/A1 * m_flow1_prev )
    # P2_next = P2_prev + dt * (10/A1 * m_flow2_prev )

    m_flow1_next = m_flow1_prev + dt * ( A1/l * (P0_prev - P1_prev - xi / rho1 * m_flow1_prev * np.abs(m_flow1_prev)) )
    m_flow2_next = m_flow2_prev + dt * ( A1/l * (P1_prev - P2_prev - xi / rho1 * m_flow2_prev * np.abs(m_flow2_prev)) )
    
    P1_next = P1_prev +  dt * (1/V_tank0) * ( beta0 * 0.4e6 / rho1 * (m_flow1_prev - m_flow2_prev))
    
    P111.append(P0_next)
    P222.append(P2_next)
    P2222.append(P1_next)
    
    P0_prev = P0_next
    m_flow2_prev = m_flow2_next
    m_flow1_prev = m_flow1_next 
    P2_prev = P2_next
    P1_prev = P1_next
    
    
plt.figure()   
plt.plot(P111)
plt.plot(P222)
# plt.plot(P2222[100000:])
plt.title(f'presssure evolution in filling of LOx of a tank of volume = {V_tank0} m^3 with P_in= X Pa \n')
plt.xlabel(r'$time$ ($10e-6$ s)')
plt.legend(['tank1','tank2', 'pipe'])
plt.ylabel(r'$Pressure$ (Pa)')

