# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 11:17:05 2023

@author: gauth
"""
import scipy as scipy
import numpy as np
import math
import matplotlib.pyplot as plt
import sys

###################
# geometric input 
###################

l_pipe_1 = 2
l_pipe_2 = 2
l_tank_1 = 10
l_tank_2 = 10

d_tank_1 = 1
d_tank_2 = 1
d_pipe_1 = 0.5
d_pipe_2 = 0.5

A_pipe_1 = (d_pipe_1/2)**2 * 3.14
A_pipe_2 = (d_pipe_2/2)**2 * 3.14
A_tank_1 = (d_tank_1/2)**2 * 3.14
A_tank_2 = (d_tank_2/2)**2 * 3.14

V_tank_1 = l_tank_1 * A_tank_1
V_tank_2 = l_tank_2 * A_tank_2
V_pipe_1 = l_pipe_1 * A_pipe_1
V_pipe_2 = l_pipe_2 * A_pipe_2

###################
# initial input 
###################

P0i = 0.6e6
P1i = 0.4e6
P2i = 0.5e6

m1i = 0
m2i = 0

t = 0
dt = 0.00005
# dt = 1e-3 converge (vir la différence entre deux division du pas si elle est minime le resultat a converge)
P0 = P0i
P1 = P1i
P2 = P2i

m1_prev = m1i
m2 = m2i

rho0i = 50
rho1i = 50
rho2i = 50

###################
# constant input 
###################

gamma = 1.4
mu = 1e-5
zeta1 = 0.5625
zeta2 = 9

P = []
Pr = []
Pre = []
Q = []
Qr = []

time = []
type_simulation = 0
# 0 : gas & 1 : liquid

def a_sound(Pressure_current,Pressure_init,rho_init):
    # if liquid, option to return something cste
    if type_simulation == 0 :
        return gamma * Pressure_current  / (rho_init * (Pressure_current/Pressure_init)**(1/gamma) )
    else : 
        return 2e9 * Pressure_current/rho_init

def rho(Pressure_current,Pressure_init,rho_init):
    # if liquid, option to return something cste
    if type_simulation == 0 : 
        return rho_init * (Pressure_current/Pressure_init)**(1/gamma)
    else :
        return rho_init

def Reynolds(m_flow, length):
    m=m_flow
    if m_flow==0:
        m=1e-12
    return m / (length*mu)
    # return 300

# while np.abs(P0-P2)>5e-5 :
while t<2:
    t = t + dt
    # sys.stdout.write("\rSimulation time advancement: {0:.2f} %".format(100* t/3))
    # sys.stdout.flush()
    # loop to calculus
    P0n = P0 + dt * ( a_sound(P0,P0i,rho0i) / V_tank_1 ) * (-1 * m1_prev)
    # P0n = P0 + dt * (-10)/A_tank_1 * m1_prev
    P1n = P1 + dt * ( (m1_prev-m2) / (V_pipe_1 + V_pipe_2) * a_sound(P1,P1i,rho1i) )
    P2n = P2 + dt * ( ( m2 / V_tank_2) * a_sound(P2,P2i,rho2i) )
    
    xi1 = 64/Reynolds(m1_prev,l_pipe_1)*1/d_pipe_1 * 0.5 * A_pipe_1**(-2) + zeta1 * 0.5 * A_pipe_1**(-2)
    m1_next = m1_prev + dt * ( 1/l_pipe_1 * ( P0*A_pipe_1 - P1 * A_pipe_1 - xi1 / rho(P1,P1i,rho1i) * A_pipe_1 * m1_prev * np.abs(m1_prev) ) ) 
    # m1_next = m1_prev + dt * A_pipe_1/l_pipe_1*(P0 - P1 - xi1/rho1i * m1_prev * np.abs(m1_prev))
    
    xi2 = 64/Reynolds(m2,l_pipe_2)*1/d_pipe_2 * 0.5 * A_pipe_2**(-2) + zeta2 * 0.5 * A_pipe_2**(-2)
    m2n = m2 + dt * ( 1/l_pipe_2 * ( P1*A_pipe_1 - P2 * A_pipe_1 - xi2 / rho(P1,P1i,rho1i) * A_pipe_2 * m2 * np.abs(m2) ) ) 
    # m2n = m2 + dt * A_pipe_2/l_pipe_2*(P1 - P2 - xi2/rho2i * m2 * np.abs(m2))

    # loop to integrate 
    P0 = P0n
    P1 = P1n
    P2 = P2n
    m1_prev = m1_next
    m2 = m2n
    
    # loop to save 
    P.append(P0n)
    Pr.append(P1n)
    Pre.append(P2n)
    Q.append(m1_prev)
    Qr.append(m2)
    time.append(t)
    
plt.figure()
plt.plot(time,P)
plt.plot(time,Pr)
plt.plot(time,Pre)
plt.title(f'Gas Oxygen model : Tank Pipe Tank \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['P0','P1', 'P2'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time,Q)
plt.plot(time,Qr)
plt.title(f'mass evolution in filling of Gas Oxygen of tanks \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['m1_prev','m2'])
plt.ylabel(r'$mass flow$ (kg/s)')