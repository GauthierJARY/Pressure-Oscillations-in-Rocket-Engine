# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 14:31:27 2023

@author: gauth

"""

"""
Comment : 
    Modele with equation for the injector 0D => 1 equation without capacitance and inductance hypothesis
    full explicit 
    done with the idea that the intjector was just a connection 
    does not converge not usefullfor the results in itself : 
        but interesting to compare with implicit 
        for the idea to neglect capacitance and inductance of the injector
"""

import scipy as scipy
import numpy as np
import math
import matplotlib.pyplot as plt
import sys

###################
# geometric input 
###################
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

###################
# initial input 
###################
P0i = 0.6e6
P1i = 0.4e6
P2i = 0.4e6
P3i = 1e5

m0i = 0
m1i = 0
m2i = 0
m3i = 0

t = 0
dt = 0.0000001
# dt = 1e-3 converge (vir la différence entre deux division du pas si elle est minime le resultat a converge)

m0_prev = m0i
m1_prev = m1i
m2_prev = m2i
m3_prev = m3i

rho0i = 1197
rho1i = 1197
rho2i = 1197
rho3i = 1197

P0 = P0i
P1 = P1i
P2 = P2i
P3 = P3i

type_simulation = 1
# 0 : gas & 1 : liquid


###################
# constant input 
###################

gamma = 1
mu = 1e-5
eta = 1e-3

zeta0 = 1 + (1-d_pipe/d_dome)
zeta1 = 0.36
zeta2 = 1 + 2.85

list_pressure_0 = []
list_pressure_1 = []
list_pressure_2 = []
list_pressure_3 = []

list_mass_flow_0 = []
list_mass_flow_1 = []
list_mass_flow_2 = []
list_mass_flow_3 = []

time = []


def a_sound(Pressure_current,Pressure_init,rho_init):
    # if liquid, option to return something cste
    if type_simulation == 0 :
        return gamma * Pressure_current  / (rho_init * (Pressure_current/Pressure_init)**(1/gamma) )
    else : 
        # return ( 2.2e9 / rho_init ) / (1+(l_pipe/0.1)*(2.2e9/50e9))
        return 2.2e9 / rho_init

def rho(Pressure_current,Pressure_init,rho_init):
    # if liquid, option to return something cste
    if type_simulation == 0 : 
        return rho_init * (Pressure_current/Pressure_init)**(1/gamma)
    else :
        return rho_init

def Reynolds(m_flow, length, section):
    m=m_flow
    if m_flow==0:
        m=1e-12
    return (m * length ) / ( section * eta ) # and not the old Re used by error : return m / (length*mu)


timecount =0
# while np.abs(P0-P2)>5e-5 :
while t<0.1:
    t = t + dt
    # timestep = 1000
    # timecount +=1
    # if timecount == timestep :
    #     sys.stdout.write("\rSimulation time advancement: {0:.3f} %\n".format(100* t/0.1))
    #     # sys.stdout.flush()
    # loop to calculus
    P0n = P0 
    P1n = P1 + dt * ( (m0_prev-m1_prev) / (V_pipe ) * a_sound(P1,P1i,rho1i) )
    P2n = P2 + dt * ( (m1_prev-m2_prev) / (V_dome ) * a_sound(P2,P2i,rho2i) )
    P3n = P3
    
    lambda0 = 64/Reynolds(m0_prev,l_pipe,A_pipe)*0/d_pipe *l_pipe
    xi0 = 0.5 * A_pipe**(-2) * ( lambda0 + zeta0 )
    m0_next = m0_prev + dt * ( (A_pipe/l_pipe) * ( P0 - P1 - ( xi0 / rho(P0,P0i,rho0i) ) * m0_prev * np.abs(m0_prev) ) ) 
    
    lambda1 = 64/Reynolds(m1_prev,l_dome,A_dome)*1/d_dome *l_dome
    xi1 = 0.5 * A_dome**(-2) * ( lambda1 + zeta1 )
    m1_next = m1_prev + dt * ( (A_pipe/l_dome) * ( P1 - P2 - ( xi1 / rho(P1,P1i,rho1i) ) * m1_prev * np.abs(m1_prev) ) ) 
    
    QINJ = A_injector * zeta2 * np.sign(P3-P2) * np.sqrt(2*rho1i*np.abs(P3-P2))
    m2_next = Nombre_injector * QINJ
    m3_next = m3_prev
    
    # loop to integrate 
    P0 = P0n
    P1 = P1n
    P2 = P2n
    P3 = P3n

    m0_prev = m0_next
    m1_prev = m1_next
    m2_prev = m2_next
    m3_prev = m3_next

    # loop to save 
    list_pressure_0.append(P0n)
    list_pressure_1.append(P1)
    list_pressure_2.append(P2n)
    list_pressure_3.append(P3n)

    list_mass_flow_0.append(m0_prev)
    list_mass_flow_1.append(m1_prev)
    list_mass_flow_2.append(m2_prev)
    list_mass_flow_3.append(m3_prev)

    time.append(t)
    
    if np.abs(t-0.1)<0.5*dt:
        aa = np.sqrt(a_sound(P1,P1i,rho1i))
        print(f'time : {t} \n P0 : {P0}\n P1 : {P1} \n P2 : {P2}\n m0 : {m0_prev} \n m1 : {m1_prev}\n m2 : {m2_prev} \n a : {aa}')
        print(f'a speed sound non corrected : {aa}')



plt.figure()
plt.plot(time[::1000],list_pressure_0[::1000])
plt.plot(time[::1000],list_pressure_1[::1000])
plt.plot(time[::1000],list_pressure_2[::1000])
plt.plot(time[::1000],list_pressure_3[::1000])

plt.title(f'Liquid Oxygen model : Pipe Dome Injector \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['P0','P1','P2','P3'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time[::1000],list_mass_flow_0[::1000])
plt.plot(time[::1000],list_mass_flow_1[::1000])
plt.plot(time[::1000],list_mass_flow_2[::1000])
plt.plot(time[::1000],list_mass_flow_3[::1000])

plt.title(f'mass evolution in Pipe Dome Injector System \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['m0','m1','m2','m3'])
plt.ylabel(r'$mass flow$ (kg/s)')



from numpy.fft import fft, fftfreq
plt.figure()
# Calcul FFT
x=list_pressure_1
N= len(x)
X = fft(x)  # Transformée de fourier
freq = fftfreq(N, dt)  # Fréquences de la transformée de Fourier
# Calcul du nombre d'échantillon
N= len(x)
# On prend la valeur absolue de l'amplitude uniquement pour les fréquences positives et normalisation
X_abs = np.abs(X[:N//2])*2.0/N
# On garde uniquement les fréquences positives
freq_pos = freq[:N//2]
m=max(X_abs)
# index=X_abs.index(m)
plt.plot(freq_pos, X_abs, label="Amplitude absolue")
plt.xlim(0, 200)  # On réduit la plage des fréquences à la zone utile
plt.grid()
plt.xlabel(r"Fréquence (Hz)")
plt.ylabel(r"Amplitude $|X(f)|$")
plt.title("Transformée de Fourier")
plt.show()

# Get the indices of maximum element in numpy array
result = np.where(X_abs == np.amax(X_abs))


