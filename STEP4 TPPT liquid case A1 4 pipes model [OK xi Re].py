# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 16:08:26 2023

@author: gauth
"""

"""
### Comment on the code : 
    CODE OK with NEW Reynolds implementation and new xi 
    4 pipes implemented with chiara equations 
    eigen values 
    speed of sound not corrected implemented but with an option to correct it 
"""

import scipy as scipy
import numpy as np
import math
import matplotlib.pyplot as plt
import sys

###################
# geometric input 
###################

l_pipe_1 = 1
l_pipe_2 = 1
l_pipe_3 = 1
l_pipe_4 = 1
l_tank_1 = 10
l_tank_2 = 10

d_tank_1 = 1
d_tank_2 = 1
d_pipe_1 = 0.5
d_pipe_2 = 0.5
d_pipe_3 = 0.5
d_pipe_4 = 0.5

A_pipe_1 = (d_pipe_1/2)**2 * 3.14
A_pipe_2 = (d_pipe_2/2)**2 * 3.14
A_pipe_3 = (d_pipe_3/2)**2 * 3.14
A_pipe_4 = (d_pipe_4/2)**2 * 3.14
A_tank_1 = (d_tank_1/2)**2 * 3.14
A_tank_2 = (d_tank_2/2)**2 * 3.14

V_tank_1 = l_tank_1 * A_tank_1
V_tank_2 = l_tank_2 * A_tank_2
V_pipe_1 = l_pipe_1 * A_pipe_1
V_pipe_2 = l_pipe_2 * A_pipe_2
V_pipe_3 = l_pipe_3 * A_pipe_3
V_pipe_4 = l_pipe_4 * A_pipe_4

###################
# initial input 
###################

P0i = 0.6e6
P01i = 0.4e6
P1i = 0.4e6
P12i = 0.4e6
P2i = 0.5e6

m1i = 0
m2i = 0
m3i = 0
m4i = 0

t = 0
dt = 0.0000005
# dt = 1e-3 converge (vir la différence entre deux division du pas si elle est minime le resultat a converge)

m1_prev = m1i
m2_prev= m2i
m3_prev = m3i
m4_prev = m4i


rho0i = 1197
rho1i = 1197
rho2i = 1197

type_simulation = 1
# 0 : gas & 1 : liquid

P0 = P0i + rho0i*10*l_tank_1
P01 = P01i
P1 = P1i
P12 = P12i
P2 = P2i + rho2i*10*l_tank_2


###################
# constant input 
###################

gamma = 1
mu = 1e-5
eta = 1e-3

zeta1 = 0.5625
zeta2 = 0
zeta3 = 0
zeta4 = 9

list_Pressure_0 = []
list_Pressure_01 = []
list_Pressure_1 = []
list_Pressure_12 = []
list_Pressure_2 = []

list_mass_flow_1 = []
list_mass_flow_2 = []
list_mass_flow_3 = []
list_mass_flow_4 = []

time = []


def a_sound(Pressure_current,Pressure_init,rho_init):
    # if liquid, option to return something cste
    if type_simulation == 0 :
        return gamma * Pressure_current  / (rho_init * (Pressure_current/Pressure_init)**(1/gamma) )
    else : 
        return ( 2.2e9 / rho_init ) / (1+(d_pipe_1/0.02)*(2.2e9/50e9))
        # return ( 2.2e9 / rho_init )

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
    # return m / (length*mu)
    return (m * length ) / ( section * eta )
    # return 300

# while np.abs(P0-P2)>5e-5 :
while t<0.1:
    t = t + dt
    # sys.stdout.write("\rSimulation time advancement: {0:.2f} %".format(100* t/3))
    # sys.stdout.flush()
    # loop to calculus
     
    P0n = P0 + dt * ( 9.81 / A_tank_1 ) * (- m1_prev)
    P01n = P01 + dt * ( (m1_prev-m2_prev) / (V_pipe_1) * a_sound(P01,P1i,rho1i) )
    P1n = P1 + dt * ( -(m3_prev-m2_prev) / (V_pipe_2) * a_sound(P1,P1i,rho1i) )
    P12n = P12 + dt * ( -(m4_prev-m3_prev) / (V_pipe_3) * a_sound(P12,P1i,rho1i) )
    P2n = P2 + dt * ( ( (m4_prev) * 9.81 / A_tank_2)  )

    lambda1 = 64/Reynolds(m1_prev,l_pipe_1,A_pipe_1)*1/d_pipe_1 *l_pipe_1
    xi1 = 0.5 * A_pipe_1**(-2) * ( lambda1 + zeta1 )
    m1_next = m1_prev + dt * ( (A_pipe_1/l_pipe_1) * ( P0 - P01 - ( xi1 / rho(P0,P0i,rho1i) ) * m1_prev * np.abs(m1_prev) ) ) 
    
    lambda2 = 64/Reynolds(m2_prev,l_pipe_2,A_pipe_2)*2/d_pipe_2 *l_pipe_2
    xi2 = 0.5 * A_pipe_2**(-2) * ( lambda2 + zeta2 )
    m2_next = m2_prev + dt * ( (A_pipe_2/l_pipe_2) * ( P01 - P1 - ( xi2 / rho(P01,P01i,rho1i) ) * m2_prev * np.abs(m2_prev) ) ) 
    
    lambda3 = 64/Reynolds(m3_prev,l_pipe_3,A_pipe_3)*3/d_pipe_3 *l_pipe_3
    xi3 = 0.5 * A_pipe_3**(-2) * ( lambda3 + zeta3 )
    m3_next = m3_prev + dt * ( (A_pipe_3/l_pipe_3) * ( P1 - P12 - ( xi3 / rho(P1,P1i,rho1i) ) * m3_prev * np.abs(m3_prev) ) ) 
    
    lambda4 = 64/Reynolds(m4_prev,l_pipe_4,A_pipe_4)*4/d_pipe_4 *l_pipe_4
    xi4 = 0.5 * A_pipe_4**(-2) * ( lambda4 + zeta4 )
    m4_next = m4_prev + dt * ( (A_pipe_4/l_pipe_4) * ( P12 - P2 - ( xi4 / rho(P12,P12i,rho1i) ) * m4_prev* np.abs(m4_prev) ) ) 

    

    # loop to integrate 
    P0 = P0n
    P01 = P01n
    P1 = P1n
    P12 = P12n
    P2 = P2n
    m1_prev = m1_next
    m2_prev= m2_next
    m3_prev= m3_next
    m4_prev= m4_next
    # m0_prev = m0_next
    # m3_prev = m3_next
    
    # loop to save 
    list_Pressure_0.append(P0n)
    list_Pressure_1.append(P1n)
    list_Pressure_2.append(P2n)
    list_Pressure_01.append(P01n)
    list_Pressure_12.append(P12n)
    list_mass_flow_1.append(m1_prev)
    list_mass_flow_2.append(m2_prev)
    list_mass_flow_3.append(m3_prev)
    list_mass_flow_4.append(m4_prev)
    time.append(t)
    if np.abs(t-0.1)<dt:
        aa = np.sqrt(a_sound(P1,P1i,rho1i))
        
        print(f'time : {t} \n P0 : {P0}\n P1 : {P1} \n P2 : {P2}\n m1 : {m1_prev} \n m2 : {m2_prev}\n m3 : {m3_prev} \n m4 : {m4_prev} \n a : {aa}')

        aa = np.sqrt(a_sound(P1,P1i,rho1i))
        print(f'a speed sound not corrected : {aa*np.sqrt(1+(l_pipe_1/0.1)*(2.2e9/50e9))}')
        print(f'a speed sound corrected : {aa}')
        V1=V_pipe_1
        lin1 = A_pipe_1/l_pipe_1*xi1/rho1i*np.abs(m1_prev)
        lin2 = A_pipe_2/l_pipe_2*xi2/rho1i*np.abs(m2_prev)
        lin3 = A_pipe_1/l_pipe_1*xi3/rho1i*np.abs(m3_prev)
        lin4 = A_pipe_1/l_pipe_1*xi4/rho1i*np.abs(m4_prev)
        A = np.array([ [0,0,0,0,0,a_sound(P01,P1i,rho1i)/V1,-a_sound(P01,P1i,rho1i)/V1,0,0],
                     [0,0,0,0,0,0,0,a_sound(P12,P1i,rho1i)/V1,-a_sound(P12,P1i,rho1i)/V1],
                     [0,0,0,0,0,0,a_sound(P1,P1i,rho1i)/V1,-a_sound(P1,P1i,rho1i)/V1,0],
                     [0,0,0,0,0,0,0,0,9.81/A_tank_2],
                     [0,0,0,0,0,-9.81/A_tank_1,0,0,0],
                     [-A_pipe_1/l_pipe_1,0,0,0,A_pipe_1/l_pipe_1,-lin1,0,0,0],
                     [A_pipe_1/l_pipe_1,0,-A_pipe_1/l_pipe_1,0,0,0,-lin2,0,0],
                     [0,-A_pipe_1/l_pipe_1,A_pipe_1/l_pipe_1,0,0,0,0,-lin3,0],
                     [0,A_pipe_1/l_pipe_1,0,-A_pipe_1/l_pipe_1,0,0,0,0,-lin4] ])
        vap = np.linalg.eigvals(A)
        print(f'eigen values are : {vap}')
        freqq = np.abs(np.imag(vap[5])/(2*3.14))
        print(f'eigen frequencies are : {freqq} Hz')
        print(f'harmonic frequencies : {aa/(8*l_pipe_1)}')


plt.figure()
plt.plot(time[::100],list_Pressure_0[::100])
plt.plot(time[::100],list_Pressure_01[::100])
plt.plot(time[::100],list_Pressure_1[::100])
plt.plot(time[::100],list_Pressure_12[::100])
plt.plot(time[::100],list_Pressure_2[::100])
plt.title(f'Liquid Oxygen model : Tank Pipe Tank \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['Tank n°1 pressure P0','Pipe pressure P01','Pipe pressure P1', 'Pipe pressure P12',' Tank n°2 pressure P2'])
plt.ylabel(r'$Pressure$ (Pa)')

plt.figure()
plt.plot(time[::1000],list_mass_flow_1[::1000])
plt.plot(time[::1000],list_mass_flow_2[::1000])
plt.plot(time[::1000],list_mass_flow_3[::1000])
plt.plot(time[::1000],list_mass_flow_4[::1000])

plt.title(f'mass evolution in filling of Liquid Oxygen of tanks \n')
plt.xlabel(r'$time$ (s)')
plt.legend(['m1','m2','m3','m4'])
plt.ylabel(r'$mass flow$ (kg/s)')



from numpy.fft import fft, fftfreq
plt.figure()
# Calcul FFT
x=list_Pressure_1
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
plt.ylim(0,500000)
plt.grid()
plt.xlabel(r"Fréquence (Hz)")
plt.ylabel(r"Amplitude $|X(f)|$")
plt.title("Transformée de Fourier")
plt.show()

# Get the indices of maximum element in numpy array
result = np.where(X_abs == np.amax(X_abs))

