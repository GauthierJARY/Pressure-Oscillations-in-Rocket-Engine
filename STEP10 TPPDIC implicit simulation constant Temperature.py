# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 11:58:33 2023

@author: gauth
"""

# =============================================================================
# Imported modules 
# =============================================================================

import scipy as scipy
import numpy as np
import math
import matplotlib.pyplot as plt
import sys
from scipy.optimize import fsolve
import time 

# =============================================================================
# Geometric parametrisation 
# =============================================================================

l_tank = 10
d_tank = 1
A_tank = (d_tank/2)**2 * 3.14
V_tank = l_tank * A_tank

l_pipe = 2
d_pipe = 0.01
A_pipe = (d_pipe/2)**2 * 3.14
V_pipe = l_pipe * A_pipe

l_dome = 0.1
d_dome = 0.15
A_dome = (d_dome/2)**2 * 3.14
V_dome = l_dome * A_dome

NI = 5 # number of injectors at the end of injection dome of each species

l_injector = 0.05
d_injector = 0.003
A_injector = (d_injector/2)**2 * 3.14
V_injector = l_injector * A_injector

l_combustor = 0.2
d_combustor = 0.08 
A_combustor = (d_combustor/2)**2 * 3.14
V_combustor =  A_combustor * l_combustor 
A_throat = A_combustor/8

# =============================================================================
# Constant values and parameters
# =============================================================================

R = 8.734
R_cc = R/(29e-3)

T3 = 3500 #K

t = 0
dt = 1e-4
time_of_simulation = 0.3

gamma = 1.3

# dynamic viscosity
eta_oxygene = 42e-6
eta_propane = 117e-6 # Pa/s

zeta1 = 0
zeta2 = (1-d_pipe/d_dome)**2
zeta4 = 0
zeta5 = (1-d_pipe/d_dome)**2
zeta6 = 30
zeta7 = 2.85 #(1-d_injector/d_combustor)**2 #2.85
zeta8 = 30
zeta9 = 2.85 #2.85

density_oxygene = 1141
density_fuel = 500

bulk_modulus_oxygene = 0.2e9
bulk_modulus_fuel = 0.36e9

material_pipe = 200e9 #Pa , Young Modulus steel
material_dome = 200e9 #Pa , steel
material_injector = 200e9 #Pa , steel
material_combustor = 117e9 #Pa , copper


g = 9.81 # gravity constant 
# =============================================================================
# Initial values of my variables
# =============================================================================

P0i = 40e5 # Pa
P1i = 0
P2i = 0
P3i = 40e5
P4i = 0
P5i = 0
P6i = 0
P7i = 0
P8i = 0
P9i = 0
m0i = 0 # kg/s
m1i = 0 # kg/s
m2i = 0
m3i = 0
m4i = 0
m5i = 0
m6i = 0
m7i = 0
m8i = 0
m9i = 0
m10i = 0
T10i = 300 # K
# =============================================================================
# Initialisation of my variables
# =============================================================================

P0 = P0i + density_oxygene*g*l_tank
P1 = P1i
P2 = P2i
P3 = P3i + density_fuel*g*l_tank
P4 = P4i
P5 = P5i
P6 = P6i
P7 = P7i
P8 = P8i
P9 = P9i
m0 = m0i
m1 = m1i
m2 = m2i
m3 = m3i
m4 = m4i
m5 = m5i
m6 = m6i
m7 = m7i
m8 = m8i
m9 = m9i
m10 = m10i
T10 = T10i

# =============================================================================
# Few functions for calculus
# =============================================================================

def a_sound(density, bulk_modulus,length, material):
    # if liquid, option to return something cste
    return ( bulk_modulus / density ) / (1+(length/0.1)*(bulk_modulus/material)) # corrected speed of sound
    # return  bulk_modulus / density

def Reynolds(m_flow, length, section, eta):
    m=m_flow
    if m_flow==0:
        m=1e-12
    return (m * length ) / ( section * eta )

def trouver_maximum_fft(array,nombre_max,array_freq, length, material):
    frequencies_vector = []
    for i in range(nombre_max):
        indice_max = np.argmax(array)
        max_value = array[indice_max]
        # print("La valeur maximale est :", max_value)
        # print("Son indice dans la liste est :", indice_max)
        # print("frequence associée:", array_freq[indice_max])
        scatter = [array_freq[indice_max],max_value]
        frequencies_vector.append(scatter)
        array[indice_max] = 0 # on le met à 0 et on réitère
    freq_harmonique_theorique = np.sqrt(a_sound(density_oxygene, bulk_modulus_oxygene, length,material))/(2*length)
    if length == l_pipe : 
        freq_harmonique_theorique=freq_harmonique_theorique/2
    # print(f"freq harmonique theorique: {freq_harmonique_theorique}\n")
    return frequencies_vector
        

def calculate_fft_and_plot(x, dt, title, length, material, density, bulk, zero_padding_factor=2):
    N = len(x)
    N_padded = N * zero_padding_factor  # Nouvelle taille après zéro padding
    X = np.fft.fft(x, n=N_padded)  # Transformée de Fourier avec zéro padding
    freq = np.fft.fftfreq(N_padded, dt)  # Fréquences de la transformée de Fourier avec zéro padding
    X_abs = np.abs(X[:N_padded//2]) * 2.0 / N_padded  # Normalisation et garde uniquement les fréquences positives
    freq_pos = freq[:N_padded//2]

    plt.figure()
    plt.semilogx(freq_pos, X_abs, label="Amplitude absolue")
    plt.ylim(0, 1e6)
    plt.grid()
    plt.xlabel(r"Fréquence (Hz)")
    plt.ylabel(r"Amplitude $|X(f)|$")
    plt.title(title)
    vec = trouver_maximum_fft(X_abs, 3, freq_pos, length, material)
    for i in vec:
        marker = np.random.choice(['o', 's', '^', 'D', 'v', 'p', '*', 'x', '+'])  # Choix aléatoire d'un symbole
        color = np.random.choice(["#FFA500",'r',"#FF69B4"])
        plt.scatter(i[0], i[1], color=color, marker=marker, label=f'max frequencies : {i[0]:.0f} Hz')
    freq_harmonique_theorique = np.sqrt(a_sound(density, bulk, length, material))/(2*length)
    plt.axvline(freq_harmonique_theorique, color='g', linestyle='--', label=f'theoric frequencies : {freq_harmonique_theorique:.0f} Hz')
    # plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1))
    plt.legend()
    plt.show()
    # print(f"For {title.lower()}, we have:")
    # trouver_maximum_fft(X_abs, 4, freq_pos, length, material)
    
    
# Function to solve implicit problem 
def function(S):
    # /!\ to properly say it, xi is defined by the time step just before 
    # (as we solve implicit we cannot add built in function in the solver)
    # we consider the time step to be sufficiently small that :
    #    P(n+1) close to P(n) HYPOTHESIS
    
    # /!\ in the report, we called zeta_geo and zeta_fluid, but here we have zeta_geo and lambda = zeta_fluid 
    # their sums is equal to xi, as defined in the report
    # there is no other impacts
    lambda1 = 64/Reynolds(m1,l_pipe,A_pipe, eta_oxygene)*1/d_pipe *l_pipe
    xi1 = 0.5 * A_pipe**(-2) * ( lambda1 + zeta1 )
    
    lambda2 = 64/Reynolds(m2,l_pipe,A_pipe, eta_oxygene)*1/d_pipe *l_pipe
    xi2 = 0.5 * A_pipe**(-2) * ( lambda2 + zeta2 )
    
    lambda4 = 64/Reynolds(m4,l_pipe,A_pipe, eta_propane)*1/d_pipe *l_pipe
    xi4 = 0.5 * A_pipe**(-2) * ( lambda4 + zeta4 )
    
    lambda5 = 64/Reynolds(m5,l_pipe,A_pipe, eta_propane)*1/d_pipe *l_pipe
    xi5 = 0.5 * A_pipe**(-2) * ( lambda5 + zeta5 )
    
    lambda6 = 64/Reynolds(m6,l_dome,A_dome, eta_oxygene)*1/d_dome *l_dome
    xi6 = 0.5 * A_dome**(-2) * ( lambda6 + zeta6 )
    
    lambda8 = 64/Reynolds(m8,l_dome,A_dome, eta_propane)*1/d_dome *l_dome
    xi8 = 0.5 * A_dome**(-2) * ( lambda8 + zeta8 )
    
    lambda7 = 64/Reynolds(m7,l_injector,A_injector, eta_oxygene)*1/d_injector *l_injector
    xi7 = 0.5 * A_injector**(-2) * ( lambda7 + zeta7 )
    
    lambda9 = 64/Reynolds(m9,l_injector,A_injector, eta_propane)*1/d_injector *l_injector
    xi9 = 0.5 * A_injector**(-2) * ( lambda9 + zeta9 )
    
    return np.array([
        (S[0]-P0)/dt + S[10]*g/A_tank ,
        (S[3]-P3)/dt + S[13]*g/A_tank , 
        (S[0]-P0)/dt - a_sound(density_oxygene, bulk_modulus_oxygene,l_pipe, material_pipe) / V_pipe * (S[10] - S[11]) ,
        (S[1]-P1)/dt - a_sound(density_oxygene, bulk_modulus_oxygene,l_pipe, material_pipe) / V_pipe * (S[11] - S[12]) ,
        (S[3]-P3)/dt - a_sound(density_fuel, bulk_modulus_fuel,l_pipe, material_pipe) / V_pipe * (S[13] - S[14]) ,
        (S[4]-P4)/dt - a_sound(density_fuel, bulk_modulus_fuel,l_pipe, material_pipe) / V_pipe * (S[14] - S[15]) ,
        (S[11] - m1)/dt - A_pipe/l_pipe * (S[0] - S[1] - xi1 / density_oxygene * S[11]*np.abs(S[11])) ,
        (S[12] - m2)/dt - A_pipe/l_pipe * (S[1] - S[2] - xi2 / density_oxygene * S[12]*np.abs(S[12])) ,
        (S[14] - m4)/dt - A_pipe/l_pipe * (S[3] - S[4] - xi4 / density_fuel * S[14]*np.abs(S[14])) ,
        (S[15] - m5)/dt - A_pipe/l_pipe * (S[4] - S[5] - xi5 / density_fuel * S[15]*np.abs(S[15])) ,
        (S[2]-P2)/dt - a_sound(density_oxygene, bulk_modulus_oxygene,l_dome, material_dome) / V_dome * (S[12] - S[16]) ,
        (S[16] - m6)/dt - A_dome/l_dome * (S[2] - S[6] - xi6 / density_oxygene * S[16]*np.abs(S[16])) ,
        (S[6]-P6)/dt - a_sound(density_oxygene, bulk_modulus_oxygene,l_injector, material_injector) / V_injector * (S[16]/NI - S[17]/NI) ,
        # (S[17] - m7)/dt - NI*A_injector/l_injector * ( S[6] - S[7] - xi7 / density_oxygene * S[17]*np.abs(S[17]) ) ,
        (S[17]/NI - m7/NI)/dt - A_injector/l_injector * ( S[6] - S[7] - xi7 / density_oxygene * S[17]*np.abs(S[17])/NI**2 ) ,
        (S[5]-P5)/dt - a_sound(density_fuel, bulk_modulus_fuel,l_dome, material_dome) / V_dome * (S[15] - S[18]) ,
        (S[18] - m8)/dt - A_dome/l_dome * (S[5] - S[8] - xi8 / density_fuel * S[18]*np.abs(S[18])) ,
        (S[8]-P8)/dt - a_sound(density_fuel, bulk_modulus_fuel,l_injector, material_injector) / V_injector * (S[18]/NI - S[19]/NI) ,
        # (S[19] - m9)/dt - NI*A_injector/l_injector * ( S[8] - S[9] - xi9 / density_oxygene * S[19]*np.abs(S[19]) ) ,
        (S[19]/NI - m9/NI)/dt - A_injector/l_injector * ( S[8] - S[9] - xi9 / density_oxygene * S[19]*np.abs(S[19])/NI**2 ) ,
        S[7] - S[9] ,
        S[20] - ( np.sqrt(gamma * (2/(gamma+1) )**((gamma+1)/(gamma-1)) ) * A_throat )/(np.sqrt(R_cc * S[21])) * S[7] ,
        # V_combustor/( R_cc * S[21] ) * (S[7] - P7)/dt + S[20] - S[19] - S[17] - 1*np.sin(95*2*np.pi*t),
        (S[7] - P7)/dt + R_cc * S[21]/V_combustor *( S[20] - S[19] - S[17] ),
        # S[7] - (538965 *( 1 + 0.15*np.sin(1043*2*np.pi*t) ) ),
        S[21] - 3000
        ])

# =============================================================================
# Plots lists
# =============================================================================
list_pressure_0 = []
list_pressure_1 = []
list_pressure_2 = []
list_pressure_3 = []
list_pressure_4 = []
list_pressure_5 = []
list_pressure_6 = []
list_pressure_7 = []
list_pressure_8 = []
list_pressure_9 = []

list_mass_flow_0 = []
list_mass_flow_1 = []
list_mass_flow_2 = []
list_mass_flow_3 = []
list_mass_flow_4 = []
list_mass_flow_5 = []
list_mass_flow_6 = []
list_mass_flow_7 = []
list_mass_flow_8 = []
list_mass_flow_9 = []
list_mass_flow_10 = []

list_temperature_10 = []

time_storage = []

start = time.time() 
print("Start computing") 

while t<time_of_simulation: 
    
    ############################################################  
    #### LOOP AND CALCULUS 
    ############################################################
    
    t = t + dt
    S_vector = [ P0, P1, P2, P3, P4, P5, P6, P7, P8, P9,
                m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10,
                T10
                ]
    root = fsolve(function, S_vector)
    P0 = root[0]
    P1 = root[1]
    P2 = root[2]
    P3 = root[3]
    P4 = root[4]
    P5 = root[5]
    P6 = root[6]
    P7 = root[7]
    P8 = root[8]
    P9 = root[9]
    m0 = root[10]
    m1 = root[11]
    m2 = root[12]
    m3 = root[13]
    m4 = root[14]
    m5 = root[15]
    m6 = root[16]
    m7 = root[17]
    m8 = root[18]
    m9 = root[19]
    m10 = root[20]
    T10 = root[21]
    
    ############################################################  
    #### eigen values with matrix
    ############################################################
    if np.abs(t-time_of_simulation)<=0.00001 :
        lambda1 = 64/Reynolds(m1,l_pipe,A_pipe, eta_oxygene)*1/d_pipe *l_pipe
        xi1 = 0.5 * A_pipe**(-2) * ( lambda1 + zeta1 )
        lambda2 = 64/Reynolds(m2,l_pipe,A_pipe, eta_oxygene)*1/d_pipe *l_pipe
        xi2 = 0.5 * A_pipe**(-2) * ( lambda2 + zeta2 )
        lambda4 = 64/Reynolds(m4,l_pipe,A_pipe, eta_propane)*1/d_pipe *l_pipe
        xi4 = 0.5 * A_pipe**(-2) * ( lambda4 + zeta4 )
        lambda5 = 64/Reynolds(m5,l_pipe,A_pipe, eta_propane)*1/d_pipe *l_pipe
        xi5 = 0.5 * A_pipe**(-2) * ( lambda5 + zeta5 )
        lambda6 = 64/Reynolds(m6,l_dome,A_dome, eta_oxygene)*1/d_dome *l_dome
        xi6 = 0.5 * A_dome**(-2) * ( lambda6 + zeta6 )
        lambda8 = 64/Reynolds(m8,l_dome,A_dome, eta_propane)*1/d_dome *l_dome
        xi8 = 0.5 * A_dome**(-2) * ( lambda8 + zeta8 )
        lambda7 = 64/Reynolds(m7,l_injector,A_injector, eta_oxygene)*1/d_injector *l_injector
        xi7 = 0.5 * A_injector**(-2) * ( lambda7 + zeta7 )
        lambda9 = 64/Reynolds(m9,l_injector,A_injector, eta_propane)*1/d_injector *l_injector
        xi9 = 0.5 * A_injector**(-2) * ( lambda9 + zeta9 )
            
        # d S / dt + A *S = b
        A = np.zeros((22,22))
        # vector size of 22x22
        # line by line
        # A[10] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        # A[10] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, o, T]
    
        A[0] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g/A_tank, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        A[1] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, a_sound(density_oxygene, bulk_modulus_oxygene,l_pipe, material_pipe) / V_pipe, -a_sound(density_oxygene, bulk_modulus_oxygene,l_pipe, material_pipe) / V_pipe, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        A[2] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, a_sound(density_oxygene, bulk_modulus_oxygene,l_dome, material_dome) / V_dome, 0, 0, 0, -a_sound(density_oxygene, bulk_modulus_oxygene,l_dome, material_dome) / V_dome, 0, 0, 0, 0, 0]
        A[3] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g/A_tank, 0, 0, 0, 0, 0, 0, 0, 0]
        A[4] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, a_sound(density_fuel, bulk_modulus_fuel,l_pipe, material_pipe) / V_pipe, -a_sound(density_fuel, bulk_modulus_fuel,l_pipe, material_pipe) / V_pipe, 0, 0, 0, 0, 0, 0]
        A[5] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, a_sound(density_fuel, bulk_modulus_fuel,l_dome, material_dome) / V_dome, 0, 0, -a_sound(density_fuel, bulk_modulus_fuel,l_dome, material_dome) / V_dome, 0, 0, 0]
        A[6] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, a_sound(density_oxygene, bulk_modulus_oxygene,l_injector, material_injector) / V_injector/NI, -a_sound(density_oxygene, bulk_modulus_oxygene,l_injector, material_injector) / V_injector/NI, 0, 0, 0, 0]
        A[7] = [0, 0, 0, 0, 0, 0, 0, -R_cc*T10/V_combustor*( np.sqrt(gamma * (2/(gamma+1) )**((gamma+1)/(gamma-1)) ) * A_throat )/(np.sqrt(R_cc * T10)), 0, 0, 0, 0, 0, 0, 0, 0, 0, -R_cc*T10/V_combustor, 0, -R_cc*T10/V_combustor, 0, 0]
        A[8] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, a_sound(density_fuel, bulk_modulus_fuel,l_injector, material_injector) / V_injector/NI, -a_sound(density_fuel, bulk_modulus_fuel,l_injector, material_injector) / V_injector/NI, 0, 0]
        A[9] = [0, 0, 0, 0, 0, 0, 0, -R_cc*T10/V_combustor*( np.sqrt(gamma * (2/(gamma+1) )**((gamma+1)/(gamma-1)) ) * A_throat )/(np.sqrt(R_cc * T10)), 0, 0, 0, 0, 0, 0, 0, 0, 0, -R_cc*T10/V_combustor, 0, -R_cc*T10/V_combustor, 0, 0]
        A[10] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        A[11] = [A_pipe/l_pipe, -A_pipe/l_pipe, 0, 0, 0, 0, 0, 0, 0, 0, 0, -A_pipe/l_pipe*xi1/density_oxygene*np.abs(m1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        A[12] = [0, A_pipe/l_pipe, -A_pipe/l_pipe, 0, 0, 0, 0, 0, 0, 0, 0, -A_pipe/l_pipe*xi2/density_oxygene*np.abs(m2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        A[13] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        A[14] = [0, 0, 0, A_pipe/l_pipe, -A_pipe/l_pipe, 0, 0, 0, 0, 0, 0, 0, 0, -A_pipe/l_pipe*xi4/density_fuel*np.abs(m4), 0, 0, 0, 0, 0, 0, 0, 0]
        A[15] = [0, 0, 0, 0, A_pipe/l_pipe, -A_pipe/l_pipe, 0, 0, 0, 0, 0, 0, 0, 0, -A_pipe/l_pipe*xi5/density_fuel*np.abs(m5), 0, 0, 0, 0, 0, 0, 0]
        A[16] = [0, 0, A_dome/l_dome, 0, 0, 0, -A_dome/l_dome, 0, 0, 0, 0, 0, 0, 0, 0, 0, -xi6/density_oxygene * A_dome/l_dome*np.abs(m6), 0, 0, 0, 0, 0]
        A[17] = [0, 0, 0, 0, 0, 0, A_injector/l_injector*NI, -A_injector/l_injector*NI, 0, 0, 0, 0, 0, 0, 0, 0, 0, -xi7/density_oxygene * A_injector/l_injector/NI*np.abs(m7), 0, 0, 0, 0]
        A[18] = [0, 0, 0, 0, 0, A_dome/l_dome, 0, 0, -A_dome/l_dome, 0, 0, 0, 0, 0, 0, 0, 0, 0, -xi8/density_oxygene * A_dome/l_dome*np.abs(m8), 0, 0, 0]
        A[19] = [0, 0, 0, 0, 0, 0, 0, 0, A_injector/l_injector*NI, -A_injector/l_injector*NI, 0, 0, 0, 0, 0, 0, 0, 0, 0, -xi9/density_fuel * A_injector/l_injector/NI*np.abs(m9), 0, 0]
        A[20] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        A[21] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        # Vecteur constant b
        b = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, T3])
    
        # Afficher la matrice A
        vap = np.linalg.eigvals(-A)
        print(f'eigen values are : {vap}')
        frequencies_matrix = []
        for i in vap : 
            f = np.imag(i)/(2*3.14)
            frequencies_matrix.append(f)
            frequencies_matrix = [abs(x) for x in frequencies_matrix]
        for j in range(len(frequencies_matrix)):
            S_name = ['P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9','m0', 'm1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10','T10']
            print(f'eigen frequencies of {S_name[j]} are : {frequencies_matrix[j]:.1f} Hz')
    A0 = A_injector * NI
    discharge_coefficient = m6/density_oxygene / A0 / np.sqrt(2/density_oxygene * (np.abs(P7 - P6)))
    ############################################################  
    #### SAVE 
    ############################################################
    
    list_pressure_0.append(P0)
    list_pressure_1.append(P1)
    list_pressure_2.append(P2)
    list_pressure_3.append(P3)
    list_pressure_4.append(P4)
    list_pressure_5.append(P5)
    list_pressure_6.append(P6)
    list_pressure_7.append(P7)
    list_pressure_8.append(P8)
    list_pressure_9.append(P9)
    
    list_mass_flow_0.append(m0)
    list_mass_flow_1.append(m1)
    list_mass_flow_2.append(m2)
    list_mass_flow_3.append(m3)
    list_mass_flow_4.append(m4)
    list_mass_flow_5.append(m5)
    list_mass_flow_6.append(m6)
    list_mass_flow_7.append(m7/NI)
    list_mass_flow_8.append(m8)
    list_mass_flow_9.append(m9/NI)
    list_mass_flow_10.append(m10)
    
    list_temperature_10.append(T10)
    
    time_storage.append(t)
    
end = time.time() 
print(f'time of simulation : {end - start}')   


# =============================================================================
# Plots and graphs 
# =============================================================================

# Set the font size for the titles
plt.rcParams['axes.titlesize'] = 7
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.labelsize'] = 9
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 7
# Set the markerscale and handlelength for the legend
plt.rcParams['legend.markerscale'] = 0.7
plt.rcParams['legend.handlelength'] = 1.0
plt.rcParams['axes.labelweight'] = 'bold'


plt.figure()
plt.plot(time_storage[::10], list_pressure_0[::10], marker='o', markersize=3, markevery=22)
plt.plot(time_storage[::10], list_pressure_1[::10], marker='s', markersize=3, markevery=20)
plt.plot(time_storage[::10], list_pressure_2[::10], marker='^', markersize=4, markevery=18)
plt.plot(time_storage[::10], list_pressure_6[::10], marker='*', markersize=3, markevery=15)
plt.plot(time_storage[::10], list_pressure_7[::10], marker='x', markersize=3, markevery=25)

plt.title(f'Pressure evolution in system Tank-Pipe-Pipe-Dome-Injectors-Combustor, Oxygene Line, with T=3000 K\n')

plt.legend(['pipe1','pipe2', 'dome', 'injector', 'combustor'])
plt.xlabel(r'Time (s)')
plt.ylabel(r'Pressure (Pa)')

plt.figure()
plt.plot(time_storage[::10], list_pressure_3[::10], marker='o', markersize=3, linestyle='-', markevery=20)
plt.plot(time_storage[::10], list_pressure_4[::10], marker='s', markersize=3, linestyle='-', markevery=22)
plt.plot(time_storage[::10], list_pressure_5[::10], marker='^', markersize=4, linestyle='-', markevery=18)
plt.plot(time_storage[::10], list_pressure_8[::10], marker='*', markersize=3, linestyle='-', markevery=15)
plt.plot(time_storage[::10], list_pressure_9[::10], marker='x', markersize=3, linestyle='-', markevery=25)


plt.title(f'Pressure evolution in system Tank-Pipe-Pipe-Dome-Injectors-Combustor, Propane Line, with T=3000 K \n')
plt.legend(['pipe3','pipe4', 'dome', 'injector', 'combustor'])
plt.xlabel(r'Time (s)')
plt.ylabel(r'Pressure (Pa)')


plt.figure()
plt.plot(time_storage[::10], list_mass_flow_0[::10], marker='o', markersize=3, markevery=28)
plt.plot(time_storage[::10], list_mass_flow_1[::10], marker='s', markersize=3, markevery=15)
plt.plot(time_storage[::10], list_mass_flow_2[::10], marker='^', markersize=3, markevery=20)
plt.plot(time_storage[::10], list_mass_flow_6[::10], marker='*', markersize=3, markevery=18)
plt.plot(time_storage[::10], list_mass_flow_7[::10], marker='x', markersize=3, markevery=20)
plt.plot(time_storage[::10], list_mass_flow_10[::10], marker='D', markersize=3, markevery=22)


plt.title(f'Mass flow evolution in system Tank-Pipe-Pipe-Dome-Injectors-Combustor, Oxygene line with T=3000 K \n')
plt.legend(['tank','pipe 1','pipe 2','dome','injector','combustor out'])
plt.ylabel(r'Mass flow (kg/s)')
plt.xlabel(r'Time (s)')

plt.figure()
plt.plot(time_storage[::10], list_mass_flow_3[::10], marker='o', markersize=3, markevery=25)
plt.plot(time_storage[::10], list_mass_flow_4[::10], marker='s', markersize=3, markevery=18)
plt.plot(time_storage[::10], list_mass_flow_5[::10], marker='^', markersize=3, markevery=15)
plt.plot(time_storage[::10], list_mass_flow_8[::10], marker='*', markersize=3, markevery=22)
plt.plot(time_storage[::10], list_mass_flow_9[::10], marker='x', markersize=3, markevery=20)
plt.plot(time_storage[::10], list_mass_flow_10[::10], marker='D', markersize=3, markevery=28)


plt.title(f'Mass flow evolution in system Tank-Pipe-Pipe-Dome-Injectors-Combustor, Propane Line; with T= 3000K\n ')
plt.legend(['tank','pipe 3','pipe 4','dome','injector','combustor out'])
plt.ylabel(r'Mass flow (kg/s)')
plt.xlabel(r'Time (s)')

print(f'Drop of pressure between injection ox and combustor is equal to : {int((1 - P7/P6)*100)} %')
print(f'Drop of pressure between injection fuel and combustor is equal to : {int((1 - P9/P8)*100)} %')



# =============================================================================
# Fast Fourier Transformation for frequencies analysis 
# =============================================================================
# print("We do FFT analysis for the oxygene line at the moment \n So, bulk oxyegene et speed of sound oxygene \n")
from numpy.fft import fft, fftfreq

# Utilisation de la fonction pour les différents cas
calculate_fft_and_plot(list_pressure_1, dt, "Transformée de Fourier de Pression dans tuyau 2\n", 2*l_pipe, material_pipe,density_oxygene,bulk_modulus_oxygene)
calculate_fft_and_plot(list_pressure_2, dt, "Transformée de Fourier de Pression dans dome injection\n", l_dome, material_injector,density_oxygene,bulk_modulus_oxygene)
calculate_fft_and_plot(list_pressure_6, dt, "Transformée de Fourier de Pression dans injecteurs\n", l_injector, material_injector,density_oxygene,bulk_modulus_oxygene)
calculate_fft_and_plot(list_pressure_7, dt, "Transformée de Fourier de Pression dans chambre de combustion\n", l_combustor, material_combustor,density_fuel,bulk_modulus_fuel)
