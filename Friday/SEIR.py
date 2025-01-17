#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:32:43 2025

@author: adminaccount
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the SEIR model
def seir_model(y, t, beta, gamma_r, gamma_d):
    S, I, R, D = y
    N = S + D + I + R  # Total population

    # Differential equations
    dS_dt = -beta * S * I / N
    dI_dt = beta * S * I / N - (gamma_r+gamma_d)*I
    dR_dt = (gamma_r)*I
    dD_dt = (gamma_d)*I
    return [dS_dt, dI_dt, dR_dt, dD_dt]

def expgrowth(I0, t, r):
    dI_dt = r*I
    return dI_dt



# Parameters
N = 100_000          # Total population
sigma = 1/4     # Incubation rate (1/average incubation period)
gamma_r = 1/7       # Recovery rate (1/average infectious period)
p_Death = 0.005   #% of people who die

doubling_time = 4

#This means that the sum of the rates needs to be 1/4, and since p_death = 0.005
gamma_r = 1/(4*(1+p_Death))
gamma_d = p_Death*gamma_r 

R_0 = 2.5

beta = R_0 * (gamma_r + gamma_d) 


# Initial conditions
I0 = 1            # Initial number of infectious individuals
D0 = 0
R0 = 0            # Initial number of recovered individuals
S0 = N - I0  - R0 -D0  # Initial number of susceptible individuals

y0 = [S0, I0, R0, D0]  # Initial state vector

# Time points (in days)
t = np.linspace(0, 60,200)  # Simulate for 160 days

# Solve the ODEs
result = odeint(seir_model, y0, t, args=(beta, gamma_r, gamma_d))
S, I, R, D = result.T  # Transpose result to get individual components

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(t, I, label='Infectious', color='red')
plt.plot(t, D, label='Dead', color='black')

plt.title('SEIR Model')
plt.xlabel('Days')
plt.ylabel('Population')
plt.legend()
plt.grid(True)
plt.show()
