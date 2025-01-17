#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 11:28:47 2025

@author: adminaccount
"""

import covasim as cv
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import lognorm


# Define your custom durations, including exp2inf
durations = {
    'exp2inf': {'dist': 'lognormal', 'par1': 3.5, 'par2': 0.3},  # Mean 3 days, std dev 1 day, exposed to infectious
    'inf2sym': {'dist': 'lognormal', 'par1': 1.0, 'par2': 0.5},  # Infectious-to-symptomatic
    'sym2sev': {'dist': 'lognormal', 'par1': 3.0, 'par2': 0.5},  # Symptomatic to severe
    'sev2crit': {'dist': 'lognormal', 'par1': 4.0, 'par2': 0.5}, # Severe to critical
    'crit2die': {'dist': 'lognormal', 'par1': 4.0, 'par2': 0.5}, # Critical to death
    'asym2rec' : {'dist': 'lognormal', 'par1': 2, 'par2':0.5  }, # Asymptomatic to recovery
    'mild2rec' : {'dist': 'lognormal', 'par1': 4, 'par2':0.5  }, # Mild to recovery
    'sev2rec' : {'dist': 'lognormal', 'par1': 5, 'par2':0.5  },  # Severe to recovery
    'crit2rec' : {'dist': 'lognormal', 'par1': 5, 'par2':0.5  }, # Critical to recovery
}



# Function to compute the PDF of a lognormal distribution
def lognormal_pdf(x, mean, sigma):
    return lognorm.pdf(x, s=sigma, scale=mean)

# Define x values for the PDF
x_values = np.linspace(0, 15, 1000)  # Define a range of x values (0 to 30 days)


# Create a figure with subplots
n_distributions = len(durations)
fig, axes = plt.subplots(n_distributions, 1, figsize=(10, 5 * n_distributions), sharex=True)

# Plot each distribution in its own panel
x_max = 15  # Define a reasonable x-axis range

for ax, (key, params) in zip(axes, durations.items()):
    y_values = lognormal_pdf(x_values, params['par1'], params['par2'])
    ax.plot(x_values, y_values, label=f"{key} (mean={params['par1']}, std={params['par2']})")
    ax.set_title(f"Distribution: {key}")
    ax.set_ylabel("Density")
    ax.grid(True)
    ax.legend()

# Common x-axis label
plt.xlabel("Days")
plt.xlim(0, x_max)  # Ensure the same x-axis across all subplots
plt.tight_layout()
plt.show()





#Parameters for covasim of relevance
pars = dict(
    pop_size = 1e5,
    pop_infected = 100,
    start_day = '2025-01-01',
    end_day = '2025-07-01',
    beta = 0.03,  #per contact infection rate
    pop_scale=100,
    dur = durations,
    rel_symp_prob = 0.7,
    rel_severe_prob = 0.25, #among symptomatics how many become severe
    rel_crit_prob = 0.4, #among severe how many become critical
    pop_type='hybrid',
    rel_death_prob= 0.5,

)


sim = cv.Sim(pars)

sim.initialize() # Create people
fig = sim.people.plot() # Show statistics of the people

sim.run()

# Extract and plot infections by contact layer
results = sim.results
daily_infections = results['new_infections']  # Total daily infections
# Assume proportionate distribution for simplicity if per-layer data is not configured

plt.figure(figsize=(10, 6))
plt.plot(daily_infections, label='Total Infections')
plt.title("Daily Infections (All Settings)")
plt.xlabel("Days")
plt.ylabel("Number of Infections")
plt.legend()
plt.grid(True)
plt.show()


msim = cv.MultiSim(sim)
msim.run(n_runs=5)
msim.plot()



#Interventions:


