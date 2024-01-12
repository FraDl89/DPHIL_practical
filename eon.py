#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 13:45:17 2024

@author: adminaccount
"""


import EoN
import networkx as nx
from matplotlib import pyplot as plt

N=10000
p=0.001
G=nx.erdos_renyi_graph(N,p) #create a barabasi-albert graph

actual_degrees = [d for v, d in G.degree()]

plt.hist(actual_degrees, bins=100)
plt.show()

tmax = 20
iterations = 10  #run 5 simulations
tau = 0.2           #transmission rate
gamma = 1.0    #recovery rate
rho = 0.005      #random fraction initially infected

for counter in range(iterations): #run simulations
    t, S, I, R = EoN.fast_SIR(G, tau, gamma, rho=rho, tmax = tmax)
    if counter == 0:
        plt.plot(t, I, color = 'k', alpha=0.3, label='Simulation')
    plt.plot(t, I, color = 'k', alpha=0.3)
    

#approximation
t, S, I, R = EoN.SIR_homogeneous_meanfield(S0=9950, I0=50, R0=0, n=(N-1)*p,tmax=25, tau=tau, gamma=gamma)
plt.plot(t, I, '-.', label = 'Homogeneous pairwise', linewidth = 5)


