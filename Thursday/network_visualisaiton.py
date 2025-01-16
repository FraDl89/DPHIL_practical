#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:46:07 2025

@author: adminaccount
"""

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy
import powerlaw
%matplotlib qt

##############################################################################
# Introduction
##############################################################################
#To generate a specific network in networkx, you start off by defining an empty graph
        
G = nx.empty_graph(n=4)

#Iterate through nodes, that have an id assigned by default
print(G.nodes())

#removing or adding nodes manually is done with these commands
#Note: this changes the structure of the network, it doesn't generate a new network
G.add_nodes_from([4,5,6])
G.remove_nodes_from([1,2,3])
#Notice how the nodes id now are different.
# You can be deliberate when choosing nodes ids
G.add_nodes_from(["Mike", "Karl"])
G.nodes()
#To add edges between nodes the command is similar:
G.add_edges_from([[0,4], ["Mike", "Karl"]])
#To view the Graph as an adjacency list (sparse array)
adjlist=nx.adjacency_matrix(G)
adjmatrix=adjlist.toarray()
#plot the graph
plt.figure()
nx.draw_networkx(G, with_labels=True)
 
#Networkx also incorporates some very well-known network models

#Regular network
G = nx.random_regular_graph(4, 30)
nx.draw_networkx(G, with_labels=True)

#Notice how every time you draw the network, you get a different layout
#This can be fixed once and for all by defining, for each node, its position


#define a position for each node. Either manually or with some pre-defined layout
pos = nx.spring_layout(G)

plt.figure()
nx.draw(G, pos=pos, node_size=25, width=0.5, alpha=0.5, edge_color='gray')  


#you can get the degree distribution from the adjacency matrix
#or you can use the built-in function, but it is not very handy

nx.degree_histogram(G)


G = nx.erdos_renyi_graph(1000, 0.24)
#Instead, you can get a view of the degree of each node:
    
deg=dict(nx.degree(G))

plt.hist(deg.values())

#scale-free networks are more realistic for many processes. One way to get a scale-free
#network is the barabasi albert graph

##############################################################################
# Scale-free network
##############################################################################
BA = nx.barabasi_albert_graph(1000, 3)
#To get the degree distribution, we iterate through the network degrees
degree_sequence = [d for n, d in BA.degree()]
degree_count = np.bincount(degree_sequence) #bin them
degrees = np.arange(len(degree_count))
# Filter out degrees with zero count
degrees = degrees[degree_count > 0]
degree_count = degree_count[degree_count > 0]

# Plot the histogram on a log-log scale
plt.figure(figsize=(8, 6))
plt.loglog(degrees, degree_count, 'o', markersize=8, label='Degree Distribution')
plt.title("Degree Distribution (Log-Log Scale)", fontsize=14)
plt.xlabel("Degree (k)", fontsize=12)
plt.ylabel("Count (P(k))", fontsize=12)
plt.grid(True, which="both", linestyle='--', linewidth=0.5)
plt.legend()
plt.show()

#This is a power-law, and indeed we can fit it with that distribution
fit = powerlaw.Fit(degree_sequence, discrete=True)
alpha = fit.power_law.alpha
xmin = fit.power_law.xmin


# Plot the empirical degree distribution
plt.figure(figsize=(8, 6))
hist, bin_edges = np.histogram(degree_sequence, bins=np.logspace(np.log10(1), np.log10(max(degree_sequence)), num=20), density=True)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
plt.loglog(bin_centers, hist, 'o', markersize=8, label='Empirical Data')

x = np.arange(xmin, max(degree_sequence) + 1)
C = hist[bin_centers >= xmin][0]  # Match the normalization at xmin
y = C * (x / xmin) ** (-alpha)  # Normalized power-law starting at xmin
plt.loglog(x, y, '--', label=f'Power-Law Fit ($\\alpha={alpha:.2f}$)', linewidth=2)
# Add labels, grid, and legend
plt.title("Degree Distribution with Power-Law Fit", fontsize=14)
plt.xlabel("Degree (k)", fontsize=12)
plt.ylabel("Frequency (P(k))", fontsize=12)
plt.grid(True, which="both", linestyle='--', linewidth=0.5)
plt.legend(fontsize=12)
plt.show()

#visualisation


# Create a spring layout for visualization
pos = nx.kamada_kawai_layout(BA, weight=None, scale=3,tol=1e-3)

# Normalize node size based on degree
node_sizes =15*np.log1p(np.array(degree_sequence))  # Scale node size for visualization


# Define degree bins and corresponding colors
bins = [0, 5, 50, float('inf')]  # Degree bins
colors = ['black', 'orange', 'red']  # Colors for each bin 

# Assign colors based on degree

bin_indexes = np.searchsorted(bins, degree_sequence)-1
node_colors = [colors[i] for i in bin_indexes]

# Plot the network
plt.figure(figsize=(18, 12))
nx.draw_networkx(BA,pos,with_labels=False, node_size=node_sizes,
    edge_color='gray',
    alpha=0.7,
    width=0.1,
    node_color=node_colors,  # Color nodes based on degree
    cmap=plt.cm.plasma  # Use a colormap to emphasize hubs
)


# Add a colorbar to indicate degree
sm = plt.cm.ScalarMappable(cmap=plt.cm.plasma, norm=plt.Normalize(vmin=min(), vmax=max(degrees.values())))
sm.set_array([])
plt.colorbar(sm, label='Node Degree')

plt.title("Network Visualization with Hubs Highlighted", fontsize=16)
plt.axis('off')
plt.show()

#Check degree degree correlation

degrees = dict(BA.degree())

# Extract degree pairs for each edge
degree_pairs = [(degrees[u], degrees[v]) for u, v in BA.edges()]

# Split the pairs into two lists

# Plot the degree-degree correlation
plt.figure(figsize=(10, 8))
plt.scatter(degree_u, degree_v, alpha=0.5, edgecolor='k', s=10)
plt.xscale('log')  # Use log scale for better visualization
plt.yscale('log')
plt.xlabel("Degree of Node u ($k_u$)", fontsize=12)
plt.ylabel("Degree of Node v ($k_v$)", fontsize=12)
plt.title("Degree-Degree Correlation", fontsize=16)
plt.grid(True, which="both", linestyle='--', linewidth=0.5)
plt.show()

assortativity = nx.degree_assortativity_coefficient(BA)

print(assortativity)

clustering = nx.average_clustering(BA)

print(clustering)


#############################################################################
# A real network
############################################################################
#download from https://github.com/RishujeetRai/DolphinsSocialNetworkAnalysis/blob/main/dolphins.gml

#contains an undirected social network of frequent associations between 62 dolphins 
#in a community living off Doubtful Sound, New Zealand, as compiled by Lusseau et al. (2003)
dolphins = "/Users/adminaccount/Downloads/dolphins.gml"

dolphin_net = nx.read_gml(dolphins)


#Exercise: analyse this network
#    1) Plot it 
#    2) look at the degree distribution
#    3) calculate assortativity
#plt.figure(figsize=(12,9))
#pos = nx.circular_layout(dolphin_net,)
#nx.draw_networkx(dolphin_net, pos=pos)













