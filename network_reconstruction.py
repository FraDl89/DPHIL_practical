#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Workflow:
 1) Read reduced bipartite network (edges file produced earlier).
 2) Plot degree CCDF buyers vs escorts.
 3) Read survey file (node type degree) and build bipartite configuration
    model(s) from those reported degrees.
 4) Compare statistics (assortativity, clustering on projections) between
    real network and config-model ensemble.

Input files:
 - REDUCED_EDGES_PATH: space-delimited lines "B<id> E<id> n_contacts first_time last_time"
 - SURVEY_PATH: tab or space delimited "node<TAB>type<TAB>degree" where node is "B<id>" or "E<id>"

Adjust paths below.
"""

import random
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite
import statistics
import warnings

# -------------------- CONFIG --------------------
REDUCED_EDGES_PATH = "/Users/adminaccount/Downloads/ia-escorts-dynamic/ia-escorts-reduced.edges"
SURVEY_PATH        = "/Users/adminaccount/Downloads/ia-escorts-dynamic/ia-escorts-survey.tsv"
SEED = 2026
N_REPS = 10            # number of config-model replications to estimate variability
STRICT_BALANCE = False # if True, raise error when stub sums mismatch; otherwise trim randomly
# ------------------------------------------------

random.seed(SEED)
np.random.seed(SEED)




# ---------- helpers ----------
def read_reduced_network(path):
    """Return bipartite Graph G built from reduced edges file."""
    G = nx.Graph()
    with open(path, "r") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#") or s.startswith("%"):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            b, e = parts[0], parts[1]
            G.add_edge(b, e)
    return G

def read_survey_degrees(path):
    """
    Read survey file of lines: node type degree
    Return two dicts: deg_buyers (node->deg), deg_escorts (node->deg)
    """
    deg_b = {}
    deg_e = {}
    with open(path, "r") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#") or s.startswith("%"):
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            node = parts[0]
            nodetype = parts[1].lower()
            try:
                deg = int(float(parts[2]))
            except Exception:
                continue
            if nodetype.startswith("b") or node.startswith("B"):
                deg_b[node] = deg
            elif nodetype.startswith("e") or node.startswith("E"):
                deg_e[node] = deg
            else:
                # fallback by prefix
                if node.startswith("B"):
                    deg_b[node] = deg
                elif node.startswith("E"):
                    deg_e[node] = deg
                else:
                    # unknown: skip
                    continue
    return deg_b, deg_e

def build_bipartite_configuration(deg_b, deg_e, seed=None, strict=False):
    """
    Build bipartite configuration-model simple graph from degree dictionaries.
    deg_b, deg_e: dict node->degree
    If sums mismatch:
      - if strict True: raise ValueError
      - else: randomly trim stubs on larger side until sums match (warn)
    Returns simple Graph H with nodes union of deg_b and deg_e (no node attributes)
    """
    rnd = random.Random(seed)
    sum_b = sum(deg_b.values())
    sum_e = sum(deg_e.values())

    if sum_b != sum_e:
        msg = f"Stub sums mismatch: buyers={sum_b} escorts={sum_e}"
        if strict:
            raise ValueError(msg)
        warnings.warn(msg + "  → balancing by random trimming of extra stubs")
        # copy deg dicts to mutate
        deg_b = dict(deg_b)
        deg_e = dict(deg_e)
        # remove stubs randomly from the larger side until equal
        if sum_b > sum_e:
            excess = sum_b - sum_e
            # create list of buyer nodes repeated by degree, choose random excess to decrement
            buyer_stubs = []
            for n,k in deg_b.items():
                buyer_stubs.extend([n]*k)
            rnd.shuffle(buyer_stubs)
            to_trim = buyer_stubs[:excess]
            for n in to_trim:
                deg_b[n] -= 1
        else:
            excess = sum_e - sum_b
            escort_stubs = []
            for n,k in deg_e.items():
                escort_stubs.extend([n]*k)
            rnd.shuffle(escort_stubs)
            to_trim = escort_stubs[:excess]
            for n in to_trim:
                deg_e[n] -= 1

    # now sums equal
    sum_b = sum(deg_b.values())
    sum_e = sum(deg_e.values())
    if sum_b != sum_e:
        raise RuntimeError("Balancing failed; unequal stub counts remain")

    # build stub lists
    buyer_stubs = []
    for n,k in deg_b.items():
        if k < 0:
            raise ValueError("Negative degree after trimming")
        buyer_stubs.extend([n]*k)
    escort_stubs = []
    for n,k in deg_e.items():
        if k < 0:
            raise ValueError("Negative degree after trimming")
        escort_stubs.extend([n]*k)

    rnd.shuffle(buyer_stubs)
    rnd.shuffle(escort_stubs)

    # pair stubs
    M = nx.MultiGraph()
    for b_node, e_node in zip(buyer_stubs, escort_stubs):
        M.add_node(b_node)
        M.add_node(e_node)
        M.add_edge(b_node, e_node)

    # collapse to simple graph
    H = nx.Graph()
    H.add_nodes_from(M.nodes())
    for u,v in M.edges():
        if H.has_edge(u,v):
            continue
        H.add_edge(u,v)
    return H

def ccdf_from_degrees(degs):
    degs = np.asarray(degs)
    degs = degs[degs > 0]
    if degs.size == 0:
        return np.array([]), np.array([])
    ks = np.sort(np.unique(degs))
    ccdf_vals = np.array([(degs >= k).mean() for k in ks])
    return ks, ccdf_vals

def compute_stats_simple(G):
    """Compute assortativity (whole graph) and average clustering on projections."""
    stats = {}
    try:
        stats['assortativity'] = nx.degree_assortativity_coefficient(G)
    except Exception:
        stats['assortativity'] = float('nan')
    # identify bipartite sides by prefix
    buyers = [n for n in G.nodes() if str(n).startswith("B")]
    escorts = [n for n in G.nodes() if str(n).startswith("E")]
    # if sides detection fails, try bipartite.sets
    if not buyers or not escorts:
        if nx.is_bipartite(G):
            bset, eset = nx.algorithms.bipartite.sets(G)
            buyers, escorts = list(bset), list(eset)
    # buyer projection clustering
    try:
        proj_b = bipartite.projected_graph(G, buyers)
        proj_e = bipartite.projected_graph(G, escorts)
        stats['clustering_buyers'] = nx.average_clustering(proj_b) if proj_b.number_of_nodes()>0 else 0.0
        stats['clustering_escorts'] = nx.average_clustering(proj_e) if proj_e.number_of_nodes()>0 else 0.0
    except Exception:
        stats['clustering_buyers'] = float('nan')
        stats['clustering_escorts'] = float('nan')
    return stats

# ---------- plotting helpers ----------
def plot_degree_ccdf(G):
    buyers = [n for n in G.nodes() if str(n).startswith("B")]
    escorts = [n for n in G.nodes() if str(n).startswith("E")]
    deg_b = np.array([G.degree(n) for n in buyers])
    deg_e = np.array([G.degree(n) for n in escorts])
    kb, cb = ccdf_from_degrees(deg_b)
    ke, ce = ccdf_from_degrees(deg_e)
    plt.figure(figsize=(6,5))
    if kb.size: plt.plot(kb, cb, 'o-', label='buyers')
    if ke.size: plt.plot(ke, ce, 'o-', label='escorts')
    plt.xscale('log'); plt.yscale('log')
    plt.xlabel('degree k'); plt.ylabel('P(K ≥ k)')
    plt.title('Degree CCDF (real network)')
    plt.legend(); plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

def barplot_stats(real_stats, conf_mean, conf_sd):
    labels = ['Assortativity', 'Clust buyers', 'Clust escorts']
    real_vals = [real_stats['assortativity'], real_stats['clustering_buyers'], real_stats['clustering_escorts']]
    conf_vals = [conf_mean['assortativity'], conf_mean['clustering_buyers'], conf_mean['clustering_escorts']]
    conf_errs = [conf_sd['assortativity'], conf_sd['clustering_buyers'], conf_sd['clustering_escorts']]

    x = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(8,4))
    ax.bar(x - width/2, real_vals, width, label='real')
    ax.bar(x + width/2, conf_vals, width, yerr=conf_errs, capsize=5, label='config (mean ± sd)')
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("value")
    ax.set_title("Real network vs configuration model")
    ax.legend()
    plt.tight_layout()
    plt.show()

# ---------- main ----------
def main():
    # 1) read real reduced network and plot CCDF
    G_real = read_reduced_network(REDUCED_EDGES_PATH)
    print("Real reduced network: nodes", G_real.number_of_nodes(), "edges", G_real.number_of_edges())
    plot_degree_ccdf(G_real)
    real_stats = compute_stats_simple(G_real)
    print("Real network stats:", real_stats)

    # 2) read survey degrees (this is the file you generated earlier)
    deg_b, deg_e = read_survey_degrees(SURVEY_PATH)
    print("Survey degrees read: buyers:", len(deg_b), " escorts:", len(deg_e))
    sum_b = sum(deg_b.values()); sum_e = sum(deg_e.values())
    print("Sum stubs: buyers", sum_b, " escorts", sum_e)

    # quick sanity: if some nodes in survey don't appear in real network, warn
    missing_buyers = [n for n in deg_b if n not in G_real.nodes()]
    missing_escorts = [n for n in deg_e if n not in G_real.nodes()]
    if missing_buyers or missing_escorts:
        warnings.warn(f"{len(missing_buyers)} buyers and {len(missing_escorts)} escorts from survey not in reduced network nodes.")

    # 3) build N_REPS config-models and compute stats
    conf_stats_list = []
    for rep in range(N_REPS):
        confG = build_bipartite_configuration(deg_b, deg_e, seed=SEED + rep, strict=STRICT_BALANCE)
        s = compute_stats_simple(confG)
        conf_stats_list.append(s)
    # compute mean and sd for each stat
    keys = ['assortativity', 'clustering_buyers', 'clustering_escorts']
    conf_mean = {k: statistics.mean([s[k] for s in conf_stats_list]) for k in keys}
    conf_sd   = {k: statistics.pstdev([s[k] for s in conf_stats_list]) for k in keys}  # population sd

    print("\nConfiguration model (ensemble) mean stats:", conf_mean)
    print("Configuration model (ensemble) sd:", conf_sd)

    # 4) plot comparison barplot
    barplot_stats(real_stats, conf_mean, conf_sd)

if __name__ == "__main__":
    main()
