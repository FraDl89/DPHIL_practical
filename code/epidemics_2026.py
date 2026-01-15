"""
SIR on ER and BA networks with explicit dt, interpolation to fixed grid,
and ODE overlay using beta_ode = beta_edge * mean_degree.

Produces:
 - individual runs (thin gray lines),
 - averaged prevalence I(t)/N (colored),
 - ODE mean-field curves (dashed).
"""
import random
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# -----------------------------
# Discrete-time SIR with dt
# -----------------------------
def run_sir_discrete(G, beta, mu, dt,
                     initial_infected=1,
                     max_time=100.0,
                     rng=None):
    """
    Discrete-time SIR with explicit dt.
    beta, mu are rates (1/time).
    Returns: t_hist, S_hist, I_hist, R_hist, node_history
    """
    if rng is None:
        rng = random.Random()

    nodes = list(G.nodes())
    n = len(nodes)
    node_index = {u: i for i, u in enumerate(nodes)}

    state = np.zeros(n, dtype=int)  # 0=S,1=I,2=R

    # seeds
    if isinstance(initial_infected, int):
        k = max(1, min(initial_infected, n))
        seeds = rng.sample(nodes, k)
    else:
        seeds = list(initial_infected)

    for u in seeds:
        state[node_index[u]] = 1

    node_history = np.full((n, 2), np.nan)  # t_infected, t_recovered
    for u in seeds:
        node_history[node_index[u], 0] = 0.0

    t = 0.0
    t_hist = [t]
    S_hist = [int(np.sum(state == 0))]
    I_hist = [int(np.sum(state == 1))]
    R_hist = [int(np.sum(state == 2))]

    while t < max_time and np.any(state == 1):
        to_infect = []
        to_recover = []

        infected_idxs = set(np.where(state == 1)[0])

        for i in range(n):
            if state[i] == 0:
                # infected neighbours count
                k_inf = 0
                for nb in G.neighbors(nodes[i]):
                    if node_index[nb] in infected_idxs:
                        k_inf += 1
                if k_inf > 0:
                    p_inf = 1.0 - np.exp(-beta * k_inf * dt)
                    if rng.random() < p_inf:
                        to_infect.append(i)
            elif state[i] == 1:
                p_rec = 1.0 - np.exp(-mu * dt)
                if rng.random() < p_rec:
                    to_recover.append(i)

        # apply updates simultaneously
        for i in to_infect:
            if state[i] == 0:
                state[i] = 1
                if np.isnan(node_history[i, 0]):
                    node_history[i, 0] = t + dt
        for i in to_recover:
            if state[i] == 1:
                state[i] = 2
                if np.isnan(node_history[i, 1]):
                    node_history[i, 1] = t + dt

        t += dt
        t_hist.append(t)
        S_hist.append(int(np.sum(state == 0)))
        I_hist.append(int(np.sum(state == 1)))
        R_hist.append(int(np.sum(state == 2)))

    return np.array(t_hist), np.array(S_hist), np.array(I_hist), np.array(R_hist), node_history

# -----------------------------
# Deterministic ODE SIR
# -----------------------------
def sir_ode(y, t, beta, mu):
    S, I, R = y
    return [-beta * S * I, beta * S * I - mu * I, mu * I]

def solve_sir_ode(beta, mu, S0, I0, R0, t_grid):
    return odeint(sir_ode, [S0, I0, R0], t_grid, args=(beta, mu))

# -----------------------------
# Run many sims, interpolate to common grid
# -----------------------------
def sims_interpolate_I(G, beta, mu, dt, t_grid,
                       n_runs=20, initial_infected=1, seed=0):
    rng_master = random.Random(seed)
    I_matrix = np.zeros((n_runs, len(t_grid)))

    for r in range(n_runs):
        rng = random.Random(rng_master.randint(0, 2**31 - 1))
        t_run, S_run, I_run, R_run, _ = run_sir_discrete(G, beta, mu, dt,
                                                        initial_infected=initial_infected,
                                                        max_time=t_grid[-1],
                                                        rng=rng)
        # interpolate I_run (counts) onto t_grid
        # t_run is monotonic; outside range we use last value (via np.interp)
        I_interp = np.interp(t_grid, t_run, I_run, left=I_run[0], right=I_run[-1])
        I_matrix[r, :] = I_interp

    I_mean = I_matrix.mean(axis=0)
    return I_matrix, I_mean

# -----------------------------
# Main experiment
# -----------------------------
def main():
    # params
    N = 1000
    beta = 0.1    # per-edge transmission rate
    mu = 0.4      # recovery rate
    dt = 0.01
    max_time = 30.0
    n_runs = 30
    initial_infected = 50
    rng_seed = 2026

    # networks: ER and BA
    p_er = 1/100
    G_er = nx.erdos_renyi_graph(N, p_er, seed=1)
    G_ba = nx.barabasi_albert_graph(N, 5, seed=1)

    networks = [("Erdos-Renyi", G_er), ("Barabasi-Albert", G_ba)]

    # common time grid
    t_grid = np.arange(0.0, max_time + dt/2, dt)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

    for ax, (name, G) in zip(axes, networks):
        print(f"Running {n_runs} sims on {name}...")
        I_matrix, I_mean = sims_interpolate_I(G, beta, mu, dt, t_grid,
                                              n_runs=n_runs,
                                              initial_infected=initial_infected,
                                              seed=rng_seed)

        # plot individual runs (subset)
        for i in range(min(12, n_runs)):
            ax.plot(t_grid, I_matrix[i, :] / N, color="gray", alpha=0.3, linewidth=0.8)

        # plot mean
        ax.plot(t_grid, I_mean / N, color="C1", lw=2.2, label="Network mean")

        # overlay ODE 
        beta_ode = 1
        S0 = (N - initial_infected) / N
        I0 = initial_infected / N
        R0 = 0.0
        sol = solve_sir_ode(beta_ode, mu, S0, I0, R0, t_grid)
        ax.plot(t_grid, sol[:, 1], linestyle="--", color="k", lw=1.8,
                label=f"ODE β={beta_ode:.2f}")

        ax.set_title(name)
        ax.set_xlabel("time")
        ax.grid(alpha=0.25)
        if ax is axes[0]:
            ax.set_ylabel("I(t) / N")
        ax.legend()

    plt.suptitle(f"SIR on ER and BA (beta={beta}, mu={mu}, dt={dt}), {n_runs} runs each")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()
    
    #Task: what is the correct mapping from a network ode per-link transmission
    #parameter and the mass-action ode \beta? Why does it fail for B-A?

if __name__ == "__main__":
    main()
