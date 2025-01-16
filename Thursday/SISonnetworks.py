import numpy as np
import networkx as nx
import random
import scipy.sparse as sp
import matplotlib.pyplot as plt

def sis_fast_simulation(G, tau, gamma, initial_infected, time_steps, seed=None):
    """
    Simulates a sparse matrix-based SIS model, recording the number of infected nodes at specified time steps.

    Parameters:
        G (networkx.Graph): The network.
        tau (float): Infection rate.
        gamma (float): Recovery rate.
        initial_infected (int): Initial number of infected nodes.
        time_steps (list): Time points at which to record the number of infected nodes.
        seed (int): Random seed for reproducibility (optional).

    Returns:
        times (list): List of recorded time points.
        infected_counts (list): List of the number of infected nodes at each recorded time point.
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    # Sparse adjacency matrix
    A = nx.to_scipy_sparse_array(G, format="csr")
    N = G.number_of_nodes()

    # Initialize states: 0 = susceptible, 1 = infected
    states = np.zeros(N, dtype=int)
    initial_infected_nodes = np.random.choice(N, initial_infected, replace=False)
    states[initial_infected_nodes] = 1

    # Precompute rates
    infection_rates = tau * A
    recovery_rates = gamma

    # Track infected counts
    infected_counts = []
    current_time = 0.0
    time_index = 1
    current_infected = np.sum(states)

    infected_counts.append(current_infected)

    while time_index < len(time_steps):
        # Calculate infection and recovery rates for each node
        
        #Exercise: instead of recalculating this product (very slow), why don't
        #we just update the rates around the nodes who changed their state
        #in the previous iteration? This will make the code go much faster on
        #larger networks 
        infected_neighbors = infection_rates.dot(states)
        infected_neighbors[states == 1] = 0  # Set infection rate to 0 for already infected nodes
        total_rates = infected_neighbors + recovery_rates * states

        # Time to next event
        total_rate = np.sum(total_rates)
        if total_rate == 0:
            break

        delta_time = np.random.exponential(1 / total_rate)
        current_time += delta_time

        # Process events until the next time step
        while time_index < len(time_steps) and current_time >= time_steps[time_index]:
            infected_counts.append(current_infected)
            time_index += 1

        # Choose which node and event occurs
        event_probs = total_rates / total_rate
        selected_node = np.random.choice(N, p=event_probs)

        if states[selected_node] == 1:  # Recovery event
            states[selected_node] = 0
            current_infected -= 1
        else:  # Infection event
            states[selected_node] = 1
            current_infected += 1

    # Record remaining time steps if simulation ends early
    while time_index < len(time_steps):
        infected_counts.append(current_infected)
        time_index += 1

    return time_steps, infected_counts




def mean_field_sis(G, tau, gamma, initial_infected, time_steps):
    """
    Computes the mean-field approximation for SIS epidemics.

    Parameters:
        G (networkx.Graph): The network.
        tau (float): Infection rate.
        gamma (float): Recovery rate.
        initial_infected (int): Initial number of infected nodes.
        time_steps (list): Time points at which to compute the approximation.

    Returns:
        times (list): List of recorded time points.
        infected_counts (list): List of the mean-field number of infected nodes at each recorded time point.
    """
    N = G.number_of_nodes()
    avg_degree = np.mean([d for _, d in G.degree()])

    # Initialize infected fraction
    infected_fraction = initial_infected / N

    infected_counts = []
    for t in time_steps:
        infected_counts.append(infected_fraction * N)
        # Mean-field SIS equation
        dI_dt = tau * avg_degree * infected_fraction * (1 - infected_fraction) - gamma * infected_fraction
        infected_fraction += dI_dt * (time_steps[1] - time_steps[0])  # Euler integration
        infected_fraction = max(0, min(1, infected_fraction))  # Keep within [0, 1]

    return time_steps, infected_counts

def run_simulation_ensemble(G, tau, gamma, initial_infected, time_steps, num_simulations, seed=None):
    """
    Runs multiple simulations of the SIS model and computes the average infected count over time.

    Parameters:
        G (networkx.Graph): The network.
        tau (float): Infection rate.
        gamma (float): Recovery rate.
        initial_infected (int): Initial number of infected nodes.
        time_steps (list): Time points at which to record the number of infected nodes.
        num_simulations (int): Number of stochastic simulations to run.
        seed (int): Random seed for reproducibility (optional).

    Returns:
        times (list): List of recorded time points.
        avg_infected_counts (list): Average number of infected nodes over all simulations.
    """
    all_infected_counts = []

    for i in range(num_simulations):
        simulation_seed = None if seed is None else seed + i
        _, infected_counts = sis_fast_simulation(G, tau, gamma, initial_infected, time_steps, seed=simulation_seed)
        all_infected_counts.append(infected_counts)

    avg_infected_counts = np.mean(all_infected_counts, axis=0)
    return time_steps, avg_infected_counts

# Example Usage
if __name__ == "__main__":
    # Generate an Erdos-Renyi network
    N = 1000  # Number of nodes
    p = 0.01  # Probability of edge creation
    G = nx.erdos_renyi_graph(N, p, seed=42)
    # Parameters
    tau = 0.1  # Infection rate
    gamma = 0.05  # Recovery rate
    initial_infected = 20  # Initial number of infected nodes
    time_steps = np.linspace(0, 20, 21)  # Record at these time steps

    # Run the simulation
    times, infected_counts = sis_fast_simulation(G, tau, gamma, initial_infected, time_steps, seed=42)

    # Plot the results
    plt.figure(figsize=(10, 6))
    plt.plot(times, infected_counts, marker="o", label="Infected Nodes")
    plt.xlabel("Time", fontsize=12)
    plt.ylabel("Number of Infected Nodes", fontsize=12)
    plt.legend()
    plt.grid()
    plt.show()

    #now run an ensemble of simulations

    num_simulations = 100 # Number of simulations

    # Run the stochastic ensemble
    times_ensemble, avg_infected_counts = run_simulation_ensemble(G, tau, gamma, initial_infected, time_steps, num_simulations, seed=42)

    # Run the mean-field approximation. Note, time-step should be really small
    time_steps_meanfield = np.linspace(0, 20, 401 )
    times_mean_field, infected_counts_mean_field = mean_field_sis(G, tau, gamma, initial_infected, time_steps_meanfield)

    # Plot the results
    plt.figure(figsize=(10, 6))
    plt.plot(times_ensemble, avg_infected_counts, marker="o", label="Stochastic Ensemble (Avg)")
    plt.plot(times_mean_field, infected_counts_mean_field, linestyle="--", label="Mean-Field Approximation")
    plt.xlabel("Time", fontsize=12)
    plt.ylabel("Number of Infected Nodes", fontsize=12)
    plt.title("SIS Model: Stochastic Ensemble vs. Mean-Field", fontsize=14)
    plt.legend()
    plt.grid()
    plt.show()
    
    
    
