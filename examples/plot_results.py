import numpy as np
import matplotlib.pyplot as plt



# I want to split it up, plot forward Euler etc. 
def plot_simulation_final_time():
    # Load data
    x_fe, u_fe = np.loadtxt("example-outputs/simulation-final-time/HeatFE.dat", unpack=True)    
    x_be, u_be = np.loadtxt("example-outputs/simulation-final-time/HeatBE.dat", unpack=True)
    x_cn, u_cn = np.loadtxt("example-outputs/simulation-final-time/HeatCN.dat", unpack=True)

    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(16, 6), sharex=True)

    # Forward Euler
    axes[0].plot(x_fe, u_fe, 'bo')
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("u(x, t)")
    axes[0].set_title("Forward Euler Method")

    # Backward Euler
    axes[1].plot(x_be, u_be, 'go')
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("u(x, t)")
    axes[1].set_title("Backward Euler Method")

    # Crank--Nicolson
    axes[2].plot(x_cn, u_cn, 'ro')
    axes[2].set_xlabel("x")
    axes[2].set_ylabel("u(x, t)")
    axes[2].set_title("Crank--Nicolson Method")

    # Show grid on all axes. 
    for ax in axes: 
        ax.grid(True)

    plt.tight_layout()
    plt.show()