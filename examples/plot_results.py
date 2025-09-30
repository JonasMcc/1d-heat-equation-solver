import numpy as np
import matplotlib.pyplot as plt

# Load data
x_fe, u_fe = np.loadtxt("example-outputs/HeatFE.dat", unpack=True)
x_be, u_be = np.loadtxt("example-outputs/HeatBE.dat", unpack=True)
x_cn, u_cn = np.loadtxt("example-outputs/HeatCN.dat", unpack=True)

# Create figure with 3 subplots
fig, axes = plt.subplots(1, 3, figsize=(16, 6), sharex=True)

# Forward Euler
axes[0].plot(x_fe, u_fe, 'bo')
axes[0].set_xlabel("x")
axes[0].set_ylabel("u(x, t)")
axes[0].set_title("Forward Euler Method")

# Backward Euler
axes[1].plot(x_be, u_be, color='green', marker='o', linestyle='-', markersize=4)
axes[1].set_xlabel("x")
axes[1].set_ylabel("u(x, t)")
axes[1].set_title("Backward Euler Method")

# Crank--Nicolson
axes[2].plot(x_cn, u_cn, color='red', marker='o', linestyle='-', markersize=4)
axes[2].set_xlabel("x")
axes[2].set_ylabel("u(x, t)")
axes[2].set_title("Crank--Nicolson Method")

# Show grid on all axes. 
for ax in axes: 
    ax.grid(True)

plt.tight_layout()
plt.show()

def plotExactSolution(): 
    # Parameters
    alpha = 1.0  # thermal diffusivity
    x = np.linspace(0, 1, 500)  # spatial domain
    times = [0, 0.1, 0.2, 0.3, 0.8, 1.0]  # times to plot

    # Plot
    plt.figure(figsize=(8,5))
    for t in times:
        u = np.sin(np.pi * x) * np.exp(-alpha * (np.pi**2) * t)
        plt.plot(x, u, label=f't={t}')

    plt.title("Heat Equation Solution: u(x,t) = sin(pi*x) exp(-alpha*pi^2*t)")
    plt.xlabel("x")
    plt.ylabel("u(x,t)")
    plt.legend()
    plt.grid(True)
    plt.show()

plotExactSolution()