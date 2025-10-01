import numpy as np
import matplotlib.pyplot as plt
import glob
import re

# Parameter
dt = 0.0001
alpha = 1.0
Nx = 50
T = 0.1

# Spatial grid
x = np.linspace(0, 1, Nx+1)
Nt = int(T / dt)

file_pattern = "example-outputs/full-simulation-FE/HeatFE_step*.dat"
files = glob.glob(file_pattern)
files_sorted = sorted(files, key=lambda f: int(re.search(r'step(\d+)\.dat', f).group(1)))

save_steps = [0, Nt // 2, Nt]

selected_files = [f for f in files_sorted if int(re.search(r'step(\d+)\.dat', f).group(1)) in save_steps]

plt.figure(figsize=(15,5))

for index, file in enumerate(selected_files):
    step_index = int(re.search(r'step(\d+)\.dat', file).group(1))
    t = step_index * dt

    # Numerical solution
    u_numeric = np.loadtxt(file)[:, 1]

    # Exact solution
    u_exact = np.exp(-alpha * np.pi**2 * t) * np.sin(np.pi * x)

    plt.subplot(1, 3, index+1)
    plt.plot(x, u_exact, 'r-', label='Exact')

    stride = 3
    plt.plot(x[::stride], u_numeric[::stride], 'bo', markersize=4, label='Numerical')
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.title(f"t = {t:.3f}")
    plt.ylim(0, 1.1)
    plt.legend()

plt.suptitle("1D Heat Equation: Forward Euler Numerical vs Exact Solution", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

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