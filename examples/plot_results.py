from typing import List
import numpy as np
import matplotlib.pyplot as plt
import glob
import re

from dataclasses import dataclass, field

@dataclass
class Parameters:
    dt: float = 0.0001
    alpha: float = 1.0
    Nx: int = 50
    T: float = 0.1

    x: np.ndarray = field(init=False)
    Nt: int = field(init=False)
    save_steps: List = field(init=False)

    def __post_init__(self):
        self.x = np.linspace(0, 1, self.Nx + 1)
        self.Nt = int(self.T / self.dt)
        self.save_steps = [0, self.Nt // 2, self.Nt]

p = Parameters()

def plot_simulation(
    solver: str,
    parameters: Parameters,
    stride: int = 3
) -> None:
    suffix: str = f"Heat{solver.upper()}_step*"
    fpattern: str = f"example-outputs/full-simulation-{solver.upper()}/{suffix}"
    files = glob.glob(fpattern)
    fsorted: List = sorted(files, key=lambda f: int(re.search(r'step(\d+)\.dat', f).group(1)))
    fselected: List = [f for f in fsorted if int(re.search(r'step(\d+)\.dat', f).group(1)) in p.save_steps]

    plt.figure(figsize=(15,5))

    for index, file in enumerate(fselected):
        step_index: int = int(re.search(r'step(\d+)\.dat', file).group(1))
        t: float = step_index * p.dt

        # Numerical solution
        u_numeric = np.loadtxt(file)[:, 1]

        # Exact solution
        u_exact = np.exp(-p.alpha * np.pi**2 * t) * np.sin(np.pi * p.x)

        plt.subplot(1, 3, index+1)
        plt.plot(p.x, u_exact, 'r-', label='Exact')
        plt.plot(p.x[::stride], u_numeric[::stride], 'bo', markersize=4, label='Numerical')
        plt.xlabel('x')
        plt.ylabel('u(x, t)')
        plt.title(f"t = {t:.3f}")
        plt.ylim(0, 1.1)
        plt.legend()
    
    solver_names = {
        "fe": "Forward Euler",
        "be": "Backward Euler",
        "cn": "Crank--Nicolson"
    }

    plt.suptitle(f"1D Heat Equation: {solver_names.get(solver)} Numerical vs Exact Solution", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()   

plot_simulation(solver="cn", parameters=p)


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