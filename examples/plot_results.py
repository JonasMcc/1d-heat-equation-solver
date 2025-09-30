import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("example-outputs/HeatCN.dat")

x = data[:, 0]
u = data[:, 1]

plt.plot(x, u, label="Solution")

plt.xlabel("x")
plt.ylabel("u(x, t)")
plt.legend()
plt.title("Crank--Nicolson for 1D Heat Equation")
plt.show()