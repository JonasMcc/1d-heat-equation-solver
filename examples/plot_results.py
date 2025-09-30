import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# Find result files relative to project root

data = np.loadtxt("example-outputs/HeatCN.dat")

x = data[:, 0]
u = data[:, 1]

plt.plot(x, u, label="Solution")

plt.xlabel("x")
plt.ylabel("u(x, t)")
plt.legend()
plt.title("Crank--Nicolson for 1D Heat Equation")
plt.show()
# Plot every 10th snapshot
#for filename in files[::10]:
#    data = np.loadtxt(filename)
#    x, u = data[:, 0], data[:, 1]#

    #step = int(filename.split("step")[-1].split(".")[0])
    #plt.plot(x, u, label=f"step {step}")

#plt.xlabel("x")
#plt.ylabel("u(x, t)")
#plt.legend()
#plt.title("1D Heat Equation")
#plt.show()