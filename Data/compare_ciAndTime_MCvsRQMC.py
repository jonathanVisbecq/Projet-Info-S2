import sys
import matplotlib.pyplot as plt
import numpy as np

directory = "/home/jonathan/Programmation/C++/Projet_C++_S2"

# Get target path
target = directory + "/Images/compare_ciAndTime_MCvsRQMC"

# Read data
f = np.genfromtxt(fname=directory + "/Data/compare_ciAndTime_MCvsRQMC.dat", delimiter=" ")

# Make figures
fig1 = plt.figure()
ax1 = fig1.gca()
ax1.plot(f[:,0],f[:,1:4])
ax1.set_xlabel("M*n - log scale")
ax1.set_xscale('log')
ax1.set_ylabel("Half confidence interval - log scale")
ax1.set_yscale('log')
ax1.legend(["Monte Carlo","Shifted QMC","Random-start Halton"],loc='best')

fig2 = plt.figure()
ax2 = fig2.gca()
ax2.plot(f[:,0],f[:,4:7])
ax2.set_xlabel("M*n - log scale")
ax2.set_xscale('log')
ax2.set_ylabel("Half confidence interval - log scale")
ax2.set_yscale('log')
ax2.legend(["Monte Carlo","Shifted QMC","Random-start Halton"],loc='best')

# Save figures
save_params = {"bbox_inches":'tight', "pad_inches":0.}

fig1.savefig(target + "__CI.png",**save_params)
fig2.savefig(target + "__time.png",**save_params)
