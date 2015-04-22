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
ax1.plot(f[:,0],f[:,1:3])
ax1.set_xlabel("M*n - log scale")
ax1.set_ylabel("Half confidence interval")
ax1.set_xscale('log')

fig2 = plt.figure()
ax2 = fig2.gca()
ax2.plot(f[:,0],f[:,3:5])
ax2.set_xlabel("M*n - log scale")
ax2.set_ylabel("Half confidence interval")
ax2.set_xscale('log')

# Save figures
save_params = {"bbox_inches":'tight', "pad_inches":0.}

fig1.savefig(target + "__CI.png",**save_params)
fig2.savefig(target + "__time.png",**save_params)
