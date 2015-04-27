import sys
import os
import matplotlib.pyplot as plt
import numpy as np

directory = "/home/jonathan/Programmation/C++/Projet_C++_S2"

# Get target path
target = directory + "/Images/compare_gaussian"

# Image saving parameters
save_params = {"bbox_inches":'tight', "pad_inches":0.}

#%% Plot stratification
# Retrieve data
strat = np.genfromtxt(fname=directory + "/Data/compare_gaussian_strat.dat", delimiter=" ")
rqmc = np.genfromtxt(fname=directory + "/Data/compare_gaussian_rqmc.dat", delimiter=" ")

idx = rqmc[:,1]==800

# Make figures
fig = plt.figure()
ax = fig.gca()
ax.plot(strat[:,2], strat[:,1], c='red')
ax.plot(rqmc[idx,4], rqmc[idx,2], c='blue')
ax.plot(rqmc[idx,3], rqmc[idx,5], c='green')
ax.set_xlabel("CI - log scale")
ax.set_xscale('log')
ax.set_ylabel("Times (s) - log scale")
ax.set_yscale('log')

# Save figures
fig.savefig(target + ".png",**save_params)
















