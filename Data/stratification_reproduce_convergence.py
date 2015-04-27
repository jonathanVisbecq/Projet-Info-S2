import sys
import os
import matplotlib.pyplot as plt
import numpy as np

directory = "/home/jonathan/Programmation/C++/Projet_C++_S2"

# Get target path
target = directory + "/Images/stratification_reproduce_convergence"

# Image saving parameters
save_params = {"bbox_inches":'tight', "pad_inches":0.}

# Path for the data
path = directory + "/Data/stratification_reproduce_convergence.dat" 

# Retrieve data
f = np.genfromtxt(fname=path, delimiter=" ")

# True optimal proportion
prop = 0.04685

# Make figures
# Second moment
fig = plt.figure()
ax = fig.gca()
ax.plot(f[:,0],f[:,1:3])
ax.axhline(y=prop,color='red')
ax.set_xlabel("Number of drawings")
ax.set_ylabel("Proportion of drawings in fifth stratum")
ax.legend(["Method a)","Method b)"],loc='best')

# Save figures
fig.savefig(target + ".png",**save_params)

