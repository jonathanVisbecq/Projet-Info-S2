import sys
import os
import matplotlib.pyplot as plt
import numpy as np

directory = "/home/jonathan/Programmation/C++/Projet_C++_S2"

# Get target path
target = directory + "/Images/stratification_reproduce_variance"

# Image saving parameters
save_params = {"bbox_inches":'tight', "pad_inches":0.}

# Path for the data
path = directory + "/Data/stratification_reproduce_variance.dat" 

# Retrieve data
f = np.genfromtxt(fname=path, delimiter=" ")

# True sigma_star
sigma_star = 0.1559335

# Make figures
# Second moment
fig_a = plt.figure()
ax_a = fig_a.gca()
ax_a.plot(f[:,0],f[:,1],c='orange')
ax_a.plot(f[:,0],f[:,3],c='black')
ax_a.axhline(y=sigma_star,color='red')
ax_a.set_xlabel("Number of drawings")
ax_a.set_ylabel("Estimation of sigma_star")
#ax_a.legend(["MC","Algo"],loc='best')

fig_b = plt.figure()
ax_b = fig_b.gca()
ax_b.plot(f[:,0],f[:,2],c='orange')
ax_b.plot(f[:,0],f[:,4],c='black')
ax_b.axhline(y=sigma_star,color='red')
ax_b.set_xlabel("Number of drawings")
ax_b.set_ylabel("Estimation of sigma_star")
#ax_b.legend(["MC","Algo"],loc='best')


# Save figures
fig_a.savefig(target + "_a.png",**save_params)
fig_b.savefig(target + "_b.png",**save_params)
