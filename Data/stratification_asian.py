import sys
import os
import matplotlib.pyplot as plt
import numpy as np

directory = "/home/jonathan/Programmation/C++/Projet_C++_S2"

# Get target path
target = directory + "/Images/stratification_asian"

# Image saving parameters
save_params = {"bbox_inches":'tight', "pad_inches":0.}

#%% Plot coordinates of u
# Retrieve data
u_call = np.genfromtxt(fname=directory + "/Data/stratification_asian_call_u.dat", delimiter=" ")
u_put = np.genfromtxt(fname=directory + "/Data/stratification_asian_put_u.dat", delimiter=" ")

# Make figures
fig = plt.figure()
ax = fig.gca()
ax.plot(range(1,len(u_call)+1), u_call,
        range(1,len(u_put)+1), u_put)
        
ax.set_xlabel("m")
ax.set_ylabel("u_m")
ax.legend(["Call","Put"],loc='best')

# Save figures
fig.savefig(target + "_u.png",**save_params)

#%% Plot stds and mean estimation on strata

f_call = np.genfromtxt(fname=directory + "/Data/stratification_asian_call_stds.dat", delimiter=" ")
f_put = np.genfromtxt(fname=directory + "/Data/stratification_asian_put_stds.dat", delimiter=" ")

# Make figures
fig_stds_call = plt.figure()
ax_stds_call = fig_stds_call.gca()
ax_stds_call.plot(range(1,len(f_call)+1), f_call[:,0])
#ax_stds_call.legend(["Method a)","Method b)"],loc='best')

fig_mean_call = plt.figure()
ax_mean_call = fig_mean_call.gca()
ax_mean_call.plot(range(1,len(f_call)+1), f_call[:,2])
#ax_mean_call.legend(["Method a)","Method b)"],loc='best')

fig_stds_put = plt.figure()
ax_stds_put = fig_stds_put.gca()
ax_stds_put.plot(range(1,len(f_put)+1), f_put[:,0])
#ax_stds_put.legend(["Method a)","Method b)"],loc='best')

fig_mean_put = plt.figure()
ax_mean_put = fig_mean_put.gca()
ax_mean_put.plot(range(1,len(f_put)+1), f_put[:,2])
#ax_mean_put.legend(["Method a)","Method b)"],loc='best')

# Save figures
fig_stds_call.savefig(target + "_call_stds.png",**save_params)
fig_mean_call.savefig(target + "_call_means.png",**save_params)
fig_stds_put.savefig(target + "_put_stds.png",**save_params)
fig_mean_put.savefig(target + "_put_means.png",**save_params)




















