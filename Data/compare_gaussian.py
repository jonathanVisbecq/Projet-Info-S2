import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model

directory = "/home/jonathan/Programmation/C++/Projet_C++_S2"

# Get target path
target = directory + "/Images/compare_gaussian"

# Image saving parameters
save_params = {"bbox_inches":'tight', "pad_inches":0.}

# Create linear regression object
regr = linear_model.LinearRegression()

#%% Plot stratification
# Retrieve data
strat = np.genfromtxt(fname=directory + "/Data/compare_gaussian_strat.dat", delimiter=" ")
rqmc = np.genfromtxt(fname=directory + "/Data/compare_gaussian_rqmc.dat", delimiter=" ")

# Make figures
fig = plt.figure()
ax = fig.gca()
ax.plot(strat[:,1], strat[:,2], c='red')

regr.fit(np.transpose(np.asarray([np.log(strat[:,1])])), np.transpose(np.asarray([np.log(strat[:,2])])))
reg_strat = regr.coef_[0,0]

ax.plot(rqmc[:,2], rqmc[:,4], c='blue')

regr.fit(np.transpose(np.asarray([np.log(rqmc[:,2])])), 
         np.transpose(np.asarray([np.log(rqmc[:,4])])))
reg_sqmc = regr.coef_[0,0]

ax.plot(rqmc[:,3], rqmc[:,5], c='green')

regr.fit(np.transpose(np.asarray([np.log(rqmc[:,3])])),
         np.transpose(np.asarray([np.log(rqmc[:,5])])))
reg_rdStart = regr.coef_[0,0]

ax.set_ylabel("CI - log scale")
ax.set_yscale('log')
ax.set_xlabel("Times (s) - log scale")
ax.set_xscale('log')

ax.legend(["{0:.3}".format(reg_strat),
           "{0:.3}".format(reg_sqmc),
           "{0:.3}".format(reg_rdStart)],loc='best')

# Save figures
fig.savefig(target + ".png",**save_params)
















