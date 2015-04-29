import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model

directory = "/home/jonathan/Programmation/C++/Projet_C++_S2"

# Get target path
target = directory + "/Images/compare_ciAndTime_MCvsRQMC"

# Read data
f = np.genfromtxt(fname=directory + "/Data/compare_ciAndTime_MCvsRQMC.dat", delimiter=" ")

# Create linear regression object
regr = linear_model.LinearRegression()




# Make figures
fig1 = plt.figure()
ax1 = fig1.gca()
ax1.plot(f[:,0],f[:,1:4])
ax1.set_xlabel("M*n")
ax1.set_xscale('log')
ax1.set_ylabel("Half confidence interval")
ax1.set_yscale('log')

regr.fit(np.transpose(np.asarray([np.log(f[:,0])])), np.transpose(np.asarray([np.log(f[:,1])])))
reg_mc = regr.coef_[0,0]
regr.fit(np.transpose(np.asarray([np.log(f[:,0])])), np.transpose(np.asarray([np.log(f[:,2])])))
reg_sqmc = regr.coef_[0,0]
regr.fit(np.transpose(np.asarray([np.log(f[:,0])])), np.transpose(np.asarray([np.log(f[:,3])])))
reg_rdStart = regr.coef_[0,0]

ax1.legend(["Monte Carlo: {0:.3}".format(reg_mc),
            "Shifted QMC: {0:.3}".format(reg_sqmc),
            "Rd-start Halton: {0:.3}".format(reg_rdStart)],loc='best')

fig2 = plt.figure()
ax2 = fig2.gca()
ax2.plot(f[:,0],f[:,4:7])
ax2.set_xlabel("M*n")
ax2.set_xscale('log')
ax2.set_ylabel("Time (s)")
ax2.set_yscale('log')

idx = f[:,0]>20000
regr.fit(np.transpose(np.asarray([np.log(f[idx,0])])), np.transpose(np.asarray([np.log(f[idx,4])])))
reg_mc = regr.coef_[0,0]
regr.fit(np.transpose(np.asarray([np.log(f[idx,0])])), np.transpose(np.asarray([np.log(f[idx,5])])))
reg_sqmc = regr.coef_[0,0]
regr.fit(np.transpose(np.asarray([np.log(f[idx,0])])), np.transpose(np.asarray([np.log(f[idx,6])])))
reg_rdStart = regr.coef_[0,0]

ax2.legend(["Monte Carlo: {0:.3}".format(reg_mc),
            "Shifted QMC: {0:.3}".format(reg_sqmc),
            "Rd-start Halton: {0:.3}".format(reg_rdStart)],loc='best')

fig3 = plt.figure()
ax3 = fig3.gca()
ax3.plot(f[:,4],f[:,1])
ax3.plot(f[:,5],f[:,2])
ax3.plot(f[:,6],f[:,3])
ax3.set_xlabel("CI - log scale")
ax3.set_xscale('log')
ax3.set_ylabel("Time (s) - log scale")
ax3.set_yscale('log')
ax3.legend(["Monte Carlo","Shifted QMC","Random-start Halton"],loc='best')

# Save figures
save_params = {"bbox_inches":'tight', "pad_inches":0.}

fig1.savefig(target + "__CI.png",**save_params)
fig2.savefig(target + "__time.png",**save_params)
fig3.savefig(target + ".png",**save_params)
