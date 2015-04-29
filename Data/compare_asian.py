import os
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from sklearn import linear_model

directory = "/home/jonathan/Programmation/C++/Projet_C++_S2"

# Get target path
target = directory + "/Images/compare_asian"

# Image saving parameters
save_params = {"bbox_inches":'tight', "pad_inches":0.}

# Create linear regression object
regr = linear_model.LinearRegression()

#%% Plot stratification
# Retrieve data

for d in [16,64]:
    for K in [45,55]:

        strat = np.genfromtxt(fname=directory + "/Data/compare_asian_strat_" + str(d) + "_" + str(K) + ".dat", delimiter=' ')
        rqmc = np.genfromtxt(fname=directory + "/Data/compare_asian_rqmc_" + str(d) + "_" + str(K) + ".dat", delimiter=' ')
                
        x = sp.linspace(1,1e5,100) / 100.
        # Make figures
        fig = plt.figure()
        ax = fig.gca()
        ax.plot(strat[:,1], strat[:,2], c='red')
        
        regr.fit(np.transpose(np.asarray([np.log(strat[:,1])])),
                 np.transpose(np.asarray([np.log(strat[:,2])])))
        reg_strat = regr.coef_[0,0]
        plt_strat, = ax.plot(x,np.exp(regr.coef_[0,0]*np.log(x)+regr.intercept_[0]),'--',c="red")
        plt_strat.set_label("{0:.3}".format(reg_strat))
        
        ax.plot(rqmc[:,2], rqmc[:,4], c='blue')
        
        regr.fit(np.transpose(np.asarray([np.log(rqmc[:,2])])),
                 np.transpose(np.asarray([np.log(rqmc[:,4])])))
        reg_sqmc = regr.coef_[0,0]
        plt_sqmc, = ax.plot(x,np.exp(regr.coef_[0,0]*np.log(x)+regr.intercept_[0]),'--',c="blue")
        plt_sqmc.set_label("{0:.3}".format(reg_sqmc))
        
        ax.plot(rqmc[:,3], rqmc[:,5], c='green')
        
        regr.fit(np.transpose(np.asarray([np.log(rqmc[:,3])])),
		 np.transpose(np.asarray([np.log(rqmc[:,5])])))
        reg_rdStart = regr.coef_[0,0]
        plt_rdStart, = ax.plot(x,np.exp(regr.coef_[0,0]*np.log(x)+regr.intercept_[0]),'--',c="green")
        plt_rdStart.set_label("{0:.3}".format(reg_rdStart))
        
        ax.set_ylabel("CI - log scale")
        ax.set_yscale('log')
        ax.set_xlabel("Times (s) - log scale")
        ax.set_xscale('log')
        
        ax.legend()
        
        # Save figures
        fig.savefig(target + "_" + str(d) + "_" + str(K) + ".png",**save_params)
        
    
    
    
    
    
    
    
    
    
    
    
    



