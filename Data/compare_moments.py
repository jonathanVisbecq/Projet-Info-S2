import sys
import os
import matplotlib.pyplot as plt
import numpy as np

directory = "/home/jonathan/Programmation/C++/Projet_C++_S2"

# Get target path
target = directory + "/Images/compare_moments"

# Image saving parameters
save_params = {"bbox_inches":'tight', "pad_inches":0.}

for func in xrange(1,6):    

    path = directory + "/Data/compare_moments_func" + str(func) + ".dat"

    if not os.path.isfile(path):
        continue     
    
    f = np.genfromtxt(fname=path, delimiter=" ")
    
    # Make figures
    # Second moment
    fig2 = plt.figure()
    ax2 = fig2.gca()
    ax2.plot(f[:,0],f[:,1:4])
    ax2.set_xlabel("M*n - log scale")
    ax2.set_xscale('log')
    ax2.set_ylabel("Variance - log scale")
    ax2.set_yscale('log')
    ax2.legend(["Monte Carlo","Shifted QMC","Random-start Halton"],loc='best')
    
    # Third moment
    fig3 = plt.figure()
    ax3 = fig3.gca()
    ax3.plot(f[:,0],f[:,4:7])
    ax3.set_xlabel("M*n - log scale")
    ax3.set_xscale('log')
    ax3.set_ylabel("Normalized Skewness - log scale")
    ax3.set_yscale('log')
    ax3.legend(["Monte Carlo","Shifted QMC","Random-start Halton"],loc='best')

    # Save figures
    fig2.savefig(target + "_func" + str(func) + "_M2.png",**save_params)
    fig3.savefig(target + "_func" + str(func) + "_M3.png",**save_params)










