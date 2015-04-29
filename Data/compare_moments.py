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
    
    fig_sqmc = plt.figure()
    ax_sqmc = fig_sqmc.gca()
    ax_sqmc.plot(f[:,0],f[:,2]**(3./2.))
    ax_sqmc.plot(f[:,0],f[:,5])
    ax_sqmc.set_xlabel("M*n")
    ax_sqmc.set_xscale('log')
    ax_sqmc.set_yscale('log')
    ax_sqmc.legend(["Std^3","Third Moment"],loc='best')
    
    fig_sqmc.savefig(target + "_func" + str(func) + "_sqmc.png",**save_params)

    fig_rdStart = plt.figure()
    ax_rdStart = fig_rdStart.gca()
    ax_rdStart.plot(f[:,0],f[:,3]**(3./2.))
    ax_rdStart.plot(f[:,0],f[:,6])
    ax_rdStart.set_xlabel("M*n")
    ax_rdStart.set_xscale('log')
    ax_rdStart.set_yscale('log')
    ax_rdStart.legend(["Std^3","Third Moment"],loc='best')
    
    fig_rdStart.savefig(target + "_func" + str(func) + "_rdStart.png",**save_params)










