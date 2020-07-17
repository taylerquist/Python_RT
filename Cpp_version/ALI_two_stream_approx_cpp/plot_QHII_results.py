import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

q_list = []

# Read the file and strip it of any extra whitespace
q_file = open('QHII_results.txt','r').read().strip()
info_file = open('info_outfile.txt').read().strip()
# Make one list containing each value from the .txt file
q_s = q_file.split("\n")
info = info_file.split("\n")
# Grab all the necessary info about the grid, number of iterations, and conversion to kpc
grid_size = 1000

cell_size_kpc = 622.08123/float(grid_size)

count = 0

#ax = plt.axes()
#ax.set_prop_cycle('color',[plt.cm.viridis(i) for i in np.linspace(0, 1, iters)])

for i in q_s:
    if (i == 'break'):
        q_arr = np.array(q_list)
        x = np.linspace(0,grid_size,len(q_arr))
        x = x*cell_size_kpc
        if count == 0:
            #plt.plot(x,q_arr,color='maroon',label='f_esc=1.0') # If varying f_esc
            plt.plot(x,q_arr,color='maroon',label='I- = 1e-23') # If varying I-
        elif count == 1:
            #plt.plot(x,q_arr,color='orangered',linestyle='-.',label='f_esc=0.1') # If varying f_esc
            plt.plot(x,q_arr,color='orangered',linestyle='-.',label='I- = 1e-24') # If varying I-
        elif count == 2:
            #plt.plot(x,q_arr,color='gold',linestyle='--',label='f_esc=0.01') # If varying f_esc
            plt.plot(x,q_arr,color='gold',linestyle='--',label='I- = 1e-25') # If varying I-
        # Reset x and y lists and counter
        q_list = []
        count += 1
    else:
        q_list.append(float(i))

# Should have plotted all the lines together 
plt.title(r'Converged QHII, 0.001$\rho$')
plt.ylabel(r'$QHII')
plt.xlabel(r'$x$ [$\mathrm{kpc}$]')
plt.ylim(0,1.1)
#plt.xlim(90.,180.)
plt.legend()
plt.show()
