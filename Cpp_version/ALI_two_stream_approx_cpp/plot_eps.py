import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

y_eps = []

# Read the file and strip it of any extra whitespace
eps_file = open('eps_out.txt','r').read().strip()
info_file = open('info_outfile.txt').read().strip()
# Make one list containing each value from the .txt file
s_eps = eps_file.split("\n")
info = info_file.split("\n")
# Grab all the necessary info about the grid, number of iterations, and conversion to kpc
grid_size = 1000

cell_size_kpc = 622.08123/float(grid_size)
#cell_size_kpc = float(cell_size_kpc)

iters = 0

for i in s_eps:
    if (i == 'break'):
        iters += 1
        
print(iters)
        
# Set the plot colors
ax = plt.axes()
ax.set_prop_cycle('color',[plt.cm.viridis(i) for i in np.linspace(0, 1, iters)])

for i in s_eps:
    if (i == 'break'):
        y_arr_eps = np.array(y_eps)
        x = np.linspace(0,grid_size,len(y_arr_eps))
        x = x*cell_size_kpc
        #plt.semilogy(x,y_arr_eps)
        plt.plot(x,y_arr_eps)
        #plt.title("log(QHII) vs distance (from galaxy) over time")
        #plt.xlabel("distance kpc")
        #plt.ylabel("log(QHII)")
        #plt.show()
        # Reset y list
        y_eps = []
    else:
        y_eps.append(float(i))

# Should have plotted all the lines together 
plt.title(r'log($\epsilon$) vs distance (from galaxy) over time')
plt.xlabel(r'distance [$\mathrm{kpc}$]')
plt.ylabel(r'log($\epsilon$)')
#plt.xlim(100,)
plt.show()

