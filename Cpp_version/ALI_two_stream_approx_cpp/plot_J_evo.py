import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

j_list = []

# Read the file and strip it of any extra whitespace
j_file = open('J_out.txt','r').read().strip()
info_file = open('info_outfile.txt').read().strip()
# Make one list containing each value from the .txt file
j_s = j_file.split("\n")
info = info_file.split("\n")
# Grab all the necessary info about the grid, number of iterations, and conversion to kpc
grid_size = 1000

cell_size_kpc = 622.08123/float(grid_size)

iters = 0

for i in j_s:
    if (i == 'break'):
        iters += 1

ax = plt.axes()
ax.set_prop_cycle('color',[plt.cm.viridis(i) for i in np.linspace(0, 1, iters)])

for i in j_s:
    if (i == 'break'):
        j_arr = np.array(j_list)
        x = np.linspace(0,grid_size,len(j_arr))
        x = x*cell_size_kpc
        plt.semilogy(x,j_arr)
        # Reset x and y lists and counter
        j_list = []
    else:
        j_list.append(float(i))

# Should have plotted all the lines together 
plt.title(r'Average Specific Intensity Evolution, grid resolution 1000')
plt.ylabel(r'$log(J) \ $ [$\mathrm{erg \ cm^{-2} \ s^{-1} \ Hz^{-1} ster^{-1}}$]')
plt.xlabel(r'$x$ [$\mathrm{kpc}$]')
#plt.ylim(1e-13,1e-3)
plt.show()

