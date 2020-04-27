import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

y_RK = []

# Read the file and strip it of any extra whitespace
RK4_file = open('ion_frac_outfile.txt','r').read().strip()
info_file = open('info_outfile.txt').read().strip()
# Make one list containing each value from the .txt file
s_RK = RK4_file.split("\n")
info = info_file.split("\n")
# Grab all the necessary info about the grid, number of iterations, and conversion to kpc
grid_size = info[0]
grid_size = int(grid_size)

n = info[1]
n = int(n)

cell_size_kpc = info[-1]
cell_size_kpc = float(cell_size_kpc)

#x = np.linspace(0,grid_size,grid_size)
#x = x*cell_size_kpc
#print(x)
# Set the plot colors
ax = plt.axes()
ax.set_prop_cycle('color',[plt.cm.viridis(i) for i in np.linspace(0, 1, n)])

# Set a counter for x-axis
index = 0
count = 0.
h = 0.1

#z_arr = np.linspace(0.,1.,70)

for i in s_RK:
    if (i == 'break'):
        #x_arr = np.array(x)
        #z_arr = x_arr/(len(x_arr)-1.)
        #z_arr = x_arr/70.
        y_arr_RK = np.array(y_RK)
        x = np.linspace(0,grid_size,len(y_arr_RK))
        x = x*cell_size_kpc
        #plt.semilogy(x,y_arr_RK)
        plt.plot(x,y_arr_RK)
        # Reset x and y lists and counter
        x = []
        y_RK = []
        count = 0.
    else:
        y_RK.append(float(i))
    count = count + h
    index += 1

# Should have plotted all the lines together 
plt.title("log(QHII) vs distance (from galaxy) over time")
plt.xlabel("distance kpc")
plt.ylabel("log(QHII)")
#plt.xlim(0.0,1.0)
plt.show()

