import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

y_RK = []

# Read the file and strip it of any extra whitespace
RK4_file = open('S_outfile.txt','r').read().strip()
info_file = open('info_outfile.txt').read().strip()
# Make one list containing each value from the .txt file
s_RK = RK4_file.split("\n")
info = info_file.split("\n")
# Grab the number of iterations S went through to divide plotting colors
grid_size = info[0]
grid_size = int(grid_size)

n = info[1]
n = int(n)

cell_size_kpc = info[-1]
cell_size_kpc = float(cell_size_kpc)
# Set the x array: the kpc distance
x = np.linspace(0.,1.,grid_size)
x = x*cell_size_kpc
# Set the plot colors
ax = plt.axes()
ax.set_prop_cycle('color',[plt.cm.plasma(i) for i in np.linspace(0, 1, n)])

# Set a counter for x-axis
index = 0
count = 0.

#z_arr = np.linspace(0.,1.,70)

for i in s_RK:
    if (i == 'break'):
        #x_arr = np.array(x)
        #z_arr = x_arr/(len(x_arr)-1.)
        #z_arr = x_arr/70.
        y_arr_RK = np.array(y_RK)
        x = np.linspace(0,grid_size,len(y_arr_RK))
        x = x*cell_size_kpc
        plt.semilogy(x,y_arr_RK)
        # Reset x and y lists and counter
        #x = []
        y_RK = []
        count = 0.
    else:
        #x.append(float(count))
        y_RK.append(float(i))
    count += 1.
    index += 1

# Should have plotted all the lines together 
plt.title("log(S) vs distance (from galaxy) over time")
plt.xlabel("distance kpc")
plt.ylabel("log(S)")
#plt.xlim(-0.1,1.1)
plt.show()

