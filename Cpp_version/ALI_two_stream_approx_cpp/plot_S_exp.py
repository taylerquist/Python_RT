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
grid_size = 1000

cell_size_kpc = 622.08123/float(grid_size)

# Find out how many iterations the code did and use it to set the color scheme
iters = 0
for i in s_RK:
    if (i == 'break'):
        iters += 1

# Set the plot colors
ax = plt.axes()
ax.set_prop_cycle('color',[plt.cm.plasma(i) for i in np.linspace(0, 1, iters)])

for i in s_RK:
    if (i == 'break'):
        y_arr_RK = np.array(y_RK)
        x = np.linspace(0,grid_size,len(y_arr_RK))
        x = x*cell_size_kpc
        plt.semilogy(x,y_arr_RK)
        # Reset y list and counter
        y_RK = []
    else:
        y_RK.append(float(i))

# Should have plotted all the lines together 
plt.title("log(S) vs distance (from galaxy) over time")
plt.xlabel(r'distance [$\mathrm{kpc}$]')
plt.ylabel(r'log(S)')
#plt.xlim(-0.1,1.1)
plt.show()

