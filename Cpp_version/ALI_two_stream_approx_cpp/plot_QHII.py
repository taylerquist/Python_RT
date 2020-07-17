import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

y_RK = []

# Read the file and strip it of any extra whitespace
RK4_file = open('QHII_out.txt','r').read().strip()
info_file = open('info_outfile.txt').read().strip()
# Make one list containing each value from the .txt file
s_RK = RK4_file.split("\n")
info = info_file.split("\n")
# Grab all the necessary info about the grid, number of iterations, and conversion to kpc
grid_size = 1000

cell_size_kpc = 622.08123/float(grid_size)
#cell_size_kpc = float(cell_size_kpc)

iters = 0

for i in s_RK:
    if (i == 'break'):
        iters += 1
        
bound = iters - 1
iters = 0        
# Set the plot colors
#ax = plt.axes()
#ax.set_prop_cycle('color',[plt.cm.viridis(i) for i in np.linspace(0, 1, iters)])

for i in s_RK:
    if (i == 'break'):
        if (iters == bound):
            y_arr_RK = np.array(y_RK)
            x = np.linspace(0,grid_size,len(y_arr_RK))
            x = x*cell_size_kpc
            # Save the array for an output file
            q_arr = y_arr_RK
            #plt.semilogy(x,y_arr_RK)
            plt.plot(x,y_arr_RK)
            #plt.title("log(QHII) vs distance (from galaxy) over time")
            #plt.xlabel("distance kpc")
            #plt.ylabel("log(QHII)")
            #plt.show()
        # Reset y list
        y_RK = []
        iters += 1
    else:
        y_RK.append(float(i))

# Should have plotted all the lines together 
plt.title("log(QHII) vs distance (from galaxy) over time")
plt.xlabel(r'distance [$\mathrm{kpc}$]')
plt.ylabel(r'log(QHII)')
#plt.xlim(100,)
plt.show()

# If the file needs to be created or overwritten
#outfile = open("QHII_results.txt","w")
# If the file is being appended
outfile = open("QHII_results.txt","a")
for val in q_arr:
    outfile.write(str(val))
    outfile.write("\n")
outfile.write("break")
outfile.write("\n")
outfile.close()
