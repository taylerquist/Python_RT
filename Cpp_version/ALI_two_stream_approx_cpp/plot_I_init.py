import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

#y = []
#y_act = []

# Read the file and strip it of any extra whitespace
#file = open('I_luminosity.txt','r').read().strip()
info_file = open('info_outfile.txt').read().strip()
# Make one list containing each value from the .txt file
#s = file.split("\n")
#count = 1.
#for i in s:
    # Read value in from the file
#    print(i)
#    y.append(j)
    # Create the actual solution to 1/r**2
#    val = 1./(count**2.)
#    y_act.append(val)
#    count += 1.
    
#print(y)
#print('---------')
# Read in some important info from the run to set x and convert to kpc
info = info_file.split("\n")

grid_size = info[0]
grid_size = int(grid_size)

cell_size_kpc = info[-1]
cell_size_kpc = float(cell_size_kpc)

#print(cell_size_kpc)

# Make the x axis
x = np.linspace(0.,(grid_size-1),grid_size)
x = x*cell_size_kpc
y = np.ones(grid_size)
y[1:] = 1/(x[1:])**2.
#print(x)
#print(y_act)

# Plot it
#plt.plot(x,y,color='m')
plt.semilogy(x,y,color='dodgerblue',linewidth=2, label='I_galaxy')
#plt.semilogy(x,y_act,color='b',linestyle='--',label='1/r**2')
plt.title(r"Initial $I^{+}$: the specific intensity for the galaxy")
plt.xlabel(r'$x$ [$\mathrm{kpc}$]')
plt.ylabel(r'$I \ $ [$\mathrm{erg \ cm^{-2} \ s^{-1} \ Hz^{-1} ster^{-1}}$]')
#plt.xlim(0.0,1.0)
#plt.legend()
plt.show()

