import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

y = []
y_act = []

# Read the file and strip it of any extra whitespace
file = open('I_luminosity.txt','r').read().strip()
info_file = open('info_outfile.txt').read().strip()
# Make one list containing each value from the .txt file
s = file.split("\n")
count = 1.
for i in s:
    # Read value in from the file
    y.append(float(i))
    # Create the actual solution to 1/r**2
    val = 1./(count**2.)
    y_act.append(val)
    count += 1.
    
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
#x = x*cell_size_kpc
#print('y_act')
#print(y_act)

# Plot it
#plt.plot(x,y,color='m')
plt.semilogy(x,y,color='orange',linewidth=3, label='Cpp')
plt.semilogy(x,y_act,color='b',linestyle='--',label='1/r**2')
plt.title("Initial I+: the specific intensity for the galaxy")
plt.xlabel("grid point number")
plt.ylabel("log(I): log(erg/(cm**2 s Hz rad**2))")
#plt.xlim(0.0,1.0)
plt.legend()
plt.show()

