import numpy as np
import matplotlib.pyplot as plt

# Initialize the arrays
del_tau = []

# Read the file and strip it of any extra whitespace
file_del_tau = open('del_tau_out.txt','r').read().strip()
f_dt = file_del_tau.split("\n")
n = len(f_dt)

# Set important values
grid_size = 1000
cell_size_kpc = 9.01567

# Set the plot colors
ax = plt.axes()
ax.set_prop_cycle('color',[plt.cm.viridis(i) for i in np.linspace(0, 1, 11)])

# Create an array of the values from the .txt files and plot it 
for i in f_dt:
    if (i == 'break'):
        del_tau_arr = np.array(del_tau) 
        x = np.linspace(0,grid_size,len(del_tau_arr))
        x = x*cell_size_kpc
        plt.plot(x,del_tau_arr)
        plt.title(r'Delta Tau vs Radius: grid resolution 1000')
        plt.xlabel(r'Radius [$\mathrm{kpc}$]')
        plt.ylabel(r'Delta Tau')
        plt.show()
        x = []
        del_tau = []
    else:
        del_tau.append(float(i))  # Save it
    
# Show the plot that has been accumulating
#plt.title(r'Delta Tau vs Radius: grid resolution 1000')
#plt.xlabel(r'Radius [$\mathrm{kpc}$]')
#plt.ylabel(r'Delta Tau')
#plt.show()
