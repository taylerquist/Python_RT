import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

I_p = []
I_m = []

# Read the file and strip it of any extra whitespace
Ip_file = open('I_p_output_IC.txt','r').read().strip()
Im_file = open('I_m_output_IC.txt','r').read().strip()
info_file = open('info_outfile.txt').read().strip()
# Make one list containing each value from the .txt file
Ip = Ip_file.split("\n")
Im = Im_file.split("\n")
#count = 0
for i in range(0,len(Ip)):
    # Read values in from the files
    I_p.append(float(Ip[i]))
    I_m.append(float(Im[i]))

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
#y = np.ones(grid_size)
#y[1:] = 1/(x[1:])**2.
#print(x)
#print(y_act)

# Plot it
fig,(ax1,ax2) = plt.subplots(2, sharex=True)
fig.suptitle('Steady State Specific Intensity: $I$')
ax1.plot(x,I_p,color='dodgerblue',label='Galaxy')
ax2.plot(x,I_m,color='indigo',label='UV Background')
ax2.set_xlabel(r'$x$ [$\mathrm{kpc}$]')
fig.text(0.01, 0.5, r'$I \ $ [$\mathrm{erg \ cm^{-2} \ s^{-1} \ Hz^{-1} ster^{-1}}$]', va='center', rotation='vertical')
#ax1.set_ylabel("I: erg/(cm**2 s Hz rad**2)")
#ax2.set_ylabel("I: erg/(cm**2 s Hz rad**2)")
#ax2.set_ylim(0.0,1.0)
#plt.xlim(0.0,1.0)
fig.set_figheight(16)
fig.set_figwidth(8)
ax1.legend()
ax2.legend()
plt.show()

