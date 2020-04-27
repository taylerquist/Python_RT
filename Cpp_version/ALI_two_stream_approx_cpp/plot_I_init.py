import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

y = []

# Read the file and strip it of any extra whitespace
file = open('I_luminosity.txt','r').read().strip()
info_file = open('info_outfile.txt').read().strip()
# Make one list containing each value from the .txt file
s = file.split("\n")
for i in s:
    y.append(float(i))
# Read in some important info from the run to set x and convert to kpc
info = info_file.split("\n")

grid_size = info[0]
grid_size = int(grid_size)

cell_size_kpc = info[-1]
cell_size_kpc = float(cell_size_kpc)
# Make the x axis
x = np.linspace(0.,1.,grid_size)
x = x*cell_size_kpc
# Plot it
#plt.plot(x,y,color='m')
plt.semilogy(x,y,color='m')
plt.title("Initial I+: the specific intensity for the galaxy")
plt.xlabel("kpc")
plt.ylabel("log(I): log(erg/(cm**2 s Hz rad**2))")
#plt.xlim(0.0,1.0)
plt.show()

