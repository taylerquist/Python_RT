import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

x = []
y_RK = []
n = 40

# Read the file and strip it of any extra whitespace
RK4_file = open('ion_frac_outfile.txt','r').read().strip()
# Make one list containing each value from the .txt file
s_RK = RK4_file.split("\n")
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
        x = len(y_RK)
        z_arr = np.linspace(0.,1.,x)
        #plt.plot(z_arr,y_arr_act,'-.',color='m')
        #plt.plot(z_arr,y_arr_RK)
        plt.semilogy(z_arr,y_arr_RK)
        # Reset x and y lists and counter
        x = []
        y_RK = []
        count = 0.
    else:
        x.append(float(count))
        y_RK.append(float(i))
    count = count + h
    index += 1

# Should have plotted all the lines together 
plt.title("QHII vs z (grid depth) over time")
plt.xlabel("z")
plt.ylabel("QHII")
plt.xlim(0.0,1.0)
plt.show()

