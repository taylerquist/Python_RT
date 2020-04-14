import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

x = []
y_act = []
y_RK = []
n = 50

# Read the file and strip it of any extra whitespace
actual_file = open('true_file.txt','r').read().strip()
RK4_file = open('RK_file.txt','r').read().strip()
# Make one list containing each value from the .txt file
s_act = actual_file.split("\n")
s_RK = RK4_file.split("\n")
# Set the plot colors
ax = plt.axes()
ax.set_prop_cycle('color',[plt.cm.viridis(i) for i in np.linspace(0, 1, n)])

# Set a counter for x-axis
index = 0
count = 0.
h = 0.1

for i in s_RK:
    if (i == 'break'):
        x_arr = np.array(x)
        z_arr = x_arr/(len(x_arr)-1.)
        y_arr_act = np.array(y_act)
        y_arr_RK = np.array(y_RK)
        #plt.plot(z_arr,y_arr_act,'-.',color='m')
        plt.plot(z_arr,y_arr_RK)
        plt.plot(z_arr,y_arr_act,'-.',color='m')
        # Reset x and y lists and counter
        x = []
        y_act = []
        y_RK = []
        count = 0.
    else:
        x.append(float(count))
        j = s_act[index]
        y_act.append(float(j))
        y_RK.append(float(i))
    count = count + h
    index += 1

# Should have plotted all the lines together 
plt.title("F vs z (grid depth) over time")
plt.xlabel("z")
plt.ylabel("F")
#plt.xlim(0.0,1.0)
plt.show()

