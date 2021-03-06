import numpy as np
import matplotlib.pyplot as plt

#--------------------------------------------
# ATTENTION this plot script should only 
# be used for spatially variable epsilon 
# and/or B 
#--------------------------------------------

# Initialize the arrays
x = []
y = []

# Set up the double x-axis plot figure
fig = plt.figure(figsize=[7,6])
#fig.set_size_inches(4.,4.)
ax1 = fig.add_subplot(1,1,1)
ax2 = ax1.twiny()

# Get the data and start adding to the plot
# Read the file and strip it of any extra whitespace
file = open('output_file.txt','r').read().strip()
# Make one list containing each value from the .txt file
s = file.split("\n")
# Initialize a counter that will create our z-array later
j = -3
# Create an array of the values from the .txt files
# As well as re-calculating the z_arr
for i in s:
    if (j == -3):
        # Epsilon for the problem will always be the first value in the .txt file 
        epsilon_l = float(i)
    elif (j == -2):
        epsilon_r = float(i)
    elif (j == -1):
        del_eps = float(i)
    elif (i == 'break'):
        # Define/Re-define plotting arrays
        x_arr = np.array(x)
        z_arr = x_arr/(len(x_arr)-1.)
        y_arr = np.array(y)
        ax1.semilogy(z_arr,y_arr,color='seagreen')
        # Reset x and y lists
        x = []
        y = []
        # Reset the placement counter considering +1 is added this iteration
        j = -1
    else:
        x.append(float(j))
        y.append(float(i))
    j += 1

# Create a dummy plot for the epsilon label at the top of the plot 
#ax2.plot(np.linspace(epsilon_l,epsilon_r,70),np.ones(70),color='w')
ax2.set_xlim(epsilon_l,epsilon_r)
ax2.set_xlabel('epsilon')
# Should have plotted all the lines together
ax1.set_title("S/B vs z, varying epsilon")
#plt.title("ALI Scheme: S/B vs z, epsilon=%1.3f" %epsilon)
ax1.set_xlabel("z")
ax1.set_ylabel("log(S/B)")
ax1.set_xlim(0.0,1.0)
ax1.set_ylim(0.1,1.1)
plt.show()
