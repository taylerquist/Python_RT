import numpy as np
import matplotlib.pyplot as plt

# Initialize the arrays
x = []
y = []

# Read the file and strip it of any extra whitespace
file = open('output_file.txt','r').read().strip()
# Make one list containing each value from the .txt file
s = file.split("\n")
# Initialize a counter that will create our z-array later
j = -1
# Create an array of the values from the .txt files
# As well as re-calculating the z_arr
for i in s:
    if (j == -1):
        # Epsilon for the problem will always be the first value in the .txt file 
        epsilon = float(i)
    elif (i == 'break'):
        # Define/Re-define plotting arrays
        x_arr = np.array(x)
        z_arr = x_arr/(len(x_arr)-1.)
        y_arr = np.array(y)
        plt.plot(z_arr,y_arr,color='k')
        # Reset x and y lists
        x = []
        y = []
        # Reset the placement counter considering +1 is added this iteration
        j = -1
    else:
        x.append(float(j))
        y.append(float(i))
    j += 1

# Should have plotted all the lines together
plt.title("S/B vs z, epsilon=%1.3f" %epsilon)
plt.xlabel("z")
plt.ylabel("S/B")
plt.ylim(0.1,1.1)
plt.xlim(0.0,1.0)
plt.show()
