import numpy as np
import matplotlib.pyplot as plt

# Read the file and strip it of any extra whitespace
  # Analytic file
a_file = open('s_sphere_analytic.txt','r').read().strip()
  # Computational file
c_file = open('s_sphere_code.txt','r').read().strip()

# Split the file into its elements
a = a_file.split("\n")
c = c_file.split("\n")

# Read the file values and save them to a list
analytic = []
comp = []
count = 0

for i in a:
    # Save the value as a float to the list
    analytic.append(float(i))
    comp.append(float(c[count]))
    count += 1

# Certain analytic values are out of the code's resolution
a_cut = []
a_cut.append(analytic[0])
a_cut.append(analytic[1])
a_cut.append(analytic[-1])

# Cut out the points that are out of the code's resolution
comp = comp[2:-1]

# Create the density range
rho = [1.,1.e-1,1.e-2,1.e-3,1.e-4,1.e-5,1.e-6,1.e-7,1.e-8]
rho_comp = rho[2:-1]
rho_cut = [1.,1.e-1,1.e-8]

# Plot it
plt.semilogx(rho_comp,comp,color='lime',marker='o',linestyle='',label='Code')
plt.semilogx(rho,analytic,color='blue',marker='o',linestyle='',label='Analytic')
plt.semilogx(rho_cut,a_cut,color='indigo',marker='o',linestyle='',label='Analytic out of range')
plt.title(r'Str$\ddot{o}$mgren Sphere Radius vs Density')
plt.xlabel(r'Density [$\mathrm{cm^{-3}}$]')
plt.ylabel(r'Radius [$\mathrm{kpc}$]')
plt.legend()
plt.show()

