import numpy as np
import matplotlib.pyplot as plt

# Read the file and strip it of any extra whitespace
  # Analytic file
a_file = open('s_sphere_analytic.txt','r').read().strip()
  # Computational file
c1_file = open('s_sphere_code_70.txt','r').read().strip()
c2_file = open('s_sphere_code_200.txt','r').read().strip()
c3_file = open('s_sphere_code_500.txt','r').read().strip()
c4_file = open('s_sphere_code_1000.txt','r').read().strip()

# Split the file into its elements
a = a_file.split("\n")
c1 = c1_file.split("\n")
c2 = c2_file.split("\n")
c3 = c3_file.split("\n")
c4 = c4_file.split("\n")

# Read the file values and save them to a list
analytic = []
comp1 = []
comp2 = []
comp3 = []
comp4 = []

count = 0

for i in a:
    # Save the value as a float to the list
    analytic.append(float(i))
    comp1.append(float(c1[count]))
    comp2.append(float(c2[count]))
    comp3.append(float(c3[count]))
    comp4.append(float(c4[count]))
    count += 1

# Certain analytic values are out of the code's resolution
a_cut = []
a_cut.append(analytic[0])
a_cut.append(analytic[1])
a_cut.append(analytic[-1])

# Cut out the points that are out of the code's resolution
comp1 = comp1[2:-1]
comp2 = comp2[1:-1]
comp3 = comp3[:-1]
comp4 = comp4[:-1]

# Create the density range
rho = [1.,1.e-1,1.e-2,1.e-3,1.e-4,1.e-5,1.e-6,1.e-7,1.e-8]
rho_comp1 = rho[2:-1]
rho_comp2 = rho[1:-1]
rho_comp34 = rho[:-1]
rho_cut = [1.,1.e-1,1.e-8]

# Plot it
plt.semilogx(rho_comp1,comp1,color='gold',marker='o',linestyle='',markersize=11,label='grid_70')
plt.semilogx(rho_comp2,comp2,color='lime',marker='o',linestyle='',markersize=9,label='grid_200')
plt.semilogx(rho_comp34,comp3,color='g',marker='o',linestyle='',markersize=7,label='grid_500')
plt.semilogx(rho_comp34,comp4,color='c',marker='o',linestyle='',markersize=5,label='grid_1000')
plt.semilogx(rho,analytic,color='blue',marker='o',linestyle='',label='Analytic')
#plt.semilogx(rho_cut,a_cut,color='indigo',marker='o',linestyle='',label='Not in range')
plt.title(r'Str$\ddot{o}$mgren Sphere Radius vs Density')
plt.xlabel(r'Density [$\mathrm{cm^{-3}}$]')
plt.ylabel(r'Radius [$\mathrm{kpc}$]')
plt.legend(loc='best')
plt.show()

plt.semilogx(rho,analytic,color='blue',marker='o',linestyle='',label='Analytic')
plt.title(r'Str$\ddot{o}$mgren Sphere Radius vs Density')
plt.xlabel(r'Density [$\mathrm{cm^{-3}}$]')
plt.ylabel(r'Radius [$\mathrm{kpc}$]')
plt.legend(loc='best')
plt.show()
