import numpy as np
import h5py
import matplotlib.pyplot as plt

hf  = h5py.File('los_130.h5','r')

# See what groups are in the file
#print(hf.keys())

# Choose a group to set the density array
rho = hf.get('density_x')
rho = np.array(rho)
l = len(rho)
#print(len(rho))
x = np.linspace(0.,l,l)

# (Manually) centered around the largest density jump)
# Plot it -- this is the section that will be taken for the density profile
#plt.plot(x[1899:1969],rho[1899:1969])
#plt.show()

#print(np.amax(rho))

# Trim down to the desired length as shown in the plot above
rho = rho[1899:1969]

#print(np.amax(rho))

# Interpolate to desired grid density
interp = True
n = 200
x_interp = np.linspace(0,len(rho)-1,n)
x_norm = np.linspace(0,len(rho)-1,70)
rho_interp = np.interp(x_interp,x_norm,rho)

# Convert to kpc for plotting
reg_con = 622.08123/69.
#interp_con = 622.08123/199.

x_interp = x_interp*reg_con
x_norm = x_norm*reg_con

#print(x_interp)
#print(x_norm)

plt.semilogy(x_interp,rho_interp, color='c',label='Interp')
plt.semilogy(x_norm,rho, color='b', linestyle='--',label='Non-interp')
plt.title('Interpolated density: grid size 200')
plt.xlabel(r'Radius [$\mathrm{kpc}$]')
plt.ylabel(r'Density [$\mathrm{h^{2} M_{sun} kpc^{3}}$]')
plt.legend()
plt.show()

# Save rho as a .txt file to read into C++ code
if (interp == True):
    new_file = open('density_profile.txt','w')
    for i in range(0,len(rho_interp)):
        string = str(rho_interp[i]) + ' \n'
        new_file.writelines(string)
else:
    new_file = open('density_profile.txt','w')
    for i in range(0,len(rho)):
        string = str(rho[i]) + ' \n'
        new_file.writelines(string)

# Close the file
new_file.close()

# Plot converted interpolation
kpc_cm = 3.086e21
M_sun_g = 1.989e3
M_H_g = 1.6726e-24
redshift = 3.
lil_h = 0.6766

density_conversion = (M_sun_g/(kpc_cm)**3.)*(1./M_H_g)*((redshift+1.)**3.)*(lil_h**2.)
#delta_x = 622.08123/(n-1);

rho_interp = rho_interp*density_conversion
rho = rho*density_conversion
plt.semilogy(x_interp,rho_interp, color='violet',label='Interp')
plt.semilogy(x_norm,rho, color='indigo', linestyle='--',label='Non-interp')
plt.title('Interpolated density: grid size 1000')
plt.xlabel(r'Radius [$\mathrm{kpc}$]')
plt.ylabel(r'Density [$\mathrm{cm^{-3}}$]')
plt.legend()
#plt.show()
