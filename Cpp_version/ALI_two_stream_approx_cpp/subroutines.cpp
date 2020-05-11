#include <math.h>
#include <stdio.h>
#include "subroutines.hpp"
#include "quad_int.h"
#include "s_interp.h"
#include "runge_kutta_solver.h"


RT::RT(void)
{
  SetConstants();
  InitializeGrid();
  printf("Done initializing grid\n");
}

void RT::PrintGrid(FILE *fp)
{
  printf("i,x,I+,I-,S,QHII,Dtau\n");
  for(int i=0;i<n;i++)
    fprintf(fp,"%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",i,xg[i]/kpc_cm,density[i],I_plus[i],I_minus[i],J[i],S[i],QHII[i],del_tau[i]);
}
void RT::SaveGrid(FILE *fp)
{
  for(int i=0;i<n;i++)
    fprintf(fp,"%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",i,xg[i]/kpc_cm,density[i],I_plus[i],I_minus[i],J[i],S[i],QHII[i],del_tau[i]);
}
void RT::PrintIntensity(FILE *fp)
{
  printf("Intensity (i,x,I+,I-,J,B,S\n");
  for(int i=0;i<n;i++)
    fprintf(fp,"%d\t%e\t%e\t%e\t%e\t%e\t%e\n",i,xg[i]/kpc_cm,I_plus[i],I_minus[i],J[i],B[i],S[i]);
}
void RT::PrintIonization(FILE *fp)
{
  printf("Ionization (i,x,QHII,epsilon,alpha,del_tau\n");
  for(int i=0;i<n;i++)
    fprintf(fp,"%d\t%e\t%e\t%e\t%e\t%e\n",i,xg[i]/kpc_cm,QHII[i],eps_grid[i],alpha[i],del_tau[i]);
}
void RT::SetConstants(void)
{
  n = 70; //number of grid cells
  sigma_nu = 7.91e-18;  // Photoionization of H1
  pi = 3.14159265358979323846;    // Value of pi
  //delta_x;                        // Physical size of grid cell at z=3
  //Delta;                          // Volumetric ratio to scale cell sizes
  kpc_cm = 3.086e21;    // Cm in one kpc
  M_sun_g = 1.989e33;   // Mass of the sun in g
  M_H_g = 1.6726e-24;   // Mass of H1 in grams
  redshift = 3.;                  // Redshift of the density profile
  lil_h = 0.6766;                 // Scaled Hubble constant
  // Unit conversion constant for the density distribution
  density_conversion = (M_sun_g/(kpc_cm*kpc_cm*kpc_cm))*(1./M_H_g)*((redshift+1.)*(redshift+1.)*(redshift+1.))*
    (lil_h*lil_h);
  delta_x = (24.4/lil_h)*(1./(redshift+1.))*(kpc_cm); // Size of one grid cell in cm
  //printf("delta_x = %e\n",delta_x);
  del_z = 1./float(n); // Grid step size

  //density_limit = 1.0e-2;
  density_limit = 100;

  sigma_T = 0.66524587158e-24; // Thomson cross section

  I_p_init = 1.0e-17; // The scaled bottom-up specific intensity; 1 --> 10^-17*Delta
  //I_p_init = 1.0e-15; // The scaled bottom-up specific intensity; 1 --> 10^-17*Delta
  //I_m_init = 1.0e-22; // The scaled top-down specific intensity; 10^-22*Delta/10^-17*Delta
  //I_p_init = 1.0e-17;
  I_m_init = 0;

}

void RT::InitializeGrid(void)
{

  FILE *fp;
  double s;

  // Set other grids, interpolation coefficients, and matrices
  //----------------------------------------------------------------

  // Setting the optical depth array
  // Set a density grid and other necessary values 
  density = new double[n]; // Density profile in cm^-3
  col_dens = new double[n];// Column density for HI in cm^-2
  I_plus = new double[n];  // The specific intensity for the vector going from bottom up
  I_minus = new double[n]; // The specific intensity for the vector going from top down
  z = new double[n];       // The grid array
  B = new double[n];       // The thermal source function (aka the Planck Function)
  S = new double[n];       // The thermal source function (aka the Planck Function)
  J = new double[n];       // The thermal source function (aka the Planck Function)

  del_tau = new double[n]; // Change in optical depth
  u_p = new double[n];     // "Plus" u interpolation coefficient array (the vector going from bottom up)
  p_p = new double[n];     // "Plus" p interpolation coefficient array (the vector going from bottom up)
  d_p = new double[n];     // "Plus" d interpolation coefficient array (the vector going from bottom up)
  u_m = new double[n];     // "Minus" u interpolation coefficient array (the vector going from top down)
  p_m = new double[n];     // "Minus" p interpolation coefficient array (the vector going from top down)
  d_m = new double[n];     // "Minus" d interpolation coefficient array (the vector going from top down)
  alpha_p = new double[n]; // "Plus" alpha interpolation coefficient array (the vector going from bottom up)
  beta_p = new double[n];  // "Plus" beta interpolation coefficient array (the vector going from bottom up)
  gamma_p = new double[n]; // "Plus" gamma interpolation coefficient array (the vector going from bottom up)
  alpha_m = new double[n]; // "Minus" alpha interpolation coefficient array (the vector going from top down)
  beta_m = new double[n];  // "Minus" beta interpolation coefficient array (the vector going from top down)
  gamma_m = new double[n]; // "Minus" gamma interpolation coefficient array (the vector going from top down)
  M_a = new double[n];     // Lower diagonal of the matrix M_star
  M_b = new double[n];     // Diagonal of the matrix M_star
  M_c = new double[n];     // Upper diagonal of the matrix M_star
  alpha = new double[n];   // Emissivity
  y_arr = new double[n];   // y array in the equation Mx = y
  I_plus = new double[n];  // The specific intensity for the vector going from bottom up
  z = new double[n];       // The grid array
  eps_grid = new double[n];
  xg = new double[n]; //grid locations

  QHII = new double[n];

  del_tau_int_init = new double[n]; // Initial I+ and I- integration values
  epsilon = 1.0;


  // Read in the density profile .txt file
  fp = fopen("density_profile.txt","r");
  for(int i=0;i<n;i++)
  {
    fscanf(fp,"%lf\n",&s);
    density[i] = s;
  }
  fclose(fp);

  double j;

  for(int i=0;i<n;i++)
  {
    j = (double) (i);
    z[i] = j*del_z;                             // Setting the grid
    xg[i] = j*delta_x;
    B[i] = 0.;                                  // Thermal source function
    //density[i] = density[i]*density_conversion; // Scaled density field 
    density[i] = 1.0e-8;

    //density[i] = 1.0e-4;

    //if(density[i]>density_limit)
    //density[i] = density_limit;

    col_dens[i] = density[i]*delta_x;           // Scaled column density
    del_tau[i] = col_dens[i]*sigma_nu;          // Change in optical depth

    //del_tau[i] = 0.;
    j = float(i) + 1.;
    Delta = pow(((j-1)/j),2.);

    //printf("j %e Delta %e\n",j,Delta);

    if(i==0)
    {
      I_plus[i] = I_p_init;                   // Bottom-up specific intensity
    }else{
      I_plus[i] = I_plus[i-1]*Delta;          // Bottom-up specific intensity
    }

    if(i==n-1)
      I_minus[i] = I_m_init;                      // Top-down specific intensity

    //printf("i%d I_plus %e I_minus %e\n",i,I_plus[i],I_minus[i]);

    //set mean intensity
    J[i] = 0.5*(I_plus[i]+I_minus[i]);          // Mean specific intensity

    //Set QHII to 0 initially
    QHII[i] = 0;

    //compute the extinction coefficient
    alpha[i] = (1.0-QHII[i])*sigma_nu*density[i];

    //compute the epsilon on the grid
    eps_grid[i] = (sigma_nu*(1. - QHII[i]))/(sigma_nu*(1. - QHII[i]) + (sigma_T*QHII[i]));

    // Source function
    S[i] = eps_grid[i] * B[i] + (1. - eps_grid[i]) * J[i];

  }

  //Set the array of third order quadratic interpolation coefficients
  quad_int(u_p,p_p,d_p,u_m,p_m,d_m,del_tau,n);
}
void RT::UpdateIntensity(void)
{
  int i;
  double j;
  double Delta;
  bool Plus;

  double S_val_plus, S_val_minus;
  //loop over the grid and update the intensity

  //first, compute I+
  I_plus[0] = I_p_init;
  for(i=1;i<n;i++)
  {
    //compute 1/r^2 attenuation
    j = float(i) + 1.;
    //Delta = pow(((j-1)/j),2.);
    //Delta = pow(0.5/(float(i)+0.5),2);
    Delta = pow(((j-1 + 0.5)/(j+0.5)),2.);


    // First solve the bottom-up interpolated value of S
    Plus = true;
    S_val_plus = S_interp(alpha,S,u_p,p_p,d_p,u_m,p_m,d_m,delta_x,i,Plus);
    
    //compute I+, with extra 1/r^2 attenuation
    I_plus[i] = exp(-del_tau[i-1]) * Delta * I_plus[i-1] + S_val_plus;
  }

  //second, compute I-
  I_minus[n-1] = I_m_init;
  for(i=n-2;i>=0;i--)
  {
    // First solve the top-down interpolated value of S
    Plus = false;
    S_val_minus = S_interp(alpha,S,u_p,p_p,d_p,u_m,p_m,d_m,delta_x,i,Plus);

    //compute I-
    I_minus[i] = exp(-del_tau[i+1]) * I_minus[i+1] + S_val_minus;
  }

  //compute J and update S
  for(i=0;i<n;i++)
  {

    //mean intensity
    J[i] = 0.5*(I_plus[i] + I_minus[i]);

    // Source function, B[i] constant, eps_grid computed in ionization routine
    S[i] = eps_grid[i] * B[i] + (1. - eps_grid[i]) * J[i];
  }

}

void RT::UpdateIonization(void)
{
  int i;

  int iq, nq;

  double h = 3.154e7 * 1e6; //RK4 step in seconds (1 M year)

  double t = 0;

  //how many QHII iterations?
  nq = 1000;

  for(iq=0;iq<nq;iq++)
  {
    //printf("iq %d h %e\n",iq,h);
    // Solve the ODE to update QHII
    RK4_solver(t, QHII, h, density, n, J, QHII);

    t += h;
  }

  //update alpha and del_tau
  for(i=0;i<n;i++)
  {
    //compute the extinction coefficient
    alpha[i] = (1.0-QHII[i])*sigma_nu*density[i];

    //compute the change in optical depth
    del_tau[i] = alpha[i]*delta_x;

    //compute the epsilon
    eps_grid[i] = (sigma_nu*(1. - QHII[i]))/(sigma_nu*(1. - QHII[i]) + (sigma_T*QHII[i]));

  }
  //Set the array of third order quadratic interpolation coefficients
  quad_int(u_p,p_p,d_p,u_m,p_m,d_m,del_tau,n);
}
