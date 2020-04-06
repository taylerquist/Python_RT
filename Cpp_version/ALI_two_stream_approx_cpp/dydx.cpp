#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

double dydx(double x, double f, double rho)
{
  double c_HII = 3.;      // Clumping fraction for the IGM at z=3, T=20,000K 
  double alpha_B;         // Recombination coefficient 
  double Y_p = 0.24;      // Fraction of Helium 
  double X_p = 1. - Y_p;  // Fraction of Hydrogen   
  double n_h_avg;         // Hydrogen number density in each cell 
  double n_dot_ion;       // Ionizing photon rate
  double t_rec_inv;       // Recombination time

  // Set the constants needed for (t_rec)**(-1)
  alpha_B = 2.6*pow(10.,-13.);
  n_h_avg = rho*X_p;
  n_dot_ion = 0.5;
  t_rec_inv = c_HII*alpha_B*(1. + (Y_p/(4.*X_p)))*n_h_avg*64.;

  // Return the spatially dependent solution
  return ((n_dot_ion/n_h_avg) - (f*t_rec_inv));
}
