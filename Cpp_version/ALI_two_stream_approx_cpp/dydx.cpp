#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

double dydx(double x, double f, double rho, double avg_I)
{

  /* 
  * Units of each input/output
  * f: This is the current value of QHII, dimensionless 
  * rho: #/cm^3
  * avg_I: (erg)(cm^-2)(s^-1)(Hz^-1)(rad^-2)
  */

  double c_HII = 3.;      // Clumping fraction for the IGM at z=3, T=20,000K 
  double alpha_B;         // Recombination coefficient, cm**3/s 
  double sigma_nu;        // Photoionization cross-section of H ions
  double Y_p = 0.24;      // Fraction of Helium 
  double X_p = 1. - Y_p;  // Fraction of Hydrogen   
  double n_h_avg;         // Hydrogen number density in each cell 
  double n_dot_ion;       // Ionizing photon rate
  double t_rec_inv;       // Recombination time

  double pi = 3.14159265358979323846; // Value of pi 
  double c_light = 3.*pow(10.,10.);   // Speed of light in cm/s
  double eps_ion = pow(10.,25.25);    // Ionizing photon production efficiency Hz/erg
  double lmfp;

  double dQp, dQn;
  double answer;

  // Set the constants needed for (t_rec)**(-1)
  alpha_B = 2.6*pow(10.,-13.);
  alpha_B = 2.6e-13;
  sigma_nu = 7.91e-18;
  //sigma_nu = 7.91*pow(10.,-18.);
  n_h_avg = rho*X_p;
  //  * avg_I: (erg)(cm^-2)(s^-1)(Hz^-1)(rad^-2)
  //n_dot_ion = avg_I*(4*pi)*c_HII*n_h_avg*sigma_nu*eps_ion;
  //  * n_dot_ion: (cm^-3)(s^-1)
  /*if(f<1)
  {
    lmfp = 1./(c_HII*(1-f)*n_h_avg*sigma_nu);
  }else{
    lmfp = 1.0e24;
  }*/    
  lmfp = 1./(c_HII*(1-f)*n_h_avg*sigma_nu);
  if(lmfp>1.0e24)
    lmfp = 1.0e24;
  //lmfp = 
  lmfp = 2.782235e+22;
  n_dot_ion = avg_I*(4*pi)*eps_ion/lmfp;


  //n_dot_ion  = (4*pi/c)*avg_I*eps_ion; //# density of ionizing photons

  t_rec_inv = c_HII*alpha_B*(1. + (Y_p/(4.*X_p)))*n_h_avg;


  // Return the spatially dependent solution
  dQp = n_dot_ion/n_h_avg;
  dQn = f*t_rec_inv;
  answer =  dQp - dQn;
  //printf("t %e avg_I %e n_dot_ion %e lmfp %e n_h_ave %e f %e t_rec_inv %e dQ %e %e answer %e\n",x,avg_I,n_dot_ion,lmfp,n_h_avg,f,t_rec_inv,dQp,dQn,answer);

  return answer;
}
