#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

// Define the third order quadratic interpolation coefficients function
void quad_int(double *up, double *pp, double *dp, double *um, double *pm, double *dm, double *delta_tau, int N)
{
  /* Computes the three quadratic interpolation coefficients for the bottom up and top down rays
   *
   * Input
   * ------
   * delta_tau: the difference in optical depth for each grid cell, double precision, size n+1=N
   * N: the size of the grid+1: n+1=N, integer, single value
   *
   * Output
   * ------
   * up, pp, dp: the three interpolation coefficients for the ray integrating from the bottom of the grid up, size n
   * um, pm, dm: the three interpolation coefficients for the ray integrating from the top of the grid down, size n */

  double *e0_p, *e1_p, *e2_p, *e0_m, *e1_m, *e2_m; // interpolation coefficients
  int k; // Cell index

  e0_p = (double *) calloc(N, sizeof(double)); // Array used to calculate interpolation coefficients
  e1_p = (double *) calloc(N, sizeof(double)); // Array used to calculate interpolation coefficients
  e2_p = (double *) calloc(N, sizeof(double)); // Array used to calculate interpolation coefficients

  e0_m = (double *) calloc(N, sizeof(double)); // Array used to calculate interpolation coefficients                                                                                                                
  e1_m = (double *) calloc(N, sizeof(double)); // Array used to calculate interpolation coefficients                                                                                                                
  e2_m = (double *) calloc(N, sizeof(double)); // Array used to calculate interpolation coefficients 

  // Loop through and set up the optical depth dependence constants for each height and ray
  // Then use these constants to set the actual interpolation coefficients for each ray
  // Optical depth dependence constants: Ch.3 section 3.8.4, equations 3.43-3.45
  // Interpolation coefficients: Ch.3, section 3.8.4, equations 3.40-3.42 
  // Note that my notation has the following equivalences: u=u, v=p, w=d 

  // This loop sets the bottom up ray (I+) interpolation coefficients
  for (k=0;k<N-1;k++)
    {
      if((delta_tau[k]>0)&&(delta_tau[k+1]>0))
      {
        e0_p[k] = 1. - exp(-1*delta_tau[k]);
        e1_p[k] = delta_tau[k] - e0_p[k];
        e2_p[k] = pow(delta_tau[k],2.) - (2.*e1_p[k]);

        dp[k] = e0_p[k] + (e2_p[k] - (2.*delta_tau[k] + delta_tau[k+1])*e1_p[k])/(delta_tau[k]*(delta_tau[k]+
                                                                                              delta_tau[k+1]));
        pp[k] = (((delta_tau[k]+delta_tau[k+1])*e1_p[k])-e2_p[k])/(delta_tau[k]*delta_tau[k+1]);
        up[k] = (e2_p[k] - delta_tau[k]*e1_p[k])/(delta_tau[k+1]*(delta_tau[k]+delta_tau[k+1]));
      }else{
        dp[k] = 0;
        pp[k] = 0;
        up[k] = 0;
      }
    }

  // This loop sets the top down (I-) interpolation coefficients
  //for (k=0;k<N;k++) //why?
  for (k=0;k<N-1;k++)
    {
      if((delta_tau[k]>0)&&(delta_tau[k+1]>0))
      {
        e0_m[k] = 1. - exp(-1*delta_tau[k+1]);
        e1_m[k] = delta_tau[k+1] - e0_m[k];
        e2_m[k] = pow(delta_tau[k+1],2.) - (2.*e1_m[k]);

        dm[k] = e0_m[k] + (e2_m[k] - (2.*delta_tau[k+1] + delta_tau[k])*e1_m[k])/(delta_tau[k+1]*(delta_tau[k]+
											    delta_tau[k+1]));
        pm[k] = (((delta_tau[k+1]+delta_tau[k])*e1_m[k])-e2_m[k])/(delta_tau[k+1]*delta_tau[k]);
        um[k] = (e2_m[k] - delta_tau[k+1]*e1_m[k])/(delta_tau[k]*(delta_tau[k+1]+delta_tau[k]));
      }else{
        dm[k] = 0;
        pm[k] = 0;
        um[k] = 0;        
      }
    }

  // Free the memory
  free(e0_p);
  free(e1_p);
  free(e2_p);
  free(e0_m);
  free(e1_m);
  free(e2_m);
}
