#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

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

  double *e0, *e1, *e2; // interpolation coefficients
  int k;                // Cell index

  e0 = (double *) calloc(N, sizeof(double)); // Array used to calculate interpolation coefficients
  e1 = (double *) calloc(N, sizeof(double)); // Array used to calculate interpolation coefficients
  e2 = (double *) calloc(N, sizeof(double)); // Array used to calculate interpolation coefficients

  // Loop through and set up the optical depth dependence for each height
  // Ch.3 section 3.8.4, equations 3.43-3.45
  for (k=0;k<N;k++)
    {
      e0[k] = 1. - exp(-delta_tau[k]);
      e1[k] = delta_tau[k] - e0[k];
      e2[k] = pow(delta_tau[k],2.) - (2.*e1[k]);
    }

  // Now compute the actual interpolation coefficients for each ray using the constants above
  // The equations for these can be found in Ch.3, section 3.8.4, equations 3.40-3.42
  // Note that my notation has the following equivalences: u=u, v=p, w=d
  for (k=0;k<N;k++)
    {
      // The coefficients for the bottom up array (I+)
      up[k] = e0[k] + (e2[k] - (2.*delta_tau[k] + delta_tau[k+1])*e1[k])/(delta_tau[k]*(delta_tau[k]+delta_tau[k+1]));
      pp[k] = (((delta_tau[k]+delta_tau[k+1])*e1[k])-e2[k])/(delta_tau[k]+delta_tau[k+1]);
      dp[k] = (e2[k] - delta_tau[k]*e1[k])/(delta_tau[k+1]*(delta_tau[k]+delta_tau[k+1]));

      // The coefficients for the top down array (I-)
      um[k] = e0[k] + (e2[k] - (2.*delta_tau[k] + delta_tau[k-1])*e1[k])/(delta_tau[k]*(delta_tau[k]+delta_tau[k-1]));
      pm[k] = (((delta_tau[k]+delta_tau[k-1])*e1[k])-e2[k])/(delta_tau[k]+delta_tau[k-1]);
      dm[k] = (e2[k] - delta_tau[k]*e1[k])/(delta_tau[k-1]*(delta_tau[k]+delta_tau[k-1]));
    }

  // Free the memory
  free(e0);
  free(e1);
  free(e2);
}
