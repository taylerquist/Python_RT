#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
void tau_fn(double *delta_tau, const double *z_arr, const double delta_z, const int N)
{
  /* Computes the difference in optical depth for each cell in the grid
   *
   * Input
   * ------
   * z_arr: the spatial grid over which the optical depth will change, double precision, size N=n+1
   * delta_z: the step size of the grid, double precision, single value
   * N: the size of the grid + 1: n+1=N, integer, single value
   *
   * Output
   * ------
   * delta_tau: an array containing the difference in optical depth at each grid for each cell, double precision,
   *            size N */

  int k;    // Cell index
  double x; // Temporary double precision variable for each loop

  for (k=0;k<N;k++)
    {
      x = z_arr[k] + (delta_z/2.);
      // The definition of the difference in optical depth for this problem is: sqrt(3.)*alpha[i]*delta_z
      // This can be found in Ch.4, section 4.4.5, equation 4.55
      // The extinction coefficient (alpha[i]) can be found in Ch.4: section 4.4.2, equation 4.32
      // alpha[i] = (10.^(5. - 6.*x[i])), where x[i] is our current cell (not grid point)
      delta_tau[k] = sqrt(3.)*pow(10.,(5.-(6.*x)))*delta_z;
    }

  // delta_tau has been found!
}
