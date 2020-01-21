#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

// Setup the source function interpolation, this includes the quadrature limiter that does not allow sources to
// become negative (aka non-physical)
void S_interp(double S_int, double *alpha_arr, double *S_arr, double *z_arr, double *up, double *pp, double *dp,
	      double *um, double *pm, double *dm, double delta_z, int x, bool plus)
{
  /* Source term interpolation function that includes a maximum quadrature limiter that does not allow any source term
   * to have a negative value.
   * This value is going to be different for the two rays going in different directions, and this is taken into
   * account here.
   *
   * Input
   * -----
   * alpha_arr: Array containing the value of the emissivity at each grid point, double precision, size n
   * S_arr: An array containing the values of the source function at every point in the grid, double precision, size n
   * z_arr: The grid, double precision, size n
   * delta_z: The step size between grid points, double precision, single value
   * x: The location within the grid that needs to be interpolated, x in z[x], integer, single value.
   * plus: This indicates which array is being dealt with; I+ --> plus = True or I- --> plus = False
   *
   * Output
   * ------
   * S_int: The value of the interpolated source function for either the I+ or I- array at grid location x,
   *        from z[x], double precision, single value */

  double Q_max;  // This is the upper value of the limiter
  double Q_curr; // This is the value of the interpolation using the coefficients ane S_arr

  /* The quadrature limiter for third order integration is discussed in Ch.3, section 3.8.4, equation 3.39
   * This equation is for the I+ vector: Q_curr = up*S[x-1] + pp*S[x] + dp*S[x+1]
   * For the maximum upper limit the equation is found in Ch.3, section 3.8.4, equation 3.37
   * Q_max = 0.5*(alpha[x] + alpha[x-1])*delta_z
   * The final value of the source term interpolation is determined by equation 3.46, section 3.8.4, Ch.3
   * S_int = max(min(Q_curr[x],Q_max[x]),0.) */

  Q_max = 0.5*(alpha_arr[x] - alpha_arr[x-1])*delta_z;

  if (plus==true)
    {
      Q_curr = up[x]*S_arr[x-1] + pp[x]*S_arr[x] + dp[x]*S_arr[x+1];
      S_int = fmax(fmin(Q_curr,Q_max),0.);
    }
  else
    {
      Q_curr = um[x]*S_arr[x+1] + pm[x]*S_arr[x] + dm[x]*S_arr[x-1];
      S_int = fmax(fmin(Q_curr,Q_max),0.);
    }

  // S_int has been found for either I+ or I- so return the function
  return;

}
