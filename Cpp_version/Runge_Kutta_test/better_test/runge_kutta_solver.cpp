#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include "dydx.h"

void RK4_solver(double x_curr, double *f_curr, double h, double *dens, int n, double *j, double *f_update)
{
  /* This is a 4th order Runge-Kutta (RK) first order, constant coefficient,
   * linear ODE solver
   * RK4: u_n+1 = u_n + k1/6 + k2/3 + k3/3 + k4/6
   *
   * General form of linear ODE: y' = g(x) - p(x)*y
   * In the case of constant coefficients g(x), p(x) = g,p
   * Input
   * -----
   * x_curr: Current value of x
   * f_curr: Current value of the function f(x) across the grid
   * h: The step at which x increases
   * dens: The array of total gas density for each cell
   * n: Number of cells
   * j: The average specific intensity
   *
   * Output
   * ------
   * f_update: The updated f array
   */

  double k1; // RK4 adjusting factors
  double k2;
  double k3;
  double k4;
  double f_k2; // Extra doubles for keeping things clean
  double f_k3;
  double f_k4;

  double x_half = x_curr + h/2.; // Needed for k calculations
  double x_p_h = x_curr + h;     // Needed for k calculations

  int i;

  // Solve for the new values of f
  for (i=0;i<n;i++)
    {
      k1 = dydx(x_curr,f_curr[i],dens[i],j[i]);
      k1 = h*k1;

      f_k2 = f_curr[i] + k1/2.;
      k2 = dydx(x_half,f_k2,dens[i],j[i]);
      k2 = h*k2;

      f_k3 = f_curr[i] + k2/2.;
      k3 = dydx(x_half,f_k3,dens[i],j[i]);
      k3 = h*k3;

      f_k4 = f_curr[i] + k3;
      k4 = dydx(x_p_h,f_k4,dens[i],j[i]);
      k4 = h*k4;

      // Solve for the new RK value
      f_update[i] = f_curr[i] + k1/6. + k2/3. + k3/3. + k4/6.;
    }

  // f_update should be filled in
}
