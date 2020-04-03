#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include "dydx.h"

double RK4_solver(double x_curr, double f_curr, double h)
{
  /* This is a 4th order Runge-Kutta (RK) first order, constant coefficient,
   * linear ODE solver
   * RK4: u_n+1 = u_n + k1/6 + k2/3 + k3/3 + k4/6
   *
   * General form of linear ODE: y' = g(x) - p(x)*y
   * In the case of constant coefficients g(x), p(x) = g,p
   * Input
   * -----
   * g_x
   * y: the current value of y to be updated
   * p_x
   * h: The step at which x increases
   *
   * Output
   * ------
   * y_update: y_n+1
   */

  double f_update; // Updated integrand value
  double k1;       // RK k value
  double k2;       // RK k value 
  double k3;       // RK k value 
  double k4;       // RK k value 

  double x_half = x_curr + h/2.; // Needed for k calculations
  double x_p_h = x_curr + h;     // Needed for k calculations
  double f_k2;
  double f_k3;
  double f_k4;

  // Set the k-values (mentioned in the equaiton in the doc string)
  k1 = dydx(x_curr,f_curr);
  k1 = h*k1;

  f_k2 = f_curr + k1/2.;
  k2 = dydx(x_half,f_k2);
  k2 = h*k2;

  f_k3 = f_curr + k2/2.;
  k3 = dydx(x_half,f_k3);
  k3 = h*k3;

  f_k4 = f_curr + k3;
  k4 = dydx(x_p_h,f_k4);
  k4 = h*k4;

  // Solve for the new RK value
  f_update = f_curr + k1/6. + k2/3. + k3/3. + k4/6.;
  //y_update = 1. + y;
  //printf("y %2.5f\n",y_update);
  return f_update;
}
