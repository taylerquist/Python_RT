#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

void RK4_solver(double g_x, double y, double p_x, double h, double y_update)
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

  double k1;     // RK k value
  double k2;     // RK k value 
  double k3;     // RK k value 
  double k4;     // RK k value 

  // Set the k-values (mentioned in the equaiton in the doc string)
  k1 = h*(g_x - (p_x*y));
  k2 = h*((g_x - (p_x*y)) + h/2.);
  k3 = h*((g_x - (p_x*y)) + h/2.);
  k4 = h*((g_x - (p_x*y)) + h);

  // Solve for the new RK value
  y_update = y + k1/6. + k2/3. + k3/3. + k4/6.;
}
