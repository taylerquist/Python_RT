#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

void RK4_solver(ph_rate,Q_n,t_rec,h,Q_new)
{
  /* This is a 4th order Runge-Kutta (RK) first order, constant coefficient,
   * linear ODE solver
   * RK4: u_n+1 = u_n + k1/6 + k2/3 + k3/3 + k4/6
   *
   * Input
   * -----
   * ph_rate: The rate at which photons are produced --> defined by dot(n_ion)/<n_h>
   * Q_n: This is the current value of the ionized fraction (we are trying to solve Q_n+1
   * t_rec: The IGM recombination time
   * h: The time step
   * N_step: Number of time steps (T_final/h)
   *
   * Output
   * ------
   * Q_new = Q_n+1, the new ionized fraction value
   */

  double k1;     // RK k value
  double k2;     // RK k value 
  double k3;     // RK k value 
  double k4;     // RK k value 

  // Set the k-values (mentioned in the equaiton in the doc string)
  k1 = h*(ph_rate - (Q_n/t_rec));
  k2 = h*((ph_rate - (Q_n/t_rec)) + h/2.);
  k3 = h*((ph_rate - (Q_n/t_rec)) + h/2.);
  k4 = h*((ph_rate - (Q_n/t_rec)) + h);

  // Solve for the new RK value
  Q_new = Q_n + k1/6. + k2/3. + k3/3. + k4/6.
}
