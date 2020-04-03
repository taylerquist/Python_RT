#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <fstream>
#include "runge_kutta_solver.h"

using std::ofstream;
using namespace std;

/* This is an accuracy test of the RK4 solver scheme for 
 * a liner ODE (which is of the same form as that seen in
 * the cosmological experiment).
 * Linear ODE general form: y' + p(x)y = g(x)
 * General form of the solution: mu*y = int(mu*g(x)*dx)
 * Where mu is the integration factor: mu = e**(int(p(x)*dx))
 * -------------------------------------------------------
 * The test case ODE: y' = 2 - (1/4)*y
 * General solution: y*e**(x/4) = 8*e**(x/4)
 * The analytical solution will be compared with the RK solver solution
 */

int main()
{
  // Allocate the variables
  float t;           // Counters
  double p = 0.25; // Value of p(x)
  double g = 2.;   // Value of g(x)
  double y;        // The variable we are evolving
  double y_new;    // Update the value of y
  double y_act;    // Analytic solution of y
  double mu;       // Integration factor that will change with each round
  double h = 0.1;  // Step size

  ofstream y_file;     // y output file
  ofstream y_act_file; // Analytic y output file
  ofstream t_file;     // t output file

  // Initialize
  y = 0.;
  y_act = 0.;
  y_new = 0.;
  mu = 0.;

  for (t=0;t<20;t=t+h)
    {
      // Analytic solution
      mu = exp(t*p);
      y_act = (8.*exp(t*p))/mu;
      // Runge-Kutta 4 solution
      RK4_solver(g, y, p, h, y_new);
      y = y_new;
      printf("y %2.5f\n",y);
    }

  return 0;
}
