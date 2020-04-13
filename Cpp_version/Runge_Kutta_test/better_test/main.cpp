#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <fstream>
#include "true_sln_slvr.h"
#include "runge_kutta_solver.h"
#include "dydx.h"

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
  double x;        // Counters
  double h = 0.001;  // Step size

  int n=5; // Size of grid
  int k;

  // For the test case these are all random coefficients
  // but they are named similarly to variables in the real problem
  double *J = new double[n];
  double *density = new double[n]; 
  double *f_new = new double[n];
  double *f_act = new double[n];

  ofstream true_sln; // True solution output file
  ofstream RK_sln;   // Solver output file

  // Initialize
  for (k=0;k<n;k++)
    {
      J[k] =  double(k);
      density[k] = 1./(1. + double(k));
      f_new[k] = 0.;
      f_act[k] = 0.;
    }

  true_sln.open("true_file.txt");
  RK_sln.open("RK_file.txt");

  for (x=0.1;x<n;x=x+h)
    {
      true_sln_slvr(x, density, n, J, f_act);
      RK4_solver(x, f_new, h, density, n, J, f_new);
      for (k=0;k<n;k++)
	{
	  true_sln << f_act[k];
	  true_sln << "\n";
	  RK_sln << f_new[k];
	  RK_sln << "\n";
	} 
      true_sln << "break\n";
      RK_sln << "break\n";
      
    }

  true_sln.close();
  RK_sln.close();

  // So far we are only re-writing the f-grid
  // Will want to eventually record each grid value for plotting

  for (k=0;k<n;k++)
    {
      printf("f_act %2.5f\n",f_act[k]);
      printf("y %2.5f\n",f_new[k]);
      printf("------------\n"); 
    }

  return 0;
}
