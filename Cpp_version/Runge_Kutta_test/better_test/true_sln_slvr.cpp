#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

void true_sln_slvr(double x_curr, double *dens, int n, double *j, double *f_update)
{
  /* This solves for the true integration value of the following quadratic
   * integrate(dens[i]*(x**2.) + j[i]*(x)) = 
   * (dens[i]/3.)*(x**3.) + (j[i]/2.)*(x**2.)
   */
  
  int z;
  
  for (z=0;z<n;z++)
    {
      f_update[z] = (dens[z]/3.)*pow(x_curr,3.) + (j[z]/2.)*pow(x_curr,2.);
    }

}
