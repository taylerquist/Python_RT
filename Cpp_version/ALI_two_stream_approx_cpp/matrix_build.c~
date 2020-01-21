#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
void matrix_build(double **lambda_s, double *m_a, double *m_b, double *m_c, double *up, double *pp, double *dp,
		  double *um, double *pm, double *dm, double eps, int n)
{
  /* This will compute the lambda_star matrix and M_star lower and upper diagonals and the diagonal.
   *
   * Input
   * -----
   * up, pp, dp: the three interpolation coefficients for the ray integrating from the bottom of the grid up, size n
   * um, pm, dm: the three interpolation coefficients for the ray integrating from the top of the grid down, size n
   * eps: Aka epsilon, the photon destruction probability, double precision, single value
   * n: size of the grid, integer, single value
   *
   * Output
   * ------
   * lambda_s: lambda_star matrix, double precision, size nxn
   * m_a = lower diagonal of the M_star matrix, double precision, size n
   * m_b = diagonal of the M_star matrix, double precision, size n
   * m_c = upper diagonal of the M_star matrix, double precision, size n */

  int k;                                        // Cell index
  double lambda_a[n], lambda_b[n], lambda_c[n]; // Temporary lower diagonal, diagonal, and upper diagonal to
  // accurately keep track of each one

  /* Loop through and set each of each tri-diagonal components for both lambda_star and M_star
   * Equations for each lambda_star tri-diagonal component can be found in Ch.4, section 4.4.5, equations 4.57-4.59
   * Note that my notation has the following equivalences: u=u, v=p, w=d
   * lambda_a = 0.5(up + dm), lambda_b = 0.5(pp + pm), lambda_c = 0.5(dp + um)
   * The equation for M_star is based off of lambda_star and can be found in Ch.4, section 4.4.4, equation 4.47
   * M_star = Identity - lambda_star * (1 - epsilon) */

  for (k=0;k<n;k++)
    {
      if (k==0)
        {
	  // Lambda_star upper, lower, and direct diagonals
	  lambda_a[k] = 0.;
	  lambda_b[k] = 0.5*(pp[k]+pm[k]);
	  lambda_c[k] = 0.5*(dp[k]+um[k]);
	  // M_star upper, lower, and direct diagonals
	  m_a[k] = 0.;
	  m_b[k] = 1. - (1. - eps)*lambda_b[k];
	  m_c[k] = -(1. - eps)*lambda_c[k];
	  // Lambda_star matrix
	  lambda_s[k][k] = lambda_b[k];
	  lambda_s[k+1][k] = lambda_c[k];
        }
      else if (k==n-1)
        {
	  // Lambda_star upper, lower, and direct diagonals
	  lambda_a[k] = 0.5*(up[k]+dm[k]);
	  lambda_b[k] = 0.5*(pp[k]+pm[k]);
	  lambda_c[k] = 0.;
	  // M_star upper, lower, and direct diagonals
	  m_a[k] = -(1. - eps)*lambda_a[k];
	  m_b[k] = 1. - (1. - eps)*lambda_b[k];
	  m_c[k] = 0.;
	  // Lambda_star matrix
	  lambda_s[k-1][k] = lambda_a[k];
	  lambda_s[k][k] = lambda_b[k];
        }
      else
        {
	  // Lambda_star upper, lower, and direct diagonals
	  lambda_a[k] = 0.5*(up[k]+dm[k]);
	  lambda_b[k] = 0.5*(pp[k]+pm[k]);
	  lambda_c[k] = 0.5*(dp[k]+um[k]);
	  // M_star upper, lower, and direct diagonals
	  m_a[k] = -(1. - eps)*lambda_a[k];
	  m_b[k] = 1. - (1. - eps)*lambda_b[k];
	  m_c[k] = -(1. - eps)*lambda_c[k];
	  // Lambda_star matrix
	  lambda_s[k-1][k] = lambda_a[k];
	  lambda_s[k][k] = lambda_b[k];
	  lambda_s[k+1][k] = lambda_c[k];
        }
    }
  // All components of lambda_star and M_star have been found
  return;
}
