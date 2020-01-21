#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

// The tri-diagonal matrix solver, the proof that this works for diagonally dominant and non-diagonally dominant
// tri-diagonal is shown in a separate file.
void tri_solver(double *a, double *b, double *c, double *y, double *X, int n)
{
  /* Compute the array X in the equation MX = y given the lower diagonal a, diagonal
   * b, and upper diagonal c of matrix M_star
   *
   * Input
   * -----
   * a: lower diagonal of matrix M_star
   * b: diagonal of matrix M_star
   * c: upper diagonal of matrix M_star
   * y: the RHS of the equation Mx = y
   * n: the size of the square matrix (nxn)
   *
   * Output
   * ------
   * X: the solution of MX = y equation */

  int i;
  double *c_prime, *y_prime;

  // Allocate the prime arrays
  c_prime = (double *) calloc(n, sizeof(double));
  y_prime = (double *) calloc(n, sizeof(double));

  /* What the algorithm is actually doing
   * for i = 0
   * ---------
   * c_prime[i] = c[i]/b[i]
   * y_prime[i] = y[i]/b[i]
   *
   * for i = 1,2,3,...,n-2
   * ---------------------
   * c_prime[i] = c[i]/(b[i] - a[i]*c_prime[i-1])
   * y_prime[i] = (y[i] - a[i]*y_prime[i-1]/(b[i] - a[i]*c_prime[i-1])
   *
   * for i = n-1
   * ------------
   * y_prime[i] = (y[i] - a[i]*y_prime[i-1]/(b[i] - a[i]*c_prime[i-1])
   * X[i] = y_prime[i]
   *
   * for i = n-2, n-3,..., 0
   * ------------------------
   * X[i] = y_prime[i] - c_prime[i]*X[i+1] */

  for (i = 0; i < n; i++) {
    if (i == 0)
      {
	c_prime[i] = c[i] / b[i];
	y_prime[i] = y[i] / b[i];
      } else if (i > 0 && i < n - 1)
      {
	c_prime[i] = c[i] / (b[i] - a[i] * c_prime[i - 1]);
	y_prime[i] = (y[i] - a[i] * y_prime[i - 1]) / (b[i] - a[i] * c_prime[i - 1]);
      } else if (i == n - 1)
      {
	y_prime[i] = (y[i] - a[i] * y_prime[i - 1]) / (b[i] - a[i] * c_prime[i - 1]);
	X[i] = y_prime[i];
      }
  }
  for (i = n - 2; i >= 0; i--)
    {
      X[i] = y_prime[i] - c_prime[i] * X[i + 1];
    }

  // Free the memory
  free(c_prime);
  free(y_prime);

  // X has been found so exit the function
  return;
}
