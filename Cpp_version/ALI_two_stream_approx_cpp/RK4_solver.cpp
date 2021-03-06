#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

void mat_mult(double **mat_a, double *vec_b, double *mat_c, int N)
{
  /* This is a matrix multiplication function that can multiply a matrix with a column vector that are compatible;
   * For example: mxn * nx1 is compatible but nxm * nx1 is not compatible.
   * mat_a * vec_b = mat_c
   *
   * Input
   * -----
   * N: the size of the grid (n)
   * mat_a: this is the first matrix on the LHS of the equation above
   * vec_b: this is the column vector on the LHS of the equation above, in C++ it must be stored as a matrix
   * N: the size of the matrix (NxN) and column vector (Nx1)
   *
   * Output
   * ------
   * mat_c: this is the resulting matrix from the multiplication */

  int i, k; // The necessary counters


  // Initialize mat_c to have 0. in all elements
  for (i = 0; i < N; ++i)
    {
      mat_c[i] = 0;
    }

  // Do the multiplication and find the resulting matrix
  for (i = 0; i < N; ++i)
    {
      for (k = 0; k < N; ++k)
        {
	  // For each element mat_c[i][j] you must add the n_a multiplications
	  mat_c[i] += mat_a[i][k] * vec_b[k];
        }
    }

  // mat_c has been found so we are done
}
