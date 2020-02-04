#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <fstream>
#include "tau_fn.h"
#include "quad_int.h"
#include "matrix_build.h"
#include "s_interp.h"
#include "tri_solver.h"
#include "mat_mult.h"

using std::ofstream;

// ALI numerical algorithm for the two stream approximation
/* Notes references
 * -----------------
 * Dullemond Radiative Transfer class notes, chapters 3 and 4:
 * http://www.ita.uni-heidelberg.de/~dullemond/lectures/radtrans_2012/Chapter_3.pdf
 * http://www.ita.uni-heidelberg.de/~dullemond/lectures/radtrans_2012/Chapter_4.pdf */

/* The following is the general tri-diagonal matrix layout, thus each diagonal component of the M_star and lambda_star
 * will be referred to as M_a, M_b,..., etc
 * | b[0] c[0]  0   ...  0  |
 * | a[1] b[1] c[1]  0  ... |
 * |  0   a[2] b[2] c[2] ...|
 * |  .    .    .    .    . |
 * |  .    .    .    .    . |
 * |  .    .    .    .    . |
 * |...a[n-2] b[n-2] c[n-2] |
 * |...  ...  a[n-1] b[n-1] |
 * Note: a[0] = 0 and c[n-1] = 0 */

int main()
{
  // Allocate the variables
  //--------------------------------
  int i,t;   // Cell indices and counters
  int n=70;  // Number of cells
  //  int N=n+1; // Number of grid points

  double del_z = 1./float(n);     // Grid step size
  double epsilon = 0.001; // Photon destruction probability should be equal to 0.1, 0.001, or 0.0001
  double j;                // Double precision counter
  double S_val_plus;       // Double precision value for the S interpolation used to solve for I+
  double S_val_minus;      // Double precision value for the S interpolation used to solve for I-

  bool Plus; // The boolean used for interpolation of the source function S

  ofstream data_file; // Going to be used to create an output file

  double *B = new double[n];       // The thermal source function (aka the Planck Function)
  double *del_tau = new double[n]; // Change in optical depth
  double *u_p = new double[n];     // "Plus" u interpolation coefficient array (the vector going from bottom up)
  double *p_p = new double[n];     // "Plus" p interpolation coefficient array (the vector going from bottom up)
  double *d_p = new double[n];     // "Plus" d interpolation coefficient array (the vector going from bottom up)
  double *u_m = new double[n];     // "Minus" u interpolation coefficient array (the vector going from top down)
  double *p_m = new double[n];     // "Minus" p interpolation coefficient array (the vector going from top down)
  double *d_m = new double[n];     // "Minus" d interpolation coefficient array (the vector going from top down)
  double *M_a = new double[n];     // Lower diagonal of the matrix M_star
  double *M_b = new double[n];     // Diagonal of the matrix M_star
  double *M_c = new double[n];     // Upper diagonal of the matrix M_star
  double *alpha = new double[n];   // Emissivity
  double *y_arr = new double[n];   // y array in the equation Mx = y
  double *I_plus = new double[n];  // The specific intensity for the vector going from bottom up
  double *I_minus = new double[n]; // The specific intensity for the vector going from top down
  double *z = new double[n];       // The grid array

  double **lmbda_s;                // This is the lambda_star matrix
  double *lmbda_s_S;              // This is the matrix multiplication product of the lambda_star matrix with S
  double *J;                      // In order to do matrix multiplication, J must be a column vector
  double *S;                      // The source function
  lmbda_s = (double **) std::malloc(n * sizeof(double *));
  lmbda_s_S = (double *) std::malloc(n * sizeof(double *));
  J = (double *) std::malloc(n * sizeof(double *));
  S = (double *) std::malloc(n * sizeof(double *));
  for (i = 0; i < n; i++) {
    lmbda_s[i] = (double *) std::malloc(n * sizeof(double));   // Lambda_star is an nxn matrix
  }

  // Initialize the variables
  //--------------------------------
  z[0] = 0.;                                    // The grid increases by step size del_z and goes from 0 to 1
  B[0] = 1.;                                    // The thermal source function in this example is constant in this example
  alpha[0] = B[0]*pow(10.,(5.-6.*z[i]));        // The emissivity, depends on z Ch.4, section 4.4.2, equation 4.32
  I_plus[0] = 1.;                               // The bottom-up specific intensity
  I_minus[0] = 0.;                              // The top-down specific intensity
  J[0] = 1.;                                 // The average specific intensity
  S[0] = epsilon*B[0] + (1.-epsilon)*J[0];// The source function Ch.4 section 4.4.4, equation 4.37

  for (i=1;i<n;i++)
    {
      j = (double) (i);
      z[i] = j*del_z;
      B[i] = 1.;
      alpha[i] = B[i]*pow(10.,(5.-6.*z[i]));
      I_plus[i] = 0.;
      I_minus[i] = 0.;
      J[i] = 0.;
      S[i] = epsilon*B[i] + (1.-epsilon)*J[i];
    }

  // Set other grids, interpolation coefficients, and matrices
  //----------------------------------------------------------------
  // Set the array of differences in optical depth for each cell by calling its function
  tau_fn(del_tau,z,del_z,n);

  // Set the array of third order quadratic interpolation coefficients
  quad_int(u_p,p_p,d_p,u_m,p_m,d_m,del_tau,n);

  // Set the tri-diagonal arrays of the lambda_star and M_star matrix and the lambda_star full matrix
  // Currently not working because of trouble with defining matrix in main and in function input
  matrix_build(lmbda_s,M_a,M_b,M_c,u_p,p_p,d_p,u_m,p_m,d_m,epsilon,n);

  // Solve for converged S
  // The 1D diffusion equation is a 2nd order PDE where lambda_star*S is the solution
  // The numerical representation used for this PDE is Central Space
  // Forward substitution + backward substitution are used to solve the PDE
  //--------------------------
  for (t = 0; t < 101; t++)
    {
      // First actually solve the matrix multiplication that will be used later to solve for the y in MX = y
      // It needs to be used in the next loop but it does not need to be re-calculated every time
      mat_mult(lmbda_s,S,lmbda_s_S,n);
      for (i = 1; i < n-1; i++)
        {
	  // First solve the bottom-up interpolated value of S
	  Plus = true;
	  S_val_plus = S_interp(alpha,S,u_p,p_p,d_p,u_m,p_m,d_m,del_z,i,Plus);
	  // First solve the top-down interpolated value of S
	  Plus = false;
	  S_val_minus = S_interp(alpha,S,u_p,p_p,d_p,u_m,p_m,d_m,del_z,i,Plus);

	  // Solve for the (interpolated) integration of the specific intensity for each ray
	  // The equations for this come from Ch.4, section 4.4.5, equations 4.53 and 4.54
	  I_plus[i] = exp(-del_tau[i-1]) * I_plus[i-1] + S_val_plus;
	  I_minus[i] = exp(-del_tau[i]) * I_minus[i+1] + S_val_minus;

	  // Update the mean specific intensity from both rays
	  // Ch.4, section 4.4.5, equation 4.56
	  J[i] = 0.5 * (I_plus[i] + I_minus[i]);

	  // Solve for each component of the y array in MX = y
	  // Ch.4, section 4.4.4, equation 4.48
	  // *** This is where things really start to differ between the Python code ***
	  y_arr[i] = epsilon*B[i] + (1. - epsilon) * (J[i] - lmbda_s_S[i]);
        }
      // Now use the tri-diagonal solver to actually update S
      tri_solver(M_a,M_b,M_c,y_arr,S,n);
    }

  data_file.open("output_file.txt");
  for (i = 0; i < n; ++i)
    {
      // Write to the file in a coherent way
      if (i == 0)
        {
	  data_file << epsilon;
	  data_file << "\n";
        }
      data_file << S[i];
      data_file << "\n";
    }
  data_file.close();

  // Free the memory
  free(lmbda_s);
  free(lmbda_s_S);
  free(J);
  free(S);
  delete[](B);
  delete[](del_tau);
  delete[](u_m);
  delete[](p_m);
  delete[](d_m);
  delete[](u_p);
  delete[](p_p);
  delete[](d_p);
  delete[](M_a);
  delete[](M_b);
  delete[](M_c);
  delete[](alpha);
  delete[](y_arr);
  delete[](I_minus);
  delete[](I_plus);
  delete[](z);

  printf("done \n");

  return 0;
}
