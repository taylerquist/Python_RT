#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include "tau_fn.h"
#include "quad_int.h"
#include "matrix_build.h"
#include "s_interp.h"
#include "tri_solver.h"
#include "mat_mult.h"

int main()
{
  // Allocate the variables
  //--------------------------------
  int i,t;     // Cell indices and counters
  int n=70;    // Number of cells
  int N=n+1;   // Number of grid points

  double del_z = 1./n;     // Grid step size
  double epsilon = 0.0001; // Photon destruction probability should be equal to 0.1, 0.001, or 0.0001
  double j;                // Double precision counter
  double S_val_plus;       // Double precision value for the S interpolation used to solve for I+
  double S_val_minus;      // Double precision value for the S interpolation used to solve for I-

  bool Plus; // The boolean used for interpolation of the source function S

  double *B = new double[n];       // The thermal source function (aka the Planck Function)
  double *del_tau = new double[N]; // Change in optical depth
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
  double *S = new double[n];       // The source function
  //double *J = new double[n];       // The mean specific intensity
  double *z = new double[n];       // The grid array

  double **lmbda_s;                // This is the lambda_star matrix
  double **lmbda_s_J;              // This is the matrix multiplication product of the lambda_star matrix with J
  double **J;              // In order to do matrix multiplication, J must be a column vector
  lmbda_s = (double **) std::malloc(n * sizeof(double *));
  lmbda_s_J = (double **) std::malloc(n * sizeof(double *));
  J = (double **) std::malloc(n * sizeof(double *));
  for(i=0;i<n;i++)
    {
      lmbda_s[i] = (double *) std::malloc(n * sizeof(double));   // Lambda_star is an nxn matrix
      lmbda_s_J[i] = (double *) std::malloc(1 * sizeof(double)); // Ends up being a column vector since nxn * nx1
      J[i] = (double *) std::malloc(1 * sizeof(double)); // Column vector, nx1
    }


  // Initialize the variables
  //--------------------------------
  z[0] = 0.;                               // The grid increases by step size del_z and goes from 0 to 1
  B[0] = 1.;                               // The thermal source function in this example is constant in this example
  alpha[0] = B[i]*pow(10.,(5.-6.*z[i]));   // The emissivity, depends on z Ch.4, section 4.4.2, equation 4.32
  I_plus[0] = 1.;                          // The bottom-up specific intensity
  I_minus[0] = 0.;                         // The top-down specific intensity
  J[0][0] = 1.;                    // The average specific intensity
  S[0] = epsilon*B[0] + (1.-epsilon)*J[0][0]; // The source function Ch.4 section 4.4.4, equation 4.37
  S_val_plus = 0.;                         // Will later be used for interpolation for the bottom-up, I+ array
  S_val_minus = 0.;                        // Will later be used for interpolation for the top-down, I-, array

  for (i=1;i<n;i++)
    {
      j = (double) (i);
      z[i] = j*del_z;
      B[i] = 1.;
      alpha[i] = B[i]*pow(10.,(5.-6.*z[i]));
      I_plus[i] = 0.;
      I_minus[i] = 0.;
      J[i][0] = 0.;
      S[i] = epsilon*B[i] + (1.-epsilon)*J[i][0];
    }

  // Set other grids, interpolation coefficients, and matrices
  //----------------------------------------------------------------
  // Set the array of differences in optical depth for each cell by calling its function
  tau_fn(del_tau,z,del_z,N);

  // Set the array of third order quadratic interpolation coefficients
  quad_int(u_p,p_p,d_p,u_m,p_m,d_m,del_tau,N);

  // Set the tri-diagonal arrays of the lambda_star and M_star matrix and the lambda_star full matrix
  // Currently not working because of trouble with defining matrix in main and in function input
  matrix_build(lmbda_s,M_a,M_b,M_c,u_p,p_p,d_p,u_m,p_m,d_m,epsilon,n);


  // Solve for converged S
  // The 1D diffusion equation is a 2nd order PDE where lambda_star*S is the solution
  // The numerical representation used for this PDE is Central Space
  // Forward substitution + backward substitution are used to solve the PDE
  //--------------------------
  for (t = 1; t < 101; t++)
    {
      // First actually solve the matrix multiplication that will be used later to solve for the y in MX = y
      // It needs to be used in the next loop but it does not need to be re-calculated every time
      mat_mult(lmbda_s,J,lmbda_s_J,n);
      for (i = 1; i < n; i++)
        {
	  // First solve the bottom-up interpolated value of S
	  Plus = true;
	  S_interp(S_val_plus,alpha,S,z,u_p,p_p,d_p,u_m,p_m,d_m,del_z,i,Plus);

	  // First solve the top-down interpolated value of S
	  Plus = false;
	  S_interp(S_val_minus,alpha,S,z,u_p,p_p,d_p,u_m,p_m,d_m,del_z,i,Plus);

	  // Solve for the (interpolated) integration of the specific intensity for each ray
	  // The equations for this come from Ch.4, section 4.4.5, equations 4.53 and 4.54
	  I_plus[i] = exp(-del_tau[i-1]) * I_plus[i-1] + S_val_plus;
	  I_minus[i] = exp(-del_tau[i]) * I_plus[i] + S_val_minus;

	  // Update the mean specific intensity from both rays
	  // Ch.4, section 4.4.5, equation 4.56
	  J[i][0] = 0.5 * (I_plus[i] + I_minus[i]);

	  // Solve for each component of the y array in MX = y
	  // Ch.4, section 4.4.4, equation 4.48
	  // Might need to do the first component outside of here...
	  y_arr[i] = epsilon*B[t] + (1. - epsilon); // * (J[t] - ls_S[t])
        }
      // Now use the tri-diagonal solver to actually update S
      tri_solver(M_a,M_b,M_c,y_arr,S,n);
      S[0] = 1.;
    }

  return 0;
}
