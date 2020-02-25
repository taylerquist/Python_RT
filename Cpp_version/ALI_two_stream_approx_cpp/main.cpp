#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <fstream>
#include "tau_fn.h"
#include "quad_int.h"
#include "matrix_build.h"
#include "matrix_build_vector.h"
#include "s_interp.h"
#include "tri_solver.h"
#include "mat_mult.h"

using std::ofstream;
using namespace std;

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


// Uncomment this if using ALI with constant epsilon  and B
//#define ALI

// Uncomment if epsilon and B are constant and ICs are not getting read in
//#define reg_IC

// Uncomment if no lambda acceleration is to be  used
//#define strict_solver

// Uncomment if the ICs are going to be read  in from a file
//#define read_in_IC

// Uncomment if the final state of the solver is to be saved (and later set as  ICs)
//#define write_outfile

// Uncomment if epsilon and/or B are varying across the grid or with time
// Automatically uses the ALI scheme to solve
#define vary_eps

int main()
{
  // Allocate the variables
  //--------------------------------
  int i,y;   // Cell indices and counters
  int n=70;  // Number of cells
  int index; // To keep track of values while reading in a file

  double del_z = 1./float(n); // Grid step size
  double j;                   // Double precision counter
  double S_val_plus;          // Double precision value for the S interpolation used to solve for I+
  double S_val_minus;         // Double precision value for the S interpolation used to solve for I-

  double s; // Used for reading values in a file

  bool Plus; // The boolean used for interpolation of the source function S

  ofstream data_file; // Going to be used to create an output file
  ofstream I_p_outfile;   // Going to be used to create an output file that saves current state of dynamics
  ofstream I_m_outfile;   // Going to be used to create an output file that saves current state of dynamics
  ifstream inFile;    // Going to be used to read in the state of dynamics

  double *B = new double[n];       // The thermal source function (aka the Planck Function)
  double *del_tau = new double[n]; // Change in optical depth
  double *u_p = new double[n];     // "Plus" u interpolation coefficient array (the vector going from bottom up)
  double *p_p = new double[n];     // "Plus" p interpolation coefficient array (the vector going from bottom up)
  double *d_p = new double[n];     // "Plus" d interpolation coefficient array (the vector going from bottom up)
  double *u_m = new double[n];     // "Minus" u interpolation coefficient array (the vector going from top down)
  double *p_m = new double[n];     // "Minus" p interpolation coefficient array (the vector going from top down)
  double *d_m = new double[n];     // "Minus" d interpolation coefficient array (the vector going from top down)
  double *alpha_p = new double[n]; // "Plus" alpha interpolation coefficient array (the vector going from bottom up)
  double *beta_p = new double[n];  // "Plus" beta interpolation coefficient array (the vector going from bottom up)
  double *gamma_p = new double[n]; // "Plus" gamma interpolation coefficient array (the vector going from bottom up)
  double *alpha_m = new double[n]; // "Minus" alpha interpolation coefficient array (the vector going from top down)
  double *beta_m = new double[n];  // "Minus" beta interpolation coefficient array (the vector going from top down)
  double *gamma_m = new double[n]; // "Minus" gamma interpolation coefficient array (the vector going from top down)
  double *M_a = new double[n];     // Lower diagonal of the matrix M_star
  double *M_b = new double[n];     // Diagonal of the matrix M_star
  double *M_c = new double[n];     // Upper diagonal of the matrix M_star
  double *alpha = new double[n];   // Emissivity
  double *y_arr = new double[n];   // y array in the equation Mx = y
  double *I_plus = new double[n];  // The specific intensity for the vector going from bottom up
  double *I_minus = new double[n]; // The specific intensity for the vector going from top down
  double *z = new double[n];       // The grid array

  double **lmbda_s;                // This is the lambda_star matrix
  double *lmbda_s_S;               // This is the matrix multiplication product of the lambda_star matrix with S
  double *J;                       // Mean specific intensity, 1xn
  double *S;                       // The source function, 1xn
  lmbda_s = (double **) std::malloc(n * sizeof(double *));
  lmbda_s_S = (double *) std::malloc(n * sizeof(double *));
  J = (double *) std::malloc(n * sizeof(double *));
  S = (double *) std::malloc(n * sizeof(double *));
  for (i = 0; i < n; i++) {
    lmbda_s[i] = (double *) std::malloc(n * sizeof(double)); // Lambda_star is an nxn matrix
  }

  // Below determines if epsilon (photon destruction probability) will be held constant
  // throughout the grid or if it spatially varies
#ifdef vary_eps

  double l_bndry_eps = 0.001;       // Left boundary for epsilon grid (epsilon at z=0)
  double r_bndry_eps = 0.1;         // Right boundary for epsilon grid (epsilon at z=1)
  double del_eps;                   // The amount epsilon changes for every grid point
  double *eps_grid = new double[n]; // Actual grid of epsilon values

  // Set the spatially varying epsilon values
  del_eps = (l_bndry_eps/r_bndry_eps)*(1./float(n-1));
  eps_grid[0] = l_bndry_eps;
  for (i=1;i<n;i++)
    {
      eps_grid[i] = l_bndry_eps/(del_eps*float(i));
    }
#else
  double epsilon; // Photon destruction probability
  epsilon = 0.1;
#endif


  // Initialize the variables

  // If the initial conditions are being read in from an output file
#ifdef read_in_IC

  // Open the first file
  inFile.open("I_p_output_IC.txt");
  if (!inFile)
    {
      cerr << "I_p_output_IC.txt does not exist!"; // Tell the user what happened
      exit(1); // Stop the system
    }

  index = 0;
  while (inFile >> s)
    {
      I_plus[index] = s; // Read in specific intensity, bottom up array
      index = index + 1;
    }
  inFile.close();

  // Open the second file
  inFile.open("I_m_output_IC.txt");
  if (!inFile)
    {
      cerr << "I_m_output_IC.txt does not exist!"; // Tell the user what happened
      exit(1); // Stop the system
    }
  index = 0;
  while (inFile >> s)
    {
      I_minus[index] = s; // Read in specific intensity, top down array
      index = index + 1;
    }
  inFile.close();

  for (i=0;i<n;i++)
    {
      j = (double) (i);
      z[i] = j*del_z;
      B[i] = 1.;
      alpha[i] = B[i]*pow(10.,(5.-6.*z[i]));
      J[i] = 0.5*(I_plus[i]+I_minus[i]);
#ifdef vary_eps
      S[i] = eps_grid[i]*B[i] + (1.-eps_grid[i])*J[i];
#else
      S[i] = epsilon*B[i] + (1.-epsilon)*J[i];
#endif
    }
#endif

#ifdef vary_eps
  for (i=0;i<n;i++)
    {
      j = (double) (i);
      z[i] = j*del_z;
      B[i] = 1.;
      alpha[i] = B[i]*pow(10.,(5.-6.*z[i]));
      I_plus[i] = 0.;
      I_minus[i] = 0.;
      J[i] = 0.5*(I_plus[i]+I_minus[i]);
      if (i==0)
	{
	  I_plus[0] = 1.;
	  J[0] = 1.;
	}
      S[i] = eps_grid[i]*B[i] + (1.-eps_grid[i])*J[i];
    }

  // If the initial conditions are not going to be manually set
#endif

#ifdef reg_IC
  z[0] = 0.;                                 // The grid increases by step size del_z and goes from 0 to 1
  B[0] = 1.;                                 // The thermal source function in this example is constant
  alpha[0] = B[0]*pow(10.,(5.-6.*z[i]));     // The emissivity, depends on z Ch.4, section 4.4.2, equation 4.32
  I_plus[0] = 1.;                            // The bottom-up specific intensity
  I_minus[0] = 0.;                           // The top-down specific intensity
  J[0] = 1.;                                 // The average specific intensity
  S[0] = epsilon*B[0] + (1.-epsilon)*J[0];   // The source function Ch.4 section 4.4.4, equation 4.37

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
#endif

  // Set other grids, interpolation coefficients, and matrices
  //----------------------------------------------------------------
  // Set the array of differences in optical depth for each cell by calling its function
  tau_fn(del_tau,z,del_z,n);

  // Set the array of third order quadratic interpolation coefficients
  quad_int(u_p,p_p,d_p,u_m,p_m,d_m,del_tau,n);

#ifdef vary_eps
  // Open the file to write to and save the left (and right) boundary epsilon value and delta epsilon
  // for plotting purposes
  data_file.open("output_file.txt");
  data_file << l_bndry_eps;
  data_file << "\n";
  data_file << r_bndry_eps;
  data_file << "\n";
  data_file << del_eps;
  data_file << "\n";

  // Solve for converged S
  // The 1D diffusion equation is a 2nd order PDE where lambda_star*
  // S is the solution
  // The numerical representation used for this PDE is Central Space
  // Forward substitution + backward substitution are used to solve the PDE
  //--------------------------
  for (y = 0; y < 102; y++)
    {
      // Set the tri-diagonal arrays of the lambda_star and M_star matrix and the lambda_star full matrix
      matrix_build_vector(lmbda_s,M_a,M_b,M_c,u_p,p_p,d_p,u_m,p_m,d_m,eps_grid,n);

      // First actually solve the matrix multiplication that will be used later to solve for the y in MX = y
      // It needs to be used in the next loop but it does not need to be re-calculated every time
      mat_mult(lmbda_s,S,lmbda_s_S,n);
      for (i = 1; i < n; i++)
        {
          // First solve the bottom-up interpolated value of S
          Plus = true;
          S_val_plus = S_interp(alpha,S,u_p,p_p,d_p,u_m,p_m,d_m,del_z,i,Plus);
          // First solve the top-down interpolated value of S
          Plus = false;
          S_val_minus = S_interp(alpha,S,u_p,p_p,d_p,u_m,p_m,d_m,del_z,i,Plus);

          // Solve for the (interpolated) integration of the specific intensity for each ray
          // The equations for this come from Ch.4, section 4.4.5, equations 4.53 and 4.54
          I_plus[i] = exp(-del_tau[i]) * I_plus[i-1] + S_val_plus;
          I_minus[i] = exp(-del_tau[i+1]) * I_minus[i+1] + S_val_minus;

          // Update the mean specific intensity from both rays
          // Ch.4, section 4.4.5, equation 4.56
          J[i] = 0.5 * (I_plus[i] + I_minus[i]);

          // Solve for each component of the y array in MX = y
          // Ch.4, section 4.4.4, equation 4.48
          y_arr[i] = eps_grid[i]*B[i] + (1. - eps_grid[i]) * (J[i] - lmbda_s_S[i]);

          // Write S to the file for plotting
          // Only write to the file after the first updated iteration
          if (y>0)
            {
              if (i==1)
                {
                  data_file << S[0];
                  data_file << "\n";
                }
	      else if (i==n-1)
                {
                  data_file << S[i];
                  data_file << "\n";
                  data_file << "break\n";
                }
              else
                {
                  data_file << S[i];
                  data_file << "\n";
                }
            }
        }
      // Now use the tri-diagonal solver to actually update S
      tri_solver(M_a,M_b,M_c,y_arr,S,n);
      S[0] = eps_grid[0] * B[0] + (1. - eps_grid[0]) * J[0];
    }
  data_file.close();
#endif

#ifdef ALI
  // Open the file to write to and save epsilon, even if it spatially varies
  data_file.open("output_file.txt");
  data_file << epsilon;
  data_file << "\n";

  // Set the tri-diagonal arrays of the lambda_star and M_star matrix and the lambda_star full matrix
  // Currently not working because of trouble with defining matrix in main and in function input
  matrix_build(lmbda_s,M_a,M_b,M_c,u_p,p_p,d_p,u_m,p_m,d_m,epsilon,n);

  // Solve for converged S
  // The 1D diffusion equation is a 2nd order PDE where lambda_star*
  // S is the solution
  // The numerical representation used for this PDE is Central Space
  // Forward substitution + backward substitution are used to solve the PDE
  //--------------------------
  for (y = 0; y < 102; y++)
    {
      // First actually solve the matrix multiplication that will be used later to solve for the y in MX = y
      // It needs to be used in the next loop but it does not need to be re-calculated every time
      mat_mult(lmbda_s,S,lmbda_s_S,n);
      for (i = 1; i < n; i++)
        {
          // First solve the bottom-up interpolated value of S
          Plus = true;
          S_val_plus = S_interp(alpha,S,u_p,p_p,d_p,u_m,p_m,d_m,del_z,i,Plus);
          // First solve the top-down interpolated value of S
          Plus = false;
          S_val_minus = S_interp(alpha,S,u_p,p_p,d_p,u_m,p_m,d_m,del_z,i,Plus);

          // Solve for the (interpolated) integration of the specific intensity for each ray
          // The equations for this come from Ch.4, section 4.4.5, equations 4.53 and 4.54
          I_plus[i] = exp(-del_tau[i]) * I_plus[i-1] + S_val_plus;
          I_minus[i] = exp(-del_tau[i+1]) * I_minus[i+1] + S_val_minus;

          // Update the mean specific intensity from both rays
          // Ch.4, section 4.4.5, equation 4.56
          J[i] = 0.5 * (I_plus[i] + I_minus[i]);

          // Solve for each component of the y array in MX = y
          // Ch.4, section 4.4.4, equation 4.48
          y_arr[i] = epsilon*B[i] + (1. - epsilon) * (J[i] - lmbda_s_S[i]);

          // Write S to the file for plotting
          // Only write to the file after the first updated iteration
          if (y>0)
            {
              if (i==1)
                {
                  data_file << S[0];
                  data_file << "\n";
                }
	      else if (i==n-1)
                {
                  data_file << S[i];
                  data_file << "\n";
                  data_file << "break\n";
                }
              else
                {
                  data_file << S[i];
                  data_file << "\n";
                }
            }
        }
      // Now use the tri-diagonal solver to actually update S
      tri_solver(M_a,M_b,M_c,y_arr,S,n);
      S[0] = epsilon * B[0] + (1. - epsilon) * J[0];
    }
  data_file.close();
#endif
#ifdef strict_solver
  data_file.open("output_file.txt");
  data_file << epsilon;
  data_file << "\n";
  for (t = 0; t < 101; t++) {
    for (i = 1; i < n; i++) {
      // First solve the bottom-up interpolated value of S
      Plus = true;
      S_val_plus = S_interp(alpha, S, u_p, p_p, d_p, u_m, p_m, d_m, del_z, i, Plus);
      // First solve the top-down interpolated value of S
      Plus = false;
      S_val_minus = S_interp(alpha, S, u_p, p_p, d_p, u_m, p_m, d_m, del_z, i, Plus);

      // Solve for the (interpolated) integration of the specific intensity for each ray
      // The equations for this come from Ch.4, section 4.4.5, equations 4.53 and 4.54
      I_plus[i] = exp(-del_tau[i]) * I_plus[i-1] + S_val_plus;
      I_minus[i] = exp(-del_tau[i+1]) * I_minus[i+1] + S_val_minus;

      // Update the mean specific intensity from both rays
      // Ch.4, section 4.4.5, equation 4.56
      J[i] = 0.5 * (I_plus[i] + I_minus[i]);
      S[i] = epsilon * B[i] + (1. - epsilon) * J[i];

      //Write S to the file for plotting
      if (i==1)
	{
	  data_file << S[0];
	  data_file << "\n";
	}
      else if (i==n-1){
	data_file << S[i];
	data_file << "\n";
	data_file << "break\n";
      }
      else{
	data_file << S[i];
	data_file << "\n";
      }
    }
  }
  data_file.close();
#endif
  for(int k=0;k<n;k++)
    {
      printf("S %2.5f\n",S[k]);
    }

  // Write final values of I+ and I- to different output files to read in later
#ifdef write_outfile
  I_p_outfile.open("I_p_output_IC.txt");
  I_m_outfile.open("I_m_output_IC.txt");
  for (i = 0; i < n; i++) {
    I_p_outfile << I_plus[i];
    I_p_outfile << "\n";
    I_m_outfile << I_minus[i];
    I_m_outfile << "\n";
  }
  I_p_outfile.close();
  I_m_outfile.close();
#endif

  // Free the memory
#ifdef vary_eps
  delete[](eps_grid);
#endif
  free(lmbda_s);
  free(lmbda_s_S);
  free(J);
  free(S);
  delete[](B);
  delete[](del_tau);
  delete[](alpha_m);
  delete[](beta_m);
  delete[](gamma_m);
  delete[](alpha_p);
  delete[](beta_p);
  delete[](gamma_p);
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
