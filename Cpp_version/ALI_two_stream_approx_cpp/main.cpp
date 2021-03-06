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
#include "runge_kutta_solver.h"
#include "dydx.h"
#include "subroutines.hpp"

using std::ofstream;
using namespace std;


#define ALT_RT
#ifdef ALT_RT

int main(int argc, char **argv)
{
  int iter;
  int niter = 100;
  int grid_size = 1000; // For saving variable outputs
  RT R; //automatically initialized grid

  //FILE * QHII_out;

  R.PrintGrid(stdout);

  //QHII_out = fopen("ion_frac_outfile.txt","w");

  //iterate
  for(iter=0;iter<niter;iter++)
  {
    //Compute I+, I-, J, and S
    R.UpdateIntensity();

    //print I+, I-, J, B, and S
    R.PrintIntensity(stdout);

    //Compute QHII, epsilon, alpha, and del_tau
    R.UpdateIonization();

    //print QHII, epsilon, alpha, and del_tau
    R.PrintIonization(stdout);

  }

  //save the final grid to a file
  R.PrintGrid(stdout);
  //R.PrintIonization(stdout);

  FILE *fp;
  char fname[200];
  sprintf(fname,"grid.txt");
  fp = fopen(fname,"w");
  R.SaveGrid(fp);
  fclose(fp);

  return 0;
}

#else  //ALT_RT

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
#define strict_solver

// Uncomment if the ICs are going to be read  in from a file
#define read_in_IC

// Uncomment if the final state of the solver is to be saved (and later set as  ICs)
#define write_outfile

// Uncomment if epsilon and/or B are varying across the grid or with time
// Automatically uses the ALI scheme to solve
#define vary_eps
//#define vary_eps_reg_IC
//#define reg_eps

// Uncomment if reading in density profile
#define read_density
//#define read_dens_set_IC

// Uncomment if evolving w.r.t. t at z=3 between a galaxy and UV background -- uses ALI automatically
// Note that read_density and vary_eps must be uncommented with this
//#define evolve
// If the strict solver, or LI scheme is desired instead
#define evolve_LI

int main()
{
  // Allocate the variables
  //--------------------------------
  int i,y; // Cell indices and counters
  int n=70;  // Number of cells
  int index; // To keep track of values while reading in a file

  double del_z = 1./float(n); // Grid step size
  double j;                   // Double precision counter
  double S_val_plus;          // Double precision value for the S interpolation used to solve for I+
  double S_val_minus;         // Double precision value for the S interpolation used to solve for I-

  double s; // Used for reading values in a file

  bool Plus; // The boolean used for interpolation of the source function S

  ofstream data_file;     // Going to be used to create an output file
  ofstream I_p_outfile;   // Going to be used to create an output file that saves current state of dynamics
  ofstream I_m_outfile;   // Going to be used to create an output file that saves current state of dynamics
  ofstream J_outfile;     // Going to be used to create an output file that saves current state of dynamics 
  ofstream rho_outfile;   // Used to created an output file that saves the converted density profile
  ofstream del_t_outfile; // Output file that saves the change in optical depth profile
  ifstream inFile;        // Going to be used to read in the state of dynamics

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
  double *eps_grid = new double[n]; // Actual grid of epsilon values  
#endif

#ifdef vary_eps_reg_IC 

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
#endif

#ifdef reg_eps
  double epsilon; // Photon destruction probability
  epsilon = 1.;
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
#ifdef vary_eps_reg_IC
      S[i] = eps_grid[i]*B[i] + (1.-eps_grid[i])*J[i];
#endif
#ifdef reg_eps
      S[i] = epsilon*B[i] + (1.-epsilon)*J[i];
#endif
    }
#endif

#ifdef vary_eps_reg_IC 
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

  // Setting the optical depth array
#ifdef read_density
  // Set a density grid and other necessary values 
  double *density = new double[n];       // Density profile  
  double *col_dens = new double[n];      // Column density for H1
  double sigma_nu = 7.91*pow(10.,-18.);  // Photoionization of H1
  double pi = 3.14159265358979323846;    // Value of pi
  double delta_x;                        // Physical size of grid cell at z=3
  double Delta;                          // Volumetric ratio to scale cell sizes
  double kpc_cm = 3.086*pow(10.,21.);    // Cm in one kpc
  double M_sun_g = 1.989*pow(10.,33.);   // Mass of the sun in g
  double M_H_g = 1.6726*pow(10.,-24.);   // Mass of H1 in grams
  double redshift = 3.;                  // Redshift of the density profile
  double lil_h = 0.6766;                 // Scaled Hubble constant
  double density_conversion;             // Unit conversion constant for the density distribution
  double I_p_init;                       // Initial I+ value
  double I_m_init;                       // Initial I- value
  double I_scale;                        // Scaling value of the specific intensity

  double *del_tau_int_init_p = new double[n]; // Initial I+ and I- integration values
  double *del_tau_int_init_m = new double[n]; // Initial I+ and I- integration values
  
  FILE * I_lum;

  // Read in the density profile .txt file
  inFile.open("density_profile.txt");
  // If there is no file with that name raise out of code
      if (!inFile)
        {
	  cerr << "density_profile.txt does not exist!"; // Tell the user what happened
	  exit(1); // Stop the system
        }

      index = 0;
      while (inFile >> s)
	{
	  density[index] = s; // Read in density profile element by element
	  index = index + 1;
	}

  inFile.close();

  density_conversion = (M_sun_g/(kpc_cm*kpc_cm*kpc_cm))*(1./M_H_g)*((redshift+1.)*(redshift+1.)*(redshift+1.))*
    (lil_h*lil_h);
  delta_x = (24.4/lil_h)*(1./(redshift+1.))*(kpc_cm); // Size of one grid cell in cm
  //Delta = (4.*pi)/(9.016*9.016*9.016);            // Volume ratio for scaling

  I_scale = 10.e-17;
  I_p_init = 10.e-17/I_scale; // The scaled bottom-up specific intensity; 1 --> 10^-17*Delta
  I_m_init = 10.e-22/I_scale; // The scaled top-down specific intensity; 10^-22*Delta/10^-17*Delta
  I_m_init = 0.;

  I_lum = fopen("I_luminosity.txt","w");
  
  for (i=0;i<n;i++)
    {
      j = (double) (i);
      z[i] = j*del_z;                             // Setting the grid
      B[i] = 0.;                                  // Thermal source function
      density[i] = density[i]*density_conversion; // Scaled density field 
      //density[i] = 1e-6;
      if (density[i] > 1.e-2)
	{
	  density[i] = 1.e-2;
	}
      col_dens[i] = density[i]*delta_x;           // Scaled column density
      del_tau[i] = col_dens[i]*sigma_nu;          // Change in optical depth
      //del_tau[i] = 0.;
      j = float(i) + 1.;
      //Delta = pow(((j-1)/j),2.);
#ifdef read_dens_set_IC 
      if (i == 0)
	{
	  I_plus[i] = I_p_init;                   // Bottom-up specific intensity
	  fprintf(I_lum, "%2.8f\n",I_plus[i]);
	}
      else 
	{
	  I_plus[i] = I_plus[0]/pow(j,2.); // Bottom-up specific intensity
	  fprintf(I_lum, "%2.8f\n",I_plus[i]);
	}
      I_minus[i] = I_m_init;                      // Top-down specific intensity
      J[i] = 0.5*(I_plus[i]+I_minus[i]);          // Mean specific intensity
      // Source function
#ifdef reg_eps
      S[i] = epsilon*B[i] + (1.-epsilon)*J[i];
#endif
#endif
      //printf("del x %2.5f\n",delta_x);
    }
  fclose(I_lum);
#else
  // Set the array of differences in optical depth for each cell by calling its function
  tau_fn(del_tau,z,del_z,n);
#endif

  // Set the array of third order quadratic interpolation coefficients
  quad_int(u_p,p_p,d_p,u_m,p_m,d_m,del_tau,n);

  //for(int k=0;k<n;k++)
  //{
  //printf("pp %2.5f\n",p_p[k]);
  //}

#ifdef evolve
  double t;
  double h=0.1;           // This is the time step by which we evolve the ODE for QHII
  double sigma_T = 0.66524587158*pow(10.,-24.); // Thomson cross section
  
  ofstream ion_frac_file; // Going to be used to create an output file for QHII
  ofstream eps_grid_file; // Saves the evolving values of epsilon across the grid

  double *QHII; // Ionized fraction of HII
  QHII = (double *) std::calloc(n,sizeof(double *));

  // Initialize the ionized fraction
  for (i = 0; i < n; i++)
    {
      //printf("Q %2.5f\n",QHII[i]);
      eps_grid[i] = (sigma_nu*(1. - QHII[i]))/(sigma_nu*(1. - QHII[i]) + (sigma_T*QHII[i]));
      S[i] = eps_grid[i] * B[i] + (1. - eps_grid[i]) * J[i];
      //printf("S_s %2.5f\n",S[i]);
    }

  // Open a file to keep track of the ionized fraction
  // Open a file to keep track of the source function
  data_file.open("S_outfile.txt");
  ion_frac_file.open("ion_frac_outfile.txt");
  eps_grid_file.open("eps_outfile.txt");
  for (i = 0; i < n; i++)
    {
      data_file << S[i];
      data_file << "\n";

      ion_frac_file << QHII[i];
      ion_frac_file << "\n";

      eps_grid_file << eps_grid[i];
      eps_grid_file << "\n";
    }
  data_file << "break\n";
  ion_frac_file << "break\n";
  eps_grid_file << "break\n";

  // Solve for converged S
  // The 1D diffusion equation is a 2nd order PDE where lambda_star*
  // S is the solution
  // The numerical representation used for this PDE is Central Space
  // Forward substitution + backward substitution are used to solve the PDE
  //--------------------------
  for (t = 0.; t < 5.; t=t+h)
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
          if (t>0)
            {
              if (i==1)
                {
                  data_file << S[0];
                  data_file << "\n";
		  ion_frac_file << QHII[0];
		  ion_frac_file << "\n";
		  eps_grid_file << eps_grid[0];
		  eps_grid_file << "\n";
                }
	      else if (i==n-1)
                {
                  data_file << S[i];
                  data_file << "\n";
                  data_file << "break\n";
		  ion_frac_file << QHII[i];
		  ion_frac_file << "\n";
		  ion_frac_file << "break\n";
		  eps_grid_file << eps_grid[i];
                  eps_grid_file << "\n";
		  eps_grid_file << "break\n";
                }
              else
                {
                  data_file << S[i];
                  data_file << "\n";
		  ion_frac_file << QHII[i];
		  ion_frac_file << "\n";
		  eps_grid_file << eps_grid[i];
                  eps_grid_file << "\n";
                }
            }
        }
      // Now use the tri-diagonal solver to actually update S
      tri_solver(M_a,M_b,M_c,y_arr,S,n);

      // Solve for the new value of epsilon by solving the ODE for QHII
      RK4_solver(t, QHII, h, density, n, J, QHII);
      for (int c = 0; c < n; c++)
	{
	  eps_grid[c] = (sigma_nu*(1. - QHII[c]))/(sigma_nu*(1. - QHII[c]) + (sigma_T*QHII[c]));
	}
      S[0] = eps_grid[0] * B[0] + (1. - eps_grid[0]) * J[0];
    }
  data_file.close();
  ion_frac_file.close();
  eps_grid_file.close();
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
  for (y = 0; y < 101; y++)
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
	  j = float(i)+1.;
	  I_plus[i] = exp(-del_tau[i]) * I_plus[i-1] * pow(((j-1)/j),2.) + S_val_plus;
          I_minus[i] = exp(-del_tau[i+1]) * I_minus[i+1] + S_val_minus;

#ifdef read_density
	  if (i == n-1)
	    {
	      I_minus[i] = I_m_init;
	    }
#endif
          // Update the mean specific intensity from both rays
          // Ch.4, section 4.4.5, equation 4.56
          J[i] = 0.5 * (I_plus[i] + I_minus[i]);

          // Solve for each component of the y array in MX = y
          // Ch.4, section 4.4.4, equation 4.48
          y_arr[i] = epsilon*B[i] + (1. - epsilon) * (J[i] - lmbda_s_S[i]);
	  
	  //printf("y %2.5f\n",(J[i] - lmbda_s_S[i]));

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
  double t;
#ifdef evolve_LI
  double h=0.1;           // This is the time step by which we evolve the ODE for QHII
  double sigma_T = 0.66524587158*pow(10.,-24.); // Thomson cross section
  
  //ofstream ion_frac_file; // Going to be used to create an output file for QHII
  //ofstream eps_grid_file; // Saves the evolving values of epsilon across the grid
  FILE * Qfile;
  FILE * eps_file;
  FILE * Sfile;

  double *QHII; // Ionized fraction of HII
  QHII = (double *) std::calloc(n,sizeof(double *));

  // Initialize the ionized fraction
  for (i = 0; i < n; i++)
    {
      //printf("Q %2.5f\n",QHII[i]);
      eps_grid[i] = (sigma_nu*(1. - QHII[i]))/(sigma_nu*(1. - QHII[i]) + (sigma_T*QHII[i]));
      S[i] = eps_grid[i] * B[i] + (1. - eps_grid[i]) * J[i]*I_scale;
      //printf("eps %2.5f\n",eps_grid[i]);
    }

  // Open a file to keep track of the ionized fraction
  // Open a file to keep track of the source function
  Qfile = fopen("ion_frac_outfile.txt","w");
  eps_file = fopen("eps_outfile.txt","w");
  Sfile = fopen("S_outfile.txt","w");
  for (i = 0; i < n; i++)
    {
      fprintf(Sfile, "%2.5f\n", S[i]);

      fprintf(Qfile, "%2.5f\n",QHII[i]);

      fprintf(eps_file, "%2.5f\n",eps_grid[i]);
    }
  fprintf(Sfile, "break\n");
  fprintf(Qfile, "break\n");
  fprintf(eps_file, "break\n");
#endif
#ifdef reg_IC
  data_file.open("output_file.txt");
  data_file << epsilon;
  data_file << "\n";
#endif

  // Convert I+, I-, and J back into actual units before starting calculations
  for (i=0;i<n;i++)
    {
      I_plus[i] = I_plus[i]*I_scale;
      I_minus[i] = I_minus[i]*I_scale;
      J[i] = J[i]*I_scale;
    }

  for (t = 0.; t < 10.; t=t+h) 
    {
      S[0] = eps_grid[0] * B[0] + (1. - eps_grid[0]) * J[0];
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
#ifdef evolve_LI
      // Integrating del_tau for the I+ and I- ray respectively
      //if (t >= 9.9)
      //{
      //del_tau_int_init_p[i] = h*(0.5*(del_tau[i+1] + del_tau[i]));
      //del_tau_int_init_m[i] = h*(0.5*(del_tau[i] + del_tau[i+1]));
	//printf("del_tau %2.5f\n",del_tau_int_init_m[i]);
      //}
      // Solve the ODE to update QHII
      RK4_solver(t, QHII, h, density, n, J, QHII);

      // Update the eps grid with the new values for QHII
      for (int c = 0; c < n; c++)
        {
          eps_grid[c] = (sigma_nu*(1. - QHII[c]))/(sigma_nu*(1. - QHII[c]) + (sigma_T*QHII[c]));
        }
      // Update the source function using the new value in the eps_grid
      S[i] = eps_grid[i] * B[i] + (1. - eps_grid[i]) * J[i];
      //S[0] = eps_grid[0] * B[0] + (1. - eps_grid[0]) * J[0];
#else
      S[i] = epsilon * B[i] + (1. - epsilon) * J[i];
      S[0] = epsilon * B[0] + (1. - epsilon) * J[0];
#endif
      //Write S to the file for plotting
      if (i==1)
	{
	  fprintf(Sfile, "%2.5f\n", S[0]);
#ifdef evolve_LI
	  fprintf(Qfile, "%2.5f\n", QHII[0]);
	  fprintf(eps_file, "%2.5f\n", eps_grid[0]);
#endif
	}
      else if (i==n-1){
	fprintf(Sfile, "%2.5f\n", S[i]);
	fprintf(Sfile, "break\n");
#ifdef evolve_LI
	fprintf(Qfile, "%2.5f\n", QHII[i]);
	fprintf(Qfile, "break\n");
	fprintf(eps_file, "%2.5f\n", eps_grid[i]);
	fprintf(eps_file, "break\n");
#endif
      }
      else{
	fprintf(Sfile, "%2.5f\n", S[i]);
#ifdef evolve_LI
	fprintf(Qfile, "%2.5f\n", QHII[i]);
	fprintf(eps_file, "%2.5f\n", eps_grid[i]);
#endif	
      }
    }
  }
  fclose(Sfile);
#ifdef evolve_LI
  fclose(Qfile);
  fclose(eps_file);
#endif
#endif
  for(int k=0;k<n;k++)
    {
      printf("S %2.5f\n",S[k]);
      //printf("del_tau %2.5f\n",del_tau[k]);
    }

  // Write final values of I+ and I- to different output files to read in later
#ifdef write_outfile
  I_p_outfile.open("I_p_output_final.txt");
  I_m_outfile.open("I_m_output_final.txt");
  J_outfile.open("J_output_final.txt");
  del_t_outfile.open("del_t_output_final.txt");
  // If density file is needed uncomment the density-related lines
  rho_outfile.open("rho_output.txt");
  for (i = 0; i < n; i++) {
    I_p_outfile << I_plus[i];
    I_p_outfile << "\n";
    I_m_outfile << I_minus[i];
    I_m_outfile << "\n";
    J_outfile << J[i];
    J_outfile << "\n";
    del_t_outfile << del_tau[i];
    del_t_outfile << "\n";
    if (i == 0)
      {
	rho_outfile << delta_x/kpc_cm;
	rho_outfile << "\n";
      }
    else
      {
	rho_outfile << density[i-1];
	rho_outfile << "\n";
      }
  }
  rho_outfile << density[n-1];
  rho_outfile << "\n";

  I_p_outfile.close();
  I_m_outfile.close();
  J_outfile.close();
  del_t_outfile.close();
  rho_outfile.close();
#endif

  // Free the memory
#ifdef vary_eps
  delete[](eps_grid);
#endif
#ifdef read_density
  delete[](density);
  delete[](col_dens);
  delete[](del_tau_int_init_p);
  delete[](del_tau_int_init_m);
#endif
#ifdef evolve
  free(QHII);
#endif
#ifdef evolve_LI
  free(QHII);
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

#endif //ALT_RT
