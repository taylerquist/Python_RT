#ifndef SUBROUTINES_H
#define SUBROUTINES_H

#include <stdio.h>

class RT {
  public:
  	int n; //number of grid cells
    double sigma_nu;// = 7.91e-18;  // Photoionization of H1
    double pi;// = 3.14159265358979323846;    // Value of pi
    double delta_x;                        // Physical size of grid cell at z=3
    double Delta;                          // Volumetric ratio to scale cell sizes
    double kpc_cm;// = 3.086e21;    // Cm in one kpc
    double M_sun_g;// = 1.989e33;   // Mass of the sun in g
    double M_H_g;// = 1.6726e-24;   // Mass of H1 in grams
    double redshift;// = 3.;                  // Redshift of the density profile
    double lil_h;// = 0.6766;                 // Scaled Hubble constant
    double density_conversion; //convert between grid (Msun/kpc^3) and cm^-3
    double density_limit; //maximum density allowed
    double del_z; //normalized grid size
    double *del_tau_int_init;
    double epsilon;
    double sigma_T;

    double I_p_init;
    double I_m_init;

    // Set a density grid and other necessary values 
    double *density;// = new double[n];       // Density profile  
    double *col_dens;// = new double[n];      // Column density for H1
    double *B;// = new double[n];       // The thermal source function (aka the Planck Function)
    double *del_tau;// = new double[n]; // Change in optical depth
    double *u_p;// = new double[n];     // "Plus" u interpolation coefficient array (the vector going from bottom up)
    double *p_p;// = new double[n];     // "Plus" p interpolation coefficient array (the vector going from bottom up)
    double *d_p;// = new double[n];     // "Plus" d interpolation coefficient array (the vector going from bottom up)
    double *u_m;// = new double[n];     // "Minus" u interpolation coefficient array (the vector going from top down)
    double *p_m;// = new double[n];     // "Minus" p interpolation coefficient array (the vector going from top down)
    double *d_m;// = new double[n];     // "Minus" d interpolation coefficient array (the vector going from top down)
    double *alpha_p;// = new double[n]; // "Plus" alpha interpolation coefficient array (the vector going from bottom up)
    double *beta_p;// = new double[n];  // "Plus" beta interpolation coefficient array (the vector going from bottom up)
    double *gamma_p;// = new double[n]; // "Plus" gamma interpolation coefficient array (the vector going from bottom up)
    double *alpha_m;// = new double[n]; // "Minus" alpha interpolation coefficient array (the vector going from top down)
    double *beta_m;// = new double[n];  // "Minus" beta interpolation coefficient array (the vector going from top down)
    double *gamma_m;// = new double[n]; // "Minus" gamma interpolation coefficient array (the vector going from top down)
    double *M_a;// = new double[n];     // Lower diagonal of the matrix M_star
    double *M_b;// = new double[n];     // Diagonal of the matrix M_star
    double *M_c;// = new double[n];     // Upper diagonal of the matrix M_star
    double *alpha;// = new double[n];   // Emissivity
    double *y_arr;// = new double[n];   // y array in the equation Mx = y
    double *I_plus;// = new double[n];  // The specific intensity for the vector going from bottom up
    double *I_minus;// = new double[n]; // The specific intensity for the vector going from top down
    double *z;// = new double[n];       // The grid array
    double *xg; //grid locations in kpc
    double *QHII;
    double *J;
    double *S;
    double *eps_grid;
	/*! \fn Constants(void)
	 *  \brief Constructor for the RT class */
	RT(void);

	void SetConstants(void);
	void InitializeGrid(void);
    void PrintGrid(FILE *fp);
    void PrintIntensity(FILE *fp);
    void PrintIonization(FILE *fp);
    void SaveGrid(FILE *fp);

    void UpdateIntensity(void);
    void UpdateIonization(void);

};
#endif // SUBROUTINES_H