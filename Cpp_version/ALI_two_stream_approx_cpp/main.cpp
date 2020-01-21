#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <cmath>

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


// Define all functions needed for the main program
//--------------------------------------------------

// Difference in optical depth function
void tau_fn(double *delta_tau, double *z_arr, double delta_z, int N)
{
    /* Computes the difference in optical depth for each cell in the grid
     *
     * Input
     * ------
     * z_arr: the spatial grid over which the optical depth will change, double precision, size N=n+1
     * delta_z: the step size of the grid, double precision, single value
     * N: the size of the grid + 1: n+1=N, integer, single value
     *
     * Output
     * ------
     * delta_tau: an array containing the difference in optical depth at each grid for each cell, double precision,
     *            size N */

    int k;    // Cell index
    double x; // Temporary double precision variable for each loop

    for (k=0;k<N;k++)
    {
        x = z_arr[k] + (delta_z/2.);
        // The definition of the difference in optical depth for this problem is: sqrt(3.)*alpha[i]*delta_z
        // This can be found in Ch.4, section 4.4.5, equation 4.55
        // The extinction coefficient (alpha[i]) can be found in Ch.4: section 4.4.2, equation 4.32
        // alpha[i] = (10.^(5. - 6.*x[i])), where x[i] is our current cell (not grid point)
        delta_tau[k] = sqrt(3.)*delta_z*pow(10.,(5.-6.*x));
    }

    // delta_tau has been found!
    return;
}

// Define the third order quadratic interpolation coefficients function
void quad_int(double *up, double *pp, double *dp, double *um, double *pm, double *dm, double *delta_tau, int N)
{
    /* Computes the three interpolation coefficients for the bottom up and top down rays
     *
     * Input
     * ------
     * delta_tau: the difference in optical depth for each grid cell, double precision, size n+1=N
     * N: the size of the grid+1: n+1=N, integer, single value
     *
     * Output
     * ------
     * up, pp, dp: the three interpolation coefficients for the ray integrating from the bottom of the grid up, size n
     * um, pm, dm: the three interpolation coefficients for the ray integrating from the top of the grid down, size n */

    double *e0, *e1, *e2; // interpolation coefficients
    int k;                // Cell index

    e0 = (double *) calloc(N, sizeof(double)); // Array used to calculate interpolation coefficients
    e1 = (double *) calloc(N, sizeof(double)); // Array used to calculate interpolation coefficients
    e2 = (double *) calloc(N, sizeof(double)); // Array used to calculate interpolation coefficients

    // Loop through and set up the optical depth dependence for each height
    // Ch.3 section 3.8.4, equations 3.43-3.45
    for (k=0;k<N;k++)
    {
        e0[k] = 1. - exp(-delta_tau[k]);
        e1[k] = delta_tau[k] - e0[k];
        e2[k] = pow(delta_tau[k],2.) - (2.*e1[k]);
    }

    // Now compute the actual interpolation coefficients for each ray using the constants above
    // The equations for these can be found in Ch.3, section 3.8.4, equations 3.40-3.42
    // Note that my notation has the following equivalences: u=u, v=p, w=d
    for (k=0;k<N;k++)
    {
        // The coefficients for the bottom up array (I+)
        up[k] = e0[k] + (e2[k] - (2.*delta_tau[k] + delta_tau[k+1])*e1[k])/(delta_tau[k]*(delta_tau[k]+delta_tau[k+1]));
        pp[k] = (((delta_tau[k]+delta_tau[k+1])*e1[k])-e2[k])/(delta_tau[k]+delta_tau[k+1]);
        dp[k] = (e2[k] - delta_tau[k]*e1[k])/(delta_tau[k+1]*(delta_tau[k]+delta_tau[k+1]));

        // The coefficients for the top down array (I-)
        um[k] = e0[k] + (e2[k] - (2.*delta_tau[k] + delta_tau[k-1])*e1[k])/(delta_tau[k]*(delta_tau[k]+delta_tau[k-1]));
        pm[k] = (((delta_tau[k]+delta_tau[k-1])*e1[k])-e2[k])/(delta_tau[k]+delta_tau[k-1]);
        dm[k] = (e2[k] - delta_tau[k]*e1[k])/(delta_tau[k-1]*(delta_tau[k]+delta_tau[k-1]));
    }

    // Free the memory
    free(e0);
    free(e1);
    free(e2);

    // Return the function
    return;
}

// Set the Lambda_star matrix and M_star matrix tri-diagonal components
void matrix_build(double **lambda_s, double *m_a, double *m_b, double *m_c,
        double *up, double *pp, double *dp,double *um, double *pm, double *dm, double eps, int n)
{
    /* This will compute the lambda_star and M lower and upper diagonals and the diagonal.
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

// Setup the source function interpolation, this includes the quadrature limiter that does not allow sources to
// become negative (aka non-physical)
void S_interp(double S_int, double *alpha_arr, double *S_arr, double *z_arr, double *up, double *pp, double *dp,
        double *um, double *pm, double *dm, double delta_z, int x, bool plus)
{
    /* Source term interpolation function that includes a maximum quadrature limiter that does not allow any source term
     * to have a negative value.
     * This value is going to be different for the two rays going in different directions, and this is taken into
     * account here.
     *
     * Input
     * -----
     * alpha_arr: Array containing the value of the emissivity at each grid point, double precision, size n
     * S_arr: An array containing the values of the source function at every point in the grid, double precision, size n
     * z_arr: The grid, double precision, size n
     * delta_z: The step size between grid points, double precision, single value
     * x: The location within the grid that needs to be interpolated, x in z[x], integer, single value.
     * plus: This indicates which array is being dealt with; I+ --> plus = True or I- --> plus = False
     *
     * Output
     * ------
     * S_int: The value of the interpolated source function for either the I+ or I- array at grid location x,
     *        from z[x], double precision, single value */

    double Q_max;  // This is the upper value of the limiter
    double Q_curr; // This is the value of the interpolation using the coefficients ane S_arr

    /* The quadrature limiter for third order integration is discussed in Ch.3, section 3.8.4, equation 3.39
     * This equation is for the I+ vector: Q_curr = up*S[x-1] + pp*S[x] + dp*S[x+1]
     * For the maximum upper limit the equation is found in Ch.3, section 3.8.4, equation 3.37
     * Q_max = 0.5*(alpha[x] + alpha[x-1])*delta_z
     * The final value of the source term interpolation is determined by equation 3.46, section 3.8.4, Ch.3
     * S_int = max(min(Q_curr[x],Q_max[x]),0.) */

    Q_max = 0.5*(alpha_arr[x] - alpha_arr[x-1])*delta_z;

    if (plus=='True')
    {
        Q_curr = up[x]*S_arr[x-1] + pp[x]*S_arr[x] + dp[x]*S_arr[x+1];
        S_int = fmax(fmin(Q_curr,Q_max),0.);
    }
    else
    {
        Q_curr = um[x]*S_arr[x+1] + pm[x]*S_arr[x] + dm[x]*S_arr[x-1];
        S_int = fmax(fmin(Q_curr,Q_max),0.);
    }

    // S_int has been found for either I+ or I- so return the function
    return;

}

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

void mat_mult(double **mat_a, double **vec_b, double **mat_c, int N)
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

    int i, j, k; // The necessary counters


    // Initialize mat_c to have 0. in all elements
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; j++)
        {
            mat_c[i][j] = 0;
        }
    }

    // Do the multiplication and find the resulting matrix
    for (i = 0; i < N; ++i)
    {
        for (j = 0; j < N; ++j)
        {
            for (k = 0; k < N; ++k)
            {
                // For each element mat_c[i][j] you must add the n_a multiplications
                mat_c[i][j] += mat_a[i][k] * vec_b[k][0];
            }
        }
    }


    // mat_c has been found so we are done
    return;
}

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