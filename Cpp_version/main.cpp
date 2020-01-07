#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <vector>


/* This script provides a correctness check for the tridiagonal matrix solver for the
 * equation: Mx = y. Where M is a tridiagonal matrix and x is being solved for. This
 * solver is then used in the main RT code. The two cases that need to be checked are
 * the diagonally dominant and non-diagonally dominant cases. */

/* The following is the layout of the M matrix
 * | b[0] c[0]  0   ...  0  |
 * | a[1] b[1] c[1]  0  ... |
 * |  0   a[2] b[2] c[2] ...|
 * |  .    .    .    .    . |
 * |  .    .    .    .    . |
 * |  .    .    .    .    . |
 * |...a[n-2] b[n-2] c[n-2] |
 * |...  ...  a[n-1] b[n-1] |
 * Note: a[0] = 0 and c[n-1] = 0 */

// Define the function that actually solve Mx = y
void tri_solver(double *a, double *b, double *c, double *y, double *X, int n)
{
    /* Compute the array x in the equation Mx = y given the lower diagonal a, diagonal
     * b, and upper diagonal c of matrix M. */

    int i;
    double *c_prime, *y_prime;

    // Allocate the prime arrays
    c_prime = (double *) calloc(n,sizeof(double));
    y_prime = (double *) calloc(n,sizeof(double));

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

    for(i=0;i<n;i++)
    {
        if(i==0)
        {
            c_prime[i] = c[i]/b[i];
            y_prime[i] = y[i]/b[i];
        }
        else if(i>0 && i<n-1)
        {
            c_prime[i] = c[i]/(b[i]-a[i]*c_prime[i-1]);
            y_prime[i] = (y[i]-a[i]*y_prime[i-1])/(b[i]-a[i]*c_prime[i-1]);
        }
        else if(i==n-1)
        {
            y_prime[i] = (y[i]-a[i]*y_prime[i-1])/(b[i]-a[i]*c_prime[i-1]);
            X[i] = y_prime[i];
        }
    }
    for(i=n-2;i>=0;i--)
    {
        X[i] = y_prime[i] - c_prime[i]*X[i+1];
    }

    // Free the memory
    free(c_prime);
    free(y_prime);

    // X has been found so exit the function
    return;
}

int main(int argc, char **argv)
{
    // Define and allocate the variables
    int j;    // Cell index
    int n=10; // For an nxn matrix we need n+1 elements
        // use n=31 for non-diagonally dominant M
        // use n=10 for diagonally dominant M
    bool flag = false; // True = non-diagonally dominant M, False = diagonally dominant M

    double *aa = new double[n];
    double *bb = new double[n];
    double *cc = new double[n];
    double *yy = new double[n];
    double *xx = new double[n];

    // Non-diagonally dominant example
    //---------------------------------
    if(flag==true)
    {
        // Initialize
        for(j=0;j<n-1;j++)
        {
            // M matrix lower diagonal, diagonal, and upper diagonal components
            aa[j] = 1.0;
            bb[j] = 2.0;
            cc[j] = 1.0;

            if(j<15)
            {
                yy[j] = (float) (j+1);
            }
            else
            {
                yy[j] = (float) (n-j);
            }
        }
        aa[0] = 0.;
        aa[n-1] = 1.;
        cc[n-1] = 0.;
        bb[n-1] = 2.;

        // Calculate solution
        tri_solver(aa,bb,cc,yy,xx,n);

        // Print calculated solution
        for(j=0;j<n;j++)
        {
            printf("j %d, x %2.1f\n",j,xx[j]);
        }
    }

    // Diagonally dominant example
    //---------------------------------
    else
    {
        // Initialize
        for(j=0;j<n-1;j++)
        {
            // M matrix lower diagonal, diagonal, and upper diagonal components
            aa[j] = 0.5;
            bb[j] = 2.0;
            cc[j] = 0.5;
            yy[j] = 4.0;
        }
        aa[0] = 0.;
        aa[n-1] = 0.5;
        cc[n-1] = 0.;
        bb[n-1] = 2.;
        yy[n-1] = 4.;

        // Calculate solution
        tri_solver(aa,bb,cc,yy,xx,n);

        // Print calculated solution
        printf("Calculated solution \n");
        for(j=0;j<n;j++)
        {
            printf("j %d, x %2.5f\n",j,xx[j]);
        }

        // Calculate analytical solution
        /* This is done by solving x = M^(-1) y
         * 1) actually build matrix M
         * 2) find the inverse
         * 3) solve for x */

        // allocate the space for M
        std::vector<std::vector<double>> M_check(n,std::vector<double>(n));

        // Assign the values to M
        // Set all elements equal to 0
        // Loop through and set lower diagonal, diagonal, and upper diagonal
        for(j=0;j<n;j++)
        {
            if(j==0)
            {
                M_check[j][j] = bb[j];
                M_check[j][j+1] = cc[j];
            }
            else if(j==n-1)
            {
                M_check[j][j] = bb[j];
                M_check[j][j-1] = aa[j];
            }
            else
            {
                M_check[j][j-1] = aa[j];
                M_check[j][j] = bb[j];
                M_check[j][j+1] = cc[j];

            }
        }

        // Need to actually do calculation but this works for now.
        printf("Actual answer \n");
        printf("[1.69059656 1.23761375 1.35894843 1.32659252 1.3346815  1.3346815 \n"
               " 1.32659252 1.35894843 1.23761375 1.69059656]");

    }

    return 0;
}