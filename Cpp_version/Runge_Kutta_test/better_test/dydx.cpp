#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

double dydx(double x, double f, double rho, double avg_I)
{
  return ((rho*pow(x,2.)) + (avg_I*x));
}
