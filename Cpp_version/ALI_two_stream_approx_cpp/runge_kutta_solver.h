#ifndef RUNGE_KUTTA_SOLVER_H
#define RUNGE_KUTTA_SOLVER_H
void RK4_solver(double x_curr, double *f_curr, double h, double *dens, int n, double *f_update);
#endif //RUNGE_KUTTA_SOLVER_H 
