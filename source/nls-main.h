/*
	Andrea Landella

	file header per main
*/

double testfun(double x);
double* testfun_vec(double* x, int DIM);

double* nlsnewton_vec(double* fun(double* x, int DIM), double* x0, double xtol, double Ftol, int DIM);