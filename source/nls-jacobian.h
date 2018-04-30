/*
	Andrea Landella

	file header per jacobian
*/
double** calc_Jac(double* fun(double* x, int DIM), double* xn, int DIM);

double calc_dFdx(double fun(double x), double xn);