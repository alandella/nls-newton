/*
	Andrea Landella

	file header per jacobian
*/

//calculate numerical jacobian, internal h = 1.e-8
double** calc_Jac(double* fun(double* x, int DIM), double* xn, int DIM);

//calculate numerical derivative, internal h = 1.e-8
double calc_dFdx(double fun(double x), double xn);