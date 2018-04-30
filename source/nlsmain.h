/*
	Andrea Landella

	file header per main
*/

// test function for scalar system
double testfun(double x);
// test function for vector system
double* testfun_vec(double* x, int DIM);

// solver for scalar system
double nlsnewton(double fun(double x), double x0, double xtol, double Ftol);
// solver for vector system
double* nlsnewton_vec(double* fun(double* x, int DIM), double* x0, double xtol, double Ftol, int DIM);

// "true" while condition evaluation
int eval_DxFx(double Dx, double xtol, double Fx, double Ftol);