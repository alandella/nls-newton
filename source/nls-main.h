/*
	Andrea Landella

	file header per main
*/

double testfun(double x);
double* testfun_vec(double* x, int DIM);

double nlsnewton(double fun(double x), double x0, double xtol, double Ftol);
double* nlsnewton_vec(double* fun(double* x, int DIM), double* x0, double xtol, double Ftol, int DIM);

int eval_DxFx(double Dx, double xtol, double Fx, double Ftol);