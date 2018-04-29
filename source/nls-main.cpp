/*
	Andrea Landella

	program file nls-newton.cpp

	We want to solve the 1-D and DIM-D NLS:

	0 = F(x)

	with a given x0, using the multidimensional Newton Method.

	In order to do that, we need two functions: one for the system, F, and the
	other for the solver, which takes as input the system function's handle.
*/

#include<stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>

#include"nls-main.h"
#include"nls-utility.h"
#include"nls-jacobian.h"

int main(){

	const int DIM = 2;

	double* ris;

	double* x0 = (double*)malloc(DIM*sizeof(double));

	x0[0] = 2.5;
	x0[1] = 5.0;

	double xtol = 1.e-9;
	double Ftol = 1.e-8;

	printf("NLS-Newton Iteration\n\n");

	ris = nlsnewton_vec(&testfun_vec, x0, xtol, Ftol, DIM);

	printf("x_opt = %e , %e\n", ris[0], ris[1]);

	double* Fris = testfun_vec(ris,DIM);

	double Fx = norm2_vec(Fris, DIM);

	printf("norm(F(x_opt)) = %e\n\n", Fx);

	system("pause");

	return 0;
}





// test function for scalar system
double testfun(double x){

	double Fx;

	Fx = pow(x, 2) - sin(x + 1);

	return Fx;
}

// solver for scalar system
/*
double nlsnewton(double testfun(double x), double x0, double xtol, double Ftol){
	...
}
*/

// test function for vector system
double* testfun_vec(double* x, int DIM){

	double* Fx;

	Fx = (double*)malloc(DIM*sizeof(double));

	//nonlinear system
	Fx[0] = -log(fabs(x[0])) + x[1];
	Fx[1] = x[0] + x[1] - 1.;

	return Fx;
}

// solver for vector system, x0 must not be the zero vector
double* nlsnewton_vec(double* fun(double* x, int DIM), double* x0, double xtol, double Ftol, int DIM){

	// we have the dummy xD point, before iteration:
	double* x0_new = (double*) malloc(DIM*sizeof(double));
	double* x1 = (double*) malloc(DIM*sizeof(double));
	//
	double* F_zero = (double*)calloc(DIM , sizeof(double));

	// we declare the function at x0 
	double* Fx0;
	double* Fx1;

	// we declare the jacobian
	double** Jx0;
	double** Jx0_inv;

	// we declare the main delta_x
	double* delta_x;

	Fx0 = fun(x0_new, DIM);

	// we evaluate the norms of function and variables
	double Dx = 10.;

	double Fx = norm2_vec(Fx0, DIM);
	
	// the first x0 is the one supplied
	for (int j = 0; j < DIM; j++){
		x0_new[j] = x0[j];
	}

	int N_iter = 0;

	double omega = 1.0;

	// entering the main calculation cycle
	//while (fmax(fabs(Dx - xtol), fabs(Fx - Ftol)) >= fmax(xtol,Ftol)*1.1){
	//while (fmax(fabs(Dx - xtol), fabs(Fx - Ftol)) >= 1.e-7){
	while ((Dx >= xtol) && (Fx >= Ftol)){
		
		N_iter++;

		//we evaluate the function at x0
		Fx0 = fun(x0_new, DIM);

		//we begin, by evaluating the jacobian at x0
		Jx0 = calc_Jac(fun, x0_new, DIM);

		//we invert the jacobian:
		Jx0_inv = calc_inv(Jx0, DIM);

		//we calculate the matrix*vector delta_x
		delta_x = calc_matvec(Jx0_inv, Fx0, DIM);

		//we obtain the new value
		for (int j = 0; j < DIM; j++){
			x1[j] = x0_new[j] - omega*delta_x[j];
		}

		//we evaluate the function at x1
		Fx1 = fun(x1, DIM);

		//we calculate the required differences:
		Dx = norm2_vec_diff(x0_new, x1, DIM);
		//
		Fx = norm2_vec(Fx1, DIM);

		printf("IT%d\tDx = %e\n\tFx = %e\n\n", N_iter, Dx, Fx);
		//getchar();

		//we advance the cycle
		for (int j = 0; j < DIM; j++){
			x0_new[j] = x1[j];
		}

		if (N_iter > 1e4){
			printf("niter exceeded\n");
			break;
		}

	}

	free(x1);
	free(F_zero);

	return x0_new;

}