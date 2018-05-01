/*
	Andrea Landella

	program file nlsnewton

	We want to solve the 1-D and DIM-D square NLS:

	0 = F(x)

	with a given x0, using the multidimensional Newton Method.

	In order to do that, we need two functions: one for the system, F, and the
	other for the solver, which takes as input the system function's handle.
*/

#include<stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include"nlssolver.h"

// test function for scalar system
double testfun(double x);

// test function for vector system
void testfun_vec(double* x, int DIM, double Fx[]);

int main(){

	/******* SCALAR SYSTEM *******/

	//we define absolute tolerances for function and variable 
	double xtol = 1.e-7;
	double Ftol = 1.e-6;

	printf("NLE-Newton Iteration\n\n");

	//initial guess
	double x0 = 2.1;

	//scalar (NLE) solver call
	double ris = nlsnewton(&testfun, x0, xtol, Ftol);

	//we print the result and the fun value at the result
	printf("x_opt = %e\n", ris);
	//
	printf("norm(F(x_opt)) = %e\n\n", testfun(ris));

	///******* VECTOR SYSTEM *******/

	//we define the system dimension
	const int DIM = 3;

	//we define absolute tolerances for function and variable norm
	double xtol_vec = 1.e-8;
	double Ftol_vec = 1.e-7;

	printf("NLS-Newton (DIM = %d) Iteration\n\n", DIM);

	//initial guess
	double x0_vec[DIM] = { 1.0, 1.0, 1.0 };

	//vectorial (NLS) solver call
	double ris_vec[DIM], Fris_vec[DIM];
	//
	nlsnewton_vec(&testfun_vec, x0_vec, xtol_vec, Ftol_vec, DIM, ris_vec);

	//we print the result and the fun norm at the result
	printf("x_opt = %e , %e , %e\n", ris_vec[0], ris_vec[1], ris_vec[2]);
	//
	testfun_vec(ris_vec, DIM, Fris_vec);
	//
	printf("norm(F(x_opt)) = %e\n\n", norm2_vec(Fris_vec, DIM));
	system("pause");

	return 0;
}

// test function for scalar system
double testfun(double x){

	double Fx;

	Fx = 5*x + log(fabs(x)) - 10000.;

	return Fx;
}

// test function for vector system
void testfun_vec(double* x, int DIM, double Fx[]){

	//nonlinear system, from 
	//http://fourier.eng.hmc.edu/e176/lectures/NM/node21.html 

	Fx[0] = 3.*x[0] - cos(x[1]*x[2]) - 3./2.;
	Fx[1] = 4.*pow(x[0], 2) - 625.*pow(x[1], 2) + 2.*x[2] - 1.;
	Fx[2] = 20.*x[2] + exp(-x[0] * x[1]) + 9.;
}