/*
	Andrea Landella

	program file nls-newton.cpp

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

#include"nlsmain.h"
#include"nlsutility.h"
#include"nlsjacobian.h"

int main(){

	/******* SCALAR SYSTEM *******/

	double xtol = 1.e-7;
	double Ftol = 1.e-6;

	printf("NLE-Newton Iteration\n\n");

	double x0 = 2.1;

	double ris = nlsnewton(&testfun, x0, xtol, Ftol);

	printf("x_opt = %e\n", ris);

	double Fx = testfun(ris);

	printf("norm(F(x_opt)) = %e\n\n", Fx);

	/******* VECTOR SYSTEM *******/

	const int DIM = 3;

	double* x0_vec = (double*)malloc(DIM*sizeof(double));
	//
	x0_vec[0] = 1.0;
	x0_vec[1] = 1.0;
	x0_vec[2] = 1.0;

	printf("NLS-Newton Iteration\n\n");

	double* ris_vec = nlsnewton_vec(&testfun_vec, x0_vec, xtol, Ftol, DIM);

	printf("x_opt = %e , %e, %e\n", ris_vec[0], ris_vec[1], ris_vec[2]);

	double* Fris_vec = testfun_vec(ris_vec, DIM);

	double Fx_vec = norm2_vec(Fris_vec, DIM);

	printf("norm(F(x_opt)) = %e\n\n", Fx_vec);

	system("pause");
	return 0;
}

// test function for scalar system
double testfun(double x){

	double Fx;

	Fx = 5*x + log(fabs(x)) - 10000.;

	return Fx;
}

// solver for scalar system
double nlsnewton(double fun(double x), double x0, double xtol, double Ftol){
	
	//we have the two iteration values and the function values
	double x0_new, Fx0;
	double x1, Fx1; 

	//we declare the jacobian, or the derivative dF/dx
	double Jx0;

	//we declare the delta_x
	double delta_x;

	//we set x0_new and evaluate the function
	x0_new = x0;
	//
	Fx0 = fun(x0_new);

	//we evaluate the norms of function and variables, first Dx is dummy
	double Dx = 10.;
	double Fx = fabs(Fx0);
	//
	if (Fx < Ftol){
		printf("Initial point is a solution!\n");
		return x0_new;
	}

	//we have the iteration number and relax parameter
	int N_iter = 0;
	//
	double omega = 1.0;

	//entering the main calculation cycle
	while ( eval_DxFx(Dx, xtol, Fx, Ftol) ){

		N_iter++;

		//we evaluate the function at x0
		Fx0 = fun(x0_new);

		//we evaluate the jacobian at x0
		Jx0 = calc_dFdx(fun, x0_new);

		//we calculate the matrix*vector term, delta_x
		delta_x = Fx0 / Jx0;

		//we obtain the new value, with relaxation
		x1 = x0_new - omega*delta_x;

		//we evaluate the function at x1
		Fx1 = fun(x1);

		//we calculate the required differences:
		Dx = fabs(x0_new - x1);
		//
		Fx = fabs(Fx1);

		printf("IT%d\tDx = %e\n\tFx = %e\n\n", N_iter, Dx, Fx);
		//getchar();

		//we advance the cycle
		x0_new = x1;

		//we have the MaxIter stop
		if (N_iter >= 1e5){
			printf("Max Iterations Exceeded\n");
			break;
		}
	}
	return x0_new;
}

// test function for vector system
double* testfun_vec(double* x, int DIM){

	double* Fx;

	Fx = (double*)malloc(DIM*sizeof(double));

	//nonlinear system, from 
	//http://fourier.eng.hmc.edu/e176/lectures/NM/node21.html 

	Fx[0] = 3.*x[0] - cos(x[1]*x[2]) - 3./2.;
	Fx[1] = 4.*pow(x[0], 2) - 625.*pow(x[1], 2) + 2.*x[2] - 1.;
	Fx[2] = 20.*x[2] + exp(-x[0] * x[1]) + 9.;

	return Fx;
}

// solver for vector system
double* nlsnewton_vec(double* fun(double* x, int DIM), double* x0, double xtol, double Ftol, int DIM){

	//we have the two iteration values and the function values
	double* x0_new = (double*) malloc(DIM*sizeof(double));
	double* x1 = (double*) malloc(DIM*sizeof(double));
	//
	double* Fx0;
	double* Fx1;

	//we declare the jacobian and its inverse
	double** Jx0;
	double** Jx0_inv;

	//we declare the delta_x
	double* delta_x;

	//we set x0_new and evaluate the function
	for (int j = 0; j < DIM; j++){
		x0_new[j] = x0[j];
	}
	//
	Fx0 = fun(x0_new, DIM);

	//we evaluate the norms of function and variables, first Dx is dummy
	double Dx = 10.;
	double Fx = norm2_vec(Fx0, DIM);
	//
	if (Fx < Ftol){
		printf("Initial point is a solution!\n");
		return x0_new;
	}

	//we have the iteration number and relax parameter
	int N_iter = 0;
	//
	double omega = 1.0;

	//entering the main calculation cycle
	while ( eval_DxFx(Dx, xtol, Fx, Ftol) ){
		
		N_iter++;

		//we evaluate the function at x0
		Fx0 = fun(x0_new, DIM);

		//we evaluate the jacobian and its inverse at x0
		Jx0 = calc_Jac(fun, x0_new, DIM);
		//
		Jx0_inv = calc_inv(Jx0, DIM);

		//we calculate the matrix*vector delta_x
		delta_x = calc_matvec(Jx0_inv, Fx0, DIM);

		//we obtain the new value, with relaxation
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

		//we have the MaxIter stop
		if (N_iter >= 1e5){
			printf("Max Iterations Exceeded\n");
			break;
		}
	}
	free(x1);
	return x0_new;
}

// "true" while condition evaluation
int eval_DxFx(double Dx, double xtol, double Fx, double Ftol){

	//return 1 -> 1 != 0 true  , continue while
	//return 0 -> 0 != 0 false , stop     while

	//define the relaxation parameter:
	double rel_par = 1.e-3;

	//evaluate the maximum of the two deltas (Dj - jtol)
	if (fabs(Dx - xtol) > fabs(Fx - Ftol)){		
		if (fabs(Dx - xtol) > (1 + rel_par)*xtol){
			return 1;
		}
		else {
			return 0;
		}
	}
	else {
		if (fabs(Fx - Ftol) > (1 + rel_par)*Ftol){
			return 1;
		}
		else{
			return 0;
		}
	}
}