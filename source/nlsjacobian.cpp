/*
	Andrea Landella

	program file nls-jacobian.cpp
*/
#include<stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>

//calculate numerical jacobian, internal h = 1.e-8
double** calc_Jac(double* fun(double* x, int DIM), double* xn, int DIM){

	//declare the jacobian matrix
	double** Jac = NULL;
	//
	Jac = (double**)malloc(DIM*sizeof(double*));
	//
	for (int j = 0; j < DIM; j++){
		Jac[j] = (double*)malloc(DIM*sizeof(double));
	}

	//declare the perturbed point xnm1 now equal to xn
	double* xnm1;
	//
	xnm1 = (double*)malloc(DIM*sizeof(double));
	//
	for (int j = 0; j < DIM; j++){
		xnm1[j] = xn[j];
	}

	double* F_xn; 
	double* F_xnm1; 	
	
	//evaluating the function at xn (unperturbed)
	F_xn = fun(xn, DIM);

	for (int i = 0; i < DIM; i++){		//row, its dF
		for (int j = 0; j < DIM; j++){	//col, its dx
			
			//perturbing and evaluating the function at xnm1 (perturbed)
			xnm1[j] += 1.e-8;
			//
			F_xnm1 = fun(xnm1, DIM);

			//evaluate the jacobian component dF[i]/dx[j]
			Jac[i][j] = (F_xn[i] - F_xnm1[i]) / (xn[j] - xnm1[j]);

			//reset perturbation
			xnm1[j] = xn[j];
		}
	}
	free(xnm1);
	return Jac;
}

//calculate numerical derivative, internal h = 1.e-8
double calc_dFdx(double fun(double x), double xn){

	//evaluate the perturbed point
	double xnm1 = xn + 1.e-8;

	//evaluate the function values
	double F_xn = fun(xn);
	double F_xnm1 = fun(xnm1);

	//return the incremental ratio
	return (F_xn - F_xnm1) / (xn - xnm1);

}


