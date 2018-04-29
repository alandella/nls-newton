/*
	Andrea Landella

	program file nls-jacobian.cpp
*/
#include<stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>

//we begin, by evaluating the jacobian at x0
double** calc_Jac(double* fun(double* x, int DIM), double* xn, int DIM){

	double** Jac = NULL;

	Jac = (double**)malloc(DIM*sizeof(double*));

	for (int j = 0; j < DIM; j++){
		Jac[j] = (double*)malloc(DIM*sizeof(double));
	}

	double* xnm1;

	xnm1 = (double*)malloc(DIM*sizeof(double));

	for (int j = 0; j < DIM; j++){
		xnm1[j] = xn[j];
	}

	//evaluate the function values
	double* F_xn; 
	double* F_xnm1; 	

	F_xn = fun(xn, DIM);

	//evaluating the jacobian components
	for (int i = 0; i < DIM; i++){		//row, its F
		for (int j = 0; j < DIM; j++){	//col, its x
			
			xnm1[j] += 1.e-8;

			F_xnm1 = fun(xnm1, DIM);

			Jac[i][j] = (F_xn[i] - F_xnm1[i]) / (xn[j] - xnm1[j]);

			xnm1[j] = xn[j];
		}
	}

	free(xnm1);
	
	return Jac;
}


