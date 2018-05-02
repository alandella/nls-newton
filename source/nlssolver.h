/*
	Andrea Landella

	file header per solver
*/
#include<stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>

/***** UTILITIES *****/

//calculate 2-norm of vector
double norm2_vec(double* v, int DIM);

//calculate 2-norm of vector difference (v - w) or (w - v), result is the same
double norm2_vec_diff(double* v, double* w, int DIM);

//calculate determinant
double calc_det(double **A, int DIM);

//calculate the hollowed matrix C
void calc_holmat(double** A, int row, int col, int DIM, double** &C);

//calculate the cofactor matrix C
void calc_cofmat(double** A, int DIM, double** &C);

//calculate transpose C
void calc_transp(double** A, int DIM, double** &C);

//calculate the matrix inverse C = A^(-1)
void calc_inv(double** A, int DIM, double** &C);

//calculate the matrix*vector product C = A*b
void calc_matvec(double** A, double* b, int DIM, double* &C);


/***** JACOBIAN *****/

//calculate numerical jacobian, internal h = 1.e-8
void calc_Jac(void fun(double* x, double Fx[]), double* xn, int DIM, double** &Jac);

//calculate numerical derivative, internal h = 1.e-8
double calc_dFdx(double fun(double x), double xn);


/***** SOLVER *****/

// solver for scalar system
double nlsnewton(double fun(double x), double x0, double xtol, double Ftol);

// solver for vector system
void nlsnewton_vec(void fun(double* x, double Fx[]), double* x0, double xtol, double Ftol, int DIM, double xN[]);

// "true" while condition evaluation
int eval_DxFx(double Dx, double xtol, double Fx, double Ftol);