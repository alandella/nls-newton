/*
	Andrea Landella

	program file solver
*/

#include<stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include"nlssolver.h"


/***** UTILITIES *****/

//calculate 2-norm of vector
double norm2_vec(double* v, int DIM){

	double f = 0.;

	for (int j = 0; j < DIM; j++){
		f += pow(v[j], 2);
	}
	return pow(f, 0.5);
}

//calculate 2-norm of vector difference (v - w) or (w - v), result is the same
double norm2_vec_diff(double* v, double* w, int DIM){

	double f = 0.;

	for (int j = 0; j < DIM; j++){
		f += pow(v[j] - w[j], 2);
	}
	return pow(f, 0.5);
}

//calculate determinant
double calc_det(double **A, int DIM) {

	int i, j, j1, j2;
	double det = 0;
	double **m = NULL;

	if (DIM < 1) { /* Error */
		printf("Error in calc_det: DIM < 1\n");
	}
	else if (DIM == 1) { /* Shouldn't get used */
		det = A[0][0];
	}
	else if (DIM == 2) {
		det = A[0][0] * A[1][1] - A[1][0] * A[0][1];
	}
	else {
		det = 0;
		for (j1 = 0; j1<DIM; j1++) {
			m = (double**)malloc((DIM - 1)*sizeof(double *)); if (!m) return NULL;
			for (i = 0; i < DIM - 1; i++){
				m[i] = (double*)malloc((DIM - 1)*sizeof(double)); if (!m[i]) return NULL;
			}
			for (i = 1; i<DIM; i++) {
				j2 = 0;
				for (j = 0; j<DIM; j++) {
					if (j == j1){
						continue;
					}
					m[i - 1][j2] = A[i][j];
					j2++;
				}
			}
			det += pow(-1.0, 1.0 + j1 + 1.0) * A[0][j1] * calc_det(m, DIM - 1);
			for (i = 0; i < DIM - 1; i++){
				free(m[i]);
			}
			free(m);
		}
	}
	return(det);
}

//calculate the hollowed matrix C
void calc_holmat(double** A, int row, int col, int DIM, double** &C){

	int jj, kk;

	for (int j = 0; j < DIM - 1; j++){
		for (int k = 0; k < DIM - 1; k++){
			if (j >= row) {
				jj = j + 1;
			}
			else {
				jj = j;
			}
			if (k >= col) {
				kk = k + 1;
			}
			else {
				kk = k;
			}
			C[j][k] = A[jj][kk];
		}
	}
}

//calculate the cofactor matrix C
void calc_cofmat(double** A, int DIM, double** &C){

	double** K = (double**)malloc(DIM*sizeof(double*)); if (!K) return;
	for (int j = 0; j < DIM; j++){
		K[j] = (double*)malloc(DIM*sizeof(double)); if (!K[j]) return;
	}

	/* Cij is the det of the ij-hollowed matrix times (-1)^(i+j) */

	for (int j = 0; j < DIM; j++){
		for (int k = 0; k < DIM; k++){

			calc_holmat(A, j, k, DIM, K);

			C[j][k] = pow(-1., j + k)*calc_det(K, DIM - 1);
		}
	}
	// free memory
	for (int j = 0; j < DIM; j++){
		free(K[j]);
	}
	free(K);
}

//calculate transpose C
void calc_transp(double** A, int DIM, double** &C){

	for (int i = 0; i < DIM; i++){
		for (int j = 0; j < DIM; j++){
			C[i][j] = A[j][i];
		}
	}
}

//calculate the matrix inverse C = A^(-1)
void calc_inv(double** A, int DIM, double** &C){

	//the inverse is 1/det * (cof)_transp
	double A_det;
	//
	double** A_cof = (double**)malloc(DIM*sizeof(double*)); if (!A_cof) return;
	double** A_tsp = (double**)malloc(DIM*sizeof(double*)); if (!A_tsp) return;
	//
	for (int j = 0; j < DIM; j++){
		A_cof[j] = (double*)malloc(DIM*sizeof(double)); if (!A_cof[j]) return;
		A_tsp[j] = (double*)malloc(DIM*sizeof(double)); if (!A_tsp[j]) return;
	}

	//we evaluate the det, with singularity flag
	A_det = calc_det(A, DIM);
	//
	if (A_det == 0.){
		printf("Warning! Matrix is singular!\n");
	}

	//we evaluate the cofactor matrix
	calc_cofmat(A, DIM, A_cof);

	//we transpose it
	calc_transp(A_cof, DIM, A_tsp);

	//we evaluate the inverse matrix:
	for (int i = 0; i < DIM; i++){
		for (int j = 0; j < DIM; j++){
			C[i][j] = 1 / A_det * A_tsp[i][j];
		}
	}
	//free memory
	for (int j = 0; j < DIM; j++){
		free(A_cof[j]);
		free(A_tsp[j]);
	}
	free(A_cof);
	free(A_tsp);
}

//we calculate the matrix*vector C = A*b
void calc_matvec(double** A, double* b, int DIM, double* &C){

	double ff;

	for (int j = 0; j < DIM; j++){
		ff = 0.;
		for (int k = 0; k < DIM; k++){
			ff += A[j][k] * b[k];
		}
		C[j] = ff;
	}
}


/***** JACOBIAN *****/

//calculate numerical jacobian, internal h = 1.e-8
void calc_Jac(void fun(double* x, int DIM, double Fx[]), double* xn, int DIM, double** &Jac){

	//declare the perturbed point xnm1 now equal to xn
	double* xnm1 = (double*)malloc(DIM*sizeof(double)); if (!xnm1) return;
	//
	for (int j = 0; j < DIM; j++){
		xnm1[j] = xn[j];
	}

	//declare the function values
	double* F_xn = (double*)malloc(DIM*sizeof(double)); if (!F_xn) return;
	double* F_xnm1 = (double*)malloc(DIM*sizeof(double)); if (!F_xnm1) return;

	//evaluating the function at xn (unperturbed)
	fun(xn, DIM, F_xn);

	for (int i = 0; i < DIM; i++){		//row, its dF
		for (int j = 0; j < DIM; j++){	//col, its dx

			//perturbing and evaluating the function at xnm1 (perturbed)
			xnm1[j] += 1.e-8;
			//
			fun(xnm1, DIM, F_xnm1);

			//evaluate the jacobian component dF[i]/dx[j]
			Jac[i][j] = (F_xn[i] - F_xnm1[i]) / (xn[j] - xnm1[j]);

			//reset perturbation
			xnm1[j] = xn[j];
		}
	}
	free(xnm1);
	free(F_xn);
	free(F_xnm1);
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


/***** SOLVER *****/

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
	while (eval_DxFx(Dx, xtol, Fx, Ftol)){

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

		//we calculate the required differences (norm of step and fval):
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

// solver for vector system
void nlsnewton_vec(void fun(double* x, int DIM, double Fx[]), double* x0, double xtol, double Ftol, int DIM, double xN[]){

	//we have the two iteration values and the function values
	double* x1 = (double*)malloc(DIM*sizeof(double)); if (!x1) return;
	//
	double* Fx0 = (double*)malloc(DIM*sizeof(double)); if (!Fx0) return;
	double* Fx1 = (double*)malloc(DIM*sizeof(double)); if (!Fx1) return;

	//we declare the jacobian and its inverse
	double** Jx0 = (double**)malloc(DIM*sizeof(double*)); if (!Jx0) return;
	double** Jx0_inv = (double**)malloc(DIM*sizeof(double*)); if (!Jx0_inv) return;
	//
	for (int j = 0; j < DIM; j++){
		Jx0[j] = (double*)malloc(DIM*sizeof(double)); if (!Jx0[j]) return;
		Jx0_inv[j] = (double*)malloc(DIM*sizeof(double)); if (!Jx0_inv[j]) return;
	}

	//we declare the delta_x
	double* delta_x = (double*)malloc(DIM*sizeof(double)); if (!delta_x) return;

	//we set x0_new and evaluate the function
	for (int j = 0; j < DIM; j++){
		xN[j] = x0[j];
	}
	//
	fun(xN, DIM, Fx0);

	//we evaluate the norms of function and variables, first Dx is dummy
	double Dx = 10.;
	double Fx = norm2_vec(Fx0, DIM);
	//
	if (Fx < Ftol){
		printf("Initial point is a solution!\n");
	}

	//we have the iteration number and relaxation parameter
	int N_iter = 0;
	//
	double omega = 1.0;

	//entering the main calculation cycle
	while (eval_DxFx(Dx, xtol, Fx, Ftol)){

		N_iter++;

		//we evaluate the function at x0
		fun(xN, DIM, Fx0);

		//we evaluate the jacobian and its inverse at x0
		calc_Jac(fun, xN, DIM, Jx0);
		//
		calc_inv(Jx0, DIM, Jx0_inv);

		//we calculate the matrix*vector delta_x
		calc_matvec(Jx0_inv, Fx0, DIM, delta_x);

		//we obtain the new value, with relaxation
		for (int j = 0; j < DIM; j++){
			x1[j] = xN[j] - omega*delta_x[j];
		}

		//we evaluate the function at x1
		fun(x1, DIM, Fx1);

		//we calculate the required differences (norm of step and fval):
		Dx = norm2_vec_diff(xN, x1, DIM);
		//
		Fx = norm2_vec(Fx1, DIM);

		printf("IT%d\tDx = %e\n\tFx = %e\n\n", N_iter, Dx, Fx);
		//getchar();

		//we advance the cycle
		for (int j = 0; j < DIM; j++){
			xN[j] = x1[j];
		}

		//we have the MaxIter stop
		if (N_iter >= 1e5){
			printf("Max Iterations Exceeded\n");
			break;
		}
	}
	// free memory
	free(x1);
	free(Fx0);
	free(Fx1);
	for (int j = 0; j < DIM; j++){
		free(Jx0[j]);
		free(Jx0_inv[j]);
	}
	free(Jx0);
	free(Jx0_inv);
	free(delta_x);
}

// "true" while condition evaluation
int eval_DxFx(double Dx, double xtol, double Fx, double Ftol){

	//return 1 -> while(1 != 0) == true  , continue while
	//return 0 -> while(0 != 0) == false , stop     while

	//we define the relaxation parameter:
	double rel_par = 1.e-3;

	//we evaluate the maximum of the two deltas (Dj - jtol)
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