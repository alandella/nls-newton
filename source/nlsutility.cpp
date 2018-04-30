/*
	Andrea Landella

	program file nls-utility.cpp
*/
#include<stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>

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
			m = (double**)malloc((DIM - 1)*sizeof(double *));
			for (i = 0; i < DIM - 1; i++){
				m[i] = (double*)malloc((DIM - 1)*sizeof(double));
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

//calculate the hollowed matrix
double** calc_holmat(double** A, int row, int col, int DIM){

	int jj, kk;

	double** C = NULL;
	//
	C = (double**)malloc((DIM - 1)*sizeof(double*));

	for (int j = 0; j < DIM - 1; j++){
		C[j] = (double*)malloc((DIM - 1)*sizeof(double));

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
	return C;
}

//calculate the cofactor matrix
double** calc_cofmat(double** A, int DIM){

	double** C = NULL;
	double** K;

	C = (double**)malloc(DIM*sizeof(double*));
	//
	for (int j = 0; j < DIM; j++){
		C[j] = (double*)malloc(DIM*sizeof(double));
	}

	for (int j = 0; j < DIM; j++){
		for (int k = 0; k < DIM; k++){
			//we write the M matrix
			K = calc_holmat(A, j, k, DIM);
			//we calculate the element
			C[j][k] = pow(-1., j + k)*calc_det(K, DIM - 1);
		}
	}
	return C;
}

//calculate transpose
double** calc_transp(double** A, int DIM){

	double** B = (double**)malloc(DIM*sizeof(double*));
	//
	for (int j = 0; j < DIM; j++){
		B[j] = (double*)malloc(DIM*sizeof(double));
	}

	for (int i = 0; i < DIM; i++){
		for (int j = 0; j < DIM; j++){
			B[i][j] = A[j][i];
		}
	}
	return B;
}

//calculate the matrix inverse
double** calc_inv(double** A, int DIM){

	//the inverse is 1/det * (cof)_transp
	double A_det;
	double** A_cof;
	double** A_tsp;
	double** A_inv = (double**)malloc(DIM*sizeof(double*));
	//
	for (int j = 0; j < DIM; j++){
		A_inv[j] = (double*)malloc(DIM*sizeof(double));
	}

	//we have the det
	A_det = calc_det(A, DIM);

	if (A_det == 0){
		printf("Warning! Matrix is singular!\n");
	}

	//we have the cofactor matrix
	A_cof = calc_cofmat(A, DIM);

	//we transpose it
	A_tsp = calc_transp(A_cof, DIM);

	//we evaluate the inverse:
	for (int i = 0; i < DIM; i++){
		for (int j = 0; j < DIM; j++){
			A_inv[i][j] = 1 / A_det * A_tsp[i][j];
		}
	}

	for (int j = 0; j < DIM; j++){
		free(A_cof[j]);
		free(A_tsp[j]);
	}
	free(A_cof);
	free(A_tsp);

	return A_inv;
}

//we calculate the matrix*vector delta_x
double* calc_matvec(double** A, double* b, int DIM){

	double* f = (double*)malloc(DIM*sizeof(double));
	double ff;

	for (int j = 0; j < DIM; j++){
		ff = 0.;
		for (int k = 0; k < DIM; k++){
			ff += A[j][k] * b[k];
		}
		f[j] = ff;
	}
	return f;
}