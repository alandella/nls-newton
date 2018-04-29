/*
Andrea Landella

file header per utility
*/

double norm2_vec(double* v, int DIM);
double norm2_vec_diff(double* v, double* w, int DIM);

//calculate determinant
double calc_det(double **A, int DIM);

//calculate the cofactor matrix
double** calc_cofmat(double** A, int DIM);

//calculate transpose
double** calc_transp(double** A, int DIM);

//calculate the matrix inverse
double** calc_inv(double** A, int DIM);

//we calculate the matrix*vector delta_x
double* calc_matvec(double** A, double* b, int DIM);

//calculate the hollowed matrix
double** calc_holmat(double** A, int row, int col, int DIM);