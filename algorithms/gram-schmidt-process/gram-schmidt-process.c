#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double norm(double *vector, int n){
	double squared_sum = 0.0;
	for (int i=0; i<n; i++){
		squared_sum += vector[i] * vector[i];
	}
	return sqrt(squared_sum);
}

double **gram_schmidt_process(double **Matrix, int n, int m){
	double** B = (double**)malloc(n * sizeof(double*));
	double vector_i[n];
	double norm_vector_i;
	double dot_product;
	double very_small_number = 1e-8;
	for (int i=0; i<n; i++){
		B[i] = (double*)malloc(m * sizeof(double));
		for (int j=0; j<m; j++){
			//copy matrix into B matrix
			B[i][j] = Matrix[i][j];
		}
	}
	for (int i=0; i<m; i++){
		// For ith column vector, do the following operation
		// Copy the ith column vector into vector_i variable.
		for(int k=0; k<n; k++){
			vector_i[k] = B[k][i];
		}
		// For subtract the components of all the orthogonal vectors found so far from column vector i
		for (int j=0; j<i; j++){
			dot_product = 0;
			// calculate the projection (dot product) of column vector i with the jth basis vector.
			for (int k=0; k<n; k++){
				dot_product += vector_i[k] * B[k][j];
			}
			// Subtract the component of vector_i along the jth vector from vector_i. B[:][j] is already normalized.
			for (int k=0; k<n; k++){
				vector_i[k] -= dot_product * B[k][j];
			}
		}
		// get the norm of the vector_i
		norm_vector_i = norm(vector_i, n);
		// Normalize the vector i if possible.
		if (norm_vector_i > very_small_number){
			for (int k=0; k<n; k++){
				B[k][i] = vector_i[k] / norm_vector_i;
			}
		}
		else{
			for (int k=0; k<n; k++)B[k][i] = 0;
		}
	}
	return B;
}

void print_matrix(double **A, int n, int m){
	printf("\n");
	for (int i=0; i<n; i++){
		for (int j=0; j<m; j++){
			printf(" %.5lf", A[i][j]);
		}
		printf("\n");
	}
}

int main(int argc, char* argv[]){
	int n, m; //row, column
	if (argc!=2){
		printf("Usage: %s INPUT_FILE\n", argv[0]);
		exit(1);
	}
	FILE *fptr = fopen(argv[1], "r");

	fscanf(fptr, "%d %d", &n, &m);

	double** Matrix = (double**)malloc(n * sizeof(double*));
	for(int i=0; i<n; i++){
		Matrix[i] = (double*)malloc(m * sizeof(double));
		for (int j=0; j<m; j++){
			fscanf(fptr, "%lf", &(Matrix[i][j]));
		}
	}

	print_matrix(Matrix, n, m);
	double **B = gram_schmidt_process(Matrix, n, m);
	print_matrix(B, n, m);
	return 0;
}
