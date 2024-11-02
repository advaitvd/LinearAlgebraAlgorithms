#include <stdio.h>
#include <stdlib.h>

#define VERY_SMALL_NUMBER (double)1e-8

struct Matrix{
    double **Matrix;
    int n, m;
};

void print_matrix(struct Matrix *M){
    printf("MatrixShape([%d, %d])\n", M->n, M->m);
    for(int i=0; i<(M->n); i++){
        for(int j=0; j<(M->m); j++){
            printf(" %lf", M->Matrix[i][j]);
        }
        printf("\n");
    }
}

struct Matrix* make_matrix(int n, int m){
    struct Matrix* M = malloc(sizeof(struct Matrix));
    M->n = n;
    M->m = m;
    M->Matrix = (double**)malloc(n * sizeof(double*));
    for (int i=0; i<n; i++){
        M->Matrix[i] = (double*)malloc(m * sizeof(double));
    }
    printf("Memory allocated successfully.\n");
    return M;
}

void free_memory(struct Matrix* M){
    free(M->Matrix);
    free(M);
    printf("Memory freed successfully.\n");
}

struct Matrix* matmul(struct Matrix* A, struct Matrix* B){
    if (A->m != B->n){
        printf("Matrices cannot be multiplied. (%d, %d) x (%d, %d) is prohibited.", A->n, A->m, B->n, B->m);
        exit(1);
    }
    struct Matrix* result = make_matrix(A->n, B->m);

    for (int i=0; i<(result->n); i++){
        for(int j=0; j<(result->m); j++){
            result->Matrix[i][j] = 0;
            for (int k=0; k<(A->m); k++){
                result->Matrix[i][j] += (A->Matrix[i][k]) * (B->Matrix[k][j]); 
            }
        }
    }
    return result;
}

struct Matrix* add(struct Matrix* A, struct Matrix* B){
    if((A->m != B->m) || (A->n != B->m)){
        printf("Matrices cannot be added. (%d, %d) + (%d, %d) is prohibited.", A->n, A->m, B->n, B->m);
        exit(1);
    }
    struct Matrix* result = make_matrix(A->n, A->m);

    for (int i=0; i<(result->n); i++){
        for (int j=0; j<(result->m); j++){
            result->Matrix[i][j] = A->Matrix[i][j] + B->Matrix[i][j];
        }
    }
    return result;
}

struct Matrix* subtract(struct Matrix* A, struct Matrix* B){
    if((A->m != B->m) || (A->n != B->m)){
        printf("Matrices cannot be subtracted. (%d, %d) - (%d, %d) is prohibited.", A->n, A->m, B->n, B->m);
        exit(1);
    }
    struct Matrix* result = make_matrix(A->n, A->m);

    for (int i=0; i<(result->n); i++){
        for (int j=0; j<(result->m); j++){
            result->Matrix[i][j] = A->Matrix[i][j] - B->Matrix[i][j];
        }
    }
    return result;
}

struct Matrix* hadamard_product(struct Matrix* A, struct Matrix* B){
    if((A->m != B->m) || (A->n != B->m)){
        printf("Matrices cannot be multiplied element-wise. (%d, %d) .* (%d, %d) is prohibited.", A->n, A->m, B->n, B->m);
        exit(1);
    }
    struct Matrix* result = make_matrix(A->n, A->m);

    for (int i=0; i<(result->n); i++){
        for (int j=0; j<(result->m); j++){
            result->Matrix[i][j] = A->Matrix[i][j] * B->Matrix[i][j];
        }
    }
    return result;
}

struct Matrix* transpose(struct Matrix* A){
    struct Matrix* result = make_matrix(A->m, A->n);
    for (int i=0; i<(A->n); i++){
        for (int j=0; j<(A->m); j++){
            result->Matrix[j][i] = A->Matrix[i][j];
        }
    }
    return result;
}

double trace(struct Matrix* A){
    if (A->n != A->m){
        printf("Trace not defined for a non square matrix.");
        exit(1);
    }
    double trace = 0;
    for (int i=0; i<(A->n); i++){
        trace += A->Matrix[i][i];
    }
    return trace;
}

struct Matrix* eye(int n){
    struct Matrix* I = make_matrix(n, n);
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            I->Matrix[i][j] = 0;
            if (i==j){
                I->Matrix[i][j] = 1;
            }
        }
    }
    return I;
}

struct Matrix* concatenate(struct Matrix* A, struct Matrix* B, int dim){
    struct Matrix* result;
    if (dim==0){
        if (A->m != B->m){
            printf("Dimension mismatch along dim=%d\n", dim);
            exit(1);
        }
        result = make_matrix((A->n + B->n), A->m);
        int i;
        for (i=0; i<A->n; i++){
            for (int j=0; j<(A->m); j++){
                result->Matrix[i][j] = A->Matrix[i][j];
            }
        }
        for (;i<(A->n + B->n); i++){
            for (int j=0; j<(B->m); j++){
                result->Matrix[i][j] = B->Matrix[i - (A->n)][j];
            }
        }
    }else if (dim==1){
        if (A->n != B->n){
            printf("Dimension mismatch along dim=%d\n", dim);
            exit(1);
        }
        result = make_matrix(A->n, A->m + B->m);
        int j;
        for(int i=0; i<(A->n); i++){
            j=0;
            for(; j<(A->m); j++){
                result->Matrix[i][j] = A->Matrix[i][j];
            }
            for(; j<(A->m + B->m); j++){
                result->Matrix[i][j] = B->Matrix[i][j - (A->m)];
            }
        }
    }else{
        printf("concatenate(Matrix* A, Matrix* B, int dim) is not defined for dim=%d\n", dim);
        exit(1);
    }
    return result;
}

void row_exchange(struct Matrix* M, int r1, int r2){
    double *temp_row;
    temp_row = M->Matrix[r1];
    M->Matrix[r1] = M->Matrix[r2];
    M->Matrix[r2] = temp_row;
}

void col_exchange(struct Matrix* M, int c1, int c2){
    double *temp;
    for (int i=0; i<(M->n); i++){
        *temp = M->Matrix[i][c1];
        M->Matrix[i][c1] = M->Matrix[i][c2];
        M->Matrix[i][c2] = *temp;
    }
}

double absolute(double val){
    return val > 0.0 ? val: -val;
}

int upper_triangular(struct Matrix* M){
    int n = M->n, m = M->m;
    double pivot;
    int pivot_column = 0, pivot_row;
    int rank = 0;
    for (int i=0; i<n && pivot_column<m; i++){
        pivot_row = i;
        pivot = M->Matrix[pivot_row][pivot_column];
        while(absolute(pivot) < VERY_SMALL_NUMBER && pivot_column < m){
            pivot_row++;
            while(absolute(pivot) < VERY_SMALL_NUMBER && pivot_row < n){
                pivot = M->Matrix[pivot_row][pivot_column];
                pivot_row++;
            }
            if (absolute(pivot) < VERY_SMALL_NUMBER){
                pivot_row = i - 1;
                pivot_column++;
            }else{
                pivot_row--;
            }
        }
        if(absolute(pivot) > VERY_SMALL_NUMBER){
            rank++;
        }else{
            break;
        }
        if (pivot_row != i)
            row_exchange(M, i, pivot_row);

        pivot_row = i;
        // make pivot = 1
        for (int j=0; j<m; j++){
            M->Matrix[i][j] /= pivot;
        }
        
        // make every element below pivot zero.
        for(int k=i+1; k<n; k++){
            double element_below_pivot = M->Matrix[k][pivot_column];
            for (int j=pivot_column; j<m; j++){
                M->Matrix[k][j] -= element_below_pivot * (M->Matrix[i][j]);
            }
        }
        pivot_column++;
    }
    return rank;
}

void eliminate(struct Matrix* M){
    int n = M->n, m = M->m;
    double pivot;
    int pivot_column = 0, pivot_row;
    for (int i=0; i<n && pivot_column<m; i++){
        pivot_row = i;
        pivot = M->Matrix[pivot_row][pivot_column];
        while(absolute(pivot) < VERY_SMALL_NUMBER && pivot_column < m){
            pivot_row++;
            while(absolute(pivot) < VERY_SMALL_NUMBER && pivot_row < n){
                pivot = M->Matrix[pivot_row][pivot_column];
                pivot_row++;
            }
            if (absolute(pivot) < VERY_SMALL_NUMBER){
                pivot_row = i - 1;
                pivot_column++;
            }else{
                pivot_row--;
            }
        }
        if(!(absolute(pivot) > VERY_SMALL_NUMBER)){
            break;
        }
        if (pivot_row != i)
            row_exchange(M, i, pivot_row);

        pivot_row = i;
        // make pivot = 1
        for (int j=0; j<m; j++){
            M->Matrix[i][j] /= pivot;
        }
        
        // make every element below pivot zero.
        for(int k=i+1; k<n; k++){
            double element_below_pivot = M->Matrix[k][pivot_column];
            for (int j=pivot_column; j<m; j++){
                M->Matrix[k][j] -= element_below_pivot * (M->Matrix[i][j]);
            }
        }
        // make every element above pivot zero
        for(int k=0; k<i; k++){
            double element_above_pivot = M->Matrix[k][pivot_column];
            for (int j=pivot_column; j<m; j++){
                M->Matrix[k][j] -= element_above_pivot * M->Matrix[i][j];
            }
        }
        pivot_column++;
    }
}

int rank(struct Matrix* M){
    struct Matrix* B = make_matrix(M->n, M->m);
    for(int i=0; i<(M->n); i++){
        for(int j=0; j<(M->m); j++){
            B->Matrix[i][j] = M->Matrix[i][j];
        }
    }
    int rank = upper_triangular(B);
    free_memory(B);
    return rank;
}

struct Matrix* inverse(struct Matrix* M){
    if (M->n != M->m){
        printf("Matrix inverse does not exist.");
        exit(1);
    }
    struct Matrix* augmented_matrix = concatenate(M, eye(M->n), 1);
    eliminate(augmented_matrix);
    struct Matrix* inv = make_matrix(M->n, M->n);
    for(int i=0; i<(M->n); i++){
        for(int j=0; j<(M->n); j++){
            inv->Matrix[i][j] = augmented_matrix->Matrix[i][j + (M->m)];
        }
    }
    return inv;
}