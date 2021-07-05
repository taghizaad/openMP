#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <synchapi.h>

float **create_empty_matrix(int row, int col) {
    float **matrix;
    matrix = (float **) malloc(sizeof(float *) * row);
    for (int i = 0; i < row; i++)
        matrix[i] = (float *) malloc(sizeof(float) * col);
    return matrix;
}

float **create_random_matrix(int row, int col) {
    float **matrix;
    matrix = (float **) malloc(sizeof(float *) * row);
    for (int i = 0; i < row; i++)
        matrix[i] = (float *) malloc(sizeof(float) * col);
    /* Initializes random number generator */
    srand(time(0));
    for (size_t i = 0; i < row; ++i)
        for (size_t j = 0; j < col; ++j)
            matrix[i][j] = (float) rand() / (float) (RAND_MAX / 5.0);
    return matrix;
}

void show_matrix(float **matrix, int row, int col) {
    printf("size of matrix is %d*%d", row, col);
    puts("");
    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j)
            printf("%-3f ", matrix[i][j]);
        puts("");
    }
}

float **add_matrix(float **mat1, float **mat2, int row1, int col1, int row2, int col2) {
    if (row1 != row2 || col1 != col2) {
        fprintf(stderr, "incompatible matrix sizes for addition\n");
        exit(EXIT_FAILURE);
    }
    float **result = create_empty_matrix(row1, col1);
    for (size_t i = 0; i < row1; ++i)
        for (size_t j = 0; j < col1; ++j)
            result[i][j] = mat1[i][j] + mat2[i][j];
    return result;
}

float **sub_matrix(float **mat1, float **mat2, int row1, int col1, int row2, int col2) {
    if (row1 != row2 || col1 != col2) {
        fprintf(stderr, "incompatible matrix sizes for addition\n");
        exit(EXIT_FAILURE);
    }
    float **result = create_empty_matrix(row1, col1);
    for (size_t i = 0; i < row1; ++i)
        for (size_t j = 0; j < col1; ++j)
            result[i][j] = mat1[i][j] - mat2[i][j];
    return result;
}

float **mul_matrix(float **mat1, float **mat2, size_t row1, size_t col1, size_t row2, size_t col2) {
    if (col1 != row2) {
        fprintf(stderr, "incompatible matrix sizes for multiplication\n");
        exit(EXIT_FAILURE);
    }
    float **result = create_empty_matrix(row1, col2);
    for (int i = 0; i < row1; ++i)
        for (int j = 0; j < col2; ++j)
            for (int k = 0; k < col1; ++k)
                result[i][j] += mat1[i][k] * mat2[k][j];
    return result;
}

/**
 * this function calculate inverse of rank-1 update of matrix A.
 * The rank-1 update of A is shown by (A+uv^T)
 * If the rank-1 update of A is invertible, then it will be calculated by:
 * (A+(u)(v^T))^(-1) = (A^-1) + ((A^-1)(u)(v^T)(A^-1))/(1+(v^T)(A^-1)(u))
 * numerator: (A^-1)(u)(v^T)(A^-1)
 * denominator: (v^T)(A^-1)(u)
 * @param A is a K*K invertible matrix
 * @param u is K*1 vector
 * @param v is K*1 vector
 * @return
 */
float *sherman_morrison(float *A_inv, float *u, float *v, int row_u, int col_u) {

/*    float *p_denominator = mul_matrix(mul_matrix(v, A_inv, col_u, row_u, row_u), u, col_u, row_u, col_u);
    float denominator = 1 + *p_denominator;
    if (denominator == 0) {
        return NULL;
    }
    float *p_numerator =
            mul_matrix(mul_matrix(mul_matrix(A_inv, u, row_u, row_u, col_u), v, row_u, col_u, row_u), A_inv, row_u,
                       row_u, row_u);
    float *rank_one_update = sub_matrix(A_inv, p_numerator, row_u * row_u);
    return rank_one_update;*/
}


/**
 * main method
 */
void main() {

/*    float *vector1;
    float *vector2;
    int rows = 2, cols = 3;
    int vec_len = rows * cols;
    vector1 = init_vector(rows, cols, 5);
    vector2 = init_vector(rows, cols, 7);
    if (vector1 != NULL && vector2 != NULL) {
        float *pDouble = add_vector(vector1, neg_vector(vector2, vec_len), vec_len);
        show_vector(pDouble, rows, cols);

    }*/

/*
    float *vec1, *vec2;
    int row1 = 3, col1 = 1, col2 = 2;
    vec1 = init_vector(row1, col1, 5);
    vec2 = init_vector(col1, col2, 7);
    show_vector(vec1, row1, col1);
    printf("---------------\n");
    show_vector(vec2, col1, col2);
    printf("---------------\n");
    float *mul_vec = mul_vector(vec1, vec2, row1, col1, col2);
    show_vector(mul_vec, row1, col2);
*/
    float **mat1 = create_random_matrix(3, 7);
    show_matrix(mat1, 3, 7);
    Sleep(5000);
    float **mat2 = create_random_matrix(7, 7);
    show_matrix(mat2, 3, 7);
    float **addMatrix = sub_matrix(mat1, mat2, 3, 7, 3, 7);
    show_matrix(addMatrix, 3, 7);


}