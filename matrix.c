#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>

/**
 *  this method shows a vector as a matrix by row*col
 * @param vector
 * @param row
 * @param col
 */
void show_vector(float *vector, int row, int col) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            printf("%f ", vector[i * col + j]);
        }
        printf("\n");
    }
}

void first_n(int number, int **vector) {
    int i;
    *vector = (int *) malloc(number * sizeof(int));
    for (i = 0; i < number; i++) {
        (*vector)[i] = i;
    }
}

/**
 *
 * @param rows
 * @param cols
 * @return  a vector of float initialized to the zero value
 */
float *init_vector(int rows, int cols, float val) {
    float *vector = (float *) calloc(rows * cols, sizeof(float));
    for (int i = 0; i < rows * cols; ++i) {
        vector[i] = val;
    }
    return vector;
}

float *add_vector(float *vec1, float *vec2, int vec_length) {
    float *result = (float *) calloc(vec_length, sizeof(float));
    for (int i = 0; i < vec_length; ++i) {
        result[i] = vec1[i] + vec2[i];
    }
    return result;
}

float *sub_vector(float *vec1, float *vec2, int vec_length) {
    float *result = (float *) calloc(vec_length, sizeof(float));
    for (int i = 0; i < vec_length; ++i) {
        result[i] = vec1[i] - vec2[i];
    }
    return result;
}

float *neg_vector(float *vec, int vec_length) {
    float *result = (float *) calloc(vec_length, sizeof(float));
    for (int i = 0; i < vec_length; ++i) {
        result[i] = -vec[i];
    }
    return result;
}

float *mul_vector(float *vec1, float *vec2, int row1, int col1, int col2) {
    float *result = (float *) calloc(row1 * col2, sizeof(float));
    for (int i = 0; i < row1; ++i) {
        for (int j = 0; j < col2; ++j) {
            float t = 0;
            for (int k = 0; k < col1; ++k) {
                t += vec1[(i * col1) + k] * vec2[(k * col2) + j];
            }
            result[(i * col2) + j] = t;
        }
    }
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

    float *p_denominator = mul_vector(mul_vector(v, A_inv, col_u, row_u, row_u), u, col_u, row_u, col_u);
    float denominator = 1 + *p_denominator;
    if (denominator == 0) {
        return NULL;
    }
    float *p_numerator =
            mul_vector(mul_vector(mul_vector(A_inv, u, row_u, row_u, col_u), v, row_u, col_u, row_u), A_inv, row_u,
                       row_u, row_u);
    float *rank_one_update = sub_vector(A_inv, p_numerator, row_u * row_u);
    return rank_one_update;
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

}