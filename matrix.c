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
        fprintf(stderr, "incompatible matrix sizes for subtraction\n");
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

float **mul_matrix_scalar(float **mat1, float scalar, size_t row1, size_t col1) {
    float **result = create_empty_matrix(row1, col1);
    for (int i = 0; i < row1; ++i)
        for (int j = 0; j < col1; ++j)
            result[i][j] = mat1[i][j] * scalar;
    return result;
}

float **transpose_matrix(float **mat, size_t row_v, size_t col_v, size_t *row_vT, size_t *col_vt){
    *row_vT = col_v;
    *col_vt = row_v;
    float **result = create_empty_matrix(row_vT, col_vt);
    for (int i = 0; i < row_v; ++i) {
        for (int j = 0; j < col_v; ++j) {
            result[j][i] = mat[i][j];
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
 * denominator: 1 + (v^T)(A^-1)(u)
 * @param A is a K*K invertible matrix
 * @param u is K*1 vector
 * @param v is K*1 vector
 * @return
 */
float **sherman_morrison(float **A_inv, float **u, float **v, size_t row_u, size_t col_u) {

    float **result = create_empty_matrix(row_u, row_u);

    float **num1 = mul_matrix(A_inv, u, row_u, row_u, row_u, col_u);
    float **num2 = mul_matrix(v, A_inv, col_u, row_u, row_u, row_u);
    float **numerator = mul_matrix(num1, num2, row_u, col_u, col_u, row_u);

    float **den1 = mul_matrix(v, A_inv, col_u, row_u, row_u, row_u);
    float **den2 = mul_matrix(den1, u, col_u, row_u, row_u, col_u);
    float denominator = 1 + **den2;

    float **second_term = mul_matrix_scalar(numerator, 1.0 / denominator, row_u, row_u);

    result = add_matrix(A_inv, second_term, row_u, row_u, row_u, row_u);

    return result;
}

float **read_matrix(size_t *rows, size_t *cols, const char *filename)
{
    if(rows == NULL || cols == NULL || filename == NULL)
        return NULL;

    *rows = 0;
    *cols = 0;

    FILE *fp = fopen(filename, "r");

    if(fp == NULL)
    {
        fprintf(stderr, "could not open %s: %s\n", filename, strerror(errno));
        return NULL;
    }

    float **matrix = NULL, **tmp;

    char line[1024];

    while(fgets(line, sizeof line, fp))
    {
        if(*cols == 0)
        {
            // determine the size of the columns based on
            // the first row
            char *scan = line;
            int dummy;
            int offset = 0;
            while(sscanf(scan, "%d%n", &dummy, &offset) == 1)
            {
                scan += offset;
                (*cols)++;
            }
        }

        tmp = realloc(matrix, (*rows + 1) * sizeof *matrix);

        if(tmp == NULL)
        {
            fclose(fp);
            return matrix; // return all you've parsed so far
        }

        matrix = tmp;

        matrix[*rows] = calloc(*cols, sizeof *matrix[*rows]);

        if(matrix[*rows] == NULL)
        {
            fclose(fp);
            if(*rows == 0) // failed in the first row, free everything
            {
                fclose(fp);
                free(matrix);
                return NULL;
            }

            return matrix; // return all you've parsed so far
        }

        int offset = 0;
        char *scan = line;
        for(size_t j = 0; j < *cols; ++j)
        {
            if(sscanf(scan, "%f%n", matrix[*rows] + j, &offset) == 1)
                scan += offset;
            else
                matrix[*rows][j] = 0; // could not read, set cell to 0
        }

        // incrementing rows
        (*rows)++;
    }

    fclose(fp);

    return matrix;
}

/**
 * main method
 */
void main() {

    size_t row_A_inv, col_A_inv, row_u, col_u, row_v, col_v, row_vT, col_vT;
    float **A_inv = read_matrix(&row_A_inv, &col_A_inv, "A_inv.dat");
    float **u = read_matrix(&row_u, &col_u, "u.dat");
    float **v =read_matrix(&row_v, &col_v, "v.dat");
    float **vT = transpose_matrix(v, row_v, col_v, &row_vT, &col_v);
    show_matrix(vT, row_vT, col_vT);

//    show_matrix(mul_matrix(u, v, row_u, col_u, row_v, col_v),row_u,col_v);

//    float **morrison = sherman_morrison(A_inv, u, v, 3, 1);
//    show_matrix(morrison, 3, 3);



}