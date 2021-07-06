#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <synchapi.h>
#include <omp.h>

float **create_empty_matrix(size_t row, size_t col) {
    float **matrix;
    matrix = (float **) malloc(sizeof(float *) * row);
    for (int i = 0; i < row; i++)
        matrix[i] = (float *) malloc(sizeof(float) * col);
    return matrix;
}

float *create_empty_vector(size_t len) {
    float *vector;
    vector = (float *) malloc(sizeof(float) * len);
    return vector;
}

float *create_random_vector(size_t len) {
    float *vector;
    vector = (float *) malloc(sizeof(float) * len);
    srand(time(0));
    for (int i = 0; i < len; ++i) {
        vector[i] = (float) rand() / (float) (RAND_MAX / 5.0);
    }
    return vector;
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

void show_matrix(float **matrix, size_t row, size_t col) {
    puts("");
    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j)
            printf("%-3f ", matrix[i][j]);
        puts("");
    }
}

void show_vector(float *matrix, size_t len) {
    puts("");
    for (size_t i = 0; i < len; ++i) {
        printf("%-3f ", matrix[i]);
        puts("");
    }
}

float **sub_matrix(float **mat1, float **mat2, int row1, int col1, int row2, int col2) {
    if (row1 != row2 || col1 != col2) {
        fprintf(stderr, "incompatible matrix sizes for subtraction\n");
        exit(EXIT_FAILURE);
    }
    float **result = create_empty_matrix(row1, col1);
    int i, j;
#pragma omp parallel shared(result, mat1, mat2) private(i, j)
    {
#pragma omp for schedule(static)
        for (i = 0; i < row1; ++i)
            for (j = 0; j < col1; ++j)
                result[i][j] = mat1[i][j] - mat2[i][j];
    }
    return result;
}

float **mul_matrix_old(float **mat1, float **mat2, size_t row1, size_t col1, size_t row2, size_t col2) {
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

float **mul_matrix_matrix(float **mat1, float **mat2, size_t row1, size_t col1, size_t row2, size_t col2) {
    if (col1 != row2) {
        fprintf(stderr, "incompatible matrix sizes for multiplication\n");
        exit(EXIT_FAILURE);
    }
    float **result = create_empty_matrix(row1, col2);
    int i, j, k;

#pragma omp parallel shared(result, mat1, mat2) private(i, j, k)
    {
#pragma omp for schedule(static)
        for (i = 0; i < row1; i++)
            for (j = 0; j < col2; j++)
                for (k = 0; k < col1; k++)
                    result[i][j] += mat1[i][k] * mat2[k][j];
    }
    return result;
}

float **mul_matrix_scalar(float **mat1, float scalar, size_t row1, size_t col1) {
    float **result = create_empty_matrix(row1, col1);
    int i, j;
#pragma omp parallel shared(result, mat1) private(i, j)
    {
#pragma omp for schedule(static)

        for (int i = 0; i < row1; ++i)
            for (int j = 0; j < col1; ++j)
                result[i][j] = mat1[i][j] * scalar;
    }
    return result;
}

float **transpose_matrix(float **mat, size_t row_v, size_t col_v, size_t *row_vT, size_t *col_vt) {
    *row_vT = col_v;
    *col_vt = row_v;
    float **result = create_empty_matrix(*row_vT, *col_vt);
    for (int i = 0; i < row_v; ++i) {
        for (int j = 0; j < col_v; ++j) {
            result[j][i] = mat[i][j];
        }
    }
    return result;
}

float *mul_matrix_vector(float **matrix, float *vector, size_t row_mat, size_t col_mat, size_t len) {
    if (col_mat != len) {
        fprintf(stderr, "incompatible matrix and vector sizes for multiplication\n");
        exit(EXIT_FAILURE);
    }
    float *result = create_empty_vector(len);
    int i, j;
#pragma omp parallel shared(matrix, result, vector) private(i, j)
    {
#pragma omp for  schedule(static)
        for (i = 0; i < row_mat; i++) {
            for (j = 0; j < col_mat; j++) {
                result[i] += matrix[i][j] * vector[j];
            }
        }
    }
    return result;
}

float *mul_vector_matrix(float *vector, float **matrix, size_t len, size_t row_mat, size_t col_mat) {
    if (len != row_mat) {
        fprintf(stderr, "incompatible matrix and vector sizes for multiplication\n");
        exit(EXIT_FAILURE);
    }
    float *result = create_empty_vector(len);
    int i, j;
#pragma omp parallel shared(matrix, result, vector) private(i, j)
    {
#pragma omp for  schedule(static)
        for (i = 0; i < col_mat; i++) {
            for (j = 0; j < len; j++) {
                result[i] += vector[j] * matrix[j][i];
            }
        }
    }
    return result;
}

float **outer_vector_vector(float *vec1, float *vec2, size_t len1, size_t len2) {
    float **result = create_empty_matrix(len1, len2);
    int i, j;
#pragma omp parallel shared(vec1, vec2, result) private(i, j)
    {
#pragma omp for  schedule(static)
        for (i = 0; i < len1; i++) {
            for (j = 0; j < len2; j++) {
                result[i][j] = vec1[i] * vec2[j];
            }
        }
    }
    return result;
}

float inner_vector_vector(float *vec1, float *vec2, size_t len1, size_t len2) {
    float result;
    int i;
#pragma omp parallel shared(vec1, vec2, result) private(i)
    {
#pragma omp for  schedule(static)
        for (i = 0; i < len1; i++) {
            result += vec1[i] * vec2[i];
        }
    }
    return result;
}

/**
 * this function calculate inverse of rank-1 update of matrix A.
 * The rank-1 update of A is shown by (A+uv^T)
 * If the rank-1 update of A is invertible, then it will be calculated by:
 * (A+(u)(v^T))^(-1) = (A^-1) - ((A^-1)(u)(v^T)(A^-1))/(1+(v^T)(A^-1)(u))
 * numerator: (A^-1)(u)(v^T)(A^-1)
 * denominator: 1 + (v^T)(A^-1)(u)
 * @param A is a K*K invertible matrix
 * @param u is K*1 vector
 * @param v is K*1 vector
 * @return
 */
float **
sherman_morrison(float **A_inv, float *u, float *v, size_t row_A_inv, size_t col_A_inv, size_t row_u, size_t col_u,
                 size_t row_vT, size_t col_vT) {

    double start_time = omp_get_wtime();

    //num1 = (A^-1)(u) --size--> n*1
    //num2 = (v^T)(A^-1) --size--> 1*n
    //numerator --size--> n*n
    float *num1 = mul_matrix_vector(A_inv, u, row_A_inv, col_A_inv, row_u * col_u);
    float *num2 = mul_vector_matrix(v, A_inv, row_vT * col_vT, row_A_inv, col_A_inv);
    float **numerator = outer_vector_vector(num1, num2, row_u, row_u);

    //den1 = (v^T)(A^-1) --size--> 1*n
    float *den1 = num2;
    float denominator = 1 + inner_vector_vector(den1, u, row_u, row_u);

    float **second_term = mul_matrix_scalar(numerator, 1.0 / denominator, row_u, row_u);


    float **result = sub_matrix(A_inv, second_term, row_u, row_u, row_u, row_u);

    double end_time = omp_get_wtime();
    printf("elapsed time in sec: %f", end_time - start_time);

    return result;

}

float **read_matrix(size_t *rows, size_t *cols, const char *filename) {
    if (rows == NULL || cols == NULL || filename == NULL)
        return NULL;

    *rows = 0;
    *cols = 0;

    FILE *fp = fopen(filename, "r");

    if (fp == NULL) {
        fprintf(stderr, "could not open %s: %s\n", filename, strerror(errno));
        return NULL;
    }

    float **matrix = NULL, **tmp;

    char line[1024];

    while (fgets(line, sizeof line, fp)) {
        if (*cols == 0) {
            // determine the size of the columns based on
            // the first row
            char *scan = line;
            int dummy;
            int offset = 0;
            while (sscanf(scan, "%d%n", &dummy, &offset) == 1) {
                scan += offset;
                (*cols)++;
            }
        }

        tmp = realloc(matrix, (*rows + 1) * sizeof *matrix);

        if (tmp == NULL) {
            fclose(fp);
            return matrix; // return all you've parsed so far
        }

        matrix = tmp;

        matrix[*rows] = calloc(*cols, sizeof *matrix[*rows]);

        if (matrix[*rows] == NULL) {
            fclose(fp);
            if (*rows == 0) // failed in the first row, free everything
            {
                fclose(fp);
                free(matrix);
                return NULL;
            }

            return matrix; // return all you've parsed so far
        }

        int offset = 0;
        char *scan = line;
        for (size_t j = 0; j < *cols; ++j) {
            if (sscanf(scan, "%f%n", matrix[*rows] + j, &offset) == 1)
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

float *mat_to_vec(float **mat, size_t row_mat, size_t col_mat, size_t len) {
    float *vector = (float *) malloc(sizeof(float) * len);
    for (size_t i = 0; i < row_mat; ++i) {
        for (size_t j = 0; j < col_mat; ++j) {
            vector[i] = mat[i][j];
        }
    }
    return vector;
}

/**
 * main method
 */
void main() {

    size_t row_A_inv, col_A_inv, row_u, col_u, row_v, col_v, row_vT, col_vT;
    float *u, *v, *vT;
    float **A_inv = read_matrix(&row_A_inv, &col_A_inv, "A_inv.dat");
    float **u_matrix = read_matrix(&row_u, &col_u, "u.dat");
    float **v_matrix = read_matrix(&row_v, &col_v, "v.dat");
//    float **vT_matrix = transpose_matrix(v, row_v, col_v, &row_vT, &col_vT);

    u = mat_to_vec(u_matrix, row_u, col_u, row_u * col_u);
    v = mat_to_vec(v_matrix, row_v, col_v, row_v * col_v);

/*    show_vector(u, row_u);
    show_vector(v, row_v);*/

    float **morrison = sherman_morrison(A_inv, u, v, row_A_inv, col_A_inv, row_u, col_u, row_v, col_v);
    show_matrix(morrison, row_A_inv, col_A_inv);


}