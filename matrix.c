#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <synchapi.h>
#include <omp.h>

double **create_empty_matrix(size_t row, size_t col) {
    double **matrix;
    matrix = (double **) malloc(sizeof(double *) * row);
    for (int i = 0; i < row; i++)
        matrix[i] = (double *) malloc(sizeof(double) * col);
    return matrix;
}

double *create_empty_vector(size_t len) {
    double *vector;
    vector = (double *) malloc(sizeof(double) * len);
    return vector;
}

double *create_random_vector(size_t len) {
    double *vector;
    vector = (double *) malloc(sizeof(double) * len);
    srand(time(NULL));
    for (int i = 0; i < len; i++) {
        vector[i] = 5.0 * rand() / RAND_MAX;
    }
    return vector;
}

double *create_random_boolean_vector(size_t len) {
    double *vector;
    vector = (double *) malloc(sizeof(double) * len);
    srand(time(NULL));
    for (int i = 0; i < len; i++) {
        vector[i] = rand() & 1;
    }
    return vector;
}

double **create_random_matrix(int row, int col) {
    double **matrix;
    matrix = (double **) malloc(sizeof(double *) * row);
    for (int i = 0; i < row; i++)
        matrix[i] = (double *) malloc(sizeof(double) * col);
    /* Initializes random number generator */
    srand(time(NULL));
    for (size_t i = 0; i < row; i++)
        for (size_t j = 0; j < col; j++)
            matrix[i][j] = 5.0 * rand() / RAND_MAX;
    return matrix;
}

void show_matrix(double **matrix, size_t row, size_t col) {
    puts("");
    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < col; j++)
            printf("%f ", matrix[i][j]);
        puts("");
    }
}

void show_vector(double *matrix, size_t len) {
    puts("");
    for (size_t i = 0; i < len; i++) {
        printf("%lf ", matrix[i]);
        puts("");
    }
}

double **sub_matrix_parallel(double **mat1, double **mat2, int row1, int col1, int row2, int col2) {
    if (row1 != row2 || col1 != col2) {
        fprintf(stderr, "incompatible matrix sizes for subtraction\n");
        exit(EXIT_FAILURE);
    }
    double **result = create_empty_matrix(row1, col1);
#pragma omp parallel for
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col1; j++)
            result[i][j] = mat1[i][j] - mat2[i][j];
    return result;
}

double **sub_matrix_sequential(double **mat1, double **mat2, int row1, int col1, int row2, int col2) {
    if (row1 != row2 || col1 != col2) {
        fprintf(stderr, "incompatible matrix sizes for subtraction\n");
        exit(EXIT_FAILURE);
    }
    double **result = create_empty_matrix(row1, col1);
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col1; j++)
            result[i][j] = mat1[i][j] - mat2[i][j];
    return result;
}

double **add_matrix_parallel(double **mat1, double **mat2, int row1, int col1, int row2, int col2) {
    if (row1 != row2 || col1 != col2) {
        fprintf(stderr, "incompatible matrix sizes for subtraction\n");
        exit(EXIT_FAILURE);
    }
    double **result = create_empty_matrix(row1, col1);
#pragma omp parallel for
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col1; j++)
            result[i][j] = mat1[i][j] + mat2[i][j];
    return result;
}

double **add_matrix_sequential(double **mat1, double **mat2, int row1, int col1, int row2, int col2) {
    if (row1 != row2 || col1 != col2) {
        fprintf(stderr, "incompatible matrix sizes for subtraction\n");
        exit(EXIT_FAILURE);
    }
    double **result = create_empty_matrix(row1, col1);
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col1; j++)
            result[i][j] = mat1[i][j] + mat2[i][j];
    return result;
}

double **mul_matrix_old(double **mat1, double **mat2, size_t row1, size_t col1, size_t row2, size_t col2) {
    if (col1 != row2) {
        fprintf(stderr, "incompatible matrix sizes for multiplication\n");
        exit(EXIT_FAILURE);
    }
    double **result = create_empty_matrix(row1, col2);
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col2; j++)
            for (int k = 0; k < col1; k++)
                result[i][j] += mat1[i][k] * mat2[k][j];
    return result;
}

double **mul_matrix_matrix(double **mat1, double **mat2, size_t row1, size_t col1, size_t row2, size_t col2) {
    if (col1 != row2) {
        fprintf(stderr, "incompatible matrix sizes for multiplication\n");
        exit(EXIT_FAILURE);
    }
    double **result = create_empty_matrix(row1, col2);

#pragma omp parallel for
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col2; j++)
            for (int k = 0; k < col1; k++)
                result[i][j] += mat1[i][k] * mat2[k][j];
    return result;
}

double **mul_matrix_scalar_parallel(double **mat1, double scalar, size_t row1, size_t col1) {
    double **result = create_empty_matrix(row1, col1);
#pragma omp parallel for
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col1; j++)
            result[i][j] = mat1[i][j] * scalar;

    return result;
}

double **mul_matrix_scalar_sequential(double **mat1, double scalar, size_t row1, size_t col1) {
    double **result = create_empty_matrix(row1, col1);
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col1; j++)
            result[i][j] = mat1[i][j] * scalar;
    return result;
}

double **transpose_matrix_sequential(double **mat, size_t row_v, size_t col_v) {

    double **result = create_empty_matrix(col_v, row_v);
    for (int i = 0; i < row_v; i++) {
        for (int j = 0; j < col_v; j++) {
            result[j][i] = mat[i][j];
        }
    }
    return result;
}

double *transpose(double *d, size_t row, size_t col) {
    double *t = create_empty_vector(row * col);
#pragma omp parallel for
    for (int i = 0; i < row; ++i)
        for (int j = 0; j < col; ++j)
            t[row * j + i] = d[col * i + j];
    return t;
}

double **transpose_matrix_parallel(double **mat, size_t row_v, size_t col_v) {

    double **result = create_empty_matrix(col_v, row_v);

#pragma omp parallel for
    for (int i = 0; i < row_v; i++)
        for (int j = 0; j < col_v; j++)
            result[j][i] = mat[i][j];
    return result;
}

double *mul_matrix_vector_parallel(double **matrix, double *vector, size_t row_mat, size_t col_mat, size_t len) {
    if (col_mat != len) {
        fprintf(stderr, "incompatible matrix and vector sizes for multiplication\n");
        exit(EXIT_FAILURE);
    }
    double *result = create_empty_vector(len);
#pragma omp parallel for
    for (int i = 0; i < row_mat; i++)
        for (int j = 0; j < col_mat; j++)
            result[i] += matrix[i][j] * vector[j];
    return result;
}

double *mul_matrix_vector_sequential(double **matrix, double *vector, size_t row_mat, size_t col_mat, size_t len) {
    double *result = create_empty_vector(len);
    for (int i = 0; i < row_mat; i++)
        for (int j = 0; j < col_mat; j++)
            result[i] += matrix[i][j] * vector[j];
    return result;
}


double *mul_vector_matrix_parallel(double **matrix, double *vector, size_t row_mat, size_t col_mat, size_t len) {
    if (len != row_mat) {
        fprintf(stderr, "incompatible matrix and vector sizes for multiplication\n");
        exit(EXIT_FAILURE);
    }
    double *result = create_empty_vector(len);
#pragma omp parallel for
    for (int i = 0; i < col_mat; i++)
        for (int j = 0; j < len; j++)
            result[i] += vector[j] * matrix[j][i];
    return result;
}

double *mul_vector_matrix_sequential(double **matrix, double *vector, size_t row_mat, size_t col_mat, size_t len) {
    if (len != row_mat) {
        fprintf(stderr, "incompatible matrix and vector sizes for multiplication\n");
        exit(EXIT_FAILURE);
    }
    double *result = create_empty_vector(len);
    for (int i = 0; i < col_mat; i++)
        for (int j = 0; j < len; j++)
            result[i] += vector[j] * matrix[j][i];
    return result;
}

double *mul_vector_matrix_parallel_with_matrix_transpose(double **matrix, double *vector, size_t row_mat, size_t col_mat,
                                                        size_t len) {
    if (len != row_mat) {
        fprintf(stderr, "incompatible matrix and vector sizes for multiplication\n");
        exit(EXIT_FAILURE);
    }
    double *result = create_empty_vector(len);
//    double **matrix_transpose = transpose_matrix_parallel(matrix, row_mat, col_mat);
#pragma omp parallel for
    for (int i = 0; i < col_mat; i++)
        for (int j = 0; j < len; j++)
            result[i] += matrix[j][i] * vector[j];
    return result;
}

double *mul_vector_matrix_sequential_with_matrix_transpose(double **matrix, double *vector, size_t row_mat, size_t col_mat,
                                                          size_t len) {
    if (len != row_mat) {
        fprintf(stderr, "incompatible matrix and vector sizes for multiplication\n");
        exit(EXIT_FAILURE);
    }
    double *result = create_empty_vector(len);
    double **matrix_transpose = transpose_matrix_sequential(matrix, row_mat, col_mat);
    for (int i = 0; i < col_mat; i++)
        for (int j = 0; j < len; j++)
            result[i] += matrix_transpose[i][j] * vector[j];
    return result;
}

double **outer_vector_vector_parallel(double *vec1, double *vec2, size_t len1, size_t len2) {
    double **result = create_empty_matrix(len1, len2);
#pragma omp parallel for
    for (int i = 0; i < len1; i++)
        for (int j = 0; j < len2; j++)
            result[i][j] = vec1[i] * vec2[j];
    return result;
}

double **outer_vector_vector_sequential(double *vec1, double *vec2, size_t len1, size_t len2) {
    double **result = create_empty_matrix(len1, len2);
    for (int i = 0; i < len1; i++)
        for (int j = 0; j < len2; j++)
            result[i][j] = vec1[i] * vec2[j];
    return result;
}

double inner_vector_vector_parallel(double *vec1, double *vec2, size_t len1) {
    double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < len1; i++)
        sum += vec1[i] * vec2[i];
    return sum;
}

double inner_vector_vector_sequential(double *vec1, double *vec2, size_t len1) {
    double sum=0.0;
    for (int i = 0; i < len1; i++) {
        sum = sum + vec1[i] * vec2[i];
    }
    return sum;
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
double **sherman_morrison_parallel(double **A_inv, double *u, double *v, size_t size) {

    //num1 = (A^-1)(u) --size--> n*1
    //num2 = (v^T)(A^-1) --size--> 1*n
    //numerator --size--> n*n
    double *num1 = mul_matrix_vector_parallel(A_inv, u, size, size, size);
    double *num2 = mul_vector_matrix_parallel_with_matrix_transpose(A_inv, v, size, size, size);
    double **numerator = outer_vector_vector_parallel(num1, num2, size, size);

    //den1 = (v^T)(A^-1) = num2 --size--> 1*n
    double *den1 = num2;
    double denominator = 1 + inner_vector_vector_parallel(den1, u, size);

    double **second_term = mul_matrix_scalar_parallel(numerator, 1.0 / denominator, size, size);


    double **result = sub_matrix_parallel(A_inv, second_term, size, size, size, size);
    return result;
}

double **sherman_morrison_sequential(double **A_inv, double *u, double *v, size_t size) {

    //num1 = (A^-1)(u) --size--> n*1
    //num2 = (v^T)(A^-1) --size--> 1*n
    //numerator --size--> n*n
    double *num1 = mul_matrix_vector_sequential(A_inv, u, size, size, size);
    double *num2 = mul_vector_matrix_sequential(A_inv, v, size, size, size);
    double **numerator = outer_vector_vector_sequential(num1, num2, size, size);

    //den1 = (v^T)(A^-1) = num2 --size--> 1*n
    double *den1 = num2;
    double denominator = 1 + inner_vector_vector_sequential(den1, u, size);

    double **second_term = mul_matrix_scalar_sequential(numerator, 1.0 / denominator, size, size);


    double **result = sub_matrix_sequential(A_inv, second_term, size, size, size, size);
    return result;
}

double **read_matrix(size_t *rows, size_t *cols, const char *filename) {
    if (rows == NULL || cols == NULL || filename == NULL)
        return NULL;

    *rows = 0;
    *cols = 0;

    FILE *fp = fopen(filename, "r");

    if (fp == NULL) {
        fprintf(stderr, "could not open %s: %s\n", filename, strerror(errno));
        return NULL;
    }

    double **matrix = NULL, **tmp;

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
        for (size_t j = 0; j < *cols; j++) {
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

double *mat_to_vec(double **mat, size_t row_mat, size_t col_mat) {
    double *vector = (double *) malloc(sizeof(double) * row_mat * col_mat);
    for (size_t i = 0; i < row_mat; i++) {
        for (size_t j = 0; j < col_mat; j++) {
            vector[i * col_mat + j] = mat[i][j];
        }
    }
    return vector;
}

double **gauss_jordan_matrix_inverse(double **A, size_t size) {
    double **I, temp;
    int i, j, k;

    I = (double **) malloc(
            size * sizeof(double *));            //memory allocation for indentity matrix I(matsize X matsize)
    for (i = 0; i < size; i++)
        I[i] = (double *) malloc(size * sizeof(double));
    for (i = 0; i < size; i++)                                  //automatically initialize the unit matrix, e.g.
        for (j = 0; j < size; j++)                              //  -       -
            if (i == j)                                        // | 1  0  0 |
                I[i][j] = 1;                                  // | 0  1  0 |
            else                                            // | 0  0  1 |
                I[i][j] = 0;                                  //  -       -

    /*---------------LoGiC starts here------------------*/      //procedure // to make the matrix A to unit matrix
    for (k = 0; k < size; k++)                                  //by some row operations,and the same row operations of
    {                                                       //Unit mat. I gives the inverse of matrix A
        temp = A[k][k];                   //'temp'
        // stores the A[k][k] value so that A[k][k]  will not change
        for (j = 0; j < size; j++)      //during the operation //A[i] //[j]/=A[k][k]  when i=j=k
        {
            A[k][j] /= temp;                                  //it performs // the following row operations to make A to unit matrix
            I[k][j] /= temp;                                  //R0=R0/A[0][0],similarly for I also R0=R0/A[0][0]
        }                                                   //R1=R1-R0*A[1][0] similarly for I
        for (i = 0; i < size; i++)                              //R2=R2-R0*A[2][0]      ,,
        {
            temp = A[i][k];                       //R1=R1/A[1][1]
            for (j = 0; j < size; j++)             //R0=R0-R1*A[0][1]
            {                                   //R2=R2-R1*A[2][1]
                if (i == k)
                    break;                      //R2=R2/A[2][2]
                A[i][j] -= A[k][j] * temp;          //R0=R0-R2*A[0][2]
                I[i][j] -= I[k][j] * temp;          //R1=R1-R2*A[1][2]
            }
        }
    }
    return I;
}

/**
 * main method
 */
void main() {

    size_t size = 64;
    size_t niter = 100000;
    double micro = .000001;
    double *v = create_random_boolean_vector(size);
//    printf("--------v-----------");
//    show_vector(v, size);
//    Sleep(1000);
    double *u = create_random_vector(size);
//    printf("--------u-----------");
//    show_vector(u, size);
//    double **uvT = outer_vector_vector_sequential(u, v, size, size);
//    Sleep(1000);
    double **A = create_random_matrix(size, size);
//    double **Anew = add_matrix_sequential(A, uvT, size, size, size, size);
//    printf("----------A-----------");
//    show_matrix(A, size, size);
    double start_Ainv_GJ = omp_get_wtime();
    double **Ainv_GJ = gauss_jordan_matrix_inverse(A, size);
    double end_Ainv_GJ = omp_get_wtime();
//    printf("----------Ainv-----------");
//    show_matrix(Ainv_GJ, size, size);
//    double start_Anewinv_GJ = omp_get_wtime();
//    double **Anewinv_GJ = gauss_jordan_matrix_inverse(Anew, size);
//    double end_Anewinv_GJ = omp_get_wtime();

    double start_Anewinv_SM_seq = omp_get_wtime();
#pragma omp parallel for
    for (int i = 0; i < niter; i++) {
        Ainv_GJ = sherman_morrison_sequential(Ainv_GJ, u, v, size);
    }
    double end_Anewinv_SM_seq = omp_get_wtime();

//    double start_Anewinv_SM_par = omp_get_wtime();
//    double **Anewinv_SM_par = sherman_morrison_parallel(Ainv_GJ, u, v, size);
//    double end_Anewinv_SM_par =omp_get_wtime();
//    printf("--------Ainv%i-----------",i);
//    show_matrix(Ainv_GJ, size,size);
    printf("Matrix size: %d * %d\n", size, size);
    printf("Niter: %d\n", niter);
    printf("Ainv_GJ time: %f\n", end_Ainv_GJ - start_Ainv_GJ);
    printf("SM_seq time: %f s\n", end_Anewinv_SM_seq - start_Anewinv_SM_seq);
    printf("Total SM_seq time: %f \xC2\xB5s\n", (end_Anewinv_SM_seq - start_Anewinv_SM_seq)/micro);
    printf("Average SM_seq time: %f \xC2\xB5s\n", (end_Anewinv_SM_seq - start_Anewinv_SM_seq)/micro/niter);






//    printf("some Anewinv values for test GJ, SMseq, SMpar: %f, %f, %f\n",
//           Anewinv_GJ[38][23],
//           Anewinv_SM_seq[38][23],
//           Anewinv_SM_par[38][23]);

}