
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


double **create_empty_matrix_double(size_t row, size_t col) {
    double **matrix;
    matrix = (double **) calloc(row, sizeof(double *));
    for (int i = 0; i < row; i++)
        matrix[i] = (double *) calloc(col, sizeof(double));
    return matrix;
}

int **create_empty_matrix_int(size_t row, size_t col) {
    int **matrix;
    matrix = (int **) calloc(row, sizeof(int *));
    for (int i = 0; i < row; i++)
        matrix[i] = (int *) calloc(col, sizeof(int));
    return matrix;
}

double **sub_matrix(double **mat1, double **mat2, int row1, int col1, int row2, int col2) {
    if (row1 != row2 || col1 != col2) {
        fprintf(stderr, "incompatible matrix sizes for subtraction\n");
        exit(EXIT_FAILURE);
    }
    double **result = create_empty_matrix_double(row1, col1);
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col1; j++)
            result[i][j] = mat1[i][j] - mat2[i][j];
    return result;
}

double **add_matrix(double **mat1, double **mat2, int row1, int col1, int row2, int col2) {
    if (row1 != row2 || col1 != col2) {
        fprintf(stderr, "incompatible matrix sizes for subtraction\n");
        exit(EXIT_FAILURE);
    }
    double **result = create_empty_matrix_double(row1, col1);
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col1; j++)
            result[i][j] = mat1[i][j] + mat2[i][j];
    return result;
}

double **mul_matrix_matrix(double **mat1, double **mat2, size_t row1, size_t col1, size_t row2, size_t col2) {
    if (col1 != row2) {
        fprintf(stderr, "incompatible matrix sizes for multiplication\n");
        exit(EXIT_FAILURE);
    }
    double **result = create_empty_matrix_double(row1, col2);

    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col2; j++)
            for (int k = 0; k < col1; k++)
                result[i][j] += mat1[i][k] * mat2[k][j];
    return result;
}

void show_matrix_double(double **matrix, size_t row, size_t col) {
    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < col; j++)
            printf("%f ", matrix[i][j]);
        puts("");
    }
}

void show_matrix_int(int **matrix, size_t row, size_t col) {
    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < col; j++)
            printf("%d ", matrix[i][j]);
        puts("");
    }
}

void show_vector_int(int *matrix, size_t len) {
    for (size_t i = 0; i < len; i++) {
        printf("%d ", matrix[i]);
        puts("");
    }
}

int readFileDouble(double **mat, int size, char *fileName) {
    int arr_row = size;
    char buffer[1024];
    char *tok;
    int column;
    double test_var;

    FILE *gain_ptr = fopen(fileName, "r");

    if (gain_ptr == NULL) {
        printf("Error Reading File\n");
        return 1; //return with failure
    }

    for (int row = 0; row < arr_row; row++) {
        fgets(buffer, 1024, gain_ptr); //buffer is declared as char buffer[1024];

        tok = strtok(buffer, "\t");    // NULL means 'continue from last token'
        //tok is declared as char *tok;
        //the numbers in my text file are seperated by tabs
        //so the tokens should be separated by \t

        column = 0;  //since fgets gets a whole line I had to separate the columns
        while (tok != NULL) {
            test_var = atof(tok);  //test variable is declared as double
            //and it is used to hold temporarily the value of tok in double
            tok = strtok(NULL, "\t");    // NULL means 'continue from last token'

            mat[row][column] = test_var;
            column++;
        }
    }
    fclose(gain_ptr);
    return 1;
}

int readFileInt(int **mat, int size, char *fileName) {
    int arr_row = size;
    char buffer[1024];
    char *tok;
    int column;
    int test_var;

    FILE *gain_ptr = fopen(fileName, "r");

    if (gain_ptr == NULL) {
        printf("Error Reading File\n");
        return 1; //return with failure
    }

    for (int row = 0; row < arr_row; row++) {
        fgets(buffer, 1024, gain_ptr); //buffer is declared as char buffer[1024];

        tok = strtok(buffer, "\t");    // NULL means 'continue from last token'
        //tok is declared as char *tok;
        //the numbers in my text file are seperated by tabs
        //so the tokens should be separated by \t

        column = 0;  //since fgets gets a whole line I had to separate the columns
        while (tok != NULL) {
            test_var = atof(tok);  //test variable is declared as double
            //and it is used to hold temporarily the value of tok in double
            tok = strtok(NULL, "\t");    // NULL means 'continue from last token'

            mat[row][column] = test_var;
            column++;
        }
    }
    fclose(gain_ptr);
    return 1;
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

int *create_empty_vector(size_t len) {
    int *vector;
    vector = (int *) calloc(len, sizeof(int));
    return vector;
}

int *findUnmatchedSwitches(int *curTT, int *baseTT) {
    int *arr = create_empty_vector(7);
    int S = 0;
    int i;
    for (i = 0; i < 6; ++i) {
        if (curTT[i] != baseTT[i]) {
            S++;
            arr[i] = 1;
        } else arr[i] = 0;
    }
    arr[i] = S;
    return arr;
}

double **calcDinv(double **AbaseInv, double **u, double **v) {
    double **unity = create_empty_matrix_double(2, 2);
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            if (i == j) {
                unity[i][j] = 1;
            }
        }
    }
//    D = eye(2)+v*A0Inv*u;
//    Dinv = inv(D);
    return gauss_jordan_matrix_inverse(add_matrix(unity,mul_matrix_matrix(mul_matrix_matrix(v, AbaseInv, 2, 16, 16, 16), u, 2, 16, 16, 2),2,2,2,2),2);
}

double **
smw(int ts, int **TT, double **Abase, double **AbaseInv, int numOfNodes, int **switchNodeMat, int gon) {

    double **u = create_empty_matrix_double(16, 2);
    double **v = create_empty_matrix_double(2, 16);
    double **A0;
    double **A0Inv;
    double **A1;
    double **A1Inv;

    A0 = Abase;
    A0Inv = AbaseInv;

    int nodeCount;
//    i represents switches
    for (int i = 0; i < 6; i++) {
        if (TT[i][ts] != TT[i][0]) {
            nodeCount = 0;
//            j represents nodes
            for (int j = 0; j < 5; ++j) {
                if (switchNodeMat[i][j] == 1) {
                    if (nodeCount == 0) {
                        u[j][0] = gon;
                        u[j][1] = -gon;
                        v[0][j] = 1;
                        nodeCount++;
                    } else {
                        u[j][0] = -gon;
                        u[j][1] = gon;
                        v[1][j] = 1;
                    }
                }
            }

//            A1 = A0 +u*v;
            A1 = add_matrix(A0, mul_matrix_matrix(u, v, 16, 2, 2, 16), 16, 16, 16, 16);
            /*update Ainv
            A1Inv = A0Inv - A0Inv * u * Dinv * v * A0Inv;*/
            A1Inv = sub_matrix(A0Inv, mul_matrix_matrix(mul_matrix_matrix(
                    mul_matrix_matrix(mul_matrix_matrix(A0Inv, u, 16, 16, 16, 2), calcDinv(A0Inv, u, v), 16, 2, 2, 2), v, 16, 2, 2,
                    16), A0Inv, 16, 16, 16, 16), 16, 16, 16, 16);
            A0 = A1;
            A0Inv = A1Inv;
        }
    }
    return A0Inv;
}


void main() {
    int size = 16;
    double **Abase = create_empty_matrix_double(size, size);
    readFileDouble(Abase, size, "../HH.txt");

    double **AbaseInv = create_empty_matrix_double(size, size);
    readFileDouble(AbaseInv, size, "../IH.txt");



    int **switchNodeMat = create_empty_matrix_int(6, 5);
    readFileInt(switchNodeMat, 6, "../switchNodeMat.txt");

    int **TT = create_empty_matrix_int(6, 64);
    readFileInt(TT, 6, "../TT.txt");

    double **Ainv = smw(3, TT, Abase, AbaseInv, 16, switchNodeMat, 1000);
    show_matrix_double(Ainv, size,size);

}



