
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void vecToMat(double *vec, int row, int col, double **matrix) {
    for (int i = 0; i < row; i++)
        matrix[i] = (double *) malloc(sizeof(double) * col);
    for (size_t i = 0; i < row; i++)
        for (size_t j = 0; j < col; j++)
            matrix[i][j] = vec[i * col + j];
}

void matToVec(double **matrix, int row, int col, double *vec) {
    for (size_t i = 0; i < row; i++)
        for (size_t j = 0; j < col; j++)
            vec[i * col + j] = matrix[i][j];
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
    {                                                       //Unit mat. I give the inverse of matrix A
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
int readFile(double *mat, char *fileName) {
    FILE *fp = fopen(fileName, "r");
    char *tok;
    double test_var;
    char chunk[128];
    // Store the chunks of text into a line buffer
    size_t len = sizeof(chunk);
    char *line = malloc(len);
    if (line == NULL) {
        perror("Unable to allocate memory for the line buffer.");
        exit(1);
    }
    // "Empty" the string
    line[0] = '\0';
    // ...
    int i = 0;
    while (fgets(chunk, sizeof(chunk), fp) != NULL) {
        // Resize the line buffer if necessary
        size_t len_used = strlen(line);
        size_t chunk_used = strlen(chunk);
        if (len - len_used < chunk_used) {
            len *= 2;
            if ((line = realloc(line, len)) == NULL) {
                perror("Unable to reallocate memory for the line buffer.");
                free(line);
                exit(1);
            }
        }
        // Copy the chunk to the end of the line buffer
        strncpy(line + len_used, chunk, len - len_used);
        len_used += chunk_used;
        // Check if line contains '\n', if yes process the line of text
        if (line[len_used - 1] == '\n') {
            tok = strtok(line, "\t");
            while (tok != NULL) {
                test_var = atof(tok);  //test variable is declared as double
                //and it is used to hold temporarily the value of tok in double
                tok = strtok(NULL, "\t");    // NULL means 'continue from last token'
                mat[i] = test_var;
                i++;
            }
//            fputs(line, stdout);
//            fputs("|*\n", stdout);
            // "Empty" the line buffer
            line[0] = '\0';
        }
    }

    fclose(fp);
    free(line);

//    printf("\n\nMax line size: %zd\n", len);
}

void show_matrix(double *matrix, size_t row, size_t col) {
    printf("sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss\n");
    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < col; j++) {
//            printf("(%d,%d)%f ", i, j, matrix[i * col + j]);
            printf("%f ", matrix[i * col + j]);
        }
        printf("\n");
    }
    printf("eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee\n");
}
void main(){
    setbuf(stdout, NULL);
    int numOfAllAffectedNodes = 21;
    double * D = (double *) malloc(sizeof(double) * numOfAllAffectedNodes * numOfAllAffectedNodes);
    readFile(D, "../smw18switch/D.txt");
    double **Dformated = (double **) malloc(sizeof(double) * numOfAllAffectedNodes);
    vecToMat(D, numOfAllAffectedNodes, numOfAllAffectedNodes, Dformated);
    double **DinvFormated = gauss_jordan_matrix_inverse(Dformated, numOfAllAffectedNodes);
    double *Dinv = (double *) malloc(sizeof(double) * numOfAllAffectedNodes * numOfAllAffectedNodes);
    matToVec(DinvFormated,numOfAllAffectedNodes,numOfAllAffectedNodes,Dinv);

    show_matrix(Dinv, numOfAllAffectedNodes, numOfAllAffectedNodes);
}