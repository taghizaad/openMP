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

void sub_matrix(int row1, int col1, int row2, int col2, double *mat1, double *mat2) {
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col1; j++)
            mat1[i * col1 + j] = mat1[i * col1 + j] - mat2[i * col1 + j];
}

void add_matrix(int row1, int col1, int row2, int col2, double *mat1, double *mat2) {
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col1; j++)
            mat1[i * col1 + j] = mat1[i * col1 + j] + mat2[i * col1 + j];
}

void mul_matrix_matrix(int row1, int col1, int row2, int col2, double *mat1, double *mat2, double *result) {
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col2; j++)
            for (int k = 0; k < col1; k++)
                result[i * col2 + j] = result[i * col2 + j] + mat1[i * col1 + k] * mat2[k * col2 + j];
}

void show_matrix(double *matrix, size_t row, size_t col) {
    printf("sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss\n");
    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < col; j++) {
//            printf("(%d,%d)%f ", i, j, matrix[i * col + j]);
            printf("%15.15e ", matrix[i * col + j]);
        }
        printf("\n");
    }
    printf("eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee\n");
}

void show_vector(double *vector, size_t len) {
    for (size_t i = 0; i < len; i++) {
        printf("%i --> %f ", i, vector[i]);
        puts("");
    }
}

void show_vector_int(int *vector, size_t len) {
    for (size_t i = 0; i < len; i++) {
        printf("%i --> %d ", i, vector[i]);
        puts("");
    }
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

void identityMatrix(int n, double *mat) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                mat[i * n + j] = 1;
            } else mat[i * n + j] = 0;
        }

    }
}

void invertMatrix(double *mat, int dim) {
    if (dim == 1) {
        if (fabs(mat[0]) > 0.000001) {
            mat[0] = 1 / mat[0];
//            show_matrix(mat, 1, 1);
        } else {
            fprintf(stderr, "division by ZERO error!");
            exit(-1);
        }
    }
    if (dim == 2) {
        double det = mat[0] * mat[3] - mat[1] * mat[2];
//        show_matrix(&det, 1, 1);
        if (fabs(det) > 0.000001) {
            double invDet = 1 / det;
            double a = mat[0];
            mat[0] = mat[3] * invDet;
            mat[1] = -mat[1] * invDet;
            mat[2] = -mat[2] * invDet;
            mat[3] = a * invDet;
        } else {
            fprintf(stderr, "division by ZERO error!");
            exit(-1);
        }
    }
}

double *chooseBase(int numOfSwitches, double *curTs, double *sigmaBase1, double *sigmaBase2) {
    int base1Score = 0, base2Score = 0;
    for (int i = 0; i < numOfSwitches; ++i) {
        if (curTs[i] == sigmaBase1[i]) {
            base1Score++;
        }
        if (curTs[i] == sigmaBase2[i]) {
            base2Score++;
        }
    }

    if (base1Score >= base2Score)
        return sigmaBase1;
    else
        return sigmaBase2;
}


void
fillAin(int ts, int rowAbaseInv, int colAbaseInv, int rowAinv, int colAinv, double AbaseInv[rowAbaseInv][colAbaseInv],
        double Ainv[rowAinv][colAinv]) {
    for (int i = 0; i < rowAbaseInv; ++i) {
        for (int j = 0; j < colAbaseInv; ++j) {
            Ainv[ts * rowAbaseInv + i][j] = AbaseInv[i][j];
        }
    }
}

void renewAbaseAndAbaseInv(int row, int col, double Abase[row][col], double Abase0[row][col],
                           double AbaseInv[row][col],
                           double AbaseInv0[row][col]) {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            Abase[i][j] = Abase0[i][j];
            AbaseInv[i][j] = AbaseInv0[i][j];
        }
    }
}

void smwModified(int numOfNodes, int numOfSwitches, int gon, double *baseArr, double *curTs, double *AbaseInv,
                 double *switchNodeMat) {

    double *u, *v, *uv, *vAbaseInv, *vAbaseInvU, *D, *DinvVAbaseInv, *AbaseInvU, *AbaseInvUDinvVAbaseInv;
    int *affectedNodes, nodeCount;


    //iterate over switches at time step ts
    for (int i = 0; i < numOfSwitches; ++i) {

        // Take the effects of i-th switch into account
        if (curTs[i] != baseArr[i]) { //if i-th switch status changed compare to base switch pattern
            affectedNodes = (int *) calloc((numOfNodes + 1), sizeof(int));
            //iterate over nodes in switchNodeMat
            for (int j = 0; j < numOfNodes; ++j) {
                //if j-th node is affected by i-th switch
                if (switchNodeMat[i * numOfNodes + j] == 1) {
                    affectedNodes[j] = j;
                    affectedNodes[numOfNodes] = affectedNodes[numOfNodes] + 1;
                } else
                    affectedNodes[j] = -1;
            }
            u = (double *) calloc(numOfNodes * affectedNodes[numOfNodes], sizeof(double));
            v = (double *) calloc(numOfNodes * affectedNodes[numOfNodes], sizeof(double));

            if (affectedNodes[numOfNodes] == 1) {
                //iterate over nodes in switchNodeMat
                for (int m = 0; m < numOfNodes; ++m) {
                    if (affectedNodes[m] != -1) {
                        u[m] = (curTs[i] - baseArr[i]) * gon;
                        v[m] = 1;
                        break;
                    }
                }
            }
            if (affectedNodes[numOfNodes] == 2) {
                nodeCount = 0;
                for (int n = 0; n < numOfNodes; ++n) {
                    if (affectedNodes[n] != -1) {
                        if (nodeCount == 0) {
                            u[n * 2] = (curTs[i] - baseArr[i]) * gon;
                            u[n * 2 + 1] = (curTs[i] - baseArr[i]) * (-gon);
                            v[n] = 1;
                            nodeCount++;
                        } else {
                            u[n * 2] = (curTs[i] - baseArr[i]) * (-gon);
                            u[n * 2 + 1] = (curTs[i] - baseArr[i]) * gon;
                            v[n + numOfNodes] = 1;
                        }
                    }
                }
            }

//            show_matrix(u, numOfNodes, affectedNodes[numOfNodes]);
//            show_matrix(v, affectedNodes[numOfNodes], numOfNodes);
//            uv = (double *) calloc(numOfNodes * numOfNodes, sizeof(double));
//            mul_matrix_matrix(numOfNodes, affectedNodes[numOfNodes], affectedNodes[numOfNodes], numOfNodes, u, v, uv);
//            show_matrix(uv, numOfNodes, numOfNodes);

//            show_matrix(AbaseInv, numOfNodes, numOfNodes);

            //vAbaseInv <-- v * AbaseInv
            vAbaseInv = (double *) malloc(sizeof(double) * affectedNodes[numOfNodes] * numOfNodes);
            mul_matrix_matrix(affectedNodes[numOfNodes], numOfNodes, numOfNodes, numOfNodes, v, AbaseInv, vAbaseInv);

            // vAbaseInvU <-- vAbaseInv * u
            vAbaseInvU = (double *) malloc(sizeof(double) * affectedNodes[numOfNodes] * affectedNodes[numOfNodes]);
            mul_matrix_matrix(affectedNodes[numOfNodes], numOfNodes, numOfNodes, affectedNodes[numOfNodes], vAbaseInv,
                              u, vAbaseInvU);

//            show_matrix(vAbaseInvU, affectedNodes[numOfNodes], affectedNodes[numOfNodes]);

//            show_matrix(AbaseInv, numOfNodes, numOfNodes);

            //D <-- I
            D = (double *) malloc(sizeof(double) * affectedNodes[numOfNodes] * affectedNodes[numOfNodes]);
            identityMatrix(affectedNodes[numOfNodes], D);

            //D <-- I + v * AbaseInv * u ==== I + vAbaseInvU
            add_matrix(affectedNodes[numOfNodes], affectedNodes[numOfNodes], affectedNodes[numOfNodes],
                       affectedNodes[numOfNodes], D, vAbaseInvU);

//            show_matrix(D, affectedNodes[numOfNodes], affectedNodes[numOfNodes]);

            //D <-- Dinv
            invertMatrix(D, affectedNodes[numOfNodes]);

//            show_matrix(D, affectedNodes[numOfNodes], affectedNodes[numOfNodes]);

            //DinvVAbaseInv <-- Dinv * vAbaseInv
            DinvVAbaseInv = (double *) malloc(sizeof(double) * affectedNodes[numOfNodes] * numOfNodes);
            mul_matrix_matrix(affectedNodes[numOfNodes], affectedNodes[numOfNodes], affectedNodes[numOfNodes],
                              numOfNodes, D, vAbaseInv, DinvVAbaseInv);

            show_matrix(DinvVAbaseInv,affectedNodes[numOfNodes],numOfNodes);


            //AbaseInvU <-- AbaseInv * u
            AbaseInvU = (double *) malloc(sizeof(double) * numOfNodes * affectedNodes[numOfNodes]);
            mul_matrix_matrix(numOfNodes, numOfNodes, numOfNodes, affectedNodes[numOfNodes], AbaseInv, u, AbaseInvU);

            //AbaseInvUDinvVAbaseInv <-- AbaseInvU * DinvVAbaseInv
            AbaseInvUDinvVAbaseInv = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
            mul_matrix_matrix(numOfNodes, affectedNodes[numOfNodes], affectedNodes[numOfNodes], numOfNodes, AbaseInvU,
                              DinvVAbaseInv, AbaseInvUDinvVAbaseInv);

//            show_matrix(AbaseInvUDinvVAbaseInv, numOfNodes, numOfNodes);

            //AbaseInv <-- AbaseInv - AbaseInv * u * Dinv * v * AbaseInv === AbaseInv - AbaseInvUDinvVAbaseInv
            sub_matrix(numOfNodes, numOfNodes, numOfNodes, numOfNodes, AbaseInv, AbaseInvUDinvVAbaseInv);

            // Deallocate allocated memory blocks
                   free(affectedNodes);
                   free(u);
                   free(v);
//                   free(uv);
                   free(vAbaseInv);
                   free(D);
                   free(vAbaseInvU);
                   free(DinvVAbaseInv);
                   free(AbaseInvU);
                   free(AbaseInvUDinvVAbaseInv);
        }
    }
}

void smwModifiedWithWholeD(int numOfNodes, int numOfSwitches, int gon, double *baseArr, double *curTs,
                           double *AbaseInv, double *switchNodeMat) {

    double *uu, *vv, *uv, *vAbaseInv, *vAbaseInvU, *D, *DinvVAbaseInv, *AbaseInvU, *AbaseInvUDinvVAbaseInv;
    int *affectedNodes, nodeCount;

    uu = (double *) calloc(numOfNodes * numOfNodes, sizeof(double));
    vv = (double *) calloc(numOfNodes * numOfNodes, sizeof(double));

    int *allAffectedNodes = (int *) calloc(numOfNodes, sizeof(int)); //all affected nodes at curTS

    //iterate over switches at time step ts
    for (int i = 0; i < numOfSwitches; ++i) {

        // Take the effects of i-th switch into account
        if (curTs[i] != baseArr[i]) { //if i-th switch status changed compare to base switch pattern
            affectedNodes = (int *) calloc(2, sizeof(int));
            nodeCount = 0;
            //iterate over nodes in switchNodeMat
            for (int j = 0; j < numOfNodes; ++j) {
                //if j-th node is affected by i-th switch
                if (switchNodeMat[i * numOfNodes + j] == 1) {
                    if (nodeCount == 0) {
                        affectedNodes[0] = j;
                        nodeCount++;
                        allAffectedNodes[j] = 1;
                    } else {
                        affectedNodes[1] = j;
                        allAffectedNodes[j] = 1;
                    }
                }
            }

            if (affectedNodes[1] == 0) { // means only 1 node affected by switch i
                uu[affectedNodes[0] * numOfNodes + affectedNodes[0]] += (curTs[i] - baseArr[i]) * gon;
                vv[affectedNodes[0] * numOfNodes + affectedNodes[0]] = 1;
            } else // means 2 nodes affected by switch i
            {
                uu[affectedNodes[0] * numOfNodes + affectedNodes[0]] += (curTs[i] - baseArr[i]) * gon;
                uu[affectedNodes[0] * numOfNodes + affectedNodes[1]] += (curTs[i] - baseArr[i]) * (-gon);
                uu[affectedNodes[1] * numOfNodes + affectedNodes[0]] += (curTs[i] - baseArr[i]) * (-gon);
                uu[affectedNodes[1] * numOfNodes + affectedNodes[1]] += (curTs[i] - baseArr[i]) * gon;

                vv[affectedNodes[0] * numOfNodes + affectedNodes[0]] = 1;
                vv[affectedNodes[1] * numOfNodes + affectedNodes[1]] = 1;
            }
        }
    }

    int numOfAllAffectedNodes = 0;
    for (int h = 0; h < numOfNodes; ++h) {
        numOfAllAffectedNodes += allAffectedNodes[h];
    }
    double *u = (double *) calloc(numOfNodes * numOfAllAffectedNodes, sizeof(double));
    double *v = (double *) calloc(numOfNodes * numOfAllAffectedNodes, sizeof(double));

    int col = 0;
    for (int t = 0; t < numOfNodes; ++t) {
        if (allAffectedNodes[t] == 1) {
            for (int q = 0; q < numOfNodes; ++q) {
                u[q * numOfAllAffectedNodes + col] = uu[q * numOfNodes + t];
                v[col * numOfNodes + q] = vv[q * numOfNodes + t];
//                vv[q * numOfAllAffectedNodes + col] = v[q * numOfNodes + t];
            }
            col++;
        }
    }

//    show_matrix(u, numOfNodes, numOfAllAffectedNodes);
//    show_matrix(v, numOfAllAffectedNodes, numOfNodes);
    uv = (double *) calloc(numOfNodes * numOfNodes, sizeof(double));
    mul_matrix_matrix(numOfNodes, numOfAllAffectedNodes, numOfAllAffectedNodes, numOfNodes, u, v, uv);
//    show_matrix(uv, numOfNodes, numOfNodes);

    //vAbaseInv <-- v * AbaseInv
    vAbaseInv = (double *) malloc(sizeof(double) * numOfAllAffectedNodes * numOfNodes);
    mul_matrix_matrix(numOfAllAffectedNodes, numOfNodes, numOfNodes, numOfNodes, v, AbaseInv, vAbaseInv);

    // vAbaseInvU <-- vAbaseInv * u
    vAbaseInvU = (double *) malloc(sizeof(double) * numOfAllAffectedNodes * numOfAllAffectedNodes);
    mul_matrix_matrix(numOfAllAffectedNodes, numOfNodes, numOfNodes, numOfAllAffectedNodes, vAbaseInv,
                      u, vAbaseInvU);

    //D <-- I
    D = (double *) malloc(sizeof(double) * numOfAllAffectedNodes * numOfAllAffectedNodes);
    identityMatrix(numOfAllAffectedNodes, D);


    //D <-- I + v * AbaseInv * u ==== I + vAbaseInvU
    add_matrix(numOfAllAffectedNodes, numOfAllAffectedNodes, numOfAllAffectedNodes,
               numOfAllAffectedNodes, D, vAbaseInvU);


    double **Dformated = (double **) malloc(sizeof(double) * numOfAllAffectedNodes);
    vecToMat(D, numOfAllAffectedNodes, numOfAllAffectedNodes, Dformated);
    double **DinvFormated = gauss_jordan_matrix_inverse(Dformated, numOfAllAffectedNodes);
    double *Dinv = (double *) malloc(sizeof(double) * numOfAllAffectedNodes * numOfAllAffectedNodes);
    matToVec(DinvFormated, numOfAllAffectedNodes, numOfAllAffectedNodes, Dinv);


    //DinvVAbaseInv <-- Dinv * vAbaseInv
    DinvVAbaseInv = (double *) malloc(sizeof(double) * numOfAllAffectedNodes * numOfNodes);
    mul_matrix_matrix(numOfAllAffectedNodes, numOfAllAffectedNodes, numOfAllAffectedNodes,
                      numOfNodes, Dinv, vAbaseInv, DinvVAbaseInv);

    //AbaseInvU <-- AbaseInv * u
    AbaseInvU = (double *) malloc(sizeof(double) * numOfNodes * numOfAllAffectedNodes);
    mul_matrix_matrix(numOfNodes, numOfNodes, numOfNodes, numOfAllAffectedNodes, AbaseInv, u, AbaseInvU);

    //AbaseInvUDinvVAbaseInv <-- AbaseInvU * DinvVAbaseInv
    AbaseInvUDinvVAbaseInv = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    mul_matrix_matrix(numOfNodes, numOfAllAffectedNodes, numOfAllAffectedNodes, numOfNodes, AbaseInvU,
                      DinvVAbaseInv, AbaseInvUDinvVAbaseInv);

    //AbaseInv <-- AbaseInv - AbaseInv * u * Dinv * v * AbaseInv === AbaseInv - AbaseInvUDinvVAbaseInv
    sub_matrix(numOfNodes, numOfNodes, numOfNodes, numOfNodes, AbaseInv, AbaseInvUDinvVAbaseInv);

    show_matrix(AbaseInv, numOfNodes, numOfNodes);

    // Deallocate allocated memory blocks
/*    free(affectedNodes);
    free(u);
    free(v);
    free(uv);
    free(vAbaseInv);
    free(D);
    free(vAbaseInvU);
    free(DinvVAbaseInv);
    free(AbaseInvU);
    free(AbaseInvUDinvVAbaseInv);*/


}


void main() {
    int numOfNodes = 21, numOfSwitches = 18, gon = 1000;
    double *Abase1, *Abase1_0, *Abase1Inv, *Abase1Inv_0, *sigmaBase1, *Abase2, *Abase2_0,
            *Abase2Inv, *Abase2Inv_0, *sigmaBase2, *switchNodeMat, *TT, *Ainv, *curTs, *baseArr;

    Abase1 = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    Abase1_0 = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    Abase1Inv = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    Abase1Inv_0 = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    sigmaBase1 = (double *) malloc(sizeof(double) * numOfSwitches);
    Abase2 = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    Abase2_0 = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    Abase2Inv = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    Abase2Inv_0 = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    sigmaBase2 = (double *) malloc(sizeof(double) * numOfSwitches);
    switchNodeMat = (double *) malloc(sizeof(double) * numOfSwitches * numOfNodes);
    TT = (double *) malloc(sizeof(double) * numOfSwitches);
//    Ainv = (double *) malloc(sizeof(double) * numOfNodes * numOfTimeSteps * numOfNodes);
    curTs = (double *) malloc(sizeof(double) * numOfSwitches);

    setbuf(stdout, NULL);

    readFile(Abase1, "../smw18switch/Abase1.txt");
    readFile(Abase1_0, "../smw18switch/Abase1.txt");

    readFile(Abase1Inv, "../smw18switch/Abase1Inv.txt");
    readFile(Abase1Inv_0, "../smw18switch/Abase1Inv.txt");

    readFile(sigmaBase1, "../smw18switch/sigmaBase1.txt");

    readFile(Abase2, "../smw18switch/Abase2.txt");
    readFile(Abase2_0, "../smw18switch/Abase2.txt");

    readFile(Abase2Inv, "../smw18switch/Abase2Inv.txt");
    readFile(Abase2Inv_0, "../smw18switch/Abase2Inv.txt");

    readFile(sigmaBase2, "../smw18switch/sigmaBase2.txt");

    readFile(switchNodeMat, "../smw18switch/switchNodeMat.txt");

    readFile(TT, "../smw18switch/TT.txt");

//    fillAin(0, numOfNodes, numOfNodes, numOfNodes * numOfTimeSteps, numOfNodes, Abase1Inv, Ainv);
//    fillAin(numOfTimeSteps-1, numOfNodes, numOfNodes, numOfNodes * numOfTimeSteps, numOfNodes, Abase2Inv, Ainv);


    int niter = 1;
    double time1, timedif;
    time1 = (double) clock() / CLOCKS_PER_SEC;

    for (int s = 0; s < numOfSwitches; ++s) {
        curTs[s] = TT[s];
    }

    //baseArr size is (numOfSwitches + 1) in which the last element shows which base selected
//    baseArr = chooseBase(numOfSwitches, curTs, sigmaBase1, sigmaBase2);
    baseArr = sigmaBase1;
    // == compares the numeric address of pointers and hence determine if they point to the same object.
    if (baseArr == sigmaBase1)
        Ainv = Abase1Inv;
    else
        Ainv = Abase2Inv;

//    show_matrix(Abase1Inv, numOfNodes, numOfNodes);
//    show_matrix(Ainv, numOfNodes, numOfNodes);
    smwModifiedWithWholeD(numOfNodes, numOfSwitches, gon, baseArr, curTs, Ainv, switchNodeMat);
//    show_matrix(Ainv, numOfNodes, numOfNodes);
//                fillAin(j, numOfNodes, numOfNodes, numOfNodes * numOfTimeSteps, numOfNodes, Abase1Inv, Ainv);
//                renewAbaseAndAbaseInv(numOfNodes, numOfNodes, Abase1, Abase1_0, Abase1Inv, Abase1Inv_0);
//        }
//    }


//    timedif = (((double) clock()) / CLOCKS_PER_SEC) - time1;

//    printf("The elapsed time for %d time-steps is %f Microseconds\n", numOfTimeSteps, timedif / niter * 1000000);
//    printf("The elapsed time for each time-steps is %f Microseconds\n", timedif / niter / numOfTimeSteps * 1000000);

}


