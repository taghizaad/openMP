#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

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

double *mul_matrix_matrix(int row1, int col1, int row2, int col2, double *mat1, double *mat2, double *result) {
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col2; j++)
            for (int k = 0; k < col1; k++)
                result[i * col2 + j] += mat1[i * col1 + k] * mat2[k * col2 + j];
    return result;
}

void show_matrix(double *matrix, size_t row, size_t col) {
    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < col; j++)
            printf("%f ", matrix[i * col + j]);
        puts("");
    }
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
        if (mat[0] != 0) {
            mat[0] = 1 / mat[0];
        } else fprintf(stderr, "division by ZERO error!");
    }
    if (dim == 2) {
        double invDet = 1 / (mat[0] * mat[3] - mat[1] * mat[2]);
        double a = mat[0];
        mat[0] = mat[3] * invDet;
        mat[1] = -mat[1] * invDet;
        mat[2] = -mat[2] * invDet;
        mat[3] = a * invDet;
    }
}

int *chooseBase(int numOfSwitches, int numOfTimeSteps, int ts, double *TT) {
    int *arrBase1 = (int *) malloc(sizeof(int) * (numOfSwitches + 1));
    int *arrBase2 = (int *) malloc(sizeof(int) * (numOfSwitches + 1));
    int base1Score = 0, base2Score = 0, i;
    for (i = 0; i < numOfSwitches; ++i) {
        if (TT[i * numOfTimeSteps + ts] != TT[i * numOfTimeSteps + 0]) {
            base1Score++;
            arrBase1[i] = i;
        } else {
            arrBase1[i] = -1;
        }
        if (TT[i * numOfTimeSteps + ts] != TT[i * numOfTimeSteps + numOfTimeSteps - 1]) {
            base2Score++;
            arrBase2[i] = i;
        } else {
            arrBase2[i] = -1;
        }
    }
    arrBase1[i] = 0;
    arrBase2[i] = numOfTimeSteps - 1;

    if (base1Score <= base2Score)
        return arrBase1;
    else
        return arrBase2;
}

void
fillAin(int ts, int rowAbaseInv, int colAbaseInv, int rowAinv, int colAinv,
        double AbaseInv[rowAbaseInv][colAbaseInv],
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

void smwModified(int ts, int numOfNodes, int numOfSwitches, int numOfTimeSteps, int gon, int *baseArr,
                 double *TT, double *Abase, double *AbaseInv, double *switchNodeMat) {

    double *u, *v, *uv, *vAbaseInv, *vAbaseInvU, *D, *DinvVAbaseInv, *AbaseInvU, *AbaseInvUDinvVAbaseInv;

    int *affectedNodes = (int *) malloc(sizeof(int) * (numOfNodes + 1));

    int nodeCount;


    //iterate over switches at time step ts
    for (int i = 0; i < numOfSwitches; ++i) {

        // Take the effects of i-th switch into account
        if (baseArr[i] != -1) { //if i-th switch status changed compare to base switch pattern
            memset(affectedNodes, 0, sizeof(int) * (numOfNodes + 1));
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
                        u[m] = (TT[i * numOfTimeSteps + ts] - TT[i * numOfTimeSteps + baseArr[numOfSwitches]]) *
                               gon;
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
                            u[n * 2] =
                                    (TT[i * numOfTimeSteps + ts] -
                                     TT[i * numOfTimeSteps + baseArr[numOfSwitches]]) * gon;
                            u[n * 2 + 1] =
                                    (TT[i * numOfTimeSteps + ts] -
                                     TT[i * numOfTimeSteps + baseArr[numOfSwitches]]) * (-gon);
                            v[n] = 1;
                            nodeCount++;
                        } else {
                            u[n * 2] =
                                    (TT[i * numOfTimeSteps + ts] -
                                     TT[i * numOfTimeSteps + baseArr[numOfSwitches]]) * (-gon);
                            u[n * 2 + 1] =
                                    (TT[i * numOfTimeSteps + ts] -
                                     TT[i * numOfTimeSteps + baseArr[numOfSwitches]]) * gon;
                            v[n + numOfNodes] = 1;
                        }
                    }
                }
            }

//            show_matrix(u, numOfNodes, affectedNodes[numOfNodes]);
//            show_matrix(v, affectedNodes[numOfNodes], numOfNodes);
            //uv = u*v
            uv = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
            mul_matrix_matrix(numOfNodes, affectedNodes[numOfNodes], affectedNodes[numOfNodes], numOfNodes, u, v, uv);


            //Abase <-- Abase + uv
            add_matrix(numOfNodes, numOfNodes, numOfNodes, numOfNodes, Abase, uv);

            //vAbaseInv <-- v * AbaseInv
            vAbaseInv = (double *) malloc(sizeof(double) * affectedNodes[numOfNodes] * numOfNodes);

            //fill vAbaseInv without matrix multiplication
            int firstRowFlag = 0;
            int secondRowFlag = 0;
            for (int k = 0; k < numOfNodes; ++k) {
                if (firstRowFlag == 1 && secondRowFlag == 1) {
                    break;
                }
                if (v[k] == 1) {
                    for (int j = 0; j < numOfNodes; ++j) {
                        vAbaseInv[j] = AbaseInv[k * numOfNodes + j];
                    }
                    if (affectedNodes[numOfNodes] == 1)
                        break;
                    firstRowFlag = 1;
                }
                if (affectedNodes[numOfNodes] == 1)
                    continue;
                if (v[numOfNodes + k] == 1) {
                    for (int j = 0; j < numOfNodes; ++j) {
                        vAbaseInv[numOfNodes + j] = AbaseInv[k * numOfNodes + j];
                    }
                    secondRowFlag = 1;
                }
            }

            //D <-- I
            D = (double *) malloc(sizeof(double) * affectedNodes[numOfNodes] * affectedNodes[numOfNodes]);
            identityMatrix(affectedNodes[numOfNodes], D);
            // vAbaseInvU <-- vAbaseInv * u
            vAbaseInvU = (double *) malloc(sizeof(double) * affectedNodes[numOfNodes] * affectedNodes[numOfNodes]);
            mul_matrix_matrix(affectedNodes[numOfNodes], numOfNodes, numOfNodes, affectedNodes[numOfNodes], vAbaseInv,
                              u, vAbaseInvU);

            //D <-- I + v * AbaseInv * u ==== I + vAbaseInvU
            add_matrix(affectedNodes[numOfNodes], affectedNodes[numOfNodes], affectedNodes[numOfNodes],
                       affectedNodes[numOfNodes], D, vAbaseInvU);

            //D <-- Dinv
            invertMatrix(D, affectedNodes[numOfNodes]);

            //DinvVAbaseInv <-- Dinv * vAbaseInv
            DinvVAbaseInv = (double *) malloc(sizeof(double) * affectedNodes[numOfNodes] * numOfNodes);
            mul_matrix_matrix(affectedNodes[numOfNodes], affectedNodes[numOfNodes], affectedNodes[numOfNodes],
                              numOfNodes, D, vAbaseInv, DinvVAbaseInv);

            //AbaseInvU <-- AbaseInv * u
            AbaseInvU = (double *) malloc(sizeof(double) * numOfNodes * affectedNodes[numOfNodes]);
            mul_matrix_matrix(numOfNodes, numOfNodes, numOfNodes, affectedNodes[numOfNodes], AbaseInv, u, AbaseInvU);

            //AbaseInvUDinvVAbaseInv <-- AbaseInvU * DinvVAbaseInv
            AbaseInvUDinvVAbaseInv = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
            mul_matrix_matrix(numOfNodes, affectedNodes[numOfNodes], affectedNodes[numOfNodes], numOfNodes, AbaseInvU,
                              DinvVAbaseInv, AbaseInvUDinvVAbaseInv);


            //AbaseInv <-- AbaseInv - AbaseInv * u * Dinv * v * AbaseInv === AbaseInv - AbaseInvUDinvVAbaseInv
            sub_matrix(numOfNodes, numOfNodes, numOfNodes, numOfNodes, AbaseInv, AbaseInvUDinvVAbaseInv);


            // Deallocate allocated memory blocks
            free(u);
            free(v);
            free(uv);
            free(vAbaseInv);
            free(D);
            free(vAbaseInvU);
            free(DinvVAbaseInv);
            free(AbaseInvU);
            free(AbaseInvUDinvVAbaseInv);
        }
    }
}

void main() {
    int numOfNodes = 21, numOfSwitches = 18, numOfTimeSteps = 3, gon = 1000, *baseArr;
    double *Abase1, *Abase1_0, *Abase1Inv, *Abase1Inv_0, *Abase2, *Abase2_0,
            *Abase2Inv, *Abase2Inv_0, *switchNodeMat, *TT, *Ainv;

    Abase1 = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    Abase1_0 = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    Abase1Inv = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    Abase1Inv_0 = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    Abase2 = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    Abase2_0 = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    Abase2Inv = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    Abase2Inv_0 = (double *) malloc(sizeof(double) * numOfNodes * numOfNodes);
    switchNodeMat = (double *) malloc(sizeof(double) * numOfSwitches * numOfNodes);
    TT = (double *) malloc(sizeof(double) * numOfSwitches * numOfTimeSteps);
    Ainv = (double *) malloc(sizeof(double) * numOfNodes * numOfTimeSteps * numOfNodes);

    setbuf(stdout, NULL);

    readFile(Abase1, "../smw18switch/Abase1.txt");
    readFile(Abase1_0, "../smw18switch/Abase1.txt");

    readFile(Abase1Inv, "../smw18switch/Abase1Inv.txt");
    readFile(Abase1Inv_0, "../smw18switch/Abase1Inv.txt");

    readFile(Abase2, "../smw18switch/Abase2.txt");
    readFile(Abase2_0, "../smw18switch/Abase2.txt");

    readFile(Abase2Inv, "../smw18switch/Abase2Inv.txt");
    readFile(Abase2Inv_0, "../smw18switch/Abase2Inv.txt");

    readFile(switchNodeMat, "../smw18switch/switchNodeMat.txt");

    readFile(TT, "../smw18switch/TT.txt");

//    fillAin(0, numOfNodes, numOfNodes, numOfNodes * numOfTimeSteps, numOfNodes, Abase1Inv, Ainv);
//    fillAin(numOfTimeSteps-1, numOfNodes, numOfNodes, numOfNodes * numOfTimeSteps, numOfNodes, Abase2Inv, Ainv);


    int niter = 1;
    double time1, timedif;
    time1 = (double) clock() / CLOCKS_PER_SEC;

//    for (int i = 0; i < niter; ++i) {
//        for (int j = numOfTimeSteps - 1; j < numOfTimeSteps; ++j) {
    int j = 1;
    //baseArr size is (numOfSwitches + 1) in which the last element shows which base selected
    baseArr = chooseBase(numOfSwitches, numOfTimeSteps, j, TT);
    if (baseArr[numOfSwitches] == 0) {
        smwModified(j, numOfNodes, numOfSwitches, numOfTimeSteps, gon, baseArr, TT,
                    Abase1, Abase1Inv, switchNodeMat);
        show_matrix(Abase1Inv, numOfNodes, numOfNodes);
//                fillAin(j, numOfNodes, numOfNodes, numOfNodes * numOfTimeSteps, numOfNodes, Abase1Inv, Ainv);
//                renewAbaseAndAbaseInv(numOfNodes, numOfNodes, Abase1, Abase1_0, Abase1Inv, Abase1Inv_0);
    } else {
        smwModified(j, numOfNodes, numOfSwitches, numOfTimeSteps, gon, baseArr, TT,
                    Abase2, Abase2Inv, switchNodeMat);
        show_matrix(Abase2Inv, numOfNodes, numOfNodes);
//                fillAin(j, numOfNodes, numOfNodes, numOfNodes * numOfTimeSteps, numOfNodes, Abase2Inv, Ainv);
//                renewAbaseAndAbaseInv(numOfNodes, numOfNodes, Abase2, Abase2_0, Abase2Inv, Abase2Inv_0);
    }
//        }
//    }


//    timedif = (((double) clock()) / CLOCKS_PER_SEC) - time1;

//    printf("The elapsed time for %d time-steps is %f Microseconds\n", numOfTimeSteps, timedif / niter * 1000000);
//    printf("The elapsed time for each time-steps is %f Microseconds\n", timedif / niter / numOfTimeSteps * 1000000);

}



