#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

void sub_matrix(int row1, int col1, int row2, int col2, double mat1[row1][col1], double mat2[row2][col2]) {
    if (row1 != row2 || col1 != col2) {
        fprintf(stderr, "incompatible matrix sizes for subtraction\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col1; j++)
            mat1[i][j] = mat1[i][j] - mat2[i][j];
}

void add_matrix(int row1, int col1, int row2, int col2, double mat1[row1][col1], double mat2[row2][col2]) {
    if (row1 != row2 || col1 != col2) {
        fprintf(stderr, "incompatible matrix sizes for subtraction\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col1; j++)
            mat1[i][j] = mat1[i][j] + mat2[i][j];
}

void mul_matrix_matrix(int row1, int col1, int row2, int col2, double mat1[row1][col1], double mat2[row2][col2],
                       double mat1mat2[row1][col2]) {
    if (col1 != row2) {
        fprintf(stderr, "incompatible matrix sizes for multiplication\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col2; j++)
            for (int k = 0; k < col1; k++)
                mat1mat2[i][j] += mat1[i][k] * mat2[k][j];
}

void show_matrix(int row, int col, double matrix[row][col]) {
    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < col; j++)
            printf("%f ", matrix[i][j]);
        puts("");
    }
}

void show_vector(double *matrix, size_t len) {
    for (size_t i = 0; i < len; i++) {
        printf("%f ", matrix[i]);
        puts("");
    }
}

int readFileDouble(int row, int col, double mat[row][col], char *fileName) {
    char buffer[1024];
    char *tok;
    double test_var;
    int column;

    FILE *gain_ptr = fopen(fileName, "r");

    if (gain_ptr == NULL) {
        printf("Error Reading File\n");
        return 1; //return with failure
    }

    for (int i = 0; i < row; i++) {
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

            mat[i][column] = test_var;
            column++;
        }
    }
    fclose(gain_ptr);
    return 1;
}

void identityMatrix(int n, double mat[n][n]) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                mat[i][j] = 1;
            } else mat[i][j] = 0;
        }

    }
}

void invertMatrix(double mat[2][2]) {
    double invDet = 1 / (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);
    double a = mat[0][0];
    mat[0][0] = mat[1][1] * invDet;
    mat[0][1] = -mat[0][1] * invDet;
    mat[1][0] = -mat[1][0] * invDet;
    mat[1][1] = a * invDet;
}

void smwModified(int ts, int numOfNodes, int numOfSwitches, int numOfTimeSteps, int numOfChangedNodes, int gon,
                 double TT[numOfSwitches][numOfTimeSteps], double Abase[numOfNodes][numOfNodes],
                 double AbaseInv[numOfNodes][numOfNodes],
                 double switchNodeMat[numOfSwitches][numOfChangedNodes]) {

    double u[numOfNodes][2];
    double v[2][numOfNodes];
    double uv[numOfNodes][numOfNodes];
    double AbaseInvU[numOfNodes][2];
    double D[2][2];
    double vAbaseInvU[2][2];
    double AbaseInvUDinvVAbaseInv[numOfNodes][numOfNodes];
    double vAbaseInv[2][numOfNodes];
    double DinvVAbaseInv[2][numOfNodes];

    int nodeCount, i, j, k;

//    i represents switches
    for (i = 0; i < 6; i++) {
        if (TT[i][ts] != TT[i][0]) {
            nodeCount = 0;
            memset(u, 0, sizeof(u));
            memset(v, 0, sizeof(v));
//            j represents nodes
            for (j = 0; j < 5; ++j) {
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

            //uv = u*v
            memset(uv, 0, sizeof(uv));
            mul_matrix_matrix(numOfNodes, 2, 2, numOfNodes, u, v, uv);

            //Abase <-- Abase + uv
            add_matrix(numOfNodes, numOfNodes, numOfNodes, numOfNodes, Abase, uv);
            /*update Ainv
            A1Inv = AbaseInv - AbaseInv * u * Dinv * v * AbaseInv;
             */

            //AbaseInvU <-- AbaseInv * u
            memset(AbaseInvU, 0, sizeof(AbaseInvU));
            mul_matrix_matrix(numOfNodes, numOfNodes, numOfNodes, 2, AbaseInv, u, AbaseInvU);

            //D <-- eye(2)
            identityMatrix(2, D);
            // vAbaseInvU <-- v * AbaseInvU
            memset(vAbaseInvU, 0, sizeof(vAbaseInvU));
            mul_matrix_matrix(2, numOfNodes, numOfNodes, 2, v, AbaseInvU, vAbaseInvU);
            //D <-- eye(2) + v * AbaseInv * u ~ eye(2) + v * AbaseInvU
            add_matrix(2, 2, 2, 2, D, vAbaseInvU);
            //D <-- Dinv
            invertMatrix(D);


            //vAbaseInv <-- v * AbaseInv
            int firstRowFlag = 0;
            int secondRowFlag = 0;
            for (k = 0; k < numOfNodes; ++k) {
                if (firstRowFlag == 1 && secondRowFlag == 1) {
                    break;
                }
                if (v[0][k] == 1) {
                    for (int j = 0; j < numOfNodes; ++j) {
                        vAbaseInv[0][j] = AbaseInv[k][j];
                    }
                    firstRowFlag = 1;
                }
                if (v[1][k] == 1) {
                    for (int j = 0; j < numOfNodes; ++j) {
                        vAbaseInv[1][j] = AbaseInv[k][j];
                    }
                    secondRowFlag = 1;
                }
            }

            //DinvVAbaseInv <-- Dinv * vAbaseInv
            memset(DinvVAbaseInv, 0, sizeof(DinvVAbaseInv));
            mul_matrix_matrix(2, 2, 2, numOfNodes, D, vAbaseInv, DinvVAbaseInv);

            //AbaseInvUDinvVAbaseInv <-- AbaseInvU * DinvVAbaseInv
            memset(AbaseInvUDinvVAbaseInv, 0, sizeof(AbaseInvUDinvVAbaseInv));
            mul_matrix_matrix(numOfNodes, 2, 2, numOfNodes, AbaseInvU, DinvVAbaseInv, AbaseInvUDinvVAbaseInv);


            //AbaseInv <-- AbaseInv - AbaseInv * u * Dinv * v * AbaseInv
            sub_matrix(numOfNodes, numOfNodes, numOfNodes, numOfNodes, AbaseInv, AbaseInvUDinvVAbaseInv);

        }
    }
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
    int i, j;
    for (i = 0; i < row; ++i) {
        for (j = 0; j < col; ++j) {
            Abase[i][j] = Abase0[i][j];
            AbaseInv[i][j] = AbaseInv0[i][j];
        }
    }
}

void main() {

    int numOfNodes = 16, numOfSwitches = 6, numOfChangedNodes = 5, numOfTimeSteps = 64, gon = 1000;
    double Abase[numOfNodes][numOfNodes], Abase0[numOfNodes][numOfNodes],
            AbaseInv[numOfNodes][numOfNodes], AbaseInv0[numOfNodes][numOfNodes],
            switchNodeMat[numOfSwitches][numOfChangedNodes], TT[numOfSwitches][numOfTimeSteps],
            Ainv[numOfNodes * numOfTimeSteps][numOfNodes], timeConsumption[numOfTimeSteps];

    readFileDouble(numOfNodes, numOfNodes, Abase, "../Abase1.txt");
    readFileDouble(numOfNodes, numOfNodes, Abase0, "../Abase1.txt");

    readFileDouble(numOfNodes, numOfNodes, AbaseInv, "../Abase1Inv.txt");
    readFileDouble(numOfNodes, numOfNodes, AbaseInv0, "../Abase1Inv.txt");

    readFileDouble(numOfSwitches, numOfChangedNodes, switchNodeMat, "../switchNodeMat.txt");

    readFileDouble(numOfSwitches, numOfTimeSteps, TT, "../TT.txt");

    fillAin(0, numOfNodes, numOfNodes, numOfNodes * numOfTimeSteps, numOfNodes, AbaseInv, Ainv);

    int niter = 100;
    double time1, timedif;
    time1 = (double) clock() / CLOCKS_PER_SEC;

    int i, j;
    for (i = 0; i < niter; ++i) {
        for (j = 1; j < numOfTimeSteps; ++j) {
            smwModified(j, numOfNodes, numOfSwitches, numOfTimeSteps, numOfChangedNodes, gon, TT, Abase, AbaseInv,
                        switchNodeMat);
            fillAin(j, numOfNodes, numOfNodes, numOfNodes * numOfTimeSteps, numOfNodes, AbaseInv, Ainv);
            renewAbaseAndAbaseInv(numOfNodes, numOfNodes, Abase, Abase0, AbaseInv, AbaseInv0);
        }
    }

//    show_matrix(numOfNodes * numOfTimeSteps, numOfNodes, Ainv);

    timedif = (((double) clock()) / CLOCKS_PER_SEC) - time1;

    printf("The elapsed time for %d time-steps is %f Microseconds\n", numOfTimeSteps, timedif / niter * 1000000);
    printf("The elapsed time for each time-steps is %f Microseconds\n", timedif / niter / numOfTimeSteps * 1000000);

}



