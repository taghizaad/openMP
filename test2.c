
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

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

void show_matrix_double(int row, int col, double matrix[row][col]) {
    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < col; j++)
            printf("%f ", matrix[i][j]);
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

void matrixInversion(double mat[2][2]) {
    double det = 1 / (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);
    int a = mat[0][0];
    mat[0][0] = mat[1][1] * det;
    mat[0][1] = -mat[0][1] * det;
    mat[1][0] = -mat[1][0] * det;
    mat[1][1] = a * det;
}


void smw(int kthTimeStep, int numOfNodes, int numOfSwitches, int numOfTimeSteps, int numOfChangedNodes, int gon,
         double TT[numOfSwitches][numOfTimeSteps], double Abase[numOfNodes][numOfNodes],
         double AbaseInv[numOfNodes][numOfNodes], double switchNodeMat[numOfSwitches][numOfChangedNodes]) {

    double u[numOfNodes][2];
    double v[2][numOfNodes];
    double uv[numOfNodes][numOfNodes];
    double AbaseInvU[numOfNodes][2];
    double D[2][2];
    double vAbaseInvU[2][2];
    double AbaseInvUDinv[numOfNodes][2];
    double AbaseInvUDinvV[numOfNodes][numOfNodes];
    double AbaseInvUDinvVAbaseInv[numOfNodes][numOfNodes];


    int nodeCount;
//    i represents switches
    for (int i = 0; i < 6; i++) {
        if (TT[i][kthTimeStep] != TT[i][0]) {
            nodeCount = 0;
            memset(u, 0, sizeof(u));
            memset(v, 0, sizeof(v));
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
            matrixInversion(D);
            //AbaseInvUDinv <-- AbaseInv * u * Dinv
            memset(AbaseInvUDinv, 0, sizeof(AbaseInvUDinv));
            mul_matrix_matrix(numOfNodes, 2, 2, 2, AbaseInvU, D, AbaseInvUDinv);
            //AbaseInvUDinvV <-- AbaseInv * u * Dinv * v
            memset(AbaseInvUDinvV, 0, sizeof(AbaseInvUDinvV));
            mul_matrix_matrix(numOfNodes, 2, 2, numOfNodes, AbaseInvUDinv, v, AbaseInvUDinvV);
            //AbaseInvUDinvVAbaseInv <--AbaseInv * u * Dinv * v * AbaseInv
            memset(AbaseInvUDinvVAbaseInv, 0, sizeof(AbaseInvUDinvVAbaseInv));
            mul_matrix_matrix(numOfNodes, numOfNodes, numOfNodes, numOfNodes, AbaseInvUDinvV, AbaseInv,
                              AbaseInvUDinvVAbaseInv);
            //AbaseInv <-- AbaseInv - AbaseInv * u * Dinv * v * AbaseInv
            sub_matrix(numOfNodes, numOfNodes, numOfNodes, numOfNodes, AbaseInv, AbaseInvUDinvVAbaseInv);
        }
    }
}


void main() {
    int numOfNodes = 16, numOfSwitches = 6, numOfChangedNodes = 5, numOfTimeSteps = 64, gon = 1000, kthTimeStep = 63;
    double Abase[numOfNodes][numOfNodes], AbaseInv[numOfNodes][numOfNodes],
            switchNodeMat[numOfSwitches][numOfChangedNodes], TT[numOfSwitches][numOfTimeSteps];
    readFileDouble(numOfNodes, numOfNodes, Abase, "../HH.txt");

    readFileDouble(numOfNodes, numOfNodes, AbaseInv, "../IH.txt");

    readFileDouble(numOfSwitches, numOfChangedNodes, switchNodeMat, "../switchNodeMat.txt");

    readFileDouble(numOfSwitches, numOfTimeSteps, TT, "../TT.txt");

    smw(kthTimeStep, numOfNodes, numOfSwitches, numOfTimeSteps, numOfChangedNodes, gon, TT, Abase, AbaseInv,
        switchNodeMat);
    show_matrix_double(numOfNodes, numOfNodes, AbaseInv);

}



