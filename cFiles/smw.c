#include "matrix.h"
#include "mex.h"
#include <math.h>

void sub_matrix(int row1, int col1, int row2, int col2, double *mat1, double *mat2);

void invertMatrix(double *mat, int dim, int sw, double *outInvArr);

void identityMatrix(int n, double *mat);

void add_matrix(int row1, int col1, int row2, int col2, double *mat1, double *mat2);

void mul_matrix_matrix(int row1, int col1, int row2, int col2, double *mat1, double *mat2, double *result);

void show_matrix(double *matrix, size_t row, size_t col);

void smw(int numOfNodes, int numOfSwitches, int gon, double *selectedBase, double *curTs, double *switchNodeMat,
         double *AbaseInv, double *outInvArr);


void show_matrix(double *matrix, size_t row, size_t col) {
    printf("sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss\n");
    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < col; j++)
            printf("%f ", matrix[i * col + j]);
        printf("\n");
    }
    printf("eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee\n");
}

void sub_matrix(int row1, int col1, int row2, int col2, double *mat1, double *mat2) {
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col1; j++)
            mat1[i * col1 + j] = mat1[i * col1 + j] - mat2[i * col1 + j];
}

void invertMatrix(double *mat, int dim, int sw, double *outInvArr) {
    if (dim == 1) {
        if (fabs(mat[0]) > 0.000001) {
            mat[0] = 1 / mat[0];
        } else {
            outInvArr[0]++;
//            printf("switch %d dim %d caused division by zero \n", sw, dim);
        }
    }
    if (dim == 2) {
        double det = mat[0] * mat[3] - mat[1] * mat[2];
        if (fabs(det) > 0.000001) {
            double invDet = 1 / det;
            double a = mat[0];
            mat[0] = mat[3] * invDet;
            mat[1] = -mat[1] * invDet;
            mat[2] = -mat[2] * invDet;
            mat[3] = a * invDet;
        } else {
            outInvArr[1]++;
//            printf("switch %d dim %d caused division by zero \n", sw, dim);
        }
    }

}

void identityMatrix(int n, double *mat) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                mat[i * n + j] = 1;
            } else
                mat[i * n + j] = 0;
        }
    }
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

void smw(int numOfNodes, int numOfSwitches, int gon, double *selectedBase, double *curTs, double *switchNodeMat,
         double *AbaseInv, double *outInvArr) {

    double *u, *v, *vAbaseInv, *vAbaseInvU, *D, *DinvVAbaseInv, *AbaseInvU, *AbaseInvUDinvVAbaseInv;
    int *affectedNodes, nodeCount;

    outInvArr[0]=0;
    outInvArr[1]=0;

    //iterate over switches at time step ts
    for (int i = 0; i < numOfSwitches; ++i) {
        // Take the effects of i-th switch into account
        if (curTs[i] != selectedBase[i]) { //if i-th switch status changed compare to base switch pattern
            affectedNodes = mxCalloc(numOfNodes + 1, sizeof(*affectedNodes));
            //iterate over nodes in switchNodeMat to see which nodes are affected by switch i
            for (int j = 0; j < numOfNodes; ++j) {
                //if j-th node is affected by i-th switch
                if (switchNodeMat[i * numOfNodes + j] == 1) {
                    affectedNodes[j] = j;
                    affectedNodes[numOfNodes] = affectedNodes[numOfNodes] + 1;
                } else
                    affectedNodes[j] = -1;
            }
            u = mxCalloc(numOfNodes * affectedNodes[numOfNodes], sizeof(*u));
            v = mxCalloc(numOfNodes * affectedNodes[numOfNodes], sizeof(*v));

            if (affectedNodes[numOfNodes] == 1) {
                //iterate over nodes in switchNodeMat
                for (int m = 0; m < numOfNodes; ++m) {
                    if (affectedNodes[m] != -1) {
                        u[m] = (curTs[i] - selectedBase[i]) * gon;
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
                            u[n * 2] = (curTs[i] - selectedBase[i]) * gon;
                            u[n * 2 + 1] = (curTs[i] - selectedBase[i]) * (-gon);
                            v[n] = 1;
                            nodeCount++;
                        } else {
                            u[n * 2] = (curTs[i] - selectedBase[i]) * (-gon);
                            u[n * 2 + 1] = (curTs[i] - selectedBase[i]) * gon;
                            v[n + numOfNodes] = 1;
                        }
                    }
                }
            }
            //vAbaseInv <-- v * AbaseInv
            vAbaseInv = mxCalloc(affectedNodes[numOfNodes] * numOfNodes, sizeof(*vAbaseInv));
            mul_matrix_matrix(affectedNodes[numOfNodes], numOfNodes, numOfNodes, numOfNodes, v, AbaseInv, vAbaseInv);

            // vAbaseInvU <-- vAbaseInv * u
            vAbaseInvU = mxCalloc(affectedNodes[numOfNodes] * affectedNodes[numOfNodes], sizeof(*vAbaseInvU));
            mul_matrix_matrix(affectedNodes[numOfNodes], numOfNodes, numOfNodes, affectedNodes[numOfNodes], vAbaseInv,
                              u, vAbaseInvU);

            //D <-- I
            D = mxCalloc(affectedNodes[numOfNodes] * affectedNodes[numOfNodes], sizeof(*D));
            identityMatrix(affectedNodes[numOfNodes], D);

            //D <-- I + vAbaseInvU
            add_matrix(affectedNodes[numOfNodes], affectedNodes[numOfNodes], affectedNodes[numOfNodes],
                       affectedNodes[numOfNodes], D, vAbaseInvU);

            //D <-- Dinv
            invertMatrix(D, affectedNodes[numOfNodes], i, outInvArr);

            //DinvVAbaseInv <-- Dinv * vAbaseInv
            DinvVAbaseInv = mxCalloc(affectedNodes[numOfNodes] * numOfNodes, sizeof(*DinvVAbaseInv));
            mul_matrix_matrix(affectedNodes[numOfNodes], affectedNodes[numOfNodes], affectedNodes[numOfNodes],
                              numOfNodes, D, vAbaseInv, DinvVAbaseInv);

            //AbaseInvU <-- AbaseInv * u
            AbaseInvU = mxCalloc(numOfNodes * affectedNodes[numOfNodes], sizeof(*AbaseInvU));
            mul_matrix_matrix(numOfNodes, numOfNodes, numOfNodes, affectedNodes[numOfNodes], AbaseInv, u, AbaseInvU);

            //AbaseInvUDinvVAbaseInv <-- AbaseInvU * DinvVAbaseInv
            AbaseInvUDinvVAbaseInv = mxCalloc(numOfNodes * numOfNodes, sizeof(*AbaseInvUDinvVAbaseInv));
            mul_matrix_matrix(numOfNodes, affectedNodes[numOfNodes], affectedNodes[numOfNodes], numOfNodes, AbaseInvU,
                              DinvVAbaseInv, AbaseInvUDinvVAbaseInv);

            //AbaseInv <-- AbaseInv - AbaseInv * u * Dinv * v * AbaseInv === AbaseInv - AbaseInvUDinvVAbaseInv
            sub_matrix(numOfNodes, numOfNodes, numOfNodes, numOfNodes, AbaseInv, AbaseInvUDinvVAbaseInv);

        }
    }
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    int numOfNodes;
    int numOfSwitches;
    int gon;
    double *selectedBase;       /* of size (numOfSwitches + 1)*/
    double *curTT;           /* of size (numOfSwitches)*/
    double *switchNodeMat; /* of size (numOfSwitches*numOfNodes)*/
    double *outMatrix;      /* of size (numOfNode*numOfNodes)*/
    double *outInvArr;      /* of size (numOfNode*numOfNodes)*/

    /* check for proper number of arguments */
    if (nrhs != 7) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "7 inputs required.");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "2 output required.");
    }

    /* get the value of the scalar input */
    numOfNodes = (int) mxGetScalar(prhs[0]);
    numOfSwitches = (int) mxGetScalar(prhs[1]);
    gon = (int) mxGetScalar(prhs[2]);

    /* create a pointer to the real data in the input matrix  */
    selectedBase = mxGetPr(prhs[3]);
    curTT = mxGetPr(prhs[4]);
    switchNodeMat = mxGetPr(prhs[5]);

    /* create the output matrix */
    plhs[0] = mxDuplicateArray(prhs[6]); //prhs[6] == AbaseInv
    plhs[1] = mxCreateDoubleMatrix(1, 2, mxREAL);
    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);
    outInvArr = mxGetPr(plhs[1]);

    /* call the computational routine */
    smw(numOfNodes, numOfSwitches, gon, selectedBase, curTT, switchNodeMat, outMatrix, outInvArr);
}
