#include "matrix.h"
#include "mex.h"
#include <math.h>

void sub_matrix(int row1, int col1, int row2, int col2, double *mat1, double *mat2, double *result);

void invertMatrix(double *mat, int dim, int sw, double *outInvArr);

void identityMatrix(int n, double *mat);

void add_matrix(int row1, int col1, int row2, int col2, double *mat1, double *mat2);

void mul_matrix_matrix(int row1, int col1, int row2, int col2, double *mat1, double *mat2, double *result);

void show_matrix(double *matrix, size_t row, size_t col);

void
smwWithWholdeD(int numOfNodes, int numOfSwitches, int gon, double *selectedBase, double *curTs, double *switchNodeMat,
               double *AbaseInv, double *outMatrix);

void vecToMat(double *vec, int row, int col, double **matrix);

void matToVec(double **matrix, int row, int col, double *vec);

double **gauss_jordan_matrix_inverse(double **A, size_t size);

double **gauss_jordan_matrix_inverse(double **A, size_t size) {
    double **I, temp;
    int i, j, k;

    I = (double **) mxCalloc(
            size, sizeof(double *));            //memory allocation for indentity matrix I(matsize X matsize)
    for (i = 0; i < size; i++)
        I[i] = (double *) mxCalloc(size, sizeof(double));
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


void vecToMat(double *vec, int row, int col, double **matrix) {
    for (int i = 0; i < row; i++)
        matrix[i] = (double *) mxCalloc(col, sizeof(double));
    for (size_t i = 0; i < row; i++)
        for (size_t j = 0; j < col; j++)
            matrix[i][j] = vec[i * col + j];
}

void matToVec(double **matrix, int row, int col, double *vec) {
    for (size_t i = 0; i < row; i++)
        for (size_t j = 0; j < col; j++)
            vec[i * col + j] = matrix[i][j];
}

void show_matrix(double *matrix, size_t row, size_t col) {
    printf("sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss\n");
    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < col; j++)
            printf("%15.15e ", matrix[i * col + j]);
        printf("\n");
    }
    printf("eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee\n");
}

void sub_matrix(int row1, int col1, int row2, int col2, double *mat1, double *mat2, double *result) {
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col1; j++)
            result[i * col1 + j] = mat1[i * col1 + j] - mat2[i * col1 + j];
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

void smwWithWholdeD(int numOfNodes, int numOfSwitches, int gon, double *baseArr, double *curTs, double *switchNodeMat,
                    double *AbaseInv, double *outMatrix) {

    double *uu, *vv, *uv, *vAbaseInv, *vAbaseInvU, *D, *DinvVAbaseInv, *AbaseInvU, *AbaseInvUDinvVAbaseInv;
    int *affectedNodes, nodeCount;

    uu = (double *) mxCalloc(numOfNodes * numOfNodes, sizeof(double));
    vv = (double *) mxCalloc(numOfNodes * numOfNodes, sizeof(double));

    int *allAffectedNodes = (int *) mxCalloc(numOfNodes, sizeof(int)); //all affected nodes at curTS

    //iterate over switches at time step ts
    for (int i = 0; i < numOfSwitches; ++i) {

        // Take the effects of i-th switch into account
        if (curTs[i] != baseArr[i]) { //if i-th switch status changed compare to base switch pattern
            affectedNodes = (int *) mxCalloc(2, sizeof(int));
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
    double *u = (double *) mxCalloc(numOfNodes * numOfAllAffectedNodes, sizeof(double));
    double *v = (double *) mxCalloc(numOfNodes * numOfAllAffectedNodes, sizeof(double));

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
    uv = (double *) mxCalloc(numOfNodes * numOfNodes, sizeof(double));
    mul_matrix_matrix(numOfNodes, numOfAllAffectedNodes, numOfAllAffectedNodes, numOfNodes, u, v, uv);
//    show_matrix(uv, numOfNodes, numOfNodes);

    //vAbaseInv <-- v * AbaseInv
    vAbaseInv = (double *) mxCalloc(numOfAllAffectedNodes * numOfNodes, sizeof(double));
    mul_matrix_matrix(numOfAllAffectedNodes, numOfNodes, numOfNodes, numOfNodes, v, AbaseInv, vAbaseInv);

    // vAbaseInvU <-- vAbaseInv * u
    vAbaseInvU = (double *) mxCalloc(numOfAllAffectedNodes * numOfAllAffectedNodes, sizeof(double));
    mul_matrix_matrix(numOfAllAffectedNodes, numOfNodes, numOfNodes, numOfAllAffectedNodes, vAbaseInv,
                      u, vAbaseInvU);

    //D <-- I
    D = (double *) mxCalloc(numOfAllAffectedNodes * numOfAllAffectedNodes, sizeof(double));
    identityMatrix(numOfAllAffectedNodes, D);


    //D <-- I + v * AbaseInv * u ==== I + vAbaseInvU
    add_matrix(numOfAllAffectedNodes, numOfAllAffectedNodes, numOfAllAffectedNodes,
               numOfAllAffectedNodes, D, vAbaseInvU);


    double **Dformated = (double **) mxCalloc(numOfAllAffectedNodes, sizeof(double));
    vecToMat(D, numOfAllAffectedNodes, numOfAllAffectedNodes, Dformated);
    double **DinvFormated = gauss_jordan_matrix_inverse(Dformated, numOfAllAffectedNodes);
    double *Dinv = (double *) mxCalloc(numOfAllAffectedNodes * numOfAllAffectedNodes, sizeof(double));
    matToVec(DinvFormated, numOfAllAffectedNodes, numOfAllAffectedNodes, Dinv);


    //DinvVAbaseInv <-- Dinv * vAbaseInv
    DinvVAbaseInv = (double *) mxCalloc(numOfAllAffectedNodes * numOfNodes, sizeof(double));
    mul_matrix_matrix(numOfAllAffectedNodes, numOfAllAffectedNodes, numOfAllAffectedNodes,
                      numOfNodes, Dinv, vAbaseInv, DinvVAbaseInv);

    //AbaseInvU <-- AbaseInv * u
    AbaseInvU = (double *) mxCalloc(numOfNodes * numOfAllAffectedNodes, sizeof(double));
    mul_matrix_matrix(numOfNodes, numOfNodes, numOfNodes, numOfAllAffectedNodes, AbaseInv, u, AbaseInvU);

    //AbaseInvUDinvVAbaseInv <-- AbaseInvU * DinvVAbaseInv
    AbaseInvUDinvVAbaseInv = (double *) mxCalloc(numOfNodes * numOfNodes, sizeof(double));
    mul_matrix_matrix(numOfNodes, numOfAllAffectedNodes, numOfAllAffectedNodes, numOfNodes, AbaseInvU,
                      DinvVAbaseInv, AbaseInvUDinvVAbaseInv);

    //AbaseInv <-- AbaseInv - AbaseInv * u * Dinv * v * AbaseInv === AbaseInv - AbaseInvUDinvVAbaseInv
    sub_matrix(numOfNodes, numOfNodes, numOfNodes, numOfNodes, AbaseInv, AbaseInvUDinvVAbaseInv, outMatrix);

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
    double *AbaseInv;   /* of size (numOfNode*numOfNodes)*/
    double *outMatrix;


    /* check for proper number of arguments */
    if (nrhs != 7) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "7 inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "1 output required.");
    }

    /* get the value of the scalar input */
    numOfNodes = (int) mxGetScalar(prhs[0]);
    numOfSwitches = (int) mxGetScalar(prhs[1]);
    gon = (int) mxGetScalar(prhs[2]);

    /* create a pointer to the real data in the input matrix  */
    selectedBase = mxGetPr(prhs[3]);
    curTT = mxGetPr(prhs[4]);
    switchNodeMat = mxGetPr(prhs[5]);
    AbaseInv = mxGetPr(prhs[6]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1, (mwSize) numOfNodes * (mwSize) numOfNodes, mxREAL);
    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    smwWithWholdeD(numOfNodes, numOfSwitches, gon, selectedBase, curTT, switchNodeMat, AbaseInv, outMatrix);
    show_matrix(outMatrix, numOfNodes, numOfNodes);
}
