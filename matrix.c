#include <stdio.h>
#include <malloc.h>

float** createArray(int m, int n)
{
    float* values = calloc(m*n, sizeof(float));
    float** rows = malloc(n*sizeof(float*));
    for (int i=0; i<n; ++i)
    {
        rows[i] = values + i*m;
    }
    return rows;
}

/*   float **mat = (float **) malloc(sizeof(float *) * rows);
   mat[i] = (float *) malloc(sizeof(float) * t);

   for (int i = 0; i < rows; i++) {
       for (int j = 0; j < cols; j++) {
           mat[i][j] = val;
       }
   }
   return mat;
}*/

void showMat(float **mat, int row,int col) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            printf("%d ",*(*(mat+i)+j));
        }
        printf("\n");
    }
}

void main() {
    float** arr = createArray(2,2);
    arr[0][0] = 12;
    arr[0][1] = 70;
    arr[1][0] = 53;
    arr[1][1] = 27;
    showMat(arr, 2,2);
}