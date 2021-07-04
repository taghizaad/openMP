#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>


void init_matrix(float *matrix, int rows, int cols) {
    int i,j;
    srand(time(0));
    for(i=0; i < rows; i++) {
        for(j=0; j < cols; j++) {
            matrix[i*cols +j] = rand();
        }
    }
}

void show_matrix(float *matrix, int row, int col){
    for(int i=0; i < row; i++) {
        for(int j=0; j < col; j++) {
            printf("%f ",matrix[i+j]);
        }
        printf("\n");
    }
}

void main() {

    float *matrix;
    int rows = 3, cols = 2;
    matrix = malloc(sizeof(float)*rows*cols);
    if(matrix != NULL) {
        init_matrix(matrix, rows, cols);
        show_matrix(matrix,rows, cols);
    }




}