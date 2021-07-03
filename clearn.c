#include <stdio.h>
#include <malloc.h>

float **createArray(int m, int n, float value) {
    float *values = calloc(m * n, sizeof(float));
    float **rows = malloc(n * sizeof(float *));
    for (int i = 0; i < n; ++i) {
        rows[i] = values + i * m;
    }
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            rows[i][j] = value;
        }
    }
    return rows;
}

int **makeMatrix(int row, int col, int value) {
    int **matrix;
    matrix = (int **) calloc(col, sizeof(int *));
    for (int i = 0; i < col; i++) {
        matrix[i] = (int *) calloc(row, sizeof(int));
    }
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            matrix[i][j] = value;
        }
    }
    return matrix;
}

void showMatrix(int **matrix, int row, int col) {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            printf("%d ", (matrix)[i][j]);
        }
        printf("\n");
    }
}

void main() {
    int row = 3, col = 2, value = 7.2;
    float **matrix;
    matrix = makeMatrix(row, col, value);
//    matrix = createArray(row, col, value);
//    showMatrix(matrix, row , col);

}