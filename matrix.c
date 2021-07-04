#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>

/**
 *  this method shows a vector as a matrix by row*col
 * @param vector
 * @param row
 * @param col
 */
void show_vector(float *vector, int row, int col) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            printf("%f ", vector[i + j]);
        }
        printf("\n");
    }
}


/**
 *
 * @param rows
 * @param cols
 * @return  a vector of float initialized to the zero value
 */
float *init_vector(int rows, int cols) {
    float *vector = (float *) calloc(rows * cols, sizeof(float));
    return vector;
}

/**
 * main method
 */
void main() {

    float *vector;
    int rows = 2, cols = 10;
    vector = init_vector(rows, cols);
    if (vector != NULL) {
        show_vector(vector, rows, cols);
    }
}