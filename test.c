#include <omp.h>
#include <stdio.h>
#include <math.h>

void myFunc();

int main() {

    double time = omp_get_wtime();
# pragma omp parallel
    {
#pragma omp for collapse(2)
        for (int i = 0; i < 10; i++) {
//            printf("Thread %d\n", omp_get_thread_num());
            for (int j = 0; j < 3; j++) {
                sqrt(4.0);
                printf("(i,j) = (%d,%d) Thread %d\n", i, j, omp_get_thread_num());
            }
        }
    }
    printf("time: %lf", omp_get_wtime() - time);
    return 0;
}

void myFunc(int i) {
    for (int j = 0; j < 3; j++) {
        sqrt(4.0);
        printf("(i,j) = (%d,%d) Thread %d\n", i, j, omp_get_thread_num());
    }
}