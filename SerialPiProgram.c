/*
 * in this program we are going to calculate the integral
 * of the following function between 0 to 1
 *  f(x) = 4/(1 + x*x)
 * with analytical solution we know this integral is equal pi
 * */

#include <stdio.h>
#include <omp.h>

static long num_steps = 100000000;
double step;

void main() {
    int i;
    double pi, x, sum = 0.0;
    step = 1.0 / (double) num_steps;
    double startTime = omp_get_wtime(); //time in seconds since a fixed point in the past

    for (i = 0; i < num_steps; i++) {
        x = (i + .5) * step;
        sum += 4.0 / (1.0 + x * x);
    }
    pi = step * sum;
    double endTime = omp_get_wtime();
    printf("pi value is: %f\n", pi);
    printf("time consumed: %f", endTime - startTime);
}



