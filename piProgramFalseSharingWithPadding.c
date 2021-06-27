/*
 * we fixed false sharing problem by padding the sum array
 * but it is not satisfactory because it required that we
 * knew the size of L1 cache line
 * If we move to another machine which has different cache
 * line size, that would be a mess!
 */
#include <stdio.h>
#include <omp.h>

static long num_steps = 100000000;
double step;
#define PAD 8   //assume 64 byte L1 cache line size
#define NUM_THREADS 2

void main() {
    int i, nthreads;
    double pi, sum[NUM_THREADS][PAD];
    step = 1.0 / (double) num_steps;
    omp_set_num_threads(NUM_THREADS);
    double startTime = omp_get_wtime();
#pragma omp parallel
    {
        int i, id, nthrds;
        double x;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        if (id == 0) nthreads = nthrds;
        for (i = id, sum[id][0] = 0.0; i < num_steps; i = i + nthreads) {
            x = (i + .5) * step;
            sum[id][0] += 4.0 / (1.0 + x * x);
        }
    }
    for (i = 0, pi = 0.0; i < nthreads; i++) {
        pi += sum[i][0] * step;
    }
    double endTime = omp_get_wtime();
    printf("nthread: %d\n", nthreads);
    printf("pi value is: %f\n", pi);
    printf("time consumed: %f", endTime - startTime);
}



