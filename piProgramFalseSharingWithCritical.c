//
// Created by h.taghizad on 6/27/2021.
//
#include <stdio.h>
#include <omp.h>

static long num_steps = 100000000;
double step;

#define NUM_THREADS 2

void main() {
    int nthreads;
    double pi = 0.0;
    step = 1.0 / (double) num_steps;
    omp_set_num_threads(NUM_THREADS);
    double startTime = omp_get_wtime();
#pragma omp parallel
    {
        int i, id, nthrds;
        /*
         * Now each thread has it's own copy of sum and because it's on stack
         * and it's not on the heap so the chances having false sharing are zilch
         * */
        double x, sum;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        if (id == 0) nthreads = nthrds;
        for (i = id, sum = 0.0; i < num_steps; i = i + nthreads) {
            x = (i + .5) * step;
            sum += 4.0 / (1.0 + x * x);
        }
        /*
         * I got to save sum into global pi. That's a potential for race condition
         * because if I have all the threads with undisciplined access like summing
         * pi willy-nilly, I could be halfway through updating one sum when I start
         * updating on the other one. It's messy so I have to have mutual exclusion
         * of that sum in the pi
         * */
        #pragma omp critical
        pi += sum * step;
    }
    double endTime = omp_get_wtime();
    printf("nthread: %d\n", nthreads);
    printf("pi value is: %f\n", pi);
    printf("time consumed: %f", endTime - startTime);
}
