/*
 * what I am using is an algorithm strategy called SPMD (single program, multiple data)
 * this is very commonly used in parallel computing.
 * the idea is I have a single program and I run multiple copies of that program
 * so each thread has its own copy
 * It's gonna use the ID of the thread and knowledge of how many threads there are
 * to adjust and change what it does
 * */

#include <stdio.h>
#include <omp.h>

static long num_steps = 100000000;
double step;

#define NUM_THREADS 2

void main() {
    int i, nthreads;
    double pi, sum[NUM_THREADS];
    step = 1.0 / (double) num_steps;
    /**
     * this runtime library routine sets the number
     * of threads. it requests a number of threads
     */
    omp_set_num_threads(NUM_THREADS);
    double startTime = omp_get_wtime();
    #pragma omp parallel
    {
        /*
         * outside the parallel region, I've lost the variables
         * inside the thread stack (following variables). at the end of the
         * parallel region, the threads go away except the master thread.
         * That means that if I sum up the value of summation inside the
         * parallel region how does it get out for me to do something with it
         * when it's done. The sum variable somehow has to be shared so that it
         * can be visible after the parallel regions done. (think about race
         * condition).
         * I used a very common trick: I promoted the scalar sum to an array
         * */
        int i, id, nthrds;
        double x;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        /*
         * I had one thread if id=0, so the master thread, saved a copy of
         * the number of threads and notice in the code nthreads was declared
         * outside the parallel region so that it would be visible after the
         * parallel region was done.
         * when you enter parallel region in OpenMp you request a number of
         * threads but an environment can choose to give you fewer threads.
         * So I actually have to check how many threads I got. That's why I
         * picked one thread, doesn't matter which one I picked, to say hey
         * will you copy the number of threads as you get we have to work with
         * inside the parallel region? will you copy that into shared variable?
         * because that way I'll know how many there were when we get out
         * */
        if (id == 0) nthreads = nthrds;
        /* I split up loop iterations between threads
         * I did the easy cheat way. It is called a cyclic distribution
         * of the loop iterations. It sometimes is called Round-robin
         * distribution
         * */
        for (i = id, sum[id] = 0.0; i < num_steps; i = i + nthrds) {
            x = (i + .5) * step;
            sum[id] += 4.0 / (1.0 + x * x);
        }
    }
    for (i = 0, pi = 0.0; i < nthreads; i++) {
        pi += sum[i] * step;
    }
    double endTime = omp_get_wtime();
    printf("nthread: %d\n", nthreads);
    printf("pi value is: %f\n", pi);
    printf("time consumed: %f", endTime - startTime);
}



