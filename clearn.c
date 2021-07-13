#include <omp.h>
#include <stdio.h>

int main() {
    omp_set_num_threads(4);
// Do this part in parallel
#pragma omp parallel
    {
        printf( "Hello, World!\n" );
    }
    return 0;
}