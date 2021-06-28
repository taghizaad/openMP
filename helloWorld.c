#include <stdio.h>
#include <omp.h> //defines any types that OpenMP uses and any function prototype

int main() {
    /**
     * it says give me bunch of threads and since I didn't tell it,
     * it gives me the default number of threads, and I don't know and
     * I don't care the system chooses what that default number is
     */
    #pragma omp parallel
    {
        /** I am using runtime library routine to set ID
         * that's gonna give me an identifier for each thread
         * and it's gonna range from 0 up to (number of threads-1)
         * so it is unique identifier for each thread or it's
         * thread ID
         */
        int ID = omp_get_thread_num(); // I am using runtime library routine to set ID
        printf("hello(%d) ", ID);
        printf("world(%d)\n ", ID);
    }
}
