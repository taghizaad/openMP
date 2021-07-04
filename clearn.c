#include <stdio.h>
#include <malloc.h>

void multiple(int a, int i, int i1);

int main(int argc, char *argv[])
{
    int a=5,b=6,*resultado;
    resultado=(int *) malloc(sizeof(int));
    multiple(a,b,&resultado);
}

int *set_up(int number_students) {
    int *vector = (int *) calloc(number_students, sizeof(int));
    return vector;
}