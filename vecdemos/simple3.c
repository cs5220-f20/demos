#include <stdalign.h>

#define SIZE 65536

struct data {
	alignas(32) double vec[SIZE];
};

void simple_loop(struct data* restrict a, struct data* restrict b)
{
    for (int i = 0; i < SIZE; ++i)
        a->vec[i] += b->vec[i];
}
