#define SIZE 65536

void simple_loop(double* restrict a, double* restrict b)
{
    for (int i = 0; i < SIZE; ++i)
        a[i] += b[i];
}
