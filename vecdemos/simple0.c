#define SIZE 65536

void simple_loop(double* a, double* b)
{
    for (int i = 0; i < SIZE; ++i)
        a[i] += b[i];
}
