#define SIZE 65536

#define GCC_ALN(var, alignment) \
	__builtin_assume_aligned(var, alignment)

void simple_loop(double* restrict a, double* restrict b)
{
    a = GCC_ALN(a, 32);
    b = GCC_ALN(b, 32);
    for (int i = 0; i < SIZE; ++i)
        a[i] += b[i];
}
