#include "lusolver/Tridiagonal_Matrix.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char const *argv[])
{
    const unsigned int N = 5;

    long double x[N], b[N];

	TridiagonalMatrix<long double> A(N);

    srand(time(0));

    A.setElement(0, 0, rand());
    A.setElement(0, 1, rand());

    x[0] = rand();
    b[0] = 0;

    for (int i = 1; i < N-1; i++)
    {
        A.setElement(i, i-1,    rand());
        A.setElement(i, i,      rand());
        A.setElement(i, i+1,    rand());

        x[i] = rand();
        b[i] = 0;
    }

    A.setElement(N-1, N-2, rand());
    A.setElement(N-1, N-1, rand());

    x[N-1] = rand();
    b[N-1] = 0;

    A.multiply(x, b);

	A.printMatrix();

    return 0;
}