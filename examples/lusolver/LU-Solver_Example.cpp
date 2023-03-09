#include "math/Data-Type.hpp"
#include "lusolver/LU-Solver.hpp"

#include <iostream>
#include <stdlib.h>
#include <time.h>

int main(int argc, char const *argv[])
{
	const unsigned int N = 10000;

	real_t x[N], b[N];

	TridiagonalMatrix A(N);

	srand(time(0));

	A.setElement(0, 0, rand());
	A.setElement(0, 1, rand());

	x[0] = rand();
	b[0] = 0;

	for (int i = 1; i < N-1; i++)
	{
		A.setElement(i, i-1, rand());
		A.setElement(i, i,   rand());
		A.setElement(i, i+1, rand());

		x[i] = rand();
		b[i] = 0;
	}

	A.setElement(N-1, N-2, rand());
	A.setElement(N-1, N-1, rand());

	x[N-1] = rand();
	b[N-1] = 0;

	A.multiply(x, b);

	LUSolver my_solver(N);

	my_solver.setEquationFirstRow(A.getElement(0, 1), A.getElement(0, 0), b[0]);

	for (int i = 1; i < N-1; i++)
	{
		my_solver.setEquationSerially(i, A.getElement(i, i+1), A.getElement(i, i), A.getElement(i, i-1), b[i]);
	}

	my_solver.setEquationLastRowSerially(A.getElement(N-1, N-1), A.getElement(N-1, N-2), b[N-1]);

	// my_solver.printMatrixEquation();

	real_t x_soln[N];

	my_solver.getSolutionSerially(x_soln);

	real_t MSE = 0;

	for (int i = 0; i < N; i++) MSE += (x[i] - x_soln[i]) * (x[i] - x_soln[i]);

	std::cout << "\nMSE : " << MSE << std::endl;

	return 0;
}