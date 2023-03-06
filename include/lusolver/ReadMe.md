# Tridiagonal Matrix Equation Solver - LU Method

This library can be used to solve the Matrix equation 
```
A.x = b
```
using Thomas algorithm / LU method.
'`A`' must be a tridiagonal matrix, i.e.
```
A[i][j] = 0

	for 
		j < i - 1
		j > i + 1
```
where, `i` and `j` are row and column indices respectively. Both indices start from `0` and go up to `n-1` for an `n x n` matrix.

## Tridiagonal Matrix

The class `TridiagonalMatrix` implements a space efficient method to store the three diagonal entries in a tridiagonal matrix. It manipulates the row and column indices to store the values in a `3n-2` flat array.
```
// Row # i
// Column # j
A[i][j] = array[ 3 * i + (j - i) ]
=> array[] = {A[0][0], A[0][1], A[1][0], A[1][1], A[1][2], A[2][1], A[2][2], A[2][3], ... }
```
It is assumed indices like `A[1][3]`, `A[2][0]` etc., i.e., where either `j < i-1` or `j > i+1`, are not accessed. Accessing such indices can result in either segmentation fault or corruption of data.

## LU Solver

The class `LUSolver` lets the user enter one row equation at a time for tridiagonal matrix equations. Then it can solve the matrix equation at one go using LU decomposition and forward substitution.

For certain applications, the calculations before setting one row equation takes considerable time. Hence, this class provides function to parallely compute the necessary values and enter one row at time.

In the LU decomposition method, the matrix `A` is decomposed into the multiplication of a lower and an upper diagonal matrices, i.e.,
```
A = L.U
```
such that
```
L[i][j] = 0

	for
		j < i-1
		j > i

L[i][j] = 1

	for j = i
```
```
U[i][j] = 0

	for
		j < i
		j > i+1
```

This gives the following constraints
```
A[i][i-1] = L[i][i-1] * U[i-1][i-1]
A[i][i]   = L[i][i-1] * U[i-1][i] + U[i][i]
A[i][i+1] = U[i][i+1]
```

Using Thomas' Algorithm
```
U[0][0] = A[0][0]
U[0][1] = A[0][1]

for iteration # i from 1 to n-1

	// U[i-1][i-1] and U[i-1][i] were already set in previous iteration

	L[i][i-1] = A[i][i-1] / U[i-1][i-1]

	U[i][i]   = A[i][i] - L[i][i-1] * U[i-1][i]
	U[i][i+1] = A[i][i+1]
```

