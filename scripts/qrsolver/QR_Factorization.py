##  @package QR_Factorization.py
#   Python program to quickly check out and understand
#   QR Factorization method for solving matrix equations
#
#   The matrices are printed at each step of factorization
#   to understand the transformation of the matrices

import numpy as np

A = np.matrix(
    [
        [1, 1, 0, 0, 0],
        [1, 2, 2, 0, 0],
        [0, 2, 3, 3, 0],
        [0, 0, 3, 4, 4],
        [0, 0, 0, 4, 5]
    ],
    dtype=float
)

print("A:\n", A)

R = np.matrix(A)
Q = np.matrix(np.identity(5))

for i in range(4):

    P = np.identity(5)
    h = np.sqrt(R[i,i]**2 + R[i+1,i]**2)
    c = R[i,i] / h
    s = R[i+1,i] / h
    P[i,i] = P[i+1,i+1] = c
    P[i+1,i] = -s
    P[i,i+1] = s

    Q = P*Q
    R = P*R

    print("Q:\n", Q)
    print("R:\n", R)

print(Q.getI()*R)


