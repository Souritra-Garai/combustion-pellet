import numpy as np

A = np.matrix(
    [
        [1, 2, 0, 0, 0],
        [3, 4, 5, 0, 0],
        [0, 6, 7, 8, 0],
        [0, 0, 9,10,11],
        [0, 0, 0,12,13]
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


