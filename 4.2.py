import numpy as np
import scipy


n = 5
A = np.zeros((n,n))
for i in range(n):
    for j in range(n):
        A[i,j] = np.math.factorial(i+j) // (np.math.factorial(i)*np.math.factorial(j))


P, L, U = scipy.linalg.lu(A)


print("Матрица Паскаля A =\n", A)
print("L =\n", L)
print("U =\n", U)
print("P =\n", P)


det_A = np.linalg.det(A)
print("det(A) =", det_A)
for i in range(2, n+1):
    det_Ai = np.math.factorial(i-1)**2
    print(f"det(A_{i}) =", det_Ai)
