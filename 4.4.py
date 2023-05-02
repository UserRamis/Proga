import numpy as np
from scipy.linalg import lu

n = 4
A = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        A[i, j] = min(i+1, j+1) / max(i+1, j+1)

# LU разложение с выбором главного элемента по строке
LU, piv = lu(A, permute_l=True)

# проверка, является ли A обратимым
if np.linalg.det(A) == 0:
    print("Матрица не является обратимой")
else:
    
    E = np.eye(n)
    X = np.linalg.solve(LU, np.dot(piv, E))
    print("A =\n", A)
    print("A^-1 =\n", X)
    print(" A^-1 тридиагональная?", np.all(X[1:, :-1] == 0))
