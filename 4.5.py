import numpy as np

#Функция для вычисления нормы матрицы
def norm(A, p):
    if p == "1":
        return np.max(np.sum(np.abs(A), axis=0))
    elif p == "E":
        return np.linalg.norm(A, 2)
    elif p == "inf":
        return np.max(np.sum(np.abs(A), axis=1))
    else:
        return None

#Функция для LU-разложения матрицы
def LU_decomposition(A):
    piv = np.arange(A.shape[0])
    LU, piv = np.linalg.lu(A)
    L = np.tril(LU, -1) + np.eye(A.shape[0])
    U = np.triu(LU)
    return L, U, piv+1

#Функция для вычисления обратной матрицы
def inverse_matrix(A):
    return np.linalg.inv(A)

#Определим размер матрицы
n = 5
A = np.fromfunction(lambda i, j: np.minimum(i+1,j+1), (n,n))

#Вычислим число обусловленности матрицы с помощью норм 1, E и бесконечности
cond_1 = norm(A, "1") * norm(inverse_matrix(A), "1")
cond_E = norm(A, "E") * norm(inverse_matrix(A), "E")
cond_inf = norm(A, "inf") * norm(inverse_matrix(A), "inf")

print("Матрица A =\n", A)
print("Обратная матрица A =\n", inverse_matrix(A))
print("Число обусловленности матрицы (норма 1) =", cond_1)
print("Число обусловленности матрицы (норма E) =", cond_E)
print("Число обусловленности матрицы (норма бесконечность) =", cond_inf)