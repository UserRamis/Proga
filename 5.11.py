import numpy as np
import time

# решение СЛАУ методом Якоби
A = np.array([[10, -1, 2], [-1, 11, -1], [2, -1, 10]])
b = np.array([6, 25, -11])
D = np.diag(A)
R = A - np.diagflat(D)
x_jacobi = np.zeros_like(b)
start_time = time.time()
for i in range(100):
    x_jacobi = (b - np.dot(R, x_jacobi)) / D# для каждой итерации выполняется матричное умножение и деление на диагональные элементы матрицы
end_time = time.time()

print("Решение методом Якоби:")
print(x_jacobi)
print("Время выполнения метода Якоби: {} секунд".format(end_time - start_time))

# решение СЛАУ методом Зейделя
A = np.array([[10, -1, 2], [-1, 11, -1], [2, -1, 10]])
b = np.array([6, 25, -11])
L = np.tril(A)
U = A - L
x_seidel = np.zeros_like(b)
start_time = time.time()
for i in range(100):
    x_seidel = np.dot(np.linalg.inv(L), b - np.dot(U, x_seidel))#numpy.linalg.inv
end_time = time.time()

print("Решение методом Зейделя:")
print(x_seidel)
print("Время выполнения метода Зейделя: {} секунд".format(end_time - start_time))
