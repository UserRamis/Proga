import numpy as np

# Матрица коэффициентов
A = np.array([[1, 2, 3],
              [4, 5, 6],
              [7, 8, 9]])

# Вектор правой части
b = np.array([10, 11, 12])

# Решение СЛАУ методом наименьших квадратов
x, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)

print(x)
