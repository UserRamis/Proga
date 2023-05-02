import numpy as np
import scipy.optimize

# Матрица 
A = np.array([[10, -1, 2],
              [-1, 11, -1],
              [2, -1, 10]])

# Вектор правой части
b = np.array([6, 25, -11])

# Начальное приближение
x0 = np.zeros_like(b)

# Функция для решения СЛАУ методом релаксации
def relax(x, A, b):
    return np.sum((np.dot(A, x) - b) ** 2)


cons = ({'type': 'eq', 'fun': lambda x: np.sum(x) - 1})

# Решение СЛАУ методом релаксации
res = scipy.optimize.minimize(relax, x0, args=(A, b), method='SLSQP', constraints=cons)

x = res.x

print(x)
