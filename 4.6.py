import numpy as np
import scipy.linalg

# циклическая трехдиагональная матрица
n = 10
h = 1/n
a = np.zeros((n,n))
for i in range(n):
    a[i,i] = 2 + h**2
    if i == 0:
        a[i,-1] = -1
        a[i,1] = -1
    elif i == n-1:
        a[i,0] = -1
        a[i,-2] = -1
    else:
        a[i,i-1] = -1
        a[i,i+1] = -1

# правая сторона
f = (1 + 4/(h**2)*np.sin(np.pi*h)**2)*np.sin(2*np.pi*np.arange(n)*h)

# Разложение LU
l, u = scipy.linalg.lu(a, permute_l=True)

# решение системы
y = scipy.linalg.solve(l, f, lower=True)
x = scipy.linalg.solve(u, y)

# точное решение
y_exact = np.sin(2*np.pi*np.arange(n)*h)

print("Solution:\n", x)
print("Error:\n", np.linalg.norm(x-y_exact, np.inf))