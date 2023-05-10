import math
import sys
import time

def residuum(a, b, x):
    n = len(a)
    res = []
    for i in range(n):
        row = 0
        for j in range(n):
            row += a[i][j] * x[j]
        res.append(row - b[i])
    norm = 0
    for i in range(n):
        if norm + res[i] ** 2 < sys.maxsize:
            norm += res[i] ** 2
    norm = math.sqrt(norm)
    return norm



def jacoby(a, b):
    n = len(b)
    start_time = time.time()
    itteration = 0
    x = [0 for i in range(n)]
    x_prev = [0 for i in range(n)]
    eps = 1e-9
    a = [[float(x) for x in row] for row in a]
    while True:
        itteration += 1
        for i in range(n):
            x_prev[i] = x[i]
            sigma = 0
            for j in range(n):
                if j != i:
                    sigma += a[i][j] * x_prev[j]
            x[i] = (b[i] - sigma) / a[i][i]
        if residuum(a, b, x) < eps:
            total_time = time.time() - start_time
            break
    return x, itteration, total_time, residuum(a, b, x)


def gauss_seidler(a, b):
    n = len(b)
    start_time = time.time()
    itteration = 0
    x = [0 for i in range(n)]
    x_prev = [0 for i in range(n)]
    eps = 1e-9
    while True:
        itteration += 1
        for i in range(n):
            x_prev[i] = x[i]
            sigma1 = sum(a[i][j] * x[j] for j in range(i))
            sigma2 = sum(a[i][j] * x_prev[j] for j in range(i+1, n))
            x[i] = (b[i] - sigma1 - sigma2) / a[i][i]
        if residuum(a, b, x) < eps:
            total_time = time.time() - start_time
            break
    return x, itteration, total_time, residuum(a, b, x)

def LU_factorization(a):
    n = len(a)
    L = [[0.0] * n for i in range(n)]
    U = [[0.0] * n for i in range(n)]
    for i in range(n):
        L[i][i] = 1.0
    for k in range(n):
        U[k][k] = a[k][k]
        for j in range(k+1, n):
            L[j][k] = a[j][k] / U[k][k]
            U[k][j] = a[k][j]
        for i in range(k+1, n):
            for j in range(k+1, n):
                a[i][j] = a[i][j] - L[i][k] * U[k][j]
    return L, U

def solve_LU(a, b):
    n = len(a)
    start_time = time.time()
    L, U = LU_factorization(a)
    y = [0.0] * n
    x = [0.0] * n
    # solve Ly = b
    for i in range(n):
        y[i] = b[i]
        for j in range(i):
            y[i] -= L[i][j] * y[j]
    # solve Ux = y
    for i in range(n-1, -1, -1):
        x[i] = y[i]
        for j in range(i+1, n):
            x[i] -= U[i][j] * x[j]
        x[i] /= U[i][i]
    total_time = time.time() - start_time
    return x, residuum(a, b, x), total_time
