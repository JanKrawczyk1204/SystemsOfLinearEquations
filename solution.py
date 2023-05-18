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
        else:
            norm = sys.maxsize
    norm = math.sqrt(norm)
    return norm


def jacoby(a, b):
    n = len(b)
    start_time = time.time()
    residuum_table = []
    itteration = 0
    x = [0 for i in range(n)]
    eps = 1e-9
    a = [[float(x) for x in row] for row in a]
    while True:
        itteration += 1
        x_prev = x[:]
        for i in range(n):
            sigma = 0
            for j in range(n):
                if j != i:
                    sigma += a[i][j] * x_prev[j]
            x[i] = (b[i] - sigma) / a[i][i]
        current_residuum = residuum(a, b, x)
        residuum_table.append(current_residuum)
        if current_residuum < eps or itteration == 200:
            total_time = time.time() - start_time
            break
    return x, itteration, total_time, current_residuum, residuum_table


def gauss_seidler(a, b):
    n = len(b)
    residuum_table = []
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
        current_residuum = residuum(a, b, x)
        residuum_table.append(current_residuum)
        if current_residuum < eps or itteration == 200:
            total_time = time.time() - start_time
            break
    return x, itteration, total_time, current_residuum, residuum_table


def LU_factorization(A):
    N = len(A)
    L = [[0.0] * N for _ in range(N)]
    U = [[0.0] * N for _ in range(N)]

    for i in range(N):
        L[i][i] = 1.0

        for j in range(i + 1):
            s1 = 0.0
            for k in range(j):
                s1 += U[k][i] * L[j][k]
            U[j][i] = A[j][i] - s1

        for j in range(i, N):
            s2 = 0.0
            for k in range(i):
                s2 += U[k][i] * L[j][k]
            L[j][i] = (A[j][i] - s2) / U[i][i]

    return L, U


def solve_LU(A, b):
    start_time = time.time()
    L, U = LU_factorization(A)
    N = len(L)
    y = [0.0 for _ in range(N)]
    x = [0.0 for _ in range(N)]

    for i in range(N):
        sum_Ly = 0.0
        for j in range(i):
            sum_Ly += L[i][j] * y[j]
        y[i] = b[i] - sum_Ly

    for i in range(N - 1, -1, -1):
        sum_Ux = 0.0
        for j in range(i + 1, N):
            sum_Ux += U[i][j] * x[j]
        x[i] = (y[i] - sum_Ux) / U[i][i]

    return x, residuum(A, b, x), time.time() - start_time



