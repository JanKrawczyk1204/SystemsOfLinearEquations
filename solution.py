import math

import numpy
def matrix_multiply(matrix1, matrix2):
    rows1 = len(matrix1)
    cols1 = len(matrix1[0])
    rows2 = len(matrix2)
    cols2 = len(matrix2[0])
    result = [[0 for j in range(cols2)] for i in range(rows1)]
    for i in range(rows1):
        for j in range(cols2):
            for k in range(cols1):
                result[i][j] += matrix1[i][k] * matrix2[k][j]
    return result


def subtract_matrix_vector(matrix, vector):
    n = len(matrix)
    result = []
    for i in range(n):
        row = []
        for j in range(n):
            row.append(matrix[i][j] - vector[j])
        result.append(row)
    return result

def add_matrix_vector(matrix, vector):
    n = len(matrix)
    result = []
    for i in range(n):
        row = []
        for j in range(n):
            row.append(matrix[i][j] + vector[j])
        result.append(row)
    return result


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
        norm += res[i] ** 2
    norm = math.sqrt(norm)
    return norm



def jacoby(a, b):
    n = len(b)
    x = [0 for i in range(n)]
    x_prev = [0 for i in range(n)]
    eps = 1e-6
    a = [[float(x) for x in row] for row in a]
    while True:
        for i in range(n):
            x_prev[i] = x[i]
            sigma = 0
            for j in range(n):
                if j != i:
                    sigma += a[i][j] * x_prev[j]
            x[i] = (b[i] - sigma) / a[i][i]
        if residuum(a, b, x) < eps:
            break
    return x


def gauss_seidler(a, b):
    n = len(b)
    x = [0 for i in range(n)]
    x_prev = [0 for i in range(n)]
    eps = 1e-6
    while True:
        for i in range(n):
            x_prev[i] = x[i]
            sigma1 = sum(a[i][j] * x[j] for j in range(i))
            sigma2 = sum(a[i][j] * x_prev[j] for j in range(i+1, n))
            x[i] = (b[i] - sigma1 - sigma2) / a[i][i]
        if residuum(a, b, x) < eps:
            break
    return x

