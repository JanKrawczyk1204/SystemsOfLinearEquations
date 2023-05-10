import math


def create_A(size, a1, a2, a3):
    array = [[0 for j in range(size)] for i in range(size)]
    for i in range(size):
        array[i][i] = a1
        if i > 0:
            array[i][i - 1] = a2
        if i < size - 1:
            array[i][i + 1] = a2
        if i > 1:
            array[i][i - 2] = a3
        if i < size - 2:
            array[i][i + 2] = a3

    return array


def create_b(size):
    array = []
    for i in range(size):
        array.append(math.sin((i+1)*9))
    return array


def print_array(array):
    for row in array:
        print(row)
