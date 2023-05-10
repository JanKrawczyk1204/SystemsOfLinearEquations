from create_matrix import create_A, create_b
from solution import jacoby, gauss_seidler, solve_LU
import matplotlib.pyplot as plt


def zadanieB():
    size = 993
    A = create_A(size, 12, -1, -1)
    b = create_b(size)
    jacoby_result, jacoby_iteration, jacoby_time = jacoby(A, b)
    gauss_result, gauss_iteration, gauss_time = gauss_seidler(A, b)
    print("Zadania A i B")
    print("Jacoby iteration = " + str(jacoby_iteration) + ", time = " + str(jacoby_time))
    print("Gauss-Seidler itteration = " + str(gauss_iteration) + ", time = " + str(gauss_time))

def zadanieC():
    size = 993
    A = create_A(size, 3, -1, -1)
    b = create_b(size)
    jacoby_result, jacoby_iteration, jacoby_time = jacoby(A, b)
    gauss_result, gauss_iteration, gauss_time = gauss_seidler(A, b)
    print("Zadania C")
    print("Jacoby iteration = " + str(jacoby_iteration) + ", time = " + str(jacoby_time))
    print("Gauss-Seidler iteration = " + str(gauss_iteration) + ", time = " + str(gauss_time))

def zadanieD():
    size = 993
    A = create_A(size, -3, -1, -1)
    b = create_b(size)
    solution, residuum, time = solve_LU(A, b)
    print("Zadanie D")
    print(residuum)
    print(time)

def zadanieE():
    size = [100, 200, 500, 1000, 1500, 2000]
    jacoby_time = []
    gauss_time = []
    lu_time = []
    for i in size:
        print(i)
        A = create_A(i, 12, -1, -1)
        b = create_b(i)
        _, _, jacoby_t = jacoby(A, b)
        _, _, gauss_t = gauss_seidler(A, b)
        _, _, lu_t = solve_LU(A, b)
        jacoby_time.append(jacoby_t)
        gauss_time.append(gauss_t)
        lu_time.append(lu_t)
    plt.plot(size, jacoby_time, label='Jacoby')
    plt.plot(size, gauss_time, label = 'Gauss-Seidler')
    plt.plot(size, lu_time, label='faktoryzacja LU')
    plt.xlabel("Rozmiar macierzy")
    plt.ylabel("Czas [s]")
    plt.title("Czas trwania algorytmów")
    plt.legend()
    plt.show()

if __name__ == '__main__':
    zadanieE()