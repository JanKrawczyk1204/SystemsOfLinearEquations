from create_matrix import create_A, create_b
from solution import jacoby, gauss_seidler, solve_LU
import matplotlib.pyplot as plt


def zadanieB():
    size = 993
    A = create_A(size, 12, -1, -1)
    b = create_b(size)
    _, jacoby_iteration, jacoby_time, _ = jacoby(A, b)
    _, gauss_iteration, gauss_time, _ = gauss_seidler(A, b)
    print("Zadania A i B")
    print("Jacoby iteration = " + str(jacoby_iteration) + ", time = " + str(jacoby_time))
    print("Gauss-Seidler itteration = " + str(gauss_iteration) + ", time = " + str(gauss_time))

def zadanieC():
    size = 993
    A = create_A(size, 3, -1, -1)
    b = create_b(size)
    _, jacoby_iteration, jacoby_time, _ = jacoby(A, b)
    _, gauss_iteration, gauss_time, _ = gauss_seidler(A, b)
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
    size = [10, 20, 50, 100, 200, 500, 600, 700, 800]
    jacoby_time = []
    gauss_time = []
    lu_time = []
    jacoby_residuum = []
    gauss_residuum = []
    lu_residuum = []
    for i in size:
        print(i)
        A = create_A(i, 12, -1, -1)
        b = create_b(i)
        _, _, jacoby_t, jacoby_r = jacoby(A, b)
        _, _, gauss_t, gauss_r = gauss_seidler(A, b)
        _, lu_r, lu_t = solve_LU(A, b)
        jacoby_time.append(jacoby_t)
        gauss_time.append(gauss_t)
        lu_time.append(lu_t)
        jacoby_residuum.append(jacoby_r)
        gauss_residuum.append(gauss_r)
        lu_residuum.append(lu_r)
    plt.plot(size, jacoby_time, label='Jacoby')
    plt.plot(size, gauss_time, label = 'Gauss-Seidler')
    plt.plot(size, lu_time, label='faktoryzacja LU')
    plt.xlabel("Rozmiar macierzy")
    plt.ylabel("Czas [s]")
    plt.title("Czas trwania algorytmów")
    plt.legend()
    plt.show()
    plt.plot(size, jacoby_residuum, label='Jacoby')
    plt.plot(size, gauss_residuum, label='Gauss-Seidler')
    plt.plot(size, lu_residuum, label='faktoryzacja LU')
    plt.xlabel("Rozmiar macierzy")
    plt.ylabel("Czas [s]")
    plt.title("Wartości normy residuum")
    plt.legend()
    plt.show()