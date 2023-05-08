from create_matrix import create_A, crate_b
from solution import jacoby, gauss_seidel

if __name__ == '__main__':
    size = 10
    A = create_A(size, 12, -1, -1)
    b = crate_b(size)
    print(jacoby(A, b))
    print("----------------")
    print(gauss_seidel(A, b))
