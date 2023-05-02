from create_matrix import create_A, crate_b,print_array

if __name__ == '__main__':
    size = 10
    A = create_A(size, 12, -1, -1)
    b = crate_b(size)
    print_array(b)
