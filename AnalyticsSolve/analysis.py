from sympy import I, shape
from sympy.matrices import Matrix, zeros
from sympy.printing import pprint

from variables import N


def load_fortran_matrix(fortran_matr_path: str, matrix_dim: int) -> Matrix:
    matrix = [[(0, 0) for _ in range(matrix_dim)] for _ in range(matrix_dim)]
    with open(fortran_matr_path, 'r') as log_file:
        lines = log_file.readlines()
        for line in lines:
            strip_line = line.replace(' ', '')
            if not strip_line.startswith('matr') or len(strip_line) == 0:
                continue

            index_part, value_part = strip_line.split('=')

            i, j = index_part.replace('matr', '').replace('[', '').replace(']', '').split(',')
            i, j = int(i), int(j)

            value = value_part.replace('(', '').replace(')', '').split(',')
            value = float(value[0]) + I * float(value[1])

            # noinspection PyTypeChecker
            matrix[i - 1][j - 1] = value
    return Matrix(matrix)


def analyse_matrices(analytics_matr: Matrix, fortran_matr_path: str) -> None:
    print("Analytics matrix:")
    pprint(analytics_matr)

    # noinspection PyBroadException
    try:
        print("Condition number:")
        print(analytics_matr.condition_number())
    except Exception:
        pass

    # noinspection PyBroadException
    try:
        fortran_matr = load_fortran_matrix(fortran_matr_path, matrix_dim=8 * (N + 1))
    except Exception:
        fortran_matr = zeros(rows=shape(analytics_matr)[0], cols=shape(analytics_matr)[1])

    print("Fortran matrix:")
    pprint(fortran_matr)

    print("Diff:")
    diff = fortran_matr - analytics_matr
    pprint(diff)

    for i in range(shape(diff)[0]):
        for j in range(shape(diff)[1]):
            if abs(diff[i, j]) > 0.9:
                print(f"Delta{(i, j)} = {diff[i, j]}")
