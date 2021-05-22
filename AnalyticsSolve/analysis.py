import numpy as np
from sympy import I, shape, simplify, expand
from sympy.matrices import Matrix, zeros
from sympy.printing import pprint

from equations import equations
from variables import N
from variables import variables, const_subs

DISCREPANCY_ITER_DELTA = 1.0


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


def calculate_discrepancy_matrix() -> Matrix:
    const_sub_eq = [eq.subs(const_subs) for eq in equations]
    const_sub_zeros = [eq.subs({v: 0.0 for v in variables}) for eq in const_sub_eq]
    columns = []
    for var in variables:
        left_vars = list(variables)
        left_vars.remove(var)

        sub_eq = [eq.subs({v: 0.0 for v in left_vars}) for eq in const_sub_eq]
        eval_eq = [eq.subs({var: DISCREPANCY_ITER_DELTA}) for eq in sub_eq]

        columns.append(
            [
                simplify(expand((e1 - e0) / DISCREPANCY_ITER_DELTA))
                for e1, e0 in zip(eval_eq, const_sub_zeros)
            ]
        )
    return Matrix(columns).T


def analyse_matrices(analytics_matr: Matrix,
                     fortran_matr_path: str,
                     calculate_jordan: bool = False) -> None:
    print("Analytics matrix:")
    pprint(analytics_matr)

    # noinspection PyBroadException
    try:
        numpy_matr = np.array(analytics_matr).astype(np.cdouble)
        cond_number = np.linalg.cond(numpy_matr)
        print(f"Condition number is {cond_number}")
    except Exception as e:
        print(f"Cannot calculate condition number: {e}")

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
                print(f"Delta{(i, j)} = {diff[i, j]} = {analytics_matr[i, j]} - {fortran_matr[i, j]}")

    if calculate_jordan:
        print("Calculating Jordan forms of both matrices...")
        PA, JA = analytics_matr.jordan_form()
        PF, JF = fortran_matr.jordan_form()

        print("Analytics matrix Jordan form: A = P * J * P^(-1)")
        print("P:")
        pprint(PA)
        print("J:")
        pprint(JA)

        print("Fortran calculated matrix Jordan form: A = P * J * P^(-1)")
        print("P:")
        pprint(PF)
        print("J:")
        pprint(JF)

        print("Deltas:")
        print("Delta P:")
        pprint(PA - PF)
        print("Delta J:")
        pprint(JA - JF)

    print("Discrepancy based matrix:")
    discrepancy_matr = calculate_discrepancy_matrix()
    pprint(discrepancy_matr)

    print("Diff:")
    diff = discrepancy_matr - analytics_matr
    pprint(diff)

    for i in range(shape(diff)[0]):
        for j in range(shape(diff)[1]):
            if abs(diff[i, j]) > 1e-3:
                print(f"Delta{(i, j)} = {diff[i, j]} = {analytics_matr[i, j]} - {discrepancy_matr[i, j]}")
