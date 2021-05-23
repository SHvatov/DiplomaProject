from typing import List

import numpy as np
from sympy import I, shape, simplify, expand, Symbol, Expr
from sympy.matrices import Matrix, zeros
from sympy.printing import pprint

from equations import equations
from variables import variables, const_subs

DISCREPANCY_ITER_DELTA = 1.0
APPROXIMATION_DELTA = complex(1.0, 0)
FORTRAN_MATRIX_PATH = "../logs/fortran_logs.txt"


def load_fortran_matrix(fortran_matr_path: str, matrix_dim: int) -> Matrix:
    """
    Loads Jacobian produced by the fortran application into memory as Matrix instance.
    Parameters
    ----------
    fortran_matr_path - path to the txt file, where matrix is stored
    matrix_dim - dimension of the matrix (all matrices are square)

    Returns
    -------
    Jacobian loaded from txt file or zero matrix if file is not found.
    """
    try:
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
    except FileNotFoundError as ex:
        print(f"Could not load fortran matrix: {ex}")
        return zeros(rows=matrix_dim, cols=matrix_dim)


def calculate_discrepancy_matrix(iter_delta: complex) -> Matrix:
    """
    Calculates the Jacobian based on the discrepancy differential.
    Parameters
    ----------
    iter_delta - value of the delta, which will be added to one of the elements on each iteration.

    Returns
    -------
    Jacobian matrix.
    """
    const_sub_eq = [eq.subs(const_subs) for eq in equations]

    def calculate_discrepancy(var: Symbol, zero_approx: List[Expr]) -> List[Expr]:
        left_vars = list(variables)
        left_vars.remove(var)

        sub_eq = [eq.subs({v: 0.0 for v in left_vars}) for eq in const_sub_eq]
        eval_eq = [eq.subs({var: iter_delta}) for eq in sub_eq]

        return [
            simplify(expand((e1 - e0) / iter_delta))
            for e1, e0 in zip(eval_eq, zero_approx)
        ]

    columns = []
    zero_discrepancy = [eq.subs({v: 0.0 for v in variables}) for eq in const_sub_eq]
    for v in variables:
        columns.append(calculate_discrepancy(v, zero_discrepancy))
    return Matrix(columns).T


def calculate_cond_number(matrix: Matrix) -> complex:
    """
    Calculates the condition number of the given matrix.
    Parameters
    ----------
    matrix - matrix, which condition number is required to be calculated.

    Returns
    -------
    Condition number of the matrix calculated using numpy or (0, 0) if something goes wrong.
    """
    # noinspection PyBroadException
    try:
        numpy_matr = np.array(matrix).astype(np.cdouble)
        cond_number = np.linalg.cond(numpy_matr)
        return complex(cond_number.real, cond_number.imag)
    except Exception as e:
        print(f"Cannot calculate condition number: {e}")
        return complex(0, 0)


def compare_matrices(lhs: Matrix, rhs: Matrix, approx_delta: complex) -> None:
    """
    Compares given matrices and prints the difference between them.
    Parameters
    ----------
    lhs - first matrix to compare
    rhs - second matrix to compare
    approx_delta - delta between elements, which considered to be a threshold value.

    Returns
    -------
    None.
    """
    diff = lhs - rhs

    print("Diff:")
    pprint(diff)

    rows, cols = shape(diff)
    for i in range(rows):
        for j in range(cols):
            if abs(diff[i, j]) > approx_delta:
                print(f"Delta{(i, j)} = {diff[i, j]} = {lhs[i, j]} - {rhs[i, j]}")


def analyse_matrices(analytics_matr: Matrix) -> None:
    """
    Performs the analysis of the matrices produced using different approaches.
    Parameters
    ----------
    analytics_matr - analytics Jacobian evaluated using the provided constants.

    Returns
    -------
    None.
    """
    print("Analytics matrix:")
    pprint(analytics_matr)

    analytics_matr_cond_number = calculate_cond_number(analytics_matr)
    print(f"Condition number is {analytics_matr_cond_number}")

    analytics_matr_dim, _ = shape(analytics_matr)
    fortran_matr = load_fortran_matrix(fortran_matr_path=FORTRAN_MATRIX_PATH,
                                       matrix_dim=analytics_matr_dim)

    print("Fortran matrix:")
    pprint(fortran_matr)

    compare_matrices(analytics_matr, fortran_matr, approx_delta=APPROXIMATION_DELTA)

    discrepancy_matr = calculate_discrepancy_matrix(iter_delta=DISCREPANCY_ITER_DELTA)

    print("Discrepancy based matrix:")
    pprint(discrepancy_matr)

    compare_matrices(analytics_matr, discrepancy_matr, approx_delta=APPROXIMATION_DELTA)
