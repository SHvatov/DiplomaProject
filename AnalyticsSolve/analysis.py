r"""
This module contains different functions that can be used to perform the analysis
of the Jacobian computed during the solution of the equations system
and to compare it with Jacobians generated using other methods: by fortran app, or
by differentiating the discrepancy function.

@author: shvatov
"""
from dataclasses import dataclass
from typing import List

import numpy as np
from sympy import I, shape, simplify, expand, Symbol, Expr
from sympy.matrices import Matrix, zeros
from sympy.printing import pprint

from equation import EquationSystem

DISCREPANCY_ITER_DELTA = 1.0
APPROXIMATION_DELTA = complex(1.0, 0)
FORTRAN_MATRIX_PATH = "../logs/fortran_logs.txt"


@dataclass
class MatrixAnalysisParams:
    """
    Basic data class, which holds the parameters for the
    analyse_matrix function.
    """
    print_matrix: bool = False
    compare_with_fortran: bool = False
    compare_with_discrepancy: bool = False
    calc_cond_number: bool = False


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


def calculate_discrepancy_matrix(iter_delta: complex, intervals_num: int) -> Matrix:
    """
    Calculates the Jacobian based on the discrepancy differential.
    Parameters
    ----------
    iter_delta - value of the delta, which will be added to one of the elements on each iteration.

    Returns
    -------
    Jacobian matrix.
    """
    equation_system = EquationSystem.acquire_equation_system(intervals_number=intervals_num)
    equations, variables, coefficients = equation_system.ordered_equations(), \
                                         equation_system.ordered_variables(), \
                                         equation_system.coefficients()
    const_sub_eq = [eq.subs(coefficients) for eq in equations]

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
        cond_number = np.linalg.cond(numpy_matr, p='fro')
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
    delta_present = False
    for i in range(rows):
        for j in range(cols):
            if abs(diff[i, j]) > approx_delta:
                print(f"Delta{(i, j)} = {diff[i, j]} = {lhs[i, j]} - {rhs[i, j]}")
                delta_present = True
    if not delta_present:
        print("No elements surpass specified delta!")


def analyse_matrix(analytics_matr: Matrix,
                   params: MatrixAnalysisParams = None) -> None:
    """
    Performs the analysis of the matrices produced using different approaches.
    Parameters
    ----------
    analytics_matr - analytics Jacobian evaluated using the provided constants.
    params - parameters of the function, see MatrixAnalysisParams for more details.
    If no parameters are provided, then immediately returns from the function.

    Returns
    -------
    None.
    """
    if params is None:
        return

    if params.print_matrix:
        print("Analytics matrix:")
        pprint(analytics_matr)

    if params.calc_cond_number:
        analytics_matr_cond_number = calculate_cond_number(analytics_matr)
        print(f"Condition number is {analytics_matr_cond_number}")

    analytics_matr_dim, _ = shape(analytics_matr)
    if params.compare_with_fortran:
        fortran_matr = load_fortran_matrix(fortran_matr_path=FORTRAN_MATRIX_PATH,
                                           matrix_dim=analytics_matr_dim)
        print("Fortran matrix:")
        pprint(fortran_matr)
        compare_matrices(analytics_matr, fortran_matr, approx_delta=APPROXIMATION_DELTA)

    if params.compare_with_discrepancy:
        intervals = analytics_matr_dim // 8 - 1
        discrepancy_matr = calculate_discrepancy_matrix(iter_delta=DISCREPANCY_ITER_DELTA,
                                                        intervals_num=intervals)

        print("Discrepancy based matrix:")
        pprint(discrepancy_matr)
        compare_matrices(analytics_matr, discrepancy_matr, approx_delta=APPROXIMATION_DELTA)
