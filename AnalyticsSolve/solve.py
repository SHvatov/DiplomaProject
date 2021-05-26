r"""
This module contains the functions, that can be used in order to solve
the original system using either Sympy, or Numpy library.

@author: shvatov
"""
import sys
from dataclasses import dataclass
from enum import Enum
from typing import List, Dict, Sequence

import numpy as np
from sympy import linear_eq_to_matrix, init_printing, Expr, Symbol, Matrix
from sympy.printing import pprint
from sympy.solvers import solve

from analysis import analyse_matrix, MatrixAnalysisParams
from plot import plot_solution

DELTA = 1e-10


class SolutionMethod(Enum):
    """
    Enum class, which defines what type of the approach will
    be used in order to solve the system. Either use sympy.solvers.solve with
    the analytics equation system, or calculate the Jacobian and solve system
    Ax = b using np.linalg.solve.
    """
    SYMPY = 1
    NUMPY = 2


@dataclass
class EquationSystemSolutionParams:
    """
    Basic data class, which holds the parameters for the
    solve function.
    """
    method: SolutionMethod = SolutionMethod.SYMPY
    verbose_output: bool = False
    check_basic_conditions: bool = False
    n: int = 4
    plot_real_part: bool = False
    analysis_params: MatrixAnalysisParams = None


def solve_sympy(equations: Sequence[Expr],
                variables: Sequence[Symbol],
                coefficients: Dict[Symbol, complex]) -> Dict[str, complex]:
    """
    Solves the given system with given variables using sympy.solve.
    Parameters
    ----------
    equations - list of expressions, that represent the equations system.
    variables - list of variables, used in the equations system.
    coefficients - list of the provided coefficients

    Returns
    -------
    Dictionary, where key is the name of the variable and value is its solution.
    """
    sub_equations = [eq.subs(coefficients) for eq in equations]
    sol = solve(sub_equations, *variables)

    result = dict()
    for k, v in sol.items():
        result[str(k)] = complex(*v.as_real_imag())
    return result


def solve_numpy(A: Matrix, b: Matrix) -> List[np.cdouble]:
    """
    Solves the given equation system represented as Ax = b using numpy.solve.
    Parameters
    ----------
    A - Jacobian matrix of coefficients of the variables.
    b - right-side vector.

    Returns
    -------
    List of complex solutions of the system.
    """
    _A = np.array(A).astype(np.cdouble)
    _b = np.array(b).astype(np.cdouble)
    return np.linalg.solve(_A, _b)


def solve_system(equations: Sequence[Expr],
                 variables: Sequence[Symbol],
                 coefficients: Dict[Symbol, complex],
                 params: EquationSystemSolutionParams = None) -> List[complex]:
    """
    Solves the given system.
    Parameters
    ----------
    equations - list of expressions that represent the system itself.
    variables - list of variables we are looking for.
    coefficients - dictionary, were key is the symbol, which represents the coefficient,
    and value is its value.
    params - function call parameters, see EquationSystemSolutionParams for more info.

    Returns
    -------
    The solution vector.
    """
    if params is None:
        return list()

    if params.verbose_output:
        init_printing(use_unicode=False, wrap_line=False)

        print("Constant coefficients / variables:")
        for coeff, val in coefficients.items():
            print(f"{coeff.name} = {val}")

        print("\nEquations:")
        for eq in equations:
            pprint(eq)

        print("\nVariables:")
        pprint(variables)

    A, b = None, None
    if params.verbose_output \
            or params.analysis_params is not None \
            or params.method == SolutionMethod.NUMPY:
        A, b = linear_eq_to_matrix(equations, *variables)
        if params.verbose_output:
            print("\nMatrix A:")
            pprint(A)

            print("\nVector b:")
            pprint(b)

    if A is not None and b is not None:
        A, b = A.subs(coefficients), b.subs(coefficients)
        if params.verbose_output:
            print("\nEvaluated matrix A:")
            pprint(A)

            print("\nEvaluated vector b:")
            pprint(b)

    if params.analysis_params is not None:
        analyse_matrix(A, params.analysis_params)

    solution = None
    if params.method == SolutionMethod.SYMPY:
        solution = solve_sympy(equations, variables, coefficients)
        if params.verbose_output:
            print("Sympy solution:")
            for k, v in solution.items():
                print(f"{k} = {v}")
        solution = solution.values()
    elif params.method == SolutionMethod.NUMPY:
        solution = solve_numpy(A, b)
        if params.verbose_output:
            print("Numpy solution:")
            for i, v in enumerate(solution):
                print(f"V{i} = {v}")

    if params.check_basic_conditions:
        N = params.n
        value_by_variable = dict()
        for i, v in enumerate(solution):
            value_by_variable[str(variables[i])] = v

        values_by_function = dict()
        for var, v in value_by_variable.items():
            func = var[0:var.index("[")]
            if func not in values_by_function.keys():
                values_by_function[func] = [v]
            else:
                values_by_function[func].append(v)

        print("\nCheck 1: Ro11 + Ro22 + Ro33 = 1")
        max_diff = 0.0
        for i in range(0, N):
            diff = complex(1, 0) \
                   - values_by_function["r11"][i] \
                   - values_by_function["r22"][i] \
                   - values_by_function["r33"][i]
            if diff.real > max_diff:
                max_diff = diff.real
        print(f"Max(1 - (Ro11 + Ro22 + Ro33)) = {max_diff}, "
              f"less than 10^(-10) - {abs(max_diff) < DELTA}")

        print("\nCheck 2: Ro11[N] = Ro22[N] = 0.5, Ro12[N] = Ro33[N] = 0")
        print(f"Ro11[N] = {values_by_function['r11'][N]}, "
              f"diff is {abs(0.5 - values_by_function['r11'][N])}"
              f"less than 10^(-10) - {abs(0.5 - values_by_function['r11'][N]) < DELTA}")
        print(f"Ro22[N] = {values_by_function['r22'][N]}, "
              f"diff is {abs(0.5 - values_by_function['r22'][N])}"
              f"less than 10^(-10) - {abs(0.5 - values_by_function['r22'][N]) < DELTA}")
        print(f"Ro33[N] = {values_by_function['r33'][N]}, "
              f"diff is {abs(values_by_function['r33'][N])}"
              f"less than 10^(-10) - {abs(values_by_function['r33'][N]) < DELTA}")
        print(f"Ro12[N] = {values_by_function['r12'][N]}, "
              f"diff is {abs(values_by_function['r12'][N])}"
              f"less than 10^(-10) - {abs(values_by_function['r12'][N]) < DELTA}")

        print("\nCheck 3: Ro(i, j)[k] = *Ro(i, j)c[k]")
        max_diff = 0.0
        for func in ["r11", "r22", "r33", "r12"]:
            func_c = func + "c"
            for i in range(0, N + 1):
                diff = values_by_function[func][i] - values_by_function[func_c][i].conjugate()
                if diff.real > max_diff.real:
                    max_diff = diff.real
        print(f"Max(Ro(i, j)[k] - *Ro(i, j)c[k]) = {max_diff}, "
              f"less than 10^(-10) - {abs(max_diff) < DELTA}")

    return solution


if __name__ == '__main__':
    from equation import EquationSystem
    from cli import prepare_parser

    parser = prepare_parser()
    args = parser.parse_args(sys.argv[1:])

    params = EquationSystemSolutionParams(method=SolutionMethod(args.method),
                                          verbose_output=args.verbose,
                                          check_basic_conditions=args.check,
                                          n=args.n,
                                          plot_real_part=args.plot,
                                          analysis_params=MatrixAnalysisParams(
                                              print_matrix=args.verbose,
                                              compare_with_fortran=args.fortran,
                                              compare_with_discrepancy=args.discrepancy,
                                              calc_cond_number=args.rcond,
                                          ))
    print(f"Params: \n{params}\n")
    equation_system = EquationSystem.acquire_equation_system(intervals_number=params.n)
    equation_system.ordered_equations()
    solution = solve_system(equations=equation_system.ordered_equations(),
                            variables=equation_system.ordered_variables(),
                            coefficients=equation_system.coefficients(),
                            params=params)

    if params.plot_real_part:
        plot_solution(solution)
