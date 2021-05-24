r"""
This module contains the functions, that can be used in order to solve
the original system using either Sympy, or Numpy library.

@author: shvatov
"""
from typing import List, Dict

import numpy as np
from sympy import linear_eq_to_matrix, init_printing, Expr, Symbol, Matrix
from sympy.printing import pprint
from sympy.solvers import solve

from analysis import analyse_matrices
from equations import equations
from variables import const_subs, variables


def solve_sympy(eqs: List[Expr], vvars: List[Symbol]) -> Dict[str, complex]:
    """
    Solves the given system with given variables using sympy.solve.
    Parameters
    ----------
    eqs - list of expressions, that represent the equations system.
    vvars - list of variables, used in the equations system.

    Returns
    -------
    Dictionary, where key is the name of the variable and value is its solution.
    """
    sub_equations = [eq.subs(const_subs) for eq in eqs]
    sol = solve(sub_equations, *vvars)

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
    numpy_A = np.array(A).astype(np.cdouble)
    numpy_B = np.array(b).astype(np.cdouble)
    return np.linalg.solve(numpy_A, numpy_B)


def main() -> None:
    init_printing(use_unicode=False, wrap_line=False)
    print("Constants:")
    pprint(const_subs)

    print("Equations:")
    for eq in equations:
        pprint(eq)

    print("Variables:")
    pprint(variables)

    A, b = linear_eq_to_matrix(equations, *variables)
    print("Matrix A:")
    pprint(A)

    print("Vector b:")
    pprint(b)

    A, b = A.subs(const_subs), b.subs(const_subs)
    print("Evaluated matrix A:")
    pprint(A)

    print("Evaluated vector b:")
    pprint(b)

    analyse_matrices(A)

    sol_s = solve_sympy(equations, variables)
    print("Sympy solution:")
    for k, v in sol_s.items():
        print(f"{k} = {v}")

    sol_n = solve_numpy(A, b)
    print("Numpy solution:")
    for v in sol_n:
        print(v)

    print("Solution delta")
    for f, s in zip(sol_s.values(), sol_n):
        diff = abs(f - s)
        print(f"{diff}, bigger than 10^-6: {diff > 1e-6}")


if __name__ == '__main__':
    main()
