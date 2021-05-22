from typing import List, Dict

import numpy as np
from sympy import linear_eq_to_matrix, init_printing, Expr, Symbol, Matrix
from sympy.printing import pprint
from sympy.solvers import solve

from analysis import analyse_matrices
from equations import equations
from fortran_utils import save_jacobian, save_rhs
from variables import const_subs, variables


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

    save_jacobian(A)
    save_rhs(b)

    # noinspection PyBroadException
    try:
        analyse_matrices(A, fortran_matr_path="../logs/fortran_logs.txt")
    except Exception as e:
        print(f"Could not perform matrix analysis: {e}")

    sol_s = solve_sympy(equations, variables)
    print("Sympy solution:")
    for k, v in sol_s.items():
        pprint(f"{k} = {v}")

    sol_n = solve_numpy(A, b)
    print("Numpy solution:")
    for v in sol_n:
        pprint(v)

    print("Solution delta")
    for f, s in zip(sol_s.values(), sol_n):
        diff = abs(f - s)
        print(f"{diff}, bigger than 10^-6: {diff > 1e-6}")


def solve_sympy(eqs: List[Expr], vvars: List[Symbol]) -> Dict[str, complex]:
    sub_equations = [eq.subs(const_subs) for eq in eqs]
    sol = solve(sub_equations, *vvars)

    result = dict()
    for k, v in sol.items():
        result[str(k)] = complex(*v.as_real_imag())
    return result


def solve_numpy(A: Matrix, b: Matrix) -> List[np.cdouble]:
    numpy_A = np.array(A).astype(np.cdouble)
    numpy_B = np.array(b).astype(np.cdouble)
    return np.linalg.solve(numpy_A, numpy_B)


if __name__ == '__main__':
    main()
