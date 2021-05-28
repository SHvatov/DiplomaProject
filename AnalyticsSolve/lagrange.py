"""
Contains the definition and implementation of the prepare_lagrange_poly.

@author: shvatov
"""
from typing import List, Optional

from sympy import Symbol, Expr, pprint, expand, simplify


def prepare_lagrange_poly(variable: Symbol,
                          points: List[complex],
                          func_values: List[complex]) -> Expr:
    """
    Prepares a Lagrange poly.
    Parameters
    ----------
    variable - symbol which represents the variable
    points - list of points
    func_values - list of the function values in the points

    Returns
    -------
    Expression, which represents the Lagrange Poly.
    """
    assert len(points) == len(func_values)
    assert len(points) > 0

    lagrange_poly: Optional[Expr] = None
    for f_index, func_value in enumerate(func_values):
        lagrange_poly_temp = None
        for p_index, point in enumerate(points):
            if p_index != f_index:
                if lagrange_poly_temp is None:
                    lagrange_poly_temp = (variable - point) / (points[f_index] - point)
                else:
                    lagrange_poly_temp *= (variable - point) / (points[f_index] - point)
        if lagrange_poly is None:
            lagrange_poly = lagrange_poly_temp * func_value
        else:
            lagrange_poly += lagrange_poly_temp * func_value
    return simplify(expand(lagrange_poly))


if __name__ == '__main__':
    pprint(prepare_lagrange_poly(Symbol("r"), [1, 2, 3, 4], [1, 4, 9, 16]))
