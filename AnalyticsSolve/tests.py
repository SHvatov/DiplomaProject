r"""
Module, which can be used as a standalone application, which performs
testing of the solutions on different test functions.

@author: shvatov
"""
from math import sin, cos
from typing import Callable, Dict, List, Tuple

from sympy import Expr, Symbol

from AnalyticsSolve.solve import solve_system, EquationSystemSolutionParams, SolutionMethod
from variables import ri_v

# Dictionary, where key is a name of the function and value is a lambda,
# which produces complex values of that function.
TEST_FUNCTIONS_1: Dict[str, Callable[[int], complex]] = {
    "r11": lambda k: complex(0.5, -1),
    "r22": lambda k: complex(0.5, 1),
    "r33": lambda k: complex(0, 0),
    "r12": lambda k: complex(0, 0),
}

TEST_FUNCTIONS: Dict[str, Callable[[int], complex]] = {
    "r11": lambda k: sin(ri_v(k)) ** 2,
    "r22": lambda k: cos(ri_v(k)) ** 2,
    "r33": lambda k: complex(-0.5, 0),
    "r12": lambda k: complex(0, 0),
}


def approximate_test_solution(equations: List[Expr],
                              variables: List[Symbol],
                              coefficients: Dict[Symbol, complex]) -> Dict[str, complex]:
    # 1. replace constants in the equation system
    print("Step 1. replace constants in the equation system...")
    const_sub_equations = [eq.subs(coefficients) for eq in equations]

    # 2. prepare test function values
    print("Step 2. prepare test function values...")
    test_ro_values = dict()
    for v in variables:
        func = v.name[0:v.name.index('[')]
        t_point = int(v.name[v.name.index('[') + 1:v.name.index(']')])

        try:
            is_conjg = func.index("c") is not None
        except ValueError:
            is_conjg = False

        if is_conjg:
            value = TEST_FUNCTIONS[func[0:len(func) - 1]](t_point).conjugate()
        else:
            value = TEST_FUNCTIONS[func](t_point)

        test_ro_values[v] = value

    # 3. prepare test variables with S(r) source addition
    print("Step 3. prepare test variables with S(r) source addition...")
    s_eq_sources = [eq.subs(test_ro_values) for eq in const_sub_equations]
    test_equations = [eq - s for eq, s in zip(const_sub_equations, s_eq_sources)]

    # 4. solve the system
    print("Step 4. solve the system...")
    params = EquationSystemSolutionParams(method=SolutionMethod.SYMPY)
    solution = solve_system(test_equations, variables, coefficients, params)

    result = dict()
    for i, v in enumerate(solution):
        result[variables[i].name] = v
    return result


if __name__ == '__main__':
    from equations import equations
    from variables import coefficients, variables

    max_diff = 0.0
    solution = approximate_test_solution(equations, variables, coefficients)

    print("Solution:")
    for func_name, sol_value in solution.items():
        f = func_name[0:func_name.index('[')]
        i = int(func_name[func_name.index('[') + 1:func_name.index(']')])

        is_conjg = False
        try:
            f = f[0:f.index("c")]
            is_conjg = True
        except ValueError:
            pass

        test_value = TEST_FUNCTIONS[f](i) if not is_conjg else TEST_FUNCTIONS[f](i).conjugate()
        diff = abs(sol_value - test_value)
        if diff > max_diff:
            max_diff = diff

        print(f"{func_name} = {sol_value}, expected: {test_value}")
        print(f"Diff = {diff}, bigger than 10^-6: {diff > 1e-6}\n")
    print(f"\nMax diff between actual and expected solutions: {max_diff}\n")

    points_by_func: Dict[str, List[Tuple[str, float, complex]]] = dict()
    for f, v in solution.items():
        f_name = f[0:f.index('[')]
        i = int(f[f.index('[') + 1:f.index(']')])
        point = (f, ri_v(i), v)
        if f_name in points_by_func:
            points_by_func[f_name].append(point)
        else:
            points_by_func[f_name] = [point]

    print("\nOutput for further plotting:")
    for f, points in points_by_func.items():
        print(f)
        for point in points:
            name, r, cval = point
            cval = cval.real if abs(cval.imag) < 1e-16 else cval
            print(f"{name} = ({r}, {cval})")
        print()
