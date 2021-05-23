from typing import Callable, Dict, List, Tuple

from sympy import init_printing
from sympy.solvers import solve

from equations import equations
from variables import const_subs, variables, ri_v

# Dictionary, where key is a name of the function and value is a lambda,
# which produces complex values of that function.
# TEST_FUNCTIONS: Dict[str, Callable[[int], complex]] = {
#     "r11": lambda k: sin(ri_v(k)) ** 2,
#     "r22": lambda k: cos(ri_v(k)) ** 2,
#     "r33": lambda k: complex(-0.5, 0),
#     "r12": lambda k: complex(0, 0),
# }
# TEST_FUNCTIONS: Dict[str, Callable[[int], complex]] = {
#     "r11": lambda k: ri_v(k) ** 2,
#     "r22": lambda k: ri_v(k) ** 2,
#     "r33": lambda k: 1 - 2 * ri_v(k) ** 2,
#     "r12": lambda k: complex(0, 0),
# }
TEST_FUNCTIONS: Dict[str, Callable[[int], complex]] = {
    "r11": lambda k: complex(0.5, -1),
    "r22": lambda k: complex(0.5, 1),
    "r33": lambda k: complex(0, 0),
    "r12": lambda k: complex(0, 0),
}


def approximate_test_solution() -> Dict[str, complex]:
    # 1. replace constants in the equation system
    print("Step 1. replace constants in the equation system...")
    const_sub_equations = [eq.subs(const_subs) for eq in equations]

    # 2. prepare test function values
    print("Step 2. prepare test function values...")
    test_ro_values = dict()
    for vv in variables:
        func = vv.name[0:vv.name.index('[')]
        t_point = int(vv.name[vv.name.index('[') + 1:vv.name.index(']')])

        try:
            is_conjg = func.index("c") is not None
        except ValueError:
            is_conjg = False

        if is_conjg:
            value = TEST_FUNCTIONS[func[0:len(func) - 1]](t_point).conjugate()
        else:
            value = TEST_FUNCTIONS[func](t_point)

        test_ro_values[vv] = value

    # 3. prepare test variables with S(r) source addition
    print("Step 3. prepare test variables with S(r) source addition...")
    s_eq_sources = [eq.subs(test_ro_values) for eq in const_sub_equations]
    test_equations = [eq - s for eq, s in zip(const_sub_equations, s_eq_sources)]

    # 4. solve the system
    print("Step 4. solve the system...")
    sol = solve(test_equations, *variables)
    result = dict()
    for k, vv in sol.items():
        result[str(k)] = complex(*vv.as_real_imag())
    return result


if __name__ == '__main__':
    init_printing(use_unicode=False, wrap_line=False)

    max_diff = 0.0
    solution = approximate_test_solution()
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
    print(f"Max diff between actual and expected solutions: {max_diff}\n")

    points_by_func: Dict[str, List[Tuple[str, float, complex]]] = dict()
    for f, v in solution.items():
        f_name = f[0:f.index('[')]
        i = int(f[f.index('[') + 1:f.index(']')])
        point = (f, ri_v(i), v)
        if f_name in points_by_func:
            points_by_func[f_name].append(point)
        else:
            points_by_func[f_name] = [point]

    print("Output for further plotting:")
    for f, points in points_by_func.items():
        print(f)
        for point in points:
            name, r, cval = point
            cval = cval.real if abs(cval.imag) < 1e-16 else cval
            print(f"{name} = ({r}, {cval})")
        print()
