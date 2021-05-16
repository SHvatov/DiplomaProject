from typing import Dict

from sympy import simplify, expand

from left_boundaries import *
from main_boundaries import *
from right_boundaries import *


def print_coefficients(boundary_name: str,
                       boundary_expr: Expr,
                       variables: Dict[str, Symbol]) -> None:
    print(f"{boundary_name}: \n {boundary_expr}")
    expanded_boundary_expr = expand(boundary_expr)

    variables_coeffs = dict()
    for var_name, var in variables.items():
        var_coeff = expanded_boundary_expr.coeff(var)
        print(f"{boundary_name} - {var_name} coefficient: \n {var_coeff}")
        variables_coeffs[var_name] = var_coeff

    independent_var_expr = expanded_boundary_expr
    for var_name, var in variables.items():
        var_coeff = variables_coeffs[var_name]
        independent_var_expr = independent_var_expr - var_coeff * var

    independent_var_expr = simplify(independent_var_expr)
    print(f"{boundary_name} - independent variables: \n {independent_var_expr}")

    print("------------------------------------------")


def print_left_boundaries() -> None:
    print("Left boundaries:")

    ro_11_0, ro_11_1, _ = get_ro_func("11", "0")
    ro_22_0, ro_22_1, _ = get_ro_func("22", "0")
    ro_33_0, ro_33_1, _ = get_ro_func("33", "0")
    ro_12_0, ro_12_1, _ = get_ro_func("12", "0")

    ro_11_0_conjg, ro_11_1_conjg, _ = get_ro_func_conjg("11", "0")
    ro_22_0_conjg, ro_22_1_conjg, _ = get_ro_func_conjg("22", "0")
    ro_33_0_conjg, ro_33_1_conjg, _ = get_ro_func_conjg("33", "0")
    ro_12_0_conjg, ro_12_1_conjg, _ = get_ro_func_conjg("12", "0")

    left_variables = {
        "Ro11[0]": ro_11_0,
        "Ro11[1]": ro_11_1,
        "Ro11*[0]": ro_11_0_conjg,
        "Ro11*[1]": ro_11_1_conjg,

        "Ro22[0]": ro_22_0,
        "Ro22[1]": ro_22_1,
        "Ro22*[0]": ro_22_0_conjg,
        "Ro22*[1]": ro_22_1_conjg,

        "Ro33[0]": ro_33_0,
        "Ro33[1]": ro_33_1,
        "Ro33*[0]": ro_33_0_conjg,
        "Ro33*[1]": ro_33_1_conjg,

        "Ro12[0]": ro_12_0,
        "Ro12[1]": ro_12_1,
        "Ro12*[0]": ro_12_0_conjg,
        "Ro12*[1]": ro_12_1_conjg,
    }

    print("Checked variables:")
    for var_name, var in left_variables.items():
        print(f"{var_name} <=> {var}")
    print("\n")

    boundary_expr = prepare_left_boundary_1()
    print_coefficients(boundary_name="LB(1)",
                       boundary_expr=boundary_expr,
                       variables=left_variables)

    boundary_expr = prepare_left_boundary_2()
    print_coefficients(boundary_name="LB(2)",
                       boundary_expr=boundary_expr,
                       variables=left_variables)

    boundary_expr = prepare_left_boundary_3()
    print_coefficients(boundary_name="LB(3)",
                       boundary_expr=boundary_expr,
                       variables=left_variables)

    boundary_expr = prepare_left_boundary_4()
    print_coefficients(boundary_name="LB(4)",
                       boundary_expr=boundary_expr,
                       variables=left_variables)

    boundary_expr = prepare_left_boundary_conjg_1()
    print_coefficients(boundary_name="LB*(1)",
                       boundary_expr=boundary_expr,
                       variables=left_variables)

    boundary_expr = prepare_left_boundary_conjg_2()
    print_coefficients(boundary_name="LB*(2)",
                       boundary_expr=boundary_expr,
                       variables=left_variables)

    boundary_expr = prepare_left_boundary_conjg_3()
    print_coefficients(boundary_name="LB*(3)",
                       boundary_expr=boundary_expr,
                       variables=left_variables)

    boundary_expr = prepare_left_boundary_conjg_4()
    print_coefficients(boundary_name="LB*(4)",
                       boundary_expr=boundary_expr,
                       variables=left_variables)


def print_right_boundaries() -> None:
    print("Right boundaries:")

    ro_11_n, _, ro_11_n_minus_1 = get_ro_func("11", "N")
    ro_22_n, _, ro_22_n_minus_1 = get_ro_func("22", "N")
    ro_33_n, _, ro_33_n_minus_1 = get_ro_func("33", "N")
    ro_12_n, _, ro_12_n_minus_1 = get_ro_func("12", "N")

    ro_11_n_conjg, _, ro_11_n_minus_1_conjg = get_ro_func_conjg("11", "N")
    ro_22_n_conjg, _, ro_22_n_minus_1_conjg = get_ro_func_conjg("22", "N")
    ro_33_n_conjg, _, ro_33_n_minus_1_conjg = get_ro_func_conjg("33", "N")
    ro_12_n_conjg, _, ro_12_n_minus_1_conjg = get_ro_func_conjg("12", "N")

    right_variables = {
        "Ro11[N]": ro_11_n,
        "Ro11[N-1]": ro_11_n_minus_1,
        "Ro11*[N]": ro_11_n_conjg,
        "Ro11*[N-1]": ro_11_n_minus_1_conjg,

        "Ro22[N]": ro_22_n,
        "Ro22[N-1]": ro_22_n_minus_1,
        "Ro22*[N]": ro_22_n_conjg,
        "Ro22*[N-1]": ro_22_n_minus_1_conjg,

        "Ro33[N]": ro_33_n,
        "Ro33[N-1]": ro_33_n_minus_1,
        "Ro33*[N]": ro_33_n_conjg,
        "Ro33*[N-1]": ro_33_n_minus_1_conjg,

        "Ro12[N]": ro_12_n,
        "Ro12[N-1]": ro_12_n_minus_1,
        "Ro12*[N]": ro_12_n_conjg,
        "Ro12*[N-1]": ro_12_n_minus_1_conjg,
    }

    print("Checked variables:")
    for var_name, var in right_variables.items():
        print(f"{var_name} <=> {var}")
    print("\n")

    boundary_expr = prepare_right_boundary_1()
    print_coefficients(boundary_name="RB(1)",
                       boundary_expr=boundary_expr,
                       variables=right_variables)

    boundary_expr = prepare_right_boundary_2()
    print_coefficients(boundary_name="RB(2)",
                       boundary_expr=boundary_expr,
                       variables=right_variables)

    boundary_expr = prepare_right_boundary_3()
    print_coefficients(boundary_name="RB(3)",
                       boundary_expr=boundary_expr,
                       variables=right_variables)

    boundary_expr = prepare_right_boundary_4()
    print_coefficients(boundary_name="RB(4)",
                       boundary_expr=boundary_expr,
                       variables=right_variables)

    boundary_expr = prepare_right_boundary_conjg_1()
    print_coefficients(boundary_name="RB*(1)",
                       boundary_expr=boundary_expr,
                       variables=right_variables)

    boundary_expr = prepare_right_boundary_conjg_2()
    print_coefficients(boundary_name="RB*(2)",
                       boundary_expr=boundary_expr,
                       variables=right_variables)

    boundary_expr = prepare_right_boundary_conjg_3()
    print_coefficients(boundary_name="RB*(3)",
                       boundary_expr=boundary_expr,
                       variables=right_variables)

    boundary_expr = prepare_right_boundary_conjg_4()
    print_coefficients(boundary_name="RB*(4)",
                       boundary_expr=boundary_expr,
                       variables=right_variables)


def print_main_boundaries() -> None:
    print("Right boundaries:")

    ro_11_i, ro_11_i_plus_1, ro_11_i_minus_1 = get_ro_func("11")
    ro_22_i, ro_22_i_plus_1, ro_22_i_minus_1 = get_ro_func("22")
    ro_33_i, ro_33_i_plus_1, ro_33_i_minus_1 = get_ro_func("33")
    ro_12_i, ro_12_i_plus_1, ro_12_i_minus_1 = get_ro_func("12")

    ro_11_i_conjg, ro_11_i_plus_1_conjg, ro_11_i_minus_1_conjg = get_ro_func_conjg("11")
    ro_22_i_conjg, ro_22_i_plus_1_conjg, ro_22_i_minus_1_conjg = get_ro_func_conjg("22")
    ro_33_i_conjg, ro_33_i_plus_1_conjg, ro_33_i_minus_1_conjg = get_ro_func_conjg("33")
    ro_12_i_conjg, ro_12_i_plus_1_conjg, ro_12_i_minus_1_conjg = get_ro_func_conjg("12")

    main_variables = {
        "Ro11[i]": ro_11_i,
        "Ro11[i+1]": ro_11_i_plus_1,
        "Ro11[i-1]": ro_11_i_minus_1,
        "Ro11*[i]": ro_11_i_conjg,
        "Ro11*[i+1]": ro_11_i_plus_1_conjg,
        "Ro11*[i-1]": ro_11_i_minus_1_conjg,

        "Ro22[i]": ro_22_i,
        "Ro22[i+1]": ro_22_i_plus_1,
        "Ro22[i-1]": ro_22_i_minus_1,
        "Ro22*[i]": ro_22_i_conjg,
        "Ro22*[i+1]": ro_22_i_plus_1_conjg,
        "Ro22*[i-1]": ro_22_i_minus_1_conjg,

        "Ro33[i]": ro_33_i,
        "Ro33[i+1]": ro_33_i_plus_1,
        "Ro33[i-1]": ro_33_i_minus_1,
        "Ro33*[i]": ro_33_i_conjg,
        "Ro33*[i+1]": ro_33_i_plus_1_conjg,
        "Ro33*[i-1]": ro_33_i_minus_1_conjg,

        "Ro12[i]": ro_12_i,
        "Ro12[i+1]": ro_12_i_plus_1,
        "Ro12[i-1]": ro_12_i_minus_1,
        "Ro12*[i]": ro_12_i_conjg,
        "Ro12*[i+1]": ro_12_i_plus_1_conjg,
        "Ro12*[i-1]": ro_12_i_minus_1_conjg,
    }

    print("Checked variables:")
    for var_name, var in main_variables.items():
        print(f"{var_name} <=> {var}")
    print("\n")

    boundary_expr = prepare_main_boundary_1()
    print_coefficients(boundary_name="MAIN(1)",
                       boundary_expr=boundary_expr,
                       variables=main_variables)

    boundary_expr = prepare_main_boundary_2()
    print_coefficients(boundary_name="MAIN(2)",
                       boundary_expr=boundary_expr,
                       variables=main_variables)

    boundary_expr = prepare_main_boundary_3()
    print_coefficients(boundary_name="MAIN(3)",
                       boundary_expr=boundary_expr,
                       variables=main_variables)

    boundary_expr = prepare_main_boundary_4()
    print_coefficients(boundary_name="MAIN(4)",
                       boundary_expr=boundary_expr,
                       variables=main_variables)

    boundary_expr = prepare_main_boundary_conjg_1()
    print_coefficients(boundary_name="MAIN*(1)",
                       boundary_expr=boundary_expr,
                       variables=main_variables)

    boundary_expr = prepare_main_boundary_conjg_2()
    print_coefficients(boundary_name="MAIN*(2)",
                       boundary_expr=boundary_expr,
                       variables=main_variables)

    boundary_expr = prepare_main_boundary_conjg_3()
    print_coefficients(boundary_name="MAIN*(3)",
                       boundary_expr=boundary_expr,
                       variables=main_variables)

    boundary_expr = prepare_main_boundary_conjg_4()
    print_coefficients(boundary_name="MAIN*(4)",
                       boundary_expr=boundary_expr,
                       variables=main_variables)


def main():
    print("******************************************")
    print_left_boundaries()
    print("******************************************")
    print_main_boundaries()
    print("******************************************")
    print_right_boundaries()
    print("******************************************")


if __name__ == '__main__':
    main()
