import sympy

from auxilary_functions import *


@lru_cache(maxsize=None)
def prepare_right_boundary_1() -> Expr:
    ro_11 = get_ro_func("11", "N")[0]
    return ro_11 - sympy.Float(0.5)


@lru_cache(maxsize=None)
def prepare_right_boundary_2() -> Expr:
    ro_22 = get_ro_func("22", "N")[0]
    return ro_22 - sympy.Float(0.5)


@lru_cache(maxsize=None)
def prepare_right_boundary_3() -> Expr:
    ro_33 = get_ro_func("33", "N")[0]
    return ro_33


@lru_cache(maxsize=None)
def prepare_right_boundary_4() -> Expr:
    ro_12 = get_ro_func("12", "N")[0]
    return ro_12


@lru_cache(maxsize=None)
def prepare_right_boundary_conjg_1() -> Expr:
    ro_11_conjg = get_ro_func_conjg("11", "N")[0]
    return ro_11_conjg - sympy.Float(0.5)


@lru_cache(maxsize=None)
def prepare_right_boundary_conjg_2() -> Expr:
    ro_22_conjg = get_ro_func_conjg("22", "N")[0]
    return ro_22_conjg - sympy.Float(0.5)


@lru_cache(maxsize=None)
def prepare_right_boundary_conjg_3() -> Expr:
    ro_33_conjg = get_ro_func_conjg("33", "N")[0]
    return ro_33_conjg


@lru_cache(maxsize=None)
def prepare_right_boundary_conjg_4() -> Expr:
    ro_12_conjg = get_ro_func_conjg("12", "N")[0]
    return ro_12_conjg
