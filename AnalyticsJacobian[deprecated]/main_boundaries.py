from auxilary_functions import *


@lru_cache(maxsize=None)
def prepare_main_boundary_1() -> Expr:
    ri, ri_plus_half, ri_minus_half = get_variable()
    ro_11, ro_11_plus_one, ro_11_minus_one = get_ro_func("11", "i")
    d_11 = get_d_const("11")
    return d_11 * (ri_plus_half * (ro_11_plus_one - ro_11) / hr
                   - ri_minus_half * (ro_11 - ro_11_minus_one) / hr) \
           - hr * ri * prepare_a1("i")


@lru_cache(maxsize=None)
def prepare_main_boundary_2() -> Expr:
    ri, ri_plus_half, ri_minus_half = get_variable()
    ro_22, ro_22_plus_one, ro_22_minus_one = get_ro_func("22", "i")
    d_22 = get_d_const("22")
    return d_22 * (ri_plus_half * (ro_22_plus_one - ro_22) / hr
                   - ri_minus_half * (ro_22 - ro_22_minus_one) / hr) \
           - hr * ri * prepare_a2("i")


@lru_cache(maxsize=None)
def prepare_main_boundary_3() -> Expr:
    ri, ri_plus_half, ri_minus_half = get_variable()
    ro_33, ro_33_plus_one, ro_33_minus_one = get_ro_func("33", "i")
    d_33 = get_d_const("33")
    return d_33 * (ri_plus_half * (ro_33_plus_one - ro_33) / hr
                   - ri_minus_half * (ro_33 - ro_33_minus_one) / hr) \
           - hr * ri * prepare_a3("i")


@lru_cache(maxsize=None)
def prepare_main_boundary_4() -> Expr:
    ri, ri_plus_half, ri_minus_half = get_variable()
    ro_12, ro_12_plus_one, ro_12_minus_one = get_ro_func("12", "i")
    d_12 = get_d_const("12")
    return d_12 * (ri_plus_half * (ro_12_plus_one - ro_12) / hr
                   - ri_minus_half * (ro_12 - ro_12_minus_one) / hr) \
           - hr * ri * prepare_a4("i")


@lru_cache(maxsize=None)
def prepare_main_boundary_conjg_1() -> Expr:
    ri, ri_plus_half, ri_minus_half = get_variable()
    ro_11_conjg, ro_11_plus_one_conjg, ro_11_minus_one_conjg = get_ro_func_conjg("11", "i")
    d_11 = get_d_const("11")
    return d_11 * (ri_plus_half * (ro_11_plus_one_conjg - ro_11_conjg) / hr
                   - ri_minus_half * (ro_11_conjg - ro_11_minus_one_conjg) / hr) \
           - hr * ri * prepare_a1_conjg("i")


@lru_cache(maxsize=None)
def prepare_main_boundary_conjg_2() -> Expr:
    ri, ri_plus_half, ri_minus_half = get_variable()
    ro_22_conjg, ro_22_plus_one_conjg, ro_22_minus_one_conjg = get_ro_func_conjg("22", "i")
    d_22 = get_d_const("22")
    return d_22 * (ri_plus_half * (ro_22_plus_one_conjg - ro_22_conjg) / hr
                   - ri_minus_half * (ro_22_conjg - ro_22_minus_one_conjg) / hr) \
           - hr * ri * prepare_a2_conjg("i")


@lru_cache(maxsize=None)
def prepare_main_boundary_conjg_3() -> Expr:
    ri, ri_plus_half, ri_minus_half = get_variable()
    ro_33_conjg, ro_33_plus_one_conjg, ro_33_minus_one_conjg = get_ro_func_conjg("33", "i")
    d_33 = get_d_const("33")
    return d_33 * (ri_plus_half * (ro_33_plus_one_conjg - ro_33_conjg) / hr
                   - ri_minus_half * (ro_33_conjg - ro_33_minus_one_conjg) / hr) \
           - hr * ri * prepare_a3_conjg("i")


@lru_cache(maxsize=None)
def prepare_main_boundary_conjg_4() -> Expr:
    ri, ri_plus_half, ri_minus_half = get_variable()
    ro_12_conjg, ro_12_plus_one_conjg, ro_12_minus_one_conjg = get_ro_func_conjg("12", "i")
    d_12 = get_d_const("12")
    return d_12 * (ri_plus_half * (ro_12_plus_one_conjg - ro_12_conjg) / hr
                   - ri_minus_half * (ro_12_conjg - ro_12_minus_one_conjg) / hr) \
           - hr * ri * prepare_a3_conjg("i")
