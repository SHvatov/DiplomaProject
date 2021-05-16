from auxilary_functions import *


@lru_cache(maxsize=None)
def prepare_left_boundary_1() -> Expr:
    _, ri_plus_half, _ = get_variable("0")
    ro_11, ro_11_plus_one, _ = get_ro_func("11", "0")
    d_11 = get_d_const("11")
    return d_11 * ri_plus_half * (ro_11_plus_one - ro_11) / hr \
           - ri_plus_half * hr / 4 * prepare_a1("0")


@lru_cache(maxsize=None)
def prepare_left_boundary_2() -> Expr:
    _, ri_plus_half, _ = get_variable("0")
    ro_22, ro_22_plus_one, _ = get_ro_func("22", "0")
    d_22 = get_d_const("22")
    return d_22 * ri_plus_half * (ro_22_plus_one - ro_22) / hr \
           - hr / 4 * prepare_a2("0") * ri_plus_half


@lru_cache(maxsize=None)
def prepare_left_boundary_3() -> Expr:
    _, ri_plus_half, _ = get_variable("0")
    ro_33, ro_33_plus_one, _ = get_ro_func("33", "0")
    d_33 = get_d_const("33")
    return d_33 * ri_plus_half * (ro_33_plus_one - ro_33) / hr \
           - hr / 4 * prepare_a3("0") * ri_plus_half


@lru_cache(maxsize=None)
def prepare_left_boundary_4() -> Expr:
    _, ri_plus_half, _ = get_variable("0")
    ro_12, ro_12_plus_one, _ = get_ro_func("12", "0")
    d_12 = get_d_const("12")
    return d_12 * ri_plus_half * (ro_12_plus_one - ro_12) / hr \
           - hr / 4 * prepare_a4("0") * ri_plus_half


@lru_cache(maxsize=None)
def prepare_left_boundary_conjg_1() -> Expr:
    _, ri_plus_half, _ = get_variable("0")
    ro_11_conjg, ro_11_plus_one_conjg, _ = get_ro_func_conjg("11", "0")
    d_11 = get_d_const("11")
    return d_11 * ri_plus_half * (ro_11_plus_one_conjg - ro_11_conjg) / hr \
           - hr / 4 * prepare_a1_conjg("0") * ri_plus_half


@lru_cache(maxsize=None)
def prepare_left_boundary_conjg_2() -> Expr:
    _, ri_plus_half, _ = get_variable("0")
    ro_22_conjg, ro_22_plus_one_conjg, _ = get_ro_func_conjg("22", "0")
    d_22 = get_d_const("22")
    return d_22 * ri_plus_half * (ro_22_plus_one_conjg - ro_22_conjg) / hr \
           - hr / 4 * prepare_a2_conjg("0") * ri_plus_half


@lru_cache(maxsize=None)
def prepare_left_boundary_conjg_3() -> Expr:
    _, ri_plus_half, _ = get_variable("0")
    ro_33_conjg, ro_33_plus_one_conjg, _ = get_ro_func_conjg("33", "0")
    d_33 = get_d_const("33")
    return d_33 * ri_plus_half * (ro_33_plus_one_conjg - ro_33_conjg) / hr \
           - hr / 4 * prepare_a3_conjg("0") * ri_plus_half


@lru_cache(maxsize=None)
def prepare_left_boundary_conjg_4() -> Expr:
    _, ri_plus_half, _ = get_variable("0")
    ro_12_conjg, ro_12_plus_one_conjg, _ = get_ro_func_conjg("12", "0")
    d_12 = get_d_const("12")
    return d_12 * ri_plus_half * (ro_12_plus_one_conjg - ro_12_conjg) / hr \
           - hr / 4 * prepare_a4_conjg("0") * ri_plus_half
