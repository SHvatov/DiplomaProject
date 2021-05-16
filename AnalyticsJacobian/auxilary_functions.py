from functions import *


@lru_cache(maxsize=None)
def prepare_a1(index: str) -> Expr:
    assert index is not None

    omega_1, omega_1_conjg = get_omega_func(number="1", index=index)

    ro_13 = get_ro13_func(index)
    ro_13_conjg = get_ro13_func_conjg(index)

    ro_11, ro_22, ro_33, _ = get_main_ro_func(index=index)
    return I * omega_1 * ro_13_conjg - I * omega_1_conjg * ro_13 \
           - gamma_31 * ro_33 + g_parallel * (ro_11 - ro_22)


@lru_cache(maxsize=None)
def prepare_a2(index: str) -> Expr:
    assert index is not None

    omega_2, omega_2_conjg = get_omega_func(number="2", index=index)

    ro_23 = get_ro23_func(index)
    ro_23_conjg = get_ro23_func_conjg(index)

    ro_11, ro_22, ro_33, _ = get_main_ro_func(index=index)
    return I * omega_2 * ro_23_conjg - I * omega_2_conjg * ro_23 \
           - gamma_32 * ro_33 + g_parallel * (ro_22 - ro_11)


@lru_cache(maxsize=None)
def prepare_a3(index: str) -> Expr:
    assert index is not None

    omega_1, omega_1_conjg = get_omega_func(number="1", index=index)
    omega_2, omega_2_conjg = get_omega_func(number="2", index=index)

    ro_13 = get_ro13_func(index)
    ro_13_conjg = get_ro13_func_conjg(index)

    ro_23 = get_ro23_func(index)
    ro_23_conjg = get_ro23_func_conjg(index)

    *_, ro_33, _ = get_main_ro_func(index=index)
    return I * omega_1_conjg * ro_13 - I * omega_1 * ro_13_conjg \
           + I * omega_2_conjg * ro_23 - I * omega_2 * ro_23_conjg \
           + (gamma_31 + gamma_32) * ro_33


@lru_cache(maxsize=None)
def prepare_a4(index: str) -> Expr:
    assert index is not None

    mjr_delta_q = get_major_delta_q_const()
    omega_1, _ = get_omega_func(number="1", index=index)
    _, omega_2_conjg = get_omega_func(number="2", index=index)

    ro_13 = get_ro13_func(index)
    ro_23_conjg = get_ro23_func_conjg(index)

    *_, ro_12 = get_main_ro_func(index=index)
    return mjr_delta_q * ro_12 - I * omega_2_conjg * ro_13 \
           + I * omega_1 * ro_23_conjg


@lru_cache(maxsize=None)
def prepare_a1_conjg(index: str) -> Expr:
    assert index is not None

    omega_1, omega_1_conjg = get_omega_func(number="1", index=index)

    ro_13 = get_ro13_func(index)
    ro_13_conjg = get_ro13_func_conjg(index)

    ro_11_conjg, ro_22_conjg, ro_33_conjg, _ = get_main_ro_func_conjg(index=index)
    return -I * omega_1_conjg * ro_13 + I * omega_1 * ro_13_conjg \
           - gamma_31 * ro_33_conjg + g_parallel * (ro_11_conjg - ro_22_conjg)


@lru_cache(maxsize=None)
def prepare_a2_conjg(index: str) -> Expr:
    assert index is not None

    omega_2, omega_2_conjg = get_omega_func(number="2", index=index)

    ro_23 = get_ro23_func(index)
    ro_23_conjg = get_ro23_func_conjg(index)

    ro_11_conjg, ro_22_conjg, ro_33_conjg, _ = get_main_ro_func_conjg(index=index)
    return -I * omega_2_conjg * ro_23 + I * omega_2 * ro_23_conjg \
           - gamma_32 * ro_33_conjg + g_parallel * (ro_22_conjg - ro_11_conjg)


@lru_cache(maxsize=None)
def prepare_a3_conjg(index: str) -> Expr:
    assert index is not None

    omega_1, omega_1_conjg = get_omega_func(number="1", index=index)
    omega_2, omega_2_conjg = get_omega_func(number="2", index=index)

    ro_13 = get_ro13_func(index)
    ro_13_conjg = get_ro13_func_conjg(index)

    ro_23 = get_ro23_func(index)
    ro_23_conjg = get_ro23_func_conjg(index)

    *_, ro_33_conjg, _ = get_main_ro_func_conjg(index=index)
    return -I * omega_1 * ro_13_conjg + I * omega_1_conjg * ro_13 \
           - I * omega_2 * ro_23_conjg + I * omega_2_conjg * ro_23 \
           + (gamma_31 + gamma_32) * ro_33_conjg


@lru_cache(maxsize=None)
def prepare_a4_conjg(index: str) -> Expr:
    assert index is not None

    mjr_delta_q_conjg = get_major_delta_q_conjg_const()
    _, omega_1_conjg = get_omega_func(number="1", index=index)
    omega_2, _ = get_omega_func(number="2", index=index)

    ro_13_conjg = get_ro13_func_conjg(index)
    ro_23 = get_ro23_func(index)

    *_, ro_12_conjg = get_main_ro_func_conjg(index=index)
    return mjr_delta_q_conjg * ro_12_conjg + I * omega_2 * ro_13_conjg \
           - I * omega_1_conjg * ro_23
