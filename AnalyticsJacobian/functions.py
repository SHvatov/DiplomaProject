from functools import lru_cache
from typing import Tuple

from sympy import Expr
from sympy.abc import I

from constants import *


@lru_cache(maxsize=None)
def prepare_symbols(func_name: str,
                    conjg: bool = False,
                    index: str = None,
                    step: str = None) -> Tuple[Symbol, Symbol, Symbol]:
    """
    Prepares a tuple of symbols, that represent function in the net.
    :param func_name: name of the function.
    :param conjg: flag, which determines, whether this function is conjugated or not; if set to true, then * is added.
    :param index: if not none, then used as an index of the function.
    :param step: step to be used in the net
    :return: tuple of sympy symbols, that represent corresponding function.
    """

    def prepare_fun_symbol(idx: str) -> str:
        return f"{func_name}[{idx}]{'_c' if conjg else ''}"

    # noinspection PyBroadException
    try:
        actual_index = int(index)
    except Exception:
        actual_index = index if index is not None else "i"

    # noinspection PyBroadException
    try:
        actual_step = float(eval(step))
    except Exception:
        actual_step = step if step is not None else "1"

    if isinstance(actual_index, (int, float)) and isinstance(actual_step, (int, float)):
        index_plus_one = str(int(actual_index) + float(actual_step))
        index_minus_one = str(int(actual_index) - float(actual_step))
    else:
        index_plus_one = f"{actual_index} + {actual_step}"
        index_minus_one = f"{actual_index} - {actual_step}"

    symbols = Symbol(prepare_fun_symbol(actual_index)), \
              Symbol(prepare_fun_symbol(index_plus_one)), \
              Symbol(prepare_fun_symbol(index_minus_one))

    assert len(symbols) == 3
    return symbols


@lru_cache(maxsize=None)
def get_variable(index: str = None) -> Tuple[Symbol, Symbol, Symbol]:
    """
    Prepares the r[i], r[i + 1/2], r[i - 1/2] variables for the further usage.
    The result of the call is cached, so that at any time same symbol objects are
    returned by this function.
    :return: tuple of sympy symbols, that represent corresponding variables in the described order.
    """
    return prepare_symbols(func_name="r", index=index, step="1/2")


@lru_cache(maxsize=None)
def get_d_const(number: str) -> Symbol:
    assert number in ["11", "22", "33", "12"]
    return Symbol(f"d_{number}")


@lru_cache(maxsize=None)
def get_major_delta_const(number: str) -> Tuple[Expr, Expr]:
    assert number in ["1", "2"]
    gamma_stroke = Symbol("gamma'")
    delta = delta_1 if number == "1" else delta_2
    return I * delta + gamma_stroke, -I * delta + gamma_stroke


@lru_cache(maxsize=None)
def get_major_delta_q_const() -> Expr:
    return I * (delta_2 - delta_1) \
           + (g_parallel + q ** 2 * get_d_const("12"))


@lru_cache(maxsize=None)
def get_major_delta_q_conjg_const() -> Expr:
    return -I * (delta_2 - delta_1) \
           + (g_parallel + q ** 2 * get_d_const("12"))


@lru_cache(maxsize=None)
def get_omega_func(number: str, index: str) -> Tuple[Symbol, Symbol]:
    assert number in ["1", "2"]
    assert index is not None
    func_name = f"omega_{number}"
    return prepare_symbols(func_name=func_name, index=index)[0], \
           prepare_symbols(func_name=func_name, conjg=True, index=index)[0]


@lru_cache(maxsize=None)
def get_ro_func(number: str, index: str = None) -> Tuple:
    assert number in ["11", "22", "33", "12"]
    func_name = f"ro_{number}"
    return prepare_symbols(func_name=func_name, index=index)


@lru_cache(maxsize=None)
def get_ro_func_conjg(number: str, index: str = None) -> Tuple:
    assert number in ["11", "22", "33", "12"]
    func_name = f"ro_{number}"
    return prepare_symbols(func_name=func_name, conjg=True, index=index)


@lru_cache(maxsize=None)
def get_ro13_func(index: str) -> Expr:
    assert index is not None
    ro_12 = get_ro_func(number="12", index=index)[0]
    ro_33 = get_ro_func(number="33", index=index)[0]
    ro_11 = get_ro_func(number="11", index=index)[0]

    omega_2 = get_omega_func(number="2", index=index)[0]
    omega_1 = get_omega_func(number="1", index=index)[0]
    mjr_delta_1 = get_major_delta_const(number="1")[0]

    return I * omega_2 * ro_12 / mjr_delta_1 \
           - I * omega_1 * (ro_33 - ro_11) / mjr_delta_1


@lru_cache(maxsize=None)
def get_ro13_func_conjg(index: str) -> Expr:
    assert index is not None
    ro_12_conjg = get_ro_func_conjg(number="12", index=index)[0]
    ro_33_conjg = get_ro_func_conjg(number="33", index=index)[0]
    ro_11_conjg = get_ro_func_conjg(number="11", index=index)[0]

    omega_2_conjg = get_omega_func(number="2", index=index)[1]
    omega_1_conjg = get_omega_func(number="1", index=index)[1]
    delta_1_conjg = get_major_delta_const(number="1")[1]

    return -I * omega_2_conjg * ro_12_conjg / delta_1_conjg \
           + I * omega_1_conjg * (ro_33_conjg - ro_11_conjg) / delta_1_conjg


@lru_cache(maxsize=None)
def get_ro23_func(index: str) -> Expr:
    assert index is not None
    ro_12_conjg = get_ro_func_conjg(number="12", index=index)[0]
    ro_33 = get_ro_func(number="33", index=index)[0]
    ro_22 = get_ro_func(number="22", index=index)[0]

    omega_1 = get_omega_func(number="1", index=index)[0]
    mjr_delta_2 = get_major_delta_const(number="2")[0]

    return I * omega_1 * ro_12_conjg / mjr_delta_2 \
           - I * omega_1 * (ro_33 - ro_22) / mjr_delta_2


@lru_cache(maxsize=None)
def get_ro23_func_conjg(index: str) -> Expr:
    assert index is not None
    ro_12 = get_ro_func(number="12", index=index)[0]
    ro_33_conjg = get_ro_func_conjg(number="33", index=index)[0]
    ro_22_conjg = get_ro_func_conjg(number="22", index=index)[0]

    omega_1_conjg = get_omega_func(number="1", index=index)[1]
    delta_2_conjg = get_major_delta_const(number="2")[1]

    return -I * omega_1_conjg * ro_12 / delta_2_conjg \
           + I * omega_1_conjg * (ro_33_conjg - ro_22_conjg) / delta_2_conjg


@lru_cache(maxsize=None)
def get_main_ro_func(index: str) -> Tuple[Symbol, Symbol, Symbol, Symbol]:
    ro_11 = get_ro_func(number="11", index=index)[0]
    ro_22 = get_ro_func(number="22", index=index)[0]
    ro_33 = get_ro_func(number="33", index=index)[0]
    ro_12 = get_ro_func(number="12", index=index)[0]
    return ro_11, ro_22, ro_33, ro_12


@lru_cache(maxsize=None)
def get_main_ro_func_conjg(index: str) -> Tuple[Symbol, Symbol, Symbol, Symbol]:
    ro_11_conjg = get_ro_func_conjg(number="11", index=index)[0]
    ro_22_conjg = get_ro_func_conjg(number="22", index=index)[0]
    ro_33_conjg = get_ro_func_conjg(number="33", index=index)[0]
    ro_12_conjg = get_ro_func_conjg(number="12", index=index)[0]
    return ro_11_conjg, ro_22_conjg, ro_33_conjg, ro_12_conjg
