from sympy import * # todo: replace with import of used function
from sympy.abc import r, I, gamma

# Global constants
RO_11 = "Ro_11"
RO_22 = "Ro_22"
RO_33 = "Ro_33"
RO_12 = "Ro_12"

OMEGA_1 = "Omega_1"
OMEGA_2 = "Omega_2"

DELTA_1 = Symbol("delta_1")
DELTA_2 = Symbol("delta_2")

# Ro11 functions
RO_11_FUNCTIONS = []
RO_11_CONJG_FUNCTIONS = []

# Ro22 functions
RO_22_FUNCTIONS = []
RO_22_CONJG_FUNCTIONS = []

# Ro33 functions
RO_33_FUNCTIONS = []
RO_33_CONJG_FUNCTIONS = []

# Ro33 functions
RO_12_FUNCTIONS = []
RO_12_CONJG_FUNCTIONS = []


def prepare_diff_scheme_function(function_name: str, index: str, conjg: bool = False) -> Function:
    return Function(f"ro${function_name}[{index}]{'*' if conjg else ''}")(r)


def prepare_ro13(index: str) -> Function:



if __name__ == '__main__':
    pass