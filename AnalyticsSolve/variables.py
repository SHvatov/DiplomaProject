r"""
This module contains the definitions of the variables and constants that are
being used in the equation system.

@author: shvatov
"""
from typing import List

from sympy import I
from sympy import Symbol
from sympy.abc import a, q, h
from sympy.functions.elementary.complexes import conjugate
from sympy.functions.elementary.exponential import exp

# size of the mesh
N = 4

# radius of the area
R = 3.3

# Rabi frequencies
C1 = Symbol("C1")
C2 = Symbol("C2")

# r[0], r[1], ..., r[N - 1], r[N]
ri = [Symbol(f"r[{i}]") for i in range(0, N + 1)]

# r[1/2], r[3/2], ..., r[(2 * N - 1) / 2]
ri_plus_half = [Symbol(f"r[{i}/2]") for i in range(1, 2 * N, 2)] + [None]
ri_minus_half = [None] + [Symbol(f"r[{i}/2]") for i in range(1, 2 * N, 2)]

# ro_**[0], ..., ro_**[N]
ro_11 = [Symbol(f"r11[{i}]") for i in range(0, N + 1)]
ro_22 = [Symbol(f"r22[{i}]") for i in range(0, N + 1)]
ro_33 = [Symbol(f"r33[{i}]") for i in range(0, N + 1)]
ro_12 = [Symbol(f"r12[{i}]") for i in range(0, N + 1)]

# conjg(ro_**[0]), ..., conjg(ro_**[N])
# Note: not using conjugated function because this functions
# are also considered to be variables
ro_11_conjg = [Symbol(f"r11c[{i}]") for i in range(0, N + 1)]
ro_22_conjg = [Symbol(f"r22c[{i}]") for i in range(0, N + 1)]
ro_33_conjg = [Symbol(f"r33c[{i}]") for i in range(0, N + 1)]
ro_12_conjg = [Symbol(f"r12c[{i}]") for i in range(0, N + 1)]

# omega_i[0], ..., omega_i[N]
omega_1 = [C1 * exp(-(ri[i] / a) ** 2) for i in range(0, N + 1)]
omega_2 = [C2 * exp(-(ri[i] / a) ** 2) for i in range(0, N + 1)]

# conjg(omega_i[0]), ..., conjg(omega_i*[N])
omega_1_conjg = [conjugate(omega_1[i]) for i in range(0, N + 1)]
omega_2_conjg = [conjugate(omega_2[i]) for i in range(0, N + 1)]

# D11, ..., D12
d_11 = Symbol("D11")
d_22 = Symbol("D22")
d_33 = Symbol("D33")
d_12 = Symbol("D12")

# gammas
gamma = Symbol("Gamma")
gamma_31 = Symbol("Gamma31")
gamma_32 = Symbol("Gamma32")

# deltas
delta_1 = Symbol("Delta1")
delta_2 = Symbol("Delta2")

# other coefficients
g_parallel = Symbol("GParallel")
g_perpendicular = Symbol("GPerpendicular")

# ro_13[0], ...
ro_13 = [
    (I * omega_2[i] * ro_12[i] - I * omega_1[i] * (ro_33[i] - ro_11[i])) / (I * delta_1 + gamma)
    for i in range(0, N + 1)
]

# conjg(ro_13[0]), ...
ro_13_conjg = [
    (
            -I * omega_2_conjg[i] * ro_12_conjg[i]
            + I * omega_1_conjg[i] * (ro_33_conjg[i] - ro_11_conjg[i])
    ) / (-I * delta_1 + gamma)
    for i in range(0, N + 1)
]

# ro_23[0], ...
ro_23 = [
    (I * omega_1[i] * ro_12_conjg[i] - I * omega_1[i] * (ro_33[i] - ro_22[i])) / (I * delta_2 + gamma)
    for i in range(0, N + 1)
]

# conjg(ro_23[0]), ...
ro_23_conjg = [
    (
            -I * omega_1_conjg[i] * ro_12[i]
            + I * omega_1_conjg[i] * (ro_33_conjg[i] - ro_22_conjg[i])
    ) / (-I * delta_2 + gamma)
    for i in range(0, N + 1)
]


def ri_v(i: int) -> float:
    return R / N * i


def rip2_v(i: int) -> float:
    return ri_v(i) + R / (2 * N)


def rim2_v(i: int) -> float:
    return ri_v(i) - R / (2 * N)


const_subs = {
    # Rabi frequencies
    C1: 3e5,
    C2: 3e5,

    # D11, ..., D12
    d_11: 10.0,
    d_22: 10.0,
    d_33: 10.0,
    d_12: 10.0,

    # gammas
    gamma: 2 * 0.875e7 + 1.9825e7,
    gamma_31: 0.875e7,
    gamma_32: 0.875e7,

    # deltas
    delta_1: 0.0,
    delta_2: 0.0,

    # other coefficients
    g_parallel: 5.0e1,
    g_perpendicular: 1.0e2,

    # other consts
    q: 4.29 / 2.9,
    a: 1.4,
    h: R / N,

    # variable r in the scheme
    **{ri[i]: ri_v(i) for i in range(0, N + 1)},
    **{ri_plus_half[i]: rip2_v(i) for i in range(0, N)},
    **{ri_minus_half[i]: rim2_v(i) for i in range(1, N + 1)},
}


def prepare_ordered_variables() -> List[Symbol]:
    """
    Prepares the variables in the specific order, so that Jacobian would have trigonal form.
    Returns
    -------
    List of the variables.
    """
    temp_vars = []
    for i in range(0, N + 1):
        temp_vars.append(ro_11[i])
        temp_vars.append(ro_22[i])
        temp_vars.append(ro_33[i])
        temp_vars.append(ro_12[i])
        temp_vars.append(ro_11_conjg[i])
        temp_vars.append(ro_22_conjg[i])
        temp_vars.append(ro_33_conjg[i])
        temp_vars.append(ro_12_conjg[i])
    return temp_vars


# Ordered ro[i] and ro*[i] variables.
variables = prepare_ordered_variables()
