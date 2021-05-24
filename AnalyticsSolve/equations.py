r"""
This module contains the definitions of the equations in the system,
that are created using integral-interpolation method beforehand.

@author: shvatov
"""
from sympy import Expr

from variables import *

left_1 = (
        d_11 * (ri_plus_half[0] * (ro_11[1] - ro_11[0]) / h)
        - 1 / 4 * ri_plus_half[0] * h * (
                I * omega_1[0] * ro_13_conjg[0]
                - I * omega_1_conjg[0] * ro_13[0]
                - gamma_31 * ro_33[0]
                + g_parallel * (ro_11[0] - ro_22[0])
        )
)

left_2 = (
        d_22 * (ri_plus_half[0] * (ro_22[1] - ro_22[0]) / h)
        - 1 / 4 * ri_plus_half[0] * h * (
                I * omega_2[0] * ro_23_conjg[0]
                - I * omega_2_conjg[0] * ro_23[0]
                - gamma_32 * ro_33[0]
                + g_parallel * (ro_22[0] - ro_11[0])
        )
)

left_3 = (
        d_33 * (ri_plus_half[0] * (ro_33[1] - ro_33[0]) / h)
        - 1 / 4 * ri_plus_half[0] * h * (
                I * omega_1_conjg[0] * ro_13[0]
                - I * omega_1[0] * ro_13_conjg[0]
                + I * omega_2_conjg[0] * ro_23[0]
                - I * omega_2[0] * ro_23_conjg[0]
                + (gamma_31 + gamma_32) * ro_33[0]
        )
)

left_4 = (
        d_12 * (ri_plus_half[0] * (ro_12[1] - ro_12[0]) / h)
        - 1 / 4 * ri_plus_half[0] * h * (
                (I * (delta_2 - delta_1) + (g_parallel + q ** 2 * d_12)) * ro_12[0]
                - I * omega_2_conjg[0] * ro_13[0]
                + I * omega_1[0] * ro_23_conjg[0]
        )
)

left_1_conjg = (
        d_11 * (ri_plus_half[0] * (ro_11_conjg[1] - ro_11_conjg[0]) / h)
        - 1 / 4 * ri_plus_half[0] * h * (
                -I * omega_1_conjg[0] * ro_13[0]
                + I * omega_1[0] * ro_13_conjg[0]
                - gamma_31 * ro_33_conjg[0]
                + g_parallel * (ro_11_conjg[0] - ro_22_conjg[0])
        )
)

left_2_conjg = (
        d_22 * (ri_plus_half[0] * (ro_22_conjg[1] - ro_22_conjg[0]) / h)
        - 1 / 4 * ri_plus_half[0] * h * (
                -I * omega_2_conjg[0] * ro_23[0]
                + I * omega_2[0] * ro_23_conjg[0]
                - gamma_32 * ro_33_conjg[0]
                + g_parallel * (ro_22_conjg[0] - ro_11_conjg[0])
        )
)

left_3_conjg = (
        d_33 * (ri_plus_half[0] * (ro_33_conjg[1] - ro_33_conjg[0]) / h)
        - 1 / 4 * ri_plus_half[0] * h * (
                -I * omega_1[0] * ro_13_conjg[0]
                + I * omega_1_conjg[0] * ro_13[0]
                - I * omega_2[0] * ro_23_conjg[0]
                + I * omega_2_conjg[0] * ro_23[0]
                + (gamma_31 + gamma_32) * ro_33_conjg[0]
        )
)

left_4_conjg = (
        d_12 * (ri_plus_half[0] * (ro_12_conjg[1] - ro_12_conjg[0]) / h)
        - 1 / 4 * ri_plus_half[0] * h * (
                (-I * (delta_2 - delta_1) + (g_parallel + q ** 2 * d_12)) * ro_12_conjg[0]
                + I * omega_2[0] * ro_13_conjg[0]
                - I * omega_1_conjg[0] * ro_23[0]
        )
)

main_1 = [
    (
            d_11 * (ri_plus_half[i] * (ro_11[i + 1] - ro_11[i]) / h)
            - d_11 * (ri_minus_half[i] * (ro_11[i] - ro_11[i - 1]) / h)
            - ri[i] * h * (
                    I * omega_1[i] * ro_13_conjg[i]
                    - I * omega_1_conjg[i] * ro_13[i]
                    - gamma_31 * ro_33[i]
                    + g_parallel * (ro_11[i] - ro_22[i])
            )
    )
    for i in range(1, N)
]

main_2 = [
    (
            d_22 * (ri_plus_half[i] * (ro_22[i + 1] - ro_22[i]) / h)
            - d_22 * (ri_minus_half[i] * (ro_22[i] - ro_22[i - 1]) / h)
            - ri[i] * h * (
                    I * omega_2[i] * ro_23_conjg[i]
                    - I * omega_2_conjg[i] * ro_23[i]
                    - gamma_32 * ro_33[i]
                    + g_parallel * (ro_22[i] - ro_11[i])
            )
    )
    for i in range(1, N)
]

main_3 = [
    (
            d_33 * (ri_plus_half[i] * (ro_33[i + 1] - ro_33[i]) / h)
            - d_33 * (ri_minus_half[i] * (ro_33[i] - ro_33[i - 1]) / h)
            - ri[i] * h * (
                    I * omega_1_conjg[i] * ro_13[i]
                    - I * omega_1[i] * ro_13_conjg[i]
                    + I * omega_2_conjg[i] * ro_23[i]
                    - I * omega_2[i] * ro_23_conjg[i]
                    + (gamma_31 + gamma_32) * ro_33[i]
            )
    )
    for i in range(1, N)
]

main_4 = [
    (
            d_12 * (ri_plus_half[i] * (ro_12[i + 1] - ro_12[i]) / h)
            - d_12 * (ri_minus_half[i] * (ro_12[i] - ro_12[i - 1]) / h)
            - ri[i] * h * (
                    (I * (delta_2 - delta_1) + (g_parallel + q ** 2 * d_12)) * ro_12[i]
                    - I * omega_2_conjg[i] * ro_13[i]
                    + I * omega_1[i] * ro_23_conjg[i]
            )
    )
    for i in range(1, N)
]

main_1_conjg = [
    (
            d_11 * (ri_plus_half[i] * (ro_11_conjg[i + 1] - ro_11_conjg[i]) / h)
            - d_11 * (ri_minus_half[i] * (ro_11_conjg[i] - ro_11_conjg[i - 1]) / h)
            - ri[i] * h * (
                    -I * omega_1_conjg[i] * ro_13[i]
                    + I * omega_1[i] * ro_13_conjg[i]
                    - gamma_31 * ro_33_conjg[i]
                    + g_parallel * (ro_11_conjg[i] - ro_22_conjg[i])
            )
    )
    for i in range(1, N)
]

main_2_conjg = [
    (
            d_22 * (ri_plus_half[i] * (ro_22_conjg[i + 1] - ro_22_conjg[i]) / h)
            - d_22 * (ri_minus_half[i] * (ro_22_conjg[i] - ro_22_conjg[i - 1]) / h)
            - ri[i] * h * (
                    -I * omega_2_conjg[i] * ro_23[i]
                    + I * omega_2[i] * ro_23_conjg[i]
                    - gamma_32 * ro_33_conjg[i]
                    + g_parallel * (ro_22_conjg[i] - ro_11_conjg[i])
            )
    )
    for i in range(1, N)
]

main_3_conjg = [
    (
            d_33 * (ri_plus_half[i] * (ro_33_conjg[i + 1] - ro_33_conjg[i]) / h)
            - d_33 * (ri_minus_half[i] * (ro_33_conjg[i] - ro_33_conjg[i - 1]) / h)
            - ri[i] * h * (
                    -I * omega_1[i] * ro_13_conjg[i]
                    + I * omega_1_conjg[i] * ro_13[i]
                    - I * omega_2[i] * ro_23_conjg[i]
                    + I * omega_2_conjg[i] * ro_23[i]
                    + (gamma_31 + gamma_32) * ro_33_conjg[i]
            )
    )
    for i in range(1, N)
]

main_4_conjg = [
    (
            d_12 * (ri_plus_half[i] * (ro_12_conjg[i + 1] - ro_12_conjg[i]) / h)
            - d_12 * (ri_minus_half[i] * (ro_12_conjg[i] - ro_12_conjg[i - 1]) / h)
            - ri[i] * h * (
                    (-I * (delta_2 - delta_1) + (g_parallel + q ** 2 * d_12)) * ro_12_conjg[i]
                    + I * omega_2[i] * ro_13_conjg[i]
                    - I * omega_1_conjg[i] * ro_23[i]
            )
    )
    for i in range(1, N)
]

right_1 = ro_11[N] - 0.5
right_2 = ro_22[N] - 0.5
right_3 = ro_33[N]
right_4 = ro_12[N]

right_1_conjg = ro_11_conjg[N] - 0.5
right_2_conjg = ro_22_conjg[N] - 0.5
right_3_conjg = ro_33_conjg[N]
right_4_conjg = ro_12_conjg[N]


def prepare_ordered_equations() -> List[Expr]:
    """
    Prepares the equations on the specific order, so that Jacobian would have trigonal form.
    Returns
    -------
    List of expressions, that represent the equation system.
    """
    temp_eqs = [
        left_1,
        left_2,
        left_3,
        left_4,
        left_1_conjg,
        left_2_conjg,
        left_3_conjg,
        left_4_conjg,
    ]

    for i in range(0, N - 1):
        temp_eqs.append(main_1[i])
        temp_eqs.append(main_2[i])
        temp_eqs.append(main_3[i])
        temp_eqs.append(main_4[i])
        temp_eqs.append(main_1_conjg[i])
        temp_eqs.append(main_2_conjg[i])
        temp_eqs.append(main_3_conjg[i])
        temp_eqs.append(main_4_conjg[i])

    temp_eqs.append(right_1)
    temp_eqs.append(right_2)
    temp_eqs.append(right_3)
    temp_eqs.append(right_4)
    temp_eqs.append(right_1_conjg)
    temp_eqs.append(right_2_conjg)
    temp_eqs.append(right_3_conjg)
    temp_eqs.append(right_4_conjg)
    return temp_eqs


# Equation system
equations = prepare_ordered_equations()
