r"""
This module contains the definition of the class EquationSystem and all related data classes,
that are used to represent the equation system in the application.

@author: shvatov
"""
from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, Callable, Any, Optional, Sequence, Collection

from sympy import Symbol, Expr, I
from sympy.functions.elementary.complexes import conjugate
from sympy.functions.elementary.exponential import exp

# Default radius of the area specified in the task
DEFAULT_RADIUS = 3.3


@dataclass
class EquationSystemCoefficient:
    """
    Defines a coefficient used in the system. Each coefficient is predefined
    and has its own symbol and its complex value.
    """
    symbol: Symbol
    value: complex

    def __hash__(self) -> int:
        return hash(self.symbol) * 31 + hash(self.value)

    @staticmethod
    def create(symbol: str, value: complex) -> Any:
        return EquationSystemCoefficient(Symbol(symbol), value)


@dataclass
class EquationSystemVariable:
    """
    Defines a variable of the equation system. Due to the usage of the integro-interpolation
    method, it is divided into multiple variables. Values are optional.
    """
    __symbols: Sequence[Optional[Symbol]]
    __values: Dict[Symbol, complex]

    def __getitem__(self, index: int) -> Symbol:
        assert index in range(len(self.__symbols))
        assert self.__symbols[index] is not None
        return self.__symbols[index]

    def values(self) -> Dict[Symbol, complex]:
        return self.__values

    def __hash__(self) -> int:
        return hash(self.__symbols) * 31 + hash(self.__values)

    @staticmethod
    def from_pattern(symbol_pattern: str,
                     elements_number: int,
                     value_generator: Callable[[int], complex] = None) -> Any:
        assert len(symbol_pattern) > 0
        assert elements_number > 0

        symbols = [Symbol(symbol_pattern.format(i)) for i in range(elements_number)]
        values = {symbols[i]: value_generator(i) for i in range(len(symbols))} if value_generator is not None else []
        return EquationSystemVariable.from_provided(symbols, values)

    @staticmethod
    def from_provided(symbols: Sequence[Optional[Symbol]], values: Dict[Symbol, complex] = None) -> Any:
        if values is None:
            values = dict()
        return EquationSystemVariable(symbols, values)


@dataclass(init=False)
class EquationSystemExpression:
    """
    Defines an expression in the equation system. Due to the usage of the integro-interpolation
    method, it is divided into multiple variables. Values are optional.
    """
    __expressions: Sequence[Expr]

    def __init__(self, expressions: Collection[Expr]):
        self.__expressions = tuple(expressions)

    def __getitem__(self, index: int) -> Expr:
        assert index in range(len(self.__expressions))
        return self.__expressions[index]

    def __hash__(self) -> int:
        return hash(self.__expressions)


@dataclass(init=False)
class EquationSystem:
    # Number of the intervals in the scheme
    N: int

    # radius of the area
    R: float

    # Size of the step in the mesh
    h: EquationSystemCoefficient

    # r[i] variables
    ri: EquationSystemVariable
    ri_plus_half: EquationSystemVariable
    ri_minus_half: EquationSystemVariable

    # ro[i] functions
    ro_11: EquationSystemVariable
    ro_22: EquationSystemVariable
    ro_33: EquationSystemVariable
    ro_12: EquationSystemVariable

    # ro[i]* functions
    ro_11_conjg: EquationSystemVariable
    ro_22_conjg: EquationSystemVariable
    ro_33_conjg: EquationSystemVariable
    ro_12_conjg: EquationSystemVariable

    # omega_i[0], ..., omega_i[N]
    omega_1: EquationSystemExpression
    omega_2: EquationSystemExpression

    # conjg(omega_i[0]), ..., conjg(omega_i*[N])
    omega_1_conjg: EquationSystemExpression
    omega_2_conjg: EquationSystemExpression

    # ro_13, ro_13*
    ro_13: EquationSystemExpression
    ro_13_conjg: EquationSystemExpression

    # ro_23, ro_23*
    ro_23: EquationSystemExpression
    ro_23_conjg: EquationSystemExpression

    # Left boundaries
    left_1: EquationSystemExpression
    left_2: EquationSystemExpression
    left_3: EquationSystemExpression
    left_4: EquationSystemExpression

    # Left boundaries - conjugated
    left_1_conjg: EquationSystemExpression
    left_2_conjg: EquationSystemExpression
    left_3_conjg: EquationSystemExpression
    left_4_conjg: EquationSystemExpression

    # Main equation system
    man_1: EquationSystemExpression
    man_2: EquationSystemExpression
    man_3: EquationSystemExpression
    man_4: EquationSystemExpression

    # Main equation system - conjugated
    man_1_conjg: EquationSystemExpression
    man_2_conjg: EquationSystemExpression
    man_3_conjg: EquationSystemExpression
    man_4_conjg: EquationSystemExpression

    # Right boundaries
    right_1: EquationSystemExpression
    right_2: EquationSystemExpression
    right_3: EquationSystemExpression
    right_4: EquationSystemExpression

    # Right boundaries - conjugated
    right_1_conjg: EquationSystemExpression
    right_2_conjg: EquationSystemExpression
    right_3_conjg: EquationSystemExpression
    right_4_conjg: EquationSystemExpression

    # Rabi frequencies
    C1: EquationSystemCoefficient = EquationSystemCoefficient.create(symbol="C1", value=3e5)
    C2: EquationSystemCoefficient = EquationSystemCoefficient.create(symbol="C2", value=3e5)

    # D11, ..., D12
    d_11: EquationSystemCoefficient = EquationSystemCoefficient.create(symbol="D11", value=10.0)
    d_22: EquationSystemCoefficient = EquationSystemCoefficient.create(symbol="D22", value=10.0)
    d_33: EquationSystemCoefficient = EquationSystemCoefficient.create(symbol="D33", value=10.0)
    d_12: EquationSystemCoefficient = EquationSystemCoefficient.create(symbol="D12", value=10.0)

    # Gamma values
    gamma: EquationSystemCoefficient = EquationSystemCoefficient.create(symbol="Gamma", value=2 * 0.875e7 + 1.9825e7)
    gamma_31: EquationSystemCoefficient = EquationSystemCoefficient.create(symbol="Gamma31", value=0.875e7)
    gamma_32: EquationSystemCoefficient = EquationSystemCoefficient.create(symbol="Gamma32", value=0.875e7)

    # Delta values
    delta_1: EquationSystemCoefficient = EquationSystemCoefficient.create(symbol="Delta1", value=0.0)
    delta_2: EquationSystemCoefficient = EquationSystemCoefficient.create(symbol="Delta2", value=0.0)

    # G coefficients
    g_parallel: EquationSystemCoefficient = EquationSystemCoefficient.create(symbol="GParallel", value=5.0e1)
    g_perpendicular: EquationSystemCoefficient = EquationSystemCoefficient.create(symbol="GPerpendicular", value=1.0e2)

    # other coefficients
    a: EquationSystemCoefficient = EquationSystemCoefficient.create(symbol="a", value=1.4)
    q: EquationSystemCoefficient = EquationSystemCoefficient.create(symbol="q", value=4.29 / 2.9)

    def __init__(self, intervals_number: int, radius: float):
        assert intervals_number > 0
        assert radius > 0

        self.N, self.R = intervals_number, radius
        self.h = EquationSystemCoefficient.create(symbol="h", value=self.R / self.N)

        ri_v_gen = lambda i: self.R / self.N * i
        self.ri = EquationSystemVariable.from_pattern(symbol_pattern="r[{}]",
                                                      elements_number=self.N + 1,
                                                      value_generator=ri_v_gen)

        ri_plus_half_v_gen = lambda i: self.R / self.N * i + self.R / (2 * self.N)
        ri_plus_half_sym = [Symbol(f"r[{i}/2]") for i in range(1, 2 * self.N, 2)] + [None]
        ri_plus_half_val = {ri_plus_half_sym[i]: ri_plus_half_v_gen(i) for i in range(0, self.N)}
        self.ri_plus_half = EquationSystemVariable.from_provided(symbols=ri_plus_half_sym,
                                                                 values=ri_plus_half_val)

        ri_minus_half_v_gen = lambda i: self.R / self.N * i - self.R / (2 * self.N)
        ri_minus_half_sym = [None] + [Symbol(f"r[{i}/2]") for i in range(1, 2 * self.N, 2)]
        ri_minus_half_val = {ri_minus_half_sym[i]: ri_minus_half_v_gen(i) for i in range(1, self.N + 1)}
        self.ri_minus_half = EquationSystemVariable.from_provided(symbols=ri_minus_half_sym,
                                                                  values=ri_minus_half_val)

        self.ro_11 = EquationSystemVariable.from_pattern(symbol_pattern="r11[{}]",
                                                         elements_number=self.N + 1)
        self.ro_22 = EquationSystemVariable.from_pattern(symbol_pattern="r22[{}]",
                                                         elements_number=self.N + 1)
        self.ro_33 = EquationSystemVariable.from_pattern(symbol_pattern="r33[{}]",
                                                         elements_number=self.N + 1)
        self.ro_12 = EquationSystemVariable.from_pattern(symbol_pattern="r12[{}]",
                                                         elements_number=self.N + 1)

        self.ro_11_conjg = EquationSystemVariable.from_pattern(symbol_pattern="r11c[{}]",
                                                               elements_number=self.N + 1)
        self.ro_22_conjg = EquationSystemVariable.from_pattern(symbol_pattern="r22c[{}]",
                                                               elements_number=self.N + 1)
        self.ro_33_conjg = EquationSystemVariable.from_pattern(symbol_pattern="r33c[{}]",
                                                               elements_number=self.N + 1)
        self.ro_12_conjg = EquationSystemVariable.from_pattern(symbol_pattern="r12c[{}]",
                                                               elements_number=self.N + 1)

        self.omega_1 = EquationSystemExpression(
            [self.C1.symbol * exp(-(self.ri[i] / self.a.symbol) ** 2) for i in range(0, self.N + 1)]
        )
        self.omega_2 = EquationSystemExpression(
            [self.C2.symbol * exp(-(self.ri[i] / self.a.symbol) ** 2) for i in range(0, self.N + 1)]
        )

        self.omega_1_conjg = EquationSystemExpression(
            [conjugate(self.omega_1[i]) for i in range(0, self.N + 1)]
        )
        self.omega_2_conjg = EquationSystemExpression(
            [conjugate(self.omega_2[i]) for i in range(0, self.N + 1)]
        )

        self.ro_13 = EquationSystemExpression(
            [
                (I * self.omega_2[i] * self.ro_12[i] - I * self.omega_1[i] *
                 (self.ro_33[i] - self.ro_11[i])) / (I * self.delta_1.symbol + self.gamma.symbol)
                for i in range(0, self.N + 1)
            ]
        )

        self.ro_13_conjg = EquationSystemExpression(
            [
                (-I * self.omega_2_conjg[i] * self.ro_12_conjg[i] + I * self.omega_1_conjg[i] *
                 (self.ro_33_conjg[i] - self.ro_11_conjg[i])) / (-I * self.delta_1.symbol + self.gamma.symbol)
                for i in range(0, self.N + 1)
            ]
        )

        self.ro_23 = EquationSystemExpression(
            [
                (I * self.omega_1[i] * self.ro_12_conjg[i] - I * self.omega_1[i] *
                 (self.ro_33[i] - self.ro_22[i])) / (I * self.delta_2.symbol + self.gamma.symbol)
                for i in range(0, self.N + 1)
            ]
        )

        self.ro_23_conjg = EquationSystemExpression(
            [
                (-I * self.omega_1_conjg[i] * self.ro_12[i] + I * self.omega_1_conjg[i] *
                 (self.ro_33_conjg[i] - self.ro_22_conjg[i])) / (-I * self.delta_2.symbol + self.gamma.symbol)
                for i in range(0, self.N + 1)
            ]
        )

        self.left_1 = EquationSystemExpression(
            [
                (
                        self.d_11.symbol * (self.ri_plus_half[0] * (self.ro_11[1] - self.ro_11[0]) / self.h.symbol)
                        - 1 / 4 * self.ri_plus_half[0] * self.h.symbol * (
                                I * self.omega_1[0] * self.ro_13_conjg[0]
                                - I * self.omega_1_conjg[0] * self.ro_13[0]
                                - self.gamma_31.symbol * self.ro_33[0]
                                + self.g_parallel.symbol * (self.ro_11[0] - self.ro_22[0])
                        )
                )
            ]
        )

        self.left_2 = EquationSystemExpression(
            [
                (
                        self.d_22.symbol * (self.ri_plus_half[0] * (self.ro_22[1] - self.ro_22[0]) / self.h.symbol)
                        - 1 / 4 * self.ri_plus_half[0] * self.h.symbol * (
                                I * self.omega_2[0] * self.ro_23_conjg[0]
                                - I * self.omega_2_conjg[0] * self.ro_23[0]
                                - self.gamma_32.symbol * self.ro_33[0]
                                + self.g_parallel.symbol * (self.ro_22[0] - self.ro_11[0])
                        )
                )
            ]
        )

        self.left_3 = EquationSystemExpression(
            [
                (
                        self.d_33.symbol * (self.ri_plus_half[0] * (self.ro_33[1] - self.ro_33[0]) / self.h.symbol)
                        - 1 / 4 * self.ri_plus_half[0] * self.h.symbol * (
                                I * self.omega_1_conjg[0] * self.ro_13[0]
                                - I * self.omega_1[0] * self.ro_13_conjg[0]
                                + I * self.omega_2_conjg[0] * self.ro_23[0]
                                - I * self.omega_2[0] * self.ro_23_conjg[0]
                                + (self.gamma_31.symbol + self.gamma_32.symbol) * self.ro_33[0]
                        )
                )
            ]
        )

        self.left_4 = EquationSystemExpression(
            [
                (
                        self.d_12.symbol * (self.ri_plus_half[0] * (self.ro_12[1] - self.ro_12[0]) / self.h.symbol)
                        - 1 / 4 * self.ri_plus_half[0] * self.h.symbol * (
                                (I * (self.delta_2.symbol - self.delta_1.symbol) + (
                                        self.g_perpendicular.symbol + self.q.symbol ** 2 * self.d_12.symbol)) *
                                self.ro_12[0]
                                - I * self.omega_2_conjg[0] * self.ro_13[0]
                                + I * self.omega_1[0] * self.ro_23_conjg[0]
                        )
                )
            ]
        )

        self.left_1_conjg = EquationSystemExpression(
            [
                (
                        self.d_11.symbol
                        * (self.ri_plus_half[0] * (self.ro_11_conjg[1] - self.ro_11_conjg[0]) / self.h.symbol)
                        - 1 / 4 * self.ri_plus_half[0] * self.h.symbol * (
                                -I * self.omega_1_conjg[0] * self.ro_13[0]
                                + I * self.omega_1[0] * self.ro_13_conjg[0]
                                - self.gamma_31.symbol * self.ro_33_conjg[0]
                                + self.g_parallel.symbol * (self.ro_11_conjg[0] - self.ro_22_conjg[0])
                        )
                )
            ]
        )

        self.left_2_conjg = EquationSystemExpression(
            [
                (
                        self.d_22.symbol *
                        (self.ri_plus_half[0] * (self.ro_22_conjg[1] - self.ro_22_conjg[0]) / self.h.symbol)
                        - 1 / 4 * self.ri_plus_half[0] * self.h.symbol * (
                                -I * self.omega_2_conjg[0] * self.ro_23[0]
                                + I * self.omega_2[0] * self.ro_23_conjg[0]
                                - self.gamma_32.symbol * self.ro_33_conjg[0]
                                + self.g_parallel.symbol * (self.ro_22_conjg[0] - self.ro_11_conjg[0])
                        )
                )
            ]
        )

        self.left_3_conjg = EquationSystemExpression(
            [
                (
                        self.d_33.symbol
                        * (self.ri_plus_half[0] * (self.ro_33_conjg[1] - self.ro_33_conjg[0]) / self.h.symbol)
                        - 1 / 4 * self.ri_plus_half[0] * self.h.symbol * (
                                -I * self.omega_1[0] * self.ro_13_conjg[0]
                                + I * self.omega_1_conjg[0] * self.ro_13[0]
                                - I * self.omega_2[0] * self.ro_23_conjg[0]
                                + I * self.omega_2_conjg[0] * self.ro_23[0]
                                + (self.gamma_31.symbol + self.gamma_32.symbol) * self.ro_33_conjg[0]
                        )
                )
            ]
        )

        self.left_4_conjg = EquationSystemExpression(
            [
                (
                        self.d_12.symbol
                        * (self.ri_plus_half[0] * (self.ro_12_conjg[1] - self.ro_12_conjg[0]) / self.h.symbol)
                        - 1 / 4 * self.ri_plus_half[0] * self.h.symbol * (
                                (-I * (self.delta_2.symbol - self.delta_1.symbol)
                                 + (self.g_perpendicular.symbol + self.q.symbol ** 2 * self.d_12.symbol))
                                * self.ro_12_conjg[0]
                                + I * self.omega_2[0] * self.ro_13_conjg[0]
                                - I * self.omega_1_conjg[0] * self.ro_23[0]
                        )
                )
            ]
        )

        self.main_1 = EquationSystemExpression(
            [
                (
                        self.d_11.symbol
                        * (self.ri_plus_half[i] * (self.ro_11[i + 1] - self.ro_11[i]) / self.h.symbol)
                        - self.d_11.symbol
                        * (self.ri_minus_half[i] * (self.ro_11[i] - self.ro_11[i - 1]) / self.h.symbol)
                        - self.ri[i] * self.h.symbol * (
                                I * self.omega_1[i] * self.ro_13_conjg[i]
                                - I * self.omega_1_conjg[i] * self.ro_13[i]
                                - self.gamma_31.symbol * self.ro_33[i]
                                + self.g_parallel.symbol * (self.ro_11[i] - self.ro_22[i])
                        )
                )
                for i in range(1, self.N)
            ]
        )

        self.main_2 = EquationSystemExpression(
            [
                (
                        self.d_22.symbol
                        * (self.ri_plus_half[i] * (self.ro_22[i + 1] - self.ro_22[i]) / self.h.symbol)
                        - self.d_22.symbol
                        * (self.ri_minus_half[i] * (self.ro_22[i] - self.ro_22[i - 1]) / self.h.symbol)
                        - self.ri[i] * self.h.symbol * (
                                I * self.omega_2[i] * self.ro_23_conjg[i]
                                - I * self.omega_2_conjg[i] * self.ro_23[i]
                                - self.gamma_32.symbol * self.ro_33[i]
                                + self.g_parallel.symbol * (self.ro_22[i] - self.ro_11[i])
                        )
                )
                for i in range(1, self.N)
            ]
        )

        self.main_3 = EquationSystemExpression(
            [
                (
                        self.d_33.symbol
                        * (self.ri_plus_half[i] * (self.ro_33[i + 1] - self.ro_33[i]) / self.h.symbol)
                        - self.d_33.symbol
                        * (self.ri_minus_half[i] * (self.ro_33[i] - self.ro_33[i - 1]) / self.h.symbol)
                        - self.ri[i] * self.h.symbol * (
                                I * self.omega_1_conjg[i] * self.ro_13[i]
                                - I * self.omega_1[i] * self.ro_13_conjg[i]
                                + I * self.omega_2_conjg[i] * self.ro_23[i]
                                - I * self.omega_2[i] * self.ro_23_conjg[i]
                                + (self.gamma_31.symbol + self.gamma_32.symbol) * self.ro_33[i]
                        )
                )
                for i in range(1, self.N)
            ]
        )

        self.main_4 = EquationSystemExpression([
            (
                    self.d_12.symbol
                    * (self.ri_plus_half[i] * (self.ro_12[i + 1] - self.ro_12[i]) / self.h.symbol)
                    - self.d_12.symbol
                    * (self.ri_minus_half[i] * (self.ro_12[i] - self.ro_12[i - 1]) / self.h.symbol)
                    - self.ri[i] * self.h.symbol * (
                            (I * (self.delta_2.symbol - self.delta_1.symbol)
                             + (self.g_perpendicular.symbol + self.q.symbol ** 2 * self.d_12.symbol))
                            * self.ro_12[i]
                            - I * self.omega_2_conjg[i] * self.ro_13[i]
                            + I * self.omega_1[i] * self.ro_23_conjg[i]
                    )
            )
            for i in range(1, self.N)
        ])

        self.main_1_conjg = EquationSystemExpression([
            (
                    self.d_11.symbol
                    * (self.ri_plus_half[i] * (self.ro_11_conjg[i + 1] - self.ro_11_conjg[i]) / self.h.symbol)
                    - self.d_11.symbol *
                    (self.ri_minus_half[i] * (self.ro_11_conjg[i] - self.ro_11_conjg[i - 1]) / self.h.symbol)
                    - self.ri[i] * self.h.symbol * (
                            -I * self.omega_1_conjg[i] * self.ro_13[i]
                            + I * self.omega_1[i] * self.ro_13_conjg[i]
                            - self.gamma_31.symbol * self.ro_33_conjg[i]
                            + self.g_parallel.symbol * (self.ro_11_conjg[i] - self.ro_22_conjg[i])
                    )
            )
            for i in range(1, self.N)
        ])

        self.main_2_conjg = EquationSystemExpression([
            (
                    self.d_22.symbol
                    * (self.ri_plus_half[i] * (self.ro_22_conjg[i + 1] - self.ro_22_conjg[i]) / self.h.symbol)
                    - self.d_22.symbol
                    * (self.ri_minus_half[i] * (self.ro_22_conjg[i] - self.ro_22_conjg[i - 1]) / self.h.symbol)
                    - self.ri[i] * self.h.symbol * (
                            -I * self.omega_2_conjg[i] * self.ro_23[i]
                            + I * self.omega_2[i] * self.ro_23_conjg[i]
                            - self.gamma_32.symbol * self.ro_33_conjg[i]
                            + self.g_parallel.symbol * (self.ro_22_conjg[i] - self.ro_11_conjg[i])
                    )
            )
            for i in range(1, self.N)
        ])

        self.main_3_conjg = EquationSystemExpression([
            (
                    self.d_33.symbol
                    * (self.ri_plus_half[i] * (self.ro_33_conjg[i + 1] - self.ro_33_conjg[i]) / self.h.symbol)
                    - self.d_33.symbol
                    * (self.ri_minus_half[i] * (self.ro_33_conjg[i] - self.ro_33_conjg[i - 1]) / self.h.symbol)
                    - self.ri[i] * self.h.symbol * (
                            -I * self.omega_1[i] * self.ro_13_conjg[i]
                            + I * self.omega_1_conjg[i] * self.ro_13[i]
                            - I * self.omega_2[i] * self.ro_23_conjg[i]
                            + I * self.omega_2_conjg[i] * self.ro_23[i]
                            + (self.gamma_31.symbol + self.gamma_32.symbol) * self.ro_33_conjg[i]
                    )
            )
            for i in range(1, self.N)
        ])

        self.main_4_conjg = EquationSystemExpression([
            (
                    self.d_12.symbol
                    * (self.ri_plus_half[i] * (self.ro_12_conjg[i + 1] - self.ro_12_conjg[i]) / self.h.symbol)
                    - self.d_12.symbol
                    * (self.ri_minus_half[i] * (self.ro_12_conjg[i] - self.ro_12_conjg[i - 1]) / self.h.symbol)
                    - self.ri[i] * self.h.symbol * (
                            (-I * (self.delta_2.symbol - self.delta_1.symbol)
                             + (self.g_perpendicular.symbol + self.q.symbol ** 2 * self.d_12.symbol))
                            * self.ro_12_conjg[i]
                            + I * self.omega_2[i] * self.ro_13_conjg[i]
                            - I * self.omega_1_conjg[i] * self.ro_23[i]
                    )
            )
            for i in range(1, self.N)
        ])

        self.right_1 = EquationSystemExpression([self.ro_11[self.N] - 0.5])
        self.right_2 = EquationSystemExpression([self.ro_22[self.N] - 0.5])
        self.right_3 = EquationSystemExpression([self.ro_33[self.N]])
        self.right_4 = EquationSystemExpression([self.ro_12[self.N]])

        self.right_1_conjg = EquationSystemExpression([self.ro_11_conjg[self.N] - 0.5])
        self.right_2_conjg = EquationSystemExpression([self.ro_22_conjg[self.N] - 0.5])
        self.right_3_conjg = EquationSystemExpression([self.ro_33_conjg[self.N]])
        self.right_4_conjg = EquationSystemExpression([self.ro_12_conjg[self.N]])

    @lru_cache(maxsize=None)
    def ordered_variables(self) -> Sequence[Symbol]:
        variables = []
        for i in range(0, self.N + 1):
            variables.append(self.ro_11[i])
            variables.append(self.ro_22[i])
            variables.append(self.ro_33[i])
            variables.append(self.ro_12[i])
            variables.append(self.ro_11_conjg[i])
            variables.append(self.ro_22_conjg[i])
            variables.append(self.ro_33_conjg[i])
            variables.append(self.ro_12_conjg[i])
        return variables

    @lru_cache(maxsize=None)
    def ordered_equations(self) -> Sequence[Expr]:
        temp_eqs = [
            self.left_1[0],
            self.left_2[0],
            self.left_3[0],
            self.left_4[0],
            self.left_1_conjg[0],
            self.left_2_conjg[0],
            self.left_3_conjg[0],
            self.left_4_conjg[0]
        ]

        for i in range(0, self.N - 1):
            temp_eqs.append(self.main_1[i])
            temp_eqs.append(self.main_2[i])
            temp_eqs.append(self.main_3[i])
            temp_eqs.append(self.main_4[i])
            temp_eqs.append(self.main_1_conjg[i])
            temp_eqs.append(self.main_2_conjg[i])
            temp_eqs.append(self.main_3_conjg[i])
            temp_eqs.append(self.main_4_conjg[i])

        temp_eqs.append(self.right_1[0])
        temp_eqs.append(self.right_2[0])
        temp_eqs.append(self.right_3[0])
        temp_eqs.append(self.right_4[0])
        temp_eqs.append(self.right_1_conjg[0])
        temp_eqs.append(self.right_2_conjg[0])
        temp_eqs.append(self.right_3_conjg[0])
        temp_eqs.append(self.right_4_conjg[0])
        return temp_eqs

    @lru_cache(maxsize=None)
    def coefficients(self) -> Dict[Symbol, complex]:
        return {
            # Rabi frequencies
            self.C1.symbol: self.C1.value,
            self.C2.symbol: self.C2.value,

            # D11, ..., D12
            self.d_11.symbol: self.d_11.value,
            self.d_22.symbol: self.d_22.value,
            self.d_33.symbol: self.d_33.value,
            self.d_12.symbol: self.d_12.value,

            # gammas
            self.gamma.symbol: self.gamma.value,
            self.gamma_31.symbol: self.gamma_31.value,
            self.gamma_32.symbol: self.gamma_32.value,

            # deltas
            self.delta_1.symbol: self.delta_1.value,
            self.delta_2.symbol: self.delta_2.value,

            # other coefficients
            self.g_parallel.symbol: self.g_parallel.value,
            self.g_perpendicular.symbol: self.g_perpendicular.value,

            # other consts
            self.q.symbol: self.q.value,
            self.a.symbol: self.a.value,
            self.h.symbol: self.h.value,

            # variable r in the scheme
            **self.ri.values(),
            **self.ri_plus_half.values(),
            **self.ri_minus_half.values()
        }

    @lru_cache(maxsize=None)
    def calculate_r(self) -> Sequence[float]:
        return [v.real for v in self.ri.values().values()]

    def __hash__(self) -> int:
        # Each equation system is uniquely identified by N and R
        return hash(self.N) * 256 + hash(self.R) * 128

    @staticmethod
    @lru_cache(maxsize=None)
    def acquire_equation_system(intervals_number: int = 4, radius: float = DEFAULT_RADIUS):
        return EquationSystem(intervals_number=intervals_number, radius=radius)
