from warnings import warn

from sympy import shape
from sympy.matrices import Matrix

DEFAULT_JACOBIAN_PATH = "../data/jacobian.txt"
DEFAULT_RHS_PATH = "../data/rhs.txt"


def save_jacobian(jacobian: Matrix, filename: str = DEFAULT_JACOBIAN_PATH) -> None:
    warn("There is no need to save the Jacobian anymore", DeprecationWarning, stacklevel=2)
    n, _ = shape(jacobian)  # matrix is square
    with open(filename, "w") as f:
        for i in range(0, n):
            column = jacobian[0:n, i]
            for elem in column:
                real, img = elem.as_real_imag()
                f.write(f"{float(real):20.9E} {float(img):20.9E}\n")


def save_rhs(rhs: Matrix, filename: str = DEFAULT_RHS_PATH) -> None:
    warn("There is no need to save the right side vector anymore", DeprecationWarning, stacklevel=2)
    with open(filename, "w") as f:
        for elem in rhs:
            real, img = elem.as_real_imag()
            f.write(f"{float(real):20.9E} {float(img):20.9E}\n")
