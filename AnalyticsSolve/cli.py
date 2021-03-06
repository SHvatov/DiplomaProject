r"""
This module contains the definition of the CLI argument parser
used in the solving application.
"""
import argparse
from typing import Any


def str2bool(v: Any) -> bool:
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def prepare_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="System solver v0.1, author shvatov")

    parser.add_argument("--method", "-m",
                        action="store",
                        default=1,
                        type=int,
                        choices=[1, 2],
                        dest="method",
                        help="method, which will be used to solve the system; either 1 (SYMPY) or 2 (NUMPY);")

    parser.add_argument("--verbose", "-vb",
                        action="store",
                        default=False,
                        type=str2bool,
                        nargs="?",
                        const=True,
                        dest="verbose",
                        help="whether output should be verbose or not;")

    parser.add_argument("--fortran-cmp", "-fc",
                        action="store",
                        default=False,
                        type=str2bool,
                        nargs="?",
                        const=True,
                        dest="fortran",
                        help="whether to compare with the matrix generated by fortran program or not;")

    parser.add_argument("--discrepancy-cmp", "-dc",
                        action="store",
                        default=False,
                        type=str2bool,
                        nargs="?",
                        const=True,
                        dest="discrepancy",
                        help="whether to compare with the matrix generated using discrepancy function or not;")

    parser.add_argument("--rcond", "-rc",
                        action="store",
                        default=False,
                        type=str2bool,
                        nargs="?",
                        const=True,
                        dest="rcond",
                        help="whether to calculate condition number or not;")

    parser.add_argument("--check", "-c",
                        action="store",
                        default=False,
                        type=str2bool,
                        nargs="?",
                        const=True,
                        dest="check",
                        help="whether to check basic conditions (right borders and ro11 + ro22 + ro33 = 1);")

    parser.add_argument("--plot", "-p",
                        action="store",
                        default=False,
                        type=str2bool,
                        nargs="?",
                        const=True,
                        dest="plot",
                        help="whether to plot functions ro11, ro22, ro33, ro44 or not")

    parser.add_argument("--intervals-number", "-n",
                        action="store",
                        default=4,
                        type=int,
                        dest="n",
                        help="number of the interval to be used in the method;")

    return parser
