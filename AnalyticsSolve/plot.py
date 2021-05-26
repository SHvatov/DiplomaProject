from typing import List, Dict, Tuple, Sequence, Callable

from matplotlib import pyplot as plt

from equation import EquationSystem


def plot_single_solution(func_name: str,
                         intervals: int,
                         indexed_points: Sequence[Tuple[int, complex]],
                         r_provider: Callable[[int], float]) -> None:
    r, real, imag = [], [], []
    for index, point in indexed_points:
        r.append(r_provider(index))
        real.append(point.real)
        imag.append(point.imag)

    figure, axis = plt.subplots(2)

    axis[0].plot(r, real)
    axis[0].set_title(f"{func_name} - Real part")

    axis[1].plot(r, imag)
    axis[1].set_title(f"{func_name} - Imaginary part")

    plt.subplots_adjust(hspace=0.5)
    plt.savefig(f"../data/{func_name}-{intervals}.png")


def plot_solution(solutions: Sequence[complex]) -> None:
    intervals = len(solutions) // 8 - 1
    eq_system = EquationSystem.acquire_equation_system(intervals_number=intervals)

    solutions_by_var: Dict[str, List[Tuple[int, complex]]] = dict()
    for variable, solution in zip(eq_system.ordered_variables(), solutions):
        var_name = str(variable.name)
        func_name = var_name[0:var_name.index("[")]
        index = int(var_name[len(func_name) + 1:var_name.index("]")])
        if func_name in solutions_by_var.keys():
            solutions_by_var[func_name].append((index, solution))
        else:
            solutions_by_var[func_name] = [(index, solution)]

    plot_single_solution("r11", intervals, solutions_by_var["r11"], lambda i: eq_system.calculate_r()[i])
    plot_single_solution("r22", intervals, solutions_by_var["r22"], lambda i: eq_system.calculate_r()[i])
    plot_single_solution("r33", intervals, solutions_by_var["r33"], lambda i: eq_system.calculate_r()[i])
    plot_single_solution("r12", intervals, solutions_by_var["r12"], lambda i: eq_system.calculate_r()[i])
