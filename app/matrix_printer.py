import sys
import math

ZERO = 0.0
PRECISION = 1e-10

def pprint_matrix(logs_file_path, matrix_dim):
    matrix = [[(0, 0) for _ in range(matrix_dim)] for _ in range(matrix_dim)]
    with open(logs_file_path, 'r') as log_file:
        lines = log_file.readlines()
        for line in lines:
            strip_line = line.replace(' ', '')
            if not strip_line.startswith('matr') or len(strip_line) == 0:
                continue
            
            index_part, value_part = strip_line.split('=')
            
            i, j = index_part.replace('matr', '').replace('[', '').replace(']', '').split(',')
            i, j = int(i), int(j)

            value = value_part.replace('(', '').replace(')', '').split(',')
            value = float(value[0]), float(value[1])

            matrix[i - 1][j - 1] = value
    
    for row in matrix:
        for elem in row:
            real, img = elem
            if math.isclose(a=real, b=ZERO, rel_tol=PRECISION) and \
                math.isclose(a=img, b=ZERO, rel_tol=PRECISION):
                print("-", end=" ")
            else:
                print("+", end=" ")

        print()


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) != 2:
        print("Inocorrect number of arguments provided! Expected 2 - path to the log file and matrix dimension.")
        print(f"Got: ${args}")
        exit(-1)
    pprint_matrix(logs_file_path=args[0], matrix_dim=int(args[1])) 