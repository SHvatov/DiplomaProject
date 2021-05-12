import sys
import math

ZERO = 0.0
PRECISION = 1e-20

def pprint_matrix(logs_file_path, matrix_dim, print_matr=False):
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
    
    if print_matr:
        # print indices
        for index in range(0, matrix_dim):
            print("{:45d}".format(index + 1), end=" ")
        print()

        # print matrix values
        index = 1
        for row in matrix:
            print("{:3d}".format(index), end=" ")
            for elem in row:
                real, img = elem
                symbol = ''
                if math.isclose(a=abs(real), b=ZERO, rel_tol=PRECISION) and \
                    math.isclose(a=abs(img), b=ZERO, rel_tol=PRECISION):
                    symbol = "-"
                else:
                    symbol = "+"
                print("({:20.3f} {:20.3f} {})".format(real, img, symbol), end=" ")
            index += 1
            print()

    for row in matrix:
        for elem in row:
            real, img = elem
            if math.isclose(a=abs(real), b=ZERO, rel_tol=PRECISION) and \
                math.isclose(a=abs(img), b=ZERO, rel_tol=PRECISION):
                print("-", end=" ")
            else:
                print("+", end=" ")

        print()


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) != 3:
        print(
            """
            Inocorrect number of arguments provided!
            Expected 3 - path to the log file, matrix dimension,  
            flag whether to print the matrix itself.
            """
        )
        print(f"Got: ${args}")
        exit(-1)
        
    pprint_matrix(
        logs_file_path=args[0],
        matrix_dim=int(args[1]), 
        print_matr=(args[2].lower() == 'true')
    ) 