gfortran ../main/*.f90 -o main.out
./main.out > log.txt
python3 matrix_printer.py ./log.txt 40 > matrix.txt
rm *.mod
rm *.out
