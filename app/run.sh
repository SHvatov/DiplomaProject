gfortran -C \
    ../main/DebugConfig.f90 \
    ../main/Constants.f90 \
    ../main/CommonFunctions.f90 \
    ../main/RoFunctions.f90 \
    ../main/VecMatrUtils.f90 \
    ../main/EquationSystemUtils.f90 \
    ../main/Discrepancy.f90 \
    ../main/CoefficientMatrix.f90 \
    ../main/System.f90 \
    ../main/BLAS-3.8.0/*.f \
    ../main/linpack/zgeco.F \
    ../main/linpack/zgefa.F \
    ../main/linpack/zgesl.F \
    ../main/Main.f90 -o main.out
./main.out > log.txt
python3 matrix_printer.py ./log.txt 40 True > matrix.txt
rm *.mod
rm *.out
