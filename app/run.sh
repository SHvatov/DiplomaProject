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
  ../main/ExternalDataAdapter.f90 \
  ../main/BLAS-3.8.0/*.f \
  ../main/linpack/zgeco.F \
  ../main/linpack/zgefa.F \
  ../main/linpack/zgesl.F \
  ../main/Main.f90 -o main.out
./main.out >../logs/fortran_logs.txt
python3 matrix_printer.py ../logs/fortran_logs.txt 48 True >../data/matrix_view.txt
rm *.mod
rm *.out
