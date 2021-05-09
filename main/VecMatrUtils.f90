module VecMatrUtils
    implicit none
contains
    ! Prints the provided vector of complex type (vec) and specified length (L).
    ! Note: 
    ! - procedure does not check the boundaries
    ! - vector is expected to start with 0 index
    subroutine printComplexVector(vec, L)
        implicit none
        complex, dimension(0:), intent(in) :: vec
        integer, intent(in) :: L

        integer :: i
        do i = 0, L - 1
            print *, "vec[", i, "] = ", vec(i)
        end do
    end subroutine printComplexVector

    ! Does the same as printComplexVector but for matricies [L * K]
    subroutine printComplexMatrix(matr, L, K)
        implicit none
        complex, dimension(0:, 0:), intent(in) :: matr
        integer, intent(in) :: L, K

        integer :: i, j
        rowloop: do i = 0, L - 1
            columnloop: do j = 0, K - 1
            print *, "matr[", i, ", ", j, "] = ", matr(i, j)
            end do columnloop
        end do rowloop
    end subroutine printComplexMatrix

    ! Analogs of the printComplexVector and printComplexMatrix, which
    ! allows to print the slices from the vectors and matricies.
    subroutine printComplexVectorSlice(vec, from, to)
        implicit none
        complex, dimension(:), intent(in) :: vec
        integer, intent(in) :: from, to

        integer :: i
        do i = from, to
            print *, "vec[", i, "] = ", vec(i)
        end do
    end subroutine printComplexVectorSlice

    subroutine printComplexMatrixSlice(matr, fromRow, toRow, fromColumn, toColumn)
        implicit none
        complex, dimension(:,:), intent(in) :: matr
        integer, intent(in) :: fromRow, toRow, fromColumn, toColumn

        integer :: i, j
        rowloop: do i = fromRow, toRow
            columnloop: do j = fromColumn, toColumn
            print *, "matr[", i, ", ", j, "] = ", matr(i, j)
            end do columnloop
        end do rowloop
    end subroutine printComplexMatrixSlice
end module VecMatrUtils