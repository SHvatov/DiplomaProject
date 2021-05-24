module CoefficientMatrix
    use Discrepancy
    use DebugConfig

    implicit none
contains
    ! Calculates the coefficient matrix A based on the initial approximation and delta 
    ! using the follwoing formula:
    !   A(*, j) = (E(ro + delta * e(j)) - E(ro)) / delta
    ! Params:
    ! - initialRoApproxMesh - matrix with dimensions (1:SYSTEM_VAR_NUM, 0:N), which contains
    ! the initial approximation of the functions Ro11, Ro22, Ro33, Ro12, Ro11*, Ro22*, Ro33*, Ro12*.
    ! - delta - DELTA value
    ! - coeffMatrix - output matrix with dimensions (1:EXTENDED_MESH_DIM, 1:EXTENDED_MESH_DIM),
    ! which will contain the coefficients of the matrix A.
    subroutine prepareCoefficientMatrix(initialRoApproxMesh, delta, coeffMatrix)
        implicit none
        ! Arguments
        complex(8), intent(in) :: delta
        complex(8), dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: initialRoApproxMesh
        complex(8), dimension(1:EXTENDED_MESH_DIM, 1:EXTENDED_MESH_DIM), intent(out) :: coeffMatrix

        ! Local variables
        integer :: j
        complex(8), dimension(1:EXTENDED_MESH_DIM) :: psiVector
        complex(8), dimension(1:EXTENDED_MESH_DIM) :: initialApproxPsiVector
        complex(8), dimension(1:SYSTEM_VAR_NUM, 0:N) :: deltaRoApproxMesh

        ! Calculate the discrepancy in the initial approximation
        call calculateDiscrepancy(initialRoApproxMesh, initialApproxPsiVector)
        if (DEBUG_MATR) then
            print *, "Discrepancy, initial approximation"
            call printComplexVectorSlice(initialApproxPsiVector, 1, EXTENDED_MESH_DIM)
        end if

        if (DEBUG_MATR) then
            print *, "Initial approximation mesh"
            call printComplexMatrixSlice(initialRoApproxMesh, 1, SYSTEM_VAR_NUM, 0, N)
        end if

        do j = 1, EXTENDED_MESH_DIM
            ! Prepare delta approximation mesh
            deltaRoApproxMesh = (0, 0)
            if (DEBUG_MATR) then
                print *, "Delta Approximation mesh before, j = ", j
                call printComplexMatrixSlice(deltaRoApproxMesh, 1, SYSTEM_VAR_NUM, 0, N)
            end if

            call prepareDeltaApproximation(initialRoApproxMesh, delta, j, deltaRoApproxMesh)
            if (DEBUG_MATR) then
                print *, "Delta Approximation mesh, j = ", j
                call printComplexMatrixSlice(deltaRoApproxMesh, 1, SYSTEM_VAR_NUM, 0, N)
            end if

            ! Calculate discrepancy
            psiVector = 0
            call calculateDiscrepancy(deltaRoApproxMesh, psiVector)
            if (DEBUG_MATR) then
                print *, "Discrepancy, delta approximation, j = ", j
                call printComplexVectorSlice(psiVector, 1, EXTENDED_MESH_DIM)
            end if

            ! Set the coeff matrix(*, i) using formula
            coeffMatrix(:, j) = (psiVector - initialApproxPsiVector) / delta
            if (DEBUG_MATR) then
                print *, "A(*, J), j = ", j
                call printComplexVectorSlice(coeffMatrix(:, j), 1, EXTENDED_MESH_DIM)
            end if
        end do

        if (OUTPUT_MATR) then
            print *, "A(*, *)"
            call printComplexMatrixSlice(coeffMatrix, 1, EXTENDED_MESH_DIM, 1, EXTENDED_MESH_DIM)
        end if
    end subroutine prepareCoefficientMatrix

    subroutine prepareDeltaApproximation(initialRoApproxMesh, delta, i, deltaRoApproxMesh)
        implicit none
        ! Arguments
        complex(8), intent(in) :: delta
        integer, intent(in) :: i
        complex(8), dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: initialRoApproxMesh
        complex(8), dimension(1:SYSTEM_VAR_NUM, 0:N), intent(out) :: deltaRoApproxMesh

        ! Local variables
        integer :: functionIterV, pointIterV, functionIndex, pointIndex

        ! Calculate the position of the point to add delta to
        pointIndex = (i - 1) / SYSTEM_VAR_NUM
        functionIndex = i - pointIndex * SYSTEM_VAR_NUM

        if (DEBUG_MATR) then
            print *, "Index [", i, "] -> [func=", functionIndex, ",point=", pointIndex, "]"
        end if

        ! Copy the matricies
        functionIter : do functionIterV = 1, SYSTEM_VAR_NUM
            pointIter : do pointIterV = 0, N
                deltaRoApproxMesh(functionIterV, pointIterV) = initialRoApproxMesh(functionIterV, pointIterV)
            end do pointIter
        end do functionIter

        ! Add delta
        deltaRoApproxMesh(functionIndex, pointIndex) = deltaRoApproxMesh(functionIndex, pointIndex) + delta
    end subroutine prepareDeltaApproximation
end module CoefficientMatrix