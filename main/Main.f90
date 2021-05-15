program main
    ! Required modules
    use DebugConfig
    use VecMatrUtils
    use CoefficientMatrix
    use Discrepancy

    implicit none

    ! Local variables
    complex :: delta = (1, 0)
    complex, dimension(1:SYSTEM_VAR_NUM, 0:N)  :: initialRoApproxMesh
    complex, dimension(1:EXTENDED_MESH_DIM, 1:EXTENDED_MESH_DIM) :: matrixA
    complex, dimension(1:EXTENDED_MESH_DIM) :: vectorB

    ! Calculate coeff matrix
    matrixA = (0, 0)
    initialRoApproxMesh = 0
    call prepareCoefficientMatrix(initialRoApproxMesh, delta, matrixA)
    if (DEBUG_MAIN) then
        print *, "A(*, *)"
        call printComplexMatrixSlice(matrixA, 1, EXTENDED_MESH_DIM, 1, EXTENDED_MESH_DIM)

        print *, "A_formatted(*, *)"
        call printComplexMatrixSliceFmt(matrixA, 1, EXTENDED_MESH_DIM, 1, EXTENDED_MESH_DIM)
    end if

    ! Calcluate left side vector as discrepancy(0)
    vectorB = (0)
    initialRoApproxMesh = 0
    call calculateDiscrepancy(initialRoApproxMesh, vectorB)
    if (DEBUG_MAIN) then
        print *, "B(*)"
        call printComplexVectorSlice(vectorB, 1, EXTENDED_MESH_DIM)

        print *, "B_formatted(*)"
        call printComplexVectorSliceFmt(vectorB, 1, EXTENDED_MESH_DIM)
    end if

    
end program main