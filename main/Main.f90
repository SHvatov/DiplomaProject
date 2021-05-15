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

    ! Util variables for z*
    complex, dimension(1:EXTENDED_MESH_DIM) :: z = 0
    integer, dimension(1:EXTENDED_MESH_DIM) :: ipvt = 0
    real :: rcond

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

    ! First call zgeco / zgefa
    ! call zgeco(matrixA, EXTENDED_MESH_DIM, EXTENDED_MESH_DIM, ipvt, rcond, z)
    ! if (DEBUG_MAIN) then
    !     print *, "R Cond"
    !     print *, rcond
    ! end if

    ! Then call zgesl
    ! call zgesl(matrixA, EXTENDED_MESH_DIM, EXTENDED_MESH_DIM, ipvt, vectorB, 0)
    ! if (DEBUG_MAIN) then
    !     print *, "Solution(*)"
    !     call printComplexVectorSlice(vectorB, 1, EXTENDED_MESH_DIM)

    !     print *, "Solution_formatted(*)"
    !     call printComplexVectorSliceFmt(vectorB, 1, EXTENDED_MESH_DIM)
    ! end if
end program main