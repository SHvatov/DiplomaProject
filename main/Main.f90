program main
    ! Required modules
    use VecMatrUtils
    use CoefficientMatrix

    implicit none

    ! Local variables
    complex :: delta = (1, 0)
    complex, dimension(1:SYSTEM_VAR_NUM, 0:N)  :: initialRoApproxMesh = 0
    complex, dimension(1:EXTENDED_MESH_DIM, 1:EXTENDED_MESH_DIM) :: coeffMatrix

    ! Calculate coeff matrix
    coeffMatrix = (0, 0)
    call prepareCoefficientMatrix(initialRoApproxMesh, delta, coeffMatrix)
end program main