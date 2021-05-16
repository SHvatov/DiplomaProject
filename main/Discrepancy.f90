module discrepancy
    use Constants
    use DebugConfig
    use EquationSystemUtils
    use VecMatrUtils

    implicit none
contains
    ! Calculates the dicrepancy vector based on the provided matrix of values of the functions 
    ! in the auxiliary mesh.
    ! Params:
    ! - roMeshMatr - matrix with dimensions (1:SYSTEM_VAR_NUM, 0:N), which contains
    ! the values of the functions Ro11, Ro22, Ro33, Ro12, Ro11*, Ro22*, Ro33*, Ro12*, in the points of the auxiliary mesh.
    ! This values are used to calculate Ro23, Ro13, Ro23*, Ro13*, which then used to calculate the dicrepancy vector.
    ! - psiVector - vector, which will contain the dicrepancy vector.
    subroutine calculateDiscrepancy(roMeshMatr, psiVector)
        implicit none
        complex(16), dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        complex(16), dimension(1:EXTENDED_MESH_DIM), intent(out) :: psiVector
        integer :: i

        ! DECLARATIONS
        ! i = 0
        complex(16), dimension(0:3) :: leftMain
        complex(16), dimension(0:3) :: leftConjg
        
        ! i = [1; Nr - 1]
        complex(16), dimension(0:3, 1:N-1) :: main
        complex(16), dimension(0:3, 1:N-1) :: mainConjg

        ! i = Nr
        complex(16), dimension(0:3) :: rightMain
        complex(16), dimension(0:3) :: rightConjg

        ! EQUATIONS
        ! i = 0
        leftMain(0) = calculateA1(roMeshMatr, 0) * Hr / 4 * riPlusHalf(0) &
            - D11 *  (Ro11Point(roMeshMatr, 1) - Ro11Point(roMeshMatr, 0)) / Hr * riPlusHalf(0)
        leftMain(1) = calculateA2(roMeshMatr, 0) * Hr / 4 * riPlusHalf(0) &
            - D22 * (Ro22Point(roMeshMatr, 1) - Ro22Point(roMeshMatr, 0)) / Hr * riPlusHalf(0)
        leftMain(2) = calculateA3(roMeshMatr, 0) * Hr / 4 * riPlusHalf(0) &
            - D33 * (Ro33Point(roMeshMatr, 1) - Ro33Point(roMeshMatr, 0)) / Hr * riPlusHalf(0)
        leftMain(3) = calculateA4(roMeshMatr, 0) * Hr / 4 * riPlusHalf(0) &
            - D12 * (Ro12Point(roMeshMatr, 1) - Ro12Point(roMeshMatr, 0)) / Hr * riPlusHalf(0)

        if (DEBUG_DISC) then
            print *, "Discrepancy, left border"
            call printComplexVector(leftMain, 4)
        end if

        leftConjg(0) = calculateA1Conjg(roMeshMatr, 0) * Hr / 4 * riPlusHalf(0) &
            - D11 * (Ro11PointConjg(roMeshMatr, 1) - Ro11PointConjg(roMeshMatr, 0)) / Hr * riPlusHalf(0)
        leftConjg(1) = calculateA2Conjg(roMeshMatr, 0) * Hr / 4 * riPlusHalf(0) &
            - D22 * (Ro22PointConjg(roMeshMatr, 1) - Ro22PointConjg(roMeshMatr, 0)) / Hr * riPlusHalf(0)
        leftConjg(2) = calculateA3Conjg(roMeshMatr, 0) * Hr / 4 * riPlusHalf(0) &
            - D33 * (Ro33PointConjg(roMeshMatr, 1) - Ro33PointConjg(roMeshMatr, 0)) / Hr * riPlusHalf(0)
        leftConjg(3) = calculateA4Conjg(roMeshMatr, 0) * Hr / 4 * riPlusHalf(0) &
            - D12 * (Ro12PointConjg(roMeshMatr, 1) - Ro12PointConjg(roMeshMatr, 0)) / Hr * riPlusHalf(0)

        if (DEBUG_DISC) then
            print *, "Discrepancy, left border, conjg"
            call printComplexVector(leftConjg, 4)
        end if

        ! i = [1; Nr - 1]
        do i = 1, (N - 1)
            ! MAIN EQUATION
            main(0, i) = calculateA1(roMeshMatr, i) * Hr * ri(i) &
                - D11 * &
                    ( &
                        riPlusHalf(i) * (Ro11Point(roMeshMatr, i + 1) - Ro11Point(roMeshMatr, i)) / Hr &
                        - riMinusHalf(i) * (Ro11Point(roMeshMatr, i) - Ro11Point(roMeshMatr, i - 1)) / Hr &
                    )
            main(1, i) = calculateA2(roMeshMatr, i) * Hr * ri(i) &
                - D22 * &
                    ( &
                        riPlusHalf(i) * (Ro22Point(roMeshMatr, i + 1) - Ro22Point(roMeshMatr, i)) / Hr &
                        - riMinusHalf(i) * (Ro22Point(roMeshMatr, i) - Ro22Point(roMeshMatr, i - 1)) / Hr &
                    )
            main(2, i) = calculateA3(roMeshMatr, i) * Hr * ri(i) &
                - D33 * &
                    ( &
                        riPlusHalf(i) * (Ro33Point(roMeshMatr, i + 1) - Ro33Point(roMeshMatr, i)) / Hr &
                        - riMinusHalf(i) * (Ro33Point(roMeshMatr, i) - Ro33Point(roMeshMatr, i - 1)) / Hr &
                    )
            main(3, i) = calculateA4(roMeshMatr, i) * Hr * ri(i) &
                - D12 * &
                    ( &
                        riPlusHalf(i) * (Ro12Point(roMeshMatr, i + 1) - Ro12Point(roMeshMatr, i)) / Hr &
                        - riMinusHalf(i) * (Ro12Point(roMeshMatr, i) - Ro12Point(roMeshMatr, i - 1)) / Hr &
                    )

            ! CONJG EQUATION
            mainConjg(0, i) = calculateA1Conjg(roMeshMatr, i) * Hr * ri(i) &
                - D11 * &
                    ( &
                        riPlusHalf(i) * (Ro11PointConjg(roMeshMatr, i + 1) - Ro11PointConjg(roMeshMatr, i)) / Hr &
                        - riMinusHalf(i) * (Ro11PointConjg(roMeshMatr, i) - Ro11PointConjg(roMeshMatr, i - 1)) / Hr &
                    )
            mainConjg(1, i) = calculateA2Conjg(roMeshMatr, i) * Hr * ri(i) &
                - D22 * &
                    ( &
                        riPlusHalf(i) * (Ro22PointConjg(roMeshMatr, i + 1) - Ro22PointConjg(roMeshMatr, i)) / Hr &
                        - riMinusHalf(i) * (Ro22PointConjg(roMeshMatr, i) - Ro22PointConjg(roMeshMatr, i - 1)) / Hr &
                    )
            mainConjg(2, i) = calculateA3Conjg(roMeshMatr, i) * Hr * ri(i) &
                - D33 * &
                    ( &
                        riPlusHalf(i) * (Ro33PointConjg(roMeshMatr, i + 1) - Ro33PointConjg(roMeshMatr, i)) / Hr &
                        - riMinusHalf(i) * (Ro33PointConjg(roMeshMatr, i) - Ro33PointConjg(roMeshMatr, i - 1)) / Hr &
                    )
            mainConjg(3, i) = calculateA4Conjg(roMeshMatr, i) * Hr * ri(i) &
                - D12 * &
                    ( &
                        riPlusHalf(i) * (Ro12PointConjg(roMeshMatr, i + 1) - Ro12PointConjg(roMeshMatr, i)) / Hr &
                        - riMinusHalf(i) * (Ro12PointConjg(roMeshMatr, i) - Ro12PointConjg(roMeshMatr, i - 1)) / Hr &
                    )
        end do

        if (DEBUG_DISC) then
            print *, "Discrepancy, main equation system"
            call printComplexMatrixSlice(main, 0, 4, 1, N - 1)

            print *, "Discrepancy, main equation system, conjg"
            call printComplexMatrixSlice(mainConjg, 0, 4, 1, N - 1)
        end if

        ! i = Nr
        rightMain(0) = 0.5 - Ro11Point(roMeshMatr, N)
        rightMain(1) = 0.5 - Ro22Point(roMeshMatr, N)
        rightMain(2) = -Ro33Point(roMeshMatr, N)
        rightMain(3) = -Ro12Point(roMeshMatr, N)

        if (DEBUG_DISC) then
            print *, "Discrepancy, right border"
            call printComplexVector(rightMain, 4)
        end if

        rightConjg(0) = 0.5 - Ro11PointConjg(roMeshMatr, N)
        rightConjg(1) = 0.5 - Ro22PointConjg(roMeshMatr, N)
        rightConjg(2) = -Ro33PointConjg(roMeshMatr, N)
        rightConjg(3) = -Ro12PointConjg(roMeshMatr, N)

        if (DEBUG_DISC) then
            print *, "Discrepancy, right border, conjg"
            call printComplexVector(rightConjg, 4)
        end if

        ! Fill ouput vector
        psiVector(1) = leftMain(0)
        psiVector(2) = leftMain(1)
        psiVector(3) = leftMain(2)
        psiVector(4) = leftMain(3)

        psiVector(5) = leftConjg(0)
        psiVector(6) = leftConjg(1)
        psiVector(7) = leftConjg(2)
        psiVector(8) = leftConjg(3)

        do i = 1, (N - 1)
            ! 8N - 7
            psiVector(8 * i + 1) = main(0, i)
            ! 8N - 6
            psiVector(8 * i + 2) = main(1, i)
            ! 8N - 5
            psiVector(8 * i + 3) = main(2, i)
            ! 8N - 4
            psiVector(8 * i + 4) = main(3, i)
            ! 8N - 3
            psiVector(8 * i + 5) = mainConjg(0, i)
            ! 8N - 2
            psiVector(8 * i + 6) = mainConjg(1, i)
            ! 8N - 1
            psiVector(8 * i + 7) = mainConjg(2, i)
            ! 8N - 0
            psiVector(8 * i + 8) = mainConjg(3, i)
        end do

        ! 8N + 1
        psiVector(EXTENDED_MESH_DIM - 7) = rightMain(0)
        psiVector(EXTENDED_MESH_DIM - 6) = rightMain(1)
        psiVector(EXTENDED_MESH_DIM - 5) = rightMain(2)
        psiVector(EXTENDED_MESH_DIM - 4) = rightMain(3)

        psiVector(EXTENDED_MESH_DIM - 3) = rightConjg(0)
        psiVector(EXTENDED_MESH_DIM - 2) = rightConjg(1)
        psiVector(EXTENDED_MESH_DIM - 1) = rightConjg(2)
        psiVector(EXTENDED_MESH_DIM) = rightConjg(3)
    end subroutine calculateDiscrepancy

    ! Helper functions, that are used to calculate left side of the equations
    function calculateA1(roMeshMatr, i) result(retval)
        implicit none
        complex(16), dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex(16) :: retval
        
        retval = & 
            ( &
                IMG_UNIT * Omega1(i) * Ro13PointConjg(roMeshMatr, i) &
                - IMG_UNIT * Omega1Conjg(i) * Ro13Point(roMeshMatr, i) &
                - Gamma31 * Ro33Point(roMeshMatr, i) &
                + GParallel * (Ro11Point(roMeshMatr, i) - Ro22Point(roMeshMatr, i)) &
            )
    end function calculateA1

    function calculateA2(roMeshMatr, i) result(retval)
        implicit none
        complex(16), dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex(16) :: retval
        
        retval = & 
            ( &
                IMG_UNIT * Omega2(i) * Ro23PointConjg(roMeshMatr, i) &
                - IMG_UNIT * Omega2Conjg(i) * Ro23Point(roMeshMatr, i) &
                - Gamma32 * Ro33Point(roMeshMatr, i) &
                + GParallel * (Ro22Point(roMeshMatr, i) - Ro11Point(roMeshMatr, i)) &
            )
    end function calculateA2

    function calculateA3(roMeshMatr, i) result(retval)
        implicit none
        complex(16), dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex(16) :: retval
        
        retval = & 
            ( &
                IMG_UNIT * Omega1Conjg(i) * Ro13Point(roMeshMatr, i) &
                - IMG_UNIT * Omega1(i) * Ro13PointConjg(roMeshMatr, i) &
                + IMG_UNIT * Omega2Conjg(i) * Ro23Point(roMeshMatr, i) &
                - IMG_UNIT * Omega2(i) * Ro23PointConjg(roMeshMatr, i) &
                + Gamma * Ro33Point(roMeshMatr, i) &
            )
    end function calculateA3

    function calculateA4(roMeshMatr, i) result(retval)
        implicit none
        complex(16), dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex(16) :: retval
        
        retval = & 
            ( &
                DeltaGamma * Ro12Point(roMeshMatr, i) &
                - IMG_UNIT * Omega2Conjg(i) * Ro13Point(roMeshMatr, i) &
                + IMG_UNIT * Omega1(i) * Ro23PointConjg(roMeshMatr, i) &
            )
    end function calculateA4

    function calculateA1Conjg(roMeshMatr, i) result(retval)
        implicit none
        complex(16), dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex(16) :: retval
        
        retval = &
            ( &
                - IMG_UNIT * Omega1Conjg(i) * Ro13Point(roMeshMatr, i) &
                + IMG_UNIT * Omega1(i) * Ro13PointConjg(roMeshMatr, i) &
                - Gamma31 * Ro33PointConjg(roMeshMatr, i) &
                + GParallel * (Ro11PointConjg(roMeshMatr, i) - Ro22PointConjg(roMeshMatr, i)) &
            )
    end function calculateA1Conjg

    function calculateA2Conjg(roMeshMatr, i) result(retval)
        implicit none
        complex(16), dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex(16) :: retval
        
        retval = & 
            ( &
                - IMG_UNIT * Omega2Conjg(i) * Ro23Point(roMeshMatr, i) &
                + IMG_UNIT * Omega2(i) * Ro23PointConjg(roMeshMatr, i) &
                - Gamma32 * Ro33PointConjg(roMeshMatr, i) &
                + GParallel * (Ro22PointConjg(roMeshMatr, i) - Ro11PointConjg(roMeshMatr, i)) &
            )
    end function calculateA2Conjg

    function calculateA3Conjg(roMeshMatr, i) result(retval)
        implicit none
        complex(16), dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex(16) :: retval
        
        retval = & 
            ( &
                - IMG_UNIT * Omega1(i) * Ro13PointConjg(roMeshMatr, i) &
                + IMG_UNIT * Omega1Conjg(i) * Ro13Point(roMeshMatr, i) &
                - IMG_UNIT * Omega2(i) * Ro23PointConjg(roMeshMatr, i) &
                + IMG_UNIT * Omega2Conjg(i) * Ro23Point(roMeshMatr, i) &
                + Gamma * Ro33PointConjg(roMeshMatr, i) &
            )
    end function calculateA3Conjg

    function calculateA4Conjg(roMeshMatr, i) result(retval)
        implicit none
        complex(16), dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex(16) :: retval
        
        retval = &
            ( &
                DeltaGamma * Ro12PointConjg(roMeshMatr, i) &
                + IMG_UNIT * Omega2(i) * Ro13PointConjg(roMeshMatr, i) &
                - IMG_UNIT * Omega1Conjg(i) * Ro23Point(roMeshMatr, i) &
            )
    end function calculateA4Conjg
end module discrepancy