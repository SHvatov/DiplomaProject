module discrepancy
    use Constants
    use Functions
    use VecMatrUtils

    implicit none
contains
    ! TODO: replace function calls with the usage of the arrays with RO values
    ! in the points of the net. Add two functions, one accepts vector 8 * (N + 1)
    ! and main one, which accepts 8 vectors of (N + 1) len.

    ! Calculates the discrepancy vector (psiVector) based 
    ! on the functions and constants specified in the corresponding modules.
    subroutine calculateDiscrepancy(psiVector)
        implicit none
        complex, dimension(1:8 * (N + 1)), intent(out) :: psiVector
        integer :: i

        ! DECLARATIONS
        ! i = 0
        complex, dimension(0:4) :: leftMain
        complex, dimension(0:4) :: leftConjg
        
        ! i = [1; Nr - 1]
        complex, dimension(0:4, 1:N-1) :: main
        complex, dimension(0:4, 1:N-1) :: mainConjg

        ! i = Nr
        complex, dimension(0:4) :: rightMain
        complex, dimension(0:4) :: rightConjg

        ! EQUATIONS
        ! i = 0
        leftMain(0) = calculateA1(0) * Hr / 4 &
            - D11 * (riPlusHalf(0) * (Ro11(1) - Ro11(0)) / Hr)
        leftMain(1) = calculateA2(0) * Hr / 4 &
            - D22 * (riPlusHalf(0) * (Ro22(1) - Ro22(0)) / Hr)
        leftMain(2) = calculateA3(0) * Hr / 4 &
            - D33 * (riPlusHalf(0) * (Ro33(1) - Ro33(0)) / Hr)
        leftMain(3) = calculateA4(0) * Hr / 4 &
            - D12 * (riPlusHalf(0) * (Ro12(1) - Ro12(0)) / Hr)

        if (DEBUG == 1) then
            print *, "Discrepancy, left border"
            call printComplexVector(leftMain, 4)
        end if

        leftConjg(0) = calculateA1Conjg(0) * Hr / 4 &
            - D11 * (riPlusHalf(0) * (Ro11Conjg(1) - Ro11Conjg(0)) / Hr)
        leftConjg(1) = calculateA2Conjg(0) * Hr / 4 &
            - D22 * (riPlusHalf(0) * (Ro22Conjg(1) - Ro22Conjg(0)) / Hr)
        leftConjg(2) = calculateA3Conjg(0) * Hr / 4 &
            - D33 * (riPlusHalf(0) * (Ro33Conjg(1) - Ro33Conjg(0)) / Hr)
        leftConjg(3) = calculateA4Conjg(0) * Hr / 4 &
            - D12 * (riPlusHalf(0) * (Ro12Conjg(1) - Ro12Conjg(0)) / Hr)

        if (DEBUG == 1) then
            print *, "Discrepancy, left border, conjg"
            call printComplexVector(leftConjg, 4)
        end if

        ! i = [1; Nr - 1]
        do i = 1, (N - 1)
            ! MAIN EQUATION
            main(0, i) = calculateA1(i) * Hr * ri(i) &
                - D11 * &
                    ( &
                        riPlusHalf(i) * (Ro11(i + 1) - Ro11(i)) / Hr &
                        - riMinusHalf(i) * (Ro11(i) - Ro11(i - 1)) / Hr &
                    )
            main(1, i) = calculateA2(i) * Hr * ri(i) &
                - D22 * &
                    ( &
                        riPlusHalf(i) * (Ro22(i + 1) - Ro22(i)) / Hr &
                        - riMinusHalf(i) * (Ro22(i) - Ro22(i - 1)) / Hr &
                    )
            main(2, i) = calculateA3(i) * Hr * ri(i) &
                - D33 * &
                    ( &
                        riPlusHalf(i) * (Ro33(i + 1) - Ro33(i)) / Hr &
                        - riMinusHalf(i) * (Ro33(i) - Ro33(i - 1)) / Hr &
                    )
            main(3, i) = calculateA4(i) * Hr * ri(i) &
                - D12 * &
                    ( &
                        riPlusHalf(i) * (Ro12(i + 1) - Ro12(i)) / Hr &
                        - riMinusHalf(i) * (Ro12(i) - Ro12(i - 1)) / Hr &
                    )

            ! CONJG EQUATION
            mainConjg(0, i) = calculateA1Conjg(i) * Hr * ri(i) &
                - D11 * &
                    ( &
                        riPlusHalf(i) * (Ro11Conjg(i + 1) - Ro11Conjg(i)) / Hr &
                        - riMinusHalf(i) * (Ro11Conjg(i) - Ro11Conjg(i - 1)) / Hr &
                    )
            mainConjg(1, i) = calculateA2Conjg(i) * Hr * ri(i) &
                - D22 * &
                    ( &
                        riPlusHalf(i) * (Ro22Conjg(i + 1) - Ro22Conjg(i)) / Hr &
                        - riMinusHalf(i) * (Ro22Conjg(i) - Ro22Conjg(i - 1)) / Hr &
                    )
            mainConjg(2, i) = calculateA3Conjg(i) * Hr * ri(i) &
                - D33 * &
                    ( &
                        riPlusHalf(i) * (Ro33Conjg(i + 1) - Ro33Conjg(i)) / Hr &
                        - riMinusHalf(i) * (Ro33Conjg(i) - Ro33Conjg(i - 1)) / Hr &
                    )
            mainConjg(3, i) = calculateA4Conjg(i) * Hr * ri(i) &
                - D12 * &
                    ( &
                        riPlusHalf(i) * (Ro12Conjg(i + 1) - Ro12Conjg(i)) / Hr &
                        - riMinusHalf(i) * (Ro12Conjg(i) - Ro12Conjg(i - 1)) / Hr &
                    )
        end do

        if (DEBUG == 1) then
            print *, "Discrepancy, main equation system"
            call printComplexMatrixSlice(main, 0, 4, 1, N - 1)

            print *, "Discrepancy, main equation system, conjg"
            call printComplexMatrixSlice(mainConjg, 0, 4, 1, N - 1)
        end if

        ! i = Nr
        rightMain(0) = 0.5 - Ro11(N)
        rightMain(1) = 0.5 - Ro22(N)
        rightMain(2) = -Ro33(N)
        rightMain(3) = -Ro12(N)

        if (DEBUG == 1) then
            print *, "Discrepancy, right border"
            call printComplexVector(rightMain, 4)
        end if

        rightConjg(0) = 0.5 - Ro11Conjg(N)
        rightConjg(1) = 0.5 - Ro22Conjg(N)
        rightConjg(2) = -Ro33Conjg(N)
        rightConjg(3) = -Ro12Conjg(N)

        if (DEBUG == 1) then
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
            psiVector(8 * i + 5) = mainConjg(1, i)
            ! 8N - 2
            psiVector(8 * i + 6) = mainConjg(2, i)
            ! 8N - 1
            psiVector(8 * i + 7) = mainConjg(3, i)
            ! 8N - 0
            psiVector(8 * i + 8) = mainConjg(4, i)
        end do

        ! 8N + 1
        psiVector(8 * (N + 1) - 7) = rightMain(0)
        psiVector(8 * (N + 1) - 6) = rightMain(1)
        psiVector(8 * (N + 1) - 5) = rightMain(2)
        psiVector(8 * (N + 1) - 4) = rightMain(3)

        psiVector(8 * (N + 1) - 3) = rightConjg(0)
        psiVector(8 * (N + 1) - 2) = rightConjg(1)
        psiVector(8 * (N + 1) - 1) = rightConjg(2)
        psiVector(8 * (N + 1)) = rightConjg(3)
    end subroutine calculateDiscrepancy

    subroutine calculateCoeffMatrix()
        implicit none
    
    end subroutine calculateCoeffMatrix

    ! Helper functions, that are used to calculate left side of the equations
    function calculateA1(i) result(retval)
        implicit none
        integer :: i
        complex :: retval
        
        retval = & 
            ( &
                IMG_UNIT * Omega1(i) * Ro13Conjg(i) &
                - IMG_UNIT * Omega1Conjg(i) * Ro13(i) &
                - Gamma31 * Ro33(i) + GParallel * (Ro11(i) - Ro22(i)) &
            )
    end function calculateA1

    function calculateA2(i) result(retval)
        implicit none
        integer :: i
        complex :: retval
        
        retval = & 
            ( &
                IMG_UNIT * Omega2(i) * Ro23Conjg(i) &
                - IMG_UNIT * Omega2Conjg(i) * Ro23(i) &
                - Gamma32 * Ro33(i) + GParallel * (Ro22(i) - Ro11(i)) &
            )
    end function calculateA2

    function calculateA3(i) result(retval)
        implicit none
        integer :: i
        complex :: retval
        
        retval = & 
            ( &
                IMG_UNIT * Omega1Conjg(i) * Ro13(i) &
                - IMG_UNIT * Omega1(i) * Ro13Conjg(i) &
                + IMG_UNIT * Omega2(i) * Ro23(i) &
                - IMG_UNIT * Omega2(i) * Ro23Conjg(i) &
                + Gamma * Ro33(i) &
            )
    end function calculateA3

    function calculateA4(i) result(retval)
        implicit none
        integer :: i
        complex :: retval
        
        retval = & 
            ( &
                DeltaGamma * Ro12(i) &
                - IMG_UNIT * Omega2Conjg(i) * Ro13(i) &
                + IMG_UNIT * Omega1(i) * Ro23Conjg(i) &
            )
    end function calculateA4

    function calculateA1Conjg(i) result(retval)
        implicit none
        integer :: i
        complex :: retval
        
        retval = &
            ( &
                - IMG_UNIT * Omega1Conjg(i) * Ro13(i) &
                + IMG_UNIT * Omega1(i) * Ro13Conjg(i) &
                - Gamma31 * Ro33Conjg(i) + GParallel * (Ro11Conjg(i) - Ro22Conjg(i)) &
            )
    end function calculateA1Conjg

    function calculateA2Conjg(i) result(retval)
        implicit none
        integer :: i
        complex :: retval
        
        retval = & 
            ( &
                - IMG_UNIT * Omega2Conjg(i) * Ro23(i) &
                + IMG_UNIT * Omega2(i) * Ro23Conjg(i) &
                - Gamma32 * Ro33Conjg(i) + GParallel * (Ro22Conjg(i) - Ro11Conjg(i)) &
            )
    end function calculateA2Conjg

    function calculateA3Conjg(i) result(retval)
        implicit none
        integer :: i
        complex :: retval
        
        retval = & 
            ( &
                - IMG_UNIT * Omega1(i) * Ro13Conjg(i) &
                + IMG_UNIT * Omega1Conjg(i) * Ro13(i) &
                - IMG_UNIT * Omega2Conjg(i) * Ro23Conjg(i) &
                + IMG_UNIT * Omega2Conjg(i) * Ro23(i) &
                + Gamma * Ro33Conjg(i) &
            )
    end function calculateA3Conjg

    function calculateA4Conjg(i) result(retval)
        implicit none
        integer :: i
        complex :: retval
        
        retval = &
            ( &
                DeltaGamma * Ro12Conjg(i) &
                + IMG_UNIT * Omega2(i) * Ro13Conjg(i) &
                - IMG_UNIT * Omega1Conjg(i) * Ro23(i) &
            )
    end function calculateA4Conjg
end module discrepancy