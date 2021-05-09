program main
    use VecMatrUtils

    implicit none

    complex, dimension(0:3) :: testVec1
    complex, dimension(0:1, 0:2) :: testMatr1

    ! Test: print vector
    testVec1(0) = 0
    testVec1(1) = 1
    testVec1(2) = 2
    testVec1(3) = 3
    print *, "Priniting testVec1"
    call printComplexVector(testVec1, 4)

    ! Test: print matrix
    testMatr1(0, 0) = 0
    testMatr1(0, 1) = 1
    testMatr1(0, 2) = 2

    testMatr1(1, 0) = 3
    testMatr1(1, 1) = 4
    testMatr1(1, 2) = 5
    print *, "Priniting testVec1"
    call printComplexMatrix(testMatr1, 2, 3)
end program main