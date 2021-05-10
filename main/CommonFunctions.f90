module CommonFunctions 
    use Constants
    
    implicit none
contains
    ! Calculates the r[i] based on the provided index (i), radius (R) 
    ! and number of the intervals (N).
    function ri(i) result (retval)
        implicit none
        integer :: i
        real :: retval
        
        retval = Hr * i
    end function ri

    ! Calculates the r[i + 1/2] based on the provided index (i), radius (R) 
    ! and number of the intervals (N).
    function riPlusHalf(i) result (retval)
        implicit none
        integer :: i
        real :: retval

        retval = ri(i) + Hr / 2
    end function riPlusHalf

    ! Calculates the r[i - 1/2] based on the provided index (i), radius (R) 
    ! and number of the intervals (N).
    function riMinusHalf(i) result (retval)
        implicit none
        integer :: i
        real :: retval

        retval = ri(i) + Hr / 2
    end function riMinusHalf

    ! Omega CommonFunctions.
    ! Returns the complex result of the Omega function in the point (r[i]).
    function Omega1(i) result(retval)
        implicit none
        integer :: i
        complex :: retval
        
        retval = C1 * e ** (ri(i) / a)
    end function Omega1

    function Omega1Conjg(i) result(retval)
        implicit none
        integer :: i
        complex :: retval
        
        retval = conjg(Omega1(i))
    end function Omega1Conjg

    function Omega2(i) result(retval)
        implicit none
        integer :: i
        complex :: retval
        
        retval = C2 * e ** (ri(i) / a)
    end function Omega2

    function Omega2Conjg(i) result(retval)
        implicit none
        integer :: i
        complex :: retval
        
        retval = conjg(Omega2(i))
    end function Omega2Conjg
end module CommonFunctions