module Functions 
    use Constants
    
    implicit none
contains
    ! Calculates the length of the insterval based on the radius (R) 
    ! and number of the intervals (N).
    function hr(R, N) result(retval)
        implicit none
        real :: R, retval
        integer :: N
        
        retval = R / N
    end function hr

    ! Calculates the r[i] based on the provided index (i), radius (R) 
    ! and number of the intervals (N).
    function ri(i, R, N) result (retval)
        implicit none
        integer :: i
        real :: R, retval
        integer :: N
        
        retval = hr(R, N) * i
    end function ri

    ! Omega functions.
    ! Returns the complex result of the Omega function in the point (ri).
    function Omega1(ri) result(retval)
        implicit none
        real :: ri
        complex :: retval
        
        retval = C1 * e ** (ri / a)
    end function Omega1

    function Omega2(ri) result(retval)
        implicit none
        real :: ri
        complex :: retval
        
        retval = C2 * e ** (ri / a)
    end function Omega2

    ! RO functions.
    ! Returns the complex result of the RO function in the point (ri).
    function Ro11(ri) result (retval)
        implicit none
        real :: ri
        complex :: retval
        
        retval = (0, 0)
    end function Ro11

    function Ro11Conjg(ri) result (retval)
        implicit none
        real :: ri
        complex :: retval
        
        retval = conjg(Ro11(ri))
    end function Ro11Conjg

    function Ro22(ri) result (retval)
        implicit none
        real :: ri
        complex :: retval
        
        retval = (0, 0)
    end function Ro22

    function Ro22Conjg(ri) result (retval)
        implicit none
        real :: ri
        complex :: retval
        
        retval = conjg(Ro22(ri))
    end function Ro22Conjg

    function Ro33(ri) result (retval)
        implicit none
        real :: ri
        complex :: retval
        
        retval = (0, 0)
    end function Ro33

    function Ro33Conjg(ri) result (retval)
        implicit none
        real :: ri
        complex :: retval
        
        retval = conjg(Ro33(ri))
    end function Ro33Conjg

    function Ro12(ri) result (retval)
        implicit none
        real :: ri
        complex :: retval
        
        retval = (0, 0)
    end function Ro12

    function Ro12Conjg(ri) result (retval)
        implicit none
        real :: ri
        complex :: retval
        
        retval = conjg(Ro12(ri))
    end function Ro12Conjg

    ! Ro13, Ro31, Ro32, Ro23 functions.
    ! Calculations are based on the equations (4) - (5) of the main Functions.
    ! Ro13 = Ro31*, Ro23 = Ro32*
    function Ro13(ri) result (retval)
        implicit none
        real :: ri
        complex :: retval
        
        retval = IMG_UNIT * Omega2(ri) * Ro12(ri) / DeltaStroke11 &
            - IMG_UNIT * Omega1(ri) * (Ro33(ri) - Ro11(ri)) / DeltaStroke11
    end function Ro13

    function Ro13Conjg(ri) result (retval)
        implicit none
        real :: ri
        complex :: retval
        
        retval = -IMG_UNIT * conjg(Omega2(ri)) * conjg(Ro12(ri)) / DeltaStroke12 &
            + IMG_UNIT * conjg(Omega1(ri)) * (conjg(Ro33(ri)) - conjg(Ro11(ri))) / DeltaStroke12
    end function Ro13Conjg

    function Ro23(ri) result (retval)
        implicit none
        real :: ri
        complex :: retval
        
        retval = IMG_UNIT * Omega1(ri) * conjg(Ro12(ri)) / DeltaStroke21 &
            - IMG_UNIT * Omega1(ri) * (Ro33(ri) - Ro11(ri)) / DeltaStroke21
    end function Ro23

    function Ro23Conjg(ri) result (retval)
        implicit none
        real :: ri
        complex :: retval
        
        retval = -IMG_UNIT * conjg(Omega1(ri)) * Ro12(ri) / DeltaStroke22 &
            + IMG_UNIT * conjg(Omega1(ri)) * (conjg(Ro33(ri)) - conjg(Ro11(ri))) / DeltaStroke22
    end function Ro23Conjg
end module Functions