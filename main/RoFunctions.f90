module RoFunctions
    use CommonFunctions

    implicit none
contains
    ! RO Functions.
    ! Returns the complex(8) result of the RO function in the point (r[i]).
    function Ro11(i) result (retval)
        implicit none
        integer :: i
        complex(8) :: retval
        
        retval = (0, 0)
    end function Ro11

    function Ro11Conjg(i) result (retval)
        implicit none
        integer :: i
        complex(8) :: retval
        
        retval = conjg(Ro11(i))
    end function Ro11Conjg

    function Ro22(i) result (retval)
        implicit none
        integer :: i
        complex(8) :: retval
        
        retval = (0, 0)
    end function Ro22

    function Ro22Conjg(i) result (retval)
        implicit none
        integer :: i
        complex(8) :: retval
        
        retval = conjg(Ro22(i))
    end function Ro22Conjg

    function Ro33(i) result (retval)
        implicit none
        integer :: i
        complex(8) :: retval
        
        retval = (0, 0)
    end function Ro33

    function Ro33Conjg(i) result (retval)
        implicit none
        integer :: i
        complex(8) :: retval
        
        retval = conjg(Ro33(i))
    end function Ro33Conjg

    function Ro12(i) result (retval)
        implicit none
        integer :: i
        complex(8) :: retval
        
        retval = (0, 0)
    end function Ro12

    function Ro12Conjg(i) result (retval)
        implicit none
        integer :: i
        complex(8) :: retval
        
        retval = conjg(Ro12(i))
    end function Ro12Conjg

    ! Ro13, Ro31, Ro32, Ro23 CommonFunctions.
    ! Calculations are based on the equations (4) - (5) of the main CommonFunctions.
    ! Ro13 = Ro31*, Ro23 = Ro32*
    function Ro13(i) result (retval)
        implicit none
        integer :: i
        complex(8) :: retval
        
        retval = IMG_UNIT * Omega2(i) * Ro12(i) / DeltaStroke11 &
            - IMG_UNIT * Omega1(i) * (Ro33(i) - Ro11(i)) / DeltaStroke11
    end function Ro13

    function Ro13Conjg(i) result (retval)
        implicit none
        integer :: i
        complex(8) :: retval
        
        retval = -IMG_UNIT * conjg(Omega2(i)) * conjg(Ro12(i)) / DeltaStroke12 &
            + IMG_UNIT * conjg(Omega1(i)) * (conjg(Ro33(i)) - conjg(Ro11(i))) / DeltaStroke12
    end function Ro13Conjg

    function Ro23(i) result (retval)
        implicit none
        integer :: i
        complex(8) :: retval
        
        retval = IMG_UNIT * Omega1(i) * conjg(Ro12(i)) / DeltaStroke21 &
            - IMG_UNIT * Omega1(i) * (Ro33(i) - Ro11(i)) / DeltaStroke21
    end function Ro23

    function Ro23Conjg(i) result (retval)
        implicit none
        integer :: i
        complex(8) :: retval
        
        retval = -IMG_UNIT * conjg(Omega1(i)) * Ro12(i) / DeltaStroke22 &
            + IMG_UNIT * conjg(Omega1(i)) * (conjg(Ro33(i)) - conjg(Ro11(i))) / DeltaStroke22
    end function Ro23Conjg
end module RoFunctions