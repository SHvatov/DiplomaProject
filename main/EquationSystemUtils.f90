module EquationSystemUtils
    use Constants
    use CommonFunctions
    
    implicit none
    
    ! Internal implementation details - indicies of the Ro elements in the matrix 
    integer, parameter, private :: RO_11 = 1, RO_22 = 2, RO_33 = 3, RO_12 = 4
    integer, parameter, private :: RO_11_CONJG = 5, RO_22_CONJG = 6, RO_33_CONJG = 7, RO_12_CONJG = 8
contains
    function Ro23Point(roMeshMatr, i) result(retval)
        implicit none
        complex, dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex :: retval
        
        retval = IMG_UNIT * Omega1(i) * Ro12PointConjg(roMeshMatr, i) / DeltaStroke21 &
            - IMG_UNIT * Omega1(i) * (Ro33Point(roMeshMatr, i) - Ro22Point(roMeshMatr, i)) / DeltaStroke21
    end function Ro23Point

    function Ro13Point(roMeshMatr, i) result(retval)
        implicit none
        complex, dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex :: retval
        
        retval = IMG_UNIT * Omega2(i) * Ro12Point(roMeshMatr, i) / DeltaStroke11 &
            - IMG_UNIT * Omega1(i) * (Ro33Point(roMeshMatr, i) - Ro11Point(roMeshMatr, i)) / DeltaStroke11
    end function Ro13Point

    function Ro11Point(roMeshMatr, i) result(retval)
        implicit none
        complex, dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex :: retval
        
        retval = roMeshMatr(RO_11, i)
    end function Ro11Point

    function Ro22Point(roMeshMatr, i) result(retval)
        implicit none
        complex, dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex :: retval
        
        retval = roMeshMatr(RO_22, i)
    end function Ro22Point

    function Ro33Point(roMeshMatr, i) result(retval)
        implicit none
        complex, dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex :: retval
        
        retval = roMeshMatr(RO_33, i)
    end function Ro33Point

    function Ro12Point(roMeshMatr, i) result(retval)
        implicit none
        complex, dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex :: retval
        
        retval = roMeshMatr(RO_12, i)
    end function Ro12Point

    function Ro11PointConjg(roMeshMatr, i) result(retval)
        implicit none
        complex, dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex :: retval
        
        retval = roMeshMatr(RO_11_CONJG, i)
    end function Ro11PointConjg

    function Ro22PointConjg(roMeshMatr, i) result(retval)
        implicit none
        complex, dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex :: retval
        
        retval = roMeshMatr(RO_22_CONJG, i)
    end function Ro22PointConjg

    function Ro33PointConjg(roMeshMatr, i) result(retval)
        implicit none
        complex, dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex :: retval
        
        retval = roMeshMatr(RO_33_CONJG, i)
    end function Ro33PointConjg

    function Ro12PointConjg(roMeshMatr, i) result(retval)
        implicit none
        complex, dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex :: retval
        
        retval = roMeshMatr(RO_12_CONJG, i)
    end function Ro12PointConjg

    function Ro23PointConjg(roMeshMatr, i) result(retval)
        implicit none
        complex, dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex :: retval
        
        retval = -IMG_UNIT * conjg(Omega1(i)) * Ro12Point(roMeshMatr, i) / DeltaStroke22 &
            + IMG_UNIT * conjg(Omega1(i)) * (Ro33PointConjg(roMeshMatr, i) - Ro22PointConjg(roMeshMatr, i)) / DeltaStroke22
    end function Ro23PointConjg

    function Ro13PointConjg(roMeshMatr, i) result(retval)
        implicit none
        complex, dimension(1:SYSTEM_VAR_NUM, 0:N), intent(in) :: roMeshMatr
        integer :: i
        complex :: retval
        
        retval = -IMG_UNIT * conjg(Omega2(i)) * Ro12PointConjg(roMeshMatr, i) / DeltaStroke12 &
            + IMG_UNIT * conjg(Omega1(i)) * (Ro33PointConjg(roMeshMatr, i) - Ro11PointConjg(roMeshMatr, i)) / DeltaStroke12
    end function Ro13PointConjg
end module EquationSystemUtils