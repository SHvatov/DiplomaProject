module ExternalDataAdapter
    use Constants

    implicit none
contains
    subroutine read_jacobian(A)
        implicit none

        complex(8), dimension(1:EXTENDED_MESH_DIM, 1:EXTENDED_MESH_DIM), intent(out) :: A

        open(unit = 24, file = "../data/jacobian.txt", status = "old", form = "formatted", action = "read")
        read(24, "(E20.9, 1X, E20.9)") A
        close(unit = 24)
    end subroutine read_jacobian

    subroutine read_rhs(b)
        implicit none

        complex(8), dimension(1:EXTENDED_MESH_DIM), intent(out) :: b

        open(unit = 24, file = "../data/rhs.txt", status = "old", form = "formatted", action = "read")
        read(24, "(E20.9, 1X, E20.9)") b
        close(unit = 24)
    end subroutine read_rhs

end module ExternalDataAdapter