module DebugConfig
    implicit none

    ! Global debug option
    logical, parameter :: DEBUG_ALL = .true.

    ! Debug discrepancy option
    logical, parameter :: DEBUG_DISC = DEBUG_ALL .and. .true.

    ! Debug matrix A calculation option
    logical, parameter :: DEBUG_MATR = DEBUG_ALL .and. .true.
    logical, parameter :: DEBUG_MATR_INDICIES = DEBUG_ALL .and. .true.
end module DebugConfig