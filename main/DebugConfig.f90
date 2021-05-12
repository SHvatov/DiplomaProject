module DebugConfig
    implicit none

    ! Global debug option
    logical, parameter :: DEBUG_ALL = .true.

    ! Debug discrepancy option
    logical, parameter :: DEBUG_DISC = DEBUG_ALL .and. .true.

    ! Debug matrix A calculation option
    logical, parameter :: DEBUG_MATR = DEBUG_ALL .and. .true.

    ! Whether to print the matr or not
    logical, parameter :: OUTPUT_MATR = .true.
end module DebugConfig