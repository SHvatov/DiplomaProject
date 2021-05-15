module DebugConfig
    implicit none

    ! Global debug option
    logical, parameter :: DEBUG_ALL = .false.

    ! Debug discrepancy option
    logical, parameter :: DEBUG_DISC = DEBUG_ALL .and. .true.

    ! Debug matrix A calculation option
    logical, parameter :: DEBUG_MATR = DEBUG_ALL .and. .true.

    ! Whether to print the matr or not
    logical, parameter :: OUTPUT_MATR = .false.

    ! Whether to print main data
    logical, parameter :: DEBUG_MAIN = .true.
end module DebugConfig