! The options are defined in here
! @author Adrian Vollmer
! @date 03.04.2014

! I know, this is not the prettiest solution, but it allows for options
! to be passed from C (and thus from Mathematica) simply as a flat
! double array.

module options

    ! Everything is public

    implicit none
    logical, save :: options_set = .false.
    logical  &
        oDEFAULT, &
        oWANT_PS, &
        oWANT_RAW, &
        oWANT_KDQ, &
        oWANT_DIFF, &
        oWANT_IM, &
        oWANT_A, &
        oLINEAR, &
        oSYMMETRY
    double precision &
        oK_LINEAR, &
        oK_TH, &
        oTOLERANCE, &
        oRKTHRESH, &
        oSCALE_A, &
        oCURV_LIMIT, &
        oSLOPE_LIMIT, &
        oOFFSET
    integer &
        oNq, &
        oDUMP_DIR, &
        oGROWTH_LENGTH, &
        oVERBOSITY, &
        oTIMESTEPS, &
        oEXTRA
    integer, parameter :: N_options = 23 ! number of options

    integer, parameter :: cPOWERLAW = 0
    integer, parameter :: cZERO = 1
    integer, parameter :: cEXPO = 2
    integer, parameter :: cCONST = 3

contains

! these compiler macros provide convenience when converting from double
#define set_integer(a) a = int(opts(i)); i = i+1 
#define set_logical(a) a = (opts(i)>0d0); i = i+1 
#define set_double(a) a = opts(i); i = i+1 


        ! convert an array to options
    subroutine set_options(opts)
        implicit none
        double precision, intent(in) :: opts(N_options)
        integer i
        
        options_set = .true.
        
        if (opts(1) > 0d0) then
            call set_default_options
        else
            i = 1
                ! this defines the order of the options in the array
            set_logical(oDEFAULT)
            set_double(oK_LINEAR)
            set_double(oK_TH)
            set_integer(oNq)
            set_double(oTOLERANCE)
            set_double(oRKTHRESH)
            set_integer(oTIMESTEPS) 
            set_logical(oLINEAR)
            set_integer(oGROWTH_LENGTH) ! must be < 500
            set_integer(oVERBOSITY)
            set_integer(oDUMP_DIR)
            set_logical(oWANT_PS)
            set_logical(oWANT_RAW)
            set_logical(oWANT_IM)
            set_logical(oWANT_DIFF)
            set_logical(oWANT_A)
            set_logical(oWANT_KDQ)
            set_double(oSCALE_A)
            set_double(oCURV_LIMIT)
            set_double(oSLOPE_LIMIT)
            set_logical(oSYMMETRY)
            set_integer(oEXTRA)
            set_double(oOFFSET)

            call validate_options
        end if

    end subroutine set_options


    subroutine set_default_options
        call set_options( (/ &
                -1d0, & !oDEFAULT
                1d-3, & !oK_LINEAR
                0.01d0, & !oK_TH
                200d0, & !oNq
                1d-2, & !oTOLERANCE
                1d-3, & !oRKTHRESH
                30d0, & !oTIMESTEPS
                -1d0, & !oLINEAR
                100d0, & !oGROWTH_LENGTH
                2d0, & !oVERBOSITY
                0d0, & !oDUMP_DIR
                -1d0, & !oWANT_PS
                -1d0, & !oWANT_RAW
                -1d0, & !oWANT_IM
                -1d0, & !oWANT_DIFF
                -1d0, & !oWANT_A
                -1d0, & !oWANT_KDQ
                1d0, & !oSCALE_A
                10d0, & !oCURV_LIMIT
                20d0, & !oSLOPE_LIMIT
                1d0, & !oSYMMETRY
                2d0, & !oEXTRA (0 = power law)
                1d-4 & ! oOFFSET
            /) )
    end subroutine set_default_options


    subroutine validate_options
        implicit none

        if (.not. ( 0d0 < oSCALE_A .and. oSCALE_A < 1d3)) &
            write(*,'(A)') "W Your oSCALE_A option is weird."

        ! TODO continue

    end subroutine validate_options

end module options
