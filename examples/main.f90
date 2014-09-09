!> @file main.f90
!! An example of how to use the trgfast library from within Fortran.
!! The two warnings about unused arguments in the background functions when
!! compiling can be ignored.
!! @author Adrian Vollmer
!! @date 20.07.2012

program trg

    use tools ! only needed for import/export functions
    use trgfast
    
    implicit none

    procedure(func_template), pointer :: pO21, pO22

    double precision, parameter :: h = 70.3d-2, OmegaM0 = .271d0,&
                OmegaL0 = 1-.271d0, bd_alpha = 1d0/(2*5d4+3), bd_mphi = 0d0

    double precision, allocatable, dimension(:,:) :: ps0, ps_out
    double precision par_real(6)
    integer par_int(12)

        ! load the linear power spectrum at present time from file
    ! call import(ps0, "pk.dat", 2)
    ! call import(ps0, "/home/vollmer/other/class/output/ref_z1_pk.dat", 2)
    call import(ps0, "/home/vollmer/projects/diss/fortran/examples/copter.ps.dat", 2)

        ! Bundle them into one array
    allocate(ps_out(4,size(ps0,2)))

        ! Set up the cosmological parameters... could be read from a
        ! config file.
    par_real(1) = 100d0
    par_real(2) = 0d0
    par_real(3) = -1.d0 !.999684! f(z)
    par_real(4) = .966d0 ! ns
    par_real(5) = 17d0 ! k_th
    par_real(6) = 1d0 ! scale_A
    par_int(1) = size(ps0,2)
    par_int(2) = 1 ! full nonlinear
    par_int(3) = 1 ! verbosity
    par_int(4) = 50 ! time evolution step size
    par_int(5) = 200 ! sampling points for Kdq
    par_int(6) = 1  ! extrapolation method
    par_int(7) = 0  ! dump Kdq?
    par_int(8) = 0  ! dump A?
    par_int(9) = 1  ! dump PS?
    par_int(10) = 0 ! dump diff?
    par_int(11) = 0 ! dump growth?
    par_int(12) = 0 ! initial conditions

        ! Make procedure pointers...
    pO21 => O21_LCDM
    pO22 => O22_LCDM
 
        ! Call the function from the TRG module
    call trgfast_init(ps0, par_real, par_int, pO21, pO22, "../output/")
        ! Request a PS
    call trgfast_get_ps(0d0, curlyH(1d0/(1+0)), ps_out)
        ! Write the portion of the resulting array that contains actual data
        ! into a file
    call export(ps_out, "../output/ps_nl_full.dat")
        ! Clean up
    call trgfast_free
    

    deallocate(ps0,ps_out)

contains
    
    function curlyH(a)
        implicit none
        double precision curlyH
        double precision, intent(in) :: a
        curlyH = a*(h/3d3)*sqrt(OmegaL0 + OmegaM0/a**3);
    end function curlyH

    function O21_LCDM(a, k)
        implicit none
        double precision O21_LCDM
        double precision, intent(in) :: a, k
        O21_LCDM = -1.5d0 * OmegaM0/(a**3*OmegaL0 + OmegaM0) + 0*k
        ! the 0*k is just to suppress warnings about unused dummy variables
    end function O21_LCDM

    function O22_LCDM(a, k)
        implicit none
        double precision O22_LCDM
        double precision, intent(in) :: a, k
        O22_LCDM = 2d0 - (-2d0*a**3*OmegaL0 + OmegaM0)/(2d0*(a**3*OmegaL0 + OmegaM0)) + 0*k
        ! the 0*k is just to suppress warnings about unused dummy variables
    end function O22_LCDM

    function O21_bransdicke(a, k)
        implicit none
        double precision O21_bransdicke
        double precision, intent(in) :: a, k
        O21_bransdicke = O21_LCDM(a,k) + 0d0*bd_alpha * a ! alpha*dphi/dlog(a)
    end function O21_bransdicke

    function O22_bransdicke(a, k)
        implicit none
        double precision O22_bransdicke
        double precision, intent(in) :: a, k
        O22_bransdicke = O22_LCDM(a,k)*(1 + 2*bd_alpha**2/(1+(bd_mphi*a/k)**2))
    end function O22_bransdicke


end program trg
