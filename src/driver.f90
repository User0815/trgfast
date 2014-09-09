! Test program that uses the trgfast libary
! @author Adrian Vollmer
! @date 15.04.2014


program driver

    use tools
    use options
    use ode
    
    implicit none


    double precision, allocatable, dimension(:,:) :: pslin
    double precision, allocatable, dimension(:,:,:) :: psnl, Omega

    double precision, parameter :: OmegaM0 = .271d0, OmegaL0 = 1-OmegaM0, zini = 100d0
    double precision opts(N_options)
    integer i


        ! prepare background function
    allocate(Omega(100,1,3))
    do i=1,size(Omega,1)
        Omega(i,1,1) = (i-1)*log(1+zini)/(size(Omega,1)-1)
        Omega(i,1,2:) = O_LCDM(Omega(i,1,1))
    end do

    opts = 0d0
    opts(1) = 1d0 ! default options
    pslin = import("examples/pk.dat", 2)
    psnl = time_evolution((/0d0, 2d0, log(1+zini) /), pslin(:,1), pslin(:,2), Omega(:,1,1), Omega(:,:,2:), opts)
    call export(psnl(:,:,2), "psout.dat")

    deallocate(Omega)
    deallocate(psnl)
    deallocate(pslin)

contains 

        ! returns Omega_21 and Omega_22
    function O_LCDM(eta)
        implicit none
        double precision O_LCDM(2), a
        double precision, intent(in) :: eta
        a = exp(eta)/(1+zini)
        O_LCDM(1) = -1.5d0 * OmegaM0/(a**3*OmegaL0 + OmegaM0) 
        O_LCDM(2) = 2d0 - (-2d0*a**3*OmegaL0 + OmegaM0)/(2d0*(a**3*OmegaL0 + OmegaM0)) 
    end function O_LCDM

end program driver
