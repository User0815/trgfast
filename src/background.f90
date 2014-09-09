! Defining the background funtions
! @author Adrian Vollmer
! @date 09.04.2014


module background

    use tools
    use options

private

    double precision, allocatable :: eta_(:), k_(:)
    double precision, allocatable :: O21_int(:,:,:), O22_int(:,:,:)

public :: Omega21_int, Omega22_int, init_background, clean_up_background

contains 
    
    subroutine init_background(eta, OmegaBulk, k)
        implicit none
        double precision, intent(in) :: eta(:), OmegaBulk(:,:,:)
        double precision, intent(in), optional :: k(:)
        double precision, allocatable :: intobj(:,:)
        integer i,n
        ! TODO warning if number of Omega values doesn't match eta_ or k
        
        allocate(eta_(size(eta)))
        eta_ = eta

        n = 1
        if (present(k)) then
            n = size(k)
            allocate(k_(n))
            k_ = k
        end if

        allocate(O21_int(5,size(eta_),n))
        allocate(O22_int(5,size(eta_),n))

        do i=1,n
            intobj = init_interpolation(eta, OmegaBulk(:,i,1), cSPLINE)
            O21_int(:,:,i) = intobj(:,:)
            deallocate(intobj)
            intobj = init_interpolation(eta, OmegaBulk(:,i,2), cSPLINE)
            O22_int(:,:,i) = intobj
            deallocate(intobj)
        end do

    end subroutine init_background


    subroutine clean_up_background()
        if (allocated(eta_)) deallocate(eta_)
        if (allocated(k_)) deallocate(k_)
        if (allocated(O21_int)) deallocate(O21_int)
        if (allocated(O22_int)) deallocate(O22_int)
    end subroutine clean_up_background


    function Omega21_int(eta, k)
        implicit none
        double precision, intent(in) :: eta, k
        integer idx, ki
        double precision Omega21_int
        idx = -1
        ki = 1
        if (allocated(k_)) ki = bsearch(k_, k)
        Omega21_int = eval_interpolation(O21_int(:,:,ki), eta, cSPLINE, idx) 
        return
    end function Omega21_int


    function Omega22_int(eta, k)
        implicit none
        double precision, intent(in) :: eta, k
        integer idx, ki
        double precision Omega22_int
        idx = -1
        ki = 1
        if (allocated(k_)) ki = bsearch(k_, k)
        Omega22_int = eval_interpolation(O22_int(:,:,ki), eta, cSPLINE, idx)
        return
    end function Omega22_int


end module background
