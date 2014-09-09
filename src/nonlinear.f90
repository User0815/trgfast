! All the functions to compute A(k).
! @author Adrian Vollmer
! @date 12.03.2014


module nonlinear

    use tools
    use options
    use iso_c_binding
    implicit none

private 
  
  double precision, allocatable :: P_ab_splines(:,:,:), logP_ab_array(:,:), &
          intg_moments_splines(:,:,:)
  double precision :: ns, exponents(3), logkmax
  integer, save :: counter = 0 ! counter for the dump files
  integer ki_max ! this variable is needed for the ODE solver to know when to stop.
                 ! also, for ki>ki_max the power spectrum will be extrapolated

public :: init_A, A, clean_up, dump_data, ki_max

contains

    subroutine init_P
        implicit none
        integer i,j
        double precision curv(3,2)

        if (.not.allocated(logP_ab_array)) then
            call message("c", msg="Power spectrum values not available", &
                    subname="init_P")
        end if

        allocate(P_ab_splines(5,size(logP_ab_array,1),3))

        do i=1,3
            P_ab_splines(:,:,i) = init_interpolation(logP_ab_array(:,1),&
                logP_ab_array(:,i+1), cSPLINE)
        end do

        do i=size(logP_ab_array,1),20,-1
            do j=1,size(curv,1)
                curv(j,:) = abs(P_ab_splines(4,i-size(curv,2)+1:i,j))
            end do
            if (stddev(curv(1,:)) < oCURV_LIMIT .and.&
                stddev(curv(2,:)) < oCURV_LIMIT .and.&
                stddev(curv(3,:)) < oCURV_LIMIT &
                .and. abs(P_ab_splines(3,i,3)) < oSLOPE_LIMIT&
                ) exit
        end do
        
            ! save the exponents at the large k end for extrapolation
            ! later. the exponent in logspace is the slope, which is the second
            ! spline coefficient, i.e. the third element
        exponents = P_ab_splines(3,i,:)
        logkmax = logP_ab_array(i,1)
        ki_max = i
            ! compute ns. P(k) ~ k^ns for small k
        ns = P_ab_splines(3,1,1) ! this is a bad idea if the power
                    ! spectrum curves upwards because of the universe is closed.
                    ! so we fix it just in case:
        if (ns < 0d0 ) ns = 1d0
    end subroutine init_P


        ! extrapolates the power spectra
    function extrapolation(P_ab, k, logk)
            implicit none
            double precision, intent(in) :: P_ab(3), logk, k
            double precision extrapolation(3)
            if (oEXTRA .eq. cPOWERLAW) then
                extrapolation = exp(P_ab) * (k/exp(logkmax))**exponents 
            elseif (oEXTRA .eq. cZERO) then
                extrapolation = 0d0
            elseif (oEXTRA .eq. cEXPO) then
                extrapolation = exp(P_ab - 4*((k - exp(logkmax))/exp(logkmax))**2)
            elseif (oEXTRA .eq. cCONST) then
                extrapolation = exp(P_ab)
            else
                call message("e", msg="Unkown extrapolation method", &
                        subname="extrapolation")
            end if
    end function extrapolation


        ! returns the interpolated (or extrapolated) value of the power spectra
    function P_ab(k)
        implicit none
        double precision, intent(in) :: k
        double precision :: P_ab(3), logk, logkmin
        integer indx, j
        
        if (.not.allocated(P_ab_splines)) then
            call message("c", msg="Power spectrum not interpolated", &
                        subname="P_ab")
        end if

        logk = log(k)
        logkmin = P_ab_splines(1,1,1)
        ! logkmax is global

        if (logk > logkmax) then
                ! extrapolation for large k
            indx = -1
            P_ab = (/( eval_interpolation(P_ab_splines(:,:,j), &
                    logkmax, cSPLINE, indx), j=1,3 )/)
            P_ab = extrapolation(P_ab, k, logk)
        else if (logk < logkmin) then
                ! extrapolation for small k
            indx = 1
            P_ab = (/( eval_interpolation(P_ab_splines(:,:,j), &
                    logkmin, cSPLINE, indx), j=1,3 )/)
            P_ab = exp(P_ab) * (k / exp(logkmin))**(ns)
        else 
            indx = -1
            P_ab = (/( eval_interpolation(P_ab_splines(:,:,j), &
                    logk, cSPLINE, indx), j=1,3 )/)
            P_ab = exp(P_ab)
        end if

        return
    end function P_ab


    subroutine init_integrated_moment()
        implicit none
        double precision, allocatable :: ps(:,:), intobj(:,:)
        integer i,extra

        extra = 30
        allocate(intg_moments_splines(6,size(P_ab_splines,2)+extra,13))
        allocate(ps(size(logP_ab_array,1)+extra,4))
        
        ps(:size(ps,1)-extra,1) =  logP_ab_array(:,1)
        ps(:size(ps,1)-extra,2:) = exp(logP_ab_array(:,2:))

        do i=size(ps,1)-extra+1,size(ps,1)
                ! extrapolate integrated moments as well
            ps(i,1) = ps(size(ps,1)-extra,1) + 10d0*(i-size(ps,1)+extra)/extra
            ps(i,2:) = P_ab(exp(ps(i,1)))
        end do

        do i=1,13
                ! the 1st argument is in log-space, the 2nd is not
            intobj = init_interpolation(ps(:,1), ps(:,i/5+2), cSPLINE)
            call integrate_splines(intobj, intg_moments_splines(:,:,i),&
                    -3+modulo(i,5)*2)
            if (allocated(intobj)) deallocate(intobj)
        end do

        deallocate(ps)

    end subroutine init_integrated_moment


    function eval_integrated_moment(i, k, indx) result(y)
        implicit none
        double precision, intent(in) :: k
        integer, intent(in) :: i
        integer, intent(inout) :: indx
        double precision y, m(13), a, expa, b, logk
        integer dummyindx,j 

        m = (/(-3+modulo(j,5)*2, j=1,13)/) + ns + 1 

        a = (intg_moments_splines(1,1,i))
        b = (intg_moments_splines(1,size(intg_moments_splines,2),i))
        logk = log(k)
        if (logk < a) then
            ! extrapolate for small k
            dummyindx = 1
            expa = exp(a)
            ! integrate it analytically (and adjust constant of integration)
            y = (k**m(i) - expa**m(i))/m(i) * exp(logP_ab_array(1,i/5+2))  + &
                    eval_int_splines(intg_moments_splines(:,:,i), expa, a, dummyindx)
        else if (logk > b) then
            ! extrapolate for large k
            dummyindx = -1
            y = eval_int_splines(intg_moments_splines(:,:,i), exp(b), b, dummyindx)
        else
            y = eval_int_splines(intg_moments_splines(:,:,i), k, logk, indx)
        end if
        
        return
    end function eval_integrated_moment


    function diff(a, b)
        implicit none
        double precision, intent(in) :: a, b
        double precision diff(13), upper, lower
        integer i, d1, d2

        if (.not.allocated(intg_moments_splines)) then
            ! call init_integrated_moment
            call message("c", msg="Moments not integrated", &
                        subname="diff")
        end if

        d1 = -1
        d2 = -1
        do i=1,13
            upper = eval_integrated_moment(i, b, d2)
            lower = eval_integrated_moment(i, a, d1)
            diff(i) = upper - lower
        end do
        return
    end function diff


    function Kdq(k, q)
        implicit none
        double precision, intent(in) :: k, q
        double precision P11q, P12q, P22q
        double precision P11k, P12k, P22k
        double precision k2,k3,k4,k5,logqk
        double precision Kdq(14), diff_(13), ps(3)

            ! Abbreviations
        k2 = k*k; k3 = k*k2; k4 = k*k3; k5 = k*k4; 
        logqk = Log(q/abs(k-q))

        ps = P_ab(k)
        P11k = ps(1)
        P12k = ps(2)
        P22k = ps(3)

        ps = P_ab(q)
        P11q = ps(1)
        P12q = ps(2)
        P22q = ps(3)

        diff_ = diff(abs(q-k), q)

        ! Include mathematica expressions. The file consists of nothing but
        !
        ! Kdq(1) = ....
        ! Kdq(2) = ....
        ! 
        ! and so on

        include "src/mathematica_kdq"

        return 
    end function Kdq


        ! this function defines the sampling points for the q-integration
    function q_sampling(k)
        implicit none
        double precision, intent(in) :: k
        double precision, allocatable :: q_sampling(:), q_left(:), q_right(:)

        double precision left, right, middle, offset

        left = k/2 ! lower limit of q-Integral
        right = 10*k ! upper limit of q-Integral
            ! the upper limit is somewhat arbitrary, but we MUST HAVE left < right
            ! TODO Experiment with this?
        if (k < oK_TH) then 
            q_sampling = alogspace(left, right, oNq)
        else
            offset = oOFFSET !how close do we get to the zero? 
            middle = k
            q_left = alogspace(offset, middle - left, oNq/2)
            q_left = middle - q_left(oNq/2:1:-1)
            q_right = alogspace(offset, right - middle, oNq-oNq/2)
            q_right = q_right + middle
            allocate(q_sampling(oNq))
            q_sampling(:oNq/2) = q_left(:)
            q_sampling(oNq/2+1:) = q_right(:)
            deallocate(q_left)
            deallocate(q_right)
        end if

        return
    end function q_sampling


    function A(k)
        implicit none
        double precision, intent(in) :: k
        double precision A(14)
        integer j

        double precision, allocatable :: q_array(:)
        double precision Kdq_array(oNq,14)
        integer qi

        q_array = q_sampling(k)
    
        do qi = 1,oNq
            Kdq_array(qi,:) = Kdq(k, q_array(qi))
        end do

        if (oWANT_KDQ) call dump_Kdq(q_array, Kdq_array, k)

        do j=1,14
            A(j) = sum( (/( (Kdq_array(qi+1,j) + Kdq_array(qi,j)) * &
                            (q_array(qi+1)-q_array(qi)), &
                        qi=1,oNq-1 )/) ) / (2*TwoPiCubed) * oSCALE_A
                ! the TwoPiCubed comes from the chosen convention of the Fourier
                ! transform. The 2 is from the trapezoidal rule.
        end do
        deallocate(q_array)
        return
    end function A


    subroutine dump_Kdq(q_, kdq_, k)
        implicit none
        character(70) filename
        double precision, intent(in) :: kdq_(:,:), q_(:), k
        double precision, allocatable :: table(:,:)
        integer, save :: kcount = 1

        if (modulo(kcount,size(logP_ab_array,1)).eq.0) then
            allocate(table(size(kdq_,1),15))
            table(:,1) = q_
            table(:,2:) = kdq_
            write(filename,'(A,i3.3,A,E8.3,A)') dump_dir()//"kdq-", counter, "-", k, ".dat"
            call export(table, filename)
            deallocate(table)
        end if
        kcount = kcount + 1
    end subroutine dump_Kdq


    subroutine dump_data
        implicit none
        character(70) filename
        double precision, allocatable :: table(:,:), k_(:)
        double precision q
        integer i,j,n,idx,qi

        if (oWANT_RAW) then
            write(filename,'(A,i3.3,A)') dump_dir()//"rawps-", counter, ".dat"
            call export(exp(logP_ab_array), filename)
        end if

        n = 4*size(logP_ab_array,1)
        k_ = alogspace(.5d0*exp(logP_ab_array(1,1)), &
                9d0*exp(logP_ab_array(size(logP_ab_array,1),1)), n)

        if (oWANT_PS) then
            write(filename,'(A,i3.3,A)') dump_dir()//"ps-", counter, ".dat"
            allocate(table(n,4))
            table(:,1) = k_
            do i=1,n
                table(i,2:) = P_ab(k_(i))
            end do
            call export(table, filename)
            deallocate(table)
        end if

        if (oWANT_IM) then
            write(filename,'(A,i3.3,A)') dump_dir()//"im-", counter, ".dat"
            allocate(table(n,14))
            table(:,1) = k_
            do i=1,n
                do j=2,14
                    idx = -1
                    table(i,j) = eval_integrated_moment(j-1, k_(i), idx)
                end do
            end do
            call export(table, filename)
            deallocate(table)
        end if

        if (oWANT_DIFF) then
            allocate(table(n,14))
            do qi=0,6
                q = .001d0 * 10**(.5d0*qi)
                write(filename,'(A,i3.3,A,E8.3,A)') dump_dir()//"diff-", counter,"-", q, ".dat"
                table(:,1) = k_
                do i=1,n
                    table(i,2:) = diff(abs(q-k_(i)),k_(i)) ! q and k switched
                end do
                call export(table, filename)
            end do
            deallocate(table)
        end if 
        deallocate(k_)

    end subroutine dump_data


        ! initialize all the arrays for this time step
    subroutine init_A(logP, opts)
        implicit none
        double precision, intent(in) :: logP(:,:)
        double precision, intent(in), optional :: opts(N_options)

        if (present(opts)) then
            call set_options(opts)
        else
            if (.not.options_set) call message("w", &
                msg="Options may not be set.", subname = "init_A")
        end if
        

        if (size(logP,2).ne.4) then
            call message("c", &
                msg="Power spectrum has the wrong dimensions", &
                subname="A") 
        end if

        if (.not.allocated(logP_ab_array)) then
            allocate(logP_ab_array(size(logP,2),4))
            logP_ab_array = logP
        end if

        if (.not.allocated(P_ab_splines)) then
            call init_P
        end if

        if (.not.allocated(intg_moments_splines)) then
            call init_integrated_moment
        end if

        counter = counter + 1
    end subroutine


    subroutine clean_up
        if (allocated(P_ab_splines)) deallocate(P_ab_splines)
        if (allocated(logP_ab_array)) deallocate(logP_ab_array)
        if (allocated(intg_moments_splines)) deallocate(intg_moments_splines)
    end subroutine clean_up


!================== wrappers for C ===================
! NOTE the function names will all be converted to lower case!

    subroutine clean_up_c() bind(c)
        call clean_up
    end subroutine clean_up_c


    subroutine init_A_c(P, k_len, opts) bind(c)
        implicit none
        integer(c_int), intent(in) :: k_len
        real(c_double), intent(in) :: P(k_len), opts(N_options)
        double precision, allocatable :: ps(:,:)

        allocate(ps(k_len/4,4))
        ps(:,1) = P(:k_len/4)
        ps(:,2) = P(k_len/4+1:k_len/4*2)
        ps(:,3) = P(k_len/4*2+1:k_len/4*3)
        ps(:,4) = P(k_len/4*3+1:)

        call init_A(ps, opts)
        deallocate(ps)

    end subroutine init_A_c


    subroutine A_c(k, res) bind(c)
        implicit none
        real(c_double), intent(in) :: k
        real(c_double), intent(out) :: res(14)
        res = A(k)
    end subroutine A_c


end module nonlinear
