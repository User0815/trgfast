! Solving the TRG ODE 
! @author Adrian Vollmer
! @date 12.03.2014


module ode

    use tools
    use options
    use nonlinear
    use background
    use iso_c_binding
    use omp_lib
    implicit none

private 
    integer counter, ki_linear
    double precision, allocatable :: k_(:), logk_(:), growth2(:,:) ! growth2 stores G^2(eta)
  
public :: time_evolution

        ! interface to the F77 routines from netlib.org
    interface
        subroutine setup(neq,tstart,ystart,tend,tol,thres,method,task,&
                         errass,hstart,work,lenwrk,mesage)
          double precision  hstart, tend, tol, tstart, thres(*), work(*), ystart(*)
          integer           lenwrk, method, neq
          logical           errass, mesage
          character*(*)     task
        end subroutine setup

        subroutine ut(f,twant,tgot,ygot,ypgot,ymax,work,uflag)
          double precision  tgot, twant, work(*), ygot(*), ymax(*), ypgot(*), &
                            tstrt, tnd, dir, hstrt, tolr
          external          f
          integer           uflag, neqn
        end subroutine ut
    end interface

contains

    subroutine dump_A(A, k)
        implicit none
        double precision, intent(in) ::  A(:,:), k(:)
        double precision, allocatable :: table(:,:)
        character(70) filename

        allocate(table(ki_max,15))
        table(:,1) = k(:ki_max)
        table(:,2:) = A(:ki_max,:)
        write(filename,'(A,i3.3,A)') dump_dir()//"A-", counter, ".dat"
        call export(table, filename)
        deallocate(table)
    end subroutine dump_A


        ! allocate global arrays, set options, etc
        ! Omega_abc = Omega_c(eta_a, k_b) where c={21,22}
    subroutine init_ode(k, O_eta_, OmegaBulk, opts, O_k_)
        implicit none
        double precision, intent(in) :: OmegaBulk(:,:,:), O_eta_(:), opts(N_options), k(:)
        double precision, intent(in), optional :: O_k_(:)

        call set_options(opts)

        if (present(O_k_)) then
            ! compiler warns here, but it's fine. the compiler doesn't
            ! understand optional arguments.
            call init_background(O_eta_, OmegaBulk, O_k_)
        else
            call init_background(O_eta_, OmegaBulk)
        end if

        ki_max = size(k)
        if (ki_max < 100) call message("w", &
                msg="Low sampling rate of the power spectrum, might not be sufficient.",&
                subname="init_ode")
        if (ki_max > 1500) call message("w",&
                msg="High sampling rate of the power spectrum, might unnecessarily increase computation time.",&
                subname="init_ode")

        allocate(k_(ki_max))
        k_ = k
        allocate(logk_(ki_max))
        logk_ = log(k)
        ki_linear = bsearch(k_, oK_LINEAR)
        allocate(growth2(oGROWTH_LENGTH,2))
    end subroutine init_ode


        ! deallocate global arrays
    subroutine clean_up_ode()
        if (allocated(k_)) deallocate(k_)
        if (allocated(logk_)) deallocate(logk_)
        if (allocated(growth2)) deallocate(growth2)
        call clean_up_background
    end subroutine clean_up_ode


        ! this is the f in x'=f(x,t)
    subroutine f_ode(eta,X,Xprime) 
        implicit none
        double precision, intent(in) ::  X(17*ki_max), eta
        double precision, intent(out) ::  Xprime(17*ki_max)
        double precision yprime_(17), Omega21, Omega22, A_(14), &
            expeta, y_(17),  k, A_all_k(ki_max,14)
        double precision, allocatable :: logP(:,:)
        integer ki
        character(200) string

        
        counter = counter + 1
        allocate(logP(ki_max,4))
        expeta = exp(eta)
        logP(:,1) = logk_(:ki_max) 
        logP(:,2) = (X(15:17*ki_max:17))
        logP(:,3) = (X(16:17*ki_max:17))
        logP(:,4) = (X(17:17*ki_max:17))

        write(string,'(A,G12.5,A,i4.4,A,i4.4)') "eta = ", eta, ", N = ", &
                    ki_max, ", counter = ", counter
            ! don't print this for the linear run
        if (ki_max>1) call message("i", msg=trim(string), subname="f_ode")

            ! compute A(k) for all k
        if (oLINEAR) then
            A_all_k = 0d0
        else
            call init_A(logP)
            A_all_k = 0d0
            ! Fortunately, the problem is "embarrassingly parallel"
!$OMP PARALLEL DO 
            do ki=ki_linear,ki_max
                A_all_k(ki,:) = A(k_(ki))
            end do
            if (oWANT_A) call dump_A(A_all_k, k_)
            call dump_data
            call clean_up
        end if

            ! compute right hand side of eq. B.4, B.5 (Pietroni 2008)
        Xprime(ki_max+1:) = 0d0
        do ki=1,ki_max
            k = k_(ki)
            Omega21 = Omega21_int(eta, k)
            Omega22 = Omega22_int(eta, k)
            y_ = X(ki*17-16:ki*17)
            y_(15:) = exp(logP(ki,2:))
            A_ = A_all_k(ki,:)
            include "src/mathematica_rhs"
            Xprime(ki*17-16:ki*17-3) = yprime_(:14)
            Xprime(ki*17-2:ki*17) = yprime_(15:)/y_(15:) 
                ! division because we're solving the ode in log-space
        end do
        deallocate(logP)
    end subroutine f_ode

    
        ! computes the derivative of the growth at eta
    function growth_factor(eta)
        implicit none
        double precision, intent(in) :: eta
        double precision growth_factor
        integer i
        i = min(bsearch(growth2(:,1), eta), size(growth2,1)-1)
        growth_factor = abs(log(growth2(i+1,2)/growth2(i,2)) / &
                (growth2(i+1,1) - growth2(i,1))) / 2
                ! factor of 2 because this is G^2(eta). Overall amplitude
                ! cancels out
        return
    end function growth_factor


        ! 4th order Runge-Kutta. Actual code from netlib.org
    subroutine rk_solver(eta_ini, eta_fin_, X, f, ps_out, sizes)
        implicit none
        double precision, intent(in) :: eta_ini, eta_fin_(:)
        double precision, intent(in) :: X(17*ki_max)
        double precision, intent(inout) :: ps_out(ki_max,4,size(eta_fin_))
        integer, intent(out) :: sizes(size(eta_fin_))
        double precision thresh(17*ki_max), work(17*ki_max*21+100), Xprime(17*ki_max), &
                        Xmax(17*ki_max), tolerance, hstart, eta_result
        logical messages, error_ass
        integer method, uflag, i
        external f

            ! options for RK4
        tolerance = oTOLERANCE
        method = 2
        hstart = 0.1d0 
        messages = .true.
        error_ass = .false.
        thresh = oRKTHRESH

        Xprime = 0d0
        call setup(17*ki_max, eta_ini, X, eta_fin_(size(eta_fin_)), tolerance, thresh, method, &
                "u", error_ass, hstart, work, size(work), messages)
        do i=1,size(eta_fin_)
            call ut(f, eta_fin_(i), eta_result, X, Xprime, Xmax, work, uflag)
            sizes(i) = ki_max
            ps_out(:ki_max,1,i) = k_(:ki_max)
            ps_out(:ki_max,2,i) = exp(X(15:17*ki_max:17))
            ps_out(:ki_max,3,i) = exp(X(16:17*ki_max:17))
            ps_out(:ki_max,4,i) = exp(X(17:17*ki_max:17))
        end do
    end subroutine rk_solver
    

        ! 4th order Runge-Kutta, non-adaptive
    subroutine rk_nonad_solver(eta_ini, eta_fin_, Xini, f, ps_out, sizes)
        implicit none
        double precision, intent(in) :: eta_ini, eta_fin_(:)
        double precision, intent(in) :: Xini(17*ki_max)
        double precision, intent(inout) :: ps_out(ki_max,4,size(eta_fin_))
        double precision, dimension(17*ki_max) :: X, Xold, k1, k2, k3, k4
        double precision eta_f, eta, etaold, h
        integer, intent(out) :: sizes(size(eta_fin_))
        integer i,steps,etai
        external f
        
        X = Xini
        if (ki_max>1) then
            steps = oTIMESTEPS
        else
            steps = 300 ! only linear run for growth function
        end if
        eta_f = eta_fin_(size(eta_fin_))
        h = (eta_f-eta_ini)/(steps)
        etai = 1
        eta = eta_ini

        do i=1,steps+1 ! one extra step just to make sure 
            etaold = eta
                ! save PS for interpolation if necessary
            if (eta+h>=eta_fin_(etai)) Xold = X
                ! apply RK expression
            call f(eta,X,k1)
            call f(eta+h/2,X+h/2*k1,k2)
            call f(eta+h/2,X+h/2*k2,k3)
            call f(eta+h,X+h*k3,k4)
            eta = eta + h
            X = X + h/6*(k1+2*(k2+k3)+k4)

            if (eta>=eta_fin_(etai)) then ! output PS at this redshift
                sizes(etai) = ki_max
                    ! linearly interpolate 
                Xold = Xold + (X-Xold)*(eta_fin_(etai) - etaold) / (eta - etaold)
                ps_out(:ki_max,1,etai) = k_(:ki_max)
                ps_out(:ki_max,2,etai) = exp(Xold(15:17*ki_max:17))
                ps_out(:ki_max,3,etai) = exp(Xold(16:17*ki_max:17))
                ps_out(:ki_max,4,etai) = exp(Xold(17:17*ki_max:17))
                etai = etai + 1
            end if
                ! Are we done computing the ps at all requested redshifts?
            if (etai > size(eta_fin_)) exit
        end do

    end subroutine rk_nonad_solver


        ! do a linear run with just one sampling point to obtain the growth
        ! function
    subroutine get_growth(eta_ini, eta_fin)
        implicit none
        double precision, intent(in) :: eta_ini, eta_fin
        double precision X(17), ps_out(1,4,size(growth2,1))
        double precision, allocatable :: eta_(:)
        integer ki_max_temp, sizes(size(growth2,1))
        logical linear_temp

        ki_max_temp = ki_max
        ki_max = 1
        linear_temp = oLINEAR
        oLINEAR = .true.
        
        X = 0d0 
        eta_ = aspace(eta_ini+0.01d0, eta_fin, size(growth2,1))
        if (oTIMESTEPS > 1) then
            call rk_nonad_solver(eta_ini, eta_, X(:17), f_ode, ps_out, sizes)
            ! call rk_solver(eta_ini, eta_, X(:17), f_ode, ps_out, sizes)
        else
            call rk_solver(eta_ini, eta_, X(:17), f_ode, ps_out, sizes)
        end if
        growth2(:,1) = eta_
        growth2(:,2) = ps_out(1,2,:) * exp(eta_)**2

        ki_max = ki_max_temp
        oLINEAR = linear_temp
        deallocate(eta_)
        return
    end subroutine get_growth


        ! evolve the power spectrum in time. this is the main routine.
    function time_evolution(eta_, k_, ps_, O_eta_, OmegaBulk, opts, growth_out, O_k_) result(ps_out)
        implicit none
        double precision, intent(in) :: eta_(:), k_(:), ps_(:), O_eta_(:), &
                        OmegaBulk(:,:,:), opts(N_options)
        double precision, intent(out), allocatable, optional ::  growth_out(:,:)
        double precision, intent(in), optional ::  O_k_(:)
        double precision X(17*size(k_)), ps_out(size(k_),4,(size(eta_)-1))
        integer sizes_(size(eta_)-1), i
        character(200) string
        
            ! initialization
        if (present(O_k_)) then
            call init_ode(k_, O_eta_, OmegaBulk, opts, O_k_)
            ! compiler warns here, but it's fine. the compiler doesn't
            ! understand optional arguments.
        else
            call init_ode(k_, O_eta_, OmegaBulk, opts)
        end if

            ! linear run
        call get_growth(eta_(1), eta_(size(eta_)))
        if (present(growth_out)) then
                ! return growth function
            allocate(growth_out(size(growth2,1),2))
                ! convert from eta to z
            growth_out(:,1) = exp(eta_(size(eta_))-growth2(:,1))-1
            growth_out(:,2) = sqrt(growth2(:,2)/growth2(size(growth2,1),2))
                ! reverse order
            growth_out(:,:) = growth_out(size(growth_out,1):1:-1,:)
        end if

            ! initial conditions
        X = 0d0 ! bispectrum vanishes
        X(15::17) = log(ps_/growth2(size(growth2,1),2))
        X(16::17) = X(15::17)
        X(17::17) = X(15::17)
        ps_out = -1d0

        counter = 0
        if (oTIMESTEPS > 0) then
            call rk_nonad_solver(eta_(1), eta_(2:), X, f_ode, ps_out, sizes_)
        else
            call rk_solver(eta_(1), eta_(2:), X, f_ode, ps_out, sizes_)
        end if

            !  return result
        do i=1,size(eta_)-1
            ps_out(:sizes_(i),2,i) = ps_out(:sizes_(i),2,i) * exp(eta_(i+1))**2
            ps_out(:sizes_(i),3,i) = ps_out(:sizes_(i),3,i) * exp(eta_(i+1))**2 / growth_factor(eta_(i+1))
            ps_out(:sizes_(i),4,i) = ps_out(:sizes_(i),4,i) * exp(eta_(i+1))**2 / growth_factor(eta_(i+1))**2
            ! extra factor due to conversion to physical power spectrum 
        end do

        
        write(string,"(A,i5.3)") "Number of evaluations of A: ", counter
        call message("i", msg=string, subname="time_evolution")
        call clean_up_ode
        return
    end function time_evolution



!=============== wrappers for C ====================
! since we're passing flat arrays for compatibility, we need to reshape them.
! After that, Fortran takes over immediately.

    subroutine time_evolution_c(eta_, eta_len, ps_in, ps_len, O_eta, O_eta_len, &
                    O_k, O_k_len, OmegaBulk, OBulk_len, ps_out, growth_out, opts) bind(c)
        implicit none
        integer(c_int), intent(in) :: OBulk_len, ps_len, eta_len, O_eta_len, O_k_len
        real(c_double), intent(in) :: eta_(eta_len), ps_in(ps_len), OmegaBulk(OBulk_len), &
                                        opts(N_options), O_eta(O_eta_len), O_k(O_k_len)
        real(c_double), intent(out) :: ps_out(2*ps_len*(eta_len-1)), growth_out(1000)
        integer i
        double precision, allocatable :: ps_out_nonflat(:,:,:), growth_nonflat(:,:)
        double precision ps_in_nonflat(ps_len/2,2), OBulk_nonflat(OBulk_len/2/O_k_len,O_k_len,2)

        OBulk_nonflat = reshape(OmegaBulk,shape(OBulk_nonflat))

        ps_in_nonflat(:,1) = ps_in(:ps_len/2)
        ps_in_nonflat(:,2) = ps_in(ps_len/2+1:)

        if (O_k_len > 1) then
            ps_out_nonflat = time_evolution(eta_, ps_in_nonflat(:,1), ps_in_nonflat(:,2), &
                            O_eta, OBulk_nonflat, opts, growth_nonflat, O_k)
        else
            ps_out_nonflat = time_evolution(eta_, ps_in_nonflat(:,1), ps_in_nonflat(:,2), &
                            O_eta, OBulk_nonflat, opts, growth_nonflat)
        end if
        growth_out = -1d0
        growth_out(:2*size(growth_nonflat,1)) = [growth_nonflat(:,1), growth_nonflat(:,2)]

            !Flatten result
        do i=1,eta_len-1
            ps_out(2*ps_len*(i-1)+1:2*ps_len*i) = [&
                            ps_out_nonflat(:,1,i), &
                            ps_out_nonflat(:,2,i), &
                            ps_out_nonflat(:,3,i), &
                            ps_out_nonflat(:,4,i)]
        end do

        deallocate(ps_out_nonflat)
        deallocate(growth_nonflat)
    end subroutine time_evolution_c


    subroutine init_ode_c(k, k_len, O_eta, O_eta_len, OmegaBulk, OBulk_len, O_k, O_k_len, opts) bind(c)
        implicit none
        integer(c_int), intent(in) :: OBulk_len, k_len, O_eta_len, O_k_len
        real(c_double), intent(in) :: O_eta(O_eta_len), OmegaBulk(OBulk_len), &
                        k(k_len), opts(N_options), O_k(O_k_len)
        double precision Omega_nonflat(OBulk_len/2/k_len,k_len,2)
        Omega_nonflat = reshape(OmegaBulk,shape(Omega_nonflat))
        if (k_len > 1) then
            call init_ode(k, O_eta, Omega_nonflat, opts, O_k)
        else
            call init_ode(k, O_eta, Omega_nonflat, opts)
        end if
    end subroutine init_ode_c


    subroutine clean_up_ode_c() bind(c)
        call clean_up_ode
    end subroutine clean_up_ode_c


    subroutine f_ode_c(eta, X, Xprime) bind(c)
        implicit none
        real(c_double), intent(in) ::  X(17*ki_max), eta
        ! integer, intent(in) :: n
        real(c_double), intent(out) ::  Xprime(17*ki_max)
        call f_ode(eta, X, Xprime)
    end subroutine f_ode_c

end module ode
