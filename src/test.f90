! This file contains test routines for the other modules.
! Adrian Vollmer


module test_mod

    use tools
    ! use netlib
    use powerspectra
    use moments
    use nonlinear
    use evolution

contains

    ! auxilary function for dverk
    ! this corresponds to the ODE y'(x)=exp(x)
    subroutine fun(n, x, y, yprime)
        implicit none
        integer n
        double precision x, y(n), yprime(n)
        yprime(1) = dexp(x)
    end subroutine fun


    subroutine test_ode
        double precision y(1),tol,xend,x,wspace(1,9),com(24)
        integer ind
        write(*,'(A)') "* Test ODE solver. Solution of y'(x)=exp(x), y(0)=1 at x=10?"
        ind=1 ! then no need to initialize com
        x=0.0d0
        xend=10.0d0
        y=(/1.d0/)
        tol=1.d-8
        ! x needs to be a variable, or else segfault
        call dverk(1,fun,x,y,xend,tol,ind,com,1,wspace)
        write(*,'(A,f15.8)') "  Number of evaluations of fun: ", com(24)
        write(*,'(A,f15.8)') "  y(10) = ", y(1)
        write(*,'(A,f15.8)') "  Should be:", dexp(10.d0)
    end subroutine test_ode
    

    subroutine test_progressbar
        integer i
        do i = 1, 10
            call progress(i) ! generate the progress bar.
            call sleep(1)
        end do
    end subroutine test_progressbar


    subroutine test_interpolation
        double precision x_(20), y_(20), objlin(2,20), objspl(5,20), x,y
        double precision, allocatable :: ps(:,:)
        integer i,j
        do i=1,20
            x_(i) = i
            y_(i) = f(x_(i))
        end do
        x = 11.1d0
        write(*,'(A,f15.8)') "* Test interpolation. f(x)=e^(x-10)², &
                &sampled at integers, at x =",x
        objlin = init_interpolation(x_, y_, cLINEAR)
        objspl = init_interpolation(x_, y_, cSPLINE)
        j = -1
        y = eval_interpolation(objlin, x, cLINEAR, j)
        write(*,"(A,f15.8)") "  Linear interpolation:", y
        j = -1
        y = eval_interpolation(objspl, x, cSPLINE, j)
        write(*,"(A,f15.8)") "  Spline interpolation:", y
        write(*,"(A,f15.8)") "  Should be: ", f(x)
        allocate(ps(2,1000))
        do i=1,size(ps,2)
            ps(1,i) = 1 + 19./size(ps,2)*i
            j = -1
            ps(2,i) = eval_interpolation(objspl, ps(1,i), cSPLINE, j)
        end do
        call export(ps, "output/interpol.dat")

        contains
            double precision function f(xx)
                double precision xx
                f = exp(-(xx-10.d0)**2)
                return
            end function f
    end subroutine test_interpolation


    subroutine test_imexport()
        double precision, allocatable, dimension(:,:) :: m
        call import(m, "scratch/z1_pk.dat",2)
        call export(m, "output/pk_out.dat", "Power spectrum test")
    end subroutine test_imexport


    subroutine test_spectra()
        double precision, dimension(4,1000) :: m
        double precision ps(3)
        integer i
        write(*,'(A)') "* Test power spectra."
        call powerspectrum(1,1.d-2,ps)        
        write(*,'(A,e15.8)') "  P11(1e-2) = ",ps(1)
        write(*,'(A,e15.8)') "  P12(1e-2) = ",ps(2)
        write(*,'(A,e15.8)') "  P22(1e-2) = ",ps(3)
        call powerspectrum(1,1.d+2,ps)        
        write(*,'(A,e15.8)') "  P11(1e+2) = ",ps(1)
        write(*,'(A,e15.8)') "  P12(1e+2) = ",ps(2)
        write(*,'(A,e15.8)') "  P11(1e+2) = ",ps(3)

        do i=1,size(m,2)
            m(1,i) = 1d1**(-7d0+i*11d0/1000)
            call powerspectrum(1, m(1,i), m(2:,i))
        end do
        call export(m, "output/test_ps.dat")
        write(*,"(A,3e15.8)") "  The slopes at k=0 are: ", ps_prefactor(1,:3)
    end subroutine test_spectra


    subroutine test_spline_integral()
        implicit none
        integer, parameter :: n = 100
        double precision  x_(n), y_(n), dummy
        double precision, allocatable :: intspl(:,:), intobj(:,:) 
        integer i,j 
        write(*,'(A)') "* Test spline integral. Integrate f(x)=exp(-x^2) from 1 to 10.&
                & Result should be: 1-1/10 = 9/10"
        do i=1,n
            x_(i) = 1d0+dble(i-1)*9d0/99
            y_(i) = exp(-x_(i)**2)
        end do
        allocate(intobj(5,size(x_)))
        intobj = init_interpolation(log(x_), y_, cSPLINE)
        allocate(intspl(6,size(x_)))
        call integrate_splines(intobj, intspl, -1)
        do i=10,10
            j = -1
            dummy = eval_int_splines(intspl, dble(i), j)
            write(*,'(A,i2.2,A,e15.8)') "  F(",i,") = ",dummy
        end do
    end subroutine test_spline_integral


    subroutine test_moments()
        double precision, allocatable ::  moms(:,:), intobj(:,:)
        double precision  temp(2,1000)
        character*100 filename
        integer i,j,l

        write(*,'(A)') "* Test moments. Integrating them and dumping them to&
                & output/integrated_moments??.dat."
        do i=1,13
            write(filename,'(A,i2.2,A)') "output/integrated_moments",i,".dat"
            call get_integrated_moment(i,1,1d0,moms)
            call export(moms,filename(:))
            write(filename,'(A,i2.2,A)') "output/momentum",i,".dat"
            dummy = get_moment(i,1,1d0,moms)
            allocate(intobj(2,size(moms,2)))
            intobj = init_interpolation(log(moms(1,:)),log(moms(2,:)), cLINEAR)
            do j=1,1000
                temp(1,j) = exp(log(moms(1,1)) + (log(moms(1,size(moms,2))) - &
                        log(moms(1,1)))/1000*j)
                l = -1
                temp(2,j) = exp(eval_interpolation(intobj, log(temp(1,j)), cLINEAR, l))
            end do 
            deallocate(intobj)
            call export(temp,filename(:))
        end do
    end subroutine test_moments


    subroutine test_pointerarray()
        implicit none
        type(tPointer_array) ptrarr(2)
        double precision, target, allocatable ::  a(:), b(:)

        write(*,'(A)') "* Test array of pointers. &
            &Output should be 2, 2, 1, 1, 1."

        allocate(a(2))
        a = (/2.d0, 2.d0/)
        allocate(ptrarr(1)%ptr(size(a)))
        ptrarr(1)%ptr = a
        deallocate(a)
        allocate(a(3))
        a = (/1.d0, 1.d0, 1.d0/)
        allocate(ptrarr(2)%ptr(size(a)))
        ptrarr(2)%ptr = a
        deallocate(a)

        write(*,*) ptrarr(1)%ptr
        write(*,*) ptrarr(2)%ptr

    end subroutine test_pointerarray


    subroutine test_parallel()
        implicit none
        integer i
!$OMP PARALLEL DO
        do i=1,2
            call sleep(10)
        end do
!$OMP END PARALLEL DO
    end subroutine test_parallel


    subroutine test_A()
        implicit none
        double precision ak(14), t1, t2, k
        double precision, allocatable :: dummy(:,:)
        integer i,j

        write(*,'(A)') "* Test A(k)."
        ! call dump_A(1) TODO
        ! do i=size(klist)-5,size(klist)
            ! k = klist(i)
            ! write(*,'(A,e15.8,A,e15.8,f8.3)') "  A(",k,") = ",A(1,1,k), t2-t1
        ! end do
    end subroutine test_A


    subroutine test_Kdq()
        implicit none
        double precision ak(14), k
        integer j

        write(*,'(A)') "* Test K(q)."
        
        ! do j=0,11
        !     k = 10.d0**(j/2.d0-3)
        !     ! ak = init_A(1,k,dump=.true.)
        ! end do
        ! call dump_Kdq(1) TODO
        ! k = klist(size(klist))*1.1d0
        ! ak = init_A(1,k,dump=.true.) TODO fix it
    end subroutine test_Kdq


    subroutine test_time_evolution()
        implicit none
        write(*,'(A)') "* Test time evolution."
#if defined(__INTEL_COMPILER)
        ! call time_evolution_oneloop(-log(aini), .false., .false.) !linear
        ! call time_evolution_oneloop(-log(aini), .true., .false.) !non-linear
        ! call time_evolution_oneloop(-log(aini), .false., .true.) !linear, smooth
#endif
        ! call time_evolution_full(-log(aini), .true.) !non-linear
        ! call time_evolution_nloop(-log(aini), 3, .true.) !non-linear
        
    end subroutine test_time_evolution


end module test_mod
	


program test

    use test_mod

    call load_PS_from_file
    call init_moments(1)
    call init_integrated_moment(1)

    ! call test_parallel
    ! call test_progressbar
    ! call test_interpolation
    ! call test_ode
    ! call test_imexport
    ! call test_spectra
    ! call test_background
    ! call test_antiderivative
    call test_spline_integral
    ! call test_moments
    ! call test_integrator
    ! call test_pointerarray
    ! call test_A
    ! call test_Kdq
    call test_time_evolution


end program test

