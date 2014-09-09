! Contains general, non-specific helper functions 
! @author Adrian Vollmer
! @date 24.07.2012

! A collection of tools. 
!
! The goal here is to keep these independent from the main program such that
! they could easily be used in other projects. No physics shall appear in this
! module. 


module tools

    use options
#ifdef __INTEL_COMPILER
    use ifport, only : system
#endif

    implicit none

        ! interface to the F77 routines from netlib.org
    interface 
        subroutine spline(n, x_, y_, b_, c_, d_)
            integer n
            double precision, dimension(n) :: x_, y_, b_, c_, d_
        end subroutine spline
    end interface

        ! some constants
    double precision, parameter :: Pi=3.141592653589793d0
    double precision, parameter :: PiHalf=Pi/2
    double precision, parameter :: TwoPiCubed = (2*pi)**3
        ! these are for interpolating functions
    integer, parameter :: cLINEAR = 0
    integer, parameter :: cSPLINE = 1


contains

    function dump_dir()
        implicit none
        character(19) dump_dir
        integer res
        write(dump_dir,'(A,i5.5,A)') "/tmp/trgfast-", oDUMP_DIR, "/"
        res = system('mkdir -p '//dump_dir) ! create directory
        ! TODO make this work for windows?
        return
    end function dump_dir


    function is_nan(x)
        double precision, intent(in) :: x
        logical is_nan
        if (x.ne.x) then
            is_nan = .true.
        else
            is_nan = .false.
        end if
        return
    end function is_nan


    double precision function mean(x)
        implicit none
        double precision, intent(in) :: x(:)
        mean = sum(x) / size(x)
        return
    end function mean


    double precision function stddev(x)
        implicit none
        double precision, intent(in) :: x(:)
        double precision m
        m = mean(x)
        stddev = sqrt( sum( (x-m)**2 )  / (size(x)) )
        return
    end function stddev

! I got this from Dr. Fortran:
! http://software.intel.com/en-us/forums/showthread.php?t=73766
! j needs to be an integer between 1 and 10. msg is an optional string.

    subroutine progress(j, msg)  
      implicit none  
      integer(kind=4)::j,k  
      character, optional, intent(in) :: msg*(*)
      character(len=18)::bar
      if (oVERBOSITY>=1) then
          bar = "???% |          | "  
          write(unit=bar(1:3),fmt="(i3)") 10*j  
          do k=1, j  
            bar(6+k:6+k)="*"  
          enddo  
          ! print the progress bar.  
          ! write(unit=6,fmt="(a1,a1,a17)") '+',char(13), bar  
          if (present(msg)) then
              write(unit=6,fmt="(a1,a20,a,$)") char(13), bar, msg
          else
              write(unit=6,fmt="(a1,a20,$)") char(13), bar  
          end if
          return  
      end if
    end subroutine progress


    function aspace(a, b, stepnumber)
        implicit none
        double precision, intent(in) :: a,b
        integer, intent(in) :: stepnumber
        double precision, allocatable :: aspace(:)
        integer i
        allocate(aspace(stepnumber))
        aspace = (/( a+i*(b-a)/(stepnumber-1), i=0,stepnumber-1)/)
        return
    end function aspace

    function alogspace(a, b, stepnumber)
        implicit none
        double precision, intent(in) :: a,b
        integer, intent(in) :: stepnumber
        double precision, allocatable :: alogspace(:)
        integer i
        allocate(alogspace(stepnumber))
        alogspace = exp((/( log(a)+i*(log(b)-log(a))/(stepnumber-1), i=0,stepnumber-1)/))
        return
    end function alogspace
    

    subroutine message(msgtype, msg, subname)
        implicit none
        character, intent(in) :: msgtype, msg*(*), subname*(*)
        if (msgtype(:)=="c") then
            write(*,'(A,A,100A)') "*** Critical error (",subname,") ",msg
            stop
        else  if (msgtype=="e") then
            write(*,'(A,A,100A)') "E (",subname,") ",msg
        else if (msgtype=="w" .and.oVERBOSITY>=1) then
            write(*,'(A,A,100A)') "W (",subname,") ",msg
        else if (msgtype=="i" .and.oVERBOSITY>=2) then
            write(*,'(A,A,100A)') "I (",subname,") ",msg
        else if (msgtype=="d" .and.oVERBOSITY>=3) then
            write(*,'(A,A,100A)') "D (",subname,") ",msg
        end if
    end subroutine message


! This subroutine writes a matrix in table form into a file. 
!
! You can leave a comment. A time stamp is saved in any case, unless bare=.true.
! The matrix must be of the form (col,row). It's compatible with import().
!
! * m the matrix you want to export
! * filename the filename
! * comment optional comment that is printed at the top of the file
! * bare unless bare is true, a comment containing time and date is printed at
!   the top of the file.

    subroutine export(m, filename, comment, bare)
        implicit none
        double precision, intent(in) :: m(:,:)
        character, intent(in) :: filename*(*)
        character, intent(in), optional ::  comment*(*)
        logical, intent(in), optional :: bare
        logical isbare
        integer i,j,s,ios, date_time(8)
        character*10 b(3)

        call date_and_time(b(1), b(2), b(3), date_time)

        isbare=.false.
        if (present(bare)) then 
            isbare=bare
        end if

        open(unit=11, file=filename, status='unknown', iostat=ios)

        if (ios == 0) then
            if (.not.isbare) then
                write(11, 900)date_time(3), date_time(2), date_time(1), &
                    date_time(5), date_time(6), date_time(7)
                write(11,'(A, A)') "# ", comment
900             format('# ', i2.2, '.', i2.2, '.', i4.4, " ", i2.2, ':', i2.2, ':', i2.2)
            end if
            s = size(m,2)
            do i=1,size(m,1)
               ! the 100 in the format statement below limits the number of possible
               ! columns.
               write(11,'(100E20.12)')(m(i,j), j=1,s)
            end do
            close(11)
        else
            if (size(m,2)==0) then
                call message("w", &
                    msg="No lines written to "//filename//".",&
                    subname="export")
            end if
        end if
        
    end subroutine export


! This subroutine imports a table of float numbers from a file that can contain
! comments marked by a '#'. 
! 
! The number of columns must be given. The result is returned in m, and each
! element can then be accessed with m(col,row). It's compatible with export().
!
! * m allocatable array to be filled with the matrix you want to import
! * filename the filename
! * columns the number of columns you want to import

    function import(filename, columns) result(m)
        implicit none
        character, intent(in) :: filename*(*)
        character(100) dummy
        integer, intent(in), optional :: columns
        integer line_count, stat, i,j
        double precision, allocatable :: mtemp(:,:), m(:,:)

        open(unit=12, file=filename, status='old', action='read', iostat=stat)
        if (stat /= 0) then
            call message("c", msg="File not found: "//filename(:), subname="import")
        end if
        open(unit=13, status = 'scratch')

        line_count=0
        do
            read(12,'(a)',iostat=stat) dummy
            dummy = adjustl(trim(dummy))
            if (stat /= 0) exit
            if (dummy(1:1).ne.'#') then ! If it's not a comment...
                write (13,'(a)') dummy
                line_count = line_count + 1
            end if
        end do

        ! Read the non-comment lines into m
        ! TODO detect no of columns and discard the argument "columns"
        allocate(mtemp(columns,line_count))
        allocate(m(line_count,columns))
        rewind(13)
        read(13,*) mtemp
        close(13)
        close(12)
            ! transpose result
        do i=1,size(m,1)
            do j=1,size(m,2)
                m(i,j) = mtemp(j,i)
            end do
        end do
        deallocate(mtemp)
    end function import

! A linear search algorithm. It returns the index of the interval that contains
! the search value.

! Riccardo was using locate and hunt from NumRep, which has unacceptable
! licensing conditions

    integer function lsearch(a_, x)
        implicit none
        double precision, intent(in) :: x, a_(:)
        integer i

        if (x<a_(1)) then
            lsearch = 0
            return
        end if

        if (x>=a_(size(a_)-1)) then
            lsearch = size(a_)-1
            return
        end if

        do i=1,size(a_)-1
            if (x<=a_(i+1)) exit
        end do
        lsearch = i
        
        return
    end function lsearch


! A recursive binary search algorithm. It returns the index of the interval
! that contains the search value. If the array is small, linear search is used.
!
! Wikipedia on binary search: 
!
! >Binary search can interact poorly with the memory hierarchy (i.e. caching),
! >because of its random-access nature.  For in-memory searching, if the span to
! >be searched is small, a linear search may have superior performance simply
! >because it exhibits better locality of reference. For external searching,
! >care must be taken or each of the first several probes will lead to a disk
! >seek. A common method is to abandon binary searching for linear searching as
! >soon as the size of the remaining span falls below a small value such as 8 or
! >16 or even more in recent computers.  The exact value depends entirely on the
! >machine running the algorithm.
!  
!
! * a_ allocated array that is to be searched 
! * x search value

    recursive integer function bsearch(a_, x) result (n)
        implicit none
        double precision, intent(in) :: x, a_(:)
        integer length, mid

        length = size(a_)
        mid = length/2 + 1

            ! switch to linear search for small arrays, because they can be
            ! cached entirely, which is faster than seemingly random accessing
        if (length <= 16) then
            n = lsearch(a_, x)
            return
        end if

        if (x<a_(mid)) then
            n = bsearch(a_(:mid),x)
            return
        else
            n = bsearch(a_(mid:),x) + mid - 1
            return
        end if
    end function bsearch


    subroutine sort(x)
        implicit none
        double precision, intent(inout) :: x(:,:)
        integer i,j,n
        double precision tmp(2)
        n = size(x,1)
        do i=1,n-1
            do j=i+1,n
                if (x(i,1) > x(j,1)) then
                   tmp = x(i,:)
                   x(i,:) = x(j,:)
                   x(j,:) = tmp
                end if
            end do
        end do
    end subroutine sort

        ! removes duplicate elements of a sorted array and returns the number of
        ! deletions
    function unique(x)
        implicit none
        double precision, intent(inout) :: x(:,:)
        integer unique, i, j
        unique = 0
        do i=size(x,1)-1,1,-1
            if (x(i,1) .eq. x(i+1,1)) then
                unique = unique + 1
                do j=i+1,size(x,1)
                    x(j-1,:) = x(j,:)
                end do
            end if
        end do
        return
    end function unique

    
! This function returns an interpolation object which can then be later used for
! evaluation. 
! 
! * x_ array of x values
! * y_ array of y values
! * int_type must take the value 0/cLINEAR (for linear interpolation)
!   or 1/cSPLINE (for natural cubic splines)

    function init_interpolation(x_, y_, int_type)
        implicit none
        double precision, allocatable :: init_interpolation(:,:)
        double precision, intent(in) :: x_(:), y_(:)
        integer, intent(in) :: int_type
        double precision, allocatable :: b_(:), c_(:), d_(:)
        integer n, i

        ! check if data is sorted
        do i = 2,size(x_)
            if (x_(i-1)>=x_(i)) then
                call message("c",msg="Data not sorted in ascending order.",&
                    subname="init_interpolation")
            end if
        end do

        n = size(x_)
        if (int_type==cSPLINE) then
            allocate(b_(n))
            allocate(c_(n))
            allocate(d_(n))
            allocate(init_interpolation(5,n))
                ! call netlib function
            call spline(n, x_, y_, b_, c_, d_)
            init_interpolation(1,:) = x_(:)
            init_interpolation(2,:) = y_(:)
            init_interpolation(3,:) = b_(:)
            init_interpolation(4,:) = c_(:)
            init_interpolation(5,:) = d_(:)
                ! check if at least this one coefficient is not NaN
                ! Hopefully this is compiler independent
            if (is_nan(init_interpolation(5,1))) then
                call message("c", msg="Interpolation failed.", &
                        subname="init_interpolation")
            end if
            return
        else if (int_type==cLINEAR) then
            allocate(init_interpolation(2,n))
            init_interpolation(1,:) = x_(:)
            init_interpolation(2,:) = y_(:)
            return
        end if
        call message("e", msg="Unknown interpolation type.",&
                subname="init_interpolation")
        return
    end function init_interpolation



! This function returns the value of an interpolated function at x when provided
! with an interpolation object computed with init_interpolation() with matching
! int_type. 
!
! Sometimes you already know in which specific interval x lies. You can speed up
! this function by suggesting an interval. If the suggested interval is
! negative, a binary search algorithm is used to find the interval which
! contains x. After returning the result, indx will contain the number of the
! interval in which x was found such that it can possibly be reused.
!
! * intobject a 2D array computed with init_interpolation
! * x value at which to evaluate the interpolated function
! * int_type specifies whether the function has been interpolated in
!   linear space or log space
! * indx suggested index

    function eval_interpolation(intobject, x, int_type, indx)
        implicit none
        double precision eval_interpolation
        double precision, intent(in) :: x, intobject(:,:)
        integer, intent(in) :: int_type
        integer, intent(inout) :: indx
        double precision xmxi
        integer n
        integer i

        eval_interpolation = 0d0

        if (indx<0) then ! don't reuse index
            n = size(intobject,2)
                ! find the x interval 
            i = bsearch(intobject(1,:),x)
                ! range check
            if (i>=n) i=n-1
            if (i<=0) i = 1
            indx = i
        else
            i = indx
        end if

        xmxi = x-intobject(1,i)
        if (int_type==cSPLINE) then
            eval_interpolation = intobject(2,i) + xmxi*(intobject(3,i) + &
                    xmxi*(intobject(4,i) + xmxi*intobject(5,i)))
            return
        else if (int_type==cLINEAR) then
            eval_interpolation = intobject(2,i) + &
                        xmxi * (intobject(2,i+1) - intobject(2,i)) &
                        / (intobject(1,i+1) - intobject(1,i))
            return
        end if
        call message("e", msg="Unknown interpolation type.", &
                subname="eval_interpolation")
        return
    end function eval_interpolation
      

! Piecewise integration of interpolating splines. 
!
! Input needs to be a spline
! object computed via init_interpolation, output is a modified spline object, m
! the moment. The result is the integral of the spline in logspace multiplied
! with x^m. It can be evaluated with eval_int_splines.
!
! * input the spline object
! * output the resulting integrated spline
! * m the moment

    subroutine integrate_splines(input, output, m)
        implicit none
        double precision, intent(in) :: input(:,:)
        double precision, intent(out) :: output(6,size(input,2))
        integer, intent(in) :: m
        integer i,mp1,d1,d2
        double precision, dimension(size(input,2)) :: x,y,b,c,d,expx

            ! some abbreviations 
        x = input(1,:)
        y = input(2,:)
        b = input(3,:)
        c = input(4,:)
        d = input(5,:)
        mp1 = m+1

        output(6,:) = (/(0d0, i=1,size(x))/)
        output(6,size(x)) = dble(mp1)
        output(1,:) = x

        if (m==-1) then
            output(2,:) = b/2
            output(3,:) = c/3
            output(4,:) = d/4
            output(5,:) = y-b*x
        else
            output(2,:) = -(mp1**2*b) + 2*mp1*c - 6*d + mp1**3*y
            output(3,:) = mp1*(mp1**2*b - 2*mp1*c + 6*d)
            output(4,:) = mp1**2*(mp1*c - 3*d)
            output(5,:) = mp1**3*d
        end if
        
            ! Integration constant TODO research it
        expx = exp(x)
        if (m>0) then 
            d1 = 1
            output(6,1) = -eval_int_splines(output,expx(1),x(1),d1)
            do i=2,size(x)-1
                d1 = i
                d2 = i-1
                output(6,i) = -eval_int_splines(output,expx(i),x(i), d1) &
                            +eval_int_splines(output,expx(i),x(i), d2)
            end do
        else
            i = size(x)
            d1 = i-1
            output(6,i-1) = -eval_int_splines(output,expx(i),x(i), d1)
            do i=size(x)-2,1,-1
                d1 = i
                d2 = i+1
                output(6,i) = -eval_int_splines(output,expx(i+1),x(i+1),d1) &
                            +eval_int_splines(output,expx(i+1),x(i+1),d2)
            end do
        end if
    end subroutine integrate_splines


! Evaluates integrated splines
!
! Integrated splines aren't technically splines anymore, in the sense that they
! aren't piecewise defined polynomials, but they are still piecewise defined
! functions - just more general. Usage is similar to eval_interpolation.
!
! * intobj a 2D array computed with integrate_splines
! * x value at which to evaluate the interpolated function
! * indx suggested index

    double precision function eval_int_splines(intobj, xin, logxin, indx)
        implicit none
        double precision, intent(in) :: xin, logxin, intobj(:,:)
        integer, intent(inout) :: indx
        double precision xminusxi,  M, x, logx
        integer i

        x = xin
        logx = logxin

        if (indx>0) then
            i = indx
        else
            i = bsearch(intobj(1,:), logx)
        end if

        if (is_nan(logx)) then
            call message("c", msg="NaN in argument encountered", subname="eval_int_splines")
        end if
            
        if (i<=0) then  ! out of range, extrapolation used 
            call message("d", msg="Index out of range, extrapolating", &
                subname="eval_int_splines") 
            i=1
            logx=intobj(1,1)
            x=exp(logx)
        else if (i>=size(intobj,2)) then
            call message("d", msg="Index out of range, extrapolating", &
                subname="eval_int_splines") 
            i=size(intobj,2)-1
            logx=intobj(1,size(intobj,2))
            x=exp(logx)
        end if

        indx = i
        xminusxi = logx-intobj(1,i)
        M = intobj(6,size(intobj,2))
            ! See dissertation for this equation
        if (M/=0d0) then
            eval_int_splines = x**M/M**4 * (&
                    intobj(2,i) + xminusxi * (intobj(3,i) + &
                    xminusxi * (intobj(4,i) + xminusxi * intobj(5,i)))) &
                    + intobj(6,i)
        else
            eval_int_splines = logx * (intobj(5,i) + logx * intobj(2,i)) &
                    + xminusxi**3 * (intobj(3,i) + xminusxi * intobj(4,i)) + intobj(6,i)
        end if

        return

    end function eval_int_splines

end module tools
