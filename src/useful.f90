module useful

   use parallel_mod
   use precision_mod
   use constants, only: tiny_number, two, one, pi, three, half
   
   implicit none

   private

   public :: abortRun, logWrite, formatTime, polynomial_interpolation, & 
   &         l_correct_step, round_off, cc2real, cc22real, getMSD2,    &
   &         gausslike_compact_edge, gausslike_compact_center,         &
   &         gausslike_compact_middle

contains

   subroutine abortRun(message)
      !
      ! This routine properly terminates a run
      !

      !-- Input variable
      character(len=*), intent(in) :: message

      !-- Local variables:
      integer :: code

      code = 32

      if ( rank == 0 ) then
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,*) '! Something went wrong, pizza will stop running now'
         write(*,*) '! See below the error message:'
         write(*,*)
         write(*,*) message
         write(*,*)
      end if

      call MPI_Abort(MPI_COMM_WORLD, code, ierr)

   end subroutine abortRun
!----------------------------------------------------------------------------
   subroutine logWrite(message, n_log_file)

      !-- Input variable:
      character(len=*), intent(in) :: message
      integer,          intent(in) :: n_log_file

       if ( rank == 0 ) then
          write(n_log_file,*) trim(message)
          write(*,*)          trim(message)
       end if

   end subroutine logWrite
!----------------------------------------------------------------------------
   subroutine formatTime(n_out, message, time_in_sec)

      !-- Input variable
      integer,          intent(in) :: n_out
      character(len=*), intent(in) :: message
      real(cp),         intent(in) :: time_in_sec

      !-- Local variable
      integer :: n_time(5)
      real(cp) :: remaining_time

      n_time(1) = int(time_in_sec/86400_cp)
      remaining_time = time_in_sec-n_time(1)*86400_cp
      n_time(2) = int(remaining_time/3600_cp)
      remaining_time = remaining_time-n_time(2)*3600_cp
      n_time(3) = int(remaining_time/60_cp)
      remaining_time = remaining_time-n_time(3)*60_cp
      n_time(4) = int(remaining_time)
      remaining_time = remaining_time-n_time(4)
      n_time(5) = int(remaining_time/1e-3_cp)

      if ( time_in_sec > 0.1_cp .and. time_in_sec < 60.0_cp ) then
         write(n_out,'(1x,A,A,f6.3,A/,1x)') message, ' ', time_in_sec, " seconds"
      else if ( time_in_sec >= 60.0_cp .and. time_in_sec < 3.6e3_cp ) then
         write(n_out,'(1x,A,A,I2,A,I2,A,I3,A/,1x)') message, ' ', n_time(3), &
         &     " m ", n_time(4), " s ", n_time(5), " ms"
      else if ( time_in_sec >= 3.6e3_cp .and. time_in_sec < 8.64e4_cp ) then
         write(n_out,'(1x,A,A,I2,A,I2,A,I3,A,I3,A/,1x)') message, ' ',    &
         &     n_time(2), " h ", n_time(3), " m ", n_time(4), " s ",      &
         &     n_time(5), " ms"
      else if ( time_in_sec >= 8.64e4_cp ) then
         write(n_out,'(1x,A,A,I2,A,I2,A,I3,A,I3,A,I3,A/,1x)') message, ' ', &
         &     n_time(1), " d ", n_time(2), " h ", n_time(3), " m ",        &
         &     n_time(4), " s ", n_time(5), " ms"
      else
         write(n_out,'(1x,A,A,es10.3,A/,1x)') message, ' ', time_in_sec, " seconds"
      end if

   end subroutine formatTime
!----------------------------------------------------------------------------
   subroutine polynomial_interpolation(xold, yold, xnew ,ynew)

      !-- Input variables
      real(cp),    intent(in) :: xold(4)
      complex(cp), intent(in) :: yold(4)
      real(cp),    intent(in) :: xnew

      !-- Output variables
      complex(cp), intent(out) :: ynew

      !-- Local variables
      real(cp) :: yold_real(4), yold_imag(4)
      real(cp) :: ynew_real, ynew_imag

      yold_real= real(yold)
      yold_imag=aimag(yold)

      call polynomial_interpolation_real(xold, yold_real, xnew, ynew_real)
      call polynomial_interpolation_real(xold, yold_imag, xnew, ynew_imag)

      ynew = cmplx(ynew_real, ynew_imag, kind=cp)

   end subroutine polynomial_interpolation
!----------------------------------------------------------------------------
   subroutine polynomial_interpolation_real(xold,yold,xnew,ynew)

      !-- Input variables:
      real(cp), intent(in) :: xold(:)
      real(cp), intent(in) :: yold(:)
      real(cp), intent(in) :: xnew

      !-- Output variables:
      real(cp), intent(out) :: ynew

      !-- Local variables:
      integer :: n_stencil
      integer :: n_st, n_st_out, n_s
      real(cp) :: diff, diff_tmp, dy
      real(cp) :: ho, hp, den, work_diff
      real(cp), allocatable :: work1(:), work2(:)

      n_stencil=size(xold)
      allocate( work1(n_stencil), work2(n_stencil) )

      n_s=1
      diff=abs(xnew-xold(1))
      do n_st=1,n_stencil
         diff_tmp=abs(xnew-xold(n_st))
         if ( diff_tmp < diff ) then
            n_s =n_st
            diff=diff_tmp
         end if
         work1(n_st)=yold(n_st)
         work2(n_st)=yold(n_st)
      end do
      ynew=yold(n_s)

      n_s=n_s-1
      do n_st_out=1,n_stencil-1
         do n_st=1,n_stencil-n_st_out
            ho       =xold(n_st)-xnew
            hp       =xold(n_st+n_st_out)-xnew
            work_diff=work1(n_st+1)-work2(n_st)
            den=ho-hp
            if ( den == 0.0_cp ) call abortRun('Stop in polynomial interpolation')
            den        =work_diff/den
            work2(n_st)=hp*den
            work1(n_st)=ho*den
         end do
         if ( 2*n_s < n_stencil-n_st_out )then
            dy=work1(n_s+1)
         else
            dy=work2(n_s)
            n_s=n_s-1
         end if
         ynew=ynew+dy
      end do

      deallocate( work1, work2 )

   end subroutine polynomial_interpolation_real
!----------------------------------------------------------------------------
   logical function l_correct_step(n,n_max,n_step,n_intervals)
      !
      ! Suppose we have a (loop) maximum of n_max steps!
      ! If n_intervals times in these steps a certain action should be carried out
      ! this can be invoked by l_correct_step=true if on input n_intervals>0
      ! and n_step=0.
      ! Alternatively the action can be invoked every n_step steps if
      ! on input n_intervals=0 and n_step>0.
      ! In both cases l_correct_step=true for n=n_max.
      !

      !-- Input variables:
      integer,  intent(in) :: n            ! current step
      integer,  intent(in) :: n_max        ! max number of steps
      integer,  intent(in) :: n_step       ! action interval
      integer,  intent(in) :: n_intervals ! number of actions

      !-- Local variables:
      integer :: n_offset     ! offset with no action
      integer :: n_steps      ! local step width


      if ( n_step /= 0 .and. n_intervals /= 0 ) then
         write(*,*) '! Error message from function l_correct_step:'
         write(*,*) '! Either n_step or n_interval have to be zero!'
         call abortRun('Stop run in l_correct_step')
      end if

      l_correct_step=.false.

      if ( n_intervals /= 0 ) then
         n_steps=n_max/n_intervals
         n_offset=n_max-n_steps*n_intervals

         if ( n > n_offset .and. mod(n-n_offset,n_steps) == 0 ) l_correct_step= .true.
      else if ( n_step /= 0 ) then
         if ( n == n_max .or. mod(n,n_step) == 0 ) l_correct_step= .true.
      end if

   end function l_correct_step
!----------------------------------------------------------------------------
   real(cp) function round_off(param)

      !-- Input variable
      real(cp), intent(in) :: param

      if ( abs(param) < tiny_number ) then
         round_off = 0.0_cp
      else
         round_off = param
      end if

   end function round_off
!----------------------------------------------------------------------------
   real(cp)  function cc2real(c,m)

      !-- Input variables:
      complex(cp), intent(in) :: c
      integer,     intent(in) :: m

      if ( m == 0 ) then
         cc2real=real(c)*real(c)
      else
         cc2real=two*(real(c)*real(c)+aimag(c)*aimag(c))
      end if

   end function cc2real
!----------------------------------------------------------------------------
   real(cp) function cc22real(c1,c2,m)

      !-- Input variables:
      complex(cp), intent(in) :: c1,c2
      integer,     intent(in) :: m

      if ( m == 0 ) then
         cc22real=real(c1)*real(c2)
      else
         cc22real=two*(real(c1)*real(c2)+aimag(c1)*aimag(c2))
      end if

   end function cc22real
!----------------------------------------------------------------------------
   subroutine getMSD2(mean,SD,x,n,dt,totalTime)
      ! This subroutine computes the mean and standard deviation according 
      ! to a method introduced by Donald Knuth (1962). I rederived his formulas
      ! for a variable time step. On output SD still needs to be normalized with 
      ! the totalTime and then you have to take the square root!!
      ! The input integer counts the number of calls. For n=1 initialisation
      ! is necessary.

      !-- Input variables:
      real(cp), intent(in) :: x         ! quantity to be averaged
      real(cp), intent(in) :: dt        ! time since last averaging step
      real(cp), intent(in) :: totalTime ! total averaging time up to now
      integer,  intent(in) :: n         ! number of calls( only n=1 needed)

      !-- Output variables:
      real(cp), intent(out) :: mean     ! Time-average
      real(cp), intent(out) :: SD       ! Standard-deviation

      !-- Local variable:
      real(cp) :: delta

      if ( n == 1) then
         mean=x
         sd=0.0_cp
      else
         delta=x-mean
         mean=mean+delta*dt/totalTime
         SD=SD+dt*delta*(x-mean)
      end if

   end subroutine getMSD2
!----------------------------------------------------------------------------
   real(cp) function gausslike_compact_center(c,L) result(gasp3)
      !
      ! This defines the central point of a compact-support Gaussian-like
      ! profile. This is adapted from a Gaspari & Cohn, QJRMS, 1999, Eq. 4.7c
      !

      !-- Input variables
      real(cp), intent(in) :: c
      real(cp), intent(in) :: L

      gasp3 = pi*L*L*L*(one-exp(-two*c/L))-two*pi*c*L*(c+L)*exp(-two*c/L)

   end function gausslike_compact_center
!----------------------------------------------------------------------------
   real(cp) function gausslike_compact_middle(r,c,L) result(gasp1)
      !
      ! This defines the first part of a compact-support Gaussian-like
      ! profile. This is adapted from a Gaspari & Cohn, QJRMS, 1999, Eq. 4.7b
      !

      !-- Input variables
      real(cp), intent(in) :: r
      real(cp), intent(in) :: c
      real(cp), intent(in) :: L

      gasp1 = pi/three*L*r*(r+three*L)*exp(-r/L)+(two*pi*L*L*(c+L)**2)/r *   &
      &       exp(-two*c/L)*(one-(one-r/(c+L))*exp(r/L))+pi*L*L*L*(exp(-r/L)-&
      &       exp((r-two*c)/L))


   end function gausslike_compact_middle
!----------------------------------------------------------------------------
   real(cp) function gausslike_compact_edge(r,c,L) result(gasp2)
      !
      ! This defines the second part of a compact-support Gaussian-like
      ! profile. This is adapted from a Gaspari & Cohn, QJRMS, 1999, Eq. 4.7a
      !

      !-- Input variables
      real(cp), intent(in) :: r
      real(cp), intent(in) :: c
      real(cp), intent(in) :: L

      gasp2 = two*pi*L/r*exp(-r/L)*(half*(r*(r+L)*(two*c-r))+one/three*( &
      &       (r-c)**3-c**3)-L*(c+L)*(r-c+L)+L*(c+L)**2*exp((r-two*c)/L))

   end function gausslike_compact_edge
!----------------------------------------------------------------------------
end module useful
