module timers_mod
   !
   ! This module is used to set several counters that are use to measure
   ! the time passed in several critical parts of the code
   !

   use precision_mod
   use parallel_mod
   use iso_fortran_env, only: output_unit
   use communications, only: my_reduce_mean
   use useful, only: formatTime

   implicit none

   private

   type, public :: timers_type
      integer :: n_r_loops
      integer :: n_m_loops
      integer :: n_m_loops_mat
      integer :: n_mpi_comms
      integer :: n_io_calls
      integer :: n_fft_calls
      integer :: n_solve_calls
      integer :: n_dct_calls
      integer :: n_lu_calls
      real(cp) :: r_loop      ! Timer of the radial loop (nonlinear terms)
      real(cp) :: m_loop      ! Timer of the m loop (linear solve/time advance)
      real(cp) :: io          ! Timer for I/O
      real(cp) :: mpi_comms   ! Timer for the main all_to_all communications
      real(cp) :: m_loop_mat  ! Timer for the m loop when matrix LU fac is needed
      real(cp) :: fft         ! Timer for FFTs
      real(cp) :: lu          ! Timer for one LU factorisation
      real(cp) :: solve       ! Timer for one linear solve
      real(cp) :: dct         ! Timer for DCTs
      real(cp) :: tot         ! Total runtime
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: write_log
   end type timers_type

   type, public :: timer_type
      integer :: n_counts
      real(cp) :: tStart
      real(cp) :: tTot
   contains
      procedure :: initialize => initialize_timer
      procedure :: finalize => finalize_timer
      procedure :: start_count
      procedure :: stop_count
   end type timer_type

contains

   subroutine initialize(this)
      !
      ! Initialize the counters with zeroes
      !

      class(timers_type) :: this

      this%n_r_loops     = 0
      this%n_m_loops     = 0
      this%n_m_loops_mat = 0
      this%n_mpi_comms   = 0
      this%n_io_calls    = 0
      this%n_fft_calls   = 0
      this%n_solve_calls = 0
      this%n_dct_calls   = 0
      this%n_lu_calls    = 0
      this%r_loop     = 0.0_cp
      this%m_loop     = 0.0_cp
      this%io         = 0.0_cp
      this%mpi_comms  = 0.0_cp
      this%m_loop_mat = 0.0_cp
      this%fft        = 0.0_cp
      this%lu         = 0.0_cp
      this%solve      = 0.0_cp
      this%dct        = 0.0_cp
      this%tot        = 0.0_cp

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this, n_steps)
      !
      ! Calculate wall time for different part of the code
      ! and average over the different ranks
      !
      class(timers_type) :: this

      !-- Input variables:
      integer, intent(in) :: n_steps

      this%io        = this%io/this%n_io_calls
      call my_reduce_mean(this%io, 0)
      this%r_loop    = this%r_loop/this%n_r_loops
      call my_reduce_mean(this%r_loop, 0)
      if ( this%n_m_loops /= 0 ) then
         this%m_loop    = this%m_loop/this%n_m_loops
         call my_reduce_mean(this%m_loop, 0)
      end if
      if ( this%n_m_loops_mat /= 0 ) then
         this%m_loop_mat=this%m_loop_mat/this%n_m_loops_mat
         call my_reduce_mean(this%m_loop_mat, 0)
      end if
      if ( n_steps /= 0 ) then
         this%tot=this%tot/n_steps
         call my_reduce_mean(this%tot, 0)
      end if
      this%mpi_comms = this%mpi_comms/this%n_mpi_comms
      call my_reduce_mean(this%mpi_comms, 0)
      if ( this%n_fft_calls /= 0 ) then
         this%fft = this%fft/this%n_fft_calls
      end if
      if ( this%n_lu_calls /= 0 ) then
         call my_reduce_mean(this%fft, 0)
         this%lu = this%lu/this%n_lu_calls
      end if
      if ( this%n_solve_calls /= 0 ) then
         call my_reduce_mean(this%lu, 0)
         this%solve = this%solve/this%n_solve_calls
      end if
      call my_reduce_mean(this%solve, 0)
      if ( this%n_dct_calls /= 0 ) then
         this%dct = this%dct/this%n_dct_calls
         call my_reduce_mean(this%dct, 0)
      end if

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine write_log(this, n)
      !
      ! Writes the timing into a log file
      !
      class(timers_type) :: this

      !-- Input variables
      integer :: n

      if ( rank == 0 ) then
         call formatTime(n, &
         &    '! Mean wall time for radial loop            :',this%r_loop)
         call formatTime(n, &
         &    '! Mean wall time for pure m loop            :',this%m_loop)
         call formatTime(n, &
         &    '! Mean wall time for m loop with matrix calc:',this%m_loop_mat)
         call formatTime(n, &
         &    '! Mean wall time for MPI communications     :',this%mpi_comms)
         call formatTime(n, &
         &    '! Mean wall time for output writting        :',this%io)
         call formatTime(n, &
         &    '! Mean wall time for one single FFT (rloop) :',this%fft)
         if ( this%n_dct_calls /= 0 ) then
            call formatTime(n, &
            &    '! Mean wall time for one single 2D-DCT      :',this%dct)
         end if
         call formatTime(n, &
         &    '! Mean wall time for one LU factor. (psi)   :',this%lu)
         call formatTime(n, &
         &    '! Mean wall time for one linear solve (psi) :',this%solve)
         call formatTime(n, &
         &    '! Mean wall time for one time step          :',this%tot)
      end if

   end subroutine write_log
!------------------------------------------------------------------------------
   subroutine initialize_timer(this)

      class(timer_type) :: this

      this%n_counts=0
      this%tStart = 0.0_cp
      this%tTot = 0.0_cp

   end subroutine initialize_timer
!-----------------------------------------------------------------------------
   subroutine finalize_timer(this, message, n_log_file)

      class(timer_type) :: this

      character(len=*), optional, intent(in) :: message
      integer, optional,          intent(in) :: n_log_file

      !-- Local variables
      integer :: n, n_out

      if ( this%n_counts > 0 ) this%tTot=this%tTot/real(this%n_counts,cp)

      call MPI_AllReduce(MPI_IN_PLACE, this%tTot, 1, MPI_DEF_REAL, MPI_SUM, &
           &             MPI_COMM_WORLD, ierr)
      this%tTot=this%tTot/real(n_procs,cp)

      if ( rank == 0 .and. present(n_log_file) .and. present(message) ) then
         if ( this%n_counts > 0 ) then
            do n=1,2
               if ( n == 1 ) then
                  n_out=n_log_file
               else if ( n == 2 ) then
                  n_out=output_unit
               end if
               call formatTime(n_out, message, this%tTot)
            end do
         end if
      end if

   end subroutine finalize_timer
!-----------------------------------------------------------------------------
   subroutine start_count(this)

      class(timer_type) :: this

      this%tStart = MPI_Wtime()

   end subroutine start_count
!-----------------------------------------------------------------------------
   subroutine stop_count(this, l_increment)

      class(timer_type) :: this
      logical, optional, intent(in) :: l_increment

      logical :: l_count
      real(cp) :: tStop

      if ( present(l_increment) ) then
         l_count=l_increment
      else
         l_count=.true.
      end if

      tStop = MPI_Wtime()

      if ( tStop > this%tStart ) then
         this%tTot=this%tTot+(tStop-this%tStart)
         if ( l_count ) this%n_counts=this%n_counts+1
      end if

   end subroutine stop_count
!-----------------------------------------------------------------------------
end module timers_mod
