module timers_mod
   !
   ! This module is used to set several counters that are use to measure
   ! the time passed in several critical parts of the code
   !

   use precision_mod
   use parallel_mod, only : rank
   use communications, only: my_reduce_mean
   use useful, only: formatTime

   implicit none

   private

   type, public :: timers_type
      integer :: n_r_loops
      integer :: n_r_loops_3D
      integer :: n_m_loops
      integer :: n_lm_loops
      integer :: n_m_loops_mat
      integer :: n_lm_loops_mat
      integer :: n_mpi_comms
      integer :: n_interp
      integer :: n_zavg
      integer :: n_zder
      integer :: n_io_calls
      integer :: n_fft_calls
      integer :: n_solve_calls
      integer :: n_dct_calls
      integer :: n_lu_calls
      real(cp) :: r_loop      ! Timer of the radial loop (nonlinear terms)
      real(cp) :: r_loop_3D   ! Timer of the 3-D radial loop (nonlinear terms)
      real(cp) :: m_loop      ! Timer of the m loop (linear solve/time advance)
      real(cp) :: lm_loop     ! Timer of the LM loop (linear solve/time advance)
      real(cp) :: io          ! Timer for I/O
      real(cp) :: interp      ! Timer for 2D -> 3D interpolations
      real(cp) :: zavg        ! Timer for zavg integrations
      real(cp) :: zder        ! Timer for zder derivations
      real(cp) :: mpi_comms   ! Timer for the main all_to_all communications
      real(cp) :: m_loop_mat  ! Timer for the m loop when matrix LU fac is needed
      real(cp) :: lm_loop_mat ! Timer for the LM loop when matrix LU fac is needed
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

contains

   subroutine initialize(this)
      !
      ! Initialize the counters with zeroes
      !

      class(timers_type) :: this

      this%n_r_loops      = 0
      this%n_r_loops_3D   = 0
      this%n_m_loops      = 0
      this%n_lm_loops     = 0
      this%n_m_loops_mat  = 0
      this%n_lm_loops_mat = 0
      this%n_mpi_comms    = 0
      this%n_io_calls     = 0
      this%n_fft_calls    = 0
      this%n_solve_calls  = 0
      this%n_dct_calls    = 0
      this%n_lu_calls     = 0
      this%n_interp       = 0
      this%n_zavg       = 0
      this%n_zder       = 0
      this%r_loop      = 0.0_cp
      this%r_loop_3D   = 0.0_cp
      this%m_loop      = 0.0_cp
      this%lm_loop     = 0.0_cp
      this%interp      = 0.0_cp
      this%zavg        = 0.0_cp
      this%zder        = 0.0_cp
      this%io          = 0.0_cp
      this%mpi_comms   = 0.0_cp
      this%m_loop_mat  = 0.0_cp
      this%lm_loop_mat = 0.0_cp
      this%fft         = 0.0_cp
      this%lu          = 0.0_cp
      this%solve       = 0.0_cp
      this%dct         = 0.0_cp
      this%tot         = 0.0_cp

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
      if ( this%n_r_loops_3D /= 0 ) then
         this%r_loop_3D = this%r_loop_3D/this%n_r_loops_3D
         call my_reduce_mean(this%r_loop_3D, 0)
      end if
      if ( this%n_m_loops /= 0 ) then
         this%m_loop    = this%m_loop/this%n_m_loops
         call my_reduce_mean(this%m_loop, 0)
      end if
      if ( this%n_lm_loops /= 0 ) then
         this%lm_loop    = this%lm_loop/this%n_lm_loops
         call my_reduce_mean(this%lm_loop, 0)
      end if
      if ( this%n_m_loops_mat /= 0 ) then
         this%m_loop_mat=this%m_loop_mat/this%n_m_loops_mat
         call my_reduce_mean(this%m_loop_mat, 0)
      end if
      if ( this%n_lm_loops_mat /= 0 ) then
         this%lm_loop_mat=this%lm_loop_mat/this%n_lm_loops_mat
         call my_reduce_mean(this%lm_loop_mat, 0)
      end if
      if ( this%n_interp /= 0 ) then
         this%interp=this%interp/this%n_interp
         call my_reduce_mean(this%interp, 0)
      end if
      if ( this%n_zavg /= 0 ) then
         this%zavg=this%zavg/this%n_zavg
         call my_reduce_mean(this%zavg, 0)
      end if
      if ( this%n_zder /= 0 ) then
         this%zder=this%zder/this%n_zder
         call my_reduce_mean(this%zder, 0)
      end if
      if ( n_steps /= 0 ) then
         this%tot=this%tot/n_steps
         call my_reduce_mean(this%tot, 0)
      end if
      this%mpi_comms = this%mpi_comms/this%n_mpi_comms
      call my_reduce_mean(this%mpi_comms, 0)
      this%fft = this%fft/this%n_fft_calls
      call my_reduce_mean(this%fft, 0)
      this%lu = this%lu/this%n_lu_calls
      call my_reduce_mean(this%lu, 0)
      this%solve = this%solve/this%n_solve_calls
      call my_reduce_mean(this%solve, 0)
      this%dct = this%dct/this%n_dct_calls
      call my_reduce_mean(this%dct, 0)

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
         &    '! Mean wall time for radial loop             :',this%r_loop)
         call formatTime(n, &
         &    '! Mean wall time for pure m loop             :',this%m_loop)
         call formatTime(n, &
         &    '! Mean wall time for m loop with matrix calc :',this%m_loop_mat)
         if ( this%n_r_loops_3D > 0 ) then
            call formatTime(n, &
            &    '! Mean wall time for 3D radial loop          :',this%r_loop_3D)
            call formatTime(n, &
            &    '! Mean wall time for pure LM loop            :',this%lm_loop)
            call formatTime(n, &
            &    '! Mean wall time for LM loop with matrix calc:',this%lm_loop_mat)
            call formatTime(n, &
            &    '! Mean wall time 2D/3D interpolations        :',this%interp)
            call formatTime(n, &
            &    '! Mean wall time zavg integrations           :',this%zavg)
            call formatTime(n, &
            &    '! Mean wall time zder derivations            :',this%zder)
         end if
         call formatTime(n, &
         &    '! Mean wall time for MPI communications      :',this%mpi_comms)
         call formatTime(n, &
         &    '! Mean wall time for output writting         :',this%io)
         call formatTime(n, &
         &    '! Mean wall time for one single FFT (rloop)  :',this%fft)
         call formatTime(n, &
         &    '! Mean wall time for one single 2D-DCT       :',this%dct)
         call formatTime(n, &
         &    '! Mean wall time for one LU factor. (psi)    :',this%lu)
         call formatTime(n, &
         &    '! Mean wall time for one linear solve (psi)  :',this%solve)
         call formatTime(n, &
         &    '! Mean wall time for one time step           :',this%tot)
      end if

   end subroutine write_log
!------------------------------------------------------------------------------
end module timers_mod
