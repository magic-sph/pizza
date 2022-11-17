module parallel_mod

   use mpimod

   implicit none

   integer :: rank, n_procs, ierr

contains

   subroutine initialize_mpi

      call MPI_init(ierr)
      call MPI_comm_rank(MPI_COMM_WORLD,rank,ierr)
      call MPI_comm_size(MPI_COMM_WORLD,n_procs, ierr)

   end subroutine initialize_mpi
!-----------------------------------------------------------------------------
   subroutine finalize_mpi

      call MPI_finalize(ierr)

   end subroutine finalize_mpi
!-----------------------------------------------------------------------------
   subroutine get_openmp_blocks(nStart, nStop)

      !--Input/Outputs variables:
      integer, intent(inout) :: nStart
      integer, intent(inout) :: nStop

      !-- Local variables
      integer :: n_threads, threadid, n_points_per_thread, n_points_left
      integer :: n_points, n_glob_start
#ifdef WITHOMP
      integer :: n_max_threads
#endif

      n_points=nStop-nStart+1
      n_glob_start=nStart

#ifdef WITHOMP
      n_threads=omp_get_num_threads()
      threadid =omp_get_thread_num()
      if ( n_points < n_threads) then
         call omp_set_num_threads(n_points)
         n_points_per_thread=1
         n_points_left=0
      else
         n_points_per_thread=n_points/n_threads
         n_points_left=n_points-n_threads*n_points_per_thread
      end if
#else
      n_threads=1
      threadid =0
      n_points_per_thread=n_points
      n_points_left=0
#endif

      !-- This is a way to reshuffle the points which are not in-balance
      !-- more regularly
      if ( threadid+1 <= n_points_left ) then
         nStart = n_glob_start+threadid*n_points_per_thread+threadid
         nStop  = nStart+n_points_per_thread
      else
         nStart = n_glob_start+threadid*n_points_per_thread+n_points_left
         nStop  = nStart+n_points_per_thread-1
      end if

#ifdef WITHOMP
      if ( n_points < n_threads) then
         n_max_threads=omp_get_max_threads()
         call omp_set_num_threads(n_max_threads)
      end if
#endif

   end subroutine get_openmp_blocks
!-----------------------------------------------------------------------------
   subroutine mpiio_setup(info)
      !
      ! This routine set ups the default MPI-IO configuration. This is based
      ! on recommandations from IDRIS "Best practices for parallel IO and
      ! MPI-IO hints"
      !

      integer, intent(out) :: info

      call MPI_Info_create(info, ierr)

      !-- Enable collective buffering
      call MPI_Info_set(info, "romio_cb_write", "automatic", ierr)
      call MPI_Info_set(info, "romio_cb_read", "automatic", ierr)

      !-- Disable data sieving (let the filesystem handles it)
      call MPI_Info_set(info, "romio_ds_write", "disable", ierr)
      call MPI_Info_set(info, "romio_ds_read", "disable", ierr)

      !-- Set the striping unit to 4M
      call MPI_Info_set(info, "striping_unit", "4194304", ierr)

      !-- Set the striping factor to 64
      !call MPI_Info_set(info, "striping_factor", "64", ierr)

      !-- Set the buffer size to 4M
      call MPI_Info_set(info,"cb_buffer_size","4194304", ierr)

   end subroutine  mpiio_setup
!------------------------------------------------------------------------------
end module parallel_mod
