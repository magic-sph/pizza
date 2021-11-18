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
end module parallel_mod
