module parallel_mod

   use mpimod
#ifdef WITH_OMP
   use omp_lib
#endif

   implicit none

   integer :: rank, n_procs, ierr, n_threads

contains

   subroutine initialize_mpi

      integer :: required_level, provided_level

#ifdef WITH_OMP
      required_level = MPI_THREAD_FUNNELED
      ! required_level = MPI_THREAD_MULTIPLE
      call MPI_init_thread(required_level, provided_level, ierr)

      if ( provided_level < required_level ) then
         call MPI_Abort(MPI_COMM_WORLD, 32, ierr)
      end if
#else
      call MPI_init(ierr)
#endif
      call MPI_comm_rank(MPI_COMM_WORLD,rank,ierr)
      call MPI_comm_size(MPI_COMM_WORLD,n_procs, ierr)

#ifdef WITH_OMP
      n_threads = omp_get_max_threads()
#else
      n_threads = 1
#endif

   end subroutine initialize_mpi
!-----------------------------------------------------------------------------
   subroutine finalize_mpi

      call MPI_finalize(ierr)

   end subroutine finalize_mpi
!-----------------------------------------------------------------------------
end module parallel_mod
