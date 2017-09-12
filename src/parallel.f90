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
end module parallel_mod
