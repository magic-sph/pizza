module precision_mod

   use mpimod

   implicit none

   private

#if (DEFAULT_PRECISION==sngl)
   integer, public, parameter :: cp=selected_real_kind(6)
   integer, public, parameter :: MPI_DEF_REAL=MPI_REAL4
   integer, public, parameter :: MPI_DEF_COMPLEX=MPI_COMPLEX8
#elif (DEFAULT_PRECISION==dble)
   integer, public, parameter :: cp=selected_real_kind(15)
   integer, public, parameter :: MPI_DEF_REAL=MPI_REAL8
   integer, public, parameter :: MPI_DEF_COMPLEX=MPI_COMPLEX16
#endif

   !-- SIZEOFs
   integer, public, parameter :: SIZEOF_DEF_COMPLEX=2*cp
   integer, public, parameter :: SIZEOF_DEF_REAL=cp

   !-- Precision for long integers
   integer, public, parameter :: lip=selected_int_kind(12)
   integer, public, parameter :: SIZEOF_INTEGER=4
   integer, public, parameter :: SIZEOF_LOGICAL=4
   integer, public, parameter :: SIZEOF_CHARACTER=1

end module precision_mod
