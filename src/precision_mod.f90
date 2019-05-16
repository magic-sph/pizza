module precision_mod

   use mpimod
   use iso_fortran_env, only: real64, real32, int32, int64

   implicit none

   private

#if (DEFAULT_PRECISION==sngl)
   integer, public, parameter :: cp=real32
   integer, public, parameter :: MPI_DEF_REAL=MPI_REAL4
   integer, public, parameter :: MPI_DEF_COMPLEX=MPI_COMPLEX8
#elif (DEFAULT_PRECISION==dble)
   integer, public, parameter :: cp=real64
   integer, public, parameter :: MPI_DEF_REAL=MPI_REAL8
   integer, public, parameter :: MPI_DEF_COMPLEX=MPI_COMPLEX16
#endif

   integer, public, parameter :: outp=real32

   !-- SIZEOFs
   integer, public, parameter :: SIZEOF_DEF_COMPLEX=2*cp
   integer, public, parameter :: SIZEOF_DEF_REAL=cp

   !-- Precision for long integers
   integer, public, parameter :: lip=int64
   integer, public, parameter :: SIZEOF_INTEGER=int32
   integer, public, parameter :: SIZEOF_LOGICAL=int32
   integer, public, parameter :: SIZEOF_CHARACTER=1

end module precision_mod
