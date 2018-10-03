module mean_sd
   !
   ! This module contains a small type that simply handles two arrays (mean and SD)
   ! This type is used for time-averaged outputs (and their standard deviations).
   !

   use mem_alloc
   use precision_mod

   implicit none

   private

   type, public :: mean_sd_type
      real(cp), allocatable :: mean(:)
      real(cp), allocatable :: SD(:)
   contains
      procedure :: initialize
      procedure :: finalize
   end type mean_sd_type

contains

   subroutine initialize(this, n_start, n_stop)
      !
      ! Memory allocation
      !
      class(mean_sd_type) :: this

      !-- Input variables:
      integer, intent(in) :: n_start
      integer, intent(in) :: n_stop

      allocate( this%mean(n_start:n_stop), this%SD(n_start:n_stop) )
      bytes_allocated=bytes_allocated+2*(n_stop-n_start+1)*SIZEOF_DEF_REAL

      this%mean(:)=0.0_cp
      this%SD(:)=0.0_cp

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation
      !

      class(mean_sd_type) :: this

      deallocate( this%SD, this%mean )

   end subroutine finalize
!------------------------------------------------------------------------------
end module mean_sd
