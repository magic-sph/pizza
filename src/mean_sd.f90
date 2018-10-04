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
      procedure :: initialize => initialize_1D
      procedure :: finalize => finalize_1D
   end type mean_sd_type

   type, public :: mean_sd_2D_type
      logical :: l_SD
      real(cp), allocatable :: mean(:,:)
      real(cp), allocatable :: SD(:,:)
   contains
      procedure :: initialize => initialize_2D
      procedure :: finalize => finalize_2D
   end type mean_sd_2D_type

contains

   subroutine initialize_1D(this, n_start, n_stop)
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

   end subroutine initialize_1D
!------------------------------------------------------------------------------
   subroutine finalize_1D(this)
      !
      ! Memory deallocation
      !
      class(mean_sd_type) :: this

      deallocate( this%SD, this%mean )

   end subroutine finalize_1D
!------------------------------------------------------------------------------
   subroutine initialize_2D(this, n_start, n_stop, n_in,l_SD)
      !
      ! Memory allocation
      !
      class(mean_sd_2D_type) :: this

      !-- Input variables:
      logical, intent(in) :: l_SD
      integer, intent(in) :: n_start
      integer, intent(in) :: n_stop
      integer, intent(in) :: n_in

      this%l_SD = l_SD

      allocate( this%mean(n_start:n_stop,n_in) )
      bytes_allocated=bytes_allocated+(n_stop-n_start+1)*n_in*SIZEOF_DEF_REAL
      this%mean(:,:)=0.0_cp
      
      if ( l_SD ) then
         allocate( this%SD(n_start:n_stop,n_in) )
         bytes_allocated=bytes_allocated+(n_stop-n_start+1)*n_in*SIZEOF_DEF_REAL
         this%SD(:,:)=0.0_cp
      end if

   end subroutine initialize_2D
!------------------------------------------------------------------------------
   subroutine finalize_2D(this)
      !
      ! Memory deallocation
      !
      class(mean_sd_2D_type) :: this

      if ( this%l_SD) deallocate(this%SD)
      deallocate( this%mean)

   end subroutine finalize_2D
!------------------------------------------------------------------------------
end module mean_sd
