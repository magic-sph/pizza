module time_array

   use precision_mod
   use constants, only: zero
   use mem_alloc, only: bytes_allocated

   implicit none

   private

   type, public :: type_tarray
      complex(cp), allocatable :: impl(:,:,:)
      complex(cp), pointer :: expl(:,:,:)
      complex(cp), allocatable :: old(:,:,:)
      logical :: l_exp
   contains
      procedure :: initialize
      procedure :: set_initial_values
      procedure :: finalize
   end type type_tarray

contains

   subroutine initialize(this, nMstart, nMstop, n_r_max, norder_imp, norder_exp, &
              &          norder_imp_lin, l_allocate_exp)

      class(type_tarray) :: this

      !-- Input variables
      integer, intent(in) :: nMstart
      integer, intent(in) :: nMstop
      integer, intent(in) :: n_r_max
      integer, intent(in) :: norder_imp
      integer, intent(in) :: norder_exp
      integer, intent(in) :: norder_imp_lin
      logical, optional, intent(in) :: l_allocate_exp ! A boolean to specify whether the explicit state has to be allocated

      !-- Local variable
      logical :: l_allocate

      if ( present(l_allocate_exp) ) then
         l_allocate = l_allocate_exp
      else
         l_allocate = .false.
      end if
      this%l_exp = l_allocate

      allocate( this%impl(nMstart:nMstop,n_r_max,norder_imp_lin-1) )
      allocate( this%old(nMstart:nMstop,n_r_max,norder_imp-1) )
      bytes_allocated = bytes_allocated + (nMstop-nMstart+1)*n_r_max*(  &
      &                 norder_imp+norder_imp_lin-2)*SIZEOF_DEF_COMPLEX

      if ( l_allocate ) then
         allocate( this%expl(nMstart:nMstop,n_r_max,norder_exp) )
         bytes_allocated = bytes_allocated + (nMstop-nMstart+1)*n_r_max*norder_exp &
         &                 *SIZEOF_DEF_COMPLEX
      end if

      call this%set_initial_values()

   end subroutine initialize
!----------------------------------------------------------------------------------
   subroutine set_initial_values(this)

      class(type_tarray) :: this

      this%old(:,:,:) =zero
      this%impl(:,:,:)=zero
      if ( this%l_exp ) this%expl(:,:,:)=zero

   end subroutine set_initial_values
!----------------------------------------------------------------------------------
   subroutine finalize(this)

      class(type_tarray) :: this

      if ( this%l_exp ) deallocate( this%expl )
      deallocate( this%old, this%impl )

   end subroutine finalize
!----------------------------------------------------------------------------------
end module time_array
