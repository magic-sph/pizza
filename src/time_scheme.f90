module time_scheme

   use precision_mod
   use namelists, only: alpha
   use constants, only: one, half
   use mem_alloc, only: bytes_allocated

   implicit none

   private

   type, public :: type_tscheme
      integer :: norder
      character(len=3) :: imp_scheme
      character(len=3) :: exp_scheme
      real(cp), allocatable :: dt(:)
      real(cp), allocatable :: wimp(:)
      real(cp), allocatable :: wexp(:)

   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_weights
      procedure :: set_dt_array
   end type type_tscheme

contains

   subroutine initialize(this)

      class(type_tscheme) :: this

      this%norder = 2
      this%imp_scheme = 'CN'
      this%exp_scheme = 'AB2'

      allocate ( this%dt(this%norder) )
      allocate ( this%wimp(this%norder) )
      allocate ( this%wexp(this%norder) )

      bytes_allocated = bytes_allocated+3*this%norder*SIZEOF_DEF_REAL

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)

      class(type_tscheme) :: this

      deallocate( this%dt, this%wimp, this%wexp )

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine set_weights(this)

      class(type_tscheme) :: this

      select case ( this%imp_scheme )
         case ('CN') 
            this%wimp(1)=alpha*this%dt(1)
            this%wimp(2)=(1-alpha)*this%dt(1)
      end select

      select case ( this%exp_scheme )
         case ('AB2') 
            this%wexp(2)=-half*this%dt(1)/this%dt(2)
            this%wexp(1)=one-this%wexp(2)
      end select

   end subroutine set_weights
!------------------------------------------------------------------------------
   subroutine set_dt_array(this, dt_new)

      class(type_tscheme) :: this
      real(cp), intent(in) :: dt_new

      this%dt   =cshift(this%dt,shift=this%norder-1)
      this%dt(1)=dt_new

   end subroutine set_dt_array
!------------------------------------------------------------------------------
end module time_scheme
