module matrix_types

   use precision_mod
   use constants, only: one
   use mem_alloc, only: bytes_allocated

   implicit none

   private

   type, public :: type_bandmat_real

      integer :: kl ! Number of lower diagonals
      integer :: ku ! Number of upper diagonals
      integer :: nbands ! Number of bands
      integer :: nlines    ! Number of lines
      real(cp), allocatable :: dat(:,:) ! data

   contains

      procedure :: initialize
      procedure :: finalize
      procedure :: mat_vec_mul

   end type type_bandmat_real

contains

   subroutine initialize(this, kl, ku, len)

      class(type_bandmat_real) :: this

      !-- Input variables
      integer, intent(in) :: kl
      integer, intent(in) :: ku
      integer, intent(in) :: len

      !-- Local variables
      integer :: n_r, n_b

      this%kl = kl
      this%ku = ku
      this%nlines = len
      this%nbands = kl+ku+1

      allocate( this%dat(this%nbands,this%nlines))

      !-- Fill it with zeros
      do n_r=1,this%nlines
         do n_b=1,this%nbands
            this%dat(n_b,n_r)=0.0_cp
         end do
      end do

      bytes_allocated = bytes_allocated+this%nbands*this%nlines*SIZEOF_DEF_REAL

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)

      class(type_bandmat_real) :: this

      deallocate( this%dat )

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine mat_vec_mul(this, vec)

      class(type_bandmat_real) :: this

      !-- Input/output variables
      complex(cp), intent(inout) :: vec(this%nlines)

      !-- Local variables:
      integer :: n_r
      real(cp) :: vecr(this%nlines), veci(this%nlines)
      real(cp) :: tmpr(this%nlines), tmpi(this%nlines)

      do n_r=1,this%nlines
         vecr(n_r) = real(vec(n_r))
         tmpr(n_r) = vecr(n_r)
         veci(n_r) =aimag(vec(n_r))
         tmpi(n_r) = veci(n_r)
      end do

      call dgbmv('N', this%nlines, this%nlines, this%kl, this%ku, one, &
           &      this%dat, this%nbands, tmpr, 1, 0.0_cp, vecr, 1)
      call dgbmv('N', this%nlines, this%nlines, this%kl, this%ku, one, &
           &      this%dat, this%nbands, tmpi, 1, 0.0_cp, veci, 1)

      do n_r=1,this%nlines
         vec(n_r) = cmplx(vecr(n_r), veci(n_r), kind=cp)
      end do

   end subroutine mat_vec_mul
!------------------------------------------------------------------------------
end module matrix_types
