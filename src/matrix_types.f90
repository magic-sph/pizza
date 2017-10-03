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

      procedure :: initialize => initialize_band
      procedure :: finalize => finalize_band
      procedure :: mat_vec_mul

   end type type_bandmat_real

   type, public :: type_bordmat_real

      integer :: kl ! Number of lower diagonals
      integer :: ku ! Number of upper diagonals
      integer :: nbands ! Number of bands
      integer :: nlines_band ! Number of lines
      integer :: ntau        ! Number of tau lines
      real(cp), allocatable :: A1(:,:) ! Upper left block
      real(cp), allocatable :: A2(:,:) ! Upper right block
      real(cp), allocatable :: A3(:,:) ! Lower left block
      real(cp), allocatable :: A4(:,:) ! Lower right block

   contains

      procedure :: initialize => initialize_bord
      procedure :: finalize => finalize_bord

   end type type_bordmat_real

contains

   subroutine initialize_band(this, kl, ku, len)

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

   end subroutine initialize_band
!------------------------------------------------------------------------------
   subroutine finalize_band(this)

      class(type_bandmat_real) :: this

      deallocate( this%dat )

   end subroutine finalize_band
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
   subroutine initialize_bord(this, kl, ku, nbounds, nlines)

      class(type_bordmat_real) :: this

      !-- Input variables
      integer, intent(in) :: kl
      integer, intent(in) :: ku
      integer, intent(in) :: nbounds
      integer, intent(in) :: nlines

      !-- Local variables
      integer :: n_r, n_b

      this%kl = kl
      this%ku = ku
      this%ntau = nbounds
      this%nlines_band = nlines-nbounds
      this%nbands = this%kl+this%ku+1

      allocate( this%A1(this%ntau, this%ntau) )
      allocate( this%A2(this%ntau, this%nlines_band) )
      allocate( this%A3(this%nlines_band, this%ntau) )
      allocate( this%A4(this%nbands+this%kl, this%nlines_band) )

      do n_r=1,nlines
         do n_b=1,this%ntau
            if ( n_r <= this%ntau ) then
               this%A1(n_b,n_r)=0.0_cp
            else
               this%A2(n_b,n_r-this%ntau)=0.0_cp
            end if
         end do
      end do

      do n_r=1,this%nlines_band
         do n_b=1,this%ntau
            this%A3(n_b,n_r)=0.0_cp
         end do
         do n_b=1,this%nbands+this%kl
            this%A4(n_b,n_r)=0.0_cp
         end do
      end do

      bytes_allocated=bytes_allocated+this%ntau*nlines*SIZEOF_DEF_REAL+&
      &               this%ntau*this%nlines_band*SIZEOF_DEF_REAL+      &
      &               (this%nbands+this%kl)*this%nlines_band*          &
      &               SIZEOF_DEF_REAL

   end subroutine initialize_bord
!------------------------------------------------------------------------------
   subroutine finalize_bord(this)

      class(type_bordmat_real) :: this

      deallocate( this%A1, this%A2, this%A3, this%A4 )

   end subroutine finalize_bord
!------------------------------------------------------------------------------
end module matrix_types
