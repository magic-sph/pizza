module matrix_types
   !
   ! This module defines several types for handling band and bordered matrices
   !

   use precision_mod
   use constants, only: one, zero
   use mem_alloc, only: bytes_allocated
   use algebra, only: prepare_bordered_mat, solve_bordered_mat, prepare_band_mat, &
       &              solve_band_mat

   implicit none

   private

   type, public :: type_bandmat_real
      integer :: kl ! Number of lower diagonals
      integer :: ku ! Number of upper diagonals
      integer :: nbands ! Number of bands
      integer :: nlines    ! Number of lines
      logical :: l_lhs_type ! L.H.S. matrix or not
      integer, allocatable :: piv(:) ! pivot
      real(cp), allocatable :: dat(:,:) ! data
   contains
      procedure :: initialize => initialize_band_real
      procedure :: finalize => finalize_band_real
      procedure :: remove_last_rows => remove_last_rows_band_real
      procedure :: remove_leading_blank_rows => remove_leading_blank_rows_band_real
      procedure :: prepare_LU => prepare_LU_band_real
      procedure :: solve_band_real_rhs_real
      procedure :: solve_band_real_rhs_complex
      procedure :: mat_real_vec_complex_mul
      procedure :: mat_real_vec_real_mul
      generic :: mat_vec_mul => mat_real_vec_complex_mul, mat_real_vec_real_mul
      generic :: solve => solve_band_real_rhs_real, solve_band_real_rhs_complex
   end type type_bandmat_real

   type, public :: type_bandmat_complex
      integer :: kl ! Number of lower diagonals
      integer :: ku ! Number of upper diagonals
      integer :: nbands ! Number of bands
      integer :: nlines    ! Number of lines
      logical :: l_lhs_type ! L.H.S. matrix or not
      integer, allocatable :: piv(:) ! pivot
      complex(cp), allocatable :: dat(:,:) ! data
   contains
      procedure :: initialize => initialize_band_complex
      procedure :: finalize => finalize_band_complex
      procedure :: remove_last_rows => remove_last_rows_band_complex
      procedure :: remove_leading_blank_rows => remove_leading_blank_rows_band_complex
      procedure :: prepare_LU => prepare_LU_band_complex
      procedure :: solve => solve_band_complex_rhs_complex
      procedure :: mat_vec_mul => mat_complex_vec_complex_mul
   end type type_bandmat_complex

   type, public :: type_bordmat_real
      integer :: kl ! Number of lower diagonals
      integer :: ku ! Number of upper diagonals
      integer :: nbands ! Number of bands
      integer :: nlines_band ! Number of lines of the banded block
      integer :: nlines      ! Number of lines 
      integer :: ntau        ! Number of tau lines
      real(cp), allocatable :: A1(:,:) ! Upper left block
      real(cp), allocatable :: A2(:,:) ! Upper right block
      real(cp), allocatable :: A3(:,:) ! Lower left block
      real(cp), allocatable :: A4(:,:) ! Lower right block
      integer, allocatable :: pivA1(:) ! Pivot for first block
      integer, allocatable :: pivA4(:) ! Pivot for fourth block
   contains
      procedure :: initialize => initialize_bord_real
      procedure :: finalize => finalize_bord_real
      procedure :: prepare_LU => prepare_LU_real
      procedure :: solve_real_mat_complex_rhs
      procedure :: solve_real_mat_real_rhs
      generic :: solve => solve_real_mat_complex_rhs, solve_real_mat_real_rhs
   end type type_bordmat_real

   type, public :: type_bordmat_complex
      integer :: kl ! Number of lower diagonals
      integer :: ku ! Number of upper diagonals
      integer :: nbands ! Number of bands
      integer :: nlines_band ! Number of lines of the banded block
      integer :: nlines      ! Number of lines 
      integer :: ntau        ! Number of tau lines
      complex(cp), allocatable :: A1(:,:) ! Upper left block
      complex(cp), allocatable :: A2(:,:) ! Upper right block
      complex(cp), allocatable :: A3(:,:) ! Lower left block
      complex(cp), allocatable :: A4(:,:) ! Lower right block
      integer, allocatable :: pivA1(:) ! Pivot for first block
      integer, allocatable :: pivA4(:) ! Pivot for fourth block
   contains
      procedure :: initialize => initialize_bord_complex
      procedure :: finalize => finalize_bord_complex
      procedure :: prepare_LU => prepare_LU_complex
      procedure :: solve => solve_complex_mat_complex_rhs
      procedure :: mat_vec_mul => bordmat_complex_vec_complex_mul
      procedure :: write => write_bordmat_complex
   end type type_bordmat_complex

   interface band_band_product
      module procedure :: band_real_band_real_product
      module procedure :: band_complex_band_real_product
   end interface band_band_product

   public :: band_band_product

contains

   subroutine initialize_band_real(this, kl, ku, len, l_lhs)

      class(type_bandmat_real) :: this

      !-- Input variables
      integer, intent(in) :: kl
      integer, intent(in) :: ku
      integer, intent(in) :: len
      logical, optional,  intent(in) :: l_lhs

      !-- Local variables
      logical :: l_lhs_loc

      if ( present(l_lhs) ) then
         l_lhs_loc = l_lhs
      else
         l_lhs_loc = .false.
      end if

      this%kl = kl
      this%ku = ku
      this%nlines = len
      this%nbands = kl+ku+1
      this%l_lhs_type = l_lhs_loc

      if ( this%l_lhs_type ) then
         allocate( this%dat(this%nbands+this%kl,this%nlines) )
         allocate( this%piv(this%nlines) )
         bytes_allocated = bytes_allocated+(this%nbands+this%kl)*this%nlines*&
         &                 SIZEOF_DEF_REAL+this%nlines*SIZEOF_INTEGER
      else
         allocate( this%dat(this%nbands,this%nlines) )
         bytes_allocated = bytes_allocated+this%nbands*this%nlines*SIZEOF_DEF_REAL
      end if

      !-- Fill it with zeros
      this%dat(:,:)=0.0_cp

   end subroutine initialize_band_real
!------------------------------------------------------------------------------
   subroutine finalize_band_real(this)

      class(type_bandmat_real) :: this

      if ( this%l_lhs_type ) deallocate( this%piv )
      deallocate( this%dat )

   end subroutine finalize_band_real
!------------------------------------------------------------------------------
   subroutine remove_last_rows_band_real(this, ncut)

      class(type_bandmat_real) :: this

      !-- Input variable
      integer, intent(in) :: ncut

      !-- Local variable
      real(cp), allocatable :: tmp(:,:)
      integer :: n_r, n_b

      allocate( tmp(this%nbands,this%nlines) )

      if ( this%l_lhs_type ) then
         do n_r=1,this%nlines
            do n_b=1,this%nbands
               tmp(n_b,n_r)=this%dat(this%kl+n_b,n_r)
            end do
         end do
      else
         do n_r=1,this%nlines
            do n_b=1,this%nbands
               tmp(n_b,n_r)=this%dat(n_b,n_r)
            end do
         end do
      end if

      deallocate( this%dat)

      this%nlines = this%nlines-ncut
      if ( this%l_lhs_type ) then
         deallocate( this%piv )
         allocate( this%piv(this%nlines) )
         allocate( this%dat(this%nbands+this%kl, this%nlines) )
         do n_r=1,this%nlines
            do n_b=1,this%nbands
               this%dat(this%kl+n_b,n_r)=tmp(n_b,n_r)
            end do
         end do
      else
         allocate( this%dat(this%nbands, this%nlines) )
         do n_r=1,this%nlines
            do n_b=1,this%nbands
               this%dat(n_b,n_r)=tmp(n_b,n_r)
            end do
         end do
      end if

      deallocate( tmp)

   end subroutine remove_last_rows_band_real
!------------------------------------------------------------------------------
   subroutine remove_leading_blank_rows_band_real(this,ncut)

      class(type_bandmat_real) :: this

      !-- Input variable
      integer, intent(in) :: ncut

      !-- Local variable
      real(cp), allocatable :: tmp(:,:)
      integer :: n_r, n_b

      if ( this%l_lhs_type ) then
         allocate( tmp(this%nbands,this%nlines) )
         do n_r=1,this%nlines
            do n_b=1,this%nbands
               tmp(n_b,n_r)=this%dat(this%kl+n_b,n_r)
            end do
         end do
         this%kl = max(this%kl-ncut,0)
         this%ku = this%ku+ncut
         deallocate( this%dat )
         allocate( this%dat(2*this%kl+this%ku+1,this%nlines) )
         this%dat(1:this%kl,:)=0.0_cp
         do n_r=1,this%nlines
            do n_b=1,this%nbands
               this%dat(this%kl+n_b,n_r)=tmp(n_b,n_r)
            end do
         end do
         deallocate( tmp )
      else
         this%kl = max(this%kl-ncut,0)
         this%ku = this%ku+ncut
      end if

   end subroutine remove_leading_blank_rows_band_real
!------------------------------------------------------------------------------
   subroutine prepare_LU_band_real(this)

      class(type_bandmat_real) :: this

      if ( this%l_lhs_type ) then
         call prepare_band_mat(this%dat,this%nlines,this%kl,this%ku,this%piv)
      end if

   end subroutine prepare_LU_band_real
!------------------------------------------------------------------------------
   subroutine solve_band_real_rhs_real(this, rhs, lenRhs)

      class(type_bandmat_real) :: this

      !-- Input variable
      integer, intent(in) :: lenRhs

      !-- Output variables
      real(cp), intent(inout) :: rhs(lenRhs)

      if ( this%l_lhs_type ) then
         call solve_band_mat(this%dat, this%nlines, this%kl, this%ku, this%piv, &
              &              rhs, lenRhs)
      end if

   end subroutine solve_band_real_rhs_real
!------------------------------------------------------------------------------
   subroutine solve_band_real_rhs_complex(this, rhs, lenRhs)

      class(type_bandmat_real) :: this

      !-- Input variable
      integer, intent(in) :: lenRhs

      !-- Output variables
      complex(cp), intent(inout) :: rhs(lenRhs)

      !-- Local variables
      integer :: n_r
      real(cp) :: rhsr(lenRhs), rhsi(lenRhs)

      if ( this%l_lhs_type ) then

         do n_r=1,lenRhs
            rhsr(n_r)= real(rhs(n_r))
            rhsi(n_r)=aimag(rhs(n_r))
         end do

         call solve_band_mat(this%dat, this%nlines, this%kl, this%ku, this%piv, &
              &              rhsr, lenRhs)
         call solve_band_mat(this%dat, this%nlines, this%kl, this%ku, this%piv, &
              &              rhsi, lenRhs)

         do n_r=1,lenRhs
            rhs(n_r)=cmplx(rhsr(n_r),rhsi(n_r),kind=cp)
         end do

      end if

   end subroutine solve_band_real_rhs_complex
!------------------------------------------------------------------------------
   subroutine mat_real_vec_complex_mul(this, vec)
      !
      ! This is a matrix-vector multiplication (real * complex)
      !

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

   end subroutine mat_real_vec_complex_mul
!------------------------------------------------------------------------------
   subroutine initialize_band_complex(this, kl, ku, len, l_lhs)

      class(type_bandmat_complex) :: this

      !-- Input variables
      integer, intent(in) :: kl
      integer, intent(in) :: ku
      integer, intent(in) :: len
      logical, optional,  intent(in) :: l_lhs

      !-- Local variables
      logical :: l_lhs_loc

      if ( present(l_lhs) ) then
         l_lhs_loc = l_lhs
      else
         l_lhs_loc = .false.
      end if

      this%kl = kl
      this%ku = ku
      this%nlines = len
      this%nbands = kl+ku+1
      this%l_lhs_type = l_lhs_loc

      if ( this%l_lhs_type ) then
         allocate( this%dat(this%nbands+this%kl,this%nlines) )
         allocate( this%piv(this%nlines) )
         bytes_allocated = bytes_allocated+(this%nbands+this%kl)*this%nlines*&
         &                 SIZEOF_DEF_REAL+this%nlines*SIZEOF_INTEGER
      else
         allocate( this%dat(this%nbands,this%nlines) )
         bytes_allocated = bytes_allocated+this%nbands*this%nlines*SIZEOF_DEF_REAL
      end if

      !-- Fill it with zeros
      this%dat(:,:)=zero

   end subroutine initialize_band_complex
!------------------------------------------------------------------------------
   subroutine finalize_band_complex(this)

      class(type_bandmat_complex) :: this

      if ( this%l_lhs_type ) deallocate( this%piv )
      deallocate( this%dat )

   end subroutine finalize_band_complex
!------------------------------------------------------------------------------
   subroutine mat_complex_vec_complex_mul(this, vec)
      !
      ! This is a matrix-vector multiplication (complex * complex)
      !

      class(type_bandmat_complex) :: this

      !-- Input/output variables
      complex(cp), intent(inout) :: vec(this%nlines)

      !-- Local variables:
      integer :: n_r
      complex(cp) :: tmp(this%nlines)

      do n_r=1,this%nlines
         tmp(n_r) = vec(n_r)
      end do

      call zgbmv('N', this%nlines, this%nlines, this%kl, this%ku, (one,0.0_cp), &
           &     this%dat, this%nbands, tmp, 1, zero, vec, 1)

   end subroutine mat_complex_vec_complex_mul
!------------------------------------------------------------------------------
   subroutine prepare_LU_band_complex(this)

      class(type_bandmat_complex) :: this

      if ( this%l_lhs_type ) then
         call prepare_band_mat(this%dat,this%nlines,this%kl,this%ku,this%piv)
      end if

   end subroutine prepare_LU_band_complex
!------------------------------------------------------------------------------
   subroutine solve_band_complex_rhs_complex(this, rhs, lenRhs)

      class(type_bandmat_complex) :: this

      !-- Input variable
      integer, intent(in) :: lenRhs

      !-- Output variables
      complex(cp), intent(inout) :: rhs(lenRhs)

      if ( this%l_lhs_type ) then
         call solve_band_mat(this%dat, this%nlines, this%kl, this%ku, this%piv, &
              &              rhs, lenRhs)
      end if

   end subroutine solve_band_complex_rhs_complex
!------------------------------------------------------------------------------
   subroutine mat_real_vec_real_mul(this, vec)
      !
      ! This is a matrix-vector multiplication (real * real)
      !

      class(type_bandmat_real) :: this

      !-- Input/output variables
      real(cp), intent(inout) :: vec(this%nlines)

      !-- Local variables:
      integer :: n_r
      real(cp) :: tmp(this%nlines)

      do n_r=1,this%nlines
         tmp(n_r) = vec(n_r)
      end do

      call dgbmv('N', this%nlines, this%nlines, this%kl, this%ku, one, &
           &      this%dat, this%nbands, tmp, 1, 0.0_cp, vec, 1)

   end subroutine mat_real_vec_real_mul
!------------------------------------------------------------------------------
   subroutine initialize_bord_real(this, kl, ku, nbounds, nlines)

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
      this%nlines = nlines
      this%nbands = this%kl+this%ku+1

      allocate( this%A1(this%ntau, this%ntau) )
      allocate( this%A2(this%ntau, this%nlines_band) )
      allocate( this%A3(this%nlines_band, this%ntau) )
      allocate( this%A4(this%nbands+this%kl, this%nlines_band) )

      allocate( this%pivA1(this%ntau) )
      allocate( this%pivA4(this%nlines_band) )

      do n_r=1,this%nlines
         do n_b=1,this%ntau
            if ( n_r <= this%ntau ) then
               this%A1(n_b,n_r)=0.0_cp
            else
               this%A2(n_b,n_r-this%ntau)=0.0_cp
            end if
         end do
      end do

      do n_b=1,this%ntau
         do n_r=1,this%nlines_band
            this%A3(n_r,n_b)=0.0_cp
         end do
      end do

      do n_r=1,this%nlines_band
         do n_b=1,this%nbands+this%kl
            this%A4(n_b,n_r)=0.0_cp
         end do
      end do

      bytes_allocated=bytes_allocated+this%ntau*this%nlines*SIZEOF_DEF_REAL+&
      &               this%ntau*this%nlines_band*SIZEOF_DEF_REAL+           &
      &               (this%nbands+this%kl)*this%nlines_band*               &
      &               SIZEOF_DEF_REAL+this%nlines*SIZEOF_INTEGER

   end subroutine initialize_bord_real
!------------------------------------------------------------------------------
   subroutine finalize_bord_real(this)

      class(type_bordmat_real) :: this

      deallocate( this%A1, this%A2, this%A3, this%A4 )
      deallocate( this%pivA1, this%pivA4)

   end subroutine finalize_bord_real
!------------------------------------------------------------------------------
   subroutine prepare_LU_real(this)

      class(type_bordmat_real) :: this

      call prepare_bordered_mat(this%A1, this%A2, this%A3, this%A4, this%ntau, &
           &                    this%nlines_band, this%kl, this%ku, this%pivA1,&
           &                    this%pivA4)

   end subroutine prepare_LU_real
!------------------------------------------------------------------------------
   subroutine solve_real_mat_complex_rhs(this, rhs, nRmax)

      class(type_bordmat_real) :: this

      !-- Input variable
      integer, intent(in) :: nRmax

      !-- In/Out variable
      complex(cp), intent(inout) :: rhs(nRmax)

      !-- Local variable
      integer :: n_r
      real(cp) :: rhsr(nRmax), rhsi(nRmax)

      do n_r=1,nRmax
         rhsr(n_r)= real(rhs(n_r))
         rhsi(n_r)=aimag(rhs(n_r))
      end do

      call solve_bordered_mat(this%A1, this%A2, this%A3, this%A4, this%ntau, &
           &                  this%nlines_band, this%kl, this%ku, this%pivA1,&
           &                  this%pivA4, rhsr, nRmax)
      call solve_bordered_mat(this%A1, this%A2, this%A3, this%A4, this%ntau, &
           &                  this%nlines_band, this%kl, this%ku, this%pivA1,&
           &                  this%pivA4, rhsi, nRmax)

      do n_r=1,nRmax
         rhs(n_r)=cmplx(rhsr(n_r),rhsi(n_r),kind=cp)
      end do

   end subroutine solve_real_mat_complex_rhs
!------------------------------------------------------------------------------
   subroutine solve_real_mat_real_rhs(this, rhs, nRmax)

      class(type_bordmat_real) :: this

      !-- Input variable
      integer, intent(in) :: nRmax

      !-- In/Out variable
      real(cp), intent(inout) :: rhs(nRmax)

      call solve_bordered_mat(this%A1, this%A2, this%A3, this%A4, this%ntau, &
           &                  this%nlines_band, this%kl, this%ku, this%pivA1,&
           &                  this%pivA4, rhs, nRmax)

   end subroutine solve_real_mat_real_rhs
!------------------------------------------------------------------------------
   subroutine initialize_bord_complex(this, kl, ku, nbounds, nlines)

      class(type_bordmat_complex) :: this

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
      this%nlines = nlines
      this%nbands = this%kl+this%ku+1

      allocate( this%A1(this%ntau, this%ntau) )
      allocate( this%A2(this%ntau, this%nlines_band) )
      allocate( this%A3(this%nlines_band, this%ntau) )
      allocate( this%A4(this%nbands+this%kl, this%nlines_band) )

      allocate( this%pivA1(this%ntau) )
      allocate( this%pivA4(this%nlines_band) )

      do n_r=1,this%nlines
         do n_b=1,this%ntau
            if ( n_r <= this%ntau ) then
               this%A1(n_b,n_r)=zero
            else
               this%A2(n_b,n_r-this%ntau)=zero
            end if
         end do
      end do

      do n_b=1,this%ntau
         do n_r=1,this%nlines_band
            this%A3(n_r,n_b)=zero
         end do
      end do

      do n_r=1,this%nlines_band
         do n_b=1,this%nbands+this%kl
            this%A4(n_b,n_r)=zero
         end do
      end do

      bytes_allocated=bytes_allocated+this%ntau*this%nlines*SIZEOF_DEF_COMPLEX+&
      &               this%ntau*this%nlines_band*SIZEOF_DEF_COMPLEX+           &
      &               (this%nbands+this%kl)*this%nlines_band*                  &
      &               SIZEOF_DEF_COMPLEX+this%nlines*SIZEOF_INTEGER

   end subroutine initialize_bord_complex
!------------------------------------------------------------------------------
   subroutine finalize_bord_complex(this)

      class(type_bordmat_complex) :: this

      deallocate( this%A1, this%A2, this%A3, this%A4 )
      deallocate( this%pivA1, this%pivA4)

   end subroutine finalize_bord_complex
!------------------------------------------------------------------------------
   subroutine prepare_LU_complex(this)

      class(type_bordmat_complex) :: this

      call prepare_bordered_mat(this%A1, this%A2, this%A3, this%A4, this%ntau, &
           &                    this%nlines_band, this%kl, this%ku, this%pivA1,&
           &                    this%pivA4)

   end subroutine prepare_LU_complex
!------------------------------------------------------------------------------
   subroutine remove_last_rows_band_complex(this, ncut)

      class(type_bandmat_complex) :: this

      !-- Input variable
      integer, intent(in) :: ncut

      !-- Local variable
      complex(cp), allocatable :: tmp(:,:)
      integer :: n_r, n_b

      allocate( tmp(this%nbands,this%nlines) )

      if ( this%l_lhs_type ) then
         do n_r=1,this%nlines
            do n_b=1,this%nbands
               tmp(n_b,n_r)=this%dat(this%kl+n_b,n_r)
            end do
         end do
      else
         do n_r=1,this%nlines
            do n_b=1,this%nbands
               tmp(n_b,n_r)=this%dat(n_b,n_r)
            end do
         end do
      end if

      deallocate( this%dat)

      this%nlines = this%nlines-ncut
      if ( this%l_lhs_type ) then
         deallocate( this%piv )
         allocate( this%piv(this%nlines) )
         allocate( this%dat(this%nbands+this%kl, this%nlines) )
         do n_r=1,this%nlines
            do n_b=1,this%nbands
               this%dat(this%kl+n_b,n_r)=tmp(n_b,n_r)
            end do
         end do
      else
         allocate( this%dat(this%nbands, this%nlines) )
         do n_r=1,this%nlines
            do n_b=1,this%nbands
               this%dat(n_b,n_r)=tmp(n_b,n_r)
            end do
         end do
      end if

      deallocate( tmp)

   end subroutine remove_last_rows_band_complex
!------------------------------------------------------------------------------
   subroutine remove_leading_blank_rows_band_complex(this,ncut)

      class(type_bandmat_complex) :: this

      !-- Input variable
      integer, intent(in) :: ncut

      !-- Local variable
      complex(cp), allocatable :: tmp(:,:)
      integer :: n_r, n_b

      if ( this%l_lhs_type ) then
         allocate( tmp(this%nbands,this%nlines) )
         do n_r=1,this%nlines
            do n_b=1,this%nbands
               tmp(n_b,n_r)=this%dat(this%kl+n_b,n_r)
            end do
         end do
         this%kl = max(this%kl-ncut,0)
         this%ku = this%ku+ncut
         deallocate( this%dat )
         allocate( this%dat(2*this%kl+this%ku+1,this%nlines) )
         this%dat(1:this%kl,:)=0.0_cp
         do n_r=1,this%nlines
            do n_b=1,this%nbands
               this%dat(this%kl+n_b,n_r)=tmp(n_b,n_r)
            end do
         end do
         deallocate( tmp )
      else
         this%kl = max(this%kl-ncut,0)
         this%ku = this%ku+ncut
      end if

   end subroutine remove_leading_blank_rows_band_complex
!------------------------------------------------------------------------------
   subroutine solve_complex_mat_complex_rhs(this, rhs, nRmax)

      class(type_bordmat_complex) :: this

      !-- Input variable
      integer, intent(in) :: nRmax

      !-- In/Out variable
      complex(cp), intent(inout) :: rhs(nRmax)

      call solve_bordered_mat(this%A1, this%A2, this%A3, this%A4, this%ntau, &
           &                  this%nlines_band, this%kl, this%ku, this%pivA1,&
           &                  this%pivA4, rhs, nRmax)

   end subroutine solve_complex_mat_complex_rhs
!------------------------------------------------------------------------------
   subroutine bordmat_complex_vec_complex_mul(this, vec)
      !
      ! This is a matrix-vector multiplication (complex * complex)
      !

      class(type_bordmat_complex) :: this

      !-- Input/output variables
      complex(cp), intent(inout) :: vec(this%nlines)

      !-- Local variables:
      integer :: n_r
      complex(cp) :: tmp1(this%ntau), tmp2(this%nlines_band)

      do n_r=1,this%ntau
         tmp1(n_r) = vec(n_r)
      end do

      do n_r=1,this%nlines_band
         tmp2(n_r) = vec(n_r+this%ntau)
      end do

      !-- A1 * vec(1:ntau)
      call zgemv('N', this%ntau, this%ntau, (one, 0.0_cp), this%A1,&
           &      this%ntau, vec(1:this%ntau), 1, (0.0_cp,0.0_cp), tmp1, 1)

      !-- A2 * vec(ntau+1:nlines)
      call zgemv('N', this%ntau, this%nlines_band, (one, 0.0_cp), &
           &     this%A2, this%ntau, vec(this%ntau+1:), 1,        &
           &     (one,0.0_cp), tmp1, 1)

      !-- A3 * vec(1:ntau)
       call zgemv('N', this%nlines_band, this%ntau, (one, 0.0_cp), this%A3,&
            &      this%nlines_band, vec(1:this%ntau), 1, (0.0_cp,0.0_cp), &
            &      tmp2, 1)

      !-- A4 * vec(ntau+1:nlines)
      call zgbmv('N', this%nlines_band, this%nlines_band, this%kl, this%ku,      &
           &    (one,0.0_cp), this%A4(this%kl+1:,:), this%nbands, vec(this%ntau+1:), 1, &
           &    (1.0_cp,0.0_cp), tmp2, 1)

      do n_r=1,this%ntau
         vec(n_r)=tmp1(n_r)
      end do
! 
      do n_r=1,this%nlines_band
         vec(n_r+this%ntau)=tmp2(n_r)
      end do

   end subroutine bordmat_complex_vec_complex_mul
!------------------------------------------------------------------------------
   subroutine write_bordmat_complex(this)
      !
      ! This subroutine write a file A_mat that contains the bordered matrix.
      ! This is a Fortran unformatted file
      !
      class(type_bordmat_complex) :: this

      !-- Local variables:
      complex(cp) :: tmp(this%nlines)
      integer :: n_r, n_col, file_handle

      open(file_handle, file='A_mat', form='unformatted')

      write(file_handle) this%nlines

      !-- Top blocks (A1 and A2)
      do n_r=1,this%ntau

         do n_col=1,this%nlines
            if ( n_col <= this%ntau ) then
               tmp(n_col)=this%A1(n_r,n_col)
            else
               tmp(n_col)=this%A2(n_r,n_col-this%ntau)
            end if
         end do
         write(file_handle) tmp

      end do

      !-- Bottom blocks (A3 and A4)
      do n_r=1,this%nlines_band

         do n_col=1,this%nlines
            if ( n_col <= this%ntau ) then
               tmp(n_col)=this%A3(n_r,n_col)
            else
               if ( this%kl+this%kl+n_r-n_col+this%ntau > 0 .and.  this%kl+this%kl+n_r-n_col+this%ntau < this%nbands+this%kl) then
                  tmp(n_col)=this%A4(this%kl+this%kl+n_r-n_col+1+this%ntau,n_col-this%ntau)
               else
                  tmp(n_col)=zero
               end if
            end if
         end do
         write(file_handle) tmp

      end do

      close(file_handle)

   end subroutine write_bordmat_complex
!------------------------------------------------------------------------------
   subroutine band_real_band_real_product(A,B,C,l_lhs)
      !
      ! This subroutine computes a matrix multiplication between
      ! two input banded matrices
      !

      !-- Input variables
      type(type_bandmat_real), intent(in) :: A
      type(type_bandmat_real), intent(in) :: B
      logical,                 intent(in) :: l_lhs

      !-- Output variables
      type(type_bandmat_real), intent(inout) :: C

      !-- Local variables
      integer :: klC, kuC, nrmax, n_r
      integer :: o_a,o_c,o_b,d_a,d_b,d_c
      integer :: row_a,row_b,row_c,row_c_offset

      if ( .not. allocated(C%dat) ) then
         klC=A%kl+B%kl
         kuC=A%ku+B%ku
         nrmax =A%nlines
         call C%initialize(klC, kuC, nrmax, l_lhs)
      end if

      if ( l_lhs ) then
         row_c_offset = klC
      else
         row_c_offset = 0
      end if

      do o_c=-min(C%ku, A%ku+B%ku),min(C%kl,A%kl+B%kl)
         do o_a=-min(C%ku,B%kl-o_c),min(A%kl,B%ku+o_c)
            o_b=o_c-o_a
            row_a=A%ku+o_a+1
            row_b=B%ku+o_b+1
            row_c=C%ku+o_c+1
            d_a  =1
            d_b  =-o_b+1
            d_c  =-o_b+1
            do n_r=max(0,-o_a,o_b),max(0,nrmax+min(0,-o_a,o_b))-1
               C%dat(row_c+row_c_offset,n_r+d_c) = C%dat(row_c+row_c_offset,n_r+d_c)+ &
               &                                   A%dat(row_a,n_r+d_a)*B%dat(row_b,n_r+d_b)
            end do
         end do
      end do

   end subroutine band_real_band_real_product
!------------------------------------------------------------------------------
   subroutine band_complex_band_real_product(A,B,C,l_lhs)
      !
      ! This subroutine computes a matrix multiplication between
      ! two input banded matrices
      !

      !-- Input variables
      type(type_bandmat_complex), intent(in) :: A
      type(type_bandmat_real),    intent(in) :: B
      logical,                    intent(in) :: l_lhs

      !-- Output variables
      type(type_bandmat_complex), intent(inout) :: C

      !-- Local variables
      integer :: klC, kuC, nrmax, n_r
      integer :: o_a,o_c,o_b,d_a,d_b,d_c
      integer :: row_a,row_b,row_c,row_c_offset

      if ( .not. allocated(C%dat) ) then
         klC=A%kl+B%kl
         kuC=A%ku+B%ku
         nrmax =A%nlines
         call C%initialize(klC, kuC, nrmax, l_lhs)
      end if

      if ( l_lhs ) then
         row_c_offset = klC
      else
         row_c_offset = 0
      end if

      do o_c=-min(C%ku, A%ku+B%ku),min(C%kl,A%kl+B%kl)
         do o_a=-min(C%ku,B%kl-o_c),min(A%kl,B%ku+o_c)
            o_b=o_c-o_a
            row_a=A%ku+o_a+1
            row_b=B%ku+o_b+1
            row_c=C%ku+o_c+1
            d_a  =1
            d_b  =-o_b+1
            d_c  =-o_b+1
            do n_r=max(0,-o_a,o_b),max(0,nrmax+min(0,-o_a,o_b))-1
               C%dat(row_c+row_c_offset,n_r+d_c) = C%dat(row_c+row_c_offset,n_r+d_c)+ &
               &                                   A%dat(row_a,n_r+d_a)*B%dat(row_b,n_r+d_b)
            end do
         end do
      end do

   end subroutine band_complex_band_real_product
!------------------------------------------------------------------------------
end module matrix_types
