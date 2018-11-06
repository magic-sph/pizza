module bordered_matrix
   !
   ! This module defines several types for handling bordered matrices
   !

   use precision_mod
   use constants, only: one, zero
   use mem_alloc, only: bytes_allocated
   use algebra, only: prepare_bordered_mat, solve_bordered_mat

   implicit none

   private

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

contains

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

      open(file_handle, file='A_mat_tau', form='unformatted', access='stream')

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
end module bordered_matrix
