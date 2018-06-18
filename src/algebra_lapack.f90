module algebra

   use precision_mod, only: cp
   use constants, only: one

   implicit none

   private

   interface prepare_band_mat
      module procedure prepare_band_mat_real
      module procedure prepare_band_mat_complex
   end interface prepare_band_mat

   interface prepare_bordered_mat
      module procedure prepare_bordered_mat_real
      module procedure prepare_bordered_mat_complex
   end interface prepare_bordered_mat

   interface prepare_full_mat
      module procedure prepare_full_mat_real
      module procedure prepare_full_mat_complex
   end interface prepare_full_mat

   interface solve_band_mat
      module procedure solve_band_mat_real
      module procedure solve_band_mat_complex
   end interface solve_band_mat

   interface solve_bordered_mat
      module procedure solve_bordered_mat_real
      module procedure solve_bordered_mat_complex
   end interface solve_bordered_mat

   interface solve_full_mat
      module procedure solve_full_mat_real_rhs_real
      module procedure solve_full_mat_real_rhs_complex
      module procedure solve_full_mat_complex_rhs_complex
   end interface solve_full_mat

   public :: prepare_full_mat, prepare_bordered_mat, solve_full_mat, &
   &         solve_bordered_mat, prepare_band_mat, solve_band_mat

contains

   subroutine solve_full_mat_real_rhs_complex(a,len_a,n,pivot,rhs)
      !
      !  This routine does the backward substitution into a lu-decomposed real 
      !  matrix a (to solve a * x = bc1) were bc1 is the right hand side  
      !  vector. On return x is stored in bc1.                            
      !                                                                     

      !-- Input variables:
      integer,  intent(in) :: n          ! dimension of problem
      integer,  intent(in) :: len_a      ! first dim of a
      integer,  intent(in) :: pivot(n)   ! pivot pointer of legth n
      real(cp), intent(in) :: a(len_a,n) ! real n X n matrix

      !-- Output variables
      complex(cp), intent(inout) :: rhs(n) ! on input RHS of problem
      real(cp) :: tmp_real(n), tmp_imag(n)
      integer :: info, n_r

#if (DEFAULT_PRECISION==sngl)
      call cgetrs('N',n,1,cmplx(a,0.0_cp,kind=cp),len_a,pivot,rhs,n,info)
#elif (DEFAULT_PRECISION==dble)
      ! call zgetrs('N',n,1,cmplx(a,0.0_cp,kind=cp),len_a,pivot,rhs,n,info)
      tmp_real(:)=real(rhs(:))
      tmp_imag(:)=aimag(rhs(:))
      call dgetrs('N',n,1,a,len_a,pivot,tmp_real,n,info)
      call dgetrs('N',n,1,a,len_a,pivot,tmp_imag,n,info)

      do n_r=1,n
         rhs(n_r)=cmplx(tmp_real(n_r), tmp_imag(n_r), kind=cp)
      end do
#endif

   end subroutine solve_full_mat_real_rhs_complex
!-----------------------------------------------------------------------------
   subroutine solve_full_mat_real_rhs_real(a,len_a,n,pivot,rhs)
      !
      !  This routine does the backward substitution into a lu-decomposed real 
      !  matrix a (to solve a * x = bc1) were bc1 is the right hand side  
      !  vector. On return x is stored in bc1.                            
      !                                                                     

      !-- Input variables:
      integer,  intent(in) :: n          ! dimension of problem
      integer,  intent(in) :: len_a      ! first dim of a
      integer,  intent(in) :: pivot(n)   ! pivot pointer of legth n
      real(cp), intent(in) :: a(len_a,n) ! real n X n matrix

      !-- Output variables
      real(cp), intent(inout) :: rhs(n) ! on input RHS of problem
      integer :: info

#if (DEFAULT_PRECISION==sngl)
      call sgetrs('N',n,1,a,len_a,pivot,rhs,n,info)
#elif (DEFAULT_PRECISION==dble)
      call dgetrs('N',n,1,a,len_a,pivot,rhs,n,info)
#endif

   end subroutine solve_full_mat_real_rhs_real
!-----------------------------------------------------------------------------
   subroutine solve_full_mat_complex_rhs_complex(a,len_a,n,pivot,rhs)
      !
      !  This routine does the backward substitution into a lu-decomposed real 
      !  matrix a (to solve a * x = bc1) were bc1 is the right hand side  
      !  vector. On return x is stored in bc1.                            
      !                                                                     

      !-- Input variables:
      integer,  intent(in) :: n          ! dimension of problem
      integer,  intent(in) :: len_a      ! first dim of a
      integer,  intent(in) :: pivot(n)   ! pivot pointer of legth n
      complex(cp), intent(in) :: a(len_a,n) ! real n X n matrix

      !-- Output variables
      complex(cp), intent(inout) :: rhs(n) ! on input RHS of problem
      integer :: info

#if (DEFAULT_PRECISION==sngl)
      call cgetrs('N',n,1,a,len_a,pivot,rhs,n,info)
#elif (DEFAULT_PRECISION==dble)
      call zgetrs('N',n,1,a,len_a,pivot,rhs,n,info)
#endif

   end subroutine solve_full_mat_complex_rhs_complex
!-----------------------------------------------------------------------------
   subroutine prepare_full_mat_real(a,len_a,n,pivot,info)
      !
      !     like the linpack routine
      !
      !     lu decomposes the real matrix a(n,n) via gaussian elimination
      !

      !-- Input variables:
      integer,  intent(in) :: len_a,n
      real(cp), intent(inout) :: a(len_a,n)

      !-- Output variables:
      integer, intent(out) :: pivot(n)   ! pivoting information
      integer, intent(out) :: info

#if (DEFAULT_PRECISION==sngl)
      call sgetrf(n,n,a(1:n,1:n),n,pivot(1:n),info)
#elif (DEFAULT_PRECISION==dble)
      call dgetrf(n,n,a,n,pivot,info)
#endif

   end subroutine prepare_full_mat_real
!-----------------------------------------------------------------------------
   subroutine prepare_full_mat_complex(a,len_a,n,pivot,info)
      !
      !     like the linpack routine
      !
      !     lu decomposes the real matrix a(n,n) via gaussian elimination
      !

      !-- Input variables:
      integer,  intent(in) :: len_a,n
      complex(cp), intent(inout) :: a(len_a,n)

      !-- Output variables:
      integer, intent(out) :: pivot(n)   ! pivoting information
      integer, intent(out) :: info

#if (DEFAULT_PRECISION==sngl)
      call cgetrf(n,n,a(1:n,1:n),n,pivot(1:n),info)
#elif (DEFAULT_PRECISION==dble)
      call zgetrf(n,n,a,n,pivot,info)
#endif

   end subroutine prepare_full_mat_complex
!-----------------------------------------------------------------------------
   subroutine prepare_band_mat_real(A,lenA,kl,ku,pivotA)

      !-- Input variables
      integer, intent(in) :: lenA
      integer, intent(in) :: kl
      integer, intent(in) :: ku

      !-- Output variables
      real(cp), intent(inout) :: A(2*kl+ku+1,lenA)
      integer,  intent(out)   :: pivotA(lenA)

      !-- Local variables
      integer :: n_bands, info

      n_bands = 2*kl+ku+1

      call dgbtrf(lenA, lenA, kl, ku, A, n_bands, pivotA, info)

   end subroutine prepare_band_mat_real
!-----------------------------------------------------------------------------
   subroutine prepare_band_mat_complex(A,lenA,kl,ku,pivotA)

      !-- Input variables
      integer, intent(in) :: lenA
      integer, intent(in) :: kl
      integer, intent(in) :: ku

      !-- Output variables
      complex(cp), intent(inout) :: A(2*kl+ku+1,lenA)
      integer,     intent(out)   :: pivotA(lenA)

      !-- Local variables
      integer :: n_bands, info

      n_bands = 2*kl+ku+1

      call zgbtrf(lenA, lenA, kl, ku, A, n_bands, pivotA, info)

   end subroutine prepare_band_mat_complex
!-----------------------------------------------------------------------------
   subroutine prepare_bordered_mat_real(A1,A2,A3,A4,n_boundaries,lenA4,kl,ku, &
              &                         pivotA1, pivotA4)

      !-- Input variables
      integer, intent(in) :: n_boundaries
      integer, intent(in) :: lenA4
      integer, intent(in) :: kl
      integer, intent(in) :: ku

      !-- Output variables
      real(cp), intent(inout) :: A1(n_boundaries,n_boundaries)
      real(cp), intent(inout) :: A2(n_boundaries,lenA4)
      real(cp), intent(inout) :: A3(lenA4,n_boundaries)
      real(cp), intent(inout) :: A4(2*kl+ku+1,lenA4)
      integer,  intent(out)   :: pivotA1(n_boundaries)
      integer,  intent(out)   :: pivotA4(lenA4)

      !-- Local variables
      integer :: n_bands_A4, info

      n_bands_A4 = 2*kl+ku+1

      !-- LU factorisation for the banded block
      call dgbtrf(lenA4, lenA4, kl, ku, A4, n_bands_A4, pivotA4, info)
      !-- Solve A4*v = A3 (on output v = A3)
      call dgbtrs('N', lenA4, kl, ku, n_boundaries, A4, n_bands_A4, pivotA4, &
           &      A3, lenA4, info)

      !-- Assemble the Schur complement of A4: A1 <- A1-A2*v
      call dgemm('N', 'N', n_boundaries, n_boundaries, lenA4, -one, A2,  &
           &     n_boundaries, A3, lenA4, one, A1,  n_boundaries)
      !-- LU factorisation of the Schur complement
      call dgetrf(n_boundaries, n_boundaries, A1, n_boundaries, pivotA1, info)

   end subroutine prepare_bordered_mat_real
!-----------------------------------------------------------------------------
   subroutine prepare_bordered_mat_complex(A1,A2,A3,A4,n_boundaries,lenA4,kl,ku, &
              &                            pivotA1, pivotA4)

      !-- Input variables
      integer, intent(in) :: n_boundaries
      integer, intent(in) :: lenA4
      integer, intent(in) :: kl
      integer, intent(in) :: ku

      !-- Output variables
      complex(cp), intent(inout) :: A1(n_boundaries,n_boundaries)
      complex(cp), intent(inout) :: A2(n_boundaries,lenA4)
      complex(cp), intent(inout) :: A3(lenA4,n_boundaries)
      complex(cp), intent(inout) :: A4(2*kl+ku+1,lenA4)
      integer,  intent(out)   :: pivotA1(n_boundaries)
      integer,  intent(out)   :: pivotA4(lenA4)

      !-- Local variables
      integer :: n_bands_A4, info

      n_bands_A4 = 2*kl+ku+1

      !-- LU factorisation for the banded block
      call zgbtrf(lenA4, lenA4, kl, ku, A4, n_bands_A4, pivotA4, info)
      !-- Solve A4*v = A3 (on output v = A3)
      call zgbtrs('N', lenA4, kl, ku, n_boundaries, A4, n_bands_A4, pivotA4, &
           &      A3, lenA4, info)

      !-- Assemble the Schur complement of A4: A1 <- A1-A2*v
      call zgemm('N', 'N', n_boundaries, n_boundaries, lenA4, -(one,0.0_cp), A2,  &
           &     n_boundaries, A3, lenA4, (one,0.0_cp), A1,  n_boundaries)
      !-- LU factorisation of the Schur complement
      call zgetrf(n_boundaries, n_boundaries, A1, n_boundaries, pivotA1, info)

   end subroutine prepare_bordered_mat_complex
!-----------------------------------------------------------------------------
   subroutine solve_band_mat_real(A, lenA, kl, ku, pivotA, rhs, lenRhs)

      !-- Input variables
      integer,  intent(in) :: kl
      integer,  intent(in) :: ku
      integer,  intent(in) :: lenA
      integer,  intent(in) :: lenRhs
      integer,  intent(in) :: pivotA(lenA)
      real(cp), intent(in) :: A(2*kl+ku+1,lenA)

      !-- Output variable
      real(cp), intent(inout) :: rhs(lenRhs)

      !-- Local variables:
      integer :: n_bands, info

      n_bands = 2*kl+ku+1

      call dgbtrs('N', lenA, kl, ku, 1, A, n_bands, pivotA, rhs(:), lenRhs, info)

   end subroutine solve_band_mat_real
!-----------------------------------------------------------------------------
   subroutine solve_band_mat_complex(A, lenA, kl, ku, pivotA, rhs, lenRhs)

      !-- Input variables
      integer,     intent(in) :: kl
      integer,     intent(in) :: ku
      integer,     intent(in) :: lenA
      integer,     intent(in) :: lenRhs
      integer,     intent(in) :: pivotA(lenA)
      complex(cp), intent(in) :: A(2*kl+ku+1,lenA)

      !-- Output variable
      complex(cp), intent(inout) :: rhs(lenRhs)

      !-- Local variables:
      integer :: n_bands, info

      n_bands = 2*kl+ku+1

      call zgbtrs('N', lenA, kl, ku, 1, A, n_bands, pivotA, rhs(:), lenA, info)

   end subroutine solve_band_mat_complex
!-----------------------------------------------------------------------------
   subroutine solve_bordered_mat_real(A1, A2, A3, A4, n_boundaries, lenA4, kl, &
              &                       ku, pivotA1, pivotA4, rhs, lenRhs)

      !-- Input variables
      integer,  intent(in) :: n_boundaries
      integer,  intent(in) :: kl
      integer,  intent(in) :: ku
      integer,  intent(in) :: lenA4
      integer,  intent(in) :: lenRhs
      integer,  intent(in) :: pivotA1(n_boundaries)
      integer,  intent(in) :: pivotA4(lenA4)
      real(cp), intent(in) :: A1(n_boundaries,n_boundaries)
      real(cp), intent(in) :: A2(n_boundaries,lenA4)
      real(cp), intent(in) :: A3(lenA4,n_boundaries)
      real(cp), intent(in) :: A4(2*kl+ku+1,lenA4)

      !-- Output variable
      real(cp), intent(inout) :: rhs(lenRhs)

      !-- Local variables:
      integer :: nStart, n_bands_A4, info

      nStart = n_boundaries+1
      n_bands_A4 = 2*kl+ku+1

      !-- Solve A4*w = rhs2
      call dgbtrs('N', lenA4, kl, ku, 1, A4, n_bands_A4, pivotA4, &
           &      rhs(nStart:), lenA4, info)

      !-- rhs1 <- rhs1-A2*rhs2
      call dgemv('N', n_boundaries, lenA4, -one, A2, n_boundaries,  &
           &     rhs(nStart:), 1, one, rhs(1:n_boundaries), 1)

      !-- Solve A1*y = rhs1
      call dgetrs('N', n_boundaries, 1, A1, n_boundaries, pivotA1, &
           &      rhs(1:n_boundaries), n_boundaries, info)

      !-- Assemble rhs2 <- rhs2-A3*rhs1
      call dgemv('N', lenA4, n_boundaries, -one, A3, lenA4, &
           &      rhs(1:n_boundaries), 1, one, rhs(nStart:), 1)

   end subroutine solve_bordered_mat_real
!-----------------------------------------------------------------------------
   subroutine solve_bordered_mat_complex(A1, A2, A3, A4, n_boundaries, lenA4, &
              &                          kl, ku, pivotA1, pivotA4, rhs, lenRhs)

      !-- Input variables
      integer,  intent(in) :: n_boundaries
      integer,  intent(in) :: kl
      integer,  intent(in) :: ku
      integer,  intent(in) :: lenA4
      integer,  intent(in) :: lenRhs
      integer,  intent(in) :: pivotA1(n_boundaries)
      integer,  intent(in) :: pivotA4(lenA4)
      complex(cp), intent(in) :: A1(n_boundaries,n_boundaries)
      complex(cp), intent(in) :: A2(n_boundaries,lenA4)
      complex(cp), intent(in) :: A3(lenA4,n_boundaries)
      complex(cp), intent(in) :: A4(2*kl+ku+1,lenA4)

      !-- Output variable
      complex(cp), intent(inout) :: rhs(lenRhs)

      !-- Local variables:
      integer :: nStart, n_bands_A4, info

      nStart = n_boundaries+1
      n_bands_A4 = 2*kl+ku+1

      !-- Solve A4*w = rhs2
      call zgbtrs('N', lenA4, kl, ku, 1, A4, n_bands_A4, pivotA4, &
           &      rhs(nStart:), lenA4, info)

      !-- rhs1 <- rhs1-A2*rhs2
      call zgemv('N', n_boundaries, lenA4, -(one,0.0_cp), A2, n_boundaries,  &
           &     rhs(nStart:), 1, (one,0.0_cp), rhs(1:n_boundaries), 1)

      !-- Solve A1*y = rhs1
      call zgetrs('N', n_boundaries, 1, A1, n_boundaries, pivotA1, &
           &      rhs(1:n_boundaries), n_boundaries, info)

      !-- Assemble rhs2 <- rhs2-A3*rhs1
      call zgemv('N', lenA4, n_boundaries, -(one,0.0_cp), A3, lenA4, &
           &      rhs(1:n_boundaries), 1, (one,0.0_cp), rhs(nStart:), 1)

   end subroutine solve_bordered_mat_complex
!-----------------------------------------------------------------------------
end module algebra
