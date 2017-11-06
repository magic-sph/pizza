module algebra

   use precision_mod, only: cp
   use constants, only: one
   use lapack95, only: getrs, getrf

   implicit none

   private

   interface prepare_bordered_mat
      module procedure prepare_bordered_mat_real
      module procedure prepare_bordered_mat_complex
   end interface prepare_bordered_mat

   interface prepare_full_mat
      module procedure :: prepare_full_mat_real
      module procedure :: prepare_full_mat_complex
   end interface prepare_full_mat

   interface solve_bordered_mat
      module procedure solve_bordered_mat_real
      module procedure solve_bordered_mat_complex
   end interface solve_bordered_mat

   interface solve_full_mat
      module procedure :: solve_full_mat_real_rhs_real
      module procedure :: solve_full_mat_real_rhs_complex
      module procedure :: solve_full_mat_complex_rhs_complex
   end interface solve_full_mat

   public :: prepare_full_mat, solve_full_mat, prepare_bordered_mat, &
   &         solve_bordered_mat

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

      call getrs(cmplx(a,0.0_cp,kind=cp),pivot,rhs)

   end subroutine solve_full_mat_real_rhs_complex
!-----------------------------------------------------------------------------
   subroutine solve_full_mat_complex_rhs_complex(a,len_a,n,pivot,rhs)
      !
      !  This routine does the backward substitution into a lu-decomposed real 
      !  matrix a (to solve a * x = bc1) were bc1 is the right hand side  
      !  vector. On return x is stored in bc1.                            
      !                                                                     

      !-- Input variables:
      integer,     intent(in) :: n          ! dimension of problem
      integer,     intent(in) :: len_a      ! first dim of a
      integer,     intent(in) :: pivot(n)   ! pivot pointer of legth n
      complex(cp), intent(in) :: a(len_a,n) ! real n X n matrix

      !-- Output variables
      complex(cp), intent(inout) :: rhs(n) ! on input RHS of problem

      call getrs(a,pivot,rhs)

   end subroutine solve_full_mat_complex_rhs_complex
!-----------------------------------------------------------------------------
   subroutine solve_full_mat_real_rhs_real(a,len_a,n,pivot,rhs)

      !-- Input variables:
      integer,  intent(in) :: n          ! dimension of problem
      integer,  intent(in) :: len_a      ! first dim of a
      integer,  intent(in) :: pivot(n)   ! pivot pointer of legth n
      real(cp), intent(in) :: a(len_a,n) ! real n X n matrix

      !-- Output variables
      real(cp), intent(inout) :: rhs(n) ! on input RHS of problem

      call getrs(a,pivot,rhs)

   end subroutine solve_full_mat_real_rhs_real
!-----------------------------------------------------------------------------
   subroutine prepare_full_mat_real(a,len_a,n,pivot,info)
      !
      !     lu decomposes the real matrix a(n,n) via gaussian elimination
      !

      !-- Input variables:
      integer,  intent(in) :: len_a,n
      real(cp), intent(inout) :: a(len_a,n)

      !-- Output variables:
      integer, intent(out) :: pivot(n)   ! pivoting information
      integer, intent(out) :: info

      call getrf(a(1:n,1:n),pivot(1:n),info)

   end subroutine prepare_full_mat_real
!-----------------------------------------------------------------------------
   subroutine prepare_full_mat_complex(a,len_a,n,pivot,info)
      !
      !     lu decomposes the real matrix a(n,n) via gaussian elimination
      !

      !-- Input variables:
      integer,     intent(in) :: len_a,n
      complex(cp), intent(inout) :: a(len_a,n)

      !-- Output variables:
      integer, intent(out) :: pivot(n)   ! pivoting information
      integer, intent(out) :: info

      call getrf(a(1:n,1:n),pivot(1:n),info)

   end subroutine prepare_full_mat_complex
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
      call getrf(A1, pivotA1, info)

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
      call getrf(A1, pivotA1, info)

   end subroutine prepare_bordered_mat_complex
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
      call getrs(A1, pivotA1, rhs(1:n_boundaries))

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
      call getrs(A1, pivotA1, rhs(1:n_boundaries))

      !-- Assemble rhs2 <- rhs2-A3*rhs1
      call zgemv('N', lenA4, n_boundaries, -(one,0.0_cp), A3, lenA4, &
           &      rhs(1:n_boundaries), 1, (one,0.0_cp), rhs(nStart:), 1)

   end subroutine solve_bordered_mat_complex
!-----------------------------------------------------------------------------
end module algebra
