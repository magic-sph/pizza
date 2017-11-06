module algebra

   use precision_mod, only: cp

   use lapack95, only: getrs, getrf

   implicit none

   private

   interface prepare_full_mat
      module procedure :: prepare_full_mat_real
      module procedure :: prepare_full_mat_complex
   end interface prepare_full_mat

   interface solve_full_mat
      module procedure :: solve_full_mat_real_rhs_real
      module procedure :: solve_full_mat_real_rhs_complex
      module procedure :: solve_full_mat_complex_rhs_complex
   end interface solve_full_mat

   public :: prepare_full_mat, solve_full_mat

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
   subroutine prepapre_full_mat_complex(a,len_a,n,pivot,info)
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
end module algebra
