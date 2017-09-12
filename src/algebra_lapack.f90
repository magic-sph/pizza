module algebra

   use precision_mod, only: cp

   implicit none

   private

   public :: cgefa, sgefa, sgesl, cgesl, rgesl

contains

   subroutine sgesl(a,len_a,n,pivot,rhs)
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
      !real(cp) :: tmp_real(n), tmp_imag(n)
      integer :: info!, n_r

#if (DEFAULT_PRECISION==sngl)
      call cgetrs('N',n,1,cmplx(a,0.0_cp,kind=cp),len_a,pivot,rhs,n,info)
#elif (DEFAULT_PRECISION==dble)
      call zgetrs('N',n,1,cmplx(a,0.0_cp,kind=cp),len_a,pivot,rhs,n,info)
      ! tmp_real(:)=real(rhs(:))
      ! tmp_imag(:)=aimag(rhs(:))
      ! call dgetrs('N',n,1,a,len_a,pivot,tmp_real,n,info)
      ! call dgetrs('N',n,1,a,len_a,pivot,tmp_imag,n,info)
! 
      ! do n_r=1,n
         ! rhs(n_r)=cmplx(tmp_real(n_r), tmp_imag(n_r), kind=cp)
      ! end do
#endif

   end subroutine sgesl
!-----------------------------------------------------------------------------
   subroutine rgesl(a,len_a,n,pivot,rhs)
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
      integer :: info, n_r

#if (DEFAULT_PRECISION==sngl)
      call sgetrs('N',n,1,a,len_a,pivot,rhs,n,info)
#elif (DEFAULT_PRECISION==dble)
      call dgetrs('N',n,1,a,len_a,pivot,rhs,n,info)
#endif

   end subroutine rgesl
!-----------------------------------------------------------------------------
   subroutine cgesl(a,len_a,n,pivot,rhs)
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
      real(cp) :: tmp_real(n), tmp_imag(n)
      integer :: info, n_r

#if (DEFAULT_PRECISION==sngl)
      call cgetrs('N',n,1,a,len_a,pivot,rhs,n,info)
#elif (DEFAULT_PRECISION==dble)
      call zgetrs('N',n,1,a,len_a,pivot,rhs,n,info)
#endif

   end subroutine cgesl
!-----------------------------------------------------------------------------
   subroutine sgefa(a,len_a,n,pivot,info)
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

   end subroutine sgefa
!-----------------------------------------------------------------------------
   subroutine cgefa(a,len_a,n,pivot,info)
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

   end subroutine cgefa
!-----------------------------------------------------------------------------
end module algebra
