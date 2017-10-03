module chebsparselib
   !
   ! This module defines several useful recurrence relations when
   ! sparse Chebyshev representation is required. This relies on
   ! the python routines provided by P. Marti et al., GGG, 2016
   !

   use precision_mod
   use constants, only: half, one, two

   implicit none

   private

   public :: rmult1, intcheb1, intcheb2, intcheb4, intcheb2rmult1, eye, &
   &         rmult2, intcheb2rmult2, intcheb1rmult1, intcheb2rmult2lapl

contains

   function eye(len_stencil) result(stencil)
      !
      ! Identity operator (validation: OK)
      !

      !-- Inputvariable
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku)            =0.0_cp
      stencil(ku+1)            =one      ! Diagonal
      stencil(ku+2:len_stencil)=0.0_cp

   end function eye
!------------------------------------------------------------------------------
   function rmult1(a, b, len_stencil) result(stencil)
      !
      ! Multiplication by r operator (validation: OK)
      !

      !-- Inputvariable
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku-1)         =0.0_cp
      stencil(ku)             =half*a ! 1st upper diagonal
      stencil(ku+1)           =b      ! Diagonal
      stencil(ku+2)           =half*a ! 1st lower diagonal
      stencil(ku+3:len_stencil)=0.0_cp

   end function rmult1
!------------------------------------------------------------------------------
   function rmult2(a, b, len_stencil) result(stencil)
      !
      ! Multiplication by r operator
      !

      !-- Inputvariable
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku-2)          =0.0_cp
      stencil(ku-1)            =0.25_cp*a*a
      stencil(ku)              =a*b ! 1st upper diagonal
      stencil(ku+1)            =half*(a*a+two*b*b) ! Diagonal
      stencil(ku+2)            =a*b ! 1st lower diagonal
      stencil(ku+3)            =0.25_cp*a*a
      stencil(ku+4:len_stencil)=0.0_cp

   end function rmult2
!------------------------------------------------------------------------------
   function intcheb1(a, idx, len_stencil) result(stencil)
      !
      ! First integral operator (validation: OK)
      !

      !-- Inputvariable
      real(cp), intent(in) :: a
      integer,  intent(in) :: idx
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku-1)          = 0.0_cp
      stencil(ku)              =-half*a/real(idx,cp)  ! 1st upper diagonal
      stencil(ku+1)            = 0.0_cp               ! Diagonal
      stencil(ku+2)            = half*a/real(idx,cp) ! 1st lower diagonal
      stencil(ku+3:len_stencil)= 0.0_cp

   end function intcheb1
!------------------------------------------------------------------------------
   function intcheb2(a, idx, len_stencil) result(stencil)
      !
      ! Second integral operator (validation: OK)
      !

      !-- Inputvariable
      real(cp), intent(in) :: a
      integer,  intent(in) :: idx
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku-2)          = 0.0_cp
      stencil(ku-1)            = a*a/4.0_cp/real(idx,cp)/real(idx+1,cp)
      stencil(ku)              = 0.0_cp
      stencil(ku+1)            =-a*a/2.0_cp/real(idx-1,cp)/real(idx+1,cp)
      stencil(ku+2)            = 0.0_cp
      stencil(ku+3)            = a*a/4.0_cp/real(idx,cp)/real(idx-1,cp)
      stencil(ku+4:len_stencil)= 0.0_cp

   end function intcheb2
!------------------------------------------------------------------------------
   function intcheb4(a, idx, len_stencil) result(stencil)
      !
      ! Fourth-integral opertor (validation: OK)
      !

      !-- Inputvariable
      real(cp), intent(in) :: a
      integer,  intent(in) :: idx
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku-4)          = 0.0_cp
      stencil(ku-3)            = a*a*a*a/16.0_cp/real(idx,cp)/real(idx+1,cp)/&
      &                          real(idx+2,cp)/real(idx+3,cp)
      stencil(ku-2)            = 0.0_cp
      stencil(ku-1)            =-a*a*a*a/4.0_cp/real(idx,cp)/real(idx-1,cp)/ &
      &                          real(idx+1,cp)/real(idx+3,cp)
      stencil(ku)              = 0.0_cp
      stencil(ku+1)            = 3.0_cp*a*a*a*a/8.0_cp/real(idx-2,cp)/     &
      &                          real(idx-1,cp)/real(idx+1,cp)/real(idx+2,cp)
      stencil(ku+2)            = 0.0_cp
      stencil(ku+3)            =-a*a*a*a/4.0_cp/real(idx,cp)/real(idx-3,cp)/ &
      &                          real(idx-1,cp)/real(idx+1,cp)
      stencil(ku+4)            = 0.0_cp
      stencil(ku+5)            = a*a*a*a/16.0_cp/real(idx,cp)/real(idx-3,cp)/ &
      &                          real(idx-2,cp)/real(idx-1,cp)
      stencil(ku+6:len_stencil)= 0.0_cp

   end function intcheb4
!------------------------------------------------------------------------------
   function intcheb1rmult1(a, b, idx, len_stencil) result(stencil)
      !
      ! Second integral of r operator
      !

      !-- Inputvariable
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
      integer,  intent(in) :: idx
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku-2)          = 0.0_cp
      stencil(ku-1)            =-0.25_cp*a*a/real(idx,cp)
      stencil(ku)              =-half*a*b/real(idx,cp)
      stencil(ku+1)            = 0.0_cp
      stencil(ku+2)            = half*a*b/real(idx,cp)
      stencil(ku+3)            = 0.25_cp*a*a/real(idx,cp)
      stencil(ku+4:len_stencil)= 0.0_cp

   end function intcheb1rmult1
!------------------------------------------------------------------------------
   function intcheb2rmult1(a, b, idx, len_stencil) result(stencil)
      !
      ! Second integral of r operator (validation: OK)
      !

      !-- Inputvariable
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
      integer,  intent(in) :: idx
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku-3)          = 0.0_cp
      stencil(ku-2)            = a*a*a/8.0_cp/real(idx,cp)/real(idx+1,cp)
      stencil(ku-1)            = a*a*b/4.0_cp/real(idx,cp)/real(idx+1,cp)
      stencil(ku)              =-a*a*a/8.0_cp/real(idx,cp)/real(idx-1,cp)
      stencil(ku+1)            =-a*a*b/2.0_cp/real(idx-1,cp)/real(idx+1,cp)
      stencil(ku+2)            =-a*a*a/8.0_cp/real(idx,cp)/real(idx+1,cp)
      stencil(ku+3)            = a*a*b/4.0_cp/real(idx,cp)/real(idx-1,cp)
      stencil(ku+4)            = a*a*a/8.0_cp/real(idx,cp)/real(idx-1,cp)
      stencil(ku+5:len_stencil)= 0.0_cp

   end function intcheb2rmult1
!------------------------------------------------------------------------------
   function intcheb2rmult2(a, b, idx, len_stencil) result(stencil)
      !
      ! Second integral of r operator
      !

      !-- Inputvariable
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
      integer,  intent(in) :: idx
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku-4)          = 0.0_cp
      stencil(ku-3)            = a*a*a*a/16.0_cp/real(idx,cp)/real(idx+1,cp)
      stencil(ku-2)            = a*a*a*b/4.0_cp/real(idx,cp)/real(idx+1,cp)
      stencil(ku-1)            =-a*a*(a*a-2.0_cp*b*b*real(idx-1,cp))/  &
      &                          8.0_cp/real(idx,cp)/real(idx-1,cp)/   &
      &                          real(idx+1,cp)
      stencil(ku)              =-a*a*a*b/4.0_cp/real(idx,cp)/real(idx-1,cp)
      stencil(ku+1)            =-a*a*(a*a+4.0_cp*b*b)/8.0_cp/real(idx-1,cp)/&
      &                          real(idx+1,cp)
      stencil(ku+2)            =-a*a*a*b/4.0_cp/real(idx,cp)/real(idx+1,cp)
      stencil(ku+3)            = a*a*(a*a+2.0_cp*b*b*real(idx+1,cp))/  &
      &                          8.0_cp/real(idx,cp)/real(idx-1,cp)/   &
      &                          real(idx+1,cp)
      stencil(ku+4)            = a*a*a*b/4.0_cp/real(idx,cp)/real(idx-1,cp)
      stencil(ku+5)            = a*a*a*a/16.0_cp/real(idx,cp)/real(idx-1,cp)
      stencil(ku+6:len_stencil)= 0.0_cp

   end function intcheb2rmult2
!------------------------------------------------------------------------------
   function intcheb2rmult2lapl(a, b, m, idx, len_stencil) result(stencil)
      !
      ! Second integral of r**2*\Delta operator (same as the combination)
      !

      !-- Inputvariable
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
      integer,  intent(in) :: m
      integer,  intent(in) :: idx
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku-2)          = 0.0_cp
      stencil(ku-1)            =-a*a*real(m-idx-2,cp)*real(m+idx+2,cp)/ &
      &                          real(4*idx*(idx+1),cp)
      stencil(ku)              = half*a*b*real(2*idx+3,cp)/real(idx,cp)
      stencil(ku+1)            = (a*a*real(m*m,cp)+a*a*real(idx*idx,cp)-   &
      &                          two*a*a+two*b*b*real(idx*idx,cp)-two*b*b)/&
      &                          real(2*(idx-1)*(idx+1),cp)
      stencil(ku+2)            = a*b*real(2*idx-3,cp)/real(2*idx,cp)
      stencil(ku+3)            = -a*a*real(m-idx+2,cp)*real(m+idx-2,cp)/ &
      &                          real(4*idx*(idx-1),cp)
      stencil(ku+4:len_stencil)= 0.0_cp

   end function intcheb2rmult2lapl
!------------------------------------------------------------------------------
end module chebsparselib
