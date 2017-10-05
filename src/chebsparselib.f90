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

   public :: rmult1, intcheb1, intcheb2, intcheb4, intcheb2rmult1, eye,  &
   &         rmult2, intcheb2rmult2, intcheb1rmult1, intcheb2rmult2lapl, &
   &         intcheb4rmult4lapl2, intcheb4rmult4lapl, intcheb4rmult4

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
   function rmult1(a, b, idx, len_stencil) result(stencil)
      !
      ! Multiplication by r operator (validation: OK)
      !

      !-- Input variables
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
      integer,  intent(in) :: idx
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variable
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku-1)         =0.0_cp
      stencil(ku)             =half*a ! 1st upper diagonal
      stencil(ku+1)           =b      ! Diagonal
      stencil(ku+2)           =half*a ! 1st lower diagonal
      stencil(ku+3:len_stencil)=0.0_cp

      stencil = mirror_stencil(idx, len_stencil)

   end function rmult1
!------------------------------------------------------------------------------
   function rmult2(a, b, idx, len_stencil) result(stencil)
      !
      ! Multiplication by r^2 operator
      !

      !-- Input variables
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
      integer,  intent(in) :: idx
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variable
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku-2)          =0.0_cp
      stencil(ku-1)            =0.25_cp*a*a
      stencil(ku)              =a*b ! 1st upper diagonal
      stencil(ku+1)            =half*(a*a+two*b*b) ! Diagonal
      stencil(ku+2)            =a*b ! 1st lower diagonal
      stencil(ku+3)            =0.25_cp*a*a
      stencil(ku+4:len_stencil)=0.0_cp

      stencil = mirror_stencil(idx, len_stencil)

   end function rmult2
!------------------------------------------------------------------------------
   function intcheb1(a, idx, len_stencil) result(stencil)
      !
      ! First integral operator (validation: OK)
      !

      !-- Input variables
      real(cp), intent(in) :: a
      integer,  intent(in) :: idx
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variable
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku-1)          = 0.0_cp
      stencil(ku)              =-half*a/real(idx,cp)  ! 1st upper diagonal
      stencil(ku+1)            = 0.0_cp               ! Diagonal
      stencil(ku+2)            = half*a/real(idx,cp) ! 1st lower diagonal
      stencil(ku+3:len_stencil)= 0.0_cp

      stencil = mirror_stencil(idx, len_stencil)

   end function intcheb1
!------------------------------------------------------------------------------
   function intcheb2(a, idx, len_stencil) result(stencil)
      !
      ! Second integral operator (validation: OK)
      !

      !-- Input variables
      real(cp), intent(in) :: a
      integer,  intent(in) :: idx
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variable
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku-2)          = 0.0_cp
      stencil(ku-1)            = a*a/4.0_cp/real(idx,cp)/real(idx+1,cp)
      stencil(ku)              = 0.0_cp
      stencil(ku+1)            =-a*a/2.0_cp/real(idx-1,cp)/real(idx+1,cp)
      stencil(ku+2)            = 0.0_cp
      stencil(ku+3)            = a*a/4.0_cp/real(idx,cp)/real(idx-1,cp)
      stencil(ku+4:len_stencil)= 0.0_cp

      stencil = mirror_stencil(idx, len_stencil)

   end function intcheb2
!------------------------------------------------------------------------------
   function intcheb4(a, idx, len_stencil) result(stencil)
      !
      ! Fourth-integral opertor (validation: OK)
      !

      !-- Input variables
      real(cp), intent(in) :: a
      integer,  intent(in) :: idx
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variable
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

      stencil = mirror_stencil(idx, len_stencil)

   end function intcheb4
!------------------------------------------------------------------------------
   function intcheb1rmult1(a, b, idx, len_stencil) result(stencil)
      !
      ! Second integral of r operator
      !

      !-- Input variables
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

      stencil = mirror_stencil(idx, len_stencil)

   end function intcheb1rmult1
!------------------------------------------------------------------------------
   function intcheb2rmult1(a, b, idx, len_stencil) result(stencil)
      !
      ! Second integral of r operator (validation: OK)
      !

      !-- Input variables
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

      stencil = mirror_stencil(idx, len_stencil)

   end function intcheb2rmult1
!------------------------------------------------------------------------------
   function intcheb2rmult2(a, b, idx, len_stencil) result(stencil)
      !
      ! Second integral of r operator
      !

      !-- Input variables
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

      stencil = mirror_stencil(idx, len_stencil)

   end function intcheb2rmult2
!------------------------------------------------------------------------------
   function intcheb2rmult2lapl(a, b, m, idx, len_stencil) result(stencil)
      !
      ! Second integral of r**2*\Delta operator (same as the combination)
      !

      !-- Input variables
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

      stencil = mirror_stencil(idx, len_stencil)

   end function intcheb2rmult2lapl
!------------------------------------------------------------------------------
   function intcheb4rmult4(a, b, n, len_stencil) result(stencil)
      !
      ! Fourth integral of r^4 operator
      !
      !-- Input variables
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
      integer,  intent(in) :: n
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku-8)          = 0.0_cp
      stencil(ku-7)            = a**8/real(256*n*(n + 1)*(n + 2)*(n + 3),cp)
      stencil(ku-6)            = a**7*b/real(32*n*(n + 1)*(n + 2)*(n + 3),cp)
      stencil(ku-5)            =-3*a**6*(a**2 - 2*b**2*n + 2*b**2)/ &
      &                         real(64*n*(n - 1)*(n + 1)*(n + 2)*(n + 3),cp)
      stencil(ku-4)            =-a**5*b*(a**2*n + 11*a**2 - 4*b**2*n + 4*b**2)/&
      &                         real(32*n*(n - 1)*(n + 1)*(n + 2)*(n + 3),cp)
      stencil(ku-3)            =-a**4*(a**4*n**2 - 19*a**4 + 12*a**2*b**2*n**2 +&
      &                          36*a**2*b**2*n - 120*a**2*b**2 - 4*b**4*n**2 + &
      &                          12*b**4*n - 8*b**4)/                           &
      &                          real(64*n*(n-2)*(n-1)*(n+1)*(n+2)*(n+3),cp)
      stencil(ku-2)            =-3*a**5*b*(a**2*n - 6*a**2 + 4*b**2*n - 8*b**2)/&
      &                          real(32*n*(n - 2)*(n - 1)*(n + 1)*(n + 2),cp)
      stencil(ku-1)            = a**4*(9*a**4*n - 33*a**4 - 6*a**2*b**2*n**2 +   &
      &                          120*a**2*b**2*n - 306*a**2*b**2 - 16*b**4*n**2 +&
      &                          80*b**4*n - 96*b**4)/                           &
      &                          real(64*n*(n-3)*(n-2)*(n-1)*(n+1)*(n+3),cp)
      stencil(ku)              = a**5*b*(3*a**2*n**2 + 15*a**2*n - 102*a**2 +  &
      &                          8*b**2*n**2 + 40*b**2*n - 192*b**2)/          &
      &                          real(32*n*(n-3)*(n-2)*(n-1)*(n+2)*(n+3),cp)
      stencil(ku+1)            = 3*a**4*(a**4*n**2 - 29*a**4 + 16*a**2*b**2*n**2 & 
      &                          - 304*a**2*b**2 + 16*b**4*n**2 - 144*b**4)/     &
      &                          real(128*(n-3)*(n-2)*(n-1)*(n+1)*(n+2)*(n+3),cp)
      stencil(ku+2)            = a**5*b*(3*a**2*n**2 - 15*a**2*n - 102*a**2 +  &
      &                          8*b**2*n**2 - 40*b**2*n - 192*b**2)/          &
      &                          real(32*n*(n-3)*(n-2)*(n+1)*(n+2)*(n+3),cp)
      stencil(ku+3)            =-a**4*(9*a**4*n + 33*a**4 + 6*a**2*b**2*n**2 +   &
      &                          120*a**2*b**2*n + 306*a**2*b**2 + 16*b**4*n**2 +&
      &                           80*b**4*n + 96*b**4)/                          &
      &                          real(64*n*(n-3)*(n-1)*(n+1)*(n+2)*(n+3),cp)
      stencil(ku+4)            =-3*a**5*b*(a**2*n + 6*a**2 + 4*b**2*n + 8*b**2)/ &
      &                          real(32*n*(n - 2)*(n - 1)*(n + 1)*(n + 2),cp)
      stencil(ku+5)            =-a**4*(a**4*n**2 - 19*a**4 + 12*a**2*b**2*n**2 - &
      &                          36*a**2*b**2*n - 120*a**2*b**2 - 4*b**4*n**2 -  &
      &                          12*b**4*n - 8*b**4)/                            &
      &                          real(64*n*(n-3)*(n-2)*(n-1)*(n+1)*(n+2),cp)
      stencil(ku+6)            =-a**5*b*(a**2*n - 11*a**2 - 4*b**2*n - 4*b**2)/ &
      &                          real(32*n*(n - 3)*(n - 2)*(n - 1)*(n + 1),cp)
      stencil(ku+7)            = 3*a**6*(a**2 + 2*b**2*n + 2*b**2)/             &
      &                          real(64*n*(n - 3)*(n - 2)*(n - 1)*(n + 1),cp)
      stencil(ku+8)            = a**7*b/real(32*n*(n - 3)*(n - 2)*(n - 1),cp)
      stencil(ku+9)            = a**8/real(256*n*(n - 3)*(n - 2)*(n - 1),cp)

      stencil(ku+10:len_stencil)= 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4rmult4
!------------------------------------------------------------------------------
   function intcheb4rmult4lapl(a, b, m, n, len_stencil) result(stencil)
      !
      ! Fourth integral of r^2\Delta operator
      !

      !-- Input variables
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
      integer,  intent(in) :: m
      integer,  intent(in) :: n
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku-6)          = 0.0_cp
      stencil(ku-5)            =-a**6*real((m-n-6)*(m+n+6),cp)/ &
      &                          real(64*n*(n+1)*(n+2)*(n+3),cp)
      stencil(ku-4)            =-a**5*b*real(2*m**2-4*n**2-41*n-105,cp)/ &
      &                          real(32*n*(n+1)*(n+2)*(n+3),cp)
      stencil(ku-3)            = a**4*(a**2*m**2*n + 5*a**2*m**2 + a**2*n**3 + &
      &                          3*a**2*n**2 - 34*a**2*n - 120*a**2 -          &
      &                          2*b**2*m**2*n + 2*b**2*m**2 + 12*b**2*n**3 +  &
      &                          90*b**2*n**2 + 114*b**2*n - 216*b**2)/        &
      &                          real(32*n*(n - 1)*(n + 1)*(n + 2)*(n + 3),cp)
      stencil(ku-2)            =a**3*b*(6*a**2*m**2 + 4*a**2*n**2 - 25*a**2*n -&
      &                         183*a**2 + 16*b**2*n**2 + 44*b**2*n - 60*b**2)/&
      &                         real(32*n*(n - 1)*(n + 1)*(n + 2),cp)
      stencil(ku-1)            =a**2*(a**4*m**2*n - 17*a**4*m**2 - a**4*n**3 - &
      &                         27*a**4*n**2 - 8*a**4*n + 372*a**4 +           &
      &                         16*a**2*b**2*m**2*n - 32*a**2*b**2*m**2 -      &
      &                         216*a**2*b**2*n**2 - 360*a**2*b**2*n +         &
      &                         1584*a**2*b**2 + 16*b**4*n**3 - 112*b**4*n +   &
      &                         96*b**4)/                                      &
      &                         real(64*n*(n - 2)*(n - 1)*(n + 1)*(n + 3),cp)
      stencil(ku)              =-a**3*b*(2*a**2*m**2*n + 16*a**2*m**2 +        &
      &                          4*a**2*n**3 + 33*a**2*n**2 - 55*a**2*n -      &
      &                          444*a**2 + 8*b**2*n**3 + 66*b**2*n**2 +       &
      &                          10*b**2*n - 348*b**2)/                        &
      &                          real(16*n*(n - 2)*(n - 1)*(n + 2)*(n + 3))
      stencil(ku+1)            =-a**2*(a**4*m**2*n**2 - 19*a**4*m**2 +         &
      &                          a**4*n**4 - 43*a**4*n**2 + 396*a**4 +         &
      &                          6*a**2*b**2*m**2*n**2 - 54*a**2*b**2*m**2 +   &
      &                          12*a**2*b**2*n**4 - 336*a**2*b**2*n**2 +      &
      &                          2052*a**2*b**2 + 8*b**4*n**4 - 104*b**4*n**2 +&
      &                          288*b**4)/                                    &
      &                         real(16*(n-3)*(n-2)*(n-1)*(n+1)*(n+2)*(n+3),cp)
      stencil(ku+2)            =-a**3*b*(2*a**2*m**2*n - 16*a**2*m**2 +        &
      &                         4*a**2*n**3 - 33*a**2*n**2 - 55*a**2*n +       &
      &                         444*a**2 + 8*b**2*n**3 - 66*b**2*n**2 +        &
      &                         10*b**2*n + 348*b**2)/                         &
      &                         real(16*n*(n - 3)*(n - 2)*(n + 1)*(n + 2),cp)
      stencil(ku+3)            =a**2*(a**4*m**2*n + 17*a**4*m**2 - a**4*n**3 + &
      &                         27*a**4*n**2 - 8*a**4*n - 372*a**4 +           &
      &                         16*a**2*b**2*m**2*n + 32*a**2*b**2*m**2 +      &
      &                         216*a**2*b**2*n**2 - 360*a**2*b**2*n -         &
      &                         1584*a**2*b**2 + 16*b**4*n**3 - 112*b**4*n -   &
      &                         96*b**4)/                                      &
      &                         real(64*n*(n - 3)*(n - 1)*(n + 1)*(n + 2),cp)
      stencil(ku+4)            =a**3*b*(6*a**2*m**2 + 4*a**2*n**2 + 25*a**2*n -&
      &                         183*a**2 + 16*b**2*n**2 - 44*b**2*n - 60*b**2)/&
      &                         real(32*n*(n - 2)*(n - 1)*(n + 1),cp)
      stencil(ku+5)            =a**4*(a**2*m**2*n - 5*a**2*m**2 + a**2*n**3 - &
      &                         3*a**2*n**2 - 34*a**2*n + 120*a**2 -          &
      &                         2*b**2*m**2*n - 2*b**2*m**2 + 12*b**2*n**3 -  &
      &                         90*b**2*n**2 + 114*b**2*n + 216*b**2)/        &
      &                         real(32*n*(n - 3)*(n - 2)*(n - 1)*(n + 1),cp)
      stencil(ku+6)            =-a**5*b*real(2*m**2-4*n**2+41*n-105,cp)/      &
      &                         real(32*n*(n - 3)*(n - 2)*(n - 1),cp)
      stencil(ku+7)            =-a**6*real((m-n+6)*(m+n-6),cp)/               &
      &                         real(64*n*(n - 3)*(n - 2)*(n - 1),cp)
      stencil(ku+8:len_stencil)= 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4rmult4lapl
!------------------------------------------------------------------------------
   function intcheb4rmult4lapl2(a, b, m, n, len_stencil) result(stencil)
      !
      ! Fourth integral of r**4*\Delta\Delta operator
      !

      !-- Input variables
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
      integer,  intent(in) :: m
      integer,  intent(in) :: n
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      integer :: ku

      ku = (len_stencil-1)/2

      stencil(1:ku-4)          = 0.0_cp
      stencil(ku-3)            = a*a*a*a*real(m-n-6,cp)*real(m-n-4,cp)* &
      &                          real(m+n+4,cp)*real(m+n+6,cp)/         &
      &                          real(16*n*(n+1)*(n+2)*(n+3),cp)
      stencil(ku-2)            =-a*a*a*b*real(2*n+9,cp)*                &
      &                           real(2*m*m-2*n*n-18*n-39,cp)/         &
      &                           real(8*n*(n+1)*(n+2),cp)
      stencil(ku-1)            =-a*a*(a*a*m*m*m*m - 6*a*a*m*m*n - 28*a*a*m*m  &
      &                            -a*a*n*n*n*n-10*a*a*n*n*n-20*a*a*n*n +     &
      &                           64*a*a*n + 192*a*a + 2*b*b*m*m*n*n +        &
      &                           4*b*b*m*m*n-6*b*b*m*m - 6*b*b*n*n*n*n -     &
      &                           60*b*b*n*n*n-173*b*b*n*n-46*b*b*n+285*b*b)/ &
      &                           real(4*n*(n-1)*(n+1)*(n+3),cp)
      stencil(ku)              =a*b*(4*a*a*m*m*n + 38*a*a*m*m + 12*a*a*n*n*n +&
      &                         54*a*a*n*n-88*a*a*n-461*a*a+16*b*b*n*n*n +    &
      &                         72*b*b*n*n+24*b*b*n-112*b*b)/                 &
      &                         real(8*n*(n-1)*(n+2),cp)
      stencil(ku+1)            = (3*a*a*a*a*m*m*m*m + 2*a*a*a*a*m*m*n*n -     &
      &                         68*a*a*a*a*m*m+3*a*a*a*a*n*n*n*n -            &
      &                         68*a*a*a*a*n*n+416*a*a*a*a+8*a*a*b*b*m*m*n*n  &
      &                         -32*a*a*b*b*m*m+24*a*a*b*b*n*n*n*n-332*a*a*   &
      &                         b*b*n*n+944*a*a*b*b+8*b*b*b*b*n*n*n*n-        &
      &                         40*b*b*b*b*n*n+32*b*b*b*b)/                   &
      &                         real(8*(n-2)*(n-1)*(n+1)*(n+2),cp)

      stencil(ku+2)            =a*b*(4*a*a*m*m*n-38*a*a*m*m+12*a*a*n*n*n -    &
      &                         54*a*a*n*n-88*a*a*n+461*a*a+16*b*b*n*n*n -    &
      &                         72*b*b*n*n+24*b*b*n+112*b*b)/                 &
      &                         real(8*n*(n-2)*(n+1),cp)
      stencil(ku+3)            =-a*a*(a*a*m*m*m*m+6*a*a*m*m*n-28*a*a*m*m -    &
      &                          a*a*n*n*n*n+10*a*a*n*n*n-20*a*a*n*n -        &
      &                          64*a*a*n+192*a*a+2*b*b*m*m*n*n-4*b*b*m*m*n - &
      &                          6*b*b*m*m-6*b*b*n*n*n*n+60*b*b*n*n*n -       &
      &                          173*b*b*n*n+46*b*b*n+285*b*b)/               &
      &                          real(4*n*(n-3)*(n-1)*(n+1),cp)
      stencil(ku+4)            =-a*a*a*b*real((2*n-9)*(2*m*m-2*n*n+18*n-39),cp) &
      &                          /real(8*n*(n-2)*(n-1))
      stencil(ku+5)            =a*a*a*a*real((m-n+4)*(m-n+6)*(m+n-6)*(m+n-4),cp)&
      &                          /real(16*n*(n-3)*(n-2)*(n-1),cp)
      stencil(ku+6:len_stencil)= 0.0_cp


      ! stencil_debug(ku-3)=a**4*(m - n - 6)*(m - n - 4)*(m + n + 4)*(m + n + 6)/ &
      ! &                   (16*n*(n + 1)*(n + 2)*(n + 3))
      ! stencil_debug(ku-2)=-a**3*b*(2*n + 9)*(2*m**2 - 2*n**2 - 18*n - 39)/ &
      ! &                   (8*n*(n + 1)*(n + 2))
      ! stencil_debug(ku-1)=-a**2*(a**2*m**4 - 6*a**2*m**2*n - 28*a**2*m**2 - &
      ! &                    a**2*n**4 - 10*a**2*n**3 - 20*a**2*n**2 + 64*a**2*n + &
      ! &                    192*a**2 + 2*b**2*m**2*n**2 + 4*b**2*m**2*n - &
      ! &                    6*b**2*m**2 - 6*b**2*n**4 - 60*b**2*n**3 - &
      ! &                    173*b**2*n**2 - 46*b**2*n + 285*b**2)/ &
      ! &                    (4*n*(n - 1)*(n + 1)*(n + 3))
      ! stencil_debug(ku)  =a*b*(4*a**2*m**2*n + 38*a**2*m**2 + 12*a**2*n**3 + &
      ! &                   54*a**2*n**2 - 88*a**2*n - 461*a**2 + 16*b**2*n**3 +&
      ! &                   72*b**2*n**2 + 24*b**2*n - 112*b**2)/(8*n*(n - 1)*(n + 2))
      ! stencil_debug(ku+1)=(3*a**4*m**4 + 2*a**4*m**2*n**2 - 68*a**4*m**2 +  &
      ! &                    3*a**4*n**4 - 68*a**4*n**2 + 416*a**4 +          &
      ! &                    8*a**2*b**2*m**2*n**2 - 32*a**2*b**2*m**2 +      &
      ! &                    24*a**2*b**2*n**4 - 332*a**2*b**2*n**2 +         &
      ! &                    944*a**2*b**2 + 8*b**4*n**4 - 40*b**4*n**2 + 32*b**4)/&
      ! &                    (8*(n - 2)*(n - 1)*(n + 1)*(n + 2))
      ! stencil_debug(ku+2)=a*b*(4*a**2*m**2*n - 38*a**2*m**2 + 12*a**2*n**3 - &
      ! &                   54*a**2*n**2 - 88*a**2*n + 461*a**2 + 16*b**2*n**3 -&
      ! &                   72*b**2*n**2 + 24*b**2*n + 112*b**2)/(8*n*(n - 2)*(n + 1))
      ! stencil_debug(ku+3)=-a**2*(a**2*m**4 + 6*a**2*m**2*n - 28*a**2*m**2 - &
      ! &                    a**2*n**4 + 10*a**2*n**3 - 20*a**2*n**2 -        &
      ! &                    64*a**2*n + 192*a**2 + 2*b**2*m**2*n**2 -        &
      ! &                    4*b**2*m**2*n - 6*b**2*m**2 - 6*b**2*n**4 +      &
      ! &                    60*b**2*n**3 - 173*b**2*n**2 + 46*b**2*n + 285*b**2)/&
      ! &                    (4*n*(n - 3)*(n - 1)*(n + 1))
      ! stencil_debug(ku+4)=-a**3*b*(2*n - 9)*(2*m**2 - 2*n**2 + 18*n - 39)/&
      ! &                    (8*n*(n - 2)*(n - 1))
      ! stencil_debug(ku+5)=a**4*(m - n + 4)*(m - n + 6)*(m + n - 6)*(m + n - 4)/ &
      ! &                   (16*n*(n - 3)*(n - 2)*(n - 1))

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4rmult4lapl2
!------------------------------------------------------------------------------
   function mirror_stencil(idx, len_stencil) result(stencil)

      !-- Input variables
      integer,  intent(in) :: idx
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variable
      integer :: n, ku

      ku = (len_stencil-1)/2

      if ( idx < ku ) then
         !print*, 'idx=', idx
         do n=1,ku-idx
            !print*, len_stencil-ku+idx-n, '<-', len_stencil-ku+idx-n, '+', len_stencil-ku+idx+n
            stencil(len_stencil-ku+idx-n) = stencil(len_stencil-ku+idx-n)+ &
            &                               stencil(len_stencil-ku+idx+n)
         end do
      end if

   end function mirror_stencil
!------------------------------------------------------------------------------
end module chebsparselib
