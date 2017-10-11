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
   &         intcheb4rmult4lapl2, intcheb4rmult4lapl, intcheb4rmult4,    &
   &         intcheb4rmult4hmult8

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
      ! Multiplication by r^2 operator (validation: OK)
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
      &                          real(4*idx,cp)/real(idx+1,cp)
      stencil(ku)              = half*a*b*real(2*idx+3,cp)/real(idx,cp)
      stencil(ku+1)            = (a*a*real(m*m,cp)+a*a*real(idx*idx,cp)-   &
      &                          two*a*a+two*b*b*real(idx*idx,cp)-two*b*b)/&
      &                          real(2*(idx-1),cp)/real(idx+1,cp)
      stencil(ku+2)            = a*b*real(2*idx-3,cp)/real(2*idx,cp)
      stencil(ku+3)            = -a*a*real(m-idx+2,cp)*real(m+idx-2,cp)/ &
      &                          real(4*idx,cp)/real(idx-1,cp)
      stencil(ku+4:len_stencil)= 0.0_cp

      stencil = mirror_stencil(idx, len_stencil)

   end function intcheb2rmult2lapl
!------------------------------------------------------------------------------
   function intcheb4rmult4(a, b, n, len_stencil) result(stencil)
      !
      ! Fourth integral of r^4 operator (validation: OK)
      !
      !-- Input variables
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
      integer,  intent(in) :: n
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      real(cp) :: dn
      integer :: ku

      ku = (len_stencil-1)/2
      dn = real(n, cp)

      stencil(1:ku-8)          = 0.0_cp
      stencil(ku-7)            = a**8/(256*dn)/(dn+1)/(dn+2)/(dn+3)
      stencil(ku-6)            = a**7*b/(32*dn)/(dn+1)/(dn+2)/(dn+3)
      stencil(ku-5)            =-3*a**6*(a**2-2*b**2*(dn-1))/(64*dn)/(dn-1)/  &
      &                          (dn+1)/(dn+2)/(dn+3)
      stencil(ku-4)            =-a**5*b*(a**2*(dn+11)-4*b**2*(dn-1))/(32*dn)/ &
      &                          (dn-1)/(dn+1)/(dn+2)/(dn+3)
      stencil(ku-3)            =-a**4*(a**4*((dn)*(dn)-19)+12*a**2*b**2*(     &
      &                          dn*(dn+3)-10)-4*b**4*((dn)*(dn-3)+2))/       &
      &                          (64*dn)/(dn-2)/(dn-1)/(dn+1)/(dn+2)/(dn+3)
      stencil(ku-2)            =-3*a**5*b*(a**2*(dn-6)+4*b**2*(dn-2))/        &
      &                          (32*dn)/(dn-2)/(dn-1)/(dn+1)/(dn+2)
      stencil(ku-1)            = a**4*(a**4*(9*dn-33)-6*a**2*b**2*            &
      &                          (dn*(dn-20)+51)-16*b**4*(dn*(dn-5)+6))/      &
      &                          (64*dn)/(dn-3)/(dn-2)/(dn-1)/(dn+1)/(dn+3)
      stencil(ku)              = a**5*b*(a**2*(3*dn*(dn+5)-102)+8*b**2*(      &
      &                          dn*(dn+5)-24))/(32*dn)/(dn-3)/(dn-2)/        &
      &                          (dn-1)/(dn+2)/(dn+3)
      stencil(ku+1)            = 3*a**4*(a**4*(dn*dn-29)+16*a**2*b**2*(dn*dn-19)+ &
      &                          16*b**4*(dn*dn-9))/(128*(dn-3))/(dn-2)/(dn-1)    &
      &                          /(dn+1)/(dn+2)/(dn+3)
      stencil(ku+2)            = a**5*b*(a**2*((3*dn)*(dn-5)-102)+8*b**2*(   &
      &                          dn*(dn-5)-24))/(32*dn)/(dn-3)/(dn-2)/       &
      &                          (dn+1)/(dn+2)/(dn+3)
      stencil(ku+3)            =-a**4*(a**4*(9*dn+33)+6*a**2*b**2*(dn*(dn+20)+51)+&
      &                          16*b**4*(dn*(dn+5)+6))/(64*dn)/(dn-3)/(dn-1)/    &
      &                          (dn+1)/(dn+2)/(dn+3)
      stencil(ku+4)            =-3*a**5*b*(a**2*(dn+6)+4*b**2*(dn+2))/       &
      &                          (32*dn)/(dn-2)/(dn-1)/(dn+1)/(dn+2)
      stencil(ku+5)            =-a**4*(a**4*(dn*dn-19)+12*a**2*b**2*(dn*(dn-3)-10)&
      &                          -4*b**4*(dn*(dn+3)+2))/(64*dn)/(dn-3)/(dn-2)/    &
      &                          (dn-1)/(dn+1)/(dn+2)
      stencil(ku+6)            =-a**5*b*(a**2*(dn-11)-4*b**2*(dn+1))/   &
      &                          (32*dn)/(dn-3)/(dn-2)/(dn-1)/(dn+1)
      stencil(ku+7)            = 3*a**6*(a**2+2*b**2*(dn+1))/           &
      &                          (64*dn)/(dn-3)/(dn-2)/(dn-1)/(dn+1)
      stencil(ku+8)            = a**7*b/(32*dn)/(dn-3)/(dn-2)/(dn-1)
      stencil(ku+9)            = a**8/(256*dn)/(dn-3)/(dn-2)/(dn-1)

      stencil(ku+10:len_stencil)= 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4rmult4
!------------------------------------------------------------------------------
   function intcheb4rmult4lapl(a, b, m, n, len_stencil) result(stencil)
      !
      ! Fourth integral of r^4\Delta operator
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
      real(cp) :: dm, dn
      integer :: ku

      dm = real(m, cp)
      dn = real(n, cp)

      ku = (len_stencil-1)/2

      stencil(1:ku-6)          = 0.0_cp
      stencil(ku-5)            =-a**6*((dm-dn-6)*(dm+dn+6))/ &
      &                          (64*dn*(dn+1)*(dn+2)*(dn+3))
      stencil(ku-4)            =-a**5*b*(2*dm**2-4*dn**2-41*dn-105)/ &
      &                          (32*dn*(dn+1)*(dn+2)*(dn+3))
      stencil(ku-3)            = a**4*(a**2*dm**2*dn + 5*a**2*dm**2 + a**2*dn**3 + &
      &                          3*a**2*dn**2 - 34*a**2*dn - 120*a**2 -            &
      &                          2*b**2*dm**2*dn + 2*b**2*dm**2 + 12*b**2*dn**3 +  &
      &                          90*b**2*dn**2 + 114*b**2*dn - 216*b**2)/          &
      &                          (32*dn*(dn - 1)*(dn + 1)*(dn + 2)*(dn + 3))
      stencil(ku-2)            =a**3*b*(6*a**2*dm**2 + 4*a**2*dn**2 - 25*a**2*dn -&
      &                         183*a**2 + 16*b**2*dn**2 + 44*b**2*dn - 60*b**2)/ &
      &                         (32*dn*(dn - 1)*(dn + 1)*(dn + 2))
      stencil(ku-1)            =a**2*(a**4*dm**2*dn - 17*a**4*dm**2 - a**4*dn**3 - &
      &                         27*a**4*dn**2 - 8*a**4*dn + 372*a**4 +             &
      &                         16*a**2*b**2*dm**2*dn - 32*a**2*b**2*dm**2 -       &
      &                         216*a**2*b**2*dn**2 - 360*a**2*b**2*dn +           &
      &                         1584*a**2*b**2 + 16*b**4*dn**3 - 112*b**4*dn +     &
      &                         96*b**4)/                                          &
      &                         (64*dn*(dn - 2)*(dn - 1)*(dn + 1)*(dn + 3))
      stencil(ku)              =-a**3*b*(2*a**2*dm**2*dn + 16*a**2*dm**2 +     &
      &                          4*a**2*dn**3 + 33*a**2*dn**2 - 55*a**2*dn -   &
      &                          444*a**2 + 8*b**2*dn**3 + 66*b**2*dn**2 +     &
      &                          10*b**2*dn - 348*b**2)/                       &
      &                          (16*dn*(dn - 2)*(dn - 1)*(dn + 2)*(dn + 3))
      stencil(ku+1)            =-a**2*(a**4*dm**2*dn**2 - 19*a**4*dm**2 +       &
      &                          a**4*dn**4 - 43*a**4*dn**2 + 396*a**4 +        &
      &                          6*a**2*b**2*dm**2*dn**2 - 54*a**2*b**2*dm**2 + &
      &                          12*a**2*b**2*dn**4 - 336*a**2*b**2*dn**2 +     &
      &                          2052*a**2*b**2 + 8*b**4*dn**4 - 104*b**4*dn**2+&
      &                          288*b**4)/                                     &
      &                         (16*(dn-3)*(dn-2)*(dn-1)*(dn+1)*(dn+2)*(dn+3))
      stencil(ku+2)            =-a**3*b*(2*a**2*dm**2*dn - 16*a**2*dm**2 +      &
      &                         4*a**2*dn**3 - 33*a**2*dn**2 - 55*a**2*dn +     &
      &                         444*a**2 + 8*b**2*dn**3 - 66*b**2*dn**2 +       &
      &                         10*b**2*dn + 348*b**2)/                         &
      &                         (16*dn*(dn - 3)*(dn - 2)*(dn + 1)*(dn + 2))
      stencil(ku+3)            =a**2*(a**4*dm**2*dn + 17*a**4*dm**2 - a**4*dn**3 +&
      &                         27*a**4*dn**2 - 8*a**4*dn - 372*a**4 +            &
      &                         16*a**2*b**2*dm**2*dn + 32*a**2*b**2*dm**2 +      &
      &                         216*a**2*b**2*dn**2 - 360*a**2*b**2*dn -          &
      &                         1584*a**2*b**2 + 16*b**4*dn**3 - 112*b**4*dn -    &
      &                         96*b**4)/                                         &
      &                         (64*dn*(dn - 3)*(dn - 1)*(dn + 1)*(dn + 2))
      stencil(ku+4)            =a**3*b*(6*a**2*dm**2 + 4*a**2*dn**2 + 25*a**2*dn -&
      &                         183*a**2 + 16*b**2*dn**2 - 44*b**2*dn - 60*b**2)/ &
      &                         (32*dn*(dn - 2)*(dn - 1)*(dn + 1))
      stencil(ku+5)            =a**4*(a**2*dm**2*dn - 5*a**2*dm**2 + a**2*dn**3 - &
      &                         3*a**2*dn**2 - 34*a**2*dn + 120*a**2 -            &
      &                         2*b**2*dm**2*dn - 2*b**2*dm**2 + 12*b**2*dn**3 -  &
      &                         90*b**2*dn**2 + 114*b**2*dn + 216*b**2)/          &
      &                         (32*dn*(dn - 3)*(dn - 2)*(dn - 1)*(dn + 1))
      stencil(ku+6)            =-a**5*b*(2*dm**2-4*dn**2+41*dn-105)/      &
      &                         (32*dn*(dn - 3)*(dn - 2)*(dn - 1))
      stencil(ku+7)            =-a**6*((dm-dn+6)*(dm+dn-6))/           &
      &                         (64*dn*(dn - 3)*(dn - 2)*(dn - 1))
      stencil(ku+8:len_stencil)= 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4rmult4lapl
!------------------------------------------------------------------------------
   function intcheb4rmult4lapl2(a, b, m, n, len_stencil) result(stencil)
      !
      ! Fourth integral of r^4*\Delta\Delta operator
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
      real(cp) :: dm, dn
      integer :: ku

      dm = real(m,cp)
      dn = real(n,cp)

      ku = (len_stencil-1)/2

      stencil(1:ku-4)          = 0.0_cp
      stencil(ku-3)            = a*a*a*a*(dm-dn-6)*(dm-dn-4)*(dm+dn+4)*(dm+dn+6)/ &
      &                          (16*dn*(dn+1)*(dn+2)*(dn+3))
      stencil(ku-2)            =-a*a*a*b*(2*dn+9)*(2*dm*dm-2*dn*dn-18*dn-39)/  &
      &                          (8*dn*(dn+1)*(dn+2))
      stencil(ku-1)            =-a*a*(a*a*dm*dm*dm*dm-6*a*a*dm*dm*dn-28*a*a*dm*dm &
      &                          -a*a*dn*dn*dn*dn-10*a*a*dn*dn*dn-20*a*a*dn*dn +  &
      &                          64*a*a*dn + 192*a*a + 2*b*b*dm*dm*dn*dn +        &
      &                          4*b*b*dm*dm*dn-6*b*b*dm*dm - 6*b*b*dn*dn*dn*dn - &
      &                          60*b*b*dn*dn*dn-173*b*b*dn*dn-46*b*b*dn+285*b*b)/&
      &                          (4*dn*(dn-1)*(dn+1)*(dn+3))
      stencil(ku)              = a*b*(4*a*a*dm*dm*dn+38*a*a*dm*dm+12*a*a*dn*dn*dn +&
      &                          54*a*a*dn*dn-88*a*a*dn-461*a*a+16*b*b*dn*dn*dn +  &
      &                          72*b*b*dn*dn+24*b*b*dn-112*b*b)/                  &
      &                          (8*dn*(dn-1)*(dn+2))
      stencil(ku+1)            = (3*a*a*a*a*dm*dm*dm*dm + 2*a*a*a*a*dm*dm*dn*dn -  &
      &                          68*a*a*a*a*dm*dm+3*a*a*a*a*dn*dn*dn*dn -          &
      &                          68*a*a*a*a*dn*dn+416*a*a*a*a+8*a*a*b*b*dm*dm*dn*dn&
      &                          -32*a*a*b*b*dm*dm+24*a*a*b*b*dn*dn*dn*dn-332*a*a* &
      &                          b*b*dn*dn+944*a*a*b*b+8*b*b*b*b*dn*dn*dn*dn-      &
      &                          40*b*b*b*b*dn*dn+32*b*b*b*b)/                     &
      &                          (8*(dn-2)*(dn-1)*(dn+1)*(dn+2))

      stencil(ku+2)            = a*b*(4*a*a*dm*dm*dn-38*a*a*dm*dm+12*a*a*dn*dn*dn -&
      &                          54*a*a*dn*dn-88*a*a*dn+461*a*a+16*b*b*dn*dn*dn -  &
      &                          72*b*b*dn*dn+24*b*b*dn+112*b*b)/                  &
      &                          (8*dn*(dn-2)*(dn+1))
      stencil(ku+3)            =-a*a*(a*a*dm*dm*dm*dm+6*a*a*dm*dm*dn-28*a*a*dm*dm- &
      &                          a*a*dn*dn*dn*dn+10*a*a*dn*dn*dn-20*a*a*dn*dn -    &
      &                          64*a*a*dn+192*a*a+2*b*b*dm*dm*dn*dn-              &
      &                          4*b*b*dm*dm*dn-6*b*b*dm*dm-6*b*b*dn*dn*dn*dn+     &
      &                          60*b*b*dn*dn*dn-173*b*b*dn*dn+46*b*b*dn+285*b*b)/ &
      &                          (4*dn*(dn-3)*(dn-1)*(dn+1))
      stencil(ku+4)            =-a*a*a*b*((2*dn-9)*(2*dm*dm-2*dn*dn+18*dn-39)) &
      &                          /(8*dn*(dn-2)*(dn-1))
      stencil(ku+5)            = a*a*a*a*((dm-dn+4)*(dm-dn+6)*(dm+dn-6)*(dm+dn-4))&
      &                          /(16*dn*(dn-3)*(dn-2)*(dn-1))
      stencil(ku+6:len_stencil)= 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4rmult4lapl2
!------------------------------------------------------------------------------
   function intcheb4rmult4hmult8(a, b, r, n, len_stencil) result(stencil)

      !-- Input variables
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
      real(cp), intent(in) :: r 
      integer,  intent(in) :: n
      integer,  intent(in) :: len_stencil

      !-- Output variable 
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      real(cp) :: dn
      integer :: ku

      dn = real(n, cp)
      ku = (len_stencil-1)/2

      stencil(1:ku-16) = 0.0_cp
      stencil(ku-15)=a**16/(65536.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku-14)=3.0_cp*a**15*b/(8192.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp) &
      &              *(dn+3.0_cp))
      stencil(ku-13)=-a**14*(-2.0_cp*a**2*dn+5.0_cp*a**2-66.0_cp*b**2*dn+66.0_cp &
      &              *b**2+4.0_cp*dn*r**2-4.0_cp*r**2)/(16384.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-12)=-a**13*b*(-21.0_cp*a**2*dn+57.0_cp*a**2-220.0_cp*b** &
      &              2*dn+220.0_cp*b**2+40.0_cp*dn*r**2-40.0_cp*r**2)/(8192.0_cp*dn*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-11)=3.0_cp*a**12*(2.0_cp*a**4*dn**2-15.0_cp*a**4*dn+27.0_cp &
      &              *a**4+132.0_cp*a**2*b**2*dn**2-660.0_cp*a**2*b**2*dn+792.0_cp*a**2* &
      &              b**2-8.0_cp*a**2*dn**2*r**2+40.0_cp*a**2*dn*r**2-48.0_cp*a**2*r**2+660.0_cp &
      &              *b**4*dn**2-1980.0_cp*b**4*dn+1320.0_cp*b**4-240.0_cp*b**2*dn**2*r**2+ &
      &              720.0_cp*b**2*dn*r**2-480.0_cp*b**2*r**2+8.0_cp*dn**2*r**4-24.0_cp*dn &
      &              *r**4+16.0_cp*r**4)/(16384.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-10)=a**11*b*(51.0_cp*a**4*dn**2-441.0_cp*a**4*dn+858.0_cp &
      &              *a**4+1100.0_cp*a**2*b**2*dn**2-5940.0_cp*a**2*b**2*dn+7480.0_cp*a**2* &
      &              b**2-200.0_cp*a**2*dn**2*r**2+1080.0_cp*a**2*dn*r**2-1360.0_cp*a**2 &
      &              *r**2+3168.0_cp*b**4*dn**2-9504.0_cp*b**4*dn+6336.0_cp*b**4-1920.0_cp &
      &              *b**2*dn**2*r**2+5760.0_cp*b**2*dn*r**2-3840.0_cp*b**2*r**2+192.0_cp &
      &              *dn**2*r**4-576.0_cp*dn*r**4+384.0_cp*r**4)/(8192.0_cp*dn*(dn-2.0_cp &
      &              )*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-9)=-a**10*(-6.0_cp*a**6*dn**3+135.0_cp*a**6*dn**2-711.0_cp &
      &              *a**6*dn+1110.0_cp*a**6-726.0_cp*a**4*b**2*dn**3+9900.0_cp*a**4*b** &
      &              2*dn**2-39666.0_cp*a**4*b**2*dn+49500.0_cp*a**4*b**2+44.0_cp*a**4*dn**3 &
      &              *r**2-600.0_cp*a**4*dn**2*r**2+2404.0_cp*a**4*dn*r**2-3000.0_cp*a** &
      &              4*r**2-7920.0_cp*a**2*b**4*dn**3+71280.0_cp*a**2*b**4*dn**2-205920.0_cp &
      &              *a**2*b**4*dn+190080.0_cp*a**2*b**4+2880.0_cp*a**2*b**2*dn**3*r**2-25920.0_cp &
      &              *a**2*b**2*dn**2*r**2+74880.0_cp*a**2*b**2*dn*r**2-69120.0_cp*a**2* &
      &              b**2*r**2-96.0_cp*a**2*dn**3*r**4+864.0_cp*a**2*dn**2*r**4-2496.0_cp &
      &              *a**2*dn*r**4+2304.0_cp*a**2*r**4-14784.0_cp*b**6*dn**3+88704.0_cp* &
      &              b**6*dn**2-162624.0_cp*b**6*dn+88704.0_cp*b**6+13440.0_cp*b**4*dn** &
      &              3*r**2-80640.0_cp*b**4*dn**2*r**2+147840.0_cp*b**4*dn*r**2-80640.0_cp &
      &              *b**4*r**2-2688.0_cp*b**2*dn**3*r**4+16128.0_cp*b**2*dn**2*r**4-29568.0_cp &
      &              *b**2*dn*r**4+16128.0_cp*b**2*r**4+64.0_cp*dn**3*r**6-384.0_cp*dn** &
      &              2*r**6+704.0_cp*dn*r**6-384.0_cp*r**6)/(16384.0_cp*dn*(dn-3.0_cp)*(dn &
      &              -2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-8)=-3.0_cp*a**9*b*(-7.0_cp*a**6*dn**3+342.0_cp*a**6*dn** &
      &              2-2117.0_cp*a**6*dn+3582.0_cp*a**6-440.0_cp*a**4*b**2*dn**3+7920.0_cp &
      &              *a**4*b**2*dn**2-35640.0_cp*a**4*b**2*dn+47520.0_cp*a**4*b**2+80.0_cp &
      &              *a**4*dn**3*r**2-1440.0_cp*a**4*dn**2*r**2+6480.0_cp*a**4*dn*r**2-8640.0_cp &
      &              *a**4*r**2-3168.0_cp*a**2*b**4*dn**3+31680.0_cp*a**2*b**4*dn**2-98208.0_cp &
      &              *a**2*b**4*dn+95040.0_cp*a**2*b**4+1920.0_cp*a**2*b**2*dn**3*r**2-19200.0_cp &
      &              *a**2*b**2*dn**2*r**2+59520.0_cp*a**2*b**2*dn*r**2-57600.0_cp*a**2* &
      &              b**2*r**2-192.0_cp*a**2*dn**3*r**4+1920.0_cp*a**2*dn**2*r**4-5952.0_cp &
      &              *a**2*dn*r**4+5760.0_cp*a**2*r**4-4224.0_cp*b**6*dn**3+25344.0_cp*b**6 &
      &              *dn**2-46464.0_cp*b**6*dn+25344.0_cp*b**6+5376.0_cp*b**4*dn**3*r**2 &
      &              -32256.0_cp*b**4*dn**2*r**2+59136.0_cp*b**4*dn*r**2-32256.0_cp*b**4 &
      &              *r**2-1792.0_cp*b**2*dn**3*r**4+10752.0_cp*b**2*dn**2*r**4-19712.0_cp &
      &              *b**2*dn*r**4+10752.0_cp*b**2*r**4+128.0_cp*dn**3*r**6-768.0_cp*dn**2* &
      &              r**6+1408.0_cp*dn*r**6-768.0_cp*r**6)/(8192.0_cp*dn*(dn-3.0_cp)*(dn &
      &              -2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-7)=a**8*(-9.0_cp*a**8*dn**3-117.0_cp*a**8*dn**2+1401.0_cp &
      &              *a**8*dn-3237.0_cp*a**8-264.0_cp*a**6*b**2*dn**3-12672.0_cp*a**6*b**2* &
      &              dn**2+100056.0_cp*a**6*b**2*dn-186912.0_cp*a**6*b**2+16.0_cp*a**6*dn**3 &
      &              *r**2+768.0_cp*a**6*dn**2*r**2-6064.0_cp*a**6*dn*r**2+11328.0_cp*a**6* &
      &              r**2+3960.0_cp*a**4*b**4*dn**3-142560.0_cp*a**4*b**4*dn**2+756360.0_cp &
      &              *a**4*b**4*dn-1092960.0_cp*a**4*b**4-1440.0_cp*a**4*b**2*dn**3*r**2 &
      &              +51840.0_cp*a**4*b**2*dn**2*r**2-275040.0_cp*a**4*b**2*dn*r**2+397440.0_cp &
      &              *a**4*b**2*r**2+48.0_cp*a**4*dn**3*r**4-1728.0_cp*a**4*dn**2*r**4+9168.0_cp &
      &              *a**4*dn*r**4-13248.0_cp*a**4*r**4+29568.0_cp*a**2*b**6*dn**3-354816.0_cp &
      &              *a**2*b**6*dn**2+1212288.0_cp*a**2*b**6*dn-1241856.0_cp*a**2*b**6-26880.0_cp &
      &              *a**2*b**4*dn**3*r**2+322560.0_cp*a**2*b**4*dn**2*r**2-1102080.0_cp &
      &              *a**2*b**4*dn*r**2+1128960.0_cp*a**2*b**4*r**2+5376.0_cp*a**2*b**2*dn**3 &
      &              *r**4-64512.0_cp*a**2*b**2*dn**2*r**4+220416.0_cp*a**2*b**2*dn*r**4 &
      &              -225792.0_cp*a**2*b**2*r**4-128.0_cp*a**2*dn**3*r**6+1536.0_cp*a**2 &
      &              *dn**2*r**6-5248.0_cp*a**2*dn*r**6+5376.0_cp*a**2*r**6+31680.0_cp*b**8 &
      &              *dn**3-190080.0_cp*b**8*dn**2+348480.0_cp*b**8*dn-190080.0_cp*b**8-53760.0_cp &
      &              *b**6*dn**3*r**2+322560.0_cp*b**6*dn**2*r**2-591360.0_cp*b**6*dn*r**2+ &
      &              322560.0_cp*b**6*r**2+26880.0_cp*b**4*dn**3*r**4-161280.0_cp*b**4*dn**2 &
      &              *r**4+295680.0_cp*b**4*dn*r**4-161280.0_cp*b**4*r**4-3840.0_cp*b**2 &
      &              *dn**3*r**6+23040.0_cp*b**2*dn**2*r**6-42240.0_cp*b**2*dn*r**6+23040.0_cp &
      &              *b**2*r**6+64.0_cp*dn**3*r**8-384.0_cp*dn**2*r**8+704.0_cp*dn*r**8-384.0_cp &
      &              *r**8)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-6)=a**7*b*(-129.0_cp*a**8*dn**3-378.0_cp*a**8*dn**2+10461.0_cp &
      &              *a**8*dn-28098.0_cp*a**8-2200.0_cp*a**6*b**2*dn**3-18480.0_cp*a**6* &
      &              b**2*dn**2+226600.0_cp*a**6*b**2*dn-480480.0_cp*a**6*b**2+400.0_cp* &
      &              a**6*dn**3*r**2+3360.0_cp*a**6*dn**2*r**2-41200.0_cp*a**6*dn*r**2+87360.0_cp &
      &              *a**6*r**2-3168.0_cp*a**4*b**4*dn**3-133056.0_cp*a**4*b**4*dn**2+915552.0_cp &
      &              *a**4*b**4*dn-1463616.0_cp*a**4*b**4+1920.0_cp*a**4*b**2*dn**3*r**2 &
      &              +80640.0_cp*a**4*b**2*dn**2*r**2-554880.0_cp*a**4*b**2*dn*r**2+887040.0_cp &
      &              *a**4*b**2*r**2-192.0_cp*a**4*dn**3*r**4-8064.0_cp*a**4*dn**2*r**4+55488.0_cp &
      &              *a**4*dn*r**4-88704.0_cp*a**4*r**4+12672.0_cp*a**2*b**6*dn**3-228096.0_cp &
      &              *a**2*b**6*dn**2+899712.0_cp*a**2*b**6*dn-988416.0_cp*a**2*b**6-16128.0_cp &
      &              *a**2*b**4*dn**3*r**2+290304.0_cp*a**2*b**4*dn**2*r**2-1145088.0_cp &
      &              *a**2*b**4*dn*r**2+1257984.0_cp*a**2*b**4*r**2+5376.0_cp*a**2*b**2*dn**3 &
      &              *r**4-96768.0_cp*a**2*b**2*dn**2*r**4+381696.0_cp*a**2*b**2*dn*r**4 &
      &              -419328.0_cp*a**2*b**2*r**4-384.0_cp*a**2*dn**3*r**6+6912.0_cp*a**2 &
      &              *dn**2*r**6-27264.0_cp*a**2*dn*r**6+29952.0_cp*a**2*r**6+14080.0_cp &
      &              *b**8*dn**3-84480.0_cp*b**8*dn**2+154880.0_cp*b**8*dn-84480.0_cp*b**8- &
      &              30720.0_cp*b**6*dn**3*r**2+184320.0_cp*b**6*dn**2*r**2-337920.0_cp* &
      &              b**6*dn*r**2+184320.0_cp*b**6*r**2+21504.0_cp*b**4*dn**3*r**4-129024.0_cp &
      &              *b**4*dn**2*r**4+236544.0_cp*b**4*dn*r**4-129024.0_cp*b**4*r**4-5120.0_cp &
      &              *b**2*dn**3*r**6+30720.0_cp*b**2*dn**2*r**6-56320.0_cp*b**2*dn*r**6 &
      &              +30720.0_cp*b**2*r**6+256.0_cp*dn**3*r**8-1536.0_cp*dn**2*r**8+2816.0_cp &
      &              *dn*r**8-1536.0_cp*r**8)/(8192.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-5)=-3.0_cp*a**6*(10.0_cp*a**10*dn**3-39.0_cp*a**10*dn**2 &
      &              -495.0_cp*a**10*dn+2106.0_cp*a**10+858.0_cp*a**8*b**2*dn**3-1452.0_cp &
      &              *a**8*b**2*dn**2-43362.0_cp*a**8*b**2*dn+143748.0_cp*a**8*b**2-52.0_cp &
      &              *a**8*dn**3*r**2+88.0_cp*a**8*dn**2*r**2+2628.0_cp*a**8*dn*r**2-8712.0_cp &
      &              *a**8*r**2+7920.0_cp*a**6*b**4*dn**3+7920.0_cp*a**6*b**4*dn**2-427680.0_cp &
      &              *a**6*b**4*dn+1077120.0_cp*a**6*b**4-2880.0_cp*a**6*b**2*dn**3*r**2 &
      &              -2880.0_cp*a**6*b**2*dn**2*r**2+155520.0_cp*a**6*b**2*dn*r**2-391680.0_cp &
      &              *a**6*b**2*r**2+96.0_cp*a**6*dn**3*r**4+96.0_cp*a**6*dn**2*r**4-5184.0_cp &
      &              *a**6*dn*r**4+13056.0_cp*a**6*r**4+14784.0_cp*a**4*b**6*dn**3+88704.0_cp &
      &              *a**4*b**6*dn**2-1020096.0_cp*a**4*b**6*dn+1862784.0_cp*a**4*b**6-13440.0_cp &
      &              *a**4*b**4*dn**3*r**2-80640.0_cp*a**4*b**4*dn**2*r**2+927360.0_cp*a**4 &
      &              *b**4*dn*r**2-1693440.0_cp*a**4*b**4*r**2+2688.0_cp*a**4*b**2*dn**3 &
      &              *r**4+16128.0_cp*a**4*b**2*dn**2*r**4-185472.0_cp*a**4*b**2*dn*r**4 &
      &              +338688.0_cp*a**4*b**2*r**4-64.0_cp*a**4*dn**3*r**6-384.0_cp*a**4*dn**2 &
      &              *r**6+4416.0_cp*a**4*dn*r**6-8064.0_cp*a**4*r**6+126720.0_cp*a**2*b**8 &
      &              *dn**2-633600.0_cp*a**2*b**8*dn+760320.0_cp*a**2*b**8-215040.0_cp*a**2 &
      &              *b**6*dn**2*r**2+1075200.0_cp*a**2*b**6*dn*r**2-1290240.0_cp*a**2*b**6 &
      &              *r**2+107520.0_cp*a**2*b**4*dn**2*r**4-537600.0_cp*a**2*b**4*dn*r** &
      &              4+645120.0_cp*a**2*b**4*r**4-15360.0_cp*a**2*b**2*dn**2*r**6+76800.0_cp &
      &              *a**2*b**2*dn*r**6-92160.0_cp*a**2*b**2*r**6+256.0_cp*a**2*dn**2*r**8- &
      &              1280.0_cp*a**2*dn*r**8+1536.0_cp*a**2*r**8-5632.0_cp*b**10*dn**3+33792.0_cp &
      &              *b**10*dn**2-61952.0_cp*b**10*dn+33792.0_cp*b**10+15360.0_cp*b**8*dn**3 &
      &              *r**2-92160.0_cp*b**8*dn**2*r**2+168960.0_cp*b**8*dn*r**2-92160.0_cp &
      &              *b**8*r**2-14336.0_cp*b**6*dn**3*r**4+86016.0_cp*b**6*dn**2*r**4-157696.0_cp &
      &              *b**6*dn*r**4+86016.0_cp*b**6*r**4+5120.0_cp*b**4*dn**3*r**6-30720.0_cp &
      &              *b**4*dn**2*r**6+56320.0_cp*b**4*dn*r**6-30720.0_cp*b**4*r**6-512.0_cp &
      &              *b**2*dn**3*r**8+3072.0_cp*b**2*dn**2*r**8-5632.0_cp*b**2*dn*r**8+3072.0_cp &
      &              *b**2*r**8)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-4)=-a**5*b*(231.0_cp*a**10*dn**3-1782.0_cp*a**10*dn**2-7359.0_cp &
      &              *a**10*dn+47718.0_cp*a**10+6380.0_cp*a**8*b**2*dn**3-33000.0_cp*a** &
      &              8*b**2*dn**2-207020.0_cp*a**8*b**2*dn+957000.0_cp*a**8*b**2-1160.0_cp &
      &              *a**8*dn**3*r**2+6000.0_cp*a**8*dn**2*r**2+37640.0_cp*a**8*dn*r**2-174000.0_cp &
      &              *a**8*r**2+34848.0_cp*a**6*b**4*dn**3-95040.0_cp*a**6*b**4*dn**2-1137312.0_cp &
      &              *a**6*b**4*dn+3706560.0_cp*a**6*b**4-21120.0_cp*a**6*b**2*dn**3*r** &
      &              2+57600.0_cp*a**6*b**2*dn**2*r**2+689280.0_cp*a**6*b**2*dn*r**2-2246400.0_cp &
      &              *a**6*b**2*r**2+2112.0_cp*a**6*dn**3*r**4-5760.0_cp*a**6*dn**2*r**4 &
      &              -68928.0_cp*a**6*dn*r**4+224640.0_cp*a**6*r**4+50688.0_cp*a**4*b**6 &
      &              *dn**3-1723392.0_cp*a**4*b**6*dn+3801600.0_cp*a**4*b**6-64512.0_cp* &
      &              a**4*b**4*dn**3*r**2+2193408.0_cp*a**4*b**4*dn*r**2-4838400.0_cp*a**4* &
      &              b**4*r**2+21504.0_cp*a**4*b**2*dn**3*r**4-731136.0_cp*a**4*b**2*dn* &
      &              r**4+1612800.0_cp*a**4*b**2*r**4-1536.0_cp*a**4*dn**3*r**6+52224.0_cp &
      &              *a**4*dn*r**6-115200.0_cp*a**4*r**6+14080.0_cp*a**2*b**8*dn**3+84480.0_cp &
      &              *a**2*b**8*dn**2-689920.0_cp*a**2*b**8*dn+929280.0_cp*a**2*b**8-30720.0_cp &
      &              *a**2*b**6*dn**3*r**2-184320.0_cp*a**2*b**6*dn**2*r**2+1505280.0_cp &
      &              *a**2*b**6*dn*r**2-2027520.0_cp*a**2*b**6*r**2+21504.0_cp*a**2*b**4 &
      &              *dn**3*r**4+129024.0_cp*a**2*b**4*dn**2*r**4-1053696.0_cp*a**2*b**4 &
      &              *dn*r**4+1419264.0_cp*a**2*b**4*r**4-5120.0_cp*a**2*b**2*dn**3*r**6 &
      &              -30720.0_cp*a**2*b**2*dn**2*r**6+250880.0_cp*a**2*b**2*dn*r**6-337920.0_cp &
      &              *a**2*b**2*r**6+256.0_cp*a**2*dn**3*r**8+1536.0_cp*a**2*dn**2*r**8-12544.0_cp &
      &              *a**2*dn*r**8+16896.0_cp*a**2*r**8-3072.0_cp*b**10*dn**3+18432.0_cp &
      &              *b**10*dn**2-33792.0_cp*b**10*dn+18432.0_cp*b**10+10240.0_cp*b**8*dn**3 &
      &              *r**2-61440.0_cp*b**8*dn**2*r**2+112640.0_cp*b**8*dn*r**2-61440.0_cp &
      &              *b**8*r**2-12288.0_cp*b**6*dn**3*r**4+73728.0_cp*b**6*dn**2*r**4-135168.0_cp &
      &              *b**6*dn*r**4+73728.0_cp*b**6*r**4+6144.0_cp*b**4*dn**3*r**6-36864.0_cp &
      &              *b**4*dn**2*r**6+67584.0_cp*b**4*dn*r**6-36864.0_cp*b**4*r**6-1024.0_cp &
      &              *b**2*dn**3*r**8+6144.0_cp*b**2*dn**2*r**8-11264.0_cp*b**2*dn*r**8+6144.0_cp &
      &              *b**2*r**8)/(8192.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-3)=-a**4*(22.0_cp*a**12*dn**3-429.0_cp*a**12*dn**2-88.0_cp &
      &              *a**12*dn+8151.0_cp*a**12+2508.0_cp*a**10*b**2*dn**3-34848.0_cp*a**10 &
      &              *b**2*dn**2-31812.0_cp*a**10*b**2*dn+618552.0_cp*a**10*b**2-152.0_cp &
      &              *a**10*dn**3*r**2+2112.0_cp*a**10*dn**2*r**2+1928.0_cp*a**10*dn*r** &
      &              2-37488.0_cp*a**10*r**2+33660.0_cp*a**8*b**4*dn**3-320760.0_cp*a**8 &
      &              *b**4*dn**2-580140.0_cp*a**8*b**4*dn+5381640.0_cp*a**8*b**4-12240.0_cp &
      &              *a**8*b**2*dn**3*r**2+116640.0_cp*a**8*b**2*dn**2*r**2+210960.0_cp* &
      &              a**8*b**2*dn*r**2-1956960.0_cp*a**8*b**2*r**2+408.0_cp*a**8*dn**3*r**4 &
      &              -3888.0_cp*a**8*dn**2*r**4-7032.0_cp*a**8*dn*r**4+65232.0_cp*a**8*r**4 &
      &              +118272.0_cp*a**6*b**6*dn**3-709632.0_cp*a**6*b**6*dn**2-2247168.0_cp &
      &              *a**6*b**6*dn+11708928.0_cp*a**6*b**6-107520.0_cp*a**6*b**4*dn**3*r**2 &
      &              +645120.0_cp*a**6*b**4*dn**2*r**2+2042880.0_cp*a**6*b**4*dn*r**2-10644480.0_cp &
      &              *a**6*b**4*r**2+21504.0_cp*a**6*b**2*dn**3*r**4-129024.0_cp*a**6*b**2* &
      &              dn**2*r**4-408576.0_cp*a**6*b**2*dn*r**4+2128896.0_cp*a**6*b**2*r** &
      &              4-512.0_cp*a**6*dn**3*r**6+3072.0_cp*a**6*dn**2*r**6+9728.0_cp*a**6 &
      &              *dn*r**6-50688.0_cp*a**6*r**6+126720.0_cp*a**4*b**8*dn**3-380160.0_cp &
      &              *a**4*b**8*dn**2-2407680.0_cp*a**4*b**8*dn+7223040.0_cp*a**4*b**8-215040.0_cp &
      &              *a**4*b**6*dn**3*r**2+645120.0_cp*a**4*b**6*dn**2*r**2+4085760.0_cp &
      &              *a**4*b**6*dn*r**2-12257280.0_cp*a**4*b**6*r**2+107520.0_cp*a**4*b**4* &
      &              dn**3*r**4-322560.0_cp*a**4*b**4*dn**2*r**4-2042880.0_cp*a**4*b**4*dn &
      &              *r**4+6128640.0_cp*a**4*b**4*r**4-15360.0_cp*a**4*b**2*dn**3*r**6+46080.0_cp &
      &              *a**4*b**2*dn**2*r**6+291840.0_cp*a**4*b**2*dn*r**6-875520.0_cp*a** &
      &              4*b**2*r**6+256.0_cp*a**4*dn**3*r**8-768.0_cp*a**4*dn**2*r**8-4864.0_cp &
      &              *a**4*dn*r**8+14592.0_cp*a**4*r**8+33792.0_cp*a**2*b**10*dn**3-642048.0_cp &
      &              *a**2*b**10*dn+1013760.0_cp*a**2*b**10-92160.0_cp*a**2*b**8*dn**3*r**2 &
      &              +1751040.0_cp*a**2*b**8*dn*r**2-2764800.0_cp*a**2*b**8*r**2+86016.0_cp &
      &              *a**2*b**6*dn**3*r**4-1634304.0_cp*a**2*b**6*dn*r**4+2580480.0_cp*a**2 &
      &              *b**6*r**4-30720.0_cp*a**2*b**4*dn**3*r**6+583680.0_cp*a**2*b**4*dn &
      &              *r**6-921600.0_cp*a**2*b**4*r**6+3072.0_cp*a**2*b**2*dn**3*r**8-58368.0_cp &
      &              *a**2*b**2*dn*r**8+92160.0_cp*a**2*b**2*r**8-1024.0_cp*b**12*dn**3+6144.0_cp &
      &              *b**12*dn**2-11264.0_cp*b**12*dn+6144.0_cp*b**12+4096.0_cp*b**10*dn**3 &
      &              *r**2-24576.0_cp*b**10*dn**2*r**2+45056.0_cp*b**10*dn*r**2-24576.0_cp &
      &              *b**10*r**2-6144.0_cp*b**8*dn**3*r**4+36864.0_cp*b**8*dn**2*r**4-67584.0_cp &
      &              *b**8*dn*r**4+36864.0_cp*b**8*r**4+4096.0_cp*b**6*dn**3*r**6-24576.0_cp &
      &              *b**6*dn**2*r**6+45056.0_cp*b**6*dn*r**6-24576.0_cp*b**6*r**6-1024.0_cp &
      &              *b**4*dn**3*r**8+6144.0_cp*b**4*dn**2*r**8-11264.0_cp*b**4*dn*r**8+6144.0_cp &
      &              *b**4*r**8)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-2)=-3.0_cp*a**5*b*(11.0_cp*a**10*dn**2-1155.0_cp*a**10*dn &
      &              +5566.0_cp*a**10+660.0_cp*a**8*b**2*dn**2-29700.0_cp*a**8*b**2*dn+122760.0_cp &
      &              *a**8*b**2-120.0_cp*a**8*dn**2*r**2+5400.0_cp*a**8*dn*r**2-22320.0_cp &
      &              *a**8*r**2+6336.0_cp*a**6*b**4*dn**2-158400.0_cp*a**6*b**4*dn+544896.0_cp &
      &              *a**6*b**4-3840.0_cp*a**6*b**2*dn**2*r**2+96000.0_cp*a**6*b**2*dn*r**2 &
      &              -330240.0_cp*a**6*b**2*r**2+384.0_cp*a**6*dn**2*r**4-9600.0_cp*a**6 &
      &              *dn*r**4+33024.0_cp*a**6*r**4+16896.0_cp*a**4*b**6*dn**2-253440.0_cp &
      &              *a**4*b**6*dn+692736.0_cp*a**4*b**6-21504.0_cp*a**4*b**4*dn**2*r**2 &
      &              +322560.0_cp*a**4*b**4*dn*r**2-881664.0_cp*a**4*b**4*r**2+7168.0_cp &
      &              *a**4*b**2*dn**2*r**4-107520.0_cp*a**4*b**2*dn*r**4+293888.0_cp*a** &
      &              4*b**2*r**4-512.0_cp*a**4*dn**2*r**6+7680.0_cp*a**4*dn*r**6-20992.0_cp &
      &              *a**4*r**6+14080.0_cp*a**2*b**8*dn**2-126720.0_cp*a**2*b**8*dn+253440.0_cp &
      &              *a**2*b**8-30720.0_cp*a**2*b**6*dn**2*r**2+276480.0_cp*a**2*b**6*dn &
      &              *r**2-552960.0_cp*a**2*b**6*r**2+21504.0_cp*a**2*b**4*dn**2*r**4-193536.0_cp &
      &              *a**2*b**4*dn*r**4+387072.0_cp*a**2*b**4*r**4-5120.0_cp*a**2*b**2*dn**2 &
      &              *r**6+46080.0_cp*a**2*b**2*dn*r**6-92160.0_cp*a**2*b**2*r**6+256.0_cp &
      &              *a**2*dn**2*r**8-2304.0_cp*a**2*dn*r**8+4608.0_cp*a**2*r**8+3072.0_cp &
      &              *b**10*dn**2-15360.0_cp*b**10*dn+18432.0_cp*b**10-10240.0_cp*b**8*dn**2 &
      &              *r**2+51200.0_cp*b**8*dn*r**2-61440.0_cp*b**8*r**2+12288.0_cp*b**6*dn**2 &
      &              *r**4-61440.0_cp*b**6*dn*r**4+73728.0_cp*b**6*r**4-6144.0_cp*b**4*dn**2 &
      &              *r**6+30720.0_cp*b**4*dn*r**6-36864.0_cp*b**4*r**6+1024.0_cp*b**2*dn**2 &
      &              *r**8-5120.0_cp*b**2*dn*r**8+6144.0_cp*b**2*r**8)/(8192.0_cp*dn*(dn &
      &              -3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku-1)=a**4*(22.0_cp*a**12*dn**2+385.0_cp*a**12*dn-3003.0_cp &
      &              *a**12+1782.0_cp*a**10*b**2*dn**2+35640.0_cp*a**10*b**2*dn-241758.0_cp &
      &              *a**10*b**2-108.0_cp*a**10*dn**2*r**2-2160.0_cp*a**10*dn*r**2+14652.0_cp &
      &              *a**10*r**2+15840.0_cp*a**8*b**4*dn**2+396000.0_cp*a**8*b**4*dn-2280960.0_cp &
      &              *a**8*b**4-5760.0_cp*a**8*b**2*dn**2*r**2-144000.0_cp*a**8*b**2*dn* &
      &              r**2+829440.0_cp*a**8*b**2*r**2+192.0_cp*a**8*dn**2*r**4+4800.0_cp* &
      &              a**8*dn*r**4-27648.0_cp*a**8*r**4+29568.0_cp*a**6*b**6*dn**2+1182720.0_cp &
      &              *a**6*b**6*dn-5588352.0_cp*a**6*b**6-26880.0_cp*a**6*b**4*dn**2*r** &
      &              2-1075200.0_cp*a**6*b**4*dn*r**2+5080320.0_cp*a**6*b**4*r**2+5376.0_cp &
      &              *a**6*b**2*dn**2*r**4+215040.0_cp*a**6*b**2*dn*r**4-1016064.0_cp*a**6* &
      &              b**2*r**4-128.0_cp*a**6*dn**2*r**6-5120.0_cp*a**6*dn*r**6+24192.0_cp &
      &              *a**6*r**6+1140480.0_cp*a**4*b**8*dn-4181760.0_cp*a**4*b**8-1935360.0_cp &
      &              *a**4*b**6*dn*r**2+7096320.0_cp*a**4*b**6*r**2+967680.0_cp*a**4*b** &
      &              4*dn*r**4-3548160.0_cp*a**4*b**4*r**4-138240.0_cp*a**4*b**2*dn*r**6 &
      &              +506880.0_cp*a**4*b**2*r**6+2304.0_cp*a**4*dn*r**8-8448.0_cp*a**4*r**8 &
      &              -16896.0_cp*a**2*b**10*dn**2+337920.0_cp*a**2*b**10*dn-861696.0_cp* &
      &              a**2*b**10+46080.0_cp*a**2*b**8*dn**2*r**2-921600.0_cp*a**2*b**8*dn &
      &              *r**2+2350080.0_cp*a**2*b**8*r**2-43008.0_cp*a**2*b**6*dn**2*r**4+860160.0_cp &
      &              *a**2*b**6*dn*r**4-2193408.0_cp*a**2*b**6*r**4+15360.0_cp*a**2*b**4 &
      &              *dn**2*r**6-307200.0_cp*a**2*b**4*dn*r**6+783360.0_cp*a**2*b**4*r** &
      &              6-1536.0_cp*a**2*b**2*dn**2*r**8+30720.0_cp*a**2*b**2*dn*r**8-78336.0_cp &
      &              *a**2*b**2*r**8-4096.0_cp*b**12*dn**2+20480.0_cp*b**12*dn-24576.0_cp &
      &              *b**12+16384.0_cp*b**10*dn**2*r**2-81920.0_cp*b**10*dn*r**2+98304.0_cp &
      &              *b**10*r**2-24576.0_cp*b**8*dn**2*r**4+122880.0_cp*b**8*dn*r**4-147456.0_cp &
      &              *b**8*r**4+16384.0_cp*b**6*dn**2*r**6-81920.0_cp*b**6*dn*r**6+98304.0_cp &
      &              *b**6*r**6-4096.0_cp*b**4*dn**2*r**8+20480.0_cp*b**4*dn*r**8-24576.0_cp &
      &              *b**4*r**8)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+3.0_cp))
      stencil(ku)=a**5*b*(297.0_cp*a**10*dn**2+1485.0_cp*a**10*dn-21978.0_cp &
      &              *a**10+7920.0_cp*a**8*b**2*dn**2+39600.0_cp*a**8*b**2*dn-506880.0_cp &
      &              *a**8*b**2-1440.0_cp*a**8*dn**2*r**2-7200.0_cp*a**8*dn*r**2+92160.0_cp &
      &              *a**8*r**2+44352.0_cp*a**6*b**4*dn**2+221760.0_cp*a**6*b**4*dn-2395008.0_cp &
      &              *a**6*b**4-26880.0_cp*a**6*b**2*dn**2*r**2-134400.0_cp*a**6*b**2*dn &
      &              *r**2+1451520.0_cp*a**6*b**2*r**2+2688.0_cp*a**6*dn**2*r**4+13440.0_cp &
      &              *a**6*dn*r**4-145152.0_cp*a**6*r**4+76032.0_cp*a**4*b**6*dn**2+380160.0_cp &
      &              *a**4*b**6*dn-3345408.0_cp*a**4*b**6-96768.0_cp*a**4*b**4*dn**2*r** &
      &              2-483840.0_cp*a**4*b**4*dn*r**2+4257792.0_cp*a**4*b**4*r**2+32256.0_cp &
      &              *a**4*b**2*dn**2*r**4+161280.0_cp*a**4*b**2*dn*r**4-1419264.0_cp*a**4* &
      &              b**2*r**4-2304.0_cp*a**4*dn**2*r**6-11520.0_cp*a**4*dn*r**6+101376.0_cp &
      &              *a**4*r**6+42240.0_cp*a**2*b**8*dn**2+211200.0_cp*a**2*b**8*dn-1436160.0_cp &
      &              *a**2*b**8-92160.0_cp*a**2*b**6*dn**2*r**2-460800.0_cp*a**2*b**6*dn &
      &              *r**2+3133440.0_cp*a**2*b**6*r**2+64512.0_cp*a**2*b**4*dn**2*r**4+322560.0_cp &
      &              *a**2*b**4*dn*r**4-2193408.0_cp*a**2*b**4*r**4-15360.0_cp*a**2*b**2 &
      &              *dn**2*r**6-76800.0_cp*a**2*b**2*dn*r**6+522240.0_cp*a**2*b**2*r**6 &
      &              +768.0_cp*a**2*dn**2*r**8+3840.0_cp*a**2*dn*r**8-26112.0_cp*a**2*r**8+ &
      &              6144.0_cp*b**10*dn**2+30720.0_cp*b**10*dn-147456.0_cp*b**10-20480.0_cp &
      &              *b**8*dn**2*r**2-102400.0_cp*b**8*dn*r**2+491520.0_cp*b**8*r**2+24576.0_cp &
      &              *b**6*dn**2*r**4+122880.0_cp*b**6*dn*r**4-589824.0_cp*b**6*r**4-12288.0_cp &
      &              *b**4*dn**2*r**6-61440.0_cp*b**4*dn*r**6+294912.0_cp*b**4*r**6+2048.0_cp &
      &              *b**2*dn**2*r**8+10240.0_cp*b**2*dn*r**8-49152.0_cp*b**2*r**8)/(8192.0_cp &
      &              *dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+1)=3.0_cp*a**4*(33.0_cp*a**12*dn**2-2277.0_cp*a**12+3168.0_cp &
      &              *a**10*b**2*dn**2-186912.0_cp*a**10*b**2-192.0_cp*a**10*dn**2*r**2+11328.0_cp &
      &              *a**10*r**2+36960.0_cp*a**8*b**4*dn**2-1811040.0_cp*a**8*b**4-13440.0_cp &
      &              *a**8*b**2*dn**2*r**2+658560.0_cp*a**8*b**2*r**2+448.0_cp*a**8*dn** &
      &              2*r**4-21952.0_cp*a**8*r**4+118272.0_cp*a**6*b**6*dn**2-4612608.0_cp &
      &              *a**6*b**6-107520.0_cp*a**6*b**4*dn**2*r**2+4193280.0_cp*a**6*b**4* &
      &              r**2+21504.0_cp*a**6*b**2*dn**2*r**4-838656.0_cp*a**6*b**2*r**4-512.0_cp &
      &              *a**6*dn**2*r**6+19968.0_cp*a**6*r**6+126720.0_cp*a**4*b**8*dn**2-3674880.0_cp &
      &              *a**4*b**8-215040.0_cp*a**4*b**6*dn**2*r**2+6236160.0_cp*a**4*b**6* &
      &              r**2+107520.0_cp*a**4*b**4*dn**2*r**4-3118080.0_cp*a**4*b**4*r**4-15360.0_cp &
      &              *a**4*b**2*dn**2*r**6+445440.0_cp*a**4*b**2*r**6+256.0_cp*a**4*dn** &
      &              2*r**8-7424.0_cp*a**4*r**8+45056.0_cp*a**2*b**10*dn**2-856064.0_cp* &
      &              a**2*b**10-122880.0_cp*a**2*b**8*dn**2*r**2+2334720.0_cp*a**2*b**8* &
      &              r**2+114688.0_cp*a**2*b**6*dn**2*r**4-2179072.0_cp*a**2*b**6*r**4-40960.0_cp &
      &              *a**2*b**4*dn**2*r**6+778240.0_cp*a**2*b**4*r**6+4096.0_cp*a**2*b** &
      &              2*dn**2*r**8-77824.0_cp*a**2*b**2*r**8+4096.0_cp*b**12*dn**2-36864.0_cp &
      &              *b**12-16384.0_cp*b**10*dn**2*r**2+147456.0_cp*b**10*r**2+24576.0_cp &
      &              *b**8*dn**2*r**4-221184.0_cp*b**8*r**4-16384.0_cp*b**6*dn**2*r**6+147456.0_cp &
      &              *b**6*r**6+4096.0_cp*b**4*dn**2*r**8-36864.0_cp*b**4*r**8)/(32768.0_cp &
      &              *(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku+2)=a**5*b*(297.0_cp*a**10*dn**2-1485.0_cp*a**10*dn-21978.0_cp &
      &              *a**10+7920.0_cp*a**8*b**2*dn**2-39600.0_cp*a**8*b**2*dn-506880.0_cp &
      &              *a**8*b**2-1440.0_cp*a**8*dn**2*r**2+7200.0_cp*a**8*dn*r**2+92160.0_cp &
      &              *a**8*r**2+44352.0_cp*a**6*b**4*dn**2-221760.0_cp*a**6*b**4*dn-2395008.0_cp &
      &              *a**6*b**4-26880.0_cp*a**6*b**2*dn**2*r**2+134400.0_cp*a**6*b**2*dn &
      &              *r**2+1451520.0_cp*a**6*b**2*r**2+2688.0_cp*a**6*dn**2*r**4-13440.0_cp &
      &              *a**6*dn*r**4-145152.0_cp*a**6*r**4+76032.0_cp*a**4*b**6*dn**2-380160.0_cp &
      &              *a**4*b**6*dn-3345408.0_cp*a**4*b**6-96768.0_cp*a**4*b**4*dn**2*r** &
      &              2+483840.0_cp*a**4*b**4*dn*r**2+4257792.0_cp*a**4*b**4*r**2+32256.0_cp &
      &              *a**4*b**2*dn**2*r**4-161280.0_cp*a**4*b**2*dn*r**4-1419264.0_cp*a**4* &
      &              b**2*r**4-2304.0_cp*a**4*dn**2*r**6+11520.0_cp*a**4*dn*r**6+101376.0_cp &
      &              *a**4*r**6+42240.0_cp*a**2*b**8*dn**2-211200.0_cp*a**2*b**8*dn-1436160.0_cp &
      &              *a**2*b**8-92160.0_cp*a**2*b**6*dn**2*r**2+460800.0_cp*a**2*b**6*dn &
      &              *r**2+3133440.0_cp*a**2*b**6*r**2+64512.0_cp*a**2*b**4*dn**2*r**4-322560.0_cp &
      &              *a**2*b**4*dn*r**4-2193408.0_cp*a**2*b**4*r**4-15360.0_cp*a**2*b**2 &
      &              *dn**2*r**6+76800.0_cp*a**2*b**2*dn*r**6+522240.0_cp*a**2*b**2*r**6 &
      &              +768.0_cp*a**2*dn**2*r**8-3840.0_cp*a**2*dn*r**8-26112.0_cp*a**2*r**8+ &
      &              6144.0_cp*b**10*dn**2-30720.0_cp*b**10*dn-147456.0_cp*b**10-20480.0_cp &
      &              *b**8*dn**2*r**2+102400.0_cp*b**8*dn*r**2+491520.0_cp*b**8*r**2+24576.0_cp &
      &              *b**6*dn**2*r**4-122880.0_cp*b**6*dn*r**4-589824.0_cp*b**6*r**4-12288.0_cp &
      &              *b**4*dn**2*r**6+61440.0_cp*b**4*dn*r**6+294912.0_cp*b**4*r**6+2048.0_cp &
      &              *b**2*dn**2*r**8-10240.0_cp*b**2*dn*r**8-49152.0_cp*b**2*r**8)/(8192.0_cp &
      &              *dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+3)=-a**4*(-22.0_cp*a**12*dn**2+385.0_cp*a**12*dn+3003.0_cp &
      &              *a**12-1782.0_cp*a**10*b**2*dn**2+35640.0_cp*a**10*b**2*dn+241758.0_cp &
      &              *a**10*b**2+108.0_cp*a**10*dn**2*r**2-2160.0_cp*a**10*dn*r**2-14652.0_cp &
      &              *a**10*r**2-15840.0_cp*a**8*b**4*dn**2+396000.0_cp*a**8*b**4*dn+2280960.0_cp &
      &              *a**8*b**4+5760.0_cp*a**8*b**2*dn**2*r**2-144000.0_cp*a**8*b**2*dn* &
      &              r**2-829440.0_cp*a**8*b**2*r**2-192.0_cp*a**8*dn**2*r**4+4800.0_cp* &
      &              a**8*dn*r**4+27648.0_cp*a**8*r**4-29568.0_cp*a**6*b**6*dn**2+1182720.0_cp &
      &              *a**6*b**6*dn+5588352.0_cp*a**6*b**6+26880.0_cp*a**6*b**4*dn**2*r** &
      &              2-1075200.0_cp*a**6*b**4*dn*r**2-5080320.0_cp*a**6*b**4*r**2-5376.0_cp &
      &              *a**6*b**2*dn**2*r**4+215040.0_cp*a**6*b**2*dn*r**4+1016064.0_cp*a**6* &
      &              b**2*r**4+128.0_cp*a**6*dn**2*r**6-5120.0_cp*a**6*dn*r**6-24192.0_cp &
      &              *a**6*r**6+1140480.0_cp*a**4*b**8*dn+4181760.0_cp*a**4*b**8-1935360.0_cp &
      &              *a**4*b**6*dn*r**2-7096320.0_cp*a**4*b**6*r**2+967680.0_cp*a**4*b** &
      &              4*dn*r**4+3548160.0_cp*a**4*b**4*r**4-138240.0_cp*a**4*b**2*dn*r**6 &
      &              -506880.0_cp*a**4*b**2*r**6+2304.0_cp*a**4*dn*r**8+8448.0_cp*a**4*r**8 &
      &              +16896.0_cp*a**2*b**10*dn**2+337920.0_cp*a**2*b**10*dn+861696.0_cp* &
      &              a**2*b**10-46080.0_cp*a**2*b**8*dn**2*r**2-921600.0_cp*a**2*b**8*dn &
      &              *r**2-2350080.0_cp*a**2*b**8*r**2+43008.0_cp*a**2*b**6*dn**2*r**4+860160.0_cp &
      &              *a**2*b**6*dn*r**4+2193408.0_cp*a**2*b**6*r**4-15360.0_cp*a**2*b**4 &
      &              *dn**2*r**6-307200.0_cp*a**2*b**4*dn*r**6-783360.0_cp*a**2*b**4*r** &
      &              6+1536.0_cp*a**2*b**2*dn**2*r**8+30720.0_cp*a**2*b**2*dn*r**8+78336.0_cp &
      &              *a**2*b**2*r**8+4096.0_cp*b**12*dn**2+20480.0_cp*b**12*dn+24576.0_cp &
      &              *b**12-16384.0_cp*b**10*dn**2*r**2-81920.0_cp*b**10*dn*r**2-98304.0_cp &
      &              *b**10*r**2+24576.0_cp*b**8*dn**2*r**4+122880.0_cp*b**8*dn*r**4+147456.0_cp &
      &              *b**8*r**4-16384.0_cp*b**6*dn**2*r**6-81920.0_cp*b**6*dn*r**6-98304.0_cp &
      &              *b**6*r**6+4096.0_cp*b**4*dn**2*r**8+20480.0_cp*b**4*dn*r**8+24576.0_cp &
      &              *b**4*r**8)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp &
      &              )*(dn+3.0_cp))
      stencil(ku+4)=-3.0_cp*a**5*b*(11.0_cp*a**10*dn**2+1155.0_cp*a**10*dn &
      &              +5566.0_cp*a**10+660.0_cp*a**8*b**2*dn**2+29700.0_cp*a**8*b**2*dn+122760.0_cp &
      &              *a**8*b**2-120.0_cp*a**8*dn**2*r**2-5400.0_cp*a**8*dn*r**2-22320.0_cp &
      &              *a**8*r**2+6336.0_cp*a**6*b**4*dn**2+158400.0_cp*a**6*b**4*dn+544896.0_cp &
      &              *a**6*b**4-3840.0_cp*a**6*b**2*dn**2*r**2-96000.0_cp*a**6*b**2*dn*r**2 &
      &              -330240.0_cp*a**6*b**2*r**2+384.0_cp*a**6*dn**2*r**4+9600.0_cp*a**6 &
      &              *dn*r**4+33024.0_cp*a**6*r**4+16896.0_cp*a**4*b**6*dn**2+253440.0_cp &
      &              *a**4*b**6*dn+692736.0_cp*a**4*b**6-21504.0_cp*a**4*b**4*dn**2*r**2 &
      &              -322560.0_cp*a**4*b**4*dn*r**2-881664.0_cp*a**4*b**4*r**2+7168.0_cp &
      &              *a**4*b**2*dn**2*r**4+107520.0_cp*a**4*b**2*dn*r**4+293888.0_cp*a** &
      &              4*b**2*r**4-512.0_cp*a**4*dn**2*r**6-7680.0_cp*a**4*dn*r**6-20992.0_cp &
      &              *a**4*r**6+14080.0_cp*a**2*b**8*dn**2+126720.0_cp*a**2*b**8*dn+253440.0_cp &
      &              *a**2*b**8-30720.0_cp*a**2*b**6*dn**2*r**2-276480.0_cp*a**2*b**6*dn &
      &              *r**2-552960.0_cp*a**2*b**6*r**2+21504.0_cp*a**2*b**4*dn**2*r**4+193536.0_cp &
      &              *a**2*b**4*dn*r**4+387072.0_cp*a**2*b**4*r**4-5120.0_cp*a**2*b**2*dn**2 &
      &              *r**6-46080.0_cp*a**2*b**2*dn*r**6-92160.0_cp*a**2*b**2*r**6+256.0_cp &
      &              *a**2*dn**2*r**8+2304.0_cp*a**2*dn*r**8+4608.0_cp*a**2*r**8+3072.0_cp &
      &              *b**10*dn**2+15360.0_cp*b**10*dn+18432.0_cp*b**10-10240.0_cp*b**8*dn**2 &
      &              *r**2-51200.0_cp*b**8*dn*r**2-61440.0_cp*b**8*r**2+12288.0_cp*b**6*dn**2 &
      &              *r**4+61440.0_cp*b**6*dn*r**4+73728.0_cp*b**6*r**4-6144.0_cp*b**4*dn**2 &
      &              *r**6-30720.0_cp*b**4*dn*r**6-36864.0_cp*b**4*r**6+1024.0_cp*b**2*dn**2 &
      &              *r**8+5120.0_cp*b**2*dn*r**8+6144.0_cp*b**2*r**8)/(8192.0_cp*dn*(dn &
      &              -2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+5)=-a**4*(22.0_cp*a**12*dn**3+429.0_cp*a**12*dn**2-88.0_cp &
      &              *a**12*dn-8151.0_cp*a**12+2508.0_cp*a**10*b**2*dn**3+34848.0_cp*a**10 &
      &              *b**2*dn**2-31812.0_cp*a**10*b**2*dn-618552.0_cp*a**10*b**2-152.0_cp &
      &              *a**10*dn**3*r**2-2112.0_cp*a**10*dn**2*r**2+1928.0_cp*a**10*dn*r** &
      &              2+37488.0_cp*a**10*r**2+33660.0_cp*a**8*b**4*dn**3+320760.0_cp*a**8 &
      &              *b**4*dn**2-580140.0_cp*a**8*b**4*dn-5381640.0_cp*a**8*b**4-12240.0_cp &
      &              *a**8*b**2*dn**3*r**2-116640.0_cp*a**8*b**2*dn**2*r**2+210960.0_cp* &
      &              a**8*b**2*dn*r**2+1956960.0_cp*a**8*b**2*r**2+408.0_cp*a**8*dn**3*r**4 &
      &              +3888.0_cp*a**8*dn**2*r**4-7032.0_cp*a**8*dn*r**4-65232.0_cp*a**8*r**4 &
      &              +118272.0_cp*a**6*b**6*dn**3+709632.0_cp*a**6*b**6*dn**2-2247168.0_cp &
      &              *a**6*b**6*dn-11708928.0_cp*a**6*b**6-107520.0_cp*a**6*b**4*dn**3*r**2 &
      &              -645120.0_cp*a**6*b**4*dn**2*r**2+2042880.0_cp*a**6*b**4*dn*r**2+10644480.0_cp &
      &              *a**6*b**4*r**2+21504.0_cp*a**6*b**2*dn**3*r**4+129024.0_cp*a**6*b**2* &
      &              dn**2*r**4-408576.0_cp*a**6*b**2*dn*r**4-2128896.0_cp*a**6*b**2*r** &
      &              4-512.0_cp*a**6*dn**3*r**6-3072.0_cp*a**6*dn**2*r**6+9728.0_cp*a**6 &
      &              *dn*r**6+50688.0_cp*a**6*r**6+126720.0_cp*a**4*b**8*dn**3+380160.0_cp &
      &              *a**4*b**8*dn**2-2407680.0_cp*a**4*b**8*dn-7223040.0_cp*a**4*b**8-215040.0_cp &
      &              *a**4*b**6*dn**3*r**2-645120.0_cp*a**4*b**6*dn**2*r**2+4085760.0_cp &
      &              *a**4*b**6*dn*r**2+12257280.0_cp*a**4*b**6*r**2+107520.0_cp*a**4*b**4* &
      &              dn**3*r**4+322560.0_cp*a**4*b**4*dn**2*r**4-2042880.0_cp*a**4*b**4*dn &
      &              *r**4-6128640.0_cp*a**4*b**4*r**4-15360.0_cp*a**4*b**2*dn**3*r**6-46080.0_cp &
      &              *a**4*b**2*dn**2*r**6+291840.0_cp*a**4*b**2*dn*r**6+875520.0_cp*a** &
      &              4*b**2*r**6+256.0_cp*a**4*dn**3*r**8+768.0_cp*a**4*dn**2*r**8-4864.0_cp &
      &              *a**4*dn*r**8-14592.0_cp*a**4*r**8+33792.0_cp*a**2*b**10*dn**3-642048.0_cp &
      &              *a**2*b**10*dn-1013760.0_cp*a**2*b**10-92160.0_cp*a**2*b**8*dn**3*r**2 &
      &              +1751040.0_cp*a**2*b**8*dn*r**2+2764800.0_cp*a**2*b**8*r**2+86016.0_cp &
      &              *a**2*b**6*dn**3*r**4-1634304.0_cp*a**2*b**6*dn*r**4-2580480.0_cp*a**2 &
      &              *b**6*r**4-30720.0_cp*a**2*b**4*dn**3*r**6+583680.0_cp*a**2*b**4*dn &
      &              *r**6+921600.0_cp*a**2*b**4*r**6+3072.0_cp*a**2*b**2*dn**3*r**8-58368.0_cp &
      &              *a**2*b**2*dn*r**8-92160.0_cp*a**2*b**2*r**8-1024.0_cp*b**12*dn**3-6144.0_cp &
      &              *b**12*dn**2-11264.0_cp*b**12*dn-6144.0_cp*b**12+4096.0_cp*b**10*dn**3 &
      &              *r**2+24576.0_cp*b**10*dn**2*r**2+45056.0_cp*b**10*dn*r**2+24576.0_cp &
      &              *b**10*r**2-6144.0_cp*b**8*dn**3*r**4-36864.0_cp*b**8*dn**2*r**4-67584.0_cp &
      &              *b**8*dn*r**4-36864.0_cp*b**8*r**4+4096.0_cp*b**6*dn**3*r**6+24576.0_cp &
      &              *b**6*dn**2*r**6+45056.0_cp*b**6*dn*r**6+24576.0_cp*b**6*r**6-1024.0_cp &
      &              *b**4*dn**3*r**8-6144.0_cp*b**4*dn**2*r**8-11264.0_cp*b**4*dn*r**8-6144.0_cp &
      &              *b**4*r**8)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+6)=-a**5*b*(231.0_cp*a**10*dn**3+1782.0_cp*a**10*dn**2-7359.0_cp &
      &              *a**10*dn-47718.0_cp*a**10+6380.0_cp*a**8*b**2*dn**3+33000.0_cp*a** &
      &              8*b**2*dn**2-207020.0_cp*a**8*b**2*dn-957000.0_cp*a**8*b**2-1160.0_cp &
      &              *a**8*dn**3*r**2-6000.0_cp*a**8*dn**2*r**2+37640.0_cp*a**8*dn*r**2+174000.0_cp &
      &              *a**8*r**2+34848.0_cp*a**6*b**4*dn**3+95040.0_cp*a**6*b**4*dn**2-1137312.0_cp &
      &              *a**6*b**4*dn-3706560.0_cp*a**6*b**4-21120.0_cp*a**6*b**2*dn**3*r** &
      &              2-57600.0_cp*a**6*b**2*dn**2*r**2+689280.0_cp*a**6*b**2*dn*r**2+2246400.0_cp &
      &              *a**6*b**2*r**2+2112.0_cp*a**6*dn**3*r**4+5760.0_cp*a**6*dn**2*r**4 &
      &              -68928.0_cp*a**6*dn*r**4-224640.0_cp*a**6*r**4+50688.0_cp*a**4*b**6 &
      &              *dn**3-1723392.0_cp*a**4*b**6*dn-3801600.0_cp*a**4*b**6-64512.0_cp* &
      &              a**4*b**4*dn**3*r**2+2193408.0_cp*a**4*b**4*dn*r**2+4838400.0_cp*a**4* &
      &              b**4*r**2+21504.0_cp*a**4*b**2*dn**3*r**4-731136.0_cp*a**4*b**2*dn* &
      &              r**4-1612800.0_cp*a**4*b**2*r**4-1536.0_cp*a**4*dn**3*r**6+52224.0_cp &
      &              *a**4*dn*r**6+115200.0_cp*a**4*r**6+14080.0_cp*a**2*b**8*dn**3-84480.0_cp &
      &              *a**2*b**8*dn**2-689920.0_cp*a**2*b**8*dn-929280.0_cp*a**2*b**8-30720.0_cp &
      &              *a**2*b**6*dn**3*r**2+184320.0_cp*a**2*b**6*dn**2*r**2+1505280.0_cp &
      &              *a**2*b**6*dn*r**2+2027520.0_cp*a**2*b**6*r**2+21504.0_cp*a**2*b**4 &
      &              *dn**3*r**4-129024.0_cp*a**2*b**4*dn**2*r**4-1053696.0_cp*a**2*b**4 &
      &              *dn*r**4-1419264.0_cp*a**2*b**4*r**4-5120.0_cp*a**2*b**2*dn**3*r**6 &
      &              +30720.0_cp*a**2*b**2*dn**2*r**6+250880.0_cp*a**2*b**2*dn*r**6+337920.0_cp &
      &              *a**2*b**2*r**6+256.0_cp*a**2*dn**3*r**8-1536.0_cp*a**2*dn**2*r**8-12544.0_cp &
      &              *a**2*dn*r**8-16896.0_cp*a**2*r**8-3072.0_cp*b**10*dn**3-18432.0_cp &
      &              *b**10*dn**2-33792.0_cp*b**10*dn-18432.0_cp*b**10+10240.0_cp*b**8*dn**3 &
      &              *r**2+61440.0_cp*b**8*dn**2*r**2+112640.0_cp*b**8*dn*r**2+61440.0_cp &
      &              *b**8*r**2-12288.0_cp*b**6*dn**3*r**4-73728.0_cp*b**6*dn**2*r**4-135168.0_cp &
      &              *b**6*dn*r**4-73728.0_cp*b**6*r**4+6144.0_cp*b**4*dn**3*r**6+36864.0_cp &
      &              *b**4*dn**2*r**6+67584.0_cp*b**4*dn*r**6+36864.0_cp*b**4*r**6-1024.0_cp &
      &              *b**2*dn**3*r**8-6144.0_cp*b**2*dn**2*r**8-11264.0_cp*b**2*dn*r**8-6144.0_cp &
      &              *b**2*r**8)/(8192.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+7)=3.0_cp*a**6*(-10.0_cp*a**10*dn**3-39.0_cp*a**10*dn**2 &
      &              +495.0_cp*a**10*dn+2106.0_cp*a**10-858.0_cp*a**8*b**2*dn**3-1452.0_cp &
      &              *a**8*b**2*dn**2+43362.0_cp*a**8*b**2*dn+143748.0_cp*a**8*b**2+52.0_cp &
      &              *a**8*dn**3*r**2+88.0_cp*a**8*dn**2*r**2-2628.0_cp*a**8*dn*r**2-8712.0_cp &
      &              *a**8*r**2-7920.0_cp*a**6*b**4*dn**3+7920.0_cp*a**6*b**4*dn**2+427680.0_cp &
      &              *a**6*b**4*dn+1077120.0_cp*a**6*b**4+2880.0_cp*a**6*b**2*dn**3*r**2 &
      &              -2880.0_cp*a**6*b**2*dn**2*r**2-155520.0_cp*a**6*b**2*dn*r**2-391680.0_cp &
      &              *a**6*b**2*r**2-96.0_cp*a**6*dn**3*r**4+96.0_cp*a**6*dn**2*r**4+5184.0_cp &
      &              *a**6*dn*r**4+13056.0_cp*a**6*r**4-14784.0_cp*a**4*b**6*dn**3+88704.0_cp &
      &              *a**4*b**6*dn**2+1020096.0_cp*a**4*b**6*dn+1862784.0_cp*a**4*b**6+13440.0_cp &
      &              *a**4*b**4*dn**3*r**2-80640.0_cp*a**4*b**4*dn**2*r**2-927360.0_cp*a**4 &
      &              *b**4*dn*r**2-1693440.0_cp*a**4*b**4*r**2-2688.0_cp*a**4*b**2*dn**3 &
      &              *r**4+16128.0_cp*a**4*b**2*dn**2*r**4+185472.0_cp*a**4*b**2*dn*r**4 &
      &              +338688.0_cp*a**4*b**2*r**4+64.0_cp*a**4*dn**3*r**6-384.0_cp*a**4*dn**2 &
      &              *r**6-4416.0_cp*a**4*dn*r**6-8064.0_cp*a**4*r**6+126720.0_cp*a**2*b**8 &
      &              *dn**2+633600.0_cp*a**2*b**8*dn+760320.0_cp*a**2*b**8-215040.0_cp*a**2 &
      &              *b**6*dn**2*r**2-1075200.0_cp*a**2*b**6*dn*r**2-1290240.0_cp*a**2*b**6 &
      &              *r**2+107520.0_cp*a**2*b**4*dn**2*r**4+537600.0_cp*a**2*b**4*dn*r** &
      &              4+645120.0_cp*a**2*b**4*r**4-15360.0_cp*a**2*b**2*dn**2*r**6-76800.0_cp &
      &              *a**2*b**2*dn*r**6-92160.0_cp*a**2*b**2*r**6+256.0_cp*a**2*dn**2*r**8+ &
      &              1280.0_cp*a**2*dn*r**8+1536.0_cp*a**2*r**8+5632.0_cp*b**10*dn**3+33792.0_cp &
      &              *b**10*dn**2+61952.0_cp*b**10*dn+33792.0_cp*b**10-15360.0_cp*b**8*dn**3 &
      &              *r**2-92160.0_cp*b**8*dn**2*r**2-168960.0_cp*b**8*dn*r**2-92160.0_cp &
      &              *b**8*r**2+14336.0_cp*b**6*dn**3*r**4+86016.0_cp*b**6*dn**2*r**4+157696.0_cp &
      &              *b**6*dn*r**4+86016.0_cp*b**6*r**4-5120.0_cp*b**4*dn**3*r**6-30720.0_cp &
      &              *b**4*dn**2*r**6-56320.0_cp*b**4*dn*r**6-30720.0_cp*b**4*r**6+512.0_cp &
      &              *b**2*dn**3*r**8+3072.0_cp*b**2*dn**2*r**8+5632.0_cp*b**2*dn*r**8+3072.0_cp &
      &              *b**2*r**8)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+8)=a**7*b*(-129.0_cp*a**8*dn**3+378.0_cp*a**8*dn**2+10461.0_cp &
      &              *a**8*dn+28098.0_cp*a**8-2200.0_cp*a**6*b**2*dn**3+18480.0_cp*a**6* &
      &              b**2*dn**2+226600.0_cp*a**6*b**2*dn+480480.0_cp*a**6*b**2+400.0_cp* &
      &              a**6*dn**3*r**2-3360.0_cp*a**6*dn**2*r**2-41200.0_cp*a**6*dn*r**2-87360.0_cp &
      &              *a**6*r**2-3168.0_cp*a**4*b**4*dn**3+133056.0_cp*a**4*b**4*dn**2+915552.0_cp &
      &              *a**4*b**4*dn+1463616.0_cp*a**4*b**4+1920.0_cp*a**4*b**2*dn**3*r**2 &
      &              -80640.0_cp*a**4*b**2*dn**2*r**2-554880.0_cp*a**4*b**2*dn*r**2-887040.0_cp &
      &              *a**4*b**2*r**2-192.0_cp*a**4*dn**3*r**4+8064.0_cp*a**4*dn**2*r**4+55488.0_cp &
      &              *a**4*dn*r**4+88704.0_cp*a**4*r**4+12672.0_cp*a**2*b**6*dn**3+228096.0_cp &
      &              *a**2*b**6*dn**2+899712.0_cp*a**2*b**6*dn+988416.0_cp*a**2*b**6-16128.0_cp &
      &              *a**2*b**4*dn**3*r**2-290304.0_cp*a**2*b**4*dn**2*r**2-1145088.0_cp &
      &              *a**2*b**4*dn*r**2-1257984.0_cp*a**2*b**4*r**2+5376.0_cp*a**2*b**2*dn**3 &
      &              *r**4+96768.0_cp*a**2*b**2*dn**2*r**4+381696.0_cp*a**2*b**2*dn*r**4 &
      &              +419328.0_cp*a**2*b**2*r**4-384.0_cp*a**2*dn**3*r**6-6912.0_cp*a**2 &
      &              *dn**2*r**6-27264.0_cp*a**2*dn*r**6-29952.0_cp*a**2*r**6+14080.0_cp &
      &              *b**8*dn**3+84480.0_cp*b**8*dn**2+154880.0_cp*b**8*dn+84480.0_cp*b**8- &
      &              30720.0_cp*b**6*dn**3*r**2-184320.0_cp*b**6*dn**2*r**2-337920.0_cp* &
      &              b**6*dn*r**2-184320.0_cp*b**6*r**2+21504.0_cp*b**4*dn**3*r**4+129024.0_cp &
      &              *b**4*dn**2*r**4+236544.0_cp*b**4*dn*r**4+129024.0_cp*b**4*r**4-5120.0_cp &
      &              *b**2*dn**3*r**6-30720.0_cp*b**2*dn**2*r**6-56320.0_cp*b**2*dn*r**6 &
      &              -30720.0_cp*b**2*r**6+256.0_cp*dn**3*r**8+1536.0_cp*dn**2*r**8+2816.0_cp &
      &              *dn*r**8+1536.0_cp*r**8)/(8192.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+9)=a**8*(-9.0_cp*a**8*dn**3+117.0_cp*a**8*dn**2+1401.0_cp &
      &              *a**8*dn+3237.0_cp*a**8-264.0_cp*a**6*b**2*dn**3+12672.0_cp*a**6*b**2* &
      &              dn**2+100056.0_cp*a**6*b**2*dn+186912.0_cp*a**6*b**2+16.0_cp*a**6*dn**3 &
      &              *r**2-768.0_cp*a**6*dn**2*r**2-6064.0_cp*a**6*dn*r**2-11328.0_cp*a**6* &
      &              r**2+3960.0_cp*a**4*b**4*dn**3+142560.0_cp*a**4*b**4*dn**2+756360.0_cp &
      &              *a**4*b**4*dn+1092960.0_cp*a**4*b**4-1440.0_cp*a**4*b**2*dn**3*r**2 &
      &              -51840.0_cp*a**4*b**2*dn**2*r**2-275040.0_cp*a**4*b**2*dn*r**2-397440.0_cp &
      &              *a**4*b**2*r**2+48.0_cp*a**4*dn**3*r**4+1728.0_cp*a**4*dn**2*r**4+9168.0_cp &
      &              *a**4*dn*r**4+13248.0_cp*a**4*r**4+29568.0_cp*a**2*b**6*dn**3+354816.0_cp &
      &              *a**2*b**6*dn**2+1212288.0_cp*a**2*b**6*dn+1241856.0_cp*a**2*b**6-26880.0_cp &
      &              *a**2*b**4*dn**3*r**2-322560.0_cp*a**2*b**4*dn**2*r**2-1102080.0_cp &
      &              *a**2*b**4*dn*r**2-1128960.0_cp*a**2*b**4*r**2+5376.0_cp*a**2*b**2*dn**3 &
      &              *r**4+64512.0_cp*a**2*b**2*dn**2*r**4+220416.0_cp*a**2*b**2*dn*r**4 &
      &              +225792.0_cp*a**2*b**2*r**4-128.0_cp*a**2*dn**3*r**6-1536.0_cp*a**2 &
      &              *dn**2*r**6-5248.0_cp*a**2*dn*r**6-5376.0_cp*a**2*r**6+31680.0_cp*b**8 &
      &              *dn**3+190080.0_cp*b**8*dn**2+348480.0_cp*b**8*dn+190080.0_cp*b**8-53760.0_cp &
      &              *b**6*dn**3*r**2-322560.0_cp*b**6*dn**2*r**2-591360.0_cp*b**6*dn*r**2- &
      &              322560.0_cp*b**6*r**2+26880.0_cp*b**4*dn**3*r**4+161280.0_cp*b**4*dn**2 &
      &              *r**4+295680.0_cp*b**4*dn*r**4+161280.0_cp*b**4*r**4-3840.0_cp*b**2 &
      &              *dn**3*r**6-23040.0_cp*b**2*dn**2*r**6-42240.0_cp*b**2*dn*r**6-23040.0_cp &
      &              *b**2*r**6+64.0_cp*dn**3*r**8+384.0_cp*dn**2*r**8+704.0_cp*dn*r**8+384.0_cp &
      &              *r**8)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+10)=-3.0_cp*a**9*b*(-7.0_cp*a**6*dn**3-342.0_cp*a**6*dn**2- &
      &              2117.0_cp*a**6*dn-3582.0_cp*a**6-440.0_cp*a**4*b**2*dn**3-7920.0_cp &
      &              *a**4*b**2*dn**2-35640.0_cp*a**4*b**2*dn-47520.0_cp*a**4*b**2+80.0_cp &
      &              *a**4*dn**3*r**2+1440.0_cp*a**4*dn**2*r**2+6480.0_cp*a**4*dn*r**2+8640.0_cp &
      &              *a**4*r**2-3168.0_cp*a**2*b**4*dn**3-31680.0_cp*a**2*b**4*dn**2-98208.0_cp &
      &              *a**2*b**4*dn-95040.0_cp*a**2*b**4+1920.0_cp*a**2*b**2*dn**3*r**2+19200.0_cp &
      &              *a**2*b**2*dn**2*r**2+59520.0_cp*a**2*b**2*dn*r**2+57600.0_cp*a**2* &
      &              b**2*r**2-192.0_cp*a**2*dn**3*r**4-1920.0_cp*a**2*dn**2*r**4-5952.0_cp &
      &              *a**2*dn*r**4-5760.0_cp*a**2*r**4-4224.0_cp*b**6*dn**3-25344.0_cp*b**6 &
      &              *dn**2-46464.0_cp*b**6*dn-25344.0_cp*b**6+5376.0_cp*b**4*dn**3*r**2 &
      &              +32256.0_cp*b**4*dn**2*r**2+59136.0_cp*b**4*dn*r**2+32256.0_cp*b**4 &
      &              *r**2-1792.0_cp*b**2*dn**3*r**4-10752.0_cp*b**2*dn**2*r**4-19712.0_cp &
      &              *b**2*dn*r**4-10752.0_cp*b**2*r**4+128.0_cp*dn**3*r**6+768.0_cp*dn**2* &
      &              r**6+1408.0_cp*dn*r**6+768.0_cp*r**6)/(8192.0_cp*dn*(dn-3.0_cp)*(dn &
      &              -2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+11)=-a**10*(-6.0_cp*a**6*dn**3-135.0_cp*a**6*dn**2-711.0_cp*a**6* &
      &              dn-1110.0_cp*a**6-726.0_cp*a**4*b**2*dn**3-9900.0_cp*a**4*b** &
      &              2*dn**2-39666.0_cp*a**4*b**2*dn-49500.0_cp*a**4*b**2+44.0_cp  &
      &              *a**4*dn**3*r**2+600.0_cp*a**4*dn**2*r**2+2404.0_cp*a**4*dn*  &
      &              r**2+3000.0_cp*a**4*r**2-7920.0_cp*a**2*b**4*dn**3-71280.0_cp*&
      &              a**2*b**4*dn**2-205920.0_cp*a**2*b**4*dn-190080.0_cp*a**2*b**4&
      &              +2880.0_cp*a**2*b**2*dn**3*r**2+25920.0_cp*a**2*b**2*dn**2*   &
      &              r**2+74880.0_cp*a**2*b**2*dn*r**2+69120.0_cp*a**2*b**2*r**2-  &
      &              96.0_cp*a**2*dn**3*r**4-864.0_cp*a**2*dn**2*r**4-2496.0_cp    &
      &              *a**2*dn*r**4-2304.0_cp*a**2*r**4-14784.0_cp*b**6*dn**3-      &
      &              88704.0_cp*b**6*dn**2-162624.0_cp*b**6*dn-88704.0_cp*b**6+    &
      &              13440.0_cp*b**4*dn**3*r**2+80640.0_cp*b**4*dn**2*r**2+        &
      &              147840.0_cp*b**4*dn*r**2+80640.0_cp*b**4*r**2-2688.0_cp*b**2* &
      &              dn**3*r**4-16128.0_cp*b**2*dn**2*r**4-29568.0_cp*b**2*dn*r**4-&
      &              16128.0_cp*b**2*r**4+64.0_cp*dn**3*r**6+384.0_cp*dn**2*r**6+  &
      &              704.0_cp*dn*r**6+384.0_cp*r**6)/(16384.0_cp*dn*(dn-3.0_cp)*   &
      &              (dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+12)=a**11*b*(51.0_cp*a**4*dn**2+441.0_cp*a**4*dn+858.0_cp*a**4+  &
      &              1100.0_cp*a**2*b**2*dn**2+5940.0_cp*a**2*b**2*dn+7480.0_cp*  &
      &              a**2*b**2-200.0_cp*a**2*dn**2*r**2-1080.0_cp*a**2*dn*r**2-   &
      &              1360.0_cp*a**2*r**2+3168.0_cp*b**4*dn**2+9504.0_cp*b**4*dn+  &
      &              6336.0_cp*b**4-1920.0_cp*b**2*dn**2*r**2-5760.0_cp*b**2*dn*  &
      &              r**2-3840.0_cp*b**2*r**2+192.0_cp*dn**2*r**4+576.0_cp*dn*r**4&
      &              +384.0_cp*r**4)/(8192.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-   &
      &              1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+13)=3.0_cp*a**12*(2.0_cp*a**4*dn**2+15.0_cp*a**4*dn+27.0_cp      &
      &              *a**4+132.0_cp*a**2*b**2*dn**2+660.0_cp*a**2*b**2*dn+        &
      &              792.0_cp*a**2*b**2-8.0_cp*a**2*dn**2*r**2-40.0_cp*a**2*      &
      &              dn*r**2-48.0_cp*a**2*r**2+660.0_cp*b**4*dn**2+1980.0_cp*     &
      &              b**4*dn+1320.0_cp*b**4-240.0_cp*b**2*dn**2*r**2-720.0_cp*    &
      &              b**2*dn*r**2-480.0_cp*b**2*r**2+8.0_cp*dn**2*r**4+24.0_cp*dn &
      &              *r**4+16.0_cp*r**4)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*( &
      &              dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+14)=-a**13*b*(-21.0_cp*a**2*dn-57.0_cp*a**2-220.0_cp*b**        &
      &              2*dn-220.0_cp*b**2+40.0_cp*dn*r**2+40.0_cp*r**2)/(8192.0_cp*&
      &              dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+15)=-a**14*(-2.0_cp*a**2*dn-5.0_cp*a**2-66.0_cp*b**2*dn-66.0_cp   &
      &              *b**2+4.0_cp*dn*r**2+4.0_cp*r**2)/(16384.0_cp*dn*(dn-3.0_cp)* &
      &              (dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+16)=3.0_cp*a**15*b/(8192.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp) &
      &              *(dn-1.0_cp))
      stencil(ku+17)=a**16/(65536.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              ))
      stencil(ku+18:len_stencil) = 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4rmult4hmult8
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
