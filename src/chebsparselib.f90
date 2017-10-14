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
   &         intcheb4rmult4hmult8, intcheb4rmult4hmult6,                 &
   &         intcheb4rmult4hmult8laplrot, intcheb4rmult4hmult8laplrot2,  &
   &         intcheb4hmult2

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
   function intcheb4hmult2(a, b, r, n, len_stencil) result(stencil)

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

      stencil(1:ku-6) = 0.0_cp
      stencil(ku-5)=-a**6/(64.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-4)=-a**5*b/(16.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-3)=a**4*(a**2*dn+5.0_cp*a**2-2.0_cp*b**2*dn+2.0_cp*b**2+2.0_cp &
      &             *dn*r**2-2.0_cp*r**2)/(32.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp)*  &
      &             (dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-2)=3.0_cp*a**5*b/(16.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku-1)=-a**4*(-a**2*dn+17.0_cp*a**2-16.0_cp*b**2*dn+32.0_cp*      &
      &             b**2+16.0_cp*dn*r**2-32.0_cp*r**2)/(64.0_cp*dn*(dn-2.0_cp)*&
      &             (dn-1.0_cp)*(dn+1.0_cp)*(dn+3.0_cp))
      stencil(ku)=-a**5*b*(dn+8.0_cp)/(8.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)* &
      &            (dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+1)=a**4*(-a**2*dn**2+19.0_cp*a**2-6.0_cp*b**2*dn**2+54.0_cp   &
      &             *b**2+6.0_cp*dn**2*r**2-54.0_cp*r**2)/(16.0_cp*(dn-3.0_cp)*&
      &             (dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+2)=-a**5*b*(dn-8.0_cp)/(8.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp &
      &             )*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+3)=-a**4*(-a**2*dn-17.0_cp*a**2-16.0_cp*b**2*dn-32.0_cp*b**2+&
      &             16.0_cp*dn*r**2+32.0_cp*r**2)/(64.0_cp*dn*(dn-3.0_cp)*    &
      &             (dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+4)=3.0_cp*a**5*b/(16.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+5)=a**4*(a**2*dn-5.0_cp*a**2-2.0_cp*b**2*dn-2.0_cp*b**2+2.0_cp &
      &             *dn*r**2+2.0_cp*r**2)/(32.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*  &
      &             (dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+6)=-a**5*b/(16.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp))
      stencil(ku+7)=-a**6/(64.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp))
      stencil(ku+8:len_stencil) = 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4hmult2
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
   function intcheb4rmult4hmult6(a, b, r, n, len_stencil) result(stencil)
      !
      ! Fourth integral of r^4*(ro^2-r^2)^3 operator
      !

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

      stencil(1:ku-14) = 0.0_cp
      stencil(ku-13)=-a**14/(16384.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku-12)=-5.0_cp*a**13*b/(4096.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp &
      &              )*(dn+3.0_cp))
      stencil(ku-11)=3.0_cp*a**12*(-a**2*dn+3.0_cp*a**2-30.0_cp*b**2*dn+30.0_cp &
      &              *b**2+2.0_cp*dn*r**2-2.0_cp*r**2)/(8192.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-10)=a**11*b*(-25.0_cp*a**2*dn+85.0_cp*a**2-240.0_cp*b**2 &
      &              *dn+240.0_cp*b**2+48.0_cp*dn*r**2-48.0_cp*r**2)/(4096.0_cp*dn*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-9)=-a**10*(11.0_cp*a**4*dn**2-117.0_cp*a**4*dn+250.0_cp* &
      &              a**4+720.0_cp*a**2*b**2*dn**2-4320.0_cp*a**2*b**2*dn+5760.0_cp*a**2 &
      &              *b**2-48.0_cp*a**2*dn**2*r**2+288.0_cp*a**2*dn*r**2-384.0_cp*a**2*r**2 &
      &              +3360.0_cp*b**4*dn**2-10080.0_cp*b**4*dn+6720.0_cp*b**4-1344.0_cp*b**2 &
      &              *dn**2*r**2+4032.0_cp*b**2*dn*r**2-2688.0_cp*b**2*r**2+48.0_cp*dn** &
      &              2*r**4-144.0_cp*dn*r**4+96.0_cp*r**4)/(16384.0_cp*dn*(dn-2.0_cp)*(dn &
      &              -1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-8)=-3.0_cp*a**9*b*(5.0_cp*a**4*dn**2-75.0_cp*a**4*dn+180.0_cp &
      &              *a**4+120.0_cp*a**2*b**2*dn**2-840.0_cp*a**2*b**2*dn+1200.0_cp*a**2 &
      &              *b**2-24.0_cp*a**2*dn**2*r**2+168.0_cp*a**2*dn*r**2-240.0_cp*a**2*r**2 &
      &              +336.0_cp*b**4*dn**2-1008.0_cp*b**4*dn+672.0_cp*b**4-224.0_cp*b**2*dn**2 &
      &              *r**2+672.0_cp*b**2*dn*r**2-448.0_cp*b**2*r**2+24.0_cp*dn**2*r**4-72.0_cp &
      &              *dn*r**4+48.0_cp*r**4)/(2048.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-7)=a**8*(a**6*dn**3+48.0_cp*a**6*dn**2-379.0_cp*a**6*dn+708.0_cp &
      &              *a**6-90.0_cp*a**4*b**2*dn**3+3240.0_cp*a**4*b**2*dn**2-17190.0_cp* &
      &              a**4*b**2*dn+24840.0_cp*a**4*b**2+6.0_cp*a**4*dn**3*r**2-216.0_cp*a**4 &
      &              *dn**2*r**2+1146.0_cp*a**4*dn*r**2-1656.0_cp*a**4*r**2-1680.0_cp*a**2* &
      &              b**4*dn**3+20160.0_cp*a**2*b**4*dn**2-68880.0_cp*a**2*b**4*dn+70560.0_cp &
      &              *a**2*b**4+672.0_cp*a**2*b**2*dn**3*r**2-8064.0_cp*a**2*b**2*dn**2* &
      &              r**2+27552.0_cp*a**2*b**2*dn*r**2-28224.0_cp*a**2*b**2*r**2-24.0_cp &
      &              *a**2*dn**3*r**4+288.0_cp*a**2*dn**2*r**4-984.0_cp*a**2*dn*r**4+1008.0_cp &
      &              *a**2*r**4-3360.0_cp*b**6*dn**3+20160.0_cp*b**6*dn**2-36960.0_cp*b**6* &
      &              dn+20160.0_cp*b**6+3360.0_cp*b**4*dn**3*r**2-20160.0_cp*b**4*dn**2* &
      &              r**2+36960.0_cp*b**4*dn*r**2-20160.0_cp*b**4*r**2-720.0_cp*b**2*dn**3* &
      &              r**4+4320.0_cp*b**2*dn**2*r**4-7920.0_cp*b**2*dn*r**4+4320.0_cp*b** &
      &              2*r**4+16.0_cp*dn**3*r**6-96.0_cp*dn**2*r**6+176.0_cp*dn*r**6-96.0_cp &
      &              *r**6)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-6)=a**7*b*(25.0_cp*a**6*dn**3+210.0_cp*a**6*dn**2-2575.0_cp &
      &              *a**6*dn+5460.0_cp*a**6+120.0_cp*a**4*b**2*dn**3+5040.0_cp*a**4*b** &
      &              2*dn**2-34680.0_cp*a**4*b**2*dn+55440.0_cp*a**4*b**2-24.0_cp*a**4*dn**3 &
      &              *r**2-1008.0_cp*a**4*dn**2*r**2+6936.0_cp*a**4*dn*r**2-11088.0_cp*a**4 &
      &              *r**2-1008.0_cp*a**2*b**4*dn**3+18144.0_cp*a**2*b**4*dn**2-71568.0_cp &
      &              *a**2*b**4*dn+78624.0_cp*a**2*b**4+672.0_cp*a**2*b**2*dn**3*r**2-12096.0_cp &
      &              *a**2*b**2*dn**2*r**2+47712.0_cp*a**2*b**2*dn*r**2-52416.0_cp*a**2* &
      &              b**2*r**2-72.0_cp*a**2*dn**3*r**4+1296.0_cp*a**2*dn**2*r**4-5112.0_cp &
      &              *a**2*dn*r**4+5616.0_cp*a**2*r**4-1920.0_cp*b**6*dn**3+11520.0_cp*b**6 &
      &              *dn**2-21120.0_cp*b**6*dn+11520.0_cp*b**6+2688.0_cp*b**4*dn**3*r**2 &
      &              -16128.0_cp*b**4*dn**2*r**2+29568.0_cp*b**4*dn*r**2-16128.0_cp*b**4 &
      &              *r**2-960.0_cp*b**2*dn**3*r**4+5760.0_cp*b**2*dn**2*r**4-10560.0_cp &
      &              *b**2*dn*r**4+5760.0_cp*b**2*r**4+64.0_cp*dn**3*r**6-384.0_cp*dn**2 &
      &              *r**6+704.0_cp*dn*r**6-384.0_cp*r**6)/(2048.0_cp*dn*(dn-3.0_cp)*(dn &
      &              -2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-5)=-3.0_cp*a**6*(-13.0_cp*a**8*dn**3+22.0_cp*a**8*dn**2+657.0_cp &
      &              *a**8*dn-2178.0_cp*a**8-720.0_cp*a**6*b**2*dn**3-720.0_cp*a**6*b**2 &
      &              *dn**2+38880.0_cp*a**6*b**2*dn-97920.0_cp*a**6*b**2+48.0_cp*a**6*dn**3 &
      &              *r**2+48.0_cp*a**6*dn**2*r**2-2592.0_cp*a**6*dn*r**2+6528.0_cp*a**6 &
      &              *r**2-3360.0_cp*a**4*b**4*dn**3-20160.0_cp*a**4*b**4*dn**2+231840.0_cp &
      &              *a**4*b**4*dn-423360.0_cp*a**4*b**4+1344.0_cp*a**4*b**2*dn**3*r**2+8064.0_cp &
      &              *a**4*b**2*dn**2*r**2-92736.0_cp*a**4*b**2*dn*r**2+169344.0_cp*a**4 &
      &              *b**2*r**2-48.0_cp*a**4*dn**3*r**4-288.0_cp*a**4*dn**2*r**4+3312.0_cp &
      &              *a**4*dn*r**4-6048.0_cp*a**4*r**4-53760.0_cp*a**2*b**6*dn**2+268800.0_cp &
      &              *a**2*b**6*dn-322560.0_cp*a**2*b**6+53760.0_cp*a**2*b**4*dn**2*r**2 &
      &              -268800.0_cp*a**2*b**4*dn*r**2+322560.0_cp*a**2*b**4*r**2-11520.0_cp &
      &              *a**2*b**2*dn**2*r**4+57600.0_cp*a**2*b**2*dn*r**4-69120.0_cp*a**2* &
      &              b**2*r**4+256.0_cp*a**2*dn**2*r**6-1280.0_cp*a**2*dn*r**6+1536.0_cp &
      &              *a**2*r**6+3840.0_cp*b**8*dn**3-23040.0_cp*b**8*dn**2+42240.0_cp*b**8* &
      &              dn-23040.0_cp*b**8-7168.0_cp*b**6*dn**3*r**2+43008.0_cp*b**6*dn**2* &
      &              r**2-78848.0_cp*b**6*dn*r**2+43008.0_cp*b**6*r**2+3840.0_cp*b**4*dn**3 &
      &              *r**4-23040.0_cp*b**4*dn**2*r**4+42240.0_cp*b**4*dn*r**4-23040.0_cp &
      &              *b**4*r**4-512.0_cp*b**2*dn**3*r**6+3072.0_cp*b**2*dn**2*r**6-5632.0_cp &
      &              *b**2*dn*r**6+3072.0_cp*b**2*r**6)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp &
      &              )*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-4)=-a**5*b*(-145.0_cp*a**8*dn**3+750.0_cp*a**8*dn**2+4705.0_cp &
      &              *a**8*dn-21750.0_cp*a**8-2640.0_cp*a**6*b**2*dn**3+7200.0_cp*a**6*b**2 &
      &              *dn**2+86160.0_cp*a**6*b**2*dn-280800.0_cp*a**6*b**2+528.0_cp*a**6*dn**3 &
      &              *r**2-1440.0_cp*a**6*dn**2*r**2-17232.0_cp*a**6*dn*r**2+56160.0_cp* &
      &              a**6*r**2-8064.0_cp*a**4*b**4*dn**3+274176.0_cp*a**4*b**4*dn-604800.0_cp &
      &              *a**4*b**4+5376.0_cp*a**4*b**2*dn**3*r**2-182784.0_cp*a**4*b**2*dn* &
      &              r**2+403200.0_cp*a**4*b**2*r**2-576.0_cp*a**4*dn**3*r**4+19584.0_cp &
      &              *a**4*dn*r**4-43200.0_cp*a**4*r**4-3840.0_cp*a**2*b**6*dn**3-23040.0_cp &
      &              *a**2*b**6*dn**2+188160.0_cp*a**2*b**6*dn-253440.0_cp*a**2*b**6+5376.0_cp &
      &              *a**2*b**4*dn**3*r**2+32256.0_cp*a**2*b**4*dn**2*r**2-263424.0_cp*a**2 &
      &              *b**4*dn*r**2+354816.0_cp*a**2*b**4*r**2-1920.0_cp*a**2*b**2*dn**3* &
      &              r**4-11520.0_cp*a**2*b**2*dn**2*r**4+94080.0_cp*a**2*b**2*dn*r**4-126720.0_cp &
      &              *a**2*b**2*r**4+128.0_cp*a**2*dn**3*r**6+768.0_cp*a**2*dn**2*r**6-6272.0_cp &
      &              *a**2*dn*r**6+8448.0_cp*a**2*r**6+1280.0_cp*b**8*dn**3-7680.0_cp*b**8* &
      &              dn**2+14080.0_cp*b**8*dn-7680.0_cp*b**8-3072.0_cp*b**6*dn**3*r**2+18432.0_cp &
      &              *b**6*dn**2*r**2-33792.0_cp*b**6*dn*r**2+18432.0_cp*b**6*r**2+2304.0_cp &
      &              *b**4*dn**3*r**4-13824.0_cp*b**4*dn**2*r**4+25344.0_cp*b**4*dn*r**4 &
      &              -13824.0_cp*b**4*r**4-512.0_cp*b**2*dn**3*r**6+3072.0_cp*b**2*dn**2 &
      &              *r**6-5632.0_cp*b**2*dn*r**6+3072.0_cp*b**2*r**6)/(4096.0_cp*dn*(dn &
      &              -3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku-3)=-a**4*(-19.0_cp*a**10*dn**3+264.0_cp*a**10*dn**2+241.0_cp &
      &              *a**10*dn-4686.0_cp*a**10-1530.0_cp*a**8*b**2*dn**3+14580.0_cp*a**8 &
      &              *b**2*dn**2+26370.0_cp*a**8*b**2*dn-244620.0_cp*a**8*b**2+102.0_cp* &
      &              a**8*dn**3*r**2-972.0_cp*a**8*dn**2*r**2-1758.0_cp*a**8*dn*r**2+16308.0_cp &
      &              *a**8*r**2-13440.0_cp*a**6*b**4*dn**3+80640.0_cp*a**6*b**4*dn**2+255360.0_cp &
      &              *a**6*b**4*dn-1330560.0_cp*a**6*b**4+5376.0_cp*a**6*b**2*dn**3*r**2 &
      &              -32256.0_cp*a**6*b**2*dn**2*r**2-102144.0_cp*a**6*b**2*dn*r**2+532224.0_cp &
      &              *a**6*b**2*r**2-192.0_cp*a**6*dn**3*r**4+1152.0_cp*a**6*dn**2*r**4+3648.0_cp &
      &              *a**6*dn*r**4-19008.0_cp*a**6*r**4-26880.0_cp*a**4*b**6*dn**3+80640.0_cp &
      &              *a**4*b**6*dn**2+510720.0_cp*a**4*b**6*dn-1532160.0_cp*a**4*b**6+26880.0_cp &
      &              *a**4*b**4*dn**3*r**2-80640.0_cp*a**4*b**4*dn**2*r**2-510720.0_cp*a**4 &
      &              *b**4*dn*r**2+1532160.0_cp*a**4*b**4*r**2-5760.0_cp*a**4*b**2*dn**3 &
      &              *r**4+17280.0_cp*a**4*b**2*dn**2*r**4+109440.0_cp*a**4*b**2*dn*r**4 &
      &              -328320.0_cp*a**4*b**2*r**4+128.0_cp*a**4*dn**3*r**6-384.0_cp*a**4*dn**2 &
      &              *r**6-2432.0_cp*a**4*dn*r**6+7296.0_cp*a**4*r**6-11520.0_cp*a**2*b**8* &
      &              dn**3+218880.0_cp*a**2*b**8*dn-345600.0_cp*a**2*b**8+21504.0_cp*a** &
      &              2*b**6*dn**3*r**2-408576.0_cp*a**2*b**6*dn*r**2+645120.0_cp*a**2*b**6* &
      &              r**2-11520.0_cp*a**2*b**4*dn**3*r**4+218880.0_cp*a**2*b**4*dn*r**4-345600.0_cp &
      &              *a**2*b**4*r**4+1536.0_cp*a**2*b**2*dn**3*r**6-29184.0_cp*a**2*b**2 &
      &              *dn*r**6+46080.0_cp*a**2*b**2*r**6+512.0_cp*b**10*dn**3-3072.0_cp*b**10 &
      &              *dn**2+5632.0_cp*b**10*dn-3072.0_cp*b**10-1536.0_cp*b**8*dn**3*r** &
      &              2+9216.0_cp*b**8*dn**2*r**2-16896.0_cp*b**8*dn*r**2+9216.0_cp*b**8* &
      &              r**2+1536.0_cp*b**6*dn**3*r**4-9216.0_cp*b**6*dn**2*r**4+16896.0_cp &
      &              *b**6*dn*r**4-9216.0_cp*b**6*r**4-512.0_cp*b**4*dn**3*r**6+3072.0_cp &
      &              *b**4*dn**2*r**6-5632.0_cp*b**4*dn*r**6+3072.0_cp*b**4*r**6)/(8192.0_cp &
      &              *dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn &
      &              +3.0_cp))
      stencil(ku-2)=-3.0_cp*a**5*b*(-15.0_cp*a**8*dn**2+675.0_cp*a**8*dn-2790.0_cp &
      &              *a**8-480.0_cp*a**6*b**2*dn**2+12000.0_cp*a**6*b**2*dn-41280.0_cp*a**6 &
      &              *b**2+96.0_cp*a**6*dn**2*r**2-2400.0_cp*a**6*dn*r**2+8256.0_cp*a**6 &
      &              *r**2-2688.0_cp*a**4*b**4*dn**2+40320.0_cp*a**4*b**4*dn-110208.0_cp &
      &              *a**4*b**4+1792.0_cp*a**4*b**2*dn**2*r**2-26880.0_cp*a**4*b**2*dn*r**2 &
      &              +73472.0_cp*a**4*b**2*r**2-192.0_cp*a**4*dn**2*r**4+2880.0_cp*a**4*dn &
      &              *r**4-7872.0_cp*a**4*r**4-3840.0_cp*a**2*b**6*dn**2+34560.0_cp*a**2 &
      &              *b**6*dn-69120.0_cp*a**2*b**6+5376.0_cp*a**2*b**4*dn**2*r**2-48384.0_cp &
      &              *a**2*b**4*dn*r**2+96768.0_cp*a**2*b**4*r**2-1920.0_cp*a**2*b**2*dn**2 &
      &              *r**4+17280.0_cp*a**2*b**2*dn*r**4-34560.0_cp*a**2*b**2*r**4+128.0_cp &
      &              *a**2*dn**2*r**6-1152.0_cp*a**2*dn*r**6+2304.0_cp*a**2*r**6-1280.0_cp &
      &              *b**8*dn**2+6400.0_cp*b**8*dn-7680.0_cp*b**8+3072.0_cp*b**6*dn**2*r**2 &
      &              -15360.0_cp*b**6*dn*r**2+18432.0_cp*b**6*r**2-2304.0_cp*b**4*dn**2* &
      &              r**4+11520.0_cp*b**4*dn*r**4-13824.0_cp*b**4*r**4+512.0_cp*b**2*dn**2* &
      &              r**6-2560.0_cp*b**2*dn*r**6+3072.0_cp*b**2*r**6)/(4096.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku-1)=a**4*(-27.0_cp*a**10*dn**2-540.0_cp*a**10*dn+3663.0_cp &
      &              *a**10-1440.0_cp*a**8*b**2*dn**2-36000.0_cp*a**8*b**2*dn+207360.0_cp &
      &              *a**8*b**2+96.0_cp*a**8*dn**2*r**2+2400.0_cp*a**8*dn*r**2-13824.0_cp &
      &              *a**8*r**2-6720.0_cp*a**6*b**4*dn**2-268800.0_cp*a**6*b**4*dn+1270080.0_cp &
      &              *a**6*b**4+2688.0_cp*a**6*b**2*dn**2*r**2+107520.0_cp*a**6*b**2*dn* &
      &              r**2-508032.0_cp*a**6*b**2*r**2-96.0_cp*a**6*dn**2*r**4-3840.0_cp*a**6 &
      &              *dn*r**4+18144.0_cp*a**6*r**4-483840.0_cp*a**4*b**6*dn+1774080.0_cp &
      &              *a**4*b**6+483840.0_cp*a**4*b**4*dn*r**2-1774080.0_cp*a**4*b**4*r** &
      &              2-103680.0_cp*a**4*b**2*dn*r**4+380160.0_cp*a**4*b**2*r**4+2304.0_cp &
      &              *a**4*dn*r**6-8448.0_cp*a**4*r**6+11520.0_cp*a**2*b**8*dn**2-230400.0_cp &
      &              *a**2*b**8*dn+587520.0_cp*a**2*b**8-21504.0_cp*a**2*b**6*dn**2*r**2 &
      &              +430080.0_cp*a**2*b**6*dn*r**2-1096704.0_cp*a**2*b**6*r**2+11520.0_cp &
      &              *a**2*b**4*dn**2*r**4-230400.0_cp*a**2*b**4*dn*r**4+587520.0_cp*a** &
      &              2*b**4*r**4-1536.0_cp*a**2*b**2*dn**2*r**6+30720.0_cp*a**2*b**2*dn* &
      &              r**6-78336.0_cp*a**2*b**2*r**6+4096.0_cp*b**10*dn**2-20480.0_cp*b**10 &
      &              *dn+24576.0_cp*b**10-12288.0_cp*b**8*dn**2*r**2+61440.0_cp*b**8*dn* &
      &              r**2-73728.0_cp*b**8*r**2+12288.0_cp*b**6*dn**2*r**4-61440.0_cp*b** &
      &              6*dn*r**4+73728.0_cp*b**6*r**4-4096.0_cp*b**4*dn**2*r**6+20480.0_cp &
      &              *b**4*dn*r**6-24576.0_cp*b**4*r**6)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp &
      &              )*(dn-1.0_cp)*(dn+1.0_cp)*(dn+3.0_cp))
      stencil(ku)=a**5*b*(-45.0_cp*a**8*dn**2-225.0_cp*a**8*dn+2880.0_cp* &
      &              a**8-840.0_cp*a**6*b**2*dn**2-4200.0_cp*a**6*b**2*dn+45360.0_cp*a** &
      &              6*b**2+168.0_cp*a**6*dn**2*r**2+840.0_cp*a**6*dn*r**2-9072.0_cp*a** &
      &              6*r**2-3024.0_cp*a**4*b**4*dn**2-15120.0_cp*a**4*b**4*dn+133056.0_cp &
      &              *a**4*b**4+2016.0_cp*a**4*b**2*dn**2*r**2+10080.0_cp*a**4*b**2*dn*r**2 &
      &              -88704.0_cp*a**4*b**2*r**2-216.0_cp*a**4*dn**2*r**4-1080.0_cp*a**4*dn &
      &              *r**4+9504.0_cp*a**4*r**4-2880.0_cp*a**2*b**6*dn**2-14400.0_cp*a**2 &
      &              *b**6*dn+97920.0_cp*a**2*b**6+4032.0_cp*a**2*b**4*dn**2*r**2+20160.0_cp &
      &              *a**2*b**4*dn*r**2-137088.0_cp*a**2*b**4*r**2-1440.0_cp*a**2*b**2*dn**2 &
      &              *r**4-7200.0_cp*a**2*b**2*dn*r**4+48960.0_cp*a**2*b**2*r**4+96.0_cp &
      &              *a**2*dn**2*r**6+480.0_cp*a**2*dn*r**6-3264.0_cp*a**2*r**6-640.0_cp &
      &              *b**8*dn**2-3200.0_cp*b**8*dn+15360.0_cp*b**8+1536.0_cp*b**6*dn**2* &
      &              r**2+7680.0_cp*b**6*dn*r**2-36864.0_cp*b**6*r**2-1152.0_cp*b**4*dn**2* &
      &              r**4-5760.0_cp*b**4*dn*r**4+27648.0_cp*b**4*r**4+256.0_cp*b**2*dn** &
      &              2*r**6+1280.0_cp*b**2*dn*r**6-6144.0_cp*b**2*r**6)/(1024.0_cp*dn*(dn &
      &              -3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+1)=3.0_cp*a**4*(-3.0_cp*a**10*dn**2+177.0_cp*a**10-210.0_cp &
      &              *a**8*b**2*dn**2+10290.0_cp*a**8*b**2+14.0_cp*a**8*dn**2*r**2-686.0_cp &
      &              *a**8*r**2-1680.0_cp*a**6*b**4*dn**2+65520.0_cp*a**6*b**4+672.0_cp* &
      &              a**6*b**2*dn**2*r**2-26208.0_cp*a**6*b**2*r**2-24.0_cp*a**6*dn**2*r**4 &
      &              +936.0_cp*a**6*r**4-3360.0_cp*a**4*b**6*dn**2+97440.0_cp*a**4*b**6+3360.0_cp &
      &              *a**4*b**4*dn**2*r**2-97440.0_cp*a**4*b**4*r**2-720.0_cp*a**4*b**2*dn**2 &
      &              *r**4+20880.0_cp*a**4*b**2*r**4+16.0_cp*a**4*dn**2*r**6-464.0_cp*a**4* &
      &              r**6-1920.0_cp*a**2*b**8*dn**2+36480.0_cp*a**2*b**8+3584.0_cp*a**2* &
      &              b**6*dn**2*r**2-68096.0_cp*a**2*b**6*r**2-1920.0_cp*a**2*b**4*dn**2 &
      &              *r**4+36480.0_cp*a**2*b**4*r**4+256.0_cp*a**2*b**2*dn**2*r**6-4864.0_cp &
      &              *a**2*b**2*r**6-256.0_cp*b**10*dn**2+2304.0_cp*b**10+768.0_cp*b**8*dn**2 &
      &              *r**2-6912.0_cp*b**8*r**2-768.0_cp*b**6*dn**2*r**4+6912.0_cp*b**6*r**4 &
      &              +256.0_cp*b**4*dn**2*r**6-2304.0_cp*b**4*r**6)/(2048.0_cp*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+2)=a**5*b*(-45.0_cp*a**8*dn**2+225.0_cp*a**8*dn+2880.0_cp &
      &              *a**8-840.0_cp*a**6*b**2*dn**2+4200.0_cp*a**6*b**2*dn+45360.0_cp*a**6* &
      &              b**2+168.0_cp*a**6*dn**2*r**2-840.0_cp*a**6*dn*r**2-9072.0_cp*a**6* &
      &              r**2-3024.0_cp*a**4*b**4*dn**2+15120.0_cp*a**4*b**4*dn+133056.0_cp* &
      &              a**4*b**4+2016.0_cp*a**4*b**2*dn**2*r**2-10080.0_cp*a**4*b**2*dn*r**2- &
      &              88704.0_cp*a**4*b**2*r**2-216.0_cp*a**4*dn**2*r**4+1080.0_cp*a**4*dn &
      &              *r**4+9504.0_cp*a**4*r**4-2880.0_cp*a**2*b**6*dn**2+14400.0_cp*a**2 &
      &              *b**6*dn+97920.0_cp*a**2*b**6+4032.0_cp*a**2*b**4*dn**2*r**2-20160.0_cp &
      &              *a**2*b**4*dn*r**2-137088.0_cp*a**2*b**4*r**2-1440.0_cp*a**2*b**2*dn**2 &
      &              *r**4+7200.0_cp*a**2*b**2*dn*r**4+48960.0_cp*a**2*b**2*r**4+96.0_cp &
      &              *a**2*dn**2*r**6-480.0_cp*a**2*dn*r**6-3264.0_cp*a**2*r**6-640.0_cp &
      &              *b**8*dn**2+3200.0_cp*b**8*dn+15360.0_cp*b**8+1536.0_cp*b**6*dn**2* &
      &              r**2-7680.0_cp*b**6*dn*r**2-36864.0_cp*b**6*r**2-1152.0_cp*b**4*dn**2* &
      &              r**4+5760.0_cp*b**4*dn*r**4+27648.0_cp*b**4*r**4+256.0_cp*b**2*dn** &
      &              2*r**6-1280.0_cp*b**2*dn*r**6-6144.0_cp*b**2*r**6)/(1024.0_cp*dn*(dn &
      &              -3.0_cp)*(dn-2.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+3)=-a**4*(27.0_cp*a**10*dn**2-540.0_cp*a**10*dn-3663.0_cp &
      &              *a**10+1440.0_cp*a**8*b**2*dn**2-36000.0_cp*a**8*b**2*dn-207360.0_cp &
      &              *a**8*b**2-96.0_cp*a**8*dn**2*r**2+2400.0_cp*a**8*dn*r**2+13824.0_cp &
      &              *a**8*r**2+6720.0_cp*a**6*b**4*dn**2-268800.0_cp*a**6*b**4*dn-1270080.0_cp &
      &              *a**6*b**4-2688.0_cp*a**6*b**2*dn**2*r**2+107520.0_cp*a**6*b**2*dn* &
      &              r**2+508032.0_cp*a**6*b**2*r**2+96.0_cp*a**6*dn**2*r**4-3840.0_cp*a**6 &
      &              *dn*r**4-18144.0_cp*a**6*r**4-483840.0_cp*a**4*b**6*dn-1774080.0_cp &
      &              *a**4*b**6+483840.0_cp*a**4*b**4*dn*r**2+1774080.0_cp*a**4*b**4*r** &
      &              2-103680.0_cp*a**4*b**2*dn*r**4-380160.0_cp*a**4*b**2*r**4+2304.0_cp &
      &              *a**4*dn*r**6+8448.0_cp*a**4*r**6-11520.0_cp*a**2*b**8*dn**2-230400.0_cp &
      &              *a**2*b**8*dn-587520.0_cp*a**2*b**8+21504.0_cp*a**2*b**6*dn**2*r**2 &
      &              +430080.0_cp*a**2*b**6*dn*r**2+1096704.0_cp*a**2*b**6*r**2-11520.0_cp &
      &              *a**2*b**4*dn**2*r**4-230400.0_cp*a**2*b**4*dn*r**4-587520.0_cp*a** &
      &              2*b**4*r**4+1536.0_cp*a**2*b**2*dn**2*r**6+30720.0_cp*a**2*b**2*dn* &
      &              r**6+78336.0_cp*a**2*b**2*r**6-4096.0_cp*b**10*dn**2-20480.0_cp*b**10 &
      &              *dn-24576.0_cp*b**10+12288.0_cp*b**8*dn**2*r**2+61440.0_cp*b**8*dn* &
      &              r**2+73728.0_cp*b**8*r**2-12288.0_cp*b**6*dn**2*r**4-61440.0_cp*b** &
      &              6*dn*r**4-73728.0_cp*b**6*r**4+4096.0_cp*b**4*dn**2*r**6+20480.0_cp &
      &              *b**4*dn*r**6+24576.0_cp*b**4*r**6)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+4)=-3.0_cp*a**5*b*(-15.0_cp*a**8*dn**2-675.0_cp*a**8*dn-2790.0_cp &
      &              *a**8-480.0_cp*a**6*b**2*dn**2-12000.0_cp*a**6*b**2*dn-41280.0_cp*a**6 &
      &              *b**2+96.0_cp*a**6*dn**2*r**2+2400.0_cp*a**6*dn*r**2+8256.0_cp*a**6 &
      &              *r**2-2688.0_cp*a**4*b**4*dn**2-40320.0_cp*a**4*b**4*dn-110208.0_cp &
      &              *a**4*b**4+1792.0_cp*a**4*b**2*dn**2*r**2+26880.0_cp*a**4*b**2*dn*r**2 &
      &              +73472.0_cp*a**4*b**2*r**2-192.0_cp*a**4*dn**2*r**4-2880.0_cp*a**4*dn &
      &              *r**4-7872.0_cp*a**4*r**4-3840.0_cp*a**2*b**6*dn**2-34560.0_cp*a**2 &
      &              *b**6*dn-69120.0_cp*a**2*b**6+5376.0_cp*a**2*b**4*dn**2*r**2+48384.0_cp &
      &              *a**2*b**4*dn*r**2+96768.0_cp*a**2*b**4*r**2-1920.0_cp*a**2*b**2*dn**2 &
      &              *r**4-17280.0_cp*a**2*b**2*dn*r**4-34560.0_cp*a**2*b**2*r**4+128.0_cp &
      &              *a**2*dn**2*r**6+1152.0_cp*a**2*dn*r**6+2304.0_cp*a**2*r**6-1280.0_cp &
      &              *b**8*dn**2-6400.0_cp*b**8*dn-7680.0_cp*b**8+3072.0_cp*b**6*dn**2*r**2 &
      &              +15360.0_cp*b**6*dn*r**2+18432.0_cp*b**6*r**2-2304.0_cp*b**4*dn**2* &
      &              r**4-11520.0_cp*b**4*dn*r**4-13824.0_cp*b**4*r**4+512.0_cp*b**2*dn**2* &
      &              r**6+2560.0_cp*b**2*dn*r**6+3072.0_cp*b**2*r**6)/(4096.0_cp*dn*(dn-2.0_cp &
      &              )*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+5)=-a**4*(-19.0_cp*a**10*dn**3-264.0_cp*a**10*dn**2+241.0_cp &
      &              *a**10*dn+4686.0_cp*a**10-1530.0_cp*a**8*b**2*dn**3-14580.0_cp*a**8 &
      &              *b**2*dn**2+26370.0_cp*a**8*b**2*dn+244620.0_cp*a**8*b**2+102.0_cp* &
      &              a**8*dn**3*r**2+972.0_cp*a**8*dn**2*r**2-1758.0_cp*a**8*dn*r**2-16308.0_cp &
      &              *a**8*r**2-13440.0_cp*a**6*b**4*dn**3-80640.0_cp*a**6*b**4*dn**2+255360.0_cp &
      &              *a**6*b**4*dn+1330560.0_cp*a**6*b**4+5376.0_cp*a**6*b**2*dn**3*r**2 &
      &              +32256.0_cp*a**6*b**2*dn**2*r**2-102144.0_cp*a**6*b**2*dn*r**2-532224.0_cp &
      &              *a**6*b**2*r**2-192.0_cp*a**6*dn**3*r**4-1152.0_cp*a**6*dn**2*r**4+3648.0_cp &
      &              *a**6*dn*r**4+19008.0_cp*a**6*r**4-26880.0_cp*a**4*b**6*dn**3-80640.0_cp &
      &              *a**4*b**6*dn**2+510720.0_cp*a**4*b**6*dn+1532160.0_cp*a**4*b**6+26880.0_cp &
      &              *a**4*b**4*dn**3*r**2+80640.0_cp*a**4*b**4*dn**2*r**2-510720.0_cp*a**4 &
      &              *b**4*dn*r**2-1532160.0_cp*a**4*b**4*r**2-5760.0_cp*a**4*b**2*dn**3 &
      &              *r**4-17280.0_cp*a**4*b**2*dn**2*r**4+109440.0_cp*a**4*b**2*dn*r**4 &
      &              +328320.0_cp*a**4*b**2*r**4+128.0_cp*a**4*dn**3*r**6+384.0_cp*a**4*dn**2 &
      &              *r**6-2432.0_cp*a**4*dn*r**6-7296.0_cp*a**4*r**6-11520.0_cp*a**2*b**8* &
      &              dn**3+218880.0_cp*a**2*b**8*dn+345600.0_cp*a**2*b**8+21504.0_cp*a** &
      &              2*b**6*dn**3*r**2-408576.0_cp*a**2*b**6*dn*r**2-645120.0_cp*a**2*b**6* &
      &              r**2-11520.0_cp*a**2*b**4*dn**3*r**4+218880.0_cp*a**2*b**4*dn*r**4+345600.0_cp &
      &              *a**2*b**4*r**4+1536.0_cp*a**2*b**2*dn**3*r**6-29184.0_cp*a**2*b**2 &
      &              *dn*r**6-46080.0_cp*a**2*b**2*r**6+512.0_cp*b**10*dn**3+3072.0_cp*b**10 &
      &              *dn**2+5632.0_cp*b**10*dn+3072.0_cp*b**10-1536.0_cp*b**8*dn**3*r** &
      &              2-9216.0_cp*b**8*dn**2*r**2-16896.0_cp*b**8*dn*r**2-9216.0_cp*b**8* &
      &              r**2+1536.0_cp*b**6*dn**3*r**4+9216.0_cp*b**6*dn**2*r**4+16896.0_cp &
      &              *b**6*dn*r**4+9216.0_cp*b**6*r**4-512.0_cp*b**4*dn**3*r**6-3072.0_cp &
      &              *b**4*dn**2*r**6-5632.0_cp*b**4*dn*r**6-3072.0_cp*b**4*r**6)/(8192.0_cp &
      &              *dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn &
      &              +3.0_cp))
      stencil(ku+6)=-a**5*b*(-145.0_cp*a**8*dn**3-750.0_cp*a**8*dn**2+4705.0_cp &
      &              *a**8*dn+21750.0_cp*a**8-2640.0_cp*a**6*b**2*dn**3-7200.0_cp*a**6*b**2 &
      &              *dn**2+86160.0_cp*a**6*b**2*dn+280800.0_cp*a**6*b**2+528.0_cp*a**6*dn**3 &
      &              *r**2+1440.0_cp*a**6*dn**2*r**2-17232.0_cp*a**6*dn*r**2-56160.0_cp* &
      &              a**6*r**2-8064.0_cp*a**4*b**4*dn**3+274176.0_cp*a**4*b**4*dn+604800.0_cp &
      &              *a**4*b**4+5376.0_cp*a**4*b**2*dn**3*r**2-182784.0_cp*a**4*b**2*dn* &
      &              r**2-403200.0_cp*a**4*b**2*r**2-576.0_cp*a**4*dn**3*r**4+19584.0_cp &
      &              *a**4*dn*r**4+43200.0_cp*a**4*r**4-3840.0_cp*a**2*b**6*dn**3+23040.0_cp &
      &              *a**2*b**6*dn**2+188160.0_cp*a**2*b**6*dn+253440.0_cp*a**2*b**6+5376.0_cp &
      &              *a**2*b**4*dn**3*r**2-32256.0_cp*a**2*b**4*dn**2*r**2-263424.0_cp*a**2 &
      &              *b**4*dn*r**2-354816.0_cp*a**2*b**4*r**2-1920.0_cp*a**2*b**2*dn**3* &
      &              r**4+11520.0_cp*a**2*b**2*dn**2*r**4+94080.0_cp*a**2*b**2*dn*r**4+126720.0_cp &
      &              *a**2*b**2*r**4+128.0_cp*a**2*dn**3*r**6-768.0_cp*a**2*dn**2*r**6-6272.0_cp &
      &              *a**2*dn*r**6-8448.0_cp*a**2*r**6+1280.0_cp*b**8*dn**3+7680.0_cp*b**8* &
      &              dn**2+14080.0_cp*b**8*dn+7680.0_cp*b**8-3072.0_cp*b**6*dn**3*r**2-18432.0_cp &
      &              *b**6*dn**2*r**2-33792.0_cp*b**6*dn*r**2-18432.0_cp*b**6*r**2+2304.0_cp &
      &              *b**4*dn**3*r**4+13824.0_cp*b**4*dn**2*r**4+25344.0_cp*b**4*dn*r**4 &
      &              +13824.0_cp*b**4*r**4-512.0_cp*b**2*dn**3*r**6-3072.0_cp*b**2*dn**2 &
      &              *r**6-5632.0_cp*b**2*dn*r**6-3072.0_cp*b**2*r**6)/(4096.0_cp*dn*(dn &
      &              -3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku+7)=3.0_cp*a**6*(13.0_cp*a**8*dn**3+22.0_cp*a**8*dn**2-657.0_cp &
      &              *a**8*dn-2178.0_cp*a**8+720.0_cp*a**6*b**2*dn**3-720.0_cp*a**6*b**2 &
      &              *dn**2-38880.0_cp*a**6*b**2*dn-97920.0_cp*a**6*b**2-48.0_cp*a**6*dn**3 &
      &              *r**2+48.0_cp*a**6*dn**2*r**2+2592.0_cp*a**6*dn*r**2+6528.0_cp*a**6 &
      &              *r**2+3360.0_cp*a**4*b**4*dn**3-20160.0_cp*a**4*b**4*dn**2-231840.0_cp &
      &              *a**4*b**4*dn-423360.0_cp*a**4*b**4-1344.0_cp*a**4*b**2*dn**3*r**2+8064.0_cp &
      &              *a**4*b**2*dn**2*r**2+92736.0_cp*a**4*b**2*dn*r**2+169344.0_cp*a**4 &
      &              *b**2*r**2+48.0_cp*a**4*dn**3*r**4-288.0_cp*a**4*dn**2*r**4-3312.0_cp &
      &              *a**4*dn*r**4-6048.0_cp*a**4*r**4-53760.0_cp*a**2*b**6*dn**2-268800.0_cp &
      &              *a**2*b**6*dn-322560.0_cp*a**2*b**6+53760.0_cp*a**2*b**4*dn**2*r**2 &
      &              +268800.0_cp*a**2*b**4*dn*r**2+322560.0_cp*a**2*b**4*r**2-11520.0_cp &
      &              *a**2*b**2*dn**2*r**4-57600.0_cp*a**2*b**2*dn*r**4-69120.0_cp*a**2* &
      &              b**2*r**4+256.0_cp*a**2*dn**2*r**6+1280.0_cp*a**2*dn*r**6+1536.0_cp &
      &              *a**2*r**6-3840.0_cp*b**8*dn**3-23040.0_cp*b**8*dn**2-42240.0_cp*b**8* &
      &              dn-23040.0_cp*b**8+7168.0_cp*b**6*dn**3*r**2+43008.0_cp*b**6*dn**2* &
      &              r**2+78848.0_cp*b**6*dn*r**2+43008.0_cp*b**6*r**2-3840.0_cp*b**4*dn**3 &
      &              *r**4-23040.0_cp*b**4*dn**2*r**4-42240.0_cp*b**4*dn*r**4-23040.0_cp &
      &              *b**4*r**4+512.0_cp*b**2*dn**3*r**6+3072.0_cp*b**2*dn**2*r**6+5632.0_cp &
      &              *b**2*dn*r**6+3072.0_cp*b**2*r**6)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp &
      &              )*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+8)=a**7*b*(25.0_cp*a**6*dn**3-210.0_cp*a**6*dn**2-2575.0_cp &
      &              *a**6*dn-5460.0_cp*a**6+120.0_cp*a**4*b**2*dn**3-5040.0_cp*a**4*b** &
      &              2*dn**2-34680.0_cp*a**4*b**2*dn-55440.0_cp*a**4*b**2-24.0_cp*a**4*dn**3 &
      &              *r**2+1008.0_cp*a**4*dn**2*r**2+6936.0_cp*a**4*dn*r**2+11088.0_cp*a**4 &
      &              *r**2-1008.0_cp*a**2*b**4*dn**3-18144.0_cp*a**2*b**4*dn**2-71568.0_cp &
      &              *a**2*b**4*dn-78624.0_cp*a**2*b**4+672.0_cp*a**2*b**2*dn**3*r**2+12096.0_cp &
      &              *a**2*b**2*dn**2*r**2+47712.0_cp*a**2*b**2*dn*r**2+52416.0_cp*a**2* &
      &              b**2*r**2-72.0_cp*a**2*dn**3*r**4-1296.0_cp*a**2*dn**2*r**4-5112.0_cp &
      &              *a**2*dn*r**4-5616.0_cp*a**2*r**4-1920.0_cp*b**6*dn**3-11520.0_cp*b**6 &
      &              *dn**2-21120.0_cp*b**6*dn-11520.0_cp*b**6+2688.0_cp*b**4*dn**3*r**2 &
      &              +16128.0_cp*b**4*dn**2*r**2+29568.0_cp*b**4*dn*r**2+16128.0_cp*b**4 &
      &              *r**2-960.0_cp*b**2*dn**3*r**4-5760.0_cp*b**2*dn**2*r**4-10560.0_cp &
      &              *b**2*dn*r**4-5760.0_cp*b**2*r**4+64.0_cp*dn**3*r**6+384.0_cp*dn**2 &
      &              *r**6+704.0_cp*dn*r**6+384.0_cp*r**6)/(2048.0_cp*dn*(dn-3.0_cp)*(dn &
      &              -2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+9)=a**8*(a**6*dn**3-48.0_cp*a**6*dn**2-379.0_cp*a**6*dn-708.0_cp &
      &              *a**6-90.0_cp*a**4*b**2*dn**3-3240.0_cp*a**4*b**2*dn**2-17190.0_cp* &
      &              a**4*b**2*dn-24840.0_cp*a**4*b**2+6.0_cp*a**4*dn**3*r**2+216.0_cp*a**4 &
      &              *dn**2*r**2+1146.0_cp*a**4*dn*r**2+1656.0_cp*a**4*r**2-1680.0_cp*a**2* &
      &              b**4*dn**3-20160.0_cp*a**2*b**4*dn**2-68880.0_cp*a**2*b**4*dn-70560.0_cp &
      &              *a**2*b**4+672.0_cp*a**2*b**2*dn**3*r**2+8064.0_cp*a**2*b**2*dn**2* &
      &              r**2+27552.0_cp*a**2*b**2*dn*r**2+28224.0_cp*a**2*b**2*r**2-24.0_cp &
      &              *a**2*dn**3*r**4-288.0_cp*a**2*dn**2*r**4-984.0_cp*a**2*dn*r**4-1008.0_cp &
      &              *a**2*r**4-3360.0_cp*b**6*dn**3-20160.0_cp*b**6*dn**2-36960.0_cp*b**6* &
      &              dn-20160.0_cp*b**6+3360.0_cp*b**4*dn**3*r**2+20160.0_cp*b**4*dn**2* &
      &              r**2+36960.0_cp*b**4*dn*r**2+20160.0_cp*b**4*r**2-720.0_cp*b**2*dn**3* &
      &              r**4-4320.0_cp*b**2*dn**2*r**4-7920.0_cp*b**2*dn*r**4-4320.0_cp*b** &
      &              2*r**4+16.0_cp*dn**3*r**6+96.0_cp*dn**2*r**6+176.0_cp*dn*r**6+96.0_cp &
      &              *r**6)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+10)=-3.0_cp*a**9*b*(5.0_cp*a**4*dn**2+75.0_cp*a**4*dn+180.0_cp &
      &              *a**4+120.0_cp*a**2*b**2*dn**2+840.0_cp*a**2*b**2*dn+1200.0_cp*a**2 &
      &              *b**2-24.0_cp*a**2*dn**2*r**2-168.0_cp*a**2*dn*r**2-240.0_cp*a**2*r**2 &
      &              +336.0_cp*b**4*dn**2+1008.0_cp*b**4*dn+672.0_cp*b**4-224.0_cp*b**2*dn**2 &
      &              *r**2-672.0_cp*b**2*dn*r**2-448.0_cp*b**2*r**2+24.0_cp*dn**2*r**4+72.0_cp &
      &              *dn*r**4+48.0_cp*r**4)/(2048.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+11)=-a**10*(11.0_cp*a**4*dn**2+117.0_cp*a**4*dn+250.0_cp &
      &              *a**4+720.0_cp*a**2*b**2*dn**2+4320.0_cp*a**2*b**2*dn+5760.0_cp*a** &
      &              2*b**2-48.0_cp*a**2*dn**2*r**2-288.0_cp*a**2*dn*r**2-384.0_cp*a**2* &
      &              r**2+3360.0_cp*b**4*dn**2+10080.0_cp*b**4*dn+6720.0_cp*b**4-1344.0_cp &
      &              *b**2*dn**2*r**2-4032.0_cp*b**2*dn*r**2-2688.0_cp*b**2*r**2+48.0_cp &
      &              *dn**2*r**4+144.0_cp*dn*r**4+96.0_cp*r**4)/(16384.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+12)=a**11*b*(-25.0_cp*a**2*dn-85.0_cp*a**2-240.0_cp*b**2 &
      &              *dn-240.0_cp*b**2+48.0_cp*dn*r**2+48.0_cp*r**2)/(4096.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+13)=3.0_cp*a**12*(-a**2*dn-3.0_cp*a**2-30.0_cp*b**2*dn-30.0_cp &
      &              *b**2+2.0_cp*dn*r**2+2.0_cp*r**2)/(8192.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp &
      &              )*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+14)=-5.0_cp*a**13*b/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp &
      &              )*(dn-1.0_cp))
      stencil(ku+15)=-a**14/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              ))
      stencil(ku+16:len_stencil) = 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4rmult4hmult6
!------------------------------------------------------------------------------
   function intcheb4rmult4hmult8(a, b, r, n, len_stencil) result(stencil)
      !
      ! Fourth integral of r^4*(ro^2-r^2)^4 operator
      !

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
      stencil(ku-15)=a**16/(65536.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-14)=3.0_cp*a**15*b/(8192.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)* &
      &              (dn+3.0_cp))
      stencil(ku-13)=-a**14*(-2.0_cp*a**2*dn+5.0_cp*a**2-66.0_cp*b**2*dn+66.0_cp  &
      &              *b**2+4.0_cp*dn*r**2-4.0_cp*r**2)/(16384.0_cp*dn*(dn-1.0_cp)*&
      &              (dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-12)=-a**13*b*(-21.0_cp*a**2*dn+57.0_cp*a**2-220.0_cp*b**2*dn+    &
      &              220.0_cp*b**2+40.0_cp*dn*r**2-40.0_cp*r**2)/(8192.0_cp*dn*   &
      &              (dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-11)=3.0_cp*a**12*(2.0_cp*a**4*dn**2-15.0_cp*a**4*dn+27.0_cp*a**4+&
      &              132.0_cp*a**2*b**2*dn**2-660.0_cp*a**2*b**2*dn+792.0_cp*a**2*&
      &              b**2-8.0_cp*a**2*dn**2*r**2+40.0_cp*a**2*dn*r**2-48.0_cp*a**2&
      &              *r**2+660.0_cp*b**4*dn**2-1980.0_cp*b**4*dn+1320.0_cp*b**4-  &
      &              240.0_cp*b**2*dn**2*r**2+720.0_cp*b**2*dn*r**2-480.0_cp*b**2*&
      &              r**2+8.0_cp*dn**2*r**4-24.0_cp*dn*r**4+16.0_cp*r**4)/        &
      &              (16384.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*          &
      &              (dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-10)=a**11*b*(51.0_cp*a**4*dn**2-441.0_cp*a**4*dn+858.0_cp*a**4+  &
      &              1100.0_cp*a**2*b**2*dn**2-5940.0_cp*a**2*b**2*dn+7480.0_cp*  &
      &              a**2*b**2-200.0_cp*a**2*dn**2*r**2+1080.0_cp*a**2*dn*r**2-   &
      &              1360.0_cp*a**2*r**2+3168.0_cp*b**4*dn**2-9504.0_cp*b**4*dn+  &
      &              6336.0_cp*b**4-1920.0_cp*b**2*dn**2*r**2+5760.0_cp*b**2*dn*  &
      &              r**2-3840.0_cp*b**2*r**2+192.0_cp*dn**2*r**4-576.0_cp*dn*    &
      &              r**4+384.0_cp*r**4)/(8192.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*   &
      &              (dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
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
   function intcheb4rmult4hmult8laplrot(a, b, r, m, n, len_stencil) result(stencil)

      !-- Input variables
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
      real(cp), intent(in) :: r 
      integer,  intent(in) :: m 
      integer,  intent(in) :: n
      integer,  intent(in) :: len_stencil

      !-- Output variable 
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      real(cp) :: dm
      real(cp) :: dn
      integer :: ku

      dm = real(m, cp)
      dn = real(n, cp)
      ku = (len_stencil-1)/2

      stencil(1:ku-14) = 0.0_cp
      stencil(ku-13)=-a**14*(dm**2-dn**2-27.0_cp*dn-182.0_cp)/(16384.0_cp &
      &              *dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-12)=-a**13*b*(5.0_cp*dm**2-6.0_cp*dn**2-151.0_cp*dn-949.0_cp &
      &              )/(4096.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-11)=a**12*(-3.0_cp*a**2*dm**2*dn+9.0_cp*a**2*dm**2+5.0_cp &
      &              *a**2*dn**3+107.0_cp*a**2*dn**2+440.0_cp*a**2*dn-1488.0_cp*a**2-90.0_cp &
      &              *b**2*dm**2*dn+90.0_cp*b**2*dm**2+132.0_cp*b**2*dn**3+2948.0_cp*b** &
      &              2*dn**2+14872.0_cp*b**2*dn-17952.0_cp*b**2+8.0_cp*dm**2*dn*r**2-8.0_cp &
      &              *dm**2*r**2-8.0_cp*dn**3*r**2-178.0_cp*dn**2*r**2-898.0_cp*dn*r**2+1084.0_cp &
      &              *r**2)/(8192.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku-10)=a**11*b*(-25.0_cp*a**2*dm**2*dn+85.0_cp*a**2*dm**2+54.0_cp &
      &              *a**2*dn**3+1051.0_cp*a**2*dn**2+3712.0_cp*a**2*dn-14465.0_cp*a**2-240.0_cp &
      &              *b**2*dm**2*dn+240.0_cp*b**2*dm**2+440.0_cp*b**2*dn**3+9020.0_cp*b**2* &
      &              dn**2+41360.0_cp*b**2*dn-50820.0_cp*b**2+64.0_cp*dm**2*dn*r**2-64.0_cp &
      &              *dm**2*r**2-80.0_cp*dn**3*r**2-1634.0_cp*dn**2*r**2-7492.0_cp*dn*r**2+ &
      &              9206.0_cp*r**2)/(4096.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn &
      &              +3.0_cp))
      stencil(ku-9)=-a**10*(11.0_cp*a**4*dm**2*dn**2-117.0_cp*a**4*dm**2*dn &
      &              +250.0_cp*a**4*dm**2-43.0_cp*a**4*dn**4-614.0_cp*a**4*dn**3+345.0_cp &
      &              *a**4*dn**2+17900.0_cp*a**4*dn-39500.0_cp*a**4+720.0_cp*a**2*b**2*dm**2 &
      &              *dn**2-4320.0_cp*a**2*b**2*dm**2*dn+5760.0_cp*a**2*b**2*dm**2-2112.0_cp &
      &              *a**2*b**2*dn**4-32736.0_cp*a**2*b**2*dn**3-33792.0_cp*a**2*b**2*dn**2 &
      &              +722304.0_cp*a**2*b**2*dn-1013760.0_cp*a**2*b**2-64.0_cp*a**2*dm**2 &
      &              *dn**2*r**2+384.0_cp*a**2*dm**2*dn*r**2-512.0_cp*a**2*dm**2*r**2+128.0_cp &
      &              *a**2*dn**4*r**2+1976.0_cp*a**2*dn**3*r**2+2040.0_cp*a**2*dn**2*r** &
      &              2-43616.0_cp*a**2*dn*r**2+61216.0_cp*a**2*r**2+3360.0_cp*b**4*dm**2 &
      &              *dn**2-10080.0_cp*b**4*dm**2*dn+6720.0_cp*b**4*dm**2-7920.0_cp*b**4 &
      &              *dn**4-132000.0_cp*b**4*dn**3-314160.0_cp*b**4*dn**2+1985280.0_cp*b**4 &
      &              *dn-1531200.0_cp*b**4-1792.0_cp*b**2*dm**2*dn**2*r**2+5376.0_cp*b** &
      &              2*dm**2*dn*r**2-3584.0_cp*b**2*dm**2*r**2+2880.0_cp*b**2*dn**4*r**2 &
      &              +47808.0_cp*b**2*dn**3*r**2+113792.0_cp*b**2*dn**2*r**2-719232.0_cp &
      &              *b**2*dn*r**2+554752.0_cp*b**2*r**2+96.0_cp*dm**2*dn**2*r**4-288.0_cp &
      &              *dm**2*dn*r**4+192.0_cp*dm**2*r**4-96.0_cp*dn**4*r**4-1584.0_cp*dn**3* &
      &              r**4-3760.0_cp*dn**2*r**4+23808.0_cp*dn*r**4-18368.0_cp*r**4)/(16384.0_cp &
      &              *dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-8)=-a**9*b*(15.0_cp*a**4*dm**2*dn**2-225.0_cp*a**4*dm**2 &
      &              *dn+540.0_cp*a**4*dm**2-102.0_cp*a**4*dn**4-1233.0_cp*a**4*dn**3+2148.0_cp &
      &              *a**4*dn**2+35127.0_cp*a**4*dn-87480.0_cp*a**4+360.0_cp*a**2*b**2*dm**2 &
      &              *dn**2-2520.0_cp*a**2*b**2*dm**2*dn+3600.0_cp*a**2*b**2*dm**2-1540.0_cp &
      &              *a**2*b**2*dn**4-20790.0_cp*a**2*b**2*dn**3-5720.0_cp*a**2*b**2*dn**2+ &
      &              436590.0_cp*a**2*b**2*dn-659340.0_cp*a**2*b**2-96.0_cp*a**2*dm**2*dn**2 &
      &              *r**2+672.0_cp*a**2*dm**2*dn*r**2-960.0_cp*a**2*dm**2*r**2+280.0_cp &
      &              *a**2*dn**4*r**2+3765.0_cp*a**2*dn**3*r**2+1034.0_cp*a**2*dn**2*r** &
      &              2-79089.0_cp*a**2*dn*r**2+119442.0_cp*a**2*r**2+1008.0_cp*b**4*dm** &
      &              2*dn**2-3024.0_cp*b**4*dm**2*dn+2016.0_cp*b**4*dm**2-3168.0_cp*b**4 &
      &              *dn**4-46992.0_cp*b**4*dn**3-88704.0_cp*b**4*dn**2+642576.0_cp*b**4 &
      &              *dn-503712.0_cp*b**4-896.0_cp*b**2*dm**2*dn**2*r**2+2688.0_cp*b**2*dm**2 &
      &              *dn*r**2-1792.0_cp*b**2*dm**2*r**2+1920.0_cp*b**2*dn**4*r**2+28368.0_cp &
      &              *b**2*dn**3*r**2+53536.0_cp*b**2*dn**2*r**2-387984.0_cp*b**2*dn*r** &
      &              2+304160.0_cp*b**2*r**2+144.0_cp*dm**2*dn**2*r**4-432.0_cp*dm**2*dn &
      &              *r**4+288.0_cp*dm**2*r**4-192.0_cp*dn**4*r**4-2820.0_cp*dn**3*r**4-5304.0_cp &
      &              *dn**2*r**4+38532.0_cp*dn*r**4-30216.0_cp*r**4)/(2048.0_cp*dn*(dn-2.0_cp &
      &              )*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-7)=a**8*(a**6*dm**2*dn**3+48.0_cp*a**6*dm**2*dn**2-379.0_cp &
      &              *a**6*dm**2*dn+708.0_cp*a**6*dm**2+25.0_cp*a**6*dn**5+124.0_cp*a**6 &
      &              *dn**4-2085.0_cp*a**6*dn**3-2932.0_cp*a**6*dn**2+57440.0_cp*a**6*dn &
      &              -109056.0_cp*a**6-90.0_cp*a**4*b**2*dm**2*dn**3+3240.0_cp*a**4*b**2 &
      &              *dm**2*dn**2-17190.0_cp*a**4*b**2*dm**2*dn+24840.0_cp*a**4*b**2*dm**2+ &
      &              1716.0_cp*a**4*b**2*dn**5+11616.0_cp*a**4*b**2*dn**4-106788.0_cp*a**4* &
      &              b**2*dn**3-304128.0_cp*a**4*b**2*dn**2+2802096.0_cp*a**4*b**2*dn-4143744.0_cp &
      &              *a**4*b**2+8.0_cp*a**4*dm**2*dn**3*r**2-288.0_cp*a**4*dm**2*dn**2*r**2 &
      &              +1528.0_cp*a**4*dm**2*dn*r**2-2208.0_cp*a**4*dm**2*r**2-104.0_cp*a**4* &
      &              dn**5*r**2-700.0_cp*a**4*dn**4*r**2+6448.0_cp*a**4*dn**3*r**2+18364.0_cp &
      &              *a**4*dn**2*r**2-169208.0_cp*a**4*dn*r**2+250224.0_cp*a**4*r**2-1680.0_cp &
      &              *a**2*b**4*dm**2*dn**3+20160.0_cp*a**2*b**4*dm**2*dn**2-68880.0_cp* &
      &              a**2*b**4*dm**2*dn+70560.0_cp*a**2*b**4*dm**2+11880.0_cp*a**2*b**4*dn**5 &
      &              +100320.0_cp*a**2*b**4*dn**4-492360.0_cp*a**2*b**4*dn**3-2476320.0_cp &
      &              *a**2*b**4*dn**2+12708960.0_cp*a**2*b**4*dn-13559040.0_cp*a**2*b**4 &
      &              +896.0_cp*a**2*b**2*dm**2*dn**3*r**2-10752.0_cp*a**2*b**2*dm**2*dn**2* &
      &              r**2+36736.0_cp*a**2*b**2*dm**2*dn*r**2-37632.0_cp*a**2*b**2*dm**2* &
      &              r**2-4320.0_cp*a**2*b**2*dn**5*r**2-36288.0_cp*a**2*b**2*dn**4*r**2 &
      &              +178400.0_cp*a**2*b**2*dn**3*r**2+897024.0_cp*a**2*b**2*dn**2*r**2-4604480.0_cp &
      &              *a**2*b**2*dn*r**2+4912512.0_cp*a**2*b**2*r**2-48.0_cp*a**2*dm**2*dn**3 &
      &              *r**4+576.0_cp*a**2*dm**2*dn**2*r**4-1968.0_cp*a**2*dm**2*dn*r**4+2016.0_cp &
      &              *a**2*dm**2*r**4+144.0_cp*a**2*dn**5*r**4+1200.0_cp*a**2*dn**4*r**4 &
      &              -5920.0_cp*a**2*dn**3*r**4-29664.0_cp*a**2*dn**2*r**4+152416.0_cp*a**2 &
      &              *dn*r**4-162624.0_cp*a**2*r**4-3360.0_cp*b**6*dm**2*dn**3+20160.0_cp &
      &              *b**6*dm**2*dn**2-36960.0_cp*b**6*dm**2*dn+20160.0_cp*b**6*dm**2+14784.0_cp &
      &              *b**6*dn**5+147840.0_cp*b**6*dn**4-310464.0_cp*b**6*dn**3-3163776.0_cp &
      &              *b**6*dn**2+8988672.0_cp*b**6*dn-5677056.0_cp*b**6+4480.0_cp*b**4*dm**2 &
      &              *dn**3*r**2-26880.0_cp*b**4*dm**2*dn**2*r**2+49280.0_cp*b**4*dm**2*dn &
      &              *r**2-26880.0_cp*b**4*dm**2*r**2-13440.0_cp*b**4*dn**5*r**2-133728.0_cp &
      &              *b**4*dn**4*r**2+281344.0_cp*b**4*dn**3*r**2+2864736.0_cp*b**4*dn** &
      &              2*r**2-8141056.0_cp*b**4*dn*r**2+5142144.0_cp*b**4*r**2-1440.0_cp*b**2 &
      &              *dm**2*dn**3*r**4+8640.0_cp*b**2*dm**2*dn**2*r**4-15840.0_cp*b**2*dm**2 &
      &              *dn*r**4+8640.0_cp*b**2*dm**2*r**4+2688.0_cp*b**2*dn**5*r**4+26544.0_cp &
      &              *b**2*dn**4*r**4-56160.0_cp*b**2*dn**3*r**4-568560.0_cp*b**2*dn**2* &
      &              r**4+1617312.0_cp*b**2*dn*r**4-1021824.0_cp*b**2*r**4+64.0_cp*dm**2 &
      &              *dn**3*r**6-384.0_cp*dm**2*dn**2*r**6+704.0_cp*dm**2*dn*r**6-384.0_cp &
      &              *dm**2*r**6-64.0_cp*dn**5*r**6-624.0_cp*dn**4*r**6+1344.0_cp*dn**3* &
      &              r**6+13296.0_cp*dn**2*r**6-37952.0_cp*dn*r**6+24000.0_cp*r**6)/(4096.0_cp &
      &              *dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn &
      &              +3.0_cp))
      stencil(ku-6)=a**7*b*(25.0_cp*a**6*dm**2*dn**3+210.0_cp*a**6*dm**2*dn**2 &
      &              -2575.0_cp*a**6*dm**2*dn+5460.0_cp*a**6*dm**2+198.0_cp*a**6*dn**5+385.0_cp &
      &              *a**6*dn**4-17303.0_cp*a**6*dn**3+2387.0_cp*a**6*dn**2+398189.0_cp* &
      &              a**6*dn-861168.0_cp*a**6+120.0_cp*a**4*b**2*dm**2*dn**3+5040.0_cp*a**4 &
      &              *b**2*dm**2*dn**2-34680.0_cp*a**4*b**2*dm**2*dn+55440.0_cp*a**4*b** &
      &              2*dm**2+4180.0_cp*a**4*b**2*dn**5+17710.0_cp*a**4*b**2*dn**4-274450.0_cp &
      &              *a**4*b**2*dn**3-327250.0_cp*a**4*b**2*dn**2+5827470.0_cp*a**4*b**2 &
      &              *dn-9577260.0_cp*a**4*b**2-32.0_cp*a**4*dm**2*dn**3*r**2-1344.0_cp* &
      &              a**4*dm**2*dn**2*r**2+9248.0_cp*a**4*dm**2*dn*r**2-14784.0_cp*a**4*dm**2 &
      &              *r**2-760.0_cp*a**4*dn**5*r**2-3199.0_cp*a**4*dn**4*r**2+49721.0_cp &
      &              *a**4*dn**3*r**2+59269.0_cp*a**4*dn**2*r**2-1055689.0_cp*a**4*dn*r**2+ &
      &              1734978.0_cp*a**4*r**2-1008.0_cp*a**2*b**4*dm**2*dn**3+18144.0_cp*a**2 &
      &              *b**4*dm**2*dn**2-71568.0_cp*a**2*b**4*dm**2*dn+78624.0_cp*a**2*b** &
      &              4*dm**2+15840.0_cp*a**2*b**4*dn**5+99792.0_cp*a**2*b**4*dn**4-704880.0_cp &
      &              *a**2*b**4*dn**3-2073456.0_cp*a**2*b**4*dn**2+13931280.0_cp*a**2*b**4* &
      &              dn-16033248.0_cp*a**2*b**4+896.0_cp*a**2*b**2*dm**2*dn**3*r**2-16128.0_cp &
      &              *a**2*b**2*dm**2*dn**2*r**2+63616.0_cp*a**2*b**2*dm**2*dn*r**2-69888.0_cp &
      &              *a**2*b**2*dm**2*r**2-9600.0_cp*a**2*b**2*dn**5*r**2-60144.0_cp*a** &
      &              2*b**2*dn**4*r**2+425744.0_cp*a**2*b**2*dn**3*r**2+1251600.0_cp*a** &
      &              2*b**2*dn**2*r**2-8412176.0_cp*a**2*b**2*dn*r**2+9681504.0_cp*a**2* &
      &              b**2*r**2-144.0_cp*a**2*dm**2*dn**3*r**4+2592.0_cp*a**2*dm**2*dn**2 &
      &              *r**4-10224.0_cp*a**2*dm**2*dn*r**4+11232.0_cp*a**2*dm**2*r**4+960.0_cp &
      &              *a**2*dn**5*r**4+5964.0_cp*a**2*dn**4*r**4-42372.0_cp*a**2*dn**3*r**4- &
      &              124116.0_cp*a**2*dn**2*r**4+835428.0_cp*a**2*dn*r**4-961560.0_cp*a**2* &
      &              r**4-1920.0_cp*b**6*dm**2*dn**3+11520.0_cp*b**6*dm**2*dn**2-21120.0_cp &
      &              *b**6*dm**2*dn+11520.0_cp*b**6*dm**2+12672.0_cp*b**6*dn**5+103488.0_cp &
      &              *b**6*dn**4-302016.0_cp*b**6*dn**3-1915584.0_cp*b**6*dn**2+5915712.0_cp &
      &              *b**6*dn-3814272.0_cp*b**6+3584.0_cp*b**4*dm**2*dn**3*r**2-21504.0_cp &
      &              *b**4*dm**2*dn**2*r**2+39424.0_cp*b**4*dm**2*dn*r**2-21504.0_cp*b** &
      &              4*dm**2*r**2-16128.0_cp*b**4*dn**5*r**2-131040.0_cp*b**4*dn**4*r**2 &
      &              +383264.0_cp*b**4*dn**3*r**2+2427936.0_cp*b**4*dn**2*r**2-7501088.0_cp &
      &              *b**4*dn*r**2+4837056.0_cp*b**4*r**2-1920.0_cp*b**2*dm**2*dn**3*r** &
      &              4+11520.0_cp*b**2*dm**2*dn**2*r**4-21120.0_cp*b**2*dm**2*dn*r**4+11520.0_cp &
      &              *b**2*dm**2*r**4+5376.0_cp*b**2*dn**5*r**4+43344.0_cp*b**2*dn**4*r**4- &
      &              127408.0_cp*b**2*dn**3*r**4-802992.0_cp*b**2*dn**2*r**4+2484016.0_cp &
      &              *b**2*dn*r**4-1602336.0_cp*b**2*r**4+256.0_cp*dm**2*dn**3*r**6-1536.0_cp &
      &              *dm**2*dn**2*r**6+2816.0_cp*dm**2*dn*r**6-1536.0_cp*dm**2*r**6-384.0_cp &
      &              *dn**5*r**6-3056.0_cp*dn**4*r**6+9104.0_cp*dn**3*r**6+56336.0_cp*dn**2 &
      &              *r**6-174992.0_cp*dn*r**6+112992.0_cp*r**6)/(2048.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-5)=-a**6*(-39.0_cp*a**8*dm**2*dn**3+66.0_cp*a**8*dm**2*dn**2 &
      &              +1971.0_cp*a**8*dm**2*dn-6534.0_cp*a**8*dm**2-121.0_cp*a**8*dn**5+671.0_cp &
      &              *a**8*dn**4+13849.0_cp*a**8*dn**3-43571.0_cp*a**8*dn**2-293040.0_cp &
      &              *a**8*dn+991188.0_cp*a**8-2160.0_cp*a**6*b**2*dm**2*dn**3-2160.0_cp &
      &              *a**6*b**2*dm**2*dn**2+116640.0_cp*a**6*b**2*dm**2*dn-293760.0_cp*a**6 &
      &              *b**2*dm**2-10560.0_cp*a**6*b**2*dn**5+19360.0_cp*a**6*b**2*dn**4+955680.0_cp &
      &              *a**6*b**2*dn**3-1640320.0_cp*a**6*b**2*dn**2-18437760.0_cp*a**6*b**2* &
      &              dn+47646720.0_cp*a**6*b**2+192.0_cp*a**6*dm**2*dn**3*r**2+192.0_cp* &
      &              a**6*dm**2*dn**2*r**2-10368.0_cp*a**6*dm**2*dn*r**2+26112.0_cp*a**6 &
      &              *dm**2*r**2+640.0_cp*a**6*dn**5*r**2-1176.0_cp*a**6*dn**4*r**2-57712.0_cp &
      &              *a**6*dn**3*r**2+99048.0_cp*a**6*dn**2*r**2+1113408.0_cp*a**6*dn*r**2- &
      &              2877216.0_cp*a**6*r**2-10080.0_cp*a**4*b**4*dm**2*dn**3-60480.0_cp* &
      &              a**4*b**4*dm**2*dn**2+695520.0_cp*a**4*b**4*dm**2*dn-1270080.0_cp*a**4 &
      &              *b**4*dm**2-102960.0_cp*a**4*b**4*dn**5-134640.0_cp*a**4*b**4*dn**4 &
      &              +6977520.0_cp*a**4*b**4*dn**3-2051280.0_cp*a**4*b**4*dn**2-120985920.0_cp &
      &              *a**4*b**4*dn+228951360.0_cp*a**4*b**4+5376.0_cp*a**4*b**2*dm**2*dn**3 &
      &              *r**2+32256.0_cp*a**4*b**2*dm**2*dn**2*r**2-370944.0_cp*a**4*b**2*dm**2 &
      &              *dn*r**2+677376.0_cp*a**4*b**2*dm**2*r**2+37440.0_cp*a**4*b**2*dn** &
      &              5*r**2+48384.0_cp*a**4*b**2*dn**4*r**2-2528448.0_cp*a**4*b**2*dn**3 &
      &              *r**2+744192.0_cp*a**4*b**2*dn**2*r**2+43834752.0_cp*a**4*b**2*dn*r**2 &
      &              -82950912.0_cp*a**4*b**2*r**2-288.0_cp*a**4*dm**2*dn**3*r**4-1728.0_cp &
      &              *a**4*dm**2*dn**2*r**4+19872.0_cp*a**4*dm**2*dn*r**4-36288.0_cp*a** &
      &              4*dm**2*r**4-1248.0_cp*a**4*dn**5*r**4-1584.0_cp*a**4*dn**4*r**4+83808.0_cp &
      &              *a**4*dn**3*r**4-24912.0_cp*a**4*dn**2*r**4-1450944.0_cp*a**4*dn*r**4+ &
      &              2745792.0_cp*a**4*r**4-161280.0_cp*a**2*b**6*dm**2*dn**2+806400.0_cp &
      &              *a**2*b**6*dm**2*dn-967680.0_cp*a**2*b**6*dm**2-236544.0_cp*a**2*b**6* &
      &              dn**5-946176.0_cp*a**2*b**6*dn**4+10881024.0_cp*a**2*b**6*dn**3+14429184.0_cp &
      &              *a**2*b**6*dn**2-167473152.0_cp*a**2*b**6*dn+212889600.0_cp*a**2*b**6+ &
      &              215040.0_cp*a**2*b**4*dm**2*dn**2*r**2-1075200.0_cp*a**2*b**4*dm**2 &
      &              *dn*r**2+1290240.0_cp*a**2*b**4*dm**2*r**2+215040.0_cp*a**2*b**4*dn**5 &
      &              *r**2+854784.0_cp*a**2*b**4*dn**4*r**2-9859584.0_cp*a**2*b**4*dn**3 &
      &              *r**2-13058304.0_cp*a**2*b**4*dn**2*r**2+151689216.0_cp*a**2*b**4*dn &
      &              *r**2-192826368.0_cp*a**2*b**4*r**2-69120.0_cp*a**2*b**2*dm**2*dn** &
      &              2*r**4+345600.0_cp*a**2*b**2*dm**2*dn*r**4-414720.0_cp*a**2*b**2*dm**2 &
      &              *r**4-43008.0_cp*a**2*b**2*dn**5*r**4-169344.0_cp*a**2*b**2*dn**4*r**4 &
      &              +1962240.0_cp*a**2*b**2*dn**3*r**4+2586240.0_cp*a**2*b**2*dn**2*r** &
      &              4-30131712.0_cp*a**2*b**2*dn*r**4+38306304.0_cp*a**2*b**2*r**4+3072.0_cp &
      &              *a**2*dm**2*dn**2*r**6-15360.0_cp*a**2*dm**2*dn*r**6+18432.0_cp*a** &
      &              2*dm**2*r**6+1024.0_cp*a**2*dn**5*r**6+3968.0_cp*a**2*dn**4*r**6-46336.0_cp &
      &              *a**2*dn**3*r**6-60032.0_cp*a**2*dn**2*r**6+706560.0_cp*a**2*dn*r** &
      &              6-898560.0_cp*a**2*r**6+11520.0_cp*b**8*dm**2*dn**3-69120.0_cp*b**8 &
      &              *dm**2*dn**2+126720.0_cp*b**8*dm**2*dn-69120.0_cp*b**8*dm**2-126720.0_cp &
      &              *b**8*dn**5-802560.0_cp*b**8*dn**4+3168000.0_cp*b**8*dn**3+12460800.0_cp &
      &              *b**8*dn**2-43591680.0_cp*b**8*dn+28892160.0_cp*b**8-28672.0_cp*b** &
      &              6*dm**2*dn**3*r**2+172032.0_cp*b**6*dm**2*dn**2*r**2-315392.0_cp*b**6* &
      &              dm**2*dn*r**2+172032.0_cp*b**6*dm**2*r**2+215040.0_cp*b**6*dn**5*r**2+ &
      &              1354752.0_cp*b**6*dn**4*r**2-5361664.0_cp*b**6*dn**3*r**2-21052416.0_cp &
      &              *b**6*dn**2*r**2+73701376.0_cp*b**6*dn*r**2-48857088.0_cp*b**6*r**2 &
      &              +23040.0_cp*b**4*dm**2*dn**3*r**4-138240.0_cp*b**4*dm**2*dn**2*r**4 &
      &              +253440.0_cp*b**4*dm**2*dn*r**4-138240.0_cp*b**4*dm**2*r**4-107520.0_cp &
      &              *b**4*dn**5*r**4-672000.0_cp*b**4*dn**4*r**4+2672640.0_cp*b**4*dn** &
      &              3*r**4+10440960.0_cp*b**4*dn**2*r**4-36618240.0_cp*b**4*dn*r**4+24284160.0_cp &
      &              *b**4*r**4-6144.0_cp*b**2*dm**2*dn**3*r**6+36864.0_cp*b**2*dm**2*dn**2 &
      &              *r**6-67584.0_cp*b**2*dm**2*dn*r**6+36864.0_cp*b**2*dm**2*r**6+15360.0_cp &
      &              *b**2*dn**5*r**6+94720.0_cp*b**2*dn**4*r**6-380928.0_cp*b**2*dn**3* &
      &              r**6-1464832.0_cp*b**2*dn**2*r**6+5164032.0_cp*b**2*dn*r**6-3428352.0_cp &
      &              *b**2*r**6+256.0_cp*dm**2*dn**3*r**8-1536.0_cp*dm**2*dn**2*r**8+2816.0_cp &
      &              *dm**2*dn*r**8-1536.0_cp*dm**2*r**8-256.0_cp*dn**5*r**8-1536.0_cp*dn**4 &
      &              *r**8+6400.0_cp*dn**3*r**8+23040.0_cp*dn**2*r**8-82944.0_cp*dn*r**8 &
      &              +55296.0_cp*r**8)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-4)=-a**5*b*(-145.0_cp*a**8*dm**2*dn**3+750.0_cp*a**8*dm**2* &
      &              dn**2+4705.0_cp*a**8*dm**2*dn-21750.0_cp*a**8*dm**2-330.0_cp*a**8*dn**5 &
      &              +4675.0_cp*a**8*dn**4+40205.0_cp*a**8*dn**3-234025.0_cp*a**8*dn**2-709115.0_cp &
      &              *a**8*dn+3377550.0_cp*a**8-2640.0_cp*a**6*b**2*dm**2*dn**3+7200.0_cp &
      &              *a**6*b**2*dm**2*dn**2+86160.0_cp*a**6*b**2*dm**2*dn-280800.0_cp*a**6* &
      &              b**2*dm**2-9240.0_cp*a**6*b**2*dn**5+69300.0_cp*a**6*b**2*dn**4+867900.0_cp &
      &              *a**6*b**2*dn**3-3217500.0_cp*a**6*b**2*dn**2-13926660.0_cp*a**6*b**2* &
      &              dn+47104200.0_cp*a**6*b**2+704.0_cp*a**6*dm**2*dn**3*r**2-1920.0_cp &
      &              *a**6*dm**2*dn**2*r**2-22976.0_cp*a**6*dm**2*dn*r**2+74880.0_cp*a** &
      &              6*dm**2*r**2+1680.0_cp*a**6*dn**5*r**2-12570.0_cp*a**6*dn**4*r**2-157246.0_cp &
      &              *a**6*dn**3*r**2+582870.0_cp*a**6*dn**2*r**2+2522974.0_cp*a**6*dn*r**2 &
      &              -8533260.0_cp*a**6*r**2-8064.0_cp*a**4*b**4*dm**2*dn**3+274176.0_cp &
      &              *a**4*b**4*dm**2*dn-604800.0_cp*a**4*b**4*dm**2-50688.0_cp*a**4*b** &
      &              4*dn**5+126720.0_cp*a**4*b**4*dn**4+3497472.0_cp*a**4*b**4*dn**3-6526080.0_cp &
      &              *a**4*b**4*dn**2-49547520.0_cp*a**4*b**4*dn+114998400.0_cp*a**4*b** &
      &              4+7168.0_cp*a**4*b**2*dm**2*dn**3*r**2-243712.0_cp*a**4*b**2*dm**2*dn &
      &              *r**2+537600.0_cp*a**4*b**2*dm**2*r**2+30720.0_cp*a**4*b**2*dn**5*r**2 &
      &              -76800.0_cp*a**4*b**2*dn**4*r**2-2112512.0_cp*a**4*b**2*dn**3*r**2+3941760.0_cp &
      &              *a**4*b**2*dn**2*r**2+29919488.0_cp*a**4*b**2*dn*r**2-69440640.0_cp &
      &              *a**4*b**2*r**2-1152.0_cp*a**4*dm**2*dn**3*r**4+39168.0_cp*a**4*dm**2* &
      &              dn*r**4-86400.0_cp*a**4*dm**2*r**4-3072.0_cp*a**4*dn**5*r**4+7680.0_cp &
      &              *a**4*dn**4*r**4+210048.0_cp*a**4*dn**3*r**4-392160.0_cp*a**4*dn**2 &
      &              *r**4-2971200.0_cp*a**4*dn*r**4+6896160.0_cp*a**4*r**4-3840.0_cp*a**2* &
      &              b**6*dm**2*dn**3-23040.0_cp*a**2*b**6*dm**2*dn**2+188160.0_cp*a**2* &
      &              b**6*dm**2*dn-253440.0_cp*a**2*b**6*dm**2-76032.0_cp*a**2*b**6*dn** &
      &              5-105600.0_cp*a**2*b**6*dn**4+3501696.0_cp*a**2*b**6*dn**3-274560.0_cp &
      &              *a**2*b**6*dn**2-42252672.0_cp*a**2*b**6*dn+61712640.0_cp*a**2*b**6 &
      &              +7168.0_cp*a**2*b**4*dm**2*dn**3*r**2+43008.0_cp*a**2*b**4*dm**2*dn**2 &
      &              *r**2-351232.0_cp*a**2*b**4*dm**2*dn*r**2+473088.0_cp*a**2*b**4*dm**2* &
      &              r**2+96768.0_cp*a**2*b**4*dn**5*r**2+133056.0_cp*a**2*b**4*dn**4*r**2- &
      &              4442816.0_cp*a**2*b**4*dn**3*r**2+353472.0_cp*a**2*b**4*dn**2*r**2+53579456.0_cp &
      &              *a**2*b**4*dn*r**2-78255744.0_cp*a**2*b**4*r**2-3840.0_cp*a**2*b**2 &
      &              *dm**2*dn**3*r**4-23040.0_cp*a**2*b**2*dm**2*dn**2*r**4+188160.0_cp &
      &              *a**2*b**2*dm**2*dn*r**4-253440.0_cp*a**2*b**2*dm**2*r**4-32256.0_cp &
      &              *a**2*b**2*dn**5*r**4-43680.0_cp*a**2*b**2*dn**4*r**4+1473568.0_cp* &
      &              a**2*b**2*dn**3*r**4-122400.0_cp*a**2*b**2*dn**2*r**4-17740576.0_cp &
      &              *a**2*b**2*dn*r**4+25913280.0_cp*a**2*b**2*r**4+512.0_cp*a**2*dm**2 &
      &              *dn**3*r**6+3072.0_cp*a**2*dm**2*dn**2*r**6-25088.0_cp*a**2*dm**2*dn &
      &              *r**6+33792.0_cp*a**2*dm**2*r**6+2304.0_cp*a**2*dn**5*r**6+3040.0_cp &
      &              *a**2*dn**4*r**6-104288.0_cp*a**2*dn**3*r**6+9824.0_cp*a**2*dn**2*r**6 &
      &              +1248608.0_cp*a**2*dn*r**6-1824576.0_cp*a**2*r**6+1280.0_cp*b**8*dm**2 &
      &              *dn**3-7680.0_cp*b**8*dm**2*dn**2+14080.0_cp*b**8*dm**2*dn-7680.0_cp &
      &              *b**8*dm**2-28160.0_cp*b**8*dn**5-126720.0_cp*b**8*dn**4+689920.0_cp &
      &              *b**8*dn**3+1562880.0_cp*b**8*dn**2-6744320.0_cp*b**8*dn+4646400.0_cp &
      &              *b**8-4096.0_cp*b**6*dm**2*dn**3*r**2+24576.0_cp*b**6*dm**2*dn**2*r**2 &
      &              -45056.0_cp*b**6*dm**2*dn*r**2+24576.0_cp*b**6*dm**2*r**2+61440.0_cp &
      &              *b**6*dn**5*r**2+274944.0_cp*b**6*dn**4*r**2-1501696.0_cp*b**6*dn** &
      &              3*r**2-3393024.0_cp*b**6*dn**2*r**2+14662144.0_cp*b**6*dn*r**2-10103808.0_cp &
      &              *b**6*r**2+4608.0_cp*b**4*dm**2*dn**3*r**4-27648.0_cp*b**4*dm**2*dn**2 &
      &              *r**4+50688.0_cp*b**4*dm**2*dn*r**4-27648.0_cp*b**4*dm**2*r**4-43008.0_cp &
      &              *b**4*dn**5*r**4-190848.0_cp*b**4*dn**4*r**4+1047936.0_cp*b**4*dn** &
      &              3*r**4+2354304.0_cp*b**4*dn**2*r**4-10202496.0_cp*b**4*dn*r**4+7034112.0_cp &
      &              *b**4*r**4-2048.0_cp*b**2*dm**2*dn**3*r**6+12288.0_cp*b**2*dm**2*dn**2 &
      &              *r**6-22528.0_cp*b**2*dm**2*dn*r**6+12288.0_cp*b**2*dm**2*r**6+10240.0_cp &
      &              *b**2*dn**5*r**6+44800.0_cp*b**2*dn**4*r**6-248576.0_cp*b**2*dn**3* &
      &              r**6-550144.0_cp*b**2*dn**2*r**6+2401024.0_cp*b**2*dn*r**6-1657344.0_cp &
      &              *b**2*r**6+256.0_cp*dm**2*dn**3*r**8-1536.0_cp*dm**2*dn**2*r**8+2816.0_cp &
      &              *dm**2*dn*r**8-1536.0_cp*dm**2*r**8-512.0_cp*dn**5*r**8-2176.0_cp*dn**4 &
      &              *r**8+12416.0_cp*dn**3*r**8+25984.0_cp*dn**2*r**8-116352.0_cp*dn*r**8+ &
      &              80640.0_cp*r**8)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)* &
      &              (dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-3)=a**4*(19.0_cp*a**10*dm**2*dn**3-264.0_cp*a**10*dm**2*dn**2 &
      &              -241.0_cp*a**10*dm**2*dn+4686.0_cp*a**10*dm**2+11.0_cp*a**10*dn**5-1254.0_cp &
      &              *a**10*dn**4-2959.0_cp*a**10*dn**3+61314.0_cp*a**10*dn**2+33440.0_cp &
      &              *a**10*dn-704352.0_cp*a**10+1530.0_cp*a**8*b**2*dm**2*dn**3-14580.0_cp &
      &              *a**8*b**2*dm**2*dn**2-26370.0_cp*a**8*b**2*dm**2*dn+244620.0_cp*a**8* &
      &              b**2*dm**2+1980.0_cp*a**8*b**2*dn**5-87120.0_cp*a**8*b**2*dn**4-302940.0_cp &
      &              *a**8*b**2*dn**3+3738240.0_cp*a**8*b**2*dn**2+3952080.0_cp*a**8*b** &
      &              2*dn-39061440.0_cp*a**8*b**2-136.0_cp*a**8*dm**2*dn**3*r**2+1296.0_cp &
      &              *a**8*dm**2*dn**2*r**2+2344.0_cp*a**8*dm**2*dn*r**2-21744.0_cp*a**8 &
      &              *dm**2*r**2-120.0_cp*a**8*dn**5*r**2+5262.0_cp*a**8*dn**4*r**2+18296.0_cp &
      &              *a**8*dn**3*r**2-225726.0_cp*a**8*dn**2*r**2-238664.0_cp*a**8*dn*r**2+ &
      &              2358792.0_cp*a**8*r**2+13440.0_cp*a**6*b**4*dm**2*dn**3-80640.0_cp* &
      &              a**6*b**4*dm**2*dn**2-255360.0_cp*a**6*b**4*dm**2*dn+1330560.0_cp*a**6 &
      &              *b**4*dm**2+31680.0_cp*a**6*b**4*dn**5-654720.0_cp*a**6*b**4*dn**4-3284160.0_cp &
      &              *a**6*b**4*dn**3+24351360.0_cp*a**6*b**4*dn**2+41775360.0_cp*a**6*b**4 &
      &              *dn-232657920.0_cp*a**6*b**4-7168.0_cp*a**6*b**2*dm**2*dn**3*r**2+43008.0_cp &
      &              *a**6*b**2*dm**2*dn**2*r**2+136192.0_cp*a**6*b**2*dm**2*dn*r**2-709632.0_cp &
      &              *a**6*b**2*dm**2*r**2-11520.0_cp*a**6*b**2*dn**5*r**2+237312.0_cp*a**6 &
      &              *b**2*dn**4*r**2+1190144.0_cp*a**6*b**2*dn**3*r**2-8822784.0_cp*a** &
      &              6*b**2*dn**2*r**2-15136256.0_cp*a**6*b**2*dn*r**2+84294144.0_cp*a** &
      &              6*b**2*r**2+384.0_cp*a**6*dm**2*dn**3*r**4-2304.0_cp*a**6*dm**2*dn**2* &
      &              r**4-7296.0_cp*a**6*dm**2*dn*r**4+38016.0_cp*a**6*dm**2*r**4+384.0_cp &
      &              *a**6*dn**5*r**4-7872.0_cp*a**6*dn**4*r**4-39424.0_cp*a**6*dn**3*r**4+ &
      &              292224.0_cp*a**6*dn**2*r**4+500992.0_cp*a**6*dn*r**4-2790144.0_cp*a**6 &
      &              *r**4+26880.0_cp*a**4*b**6*dm**2*dn**3-80640.0_cp*a**4*b**6*dm**2*dn**2 &
      &              -510720.0_cp*a**4*b**6*dm**2*dn+1532160.0_cp*a**4*b**6*dm**2+118272.0_cp &
      &              *a**4*b**6*dn**5-1064448.0_cp*a**4*b**6*dn**4-8397312.0_cp*a**4*b** &
      &              6*dn**3+34417152.0_cp*a**4*b**6*dn**2+95563776.0_cp*a**4*b**6*dn-312238080.0_cp &
      &              *a**4*b**6-35840.0_cp*a**4*b**4*dm**2*dn**3*r**2+107520.0_cp*a**4*b**4 &
      &              *dm**2*dn**2*r**2+680960.0_cp*a**4*b**4*dm**2*dn*r**2-2042880.0_cp* &
      &              a**4*b**4*dm**2*r**2-107520.0_cp*a**4*b**4*dn**5*r**2+964992.0_cp*a**4 &
      &              *b**4*dn**4*r**2+7608832.0_cp*a**4*b**4*dn**3*r**2-31178112.0_cp*a**4* &
      &              b**4*dn**2*r**2-86560768.0_cp*a**4*b**4*dn*r**2+282809856.0_cp*a**4 &
      &              *b**4*r**2+11520.0_cp*a**4*b**2*dm**2*dn**3*r**4-34560.0_cp*a**4*b**2* &
      &              dm**2*dn**2*r**4-218880.0_cp*a**4*b**2*dm**2*dn*r**4+656640.0_cp*a**4* &
      &              b**2*dm**2*r**4+21504.0_cp*a**4*b**2*dn**5*r**4-192192.0_cp*a**4*b**2* &
      &              dn**4*r**4-1512960.0_cp*a**4*b**2*dn**3*r**4+6198720.0_cp*a**4*b**2 &
      &              *dn**2*r**4+17193216.0_cp*a**4*b**2*dn*r**4-56176128.0_cp*a**4*b**2 &
      &              *r**4-512.0_cp*a**4*dm**2*dn**3*r**6+1536.0_cp*a**4*dm**2*dn**2*r** &
      &              6+9728.0_cp*a**4*dm**2*dn*r**6-29184.0_cp*a**4*dm**2*r**6-512.0_cp* &
      &              a**4*dn**5*r**6+4544.0_cp*a**4*dn**4*r**6+35584.0_cp*a**4*dn**3*r** &
      &              6-145856.0_cp*a**4*dn**2*r**6-402944.0_cp*a**4*dn*r**6+1317120.0_cp &
      &              *a**4*r**6+11520.0_cp*a**2*b**8*dm**2*dn**3-218880.0_cp*a**2*b**8*dm**2 &
      &              *dn+345600.0_cp*a**2*b**8*dm**2+126720.0_cp*a**2*b**8*dn**5-253440.0_cp &
      &              *a**2*b**8*dn**4-5702400.0_cp*a**2*b**8*dn**3+8870400.0_cp*a**2*b** &
      &              8*dn**2+53729280.0_cp*a**2*b**8*dn-97320960.0_cp*a**2*b**8-28672.0_cp &
      &              *a**2*b**6*dm**2*dn**3*r**2+544768.0_cp*a**2*b**6*dm**2*dn*r**2-860160.0_cp &
      &              *a**2*b**6*dm**2*r**2-215040.0_cp*a**2*b**6*dn**5*r**2+430080.0_cp* &
      &              a**2*b**6*dn**4*r**2+9648128.0_cp*a**2*b**6*dn**3*r**2-15009792.0_cp &
      &              *a**2*b**6*dn**2*r**2-90847232.0_cp*a**2*b**6*dn*r**2+164548608.0_cp &
      &              *a**2*b**6*r**2+23040.0_cp*a**2*b**4*dm**2*dn**3*r**4-437760.0_cp*a**2 &
      &              *b**4*dm**2*dn*r**4+691200.0_cp*a**2*b**4*dm**2*r**4+107520.0_cp*a**2* &
      &              b**4*dn**5*r**4-215040.0_cp*a**2*b**4*dn**4*r**4-4800000.0_cp*a**2* &
      &              b**4*dn**3*r**4+7472640.0_cp*a**2*b**4*dn**2*r**4+45127680.0_cp*a** &
      &              2*b**4*dn*r**4-81745920.0_cp*a**2*b**4*r**4-6144.0_cp*a**2*b**2*dm**2* &
      &              dn**3*r**6+116736.0_cp*a**2*b**2*dm**2*dn*r**6-184320.0_cp*a**2*b** &
      &              2*dm**2*r**6-15360.0_cp*a**2*b**2*dn**5*r**6+30720.0_cp*a**2*b**2*dn**4 &
      &              *r**6+678912.0_cp*a**2*b**2*dn**3*r**6-1059840.0_cp*a**2*b**2*dn**2 &
      &              *r**6-6355968.0_cp*a**2*b**2*dn*r**6+11520000.0_cp*a**2*b**2*r**6+256.0_cp &
      &              *a**2*dm**2*dn**3*r**8-4864.0_cp*a**2*dm**2*dn*r**8+7680.0_cp*a**2*dm**2 &
      &              *r**8+256.0_cp*a**2*dn**5*r**8-512.0_cp*a**2*dn**4*r**8-11008.0_cp* &
      &              a**2*dn**3*r**8+17408.0_cp*a**2*dn**2*r**8+101376.0_cp*a**2*dn*r**8 &
      &              -184320.0_cp*a**2*r**8-512.0_cp*b**10*dm**2*dn**3+3072.0_cp*b**10*dm**2 &
      &              *dn**2-5632.0_cp*b**10*dm**2*dn+3072.0_cp*b**10*dm**2+33792.0_cp*b**10 &
      &              *dn**5+90112.0_cp*b**10*dn**4-754688.0_cp*b**10*dn**3-765952.0_cp*b**10 &
      &              *dn**2+5181440.0_cp*b**10*dn-3784704.0_cp*b**10+2048.0_cp*b**8*dm**2* &
      &              dn**3*r**2-12288.0_cp*b**8*dm**2*dn**2*r**2+22528.0_cp*b**8*dm**2*dn &
      &              *r**2-12288.0_cp*b**8*dm**2*r**2-92160.0_cp*b**8*dn**5*r**2-244224.0_cp &
      &              *b**8*dn**4*r**2+2054144.0_cp*b**8*dn**3*r**2+2075136.0_cp*b**8*dn**2* &
      &              r**2-14084096.0_cp*b**8*dn*r**2+10291200.0_cp*b**8*r**2-3072.0_cp*b**6 &
      &              *dm**2*dn**3*r**4+18432.0_cp*b**6*dm**2*dn**2*r**4-33792.0_cp*b**6*dm**2 &
      &              *dn*r**4+18432.0_cp*b**6*dm**2*r**4+86016.0_cp*b**6*dn**5*r**4+225792.0_cp &
      &              *b**6*dn**4*r**4-1911808.0_cp*b**6*dn**3*r**4-1915392.0_cp*b**6*dn**2* &
      &              r**4+13075456.0_cp*b**6*dn*r**4-9560064.0_cp*b**6*r**4+2048.0_cp*b**4* &
      &              dm**2*dn**3*r**6-12288.0_cp*b**4*dm**2*dn**2*r**6+22528.0_cp*b**4*dm**2 &
      &              *dn*r**6-12288.0_cp*b**4*dm**2*r**6-30720.0_cp*b**4*dn**5*r**6-79360.0_cp &
      &              *b**4*dn**4*r**6+679936.0_cp*b**4*dn**3*r**6+669184.0_cp*b**4*dn**2 &
      &              *r**6-4624384.0_cp*b**4*dn*r**6+3385344.0_cp*b**4*r**6-512.0_cp*b** &
      &              2*dm**2*dn**3*r**8+3072.0_cp*b**2*dm**2*dn**2*r**8-5632.0_cp*b**2*dm**2 &
      &              *dn*r**8+3072.0_cp*b**2*dm**2*r**8+3072.0_cp*b**2*dn**5*r**8+7680.0_cp &
      &              *b**2*dn**4*r**8-67584.0_cp*b**2*dn**3*r**8-62976.0_cp*b**2*dn**2*r**8 &
      &              +451584.0_cp*b**2*dn*r**8-331776.0_cp*b**2*r**8)/(8192.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-2)=a**3*b*(45.0_cp*a**10*dm**2*dn**2-2025.0_cp*a**10*dm**2* &
      &              dn+8370.0_cp*a**10*dm**2-198.0_cp*a**10*dn**4-8811.0_cp*a**10*dn**3 &
      &              +36036.0_cp*a**10*dn**2+306999.0_cp*a**10*dn-1288386.0_cp*a**10+1440.0_cp &
      &              *a**8*b**2*dm**2*dn**2-36000.0_cp*a**8*b**2*dm**2*dn+123840.0_cp*a**8* &
      &              b**2*dm**2-2640.0_cp*a**8*b**2*dn**4-191400.0_cp*a**8*b**2*dn**3+580800.0_cp &
      &              *a**8*b**2*dn**2+5814600.0_cp*a**8*b**2*dn-20457360.0_cp*a**8*b**2-384.0_cp &
      &              *a**8*dm**2*dn**2*r**2+9600.0_cp*a**8*dm**2*dn*r**2-33024.0_cp*a**8 &
      &              *dm**2*r**2+480.0_cp*a**8*dn**4*r**2+34668.0_cp*a**8*dn**3*r**2-105192.0_cp &
      &              *a**8*dn**2*r**2-1053372.0_cp*a**8*dn*r**2+3706008.0_cp*a**8*r**2+8064.0_cp &
      &              *a**6*b**4*dm**2*dn**2-120960.0_cp*a**6*b**4*dm**2*dn+330624.0_cp*a**6 &
      &              *b**4*dm**2-827904.0_cp*a**6*b**4*dn**3+1596672.0_cp*a**6*b**4*dn** &
      &              2+21495936.0_cp*a**6*b**4*dn-60939648.0_cp*a**6*b**4-7168.0_cp*a**6 &
      &              *b**2*dm**2*dn**2*r**2+107520.0_cp*a**6*b**2*dm**2*dn*r**2-293888.0_cp &
      &              *a**6*b**2*dm**2*r**2+499968.0_cp*a**6*b**2*dn**3*r**2-964096.0_cp* &
      &              a**6*b**2*dn**2*r**2-12980352.0_cp*a**6*b**2*dn*r**2+36797824.0_cp* &
      &              a**6*b**2*r**2+1152.0_cp*a**6*dm**2*dn**2*r**4-17280.0_cp*a**6*dm** &
      &              2*dn*r**4+47232.0_cp*a**6*dm**2*r**4-49728.0_cp*a**6*dn**3*r**4+96000.0_cp &
      &              *a**6*dn**2*r**4+1288992.0_cp*a**6*dn*r**4-3654240.0_cp*a**6*r**4+11520.0_cp &
      &              *a**4*b**6*dm**2*dn**2-103680.0_cp*a**4*b**6*dm**2*dn+207360.0_cp*a**4 &
      &              *b**6*dm**2+25344.0_cp*a**4*b**6*dn**4-1001088.0_cp*a**4*b**6*dn**3 &
      &              +709632.0_cp*a**4*b**6*dn**2+21681792.0_cp*a**4*b**6*dn-46455552.0_cp &
      &              *a**4*b**6-21504.0_cp*a**4*b**4*dm**2*dn**2*r**2+193536.0_cp*a**4*b**4 &
      &              *dm**2*dn*r**2-387072.0_cp*a**4*b**4*dm**2*r**2-32256.0_cp*a**4*b** &
      &              4*dn**4*r**2+1270080.0_cp*a**4*b**4*dn**3*r**2-900480.0_cp*a**4*b** &
      &              4*dn**2*r**2-27494208.0_cp*a**4*b**4*dn*r**2+58907520.0_cp*a**4*b** &
      &              4*r**2+11520.0_cp*a**4*b**2*dm**2*dn**2*r**4-103680.0_cp*a**4*b**2*dm**2 &
      &              *dn*r**4+207360.0_cp*a**4*b**2*dm**2*r**4+10752.0_cp*a**4*b**2*dn** &
      &              4*r**4-421344.0_cp*a**4*b**2*dn**3*r**4+300096.0_cp*a**4*b**2*dn**2 &
      &              *r**4+9102816.0_cp*a**4*b**2*dn*r**4-19503936.0_cp*a**4*b**2*r**4-1536.0_cp &
      &              *a**4*dm**2*dn**2*r**6+13824.0_cp*a**4*dm**2*dn*r**6-27648.0_cp*a** &
      &              4*dm**2*r**6-768.0_cp*a**4*dn**4*r**6+29856.0_cp*a**4*dn**3*r**6-21696.0_cp &
      &              *a**4*dn**2*r**6-640416.0_cp*a**4*dn*r**6+1372608.0_cp*a**4*r**6+3840.0_cp &
      &              *a**2*b**8*dm**2*dn**2-19200.0_cp*a**2*b**8*dm**2*dn+23040.0_cp*a** &
      &              2*b**8*dm**2+28160.0_cp*a**2*b**8*dn**4-323840.0_cp*a**2*b**8*dn**3 &
      &              -281600.0_cp*a**2*b**8*dn**2+5730560.0_cp*a**2*b**8*dn-8194560.0_cp &
      &              *a**2*b**8-12288.0_cp*a**2*b**6*dm**2*dn**2*r**2+61440.0_cp*a**2*b**6* &
      &              dm**2*dn*r**2-73728.0_cp*a**2*b**6*dm**2*r**2-61440.0_cp*a**2*b**6*dn**4 &
      &              *r**2+705024.0_cp*a**2*b**6*dn**3*r**2+611328.0_cp*a**2*b**6*dn**2* &
      &              r**2-12458496.0_cp*a**2*b**6*dn*r**2+17814528.0_cp*a**2*b**6*r**2+13824.0_cp &
      &              *a**2*b**4*dm**2*dn**2*r**4-69120.0_cp*a**2*b**4*dm**2*dn*r**4+82944.0_cp &
      &              *a**2*b**4*dm**2*r**4+43008.0_cp*a**2*b**4*dn**4*r**4-491904.0_cp*a**2 &
      &              *b**4*dn**3*r**4-423168.0_cp*a**2*b**4*dn**2*r**4+8666496.0_cp*a**2 &
      &              *b**4*dn*r**4-12393216.0_cp*a**2*b**4*r**4-6144.0_cp*a**2*b**2*dm** &
      &              2*dn**2*r**6+30720.0_cp*a**2*b**2*dm**2*dn*r**6-36864.0_cp*a**2*b** &
      &              2*dm**2*r**6-10240.0_cp*a**2*b**2*dn**4*r**6+116480.0_cp*a**2*b**2*dn**3 &
      &              *r**6+97792.0_cp*a**2*b**2*dn**2*r**6-2036480.0_cp*a**2*b**2*dn*r** &
      &              6+2913792.0_cp*a**2*b**2*r**6+768.0_cp*a**2*dm**2*dn**2*r**8-3840.0_cp &
      &              *a**2*dm**2*dn*r**8+4608.0_cp*a**2*dm**2*r**8+512.0_cp*a**2*dn**4*r**8 &
      &              -5760.0_cp*a**2*dn**3*r**8-4352.0_cp*a**2*dn**2*r**8+97920.0_cp*a** &
      &              2*dn*r**8-140544.0_cp*a**2*r**8+6144.0_cp*b**10*dn**4-13312.0_cp*b**10 &
      &              *dn**3-73728.0_cp*b**10*dn**2+222208.0_cp*b**10*dn-141312.0_cp*b**10 &
      &              -20480.0_cp*b**8*dn**4*r**2+44544.0_cp*b**8*dn**3*r**2+244736.0_cp* &
      &              b**8*dn**2*r**2-738816.0_cp*b**8*dn*r**2+470016.0_cp*b**8*r**2+24576.0_cp &
      &              *b**6*dn**4*r**4-53760.0_cp*b**6*dn**3*r**4-291840.0_cp*b**6*dn**2* &
      &              r**4+883200.0_cp*b**6*dn*r**4-562176.0_cp*b**6*r**4-12288.0_cp*b**4 &
      &              *dn**4*r**6+27136.0_cp*b**4*dn**3*r**6+144384.0_cp*b**4*dn**2*r**6-438784.0_cp &
      &              *b**4*dn*r**6+279552.0_cp*b**4*r**6+2048.0_cp*b**2*dn**4*r**8-4608.0_cp &
      &              *b**2*dn**3*r**8-23552.0_cp*b**2*dn**2*r**8+72192.0_cp*b**2*dn*r**8 &
      &              -46080.0_cp*b**2*r**8)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku-1)=a**2*(-27.0_cp*a**12*dm**2*dn**2-540.0_cp*a**12*dm**2 &
      &              *dn+3663.0_cp*a**12*dm**2-165.0_cp*a**12*dn**4-2277.0_cp*a**12*dn** &
      &              3+20559.0_cp*a**12*dn**2+79893.0_cp*a**12*dn-547866.0_cp*a**12-1440.0_cp &
      &              *a**10*b**2*dm**2*dn**2-36000.0_cp*a**10*b**2*dm**2*dn+207360.0_cp* &
      &              a**10*b**2*dm**2-12672.0_cp*a**10*b**2*dn**4-183744.0_cp*a**10*b**2 &
      &              *dn**3+1387584.0_cp*a**10*b**2*dn**2+5613696.0_cp*a**10*b**2*dn-32845824.0_cp &
      &              *a**10*b**2+128.0_cp*a**10*dm**2*dn**2*r**2+3200.0_cp*a**10*dm**2*dn &
      &              *r**2-18432.0_cp*a**10*dm**2*r**2+768.0_cp*a**10*dn**4*r**2+11088.0_cp &
      &              *a**10*dn**3*r**2-83776.0_cp*a**10*dn**2*r**2-338992.0_cp*a**10*dn* &
      &              r**2+1983456.0_cp*a**10*r**2-6720.0_cp*a**8*b**4*dm**2*dn**2-268800.0_cp &
      &              *a**8*b**4*dm**2*dn+1270080.0_cp*a**8*b**4*dm**2-110880.0_cp*a**8*b**4 &
      &              *dn**4-1737120.0_cp*a**8*b**4*dn**3+10533600.0_cp*a**8*b**4*dn**2+45202080.0_cp &
      &              *a**8*b**4*dn-218877120.0_cp*a**8*b**4+3584.0_cp*a**8*b**2*dm**2*dn**2 &
      &              *r**2+143360.0_cp*a**8*b**2*dm**2*dn*r**2-677376.0_cp*a**8*b**2*dm**2* &
      &              r**2+40320.0_cp*a**8*b**2*dn**4*r**2+628992.0_cp*a**8*b**2*dn**3*r**2- &
      &              3816064.0_cp*a**8*b**2*dn**2*r**2-16377088.0_cp*a**8*b**2*dn*r**2+79301376.0_cp &
      &              *a**8*b**2*r**2-192.0_cp*a**8*dm**2*dn**2*r**4-7680.0_cp*a**8*dm**2 &
      &              *dn*r**4+36288.0_cp*a**8*dm**2*r**4-1344.0_cp*a**8*dn**4*r**4-20832.0_cp &
      &              *a**8*dn**3*r**4+126464.0_cp*a**8*dn**2*r**4+542048.0_cp*a**8*dn*r**4- &
      &              2624832.0_cp*a**8*r**4-483840.0_cp*a**6*b**6*dm**2*dn+1774080.0_cp* &
      &              a**6*b**6*dm**2-236544.0_cp*a**6*b**6*dn**4-4257792.0_cp*a**6*b**6*dn**3 &
      &              +19396608.0_cp*a**6*b**6*dn**2+91542528.0_cp*a**6*b**6*dn-350558208.0_cp &
      &              *a**6*b**6+645120.0_cp*a**6*b**4*dm**2*dn*r**2-2365440.0_cp*a**6*b**4* &
      &              dm**2*r**2+215040.0_cp*a**6*b**4*dn**4*r**2+3854592.0_cp*a**6*b**4*dn**3 &
      &              *r**2-17568768.0_cp*a**6*b**4*dn**2*r**2-82914048.0_cp*a**6*b**4*dn &
      &              *r**2+317517312.0_cp*a**6*b**4*r**2-207360.0_cp*a**6*b**2*dm**2*dn* &
      &              r**4+760320.0_cp*a**6*b**2*dm**2*r**4-43008.0_cp*a**6*b**2*dn**4*r**4- &
      &              766080.0_cp*a**6*b**2*dn**3*r**4+3494400.0_cp*a**6*b**2*dn**2*r**4+16467840.0_cp &
      &              *a**6*b**2*dn*r**4-63067392.0_cp*a**6*b**2*r**4+9216.0_cp*a**6*dm** &
      &              2*dn*r**6-33792.0_cp*a**6*dm**2*r**6+1024.0_cp*a**6*dn**4*r**6+18048.0_cp &
      &              *a**6*dn**3*r**6-82432.0_cp*a**6*dn**2*r**6-385920.0_cp*a**6*dn*r** &
      &              6+1478400.0_cp*a**6*r**6+11520.0_cp*a**4*b**8*dm**2*dn**2-230400.0_cp &
      &              *a**4*b**8*dm**2*dn+587520.0_cp*a**4*b**8*dm**2-126720.0_cp*a**4*b**8* &
      &              dn**4-3168000.0_cp*a**4*b**8*dn**3+9504000.0_cp*a**4*b**8*dn**2+53856000.0_cp &
      &              *a**4*b**8*dn-151303680.0_cp*a**4*b**8-28672.0_cp*a**4*b**6*dm**2*dn**2 &
      &              *r**2+573440.0_cp*a**4*b**6*dm**2*dn*r**2-1462272.0_cp*a**4*b**6*dm**2 &
      &              *r**2+215040.0_cp*a**4*b**6*dn**4*r**2+5354496.0_cp*a**4*b**6*dn**3 &
      &              *r**2-16070656.0_cp*a**4*b**6*dn**2*r**2-91055104.0_cp*a**4*b**6*dn &
      &              *r**2+255811584.0_cp*a**4*b**6*r**2+23040.0_cp*a**4*b**4*dm**2*dn** &
      &              2*r**4-460800.0_cp*a**4*b**4*dm**2*dn*r**4+1175040.0_cp*a**4*b**4*dm**2 &
      &              *r**4-107520.0_cp*a**4*b**4*dn**4*r**4-2661120.0_cp*a**4*b**4*dn**3 &
      &              *r**4+7994880.0_cp*a**4*b**4*dn**2*r**4+45223680.0_cp*a**4*b**4*dn* &
      &              r**4-127065600.0_cp*a**4*b**4*r**4-6144.0_cp*a**4*b**2*dm**2*dn**2* &
      &              r**6+122880.0_cp*a**4*b**2*dm**2*dn*r**6-313344.0_cp*a**4*b**2*dm** &
      &              2*r**6+15360.0_cp*a**4*b**2*dn**4*r**6+376320.0_cp*a**4*b**2*dn**3* &
      &              r**6-1133568.0_cp*a**4*b**2*dn**2*r**6-6366720.0_cp*a**4*b**2*dn*r**6+ &
      &              17897472.0_cp*a**4*b**2*r**6+256.0_cp*a**4*dm**2*dn**2*r**8-5120.0_cp &
      &              *a**4*dm**2*dn*r**8+13056.0_cp*a**4*dm**2*r**8-256.0_cp*a**4*dn**4* &
      &              r**8-6144.0_cp*a**4*dn**3*r**8+18688.0_cp*a**4*dn**2*r**8+101376.0_cp &
      &              *a**4*dn*r**8-285696.0_cp*a**4*r**8+4096.0_cp*a**2*b**10*dm**2*dn** &
      &              2-20480.0_cp*a**2*b**10*dm**2*dn+24576.0_cp*a**2*b**10*dm**2-630784.0_cp &
      &              *a**2*b**10*dn**3+811008.0_cp*a**2*b**10*dn**2+7929856.0_cp*a**2*b**10 &
      &              *dn-14057472.0_cp*a**2*b**10-16384.0_cp*a**2*b**8*dm**2*dn**2*r**2+81920.0_cp &
      &              *a**2*b**8*dm**2*dn*r**2-98304.0_cp*a**2*b**8*dm**2*r**2+1714176.0_cp &
      &              *a**2*b**8*dn**3*r**2-2203648.0_cp*a**2*b**8*dn**2*r**2-21551104.0_cp &
      &              *a**2*b**8*dn*r**2+38203392.0_cp*a**2*b**8*r**2+24576.0_cp*a**2*b** &
      &              6*dm**2*dn**2*r**4-122880.0_cp*a**2*b**6*dm**2*dn*r**4+147456.0_cp* &
      &              a**2*b**6*dm**2*r**4-1591296.0_cp*a**2*b**6*dn**3*r**4+2048000.0_cp &
      &              *a**2*b**6*dn**2*r**4+19994624.0_cp*a**2*b**6*dn*r**4-35450880.0_cp &
      &              *a**2*b**6*r**4-16384.0_cp*a**2*b**4*dm**2*dn**2*r**6+81920.0_cp*a**2* &
      &              b**4*dm**2*dn*r**6-98304.0_cp*a**2*b**4*dm**2*r**6+563200.0_cp*a**2 &
      &              *b**4*dn**3*r**6-729088.0_cp*a**2*b**4*dn**2*r**6-7055360.0_cp*a**2 &
      &              *b**4*dn*r**6+12521472.0_cp*a**2*b**4*r**6+4096.0_cp*a**2*b**2*dm** &
      &              2*dn**2*r**8-20480.0_cp*a**2*b**2*dm**2*dn*r**8+24576.0_cp*a**2*b** &
      &              2*dm**2*r**8-55296.0_cp*a**2*b**2*dn**3*r**8+73728.0_cp*a**2*b**2*dn**2 &
      &              *r**8+681984.0_cp*a**2*b**2*dn*r**8-1216512.0_cp*a**2*b**2*r**8+4096.0_cp &
      &              *b**12*dn**4-12288.0_cp*b**12*dn**3-28672.0_cp*b**12*dn**2+110592.0_cp &
      &              *b**12*dn-73728.0_cp*b**12-16384.0_cp*b**10*dn**4*r**2+49152.0_cp*b**10 &
      &              *dn**3*r**2+114688.0_cp*b**10*dn**2*r**2-442368.0_cp*b**10*dn*r**2 &
      &              +294912.0_cp*b**10*r**2+24576.0_cp*b**8*dn**4*r**4-73728.0_cp*b**8*dn**3 &
      &              *r**4-172032.0_cp*b**8*dn**2*r**4+663552.0_cp*b**8*dn*r**4-442368.0_cp &
      &              *b**8*r**4-16384.0_cp*b**6*dn**4*r**6+49152.0_cp*b**6*dn**3*r**6+114688.0_cp &
      &              *b**6*dn**2*r**6-442368.0_cp*b**6*dn*r**6+294912.0_cp*b**6*r**6+4096.0_cp &
      &              *b**4*dn**4*r**8-12288.0_cp*b**4*dn**3*r**8-28672.0_cp*b**4*dn**2*r**8 &
      &              +110592.0_cp*b**4*dn*r**8-73728.0_cp*b**4*r**8)/(16384.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+3.0_cp))
      stencil(ku)=-a**3*b*(45.0_cp*a**10*dm**2*dn**2+225.0_cp*a**10*dm**2 &
      &              *dn-2880.0_cp*a**10*dm**2+198.0_cp*a**10*dn**4+1089.0_cp*a**10*dn** &
      &              3-20592.0_cp*a**10*dn**2-34551.0_cp*a**10*dn+441540.0_cp*a**10+840.0_cp &
      &              *a**8*b**2*dm**2*dn**2+4200.0_cp*a**8*b**2*dm**2*dn-45360.0_cp*a**8 &
      &              *b**2*dm**2+4620.0_cp*a**8*b**2*dn**4+25410.0_cp*a**8*b**2*dn**3-406560.0_cp &
      &              *a**8*b**2*dn**2-690690.0_cp*a**8*b**2*dn+7442820.0_cp*a**8*b**2-224.0_cp &
      &              *a**8*dm**2*dn**2*r**2-1120.0_cp*a**8*dm**2*dn*r**2+12096.0_cp*a**8 &
      &              *dm**2*r**2-840.0_cp*a**8*dn**4*r**2-4599.0_cp*a**8*dn**3*r**2+73654.0_cp &
      &              *a**8*dn**2*r**2+125111.0_cp*a**8*dn*r**2-1348326.0_cp*a**8*r**2+3024.0_cp &
      &              *a**6*b**4*dm**2*dn**2+15120.0_cp*a**6*b**4*dm**2*dn-133056.0_cp*a**6* &
      &              b**4*dm**2+22176.0_cp*a**6*b**4*dn**4+121968.0_cp*a**6*b**4*dn**3-1596672.0_cp &
      &              *a**6*b**4*dn**2-2760912.0_cp*a**6*b**4*dn+24216192.0_cp*a**6*b**4-2688.0_cp &
      &              *a**6*b**2*dm**2*dn**2*r**2-13440.0_cp*a**6*b**2*dm**2*dn*r**2+118272.0_cp &
      &              *a**6*b**2*dm**2*r**2-13440.0_cp*a**6*b**2*dn**4*r**2-73584.0_cp*a**6* &
      &              b**2*dn**3*r**2+964320.0_cp*a**6*b**2*dn**2*r**2+1666896.0_cp*a**6* &
      &              b**2*dn*r**2-14622720.0_cp*a**6*b**2*r**2+432.0_cp*a**6*dm**2*dn**2 &
      &              *r**4+2160.0_cp*a**6*dm**2*dn*r**4-19008.0_cp*a**6*dm**2*r**4+1344.0_cp &
      &              *a**6*dn**4*r**4+7308.0_cp*a**6*dn**3*r**4-95880.0_cp*a**6*dn**2*r**4- &
      &              165492.0_cp*a**6*dn*r**4+1452096.0_cp*a**6*r**4+2880.0_cp*a**4*b**6 &
      &              *dm**2*dn**2+14400.0_cp*a**4*b**6*dm**2*dn-97920.0_cp*a**4*b**6*dm**2+ &
      &              31680.0_cp*a**4*b**6*dn**4+174240.0_cp*a**4*b**6*dn**3-1774080.0_cp &
      &              *a**4*b**6*dn**2-3152160.0_cp*a**4*b**6*dn+21320640.0_cp*a**4*b**6-5376.0_cp &
      &              *a**4*b**4*dm**2*dn**2*r**2-26880.0_cp*a**4*b**4*dm**2*dn*r**2+182784.0_cp &
      &              *a**4*b**4*dm**2*r**2-40320.0_cp*a**4*b**4*dn**4*r**2-220752.0_cp*a**4 &
      &              *b**4*dn**3*r**2+2250528.0_cp*a**4*b**4*dn**2*r**2+3996048.0_cp*a** &
      &              4*b**4*dn*r**2-27035232.0_cp*a**4*b**4*r**2+2880.0_cp*a**4*b**2*dm**2* &
      &              dn**2*r**4+14400.0_cp*a**4*b**2*dm**2*dn*r**4-97920.0_cp*a**4*b**2*dm**2 &
      &              *r**4+13440.0_cp*a**4*b**2*dn**4*r**4+73080.0_cp*a**4*b**2*dn**3*r**4- &
      &              746160.0_cp*a**4*b**2*dn**2*r**4-1322520.0_cp*a**4*b**2*dn*r**4+8950800.0_cp &
      &              *a**4*b**2*r**4-384.0_cp*a**4*dm**2*dn**2*r**6-1920.0_cp*a**4*dm**2 &
      &              *dn*r**6+13056.0_cp*a**4*dm**2*r**6-960.0_cp*a**4*dn**4*r**6-5160.0_cp &
      &              *a**4*dn**3*r**6+52752.0_cp*a**4*dn**2*r**6+93000.0_cp*a**4*dn*r**6 &
      &              -629808.0_cp*a**4*r**6+640.0_cp*a**2*b**8*dm**2*dn**2+3200.0_cp*a** &
      &              2*b**8*dm**2*dn-15360.0_cp*a**2*b**8*dm**2+14080.0_cp*a**2*b**8*dn**4+ &
      &              77440.0_cp*a**2*b**8*dn**3-563200.0_cp*a**2*b**8*dn**2-1048960.0_cp &
      &              *a**2*b**8*dn+4984320.0_cp*a**2*b**8-2048.0_cp*a**2*b**6*dm**2*dn** &
      &              2*r**2-10240.0_cp*a**2*b**6*dm**2*dn*r**2+49152.0_cp*a**2*b**6*dm** &
      &              2*r**2-30720.0_cp*a**2*b**6*dn**4*r**2-168192.0_cp*a**2*b**6*dn**3* &
      &              r**2+1225216.0_cp*a**2*b**6*dn**2*r**2+2279168.0_cp*a**2*b**6*dn*r**2- &
      &              10834944.0_cp*a**2*b**6*r**2+2304.0_cp*a**2*b**4*dm**2*dn**2*r**4+11520.0_cp &
      &              *a**2*b**4*dm**2*dn*r**4-55296.0_cp*a**2*b**4*dm**2*r**4+21504.0_cp &
      &              *a**2*b**4*dn**4*r**4+116928.0_cp*a**2*b**4*dn**3*r**4-853632.0_cp* &
      &              a**2*b**4*dn**2*r**4-1584192.0_cp*a**2*b**4*dn*r**4+7536384.0_cp*a**2* &
      &              b**4*r**4-1024.0_cp*a**2*b**2*dm**2*dn**2*r**6-5120.0_cp*a**2*b**2*dm**2 &
      &              *dn*r**6+24576.0_cp*a**2*b**2*dm**2*r**6-5120.0_cp*a**2*b**2*dn**4* &
      &              r**6-27520.0_cp*a**2*b**2*dn**3*r**6+201472.0_cp*a**2*b**2*dn**2*r**6+ &
      &              371840.0_cp*a**2*b**2*dn*r**6-1771008.0_cp*a**2*b**2*r**6+128.0_cp* &
      &              a**2*dm**2*dn**2*r**8+640.0_cp*a**2*dm**2*dn*r**8-3072.0_cp*a**2*dm**2 &
      &              *r**8+256.0_cp*a**2*dn**4*r**8+1344.0_cp*a**2*dn**3*r**8-9856.0_cp* &
      &              a**2*dn**2*r**8-17856.0_cp*a**2*dn*r**8+85248.0_cp*a**2*r**8+1536.0_cp &
      &              *b**10*dn**4+8448.0_cp*b**10*dn**3-36864.0_cp*b**10*dn**2-76032.0_cp &
      &              *b**10*dn+207360.0_cp*b**10-5120.0_cp*b**8*dn**4*r**2-28032.0_cp*b**8* &
      &              dn**3*r**2+122624.0_cp*b**8*dn**2*r**2+252288.0_cp*b**8*dn*r**2-688896.0_cp &
      &              *b**8*r**2+6144.0_cp*b**6*dn**4*r**4+33408.0_cp*b**6*dn**3*r**4-146688.0_cp &
      &              *b**6*dn**2*r**4-300672.0_cp*b**6*dn*r**4+822528.0_cp*b**6*r**4-3072.0_cp &
      &              *b**4*dn**4*r**6-16512.0_cp*b**4*dn**3*r**6+72960.0_cp*b**4*dn**2*r**6 &
      &              +148608.0_cp*b**4*dn*r**6-407808.0_cp*b**4*r**6+512.0_cp*b**2*dn**4 &
      &              *r**8+2688.0_cp*b**2*dn**3*r**8-12032.0_cp*b**2*dn**2*r**8-24192.0_cp &
      &              *b**2*dn*r**8+66816.0_cp*b**2*r**8)/(1024.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp &
      &              )*(dn-1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+1)=-a**2*(9.0_cp*a**12*dm**2*dn**2-531.0_cp*a**12*dm**2+33.0_cp &
      &              *a**12*dn**4-3597.0_cp*a**12*dn**2+79200.0_cp*a**12+630.0_cp*a**10* &
      &              b**2*dm**2*dn**2-30870.0_cp*a**10*b**2*dm**2+2772.0_cp*a**10*b**2*dn**4 &
      &              -257796.0_cp*a**10*b**2*dn**2+4867632.0_cp*a**10*b**2-56.0_cp*a**10 &
      &              *dm**2*dn**2*r**2+2744.0_cp*a**10*dm**2*r**2-168.0_cp*a**10*dn**4*r**2 &
      &              +15568.0_cp*a**10*dn**2*r**2-293944.0_cp*a**10*r**2+5040.0_cp*a**8* &
      &              b**4*dm**2*dn**2-196560.0_cp*a**8*b**4*dm**2+27720.0_cp*a**8*b**4*dn**4 &
      &              -2134440.0_cp*a**8*b**4*dn**2+33596640.0_cp*a**8*b**4-2688.0_cp*a** &
      &              8*b**2*dm**2*dn**2*r**2+104832.0_cp*a**8*b**2*dm**2*r**2-10080.0_cp &
      &              *a**8*b**2*dn**4*r**2+773472.0_cp*a**8*b**2*dn**2*r**2-12172608.0_cp &
      &              *a**8*b**2*r**2+144.0_cp*a**8*dm**2*dn**2*r**4-5616.0_cp*a**8*dm**2 &
      &              *r**4+336.0_cp*a**8*dn**4*r**4-25632.0_cp*a**8*dn**2*r**4+402912.0_cp &
      &              *a**8*r**4+10080.0_cp*a**6*b**6*dm**2*dn**2-292320.0_cp*a**6*b**6*dm**2 &
      &              +73920.0_cp*a**6*b**6*dn**4-4509120.0_cp*a**6*b**6*dn**2+56770560.0_cp &
      &              *a**6*b**6-13440.0_cp*a**6*b**4*dm**2*dn**2*r**2+389760.0_cp*a**6*b**4 &
      &              *dm**2*r**2-67200.0_cp*a**6*b**4*dn**4*r**2+4085760.0_cp*a**6*b**4*dn**2 &
      &              *r**2-51421440.0_cp*a**6*b**4*r**2+4320.0_cp*a**6*b**2*dm**2*dn**2* &
      &              r**4-125280.0_cp*a**6*b**2*dm**2*r**4+13440.0_cp*a**6*b**2*dn**4*r**4- &
      &              812640.0_cp*a**6*b**2*dn**2*r**4+10213920.0_cp*a**6*b**2*r**4-192.0_cp &
      &              *a**6*dm**2*dn**2*r**6+5568.0_cp*a**6*dm**2*r**6-320.0_cp*a**6*dn** &
      &              4*r**6+19136.0_cp*a**6*dn**2*r**6-239424.0_cp*a**6*r**6+5760.0_cp*a**4 &
      &              *b**8*dm**2*dn**2-109440.0_cp*a**4*b**8*dm**2+63360.0_cp*a**4*b**8*dn**4 &
      &              -2851200.0_cp*a**4*b**8*dn**2+26864640.0_cp*a**4*b**8-14336.0_cp*a**4* &
      &              b**6*dm**2*dn**2*r**2+272384.0_cp*a**4*b**6*dm**2*r**2-107520.0_cp* &
      &              a**4*b**6*dn**4*r**2+4824064.0_cp*a**4*b**6*dn**2*r**2-45423616.0_cp &
      &              *a**4*b**6*r**2+11520.0_cp*a**4*b**4*dm**2*dn**2*r**4-218880.0_cp*a**4 &
      &              *b**4*dm**2*r**4+53760.0_cp*a**4*b**4*dn**4*r**4-2400000.0_cp*a**4* &
      &              b**4*dn**2*r**4+22563840.0_cp*a**4*b**4*r**4-3072.0_cp*a**4*b**2*dm**2 &
      &              *dn**2*r**6+58368.0_cp*a**4*b**2*dm**2*r**6-7680.0_cp*a**4*b**2*dn**4* &
      &              r**6+339456.0_cp*a**4*b**2*dn**2*r**6-3177984.0_cp*a**4*b**2*r**6+128.0_cp &
      &              *a**4*dm**2*dn**2*r**8-2432.0_cp*a**4*dm**2*r**8+128.0_cp*a**4*dn** &
      &              4*r**8-5504.0_cp*a**4*dn**2*r**8+50688.0_cp*a**4*r**8+768.0_cp*a**2 &
      &              *b**10*dm**2*dn**2-6912.0_cp*a**2*b**10*dm**2+16896.0_cp*a**2*b**10 &
      &              *dn**4-489984.0_cp*a**2*b**10*dn**2+3041280.0_cp*a**2*b**10-3072.0_cp &
      &              *a**2*b**8*dm**2*dn**2*r**2+27648.0_cp*a**2*b**8*dm**2*r**2-46080.0_cp &
      &              *a**2*b**8*dn**4*r**2+1333248.0_cp*a**2*b**8*dn**2*r**2-8266752.0_cp &
      &              *a**2*b**8*r**2+4608.0_cp*a**2*b**6*dm**2*dn**2*r**4-41472.0_cp*a** &
      &              2*b**6*dm**2*r**4+43008.0_cp*a**2*b**6*dn**4*r**4-1239552.0_cp*a**2 &
      &              *b**6*dn**2*r**4+7672320.0_cp*a**2*b**6*r**4-3072.0_cp*a**2*b**4*dm**2 &
      &              *dn**2*r**6+27648.0_cp*a**2*b**4*dm**2*r**6-15360.0_cp*a**2*b**4*dn**4 &
      &              *r**6+439296.0_cp*a**2*b**4*dn**2*r**6-2709504.0_cp*a**2*b**4*r**6+768.0_cp &
      &              *a**2*b**2*dm**2*dn**2*r**8-6912.0_cp*a**2*b**2*dm**2*r**8+1536.0_cp &
      &              *a**2*b**2*dn**4*r**8-43008.0_cp*a**2*b**2*dn**2*r**8+262656.0_cp*a**2 &
      &              *b**2*r**8+1024.0_cp*b**12*dn**4-13312.0_cp*b**12*dn**2+36864.0_cp* &
      &              b**12-4096.0_cp*b**10*dn**4*r**2+53248.0_cp*b**10*dn**2*r**2-147456.0_cp &
      &              *b**10*r**2+6144.0_cp*b**8*dn**4*r**4-79872.0_cp*b**8*dn**2*r**4+221184.0_cp &
      &              *b**8*r**4-4096.0_cp*b**6*dn**4*r**6+53248.0_cp*b**6*dn**2*r**6-147456.0_cp &
      &              *b**6*r**6+1024.0_cp*b**4*dn**4*r**8-13312.0_cp*b**4*dn**2*r**8+36864.0_cp &
      &              *b**4*r**8)/(2048.0_cp*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+2)=-a**3*b*(45.0_cp*a**10*dm**2*dn**2-225.0_cp*a**10*dm**2* &
      &              dn-2880.0_cp*a**10*dm**2+198.0_cp*a**10*dn**4-1089.0_cp*a**10*dn**3 &
      &              -20592.0_cp*a**10*dn**2+34551.0_cp*a**10*dn+441540.0_cp*a**10+840.0_cp &
      &              *a**8*b**2*dm**2*dn**2-4200.0_cp*a**8*b**2*dm**2*dn-45360.0_cp*a**8 &
      &              *b**2*dm**2+4620.0_cp*a**8*b**2*dn**4-25410.0_cp*a**8*b**2*dn**3-406560.0_cp &
      &              *a**8*b**2*dn**2+690690.0_cp*a**8*b**2*dn+7442820.0_cp*a**8*b**2-224.0_cp &
      &              *a**8*dm**2*dn**2*r**2+1120.0_cp*a**8*dm**2*dn*r**2+12096.0_cp*a**8 &
      &              *dm**2*r**2-840.0_cp*a**8*dn**4*r**2+4599.0_cp*a**8*dn**3*r**2+73654.0_cp &
      &              *a**8*dn**2*r**2-125111.0_cp*a**8*dn*r**2-1348326.0_cp*a**8*r**2+3024.0_cp &
      &              *a**6*b**4*dm**2*dn**2-15120.0_cp*a**6*b**4*dm**2*dn-133056.0_cp*a**6* &
      &              b**4*dm**2+22176.0_cp*a**6*b**4*dn**4-121968.0_cp*a**6*b**4*dn**3-1596672.0_cp &
      &              *a**6*b**4*dn**2+2760912.0_cp*a**6*b**4*dn+24216192.0_cp*a**6*b**4-2688.0_cp &
      &              *a**6*b**2*dm**2*dn**2*r**2+13440.0_cp*a**6*b**2*dm**2*dn*r**2+118272.0_cp &
      &              *a**6*b**2*dm**2*r**2-13440.0_cp*a**6*b**2*dn**4*r**2+73584.0_cp*a**6* &
      &              b**2*dn**3*r**2+964320.0_cp*a**6*b**2*dn**2*r**2-1666896.0_cp*a**6* &
      &              b**2*dn*r**2-14622720.0_cp*a**6*b**2*r**2+432.0_cp*a**6*dm**2*dn**2 &
      &              *r**4-2160.0_cp*a**6*dm**2*dn*r**4-19008.0_cp*a**6*dm**2*r**4+1344.0_cp &
      &              *a**6*dn**4*r**4-7308.0_cp*a**6*dn**3*r**4-95880.0_cp*a**6*dn**2*r**4+ &
      &              165492.0_cp*a**6*dn*r**4+1452096.0_cp*a**6*r**4+2880.0_cp*a**4*b**6 &
      &              *dm**2*dn**2-14400.0_cp*a**4*b**6*dm**2*dn-97920.0_cp*a**4*b**6*dm**2+ &
      &              31680.0_cp*a**4*b**6*dn**4-174240.0_cp*a**4*b**6*dn**3-1774080.0_cp &
      &              *a**4*b**6*dn**2+3152160.0_cp*a**4*b**6*dn+21320640.0_cp*a**4*b**6-5376.0_cp &
      &              *a**4*b**4*dm**2*dn**2*r**2+26880.0_cp*a**4*b**4*dm**2*dn*r**2+182784.0_cp &
      &              *a**4*b**4*dm**2*r**2-40320.0_cp*a**4*b**4*dn**4*r**2+220752.0_cp*a**4 &
      &              *b**4*dn**3*r**2+2250528.0_cp*a**4*b**4*dn**2*r**2-3996048.0_cp*a** &
      &              4*b**4*dn*r**2-27035232.0_cp*a**4*b**4*r**2+2880.0_cp*a**4*b**2*dm**2* &
      &              dn**2*r**4-14400.0_cp*a**4*b**2*dm**2*dn*r**4-97920.0_cp*a**4*b**2*dm**2 &
      &              *r**4+13440.0_cp*a**4*b**2*dn**4*r**4-73080.0_cp*a**4*b**2*dn**3*r**4- &
      &              746160.0_cp*a**4*b**2*dn**2*r**4+1322520.0_cp*a**4*b**2*dn*r**4+8950800.0_cp &
      &              *a**4*b**2*r**4-384.0_cp*a**4*dm**2*dn**2*r**6+1920.0_cp*a**4*dm**2 &
      &              *dn*r**6+13056.0_cp*a**4*dm**2*r**6-960.0_cp*a**4*dn**4*r**6+5160.0_cp &
      &              *a**4*dn**3*r**6+52752.0_cp*a**4*dn**2*r**6-93000.0_cp*a**4*dn*r**6 &
      &              -629808.0_cp*a**4*r**6+640.0_cp*a**2*b**8*dm**2*dn**2-3200.0_cp*a** &
      &              2*b**8*dm**2*dn-15360.0_cp*a**2*b**8*dm**2+14080.0_cp*a**2*b**8*dn**4- &
      &              77440.0_cp*a**2*b**8*dn**3-563200.0_cp*a**2*b**8*dn**2+1048960.0_cp &
      &              *a**2*b**8*dn+4984320.0_cp*a**2*b**8-2048.0_cp*a**2*b**6*dm**2*dn** &
      &              2*r**2+10240.0_cp*a**2*b**6*dm**2*dn*r**2+49152.0_cp*a**2*b**6*dm** &
      &              2*r**2-30720.0_cp*a**2*b**6*dn**4*r**2+168192.0_cp*a**2*b**6*dn**3* &
      &              r**2+1225216.0_cp*a**2*b**6*dn**2*r**2-2279168.0_cp*a**2*b**6*dn*r**2- &
      &              10834944.0_cp*a**2*b**6*r**2+2304.0_cp*a**2*b**4*dm**2*dn**2*r**4-11520.0_cp &
      &              *a**2*b**4*dm**2*dn*r**4-55296.0_cp*a**2*b**4*dm**2*r**4+21504.0_cp &
      &              *a**2*b**4*dn**4*r**4-116928.0_cp*a**2*b**4*dn**3*r**4-853632.0_cp* &
      &              a**2*b**4*dn**2*r**4+1584192.0_cp*a**2*b**4*dn*r**4+7536384.0_cp*a**2* &
      &              b**4*r**4-1024.0_cp*a**2*b**2*dm**2*dn**2*r**6+5120.0_cp*a**2*b**2*dm**2 &
      &              *dn*r**6+24576.0_cp*a**2*b**2*dm**2*r**6-5120.0_cp*a**2*b**2*dn**4* &
      &              r**6+27520.0_cp*a**2*b**2*dn**3*r**6+201472.0_cp*a**2*b**2*dn**2*r**6- &
      &              371840.0_cp*a**2*b**2*dn*r**6-1771008.0_cp*a**2*b**2*r**6+128.0_cp* &
      &              a**2*dm**2*dn**2*r**8-640.0_cp*a**2*dm**2*dn*r**8-3072.0_cp*a**2*dm**2 &
      &              *r**8+256.0_cp*a**2*dn**4*r**8-1344.0_cp*a**2*dn**3*r**8-9856.0_cp* &
      &              a**2*dn**2*r**8+17856.0_cp*a**2*dn*r**8+85248.0_cp*a**2*r**8+1536.0_cp &
      &              *b**10*dn**4-8448.0_cp*b**10*dn**3-36864.0_cp*b**10*dn**2+76032.0_cp &
      &              *b**10*dn+207360.0_cp*b**10-5120.0_cp*b**8*dn**4*r**2+28032.0_cp*b**8* &
      &              dn**3*r**2+122624.0_cp*b**8*dn**2*r**2-252288.0_cp*b**8*dn*r**2-688896.0_cp &
      &              *b**8*r**2+6144.0_cp*b**6*dn**4*r**4-33408.0_cp*b**6*dn**3*r**4-146688.0_cp &
      &              *b**6*dn**2*r**4+300672.0_cp*b**6*dn*r**4+822528.0_cp*b**6*r**4-3072.0_cp &
      &              *b**4*dn**4*r**6+16512.0_cp*b**4*dn**3*r**6+72960.0_cp*b**4*dn**2*r**6 &
      &              -148608.0_cp*b**4*dn*r**6-407808.0_cp*b**4*r**6+512.0_cp*b**2*dn**4 &
      &              *r**8-2688.0_cp*b**2*dn**3*r**8-12032.0_cp*b**2*dn**2*r**8+24192.0_cp &
      &              *b**2*dn*r**8+66816.0_cp*b**2*r**8)/(1024.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+3)=a**2*(-27.0_cp*a**12*dm**2*dn**2+540.0_cp*a**12*dm**2 &
      &              *dn+3663.0_cp*a**12*dm**2-165.0_cp*a**12*dn**4+2277.0_cp*a**12*dn** &
      &              3+20559.0_cp*a**12*dn**2-79893.0_cp*a**12*dn-547866.0_cp*a**12-1440.0_cp &
      &              *a**10*b**2*dm**2*dn**2+36000.0_cp*a**10*b**2*dm**2*dn+207360.0_cp* &
      &              a**10*b**2*dm**2-12672.0_cp*a**10*b**2*dn**4+183744.0_cp*a**10*b**2 &
      &              *dn**3+1387584.0_cp*a**10*b**2*dn**2-5613696.0_cp*a**10*b**2*dn-32845824.0_cp &
      &              *a**10*b**2+128.0_cp*a**10*dm**2*dn**2*r**2-3200.0_cp*a**10*dm**2*dn &
      &              *r**2-18432.0_cp*a**10*dm**2*r**2+768.0_cp*a**10*dn**4*r**2-11088.0_cp &
      &              *a**10*dn**3*r**2-83776.0_cp*a**10*dn**2*r**2+338992.0_cp*a**10*dn* &
      &              r**2+1983456.0_cp*a**10*r**2-6720.0_cp*a**8*b**4*dm**2*dn**2+268800.0_cp &
      &              *a**8*b**4*dm**2*dn+1270080.0_cp*a**8*b**4*dm**2-110880.0_cp*a**8*b**4 &
      &              *dn**4+1737120.0_cp*a**8*b**4*dn**3+10533600.0_cp*a**8*b**4*dn**2-45202080.0_cp &
      &              *a**8*b**4*dn-218877120.0_cp*a**8*b**4+3584.0_cp*a**8*b**2*dm**2*dn**2 &
      &              *r**2-143360.0_cp*a**8*b**2*dm**2*dn*r**2-677376.0_cp*a**8*b**2*dm**2* &
      &              r**2+40320.0_cp*a**8*b**2*dn**4*r**2-628992.0_cp*a**8*b**2*dn**3*r**2- &
      &              3816064.0_cp*a**8*b**2*dn**2*r**2+16377088.0_cp*a**8*b**2*dn*r**2+79301376.0_cp &
      &              *a**8*b**2*r**2-192.0_cp*a**8*dm**2*dn**2*r**4+7680.0_cp*a**8*dm**2 &
      &              *dn*r**4+36288.0_cp*a**8*dm**2*r**4-1344.0_cp*a**8*dn**4*r**4+20832.0_cp &
      &              *a**8*dn**3*r**4+126464.0_cp*a**8*dn**2*r**4-542048.0_cp*a**8*dn*r**4- &
      &              2624832.0_cp*a**8*r**4+483840.0_cp*a**6*b**6*dm**2*dn+1774080.0_cp* &
      &              a**6*b**6*dm**2-236544.0_cp*a**6*b**6*dn**4+4257792.0_cp*a**6*b**6*dn**3 &
      &              +19396608.0_cp*a**6*b**6*dn**2-91542528.0_cp*a**6*b**6*dn-350558208.0_cp &
      &              *a**6*b**6-645120.0_cp*a**6*b**4*dm**2*dn*r**2-2365440.0_cp*a**6*b**4* &
      &              dm**2*r**2+215040.0_cp*a**6*b**4*dn**4*r**2-3854592.0_cp*a**6*b**4*dn**3 &
      &              *r**2-17568768.0_cp*a**6*b**4*dn**2*r**2+82914048.0_cp*a**6*b**4*dn &
      &              *r**2+317517312.0_cp*a**6*b**4*r**2+207360.0_cp*a**6*b**2*dm**2*dn* &
      &              r**4+760320.0_cp*a**6*b**2*dm**2*r**4-43008.0_cp*a**6*b**2*dn**4*r**4+ &
      &              766080.0_cp*a**6*b**2*dn**3*r**4+3494400.0_cp*a**6*b**2*dn**2*r**4-16467840.0_cp &
      &              *a**6*b**2*dn*r**4-63067392.0_cp*a**6*b**2*r**4-9216.0_cp*a**6*dm** &
      &              2*dn*r**6-33792.0_cp*a**6*dm**2*r**6+1024.0_cp*a**6*dn**4*r**6-18048.0_cp &
      &              *a**6*dn**3*r**6-82432.0_cp*a**6*dn**2*r**6+385920.0_cp*a**6*dn*r** &
      &              6+1478400.0_cp*a**6*r**6+11520.0_cp*a**4*b**8*dm**2*dn**2+230400.0_cp &
      &              *a**4*b**8*dm**2*dn+587520.0_cp*a**4*b**8*dm**2-126720.0_cp*a**4*b**8* &
      &              dn**4+3168000.0_cp*a**4*b**8*dn**3+9504000.0_cp*a**4*b**8*dn**2-53856000.0_cp &
      &              *a**4*b**8*dn-151303680.0_cp*a**4*b**8-28672.0_cp*a**4*b**6*dm**2*dn**2 &
      &              *r**2-573440.0_cp*a**4*b**6*dm**2*dn*r**2-1462272.0_cp*a**4*b**6*dm**2 &
      &              *r**2+215040.0_cp*a**4*b**6*dn**4*r**2-5354496.0_cp*a**4*b**6*dn**3 &
      &              *r**2-16070656.0_cp*a**4*b**6*dn**2*r**2+91055104.0_cp*a**4*b**6*dn &
      &              *r**2+255811584.0_cp*a**4*b**6*r**2+23040.0_cp*a**4*b**4*dm**2*dn** &
      &              2*r**4+460800.0_cp*a**4*b**4*dm**2*dn*r**4+1175040.0_cp*a**4*b**4*dm**2 &
      &              *r**4-107520.0_cp*a**4*b**4*dn**4*r**4+2661120.0_cp*a**4*b**4*dn**3 &
      &              *r**4+7994880.0_cp*a**4*b**4*dn**2*r**4-45223680.0_cp*a**4*b**4*dn* &
      &              r**4-127065600.0_cp*a**4*b**4*r**4-6144.0_cp*a**4*b**2*dm**2*dn**2* &
      &              r**6-122880.0_cp*a**4*b**2*dm**2*dn*r**6-313344.0_cp*a**4*b**2*dm** &
      &              2*r**6+15360.0_cp*a**4*b**2*dn**4*r**6-376320.0_cp*a**4*b**2*dn**3* &
      &              r**6-1133568.0_cp*a**4*b**2*dn**2*r**6+6366720.0_cp*a**4*b**2*dn*r**6+ &
      &              17897472.0_cp*a**4*b**2*r**6+256.0_cp*a**4*dm**2*dn**2*r**8+5120.0_cp &
      &              *a**4*dm**2*dn*r**8+13056.0_cp*a**4*dm**2*r**8-256.0_cp*a**4*dn**4* &
      &              r**8+6144.0_cp*a**4*dn**3*r**8+18688.0_cp*a**4*dn**2*r**8-101376.0_cp &
      &              *a**4*dn*r**8-285696.0_cp*a**4*r**8+4096.0_cp*a**2*b**10*dm**2*dn** &
      &              2+20480.0_cp*a**2*b**10*dm**2*dn+24576.0_cp*a**2*b**10*dm**2+630784.0_cp &
      &              *a**2*b**10*dn**3+811008.0_cp*a**2*b**10*dn**2-7929856.0_cp*a**2*b**10 &
      &              *dn-14057472.0_cp*a**2*b**10-16384.0_cp*a**2*b**8*dm**2*dn**2*r**2-81920.0_cp &
      &              *a**2*b**8*dm**2*dn*r**2-98304.0_cp*a**2*b**8*dm**2*r**2-1714176.0_cp &
      &              *a**2*b**8*dn**3*r**2-2203648.0_cp*a**2*b**8*dn**2*r**2+21551104.0_cp &
      &              *a**2*b**8*dn*r**2+38203392.0_cp*a**2*b**8*r**2+24576.0_cp*a**2*b** &
      &              6*dm**2*dn**2*r**4+122880.0_cp*a**2*b**6*dm**2*dn*r**4+147456.0_cp* &
      &              a**2*b**6*dm**2*r**4+1591296.0_cp*a**2*b**6*dn**3*r**4+2048000.0_cp &
      &              *a**2*b**6*dn**2*r**4-19994624.0_cp*a**2*b**6*dn*r**4-35450880.0_cp &
      &              *a**2*b**6*r**4-16384.0_cp*a**2*b**4*dm**2*dn**2*r**6-81920.0_cp*a**2* &
      &              b**4*dm**2*dn*r**6-98304.0_cp*a**2*b**4*dm**2*r**6-563200.0_cp*a**2 &
      &              *b**4*dn**3*r**6-729088.0_cp*a**2*b**4*dn**2*r**6+7055360.0_cp*a**2 &
      &              *b**4*dn*r**6+12521472.0_cp*a**2*b**4*r**6+4096.0_cp*a**2*b**2*dm** &
      &              2*dn**2*r**8+20480.0_cp*a**2*b**2*dm**2*dn*r**8+24576.0_cp*a**2*b** &
      &              2*dm**2*r**8+55296.0_cp*a**2*b**2*dn**3*r**8+73728.0_cp*a**2*b**2*dn**2 &
      &              *r**8-681984.0_cp*a**2*b**2*dn*r**8-1216512.0_cp*a**2*b**2*r**8+4096.0_cp &
      &              *b**12*dn**4+12288.0_cp*b**12*dn**3-28672.0_cp*b**12*dn**2-110592.0_cp &
      &              *b**12*dn-73728.0_cp*b**12-16384.0_cp*b**10*dn**4*r**2-49152.0_cp*b**10 &
      &              *dn**3*r**2+114688.0_cp*b**10*dn**2*r**2+442368.0_cp*b**10*dn*r**2 &
      &              +294912.0_cp*b**10*r**2+24576.0_cp*b**8*dn**4*r**4+73728.0_cp*b**8*dn**3 &
      &              *r**4-172032.0_cp*b**8*dn**2*r**4-663552.0_cp*b**8*dn*r**4-442368.0_cp &
      &              *b**8*r**4-16384.0_cp*b**6*dn**4*r**6-49152.0_cp*b**6*dn**3*r**6+114688.0_cp &
      &              *b**6*dn**2*r**6+442368.0_cp*b**6*dn*r**6+294912.0_cp*b**6*r**6+4096.0_cp &
      &              *b**4*dn**4*r**8+12288.0_cp*b**4*dn**3*r**8-28672.0_cp*b**4*dn**2*r**8 &
      &              -110592.0_cp*b**4*dn*r**8-73728.0_cp*b**4*r**8)/(16384.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+4)=a**3*b*(45.0_cp*a**10*dm**2*dn**2+2025.0_cp*a**10*dm**2* &
      &              dn+8370.0_cp*a**10*dm**2-198.0_cp*a**10*dn**4+8811.0_cp*a**10*dn**3 &
      &              +36036.0_cp*a**10*dn**2-306999.0_cp*a**10*dn-1288386.0_cp*a**10+1440.0_cp &
      &              *a**8*b**2*dm**2*dn**2+36000.0_cp*a**8*b**2*dm**2*dn+123840.0_cp*a**8* &
      &              b**2*dm**2-2640.0_cp*a**8*b**2*dn**4+191400.0_cp*a**8*b**2*dn**3+580800.0_cp &
      &              *a**8*b**2*dn**2-5814600.0_cp*a**8*b**2*dn-20457360.0_cp*a**8*b**2-384.0_cp &
      &              *a**8*dm**2*dn**2*r**2-9600.0_cp*a**8*dm**2*dn*r**2-33024.0_cp*a**8 &
      &              *dm**2*r**2+480.0_cp*a**8*dn**4*r**2-34668.0_cp*a**8*dn**3*r**2-105192.0_cp &
      &              *a**8*dn**2*r**2+1053372.0_cp*a**8*dn*r**2+3706008.0_cp*a**8*r**2+8064.0_cp &
      &              *a**6*b**4*dm**2*dn**2+120960.0_cp*a**6*b**4*dm**2*dn+330624.0_cp*a**6 &
      &              *b**4*dm**2+827904.0_cp*a**6*b**4*dn**3+1596672.0_cp*a**6*b**4*dn** &
      &              2-21495936.0_cp*a**6*b**4*dn-60939648.0_cp*a**6*b**4-7168.0_cp*a**6 &
      &              *b**2*dm**2*dn**2*r**2-107520.0_cp*a**6*b**2*dm**2*dn*r**2-293888.0_cp &
      &              *a**6*b**2*dm**2*r**2-499968.0_cp*a**6*b**2*dn**3*r**2-964096.0_cp* &
      &              a**6*b**2*dn**2*r**2+12980352.0_cp*a**6*b**2*dn*r**2+36797824.0_cp* &
      &              a**6*b**2*r**2+1152.0_cp*a**6*dm**2*dn**2*r**4+17280.0_cp*a**6*dm** &
      &              2*dn*r**4+47232.0_cp*a**6*dm**2*r**4+49728.0_cp*a**6*dn**3*r**4+96000.0_cp &
      &              *a**6*dn**2*r**4-1288992.0_cp*a**6*dn*r**4-3654240.0_cp*a**6*r**4+11520.0_cp &
      &              *a**4*b**6*dm**2*dn**2+103680.0_cp*a**4*b**6*dm**2*dn+207360.0_cp*a**4 &
      &              *b**6*dm**2+25344.0_cp*a**4*b**6*dn**4+1001088.0_cp*a**4*b**6*dn**3 &
      &              +709632.0_cp*a**4*b**6*dn**2-21681792.0_cp*a**4*b**6*dn-46455552.0_cp &
      &              *a**4*b**6-21504.0_cp*a**4*b**4*dm**2*dn**2*r**2-193536.0_cp*a**4*b**4 &
      &              *dm**2*dn*r**2-387072.0_cp*a**4*b**4*dm**2*r**2-32256.0_cp*a**4*b** &
      &              4*dn**4*r**2-1270080.0_cp*a**4*b**4*dn**3*r**2-900480.0_cp*a**4*b** &
      &              4*dn**2*r**2+27494208.0_cp*a**4*b**4*dn*r**2+58907520.0_cp*a**4*b** &
      &              4*r**2+11520.0_cp*a**4*b**2*dm**2*dn**2*r**4+103680.0_cp*a**4*b**2*dm**2 &
      &              *dn*r**4+207360.0_cp*a**4*b**2*dm**2*r**4+10752.0_cp*a**4*b**2*dn** &
      &              4*r**4+421344.0_cp*a**4*b**2*dn**3*r**4+300096.0_cp*a**4*b**2*dn**2 &
      &              *r**4-9102816.0_cp*a**4*b**2*dn*r**4-19503936.0_cp*a**4*b**2*r**4-1536.0_cp &
      &              *a**4*dm**2*dn**2*r**6-13824.0_cp*a**4*dm**2*dn*r**6-27648.0_cp*a** &
      &              4*dm**2*r**6-768.0_cp*a**4*dn**4*r**6-29856.0_cp*a**4*dn**3*r**6-21696.0_cp &
      &              *a**4*dn**2*r**6+640416.0_cp*a**4*dn*r**6+1372608.0_cp*a**4*r**6+3840.0_cp &
      &              *a**2*b**8*dm**2*dn**2+19200.0_cp*a**2*b**8*dm**2*dn+23040.0_cp*a** &
      &              2*b**8*dm**2+28160.0_cp*a**2*b**8*dn**4+323840.0_cp*a**2*b**8*dn**3 &
      &              -281600.0_cp*a**2*b**8*dn**2-5730560.0_cp*a**2*b**8*dn-8194560.0_cp &
      &              *a**2*b**8-12288.0_cp*a**2*b**6*dm**2*dn**2*r**2-61440.0_cp*a**2*b**6* &
      &              dm**2*dn*r**2-73728.0_cp*a**2*b**6*dm**2*r**2-61440.0_cp*a**2*b**6*dn**4 &
      &              *r**2-705024.0_cp*a**2*b**6*dn**3*r**2+611328.0_cp*a**2*b**6*dn**2* &
      &              r**2+12458496.0_cp*a**2*b**6*dn*r**2+17814528.0_cp*a**2*b**6*r**2+13824.0_cp &
      &              *a**2*b**4*dm**2*dn**2*r**4+69120.0_cp*a**2*b**4*dm**2*dn*r**4+82944.0_cp &
      &              *a**2*b**4*dm**2*r**4+43008.0_cp*a**2*b**4*dn**4*r**4+491904.0_cp*a**2 &
      &              *b**4*dn**3*r**4-423168.0_cp*a**2*b**4*dn**2*r**4-8666496.0_cp*a**2 &
      &              *b**4*dn*r**4-12393216.0_cp*a**2*b**4*r**4-6144.0_cp*a**2*b**2*dm** &
      &              2*dn**2*r**6-30720.0_cp*a**2*b**2*dm**2*dn*r**6-36864.0_cp*a**2*b** &
      &              2*dm**2*r**6-10240.0_cp*a**2*b**2*dn**4*r**6-116480.0_cp*a**2*b**2*dn**3 &
      &              *r**6+97792.0_cp*a**2*b**2*dn**2*r**6+2036480.0_cp*a**2*b**2*dn*r** &
      &              6+2913792.0_cp*a**2*b**2*r**6+768.0_cp*a**2*dm**2*dn**2*r**8+3840.0_cp &
      &              *a**2*dm**2*dn*r**8+4608.0_cp*a**2*dm**2*r**8+512.0_cp*a**2*dn**4*r**8 &
      &              +5760.0_cp*a**2*dn**3*r**8-4352.0_cp*a**2*dn**2*r**8-97920.0_cp*a** &
      &              2*dn*r**8-140544.0_cp*a**2*r**8+6144.0_cp*b**10*dn**4+13312.0_cp*b**10 &
      &              *dn**3-73728.0_cp*b**10*dn**2-222208.0_cp*b**10*dn-141312.0_cp*b**10 &
      &              -20480.0_cp*b**8*dn**4*r**2-44544.0_cp*b**8*dn**3*r**2+244736.0_cp* &
      &              b**8*dn**2*r**2+738816.0_cp*b**8*dn*r**2+470016.0_cp*b**8*r**2+24576.0_cp &
      &              *b**6*dn**4*r**4+53760.0_cp*b**6*dn**3*r**4-291840.0_cp*b**6*dn**2* &
      &              r**4-883200.0_cp*b**6*dn*r**4-562176.0_cp*b**6*r**4-12288.0_cp*b**4 &
      &              *dn**4*r**6-27136.0_cp*b**4*dn**3*r**6+144384.0_cp*b**4*dn**2*r**6+438784.0_cp &
      &              *b**4*dn*r**6+279552.0_cp*b**4*r**6+2048.0_cp*b**2*dn**4*r**8+4608.0_cp &
      &              *b**2*dn**3*r**8-23552.0_cp*b**2*dn**2*r**8-72192.0_cp*b**2*dn*r**8 &
      &              -46080.0_cp*b**2*r**8)/(4096.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+5)=a**4*(19.0_cp*a**10*dm**2*dn**3+264.0_cp*a**10*dm**2*dn**2 &
      &              -241.0_cp*a**10*dm**2*dn-4686.0_cp*a**10*dm**2+11.0_cp*a**10*dn**5+1254.0_cp &
      &              *a**10*dn**4-2959.0_cp*a**10*dn**3-61314.0_cp*a**10*dn**2+33440.0_cp &
      &              *a**10*dn+704352.0_cp*a**10+1530.0_cp*a**8*b**2*dm**2*dn**3+14580.0_cp &
      &              *a**8*b**2*dm**2*dn**2-26370.0_cp*a**8*b**2*dm**2*dn-244620.0_cp*a**8* &
      &              b**2*dm**2+1980.0_cp*a**8*b**2*dn**5+87120.0_cp*a**8*b**2*dn**4-302940.0_cp &
      &              *a**8*b**2*dn**3-3738240.0_cp*a**8*b**2*dn**2+3952080.0_cp*a**8*b** &
      &              2*dn+39061440.0_cp*a**8*b**2-136.0_cp*a**8*dm**2*dn**3*r**2-1296.0_cp &
      &              *a**8*dm**2*dn**2*r**2+2344.0_cp*a**8*dm**2*dn*r**2+21744.0_cp*a**8 &
      &              *dm**2*r**2-120.0_cp*a**8*dn**5*r**2-5262.0_cp*a**8*dn**4*r**2+18296.0_cp &
      &              *a**8*dn**3*r**2+225726.0_cp*a**8*dn**2*r**2-238664.0_cp*a**8*dn*r**2- &
      &              2358792.0_cp*a**8*r**2+13440.0_cp*a**6*b**4*dm**2*dn**3+80640.0_cp* &
      &              a**6*b**4*dm**2*dn**2-255360.0_cp*a**6*b**4*dm**2*dn-1330560.0_cp*a**6 &
      &              *b**4*dm**2+31680.0_cp*a**6*b**4*dn**5+654720.0_cp*a**6*b**4*dn**4-3284160.0_cp &
      &              *a**6*b**4*dn**3-24351360.0_cp*a**6*b**4*dn**2+41775360.0_cp*a**6*b**4 &
      &              *dn+232657920.0_cp*a**6*b**4-7168.0_cp*a**6*b**2*dm**2*dn**3*r**2-43008.0_cp &
      &              *a**6*b**2*dm**2*dn**2*r**2+136192.0_cp*a**6*b**2*dm**2*dn*r**2+709632.0_cp &
      &              *a**6*b**2*dm**2*r**2-11520.0_cp*a**6*b**2*dn**5*r**2-237312.0_cp*a**6 &
      &              *b**2*dn**4*r**2+1190144.0_cp*a**6*b**2*dn**3*r**2+8822784.0_cp*a** &
      &              6*b**2*dn**2*r**2-15136256.0_cp*a**6*b**2*dn*r**2-84294144.0_cp*a** &
      &              6*b**2*r**2+384.0_cp*a**6*dm**2*dn**3*r**4+2304.0_cp*a**6*dm**2*dn**2* &
      &              r**4-7296.0_cp*a**6*dm**2*dn*r**4-38016.0_cp*a**6*dm**2*r**4+384.0_cp &
      &              *a**6*dn**5*r**4+7872.0_cp*a**6*dn**4*r**4-39424.0_cp*a**6*dn**3*r**4- &
      &              292224.0_cp*a**6*dn**2*r**4+500992.0_cp*a**6*dn*r**4+2790144.0_cp*a**6 &
      &              *r**4+26880.0_cp*a**4*b**6*dm**2*dn**3+80640.0_cp*a**4*b**6*dm**2*dn**2 &
      &              -510720.0_cp*a**4*b**6*dm**2*dn-1532160.0_cp*a**4*b**6*dm**2+118272.0_cp &
      &              *a**4*b**6*dn**5+1064448.0_cp*a**4*b**6*dn**4-8397312.0_cp*a**4*b** &
      &              6*dn**3-34417152.0_cp*a**4*b**6*dn**2+95563776.0_cp*a**4*b**6*dn+312238080.0_cp &
      &              *a**4*b**6-35840.0_cp*a**4*b**4*dm**2*dn**3*r**2-107520.0_cp*a**4*b**4 &
      &              *dm**2*dn**2*r**2+680960.0_cp*a**4*b**4*dm**2*dn*r**2+2042880.0_cp* &
      &              a**4*b**4*dm**2*r**2-107520.0_cp*a**4*b**4*dn**5*r**2-964992.0_cp*a**4 &
      &              *b**4*dn**4*r**2+7608832.0_cp*a**4*b**4*dn**3*r**2+31178112.0_cp*a**4* &
      &              b**4*dn**2*r**2-86560768.0_cp*a**4*b**4*dn*r**2-282809856.0_cp*a**4 &
      &              *b**4*r**2+11520.0_cp*a**4*b**2*dm**2*dn**3*r**4+34560.0_cp*a**4*b**2* &
      &              dm**2*dn**2*r**4-218880.0_cp*a**4*b**2*dm**2*dn*r**4-656640.0_cp*a**4* &
      &              b**2*dm**2*r**4+21504.0_cp*a**4*b**2*dn**5*r**4+192192.0_cp*a**4*b**2* &
      &              dn**4*r**4-1512960.0_cp*a**4*b**2*dn**3*r**4-6198720.0_cp*a**4*b**2 &
      &              *dn**2*r**4+17193216.0_cp*a**4*b**2*dn*r**4+56176128.0_cp*a**4*b**2 &
      &              *r**4-512.0_cp*a**4*dm**2*dn**3*r**6-1536.0_cp*a**4*dm**2*dn**2*r** &
      &              6+9728.0_cp*a**4*dm**2*dn*r**6+29184.0_cp*a**4*dm**2*r**6-512.0_cp* &
      &              a**4*dn**5*r**6-4544.0_cp*a**4*dn**4*r**6+35584.0_cp*a**4*dn**3*r** &
      &              6+145856.0_cp*a**4*dn**2*r**6-402944.0_cp*a**4*dn*r**6-1317120.0_cp &
      &              *a**4*r**6+11520.0_cp*a**2*b**8*dm**2*dn**3-218880.0_cp*a**2*b**8*dm**2 &
      &              *dn-345600.0_cp*a**2*b**8*dm**2+126720.0_cp*a**2*b**8*dn**5+253440.0_cp &
      &              *a**2*b**8*dn**4-5702400.0_cp*a**2*b**8*dn**3-8870400.0_cp*a**2*b** &
      &              8*dn**2+53729280.0_cp*a**2*b**8*dn+97320960.0_cp*a**2*b**8-28672.0_cp &
      &              *a**2*b**6*dm**2*dn**3*r**2+544768.0_cp*a**2*b**6*dm**2*dn*r**2+860160.0_cp &
      &              *a**2*b**6*dm**2*r**2-215040.0_cp*a**2*b**6*dn**5*r**2-430080.0_cp* &
      &              a**2*b**6*dn**4*r**2+9648128.0_cp*a**2*b**6*dn**3*r**2+15009792.0_cp &
      &              *a**2*b**6*dn**2*r**2-90847232.0_cp*a**2*b**6*dn*r**2-164548608.0_cp &
      &              *a**2*b**6*r**2+23040.0_cp*a**2*b**4*dm**2*dn**3*r**4-437760.0_cp*a**2 &
      &              *b**4*dm**2*dn*r**4-691200.0_cp*a**2*b**4*dm**2*r**4+107520.0_cp*a**2* &
      &              b**4*dn**5*r**4+215040.0_cp*a**2*b**4*dn**4*r**4-4800000.0_cp*a**2* &
      &              b**4*dn**3*r**4-7472640.0_cp*a**2*b**4*dn**2*r**4+45127680.0_cp*a** &
      &              2*b**4*dn*r**4+81745920.0_cp*a**2*b**4*r**4-6144.0_cp*a**2*b**2*dm**2* &
      &              dn**3*r**6+116736.0_cp*a**2*b**2*dm**2*dn*r**6+184320.0_cp*a**2*b** &
      &              2*dm**2*r**6-15360.0_cp*a**2*b**2*dn**5*r**6-30720.0_cp*a**2*b**2*dn**4 &
      &              *r**6+678912.0_cp*a**2*b**2*dn**3*r**6+1059840.0_cp*a**2*b**2*dn**2 &
      &              *r**6-6355968.0_cp*a**2*b**2*dn*r**6-11520000.0_cp*a**2*b**2*r**6+256.0_cp &
      &              *a**2*dm**2*dn**3*r**8-4864.0_cp*a**2*dm**2*dn*r**8-7680.0_cp*a**2*dm**2 &
      &              *r**8+256.0_cp*a**2*dn**5*r**8+512.0_cp*a**2*dn**4*r**8-11008.0_cp* &
      &              a**2*dn**3*r**8-17408.0_cp*a**2*dn**2*r**8+101376.0_cp*a**2*dn*r**8 &
      &              +184320.0_cp*a**2*r**8-512.0_cp*b**10*dm**2*dn**3-3072.0_cp*b**10*dm**2 &
      &              *dn**2-5632.0_cp*b**10*dm**2*dn-3072.0_cp*b**10*dm**2+33792.0_cp*b**10 &
      &              *dn**5-90112.0_cp*b**10*dn**4-754688.0_cp*b**10*dn**3+765952.0_cp*b**10 &
      &              *dn**2+5181440.0_cp*b**10*dn+3784704.0_cp*b**10+2048.0_cp*b**8*dm**2* &
      &              dn**3*r**2+12288.0_cp*b**8*dm**2*dn**2*r**2+22528.0_cp*b**8*dm**2*dn &
      &              *r**2+12288.0_cp*b**8*dm**2*r**2-92160.0_cp*b**8*dn**5*r**2+244224.0_cp &
      &              *b**8*dn**4*r**2+2054144.0_cp*b**8*dn**3*r**2-2075136.0_cp*b**8*dn**2* &
      &              r**2-14084096.0_cp*b**8*dn*r**2-10291200.0_cp*b**8*r**2-3072.0_cp*b**6 &
      &              *dm**2*dn**3*r**4-18432.0_cp*b**6*dm**2*dn**2*r**4-33792.0_cp*b**6*dm**2 &
      &              *dn*r**4-18432.0_cp*b**6*dm**2*r**4+86016.0_cp*b**6*dn**5*r**4-225792.0_cp &
      &              *b**6*dn**4*r**4-1911808.0_cp*b**6*dn**3*r**4+1915392.0_cp*b**6*dn**2* &
      &              r**4+13075456.0_cp*b**6*dn*r**4+9560064.0_cp*b**6*r**4+2048.0_cp*b**4* &
      &              dm**2*dn**3*r**6+12288.0_cp*b**4*dm**2*dn**2*r**6+22528.0_cp*b**4*dm**2 &
      &              *dn*r**6+12288.0_cp*b**4*dm**2*r**6-30720.0_cp*b**4*dn**5*r**6+79360.0_cp &
      &              *b**4*dn**4*r**6+679936.0_cp*b**4*dn**3*r**6-669184.0_cp*b**4*dn**2 &
      &              *r**6-4624384.0_cp*b**4*dn*r**6-3385344.0_cp*b**4*r**6-512.0_cp*b** &
      &              2*dm**2*dn**3*r**8-3072.0_cp*b**2*dm**2*dn**2*r**8-5632.0_cp*b**2*dm**2 &
      &              *dn*r**8-3072.0_cp*b**2*dm**2*r**8+3072.0_cp*b**2*dn**5*r**8-7680.0_cp &
      &              *b**2*dn**4*r**8-67584.0_cp*b**2*dn**3*r**8+62976.0_cp*b**2*dn**2*r**8 &
      &              +451584.0_cp*b**2*dn*r**8+331776.0_cp*b**2*r**8)/(8192.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+6)=-a**5*b*(-145.0_cp*a**8*dm**2*dn**3-750.0_cp*a**8*dm**2* &
      &              dn**2+4705.0_cp*a**8*dm**2*dn+21750.0_cp*a**8*dm**2-330.0_cp*a**8*dn**5 &
      &              -4675.0_cp*a**8*dn**4+40205.0_cp*a**8*dn**3+234025.0_cp*a**8*dn**2-709115.0_cp &
      &              *a**8*dn-3377550.0_cp*a**8-2640.0_cp*a**6*b**2*dm**2*dn**3-7200.0_cp &
      &              *a**6*b**2*dm**2*dn**2+86160.0_cp*a**6*b**2*dm**2*dn+280800.0_cp*a**6* &
      &              b**2*dm**2-9240.0_cp*a**6*b**2*dn**5-69300.0_cp*a**6*b**2*dn**4+867900.0_cp &
      &              *a**6*b**2*dn**3+3217500.0_cp*a**6*b**2*dn**2-13926660.0_cp*a**6*b**2* &
      &              dn-47104200.0_cp*a**6*b**2+704.0_cp*a**6*dm**2*dn**3*r**2+1920.0_cp &
      &              *a**6*dm**2*dn**2*r**2-22976.0_cp*a**6*dm**2*dn*r**2-74880.0_cp*a** &
      &              6*dm**2*r**2+1680.0_cp*a**6*dn**5*r**2+12570.0_cp*a**6*dn**4*r**2-157246.0_cp &
      &              *a**6*dn**3*r**2-582870.0_cp*a**6*dn**2*r**2+2522974.0_cp*a**6*dn*r**2 &
      &              +8533260.0_cp*a**6*r**2-8064.0_cp*a**4*b**4*dm**2*dn**3+274176.0_cp &
      &              *a**4*b**4*dm**2*dn+604800.0_cp*a**4*b**4*dm**2-50688.0_cp*a**4*b** &
      &              4*dn**5-126720.0_cp*a**4*b**4*dn**4+3497472.0_cp*a**4*b**4*dn**3+6526080.0_cp &
      &              *a**4*b**4*dn**2-49547520.0_cp*a**4*b**4*dn-114998400.0_cp*a**4*b** &
      &              4+7168.0_cp*a**4*b**2*dm**2*dn**3*r**2-243712.0_cp*a**4*b**2*dm**2*dn &
      &              *r**2-537600.0_cp*a**4*b**2*dm**2*r**2+30720.0_cp*a**4*b**2*dn**5*r**2 &
      &              +76800.0_cp*a**4*b**2*dn**4*r**2-2112512.0_cp*a**4*b**2*dn**3*r**2-3941760.0_cp &
      &              *a**4*b**2*dn**2*r**2+29919488.0_cp*a**4*b**2*dn*r**2+69440640.0_cp &
      &              *a**4*b**2*r**2-1152.0_cp*a**4*dm**2*dn**3*r**4+39168.0_cp*a**4*dm**2* &
      &              dn*r**4+86400.0_cp*a**4*dm**2*r**4-3072.0_cp*a**4*dn**5*r**4-7680.0_cp &
      &              *a**4*dn**4*r**4+210048.0_cp*a**4*dn**3*r**4+392160.0_cp*a**4*dn**2 &
      &              *r**4-2971200.0_cp*a**4*dn*r**4-6896160.0_cp*a**4*r**4-3840.0_cp*a**2* &
      &              b**6*dm**2*dn**3+23040.0_cp*a**2*b**6*dm**2*dn**2+188160.0_cp*a**2* &
      &              b**6*dm**2*dn+253440.0_cp*a**2*b**6*dm**2-76032.0_cp*a**2*b**6*dn** &
      &              5+105600.0_cp*a**2*b**6*dn**4+3501696.0_cp*a**2*b**6*dn**3+274560.0_cp &
      &              *a**2*b**6*dn**2-42252672.0_cp*a**2*b**6*dn-61712640.0_cp*a**2*b**6 &
      &              +7168.0_cp*a**2*b**4*dm**2*dn**3*r**2-43008.0_cp*a**2*b**4*dm**2*dn**2 &
      &              *r**2-351232.0_cp*a**2*b**4*dm**2*dn*r**2-473088.0_cp*a**2*b**4*dm**2* &
      &              r**2+96768.0_cp*a**2*b**4*dn**5*r**2-133056.0_cp*a**2*b**4*dn**4*r**2- &
      &              4442816.0_cp*a**2*b**4*dn**3*r**2-353472.0_cp*a**2*b**4*dn**2*r**2+53579456.0_cp &
      &              *a**2*b**4*dn*r**2+78255744.0_cp*a**2*b**4*r**2-3840.0_cp*a**2*b**2 &
      &              *dm**2*dn**3*r**4+23040.0_cp*a**2*b**2*dm**2*dn**2*r**4+188160.0_cp &
      &              *a**2*b**2*dm**2*dn*r**4+253440.0_cp*a**2*b**2*dm**2*r**4-32256.0_cp &
      &              *a**2*b**2*dn**5*r**4+43680.0_cp*a**2*b**2*dn**4*r**4+1473568.0_cp* &
      &              a**2*b**2*dn**3*r**4+122400.0_cp*a**2*b**2*dn**2*r**4-17740576.0_cp &
      &              *a**2*b**2*dn*r**4-25913280.0_cp*a**2*b**2*r**4+512.0_cp*a**2*dm**2 &
      &              *dn**3*r**6-3072.0_cp*a**2*dm**2*dn**2*r**6-25088.0_cp*a**2*dm**2*dn &
      &              *r**6-33792.0_cp*a**2*dm**2*r**6+2304.0_cp*a**2*dn**5*r**6-3040.0_cp &
      &              *a**2*dn**4*r**6-104288.0_cp*a**2*dn**3*r**6-9824.0_cp*a**2*dn**2*r**6 &
      &              +1248608.0_cp*a**2*dn*r**6+1824576.0_cp*a**2*r**6+1280.0_cp*b**8*dm**2 &
      &              *dn**3+7680.0_cp*b**8*dm**2*dn**2+14080.0_cp*b**8*dm**2*dn+7680.0_cp &
      &              *b**8*dm**2-28160.0_cp*b**8*dn**5+126720.0_cp*b**8*dn**4+689920.0_cp &
      &              *b**8*dn**3-1562880.0_cp*b**8*dn**2-6744320.0_cp*b**8*dn-4646400.0_cp &
      &              *b**8-4096.0_cp*b**6*dm**2*dn**3*r**2-24576.0_cp*b**6*dm**2*dn**2*r**2 &
      &              -45056.0_cp*b**6*dm**2*dn*r**2-24576.0_cp*b**6*dm**2*r**2+61440.0_cp &
      &              *b**6*dn**5*r**2-274944.0_cp*b**6*dn**4*r**2-1501696.0_cp*b**6*dn** &
      &              3*r**2+3393024.0_cp*b**6*dn**2*r**2+14662144.0_cp*b**6*dn*r**2+10103808.0_cp &
      &              *b**6*r**2+4608.0_cp*b**4*dm**2*dn**3*r**4+27648.0_cp*b**4*dm**2*dn**2 &
      &              *r**4+50688.0_cp*b**4*dm**2*dn*r**4+27648.0_cp*b**4*dm**2*r**4-43008.0_cp &
      &              *b**4*dn**5*r**4+190848.0_cp*b**4*dn**4*r**4+1047936.0_cp*b**4*dn** &
      &              3*r**4-2354304.0_cp*b**4*dn**2*r**4-10202496.0_cp*b**4*dn*r**4-7034112.0_cp &
      &              *b**4*r**4-2048.0_cp*b**2*dm**2*dn**3*r**6-12288.0_cp*b**2*dm**2*dn**2 &
      &              *r**6-22528.0_cp*b**2*dm**2*dn*r**6-12288.0_cp*b**2*dm**2*r**6+10240.0_cp &
      &              *b**2*dn**5*r**6-44800.0_cp*b**2*dn**4*r**6-248576.0_cp*b**2*dn**3* &
      &              r**6+550144.0_cp*b**2*dn**2*r**6+2401024.0_cp*b**2*dn*r**6+1657344.0_cp &
      &              *b**2*r**6+256.0_cp*dm**2*dn**3*r**8+1536.0_cp*dm**2*dn**2*r**8+2816.0_cp &
      &              *dm**2*dn*r**8+1536.0_cp*dm**2*r**8-512.0_cp*dn**5*r**8+2176.0_cp*dn**4 &
      &              *r**8+12416.0_cp*dn**3*r**8-25984.0_cp*dn**2*r**8-116352.0_cp*dn*r**8- &
      &              80640.0_cp*r**8)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)* &
      &              (dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+7)=-a**6*(-39.0_cp*a**8*dm**2*dn**3-66.0_cp*a**8*dm**2*dn**2 &
      &              +1971.0_cp*a**8*dm**2*dn+6534.0_cp*a**8*dm**2-121.0_cp*a**8*dn**5-671.0_cp &
      &              *a**8*dn**4+13849.0_cp*a**8*dn**3+43571.0_cp*a**8*dn**2-293040.0_cp &
      &              *a**8*dn-991188.0_cp*a**8-2160.0_cp*a**6*b**2*dm**2*dn**3+2160.0_cp &
      &              *a**6*b**2*dm**2*dn**2+116640.0_cp*a**6*b**2*dm**2*dn+293760.0_cp*a**6 &
      &              *b**2*dm**2-10560.0_cp*a**6*b**2*dn**5-19360.0_cp*a**6*b**2*dn**4+955680.0_cp &
      &              *a**6*b**2*dn**3+1640320.0_cp*a**6*b**2*dn**2-18437760.0_cp*a**6*b**2* &
      &              dn-47646720.0_cp*a**6*b**2+192.0_cp*a**6*dm**2*dn**3*r**2-192.0_cp* &
      &              a**6*dm**2*dn**2*r**2-10368.0_cp*a**6*dm**2*dn*r**2-26112.0_cp*a**6 &
      &              *dm**2*r**2+640.0_cp*a**6*dn**5*r**2+1176.0_cp*a**6*dn**4*r**2-57712.0_cp &
      &              *a**6*dn**3*r**2-99048.0_cp*a**6*dn**2*r**2+1113408.0_cp*a**6*dn*r**2+ &
      &              2877216.0_cp*a**6*r**2-10080.0_cp*a**4*b**4*dm**2*dn**3+60480.0_cp* &
      &              a**4*b**4*dm**2*dn**2+695520.0_cp*a**4*b**4*dm**2*dn+1270080.0_cp*a**4 &
      &              *b**4*dm**2-102960.0_cp*a**4*b**4*dn**5+134640.0_cp*a**4*b**4*dn**4 &
      &              +6977520.0_cp*a**4*b**4*dn**3+2051280.0_cp*a**4*b**4*dn**2-120985920.0_cp &
      &              *a**4*b**4*dn-228951360.0_cp*a**4*b**4+5376.0_cp*a**4*b**2*dm**2*dn**3 &
      &              *r**2-32256.0_cp*a**4*b**2*dm**2*dn**2*r**2-370944.0_cp*a**4*b**2*dm**2 &
      &              *dn*r**2-677376.0_cp*a**4*b**2*dm**2*r**2+37440.0_cp*a**4*b**2*dn** &
      &              5*r**2-48384.0_cp*a**4*b**2*dn**4*r**2-2528448.0_cp*a**4*b**2*dn**3 &
      &              *r**2-744192.0_cp*a**4*b**2*dn**2*r**2+43834752.0_cp*a**4*b**2*dn*r**2 &
      &              +82950912.0_cp*a**4*b**2*r**2-288.0_cp*a**4*dm**2*dn**3*r**4+1728.0_cp &
      &              *a**4*dm**2*dn**2*r**4+19872.0_cp*a**4*dm**2*dn*r**4+36288.0_cp*a** &
      &              4*dm**2*r**4-1248.0_cp*a**4*dn**5*r**4+1584.0_cp*a**4*dn**4*r**4+83808.0_cp &
      &              *a**4*dn**3*r**4+24912.0_cp*a**4*dn**2*r**4-1450944.0_cp*a**4*dn*r**4- &
      &              2745792.0_cp*a**4*r**4+161280.0_cp*a**2*b**6*dm**2*dn**2+806400.0_cp &
      &              *a**2*b**6*dm**2*dn+967680.0_cp*a**2*b**6*dm**2-236544.0_cp*a**2*b**6* &
      &              dn**5+946176.0_cp*a**2*b**6*dn**4+10881024.0_cp*a**2*b**6*dn**3-14429184.0_cp &
      &              *a**2*b**6*dn**2-167473152.0_cp*a**2*b**6*dn-212889600.0_cp*a**2*b**6- &
      &              215040.0_cp*a**2*b**4*dm**2*dn**2*r**2-1075200.0_cp*a**2*b**4*dm**2 &
      &              *dn*r**2-1290240.0_cp*a**2*b**4*dm**2*r**2+215040.0_cp*a**2*b**4*dn**5 &
      &              *r**2-854784.0_cp*a**2*b**4*dn**4*r**2-9859584.0_cp*a**2*b**4*dn**3 &
      &              *r**2+13058304.0_cp*a**2*b**4*dn**2*r**2+151689216.0_cp*a**2*b**4*dn &
      &              *r**2+192826368.0_cp*a**2*b**4*r**2+69120.0_cp*a**2*b**2*dm**2*dn** &
      &              2*r**4+345600.0_cp*a**2*b**2*dm**2*dn*r**4+414720.0_cp*a**2*b**2*dm**2 &
      &              *r**4-43008.0_cp*a**2*b**2*dn**5*r**4+169344.0_cp*a**2*b**2*dn**4*r**4 &
      &              +1962240.0_cp*a**2*b**2*dn**3*r**4-2586240.0_cp*a**2*b**2*dn**2*r** &
      &              4-30131712.0_cp*a**2*b**2*dn*r**4-38306304.0_cp*a**2*b**2*r**4-3072.0_cp &
      &              *a**2*dm**2*dn**2*r**6-15360.0_cp*a**2*dm**2*dn*r**6-18432.0_cp*a** &
      &              2*dm**2*r**6+1024.0_cp*a**2*dn**5*r**6-3968.0_cp*a**2*dn**4*r**6-46336.0_cp &
      &              *a**2*dn**3*r**6+60032.0_cp*a**2*dn**2*r**6+706560.0_cp*a**2*dn*r** &
      &              6+898560.0_cp*a**2*r**6+11520.0_cp*b**8*dm**2*dn**3+69120.0_cp*b**8 &
      &              *dm**2*dn**2+126720.0_cp*b**8*dm**2*dn+69120.0_cp*b**8*dm**2-126720.0_cp &
      &              *b**8*dn**5+802560.0_cp*b**8*dn**4+3168000.0_cp*b**8*dn**3-12460800.0_cp &
      &              *b**8*dn**2-43591680.0_cp*b**8*dn-28892160.0_cp*b**8-28672.0_cp*b** &
      &              6*dm**2*dn**3*r**2-172032.0_cp*b**6*dm**2*dn**2*r**2-315392.0_cp*b**6* &
      &              dm**2*dn*r**2-172032.0_cp*b**6*dm**2*r**2+215040.0_cp*b**6*dn**5*r**2- &
      &              1354752.0_cp*b**6*dn**4*r**2-5361664.0_cp*b**6*dn**3*r**2+21052416.0_cp &
      &              *b**6*dn**2*r**2+73701376.0_cp*b**6*dn*r**2+48857088.0_cp*b**6*r**2 &
      &              +23040.0_cp*b**4*dm**2*dn**3*r**4+138240.0_cp*b**4*dm**2*dn**2*r**4 &
      &              +253440.0_cp*b**4*dm**2*dn*r**4+138240.0_cp*b**4*dm**2*r**4-107520.0_cp &
      &              *b**4*dn**5*r**4+672000.0_cp*b**4*dn**4*r**4+2672640.0_cp*b**4*dn** &
      &              3*r**4-10440960.0_cp*b**4*dn**2*r**4-36618240.0_cp*b**4*dn*r**4-24284160.0_cp &
      &              *b**4*r**4-6144.0_cp*b**2*dm**2*dn**3*r**6-36864.0_cp*b**2*dm**2*dn**2 &
      &              *r**6-67584.0_cp*b**2*dm**2*dn*r**6-36864.0_cp*b**2*dm**2*r**6+15360.0_cp &
      &              *b**2*dn**5*r**6-94720.0_cp*b**2*dn**4*r**6-380928.0_cp*b**2*dn**3* &
      &              r**6+1464832.0_cp*b**2*dn**2*r**6+5164032.0_cp*b**2*dn*r**6+3428352.0_cp &
      &              *b**2*r**6+256.0_cp*dm**2*dn**3*r**8+1536.0_cp*dm**2*dn**2*r**8+2816.0_cp &
      &              *dm**2*dn*r**8+1536.0_cp*dm**2*r**8-256.0_cp*dn**5*r**8+1536.0_cp*dn**4 &
      &              *r**8+6400.0_cp*dn**3*r**8-23040.0_cp*dn**2*r**8-82944.0_cp*dn*r**8 &
      &              -55296.0_cp*r**8)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+8)=a**7*b*(25.0_cp*a**6*dm**2*dn**3-210.0_cp*a**6*dm**2*dn**2 &
      &              -2575.0_cp*a**6*dm**2*dn-5460.0_cp*a**6*dm**2+198.0_cp*a**6*dn**5-385.0_cp &
      &              *a**6*dn**4-17303.0_cp*a**6*dn**3-2387.0_cp*a**6*dn**2+398189.0_cp* &
      &              a**6*dn+861168.0_cp*a**6+120.0_cp*a**4*b**2*dm**2*dn**3-5040.0_cp*a**4 &
      &              *b**2*dm**2*dn**2-34680.0_cp*a**4*b**2*dm**2*dn-55440.0_cp*a**4*b** &
      &              2*dm**2+4180.0_cp*a**4*b**2*dn**5-17710.0_cp*a**4*b**2*dn**4-274450.0_cp &
      &              *a**4*b**2*dn**3+327250.0_cp*a**4*b**2*dn**2+5827470.0_cp*a**4*b**2 &
      &              *dn+9577260.0_cp*a**4*b**2-32.0_cp*a**4*dm**2*dn**3*r**2+1344.0_cp* &
      &              a**4*dm**2*dn**2*r**2+9248.0_cp*a**4*dm**2*dn*r**2+14784.0_cp*a**4*dm**2 &
      &              *r**2-760.0_cp*a**4*dn**5*r**2+3199.0_cp*a**4*dn**4*r**2+49721.0_cp &
      &              *a**4*dn**3*r**2-59269.0_cp*a**4*dn**2*r**2-1055689.0_cp*a**4*dn*r**2- &
      &              1734978.0_cp*a**4*r**2-1008.0_cp*a**2*b**4*dm**2*dn**3-18144.0_cp*a**2 &
      &              *b**4*dm**2*dn**2-71568.0_cp*a**2*b**4*dm**2*dn-78624.0_cp*a**2*b** &
      &              4*dm**2+15840.0_cp*a**2*b**4*dn**5-99792.0_cp*a**2*b**4*dn**4-704880.0_cp &
      &              *a**2*b**4*dn**3+2073456.0_cp*a**2*b**4*dn**2+13931280.0_cp*a**2*b**4* &
      &              dn+16033248.0_cp*a**2*b**4+896.0_cp*a**2*b**2*dm**2*dn**3*r**2+16128.0_cp &
      &              *a**2*b**2*dm**2*dn**2*r**2+63616.0_cp*a**2*b**2*dm**2*dn*r**2+69888.0_cp &
      &              *a**2*b**2*dm**2*r**2-9600.0_cp*a**2*b**2*dn**5*r**2+60144.0_cp*a** &
      &              2*b**2*dn**4*r**2+425744.0_cp*a**2*b**2*dn**3*r**2-1251600.0_cp*a** &
      &              2*b**2*dn**2*r**2-8412176.0_cp*a**2*b**2*dn*r**2-9681504.0_cp*a**2* &
      &              b**2*r**2-144.0_cp*a**2*dm**2*dn**3*r**4-2592.0_cp*a**2*dm**2*dn**2 &
      &              *r**4-10224.0_cp*a**2*dm**2*dn*r**4-11232.0_cp*a**2*dm**2*r**4+960.0_cp &
      &              *a**2*dn**5*r**4-5964.0_cp*a**2*dn**4*r**4-42372.0_cp*a**2*dn**3*r**4+ &
      &              124116.0_cp*a**2*dn**2*r**4+835428.0_cp*a**2*dn*r**4+961560.0_cp*a**2* &
      &              r**4-1920.0_cp*b**6*dm**2*dn**3-11520.0_cp*b**6*dm**2*dn**2-21120.0_cp &
      &              *b**6*dm**2*dn-11520.0_cp*b**6*dm**2+12672.0_cp*b**6*dn**5-103488.0_cp &
      &              *b**6*dn**4-302016.0_cp*b**6*dn**3+1915584.0_cp*b**6*dn**2+5915712.0_cp &
      &              *b**6*dn+3814272.0_cp*b**6+3584.0_cp*b**4*dm**2*dn**3*r**2+21504.0_cp &
      &              *b**4*dm**2*dn**2*r**2+39424.0_cp*b**4*dm**2*dn*r**2+21504.0_cp*b** &
      &              4*dm**2*r**2-16128.0_cp*b**4*dn**5*r**2+131040.0_cp*b**4*dn**4*r**2 &
      &              +383264.0_cp*b**4*dn**3*r**2-2427936.0_cp*b**4*dn**2*r**2-7501088.0_cp &
      &              *b**4*dn*r**2-4837056.0_cp*b**4*r**2-1920.0_cp*b**2*dm**2*dn**3*r** &
      &              4-11520.0_cp*b**2*dm**2*dn**2*r**4-21120.0_cp*b**2*dm**2*dn*r**4-11520.0_cp &
      &              *b**2*dm**2*r**4+5376.0_cp*b**2*dn**5*r**4-43344.0_cp*b**2*dn**4*r**4- &
      &              127408.0_cp*b**2*dn**3*r**4+802992.0_cp*b**2*dn**2*r**4+2484016.0_cp &
      &              *b**2*dn*r**4+1602336.0_cp*b**2*r**4+256.0_cp*dm**2*dn**3*r**6+1536.0_cp &
      &              *dm**2*dn**2*r**6+2816.0_cp*dm**2*dn*r**6+1536.0_cp*dm**2*r**6-384.0_cp &
      &              *dn**5*r**6+3056.0_cp*dn**4*r**6+9104.0_cp*dn**3*r**6-56336.0_cp*dn**2 &
      &              *r**6-174992.0_cp*dn*r**6-112992.0_cp*r**6)/(2048.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+9)=a**8*(a**6*dm**2*dn**3-48.0_cp*a**6*dm**2*dn**2-379.0_cp &
      &              *a**6*dm**2*dn-708.0_cp*a**6*dm**2+25.0_cp*a**6*dn**5-124.0_cp*a**6 &
      &              *dn**4-2085.0_cp*a**6*dn**3+2932.0_cp*a**6*dn**2+57440.0_cp*a**6*dn &
      &              +109056.0_cp*a**6-90.0_cp*a**4*b**2*dm**2*dn**3-3240.0_cp*a**4*b**2 &
      &              *dm**2*dn**2-17190.0_cp*a**4*b**2*dm**2*dn-24840.0_cp*a**4*b**2*dm**2+ &
      &              1716.0_cp*a**4*b**2*dn**5-11616.0_cp*a**4*b**2*dn**4-106788.0_cp*a**4* &
      &              b**2*dn**3+304128.0_cp*a**4*b**2*dn**2+2802096.0_cp*a**4*b**2*dn+4143744.0_cp &
      &              *a**4*b**2+8.0_cp*a**4*dm**2*dn**3*r**2+288.0_cp*a**4*dm**2*dn**2*r**2 &
      &              +1528.0_cp*a**4*dm**2*dn*r**2+2208.0_cp*a**4*dm**2*r**2-104.0_cp*a**4* &
      &              dn**5*r**2+700.0_cp*a**4*dn**4*r**2+6448.0_cp*a**4*dn**3*r**2-18364.0_cp &
      &              *a**4*dn**2*r**2-169208.0_cp*a**4*dn*r**2-250224.0_cp*a**4*r**2-1680.0_cp &
      &              *a**2*b**4*dm**2*dn**3-20160.0_cp*a**2*b**4*dm**2*dn**2-68880.0_cp* &
      &              a**2*b**4*dm**2*dn-70560.0_cp*a**2*b**4*dm**2+11880.0_cp*a**2*b**4*dn**5 &
      &              -100320.0_cp*a**2*b**4*dn**4-492360.0_cp*a**2*b**4*dn**3+2476320.0_cp &
      &              *a**2*b**4*dn**2+12708960.0_cp*a**2*b**4*dn+13559040.0_cp*a**2*b**4 &
      &              +896.0_cp*a**2*b**2*dm**2*dn**3*r**2+10752.0_cp*a**2*b**2*dm**2*dn**2* &
      &              r**2+36736.0_cp*a**2*b**2*dm**2*dn*r**2+37632.0_cp*a**2*b**2*dm**2* &
      &              r**2-4320.0_cp*a**2*b**2*dn**5*r**2+36288.0_cp*a**2*b**2*dn**4*r**2 &
      &              +178400.0_cp*a**2*b**2*dn**3*r**2-897024.0_cp*a**2*b**2*dn**2*r**2-4604480.0_cp &
      &              *a**2*b**2*dn*r**2-4912512.0_cp*a**2*b**2*r**2-48.0_cp*a**2*dm**2*dn**3 &
      &              *r**4-576.0_cp*a**2*dm**2*dn**2*r**4-1968.0_cp*a**2*dm**2*dn*r**4-2016.0_cp &
      &              *a**2*dm**2*r**4+144.0_cp*a**2*dn**5*r**4-1200.0_cp*a**2*dn**4*r**4 &
      &              -5920.0_cp*a**2*dn**3*r**4+29664.0_cp*a**2*dn**2*r**4+152416.0_cp*a**2 &
      &              *dn*r**4+162624.0_cp*a**2*r**4-3360.0_cp*b**6*dm**2*dn**3-20160.0_cp &
      &              *b**6*dm**2*dn**2-36960.0_cp*b**6*dm**2*dn-20160.0_cp*b**6*dm**2+14784.0_cp &
      &              *b**6*dn**5-147840.0_cp*b**6*dn**4-310464.0_cp*b**6*dn**3+3163776.0_cp &
      &              *b**6*dn**2+8988672.0_cp*b**6*dn+5677056.0_cp*b**6+4480.0_cp*b**4*dm**2 &
      &              *dn**3*r**2+26880.0_cp*b**4*dm**2*dn**2*r**2+49280.0_cp*b**4*dm**2*dn &
      &              *r**2+26880.0_cp*b**4*dm**2*r**2-13440.0_cp*b**4*dn**5*r**2+133728.0_cp &
      &              *b**4*dn**4*r**2+281344.0_cp*b**4*dn**3*r**2-2864736.0_cp*b**4*dn** &
      &              2*r**2-8141056.0_cp*b**4*dn*r**2-5142144.0_cp*b**4*r**2-1440.0_cp*b**2 &
      &              *dm**2*dn**3*r**4-8640.0_cp*b**2*dm**2*dn**2*r**4-15840.0_cp*b**2*dm**2 &
      &              *dn*r**4-8640.0_cp*b**2*dm**2*r**4+2688.0_cp*b**2*dn**5*r**4-26544.0_cp &
      &              *b**2*dn**4*r**4-56160.0_cp*b**2*dn**3*r**4+568560.0_cp*b**2*dn**2* &
      &              r**4+1617312.0_cp*b**2*dn*r**4+1021824.0_cp*b**2*r**4+64.0_cp*dm**2 &
      &              *dn**3*r**6+384.0_cp*dm**2*dn**2*r**6+704.0_cp*dm**2*dn*r**6+384.0_cp &
      &              *dm**2*r**6-64.0_cp*dn**5*r**6+624.0_cp*dn**4*r**6+1344.0_cp*dn**3* &
      &              r**6-13296.0_cp*dn**2*r**6-37952.0_cp*dn*r**6-24000.0_cp*r**6)/(4096.0_cp &
      &              *dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn &
      &              +3.0_cp))
      stencil(ku+10)=-a**9*b*(15.0_cp*a**4*dm**2*dn**2+225.0_cp*a**4*dm** &
      &              2*dn+540.0_cp*a**4*dm**2-102.0_cp*a**4*dn**4+1233.0_cp*a**4*dn**3+2148.0_cp &
      &              *a**4*dn**2-35127.0_cp*a**4*dn-87480.0_cp*a**4+360.0_cp*a**2*b**2*dm**2 &
      &              *dn**2+2520.0_cp*a**2*b**2*dm**2*dn+3600.0_cp*a**2*b**2*dm**2-1540.0_cp &
      &              *a**2*b**2*dn**4+20790.0_cp*a**2*b**2*dn**3-5720.0_cp*a**2*b**2*dn**2- &
      &              436590.0_cp*a**2*b**2*dn-659340.0_cp*a**2*b**2-96.0_cp*a**2*dm**2*dn**2 &
      &              *r**2-672.0_cp*a**2*dm**2*dn*r**2-960.0_cp*a**2*dm**2*r**2+280.0_cp &
      &              *a**2*dn**4*r**2-3765.0_cp*a**2*dn**3*r**2+1034.0_cp*a**2*dn**2*r** &
      &              2+79089.0_cp*a**2*dn*r**2+119442.0_cp*a**2*r**2+1008.0_cp*b**4*dm** &
      &              2*dn**2+3024.0_cp*b**4*dm**2*dn+2016.0_cp*b**4*dm**2-3168.0_cp*b**4 &
      &              *dn**4+46992.0_cp*b**4*dn**3-88704.0_cp*b**4*dn**2-642576.0_cp*b**4 &
      &              *dn-503712.0_cp*b**4-896.0_cp*b**2*dm**2*dn**2*r**2-2688.0_cp*b**2*dm**2 &
      &              *dn*r**2-1792.0_cp*b**2*dm**2*r**2+1920.0_cp*b**2*dn**4*r**2-28368.0_cp &
      &              *b**2*dn**3*r**2+53536.0_cp*b**2*dn**2*r**2+387984.0_cp*b**2*dn*r** &
      &              2+304160.0_cp*b**2*r**2+144.0_cp*dm**2*dn**2*r**4+432.0_cp*dm**2*dn &
      &              *r**4+288.0_cp*dm**2*r**4-192.0_cp*dn**4*r**4+2820.0_cp*dn**3*r**4-5304.0_cp &
      &              *dn**2*r**4-38532.0_cp*dn*r**4-30216.0_cp*r**4)/(2048.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+11)=-a**10*(11.0_cp*a**4*dm**2*dn**2+117.0_cp*a**4*dm**2 &
      &              *dn+250.0_cp*a**4*dm**2-43.0_cp*a**4*dn**4+614.0_cp*a**4*dn**3+345.0_cp &
      &              *a**4*dn**2-17900.0_cp*a**4*dn-39500.0_cp*a**4+720.0_cp*a**2*b**2*dm**2 &
      &              *dn**2+4320.0_cp*a**2*b**2*dm**2*dn+5760.0_cp*a**2*b**2*dm**2-2112.0_cp &
      &              *a**2*b**2*dn**4+32736.0_cp*a**2*b**2*dn**3-33792.0_cp*a**2*b**2*dn**2 &
      &              -722304.0_cp*a**2*b**2*dn-1013760.0_cp*a**2*b**2-64.0_cp*a**2*dm**2 &
      &              *dn**2*r**2-384.0_cp*a**2*dm**2*dn*r**2-512.0_cp*a**2*dm**2*r**2+128.0_cp &
      &              *a**2*dn**4*r**2-1976.0_cp*a**2*dn**3*r**2+2040.0_cp*a**2*dn**2*r** &
      &              2+43616.0_cp*a**2*dn*r**2+61216.0_cp*a**2*r**2+3360.0_cp*b**4*dm**2 &
      &              *dn**2+10080.0_cp*b**4*dm**2*dn+6720.0_cp*b**4*dm**2-7920.0_cp*b**4 &
      &              *dn**4+132000.0_cp*b**4*dn**3-314160.0_cp*b**4*dn**2-1985280.0_cp*b**4 &
      &              *dn-1531200.0_cp*b**4-1792.0_cp*b**2*dm**2*dn**2*r**2-5376.0_cp*b** &
      &              2*dm**2*dn*r**2-3584.0_cp*b**2*dm**2*r**2+2880.0_cp*b**2*dn**4*r**2 &
      &              -47808.0_cp*b**2*dn**3*r**2+113792.0_cp*b**2*dn**2*r**2+719232.0_cp &
      &              *b**2*dn*r**2+554752.0_cp*b**2*r**2+96.0_cp*dm**2*dn**2*r**4+288.0_cp &
      &              *dm**2*dn*r**4+192.0_cp*dm**2*r**4-96.0_cp*dn**4*r**4+1584.0_cp*dn**3* &
      &              r**4-3760.0_cp*dn**2*r**4-23808.0_cp*dn*r**4-18368.0_cp*r**4)/(16384.0_cp &
      &              *dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+12)=a**11*b*(-25.0_cp*a**2*dm**2*dn-85.0_cp*a**2*dm**2+54.0_cp &
      &              *a**2*dn**3-1051.0_cp*a**2*dn**2+3712.0_cp*a**2*dn+14465.0_cp*a**2-240.0_cp &
      &              *b**2*dm**2*dn-240.0_cp*b**2*dm**2+440.0_cp*b**2*dn**3-9020.0_cp*b**2* &
      &              dn**2+41360.0_cp*b**2*dn+50820.0_cp*b**2+64.0_cp*dm**2*dn*r**2+64.0_cp &
      &              *dm**2*r**2-80.0_cp*dn**3*r**2+1634.0_cp*dn**2*r**2-7492.0_cp*dn*r**2- &
      &              9206.0_cp*r**2)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn &
      &              +1.0_cp))
      stencil(ku+13)=a**12*(-3.0_cp*a**2*dm**2*dn-9.0_cp*a**2*dm**2+5.0_cp &
      &              *a**2*dn**3-107.0_cp*a**2*dn**2+440.0_cp*a**2*dn+1488.0_cp*a**2-90.0_cp &
      &              *b**2*dm**2*dn-90.0_cp*b**2*dm**2+132.0_cp*b**2*dn**3-2948.0_cp*b** &
      &              2*dn**2+14872.0_cp*b**2*dn+17952.0_cp*b**2+8.0_cp*dm**2*dn*r**2+8.0_cp &
      &              *dm**2*r**2-8.0_cp*dn**3*r**2+178.0_cp*dn**2*r**2-898.0_cp*dn*r**2-1084.0_cp &
      &              *r**2)/(8192.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              ))
      stencil(ku+14)=-a**13*b*(5.0_cp*dm**2-6.0_cp*dn**2+151.0_cp*dn-949.0_cp &
      &              )/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp))
      stencil(ku+15)=-a**14*(dm**2-dn**2+27.0_cp*dn-182.0_cp)/(16384.0_cp &
      &              *dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp))
      stencil(ku+16:len_stencil) = 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4rmult4hmult8laplrot
!------------------------------------------------------------------------------
   function intcheb4rmult4hmult8laplrot2(a, b, r, m, n, len_stencil) result(stencil)

      !-- Input variables
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
      real(cp), intent(in) :: r 
      integer,  intent(in) :: m 
      integer,  intent(in) :: n
      integer,  intent(in) :: len_stencil

      !-- Output variable 
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      real(cp) :: dm
      real(cp) :: dn
      integer :: ku

      dm = real(m, cp)
      dn = real(n, cp)
      ku = (len_stencil-1)/2

      stencil(1:ku-12) = 0.0_cp
      stencil(ku-11)=a**12*(dm-dn-14.0_cp)*(dm+dn+14.0_cp)*(dm**2-dn**2-23.0_cp &
      &              *dn-132.0_cp)/(4096.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-10)=a**11*b*(8.0_cp*dm**4-20.0_cp*dm**2*dn**2-469.0_cp*dm**2 &
      &              *dn-2771.0_cp*dm**2+12.0_cp*dn**4+567.0_cp*dn**3+10012.0_cp*dn**2+78299.0_cp &
      &              *dn+228822.0_cp)/(2048.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-9)=-a**10*(-2.0_cp*a**2*dm**4*dn+8.0_cp*a**2*dm**4+8.0_cp &
      &              *a**2*dm**2*dn**3+153.0_cp*a**2*dm**2*dn**2+507.0_cp*a**2*dm**2*dn-2342.0_cp &
      &              *a**2*dm**2-6.0_cp*a**2*dn**5-255.0_cp*a**2*dn**4-3919.0_cp*a**2*dn**3 &
      &              -24204.0_cp*a**2*dn**2-28316.0_cp*a**2*dn+168240.0_cp*a**2-56.0_cp* &
      &              b**2*dm**4*dn+56.0_cp*b**2*dm**4+180.0_cp*b**2*dm**2*dn**3+3672.0_cp &
      &              *b**2*dm**2*dn**2+16892.0_cp*b**2*dm**2*dn-20744.0_cp*b**2*dm**2-132.0_cp &
      &              *b**2*dn**5-5610.0_cp*b**2*dn**4-87600.0_cp*b**2*dn**3-578706.0_cp* &
      &              b**2*dn**2-1136232.0_cp*b**2*dn+1808280.0_cp*b**2+8.0_cp*dm**4*dn*r**2 &
      &              -8.0_cp*dm**4*r**2-16.0_cp*dm**2*dn**3*r**2-330.0_cp*dm**2*dn**2*r**2- &
      &              1550.0_cp*dm**2*dn*r**2+1896.0_cp*dm**2*r**2+8.0_cp*dn**5*r**2+338.0_cp &
      &              *dn**4*r**2+5254.0_cp*dn**3*r**2+34632.0_cp*dn**2*r**2+68152.0_cp*dn &
      &              *r**2-108384.0_cp*r**2)/(2048.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp &
      &              )*(dn+3.0_cp))
      stencil(ku-8)=-a**9*b*(-24.0_cp*a**2*dm**4*dn+120.0_cp*a**2*dm**4+140.0_cp &
      &              *a**2*dm**2*dn**3+2365.0_cp*a**2*dm**2*dn**2+6048.0_cp*a**2*dm**2*dn &
      &              -36417.0_cp*a**2*dm**2-132.0_cp*a**2*dn**5-5115.0_cp*a**2*dn**4-71305.0_cp &
      &              *a**2*dn**3-392435.0_cp*a**2*dn**2-318219.0_cp*a**2*dn+2706966.0_cp &
      &              *a**2-224.0_cp*b**2*dm**4*dn+224.0_cp*b**2*dm**4+960.0_cp*b**2*dm** &
      &              2*dn**3+17616.0_cp*b**2*dm**2*dn**2+71744.0_cp*b**2*dm**2*dn-90320.0_cp &
      &              *b**2*dm**2-880.0_cp*b**2*dn**5-34100.0_cp*b**2*dn**4-484580.0_cp*b**2 &
      &              *dn**3-2898148.0_cp*b**2*dn**2-4984044.0_cp*b**2*dn+8401752.0_cp*b**2+ &
      &              96.0_cp*dm**4*dn*r**2-96.0_cp*dm**4*r**2-256.0_cp*dm**2*dn**3*r**2-4748.0_cp &
      &              *dm**2*dn**2*r**2-19728.0_cp*dm**2*dn*r**2+24732.0_cp*dm**2*r**2+160.0_cp &
      &              *dn**5*r**2+6164.0_cp*dn**4*r**2+87204.0_cp*dn**3*r**2+520360.0_cp*dn**2 &
      &              *r**2+896540.0_cp*dn*r**2-1510428.0_cp*r**2)/(2048.0_cp*dn*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-7)=a**8*(a**4*dm**4*dn**2-33.0_cp*a**4*dm**4*dn+92.0_cp* &
      &              a**4*dm**4-26.0_cp*a**4*dm**2*dn**4-280.0_cp*a**4*dm**2*dn**3+904.0_cp &
      &              *a**4*dm**2*dn**2+8842.0_cp*a**4*dm**2*dn-25664.0_cp*a**4*dm**2+33.0_cp &
      &              *a**4*dn**6+1089.0_cp*a**4*dn**5+11788.0_cp*a**4*dn**4+31790.0_cp*a**4 &
      &              *dn**3-176752.0_cp*a**4*dn**2-572264.0_cp*a**4*dn+1760448.0_cp*a**4 &
      &              +112.0_cp*a**2*b**2*dm**4*dn**2-1008.0_cp*a**2*b**2*dm**4*dn+1568.0_cp &
      &              *a**2*b**2*dm**4-1080.0_cp*a**2*b**2*dm**2*dn**4-13608.0_cp*a**2*b**2* &
      &              dm**2*dn**3+4736.0_cp*a**2*b**2*dm**2*dn**2+302112.0_cp*a**2*b**2*dm**2 &
      &              *dn-497024.0_cp*a**2*b**2*dm**2+1320.0_cp*a**2*b**2*dn**6+43560.0_cp &
      &              *a**2*b**2*dn**5+485340.0_cp*a**2*b**2*dn**4+1627524.0_cp*a**2*b**2 &
      &              *dn**3-4433232.0_cp*a**2*b**2*dn**2-21494160.0_cp*a**2*b**2*dn+38457216.0_cp &
      &              *a**2*b**2-16.0_cp*a**2*dm**4*dn**2*r**2+144.0_cp*a**2*dm**4*dn*r** &
      &              2-224.0_cp*a**2*dm**4*r**2+96.0_cp*a**2*dm**2*dn**4*r**2+1224.0_cp* &
      &              a**2*dm**2*dn**3*r**2-360.0_cp*a**2*dm**2*dn**2*r**2-27720.0_cp*a** &
      &              2*dm**2*dn*r**2+45552.0_cp*a**2*dm**2*r**2-80.0_cp*a**2*dn**6*r**2-2624.0_cp &
      &              *a**2*dn**5*r**2-29096.0_cp*a**2*dn**4*r**2-97264.0_cp*a**2*dn**3*r**2 &
      &              +265120.0_cp*a**2*dn**2*r**2+1289280.0_cp*a**2*dn*r**2-2306304.0_cp &
      &              *a**2*r**2+560.0_cp*b**4*dm**4*dn**2-1680.0_cp*b**4*dm**4*dn+1120.0_cp &
      &              *b**4*dm**4-3360.0_cp*b**4*dm**2*dn**4-48048.0_cp*b**4*dm**2*dn**3-84560.0_cp &
      &              *b**4*dm**2*dn**2+640416.0_cp*b**4*dm**2*dn-504448.0_cp*b**4*dm**2+3960.0_cp &
      &              *b**4*dn**6+130680.0_cp*b**4*dn**5+1497480.0_cp*b**4*dn**4+5950344.0_cp &
      &              *b**4*dn**3-5013600.0_cp*b**4*dn**2-53374368.0_cp*b**4*dn+50805504.0_cp &
      &              *b**4-480.0_cp*b**2*dm**4*dn**2*r**2+1440.0_cp*b**2*dm**4*dn*r**2-960.0_cp &
      &              *b**2*dm**4*r**2+1792.0_cp*b**2*dm**2*dn**4*r**2+25928.0_cp*b**2*dm**2 &
      &              *dn**3*r**2+47576.0_cp*b**2*dm**2*dn**2*r**2-351104.0_cp*b**2*dm**2 &
      &              *dn*r**2+275808.0_cp*b**2*dm**2*r**2-1440.0_cp*b**2*dn**6*r**2-47232.0_cp &
      &              *b**2*dn**5*r**2-538624.0_cp*b**2*dn**4*r**2-2134040.0_cp*b**2*dn** &
      &              3*r**2+1793896.0_cp*b**2*dn**2*r**2+19190864.0_cp*b**2*dn*r**2-18263424.0_cp &
      &              *b**2*r**2+48.0_cp*dm**4*dn**2*r**4-144.0_cp*dm**4*dn*r**4+96.0_cp*dm**4 &
      &              *r**4-96.0_cp*dm**2*dn**4*r**4-1416.0_cp*dm**2*dn**3*r**4-2792.0_cp &
      &              *dm**2*dn**2*r**4+19728.0_cp*dm**2*dn*r**4-15424.0_cp*dm**2*r**4+48.0_cp &
      &              *dn**6*r**4+1560.0_cp*dn**5*r**4+17640.0_cp*dn**4*r**4+69232.0_cp*dn**3 &
      &              *r**4-60000.0_cp*dn**2*r**4-620608.0_cp*dn*r**4+592128.0_cp*r**4)/(2048.0_cp &
      &              *dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-6)=a**7*b*(-8.0_cp*a**4*dm**4*dn**2-360.0_cp*a**4*dm**4*dn &
      &              +1232.0_cp*a**4*dm**4-380.0_cp*a**4*dm**2*dn**4-3023.0_cp*a**4*dm** &
      &              2*dn**3+18448.0_cp*a**4*dm**2*dn**2+98177.0_cp*a**4*dm**2*dn-354662.0_cp &
      &              *a**4*dm**2+660.0_cp*a**4*dn**6+19305.0_cp*a**4*dn**5+179965.0_cp*a**4 &
      &              *dn**4+323119.0_cp*a**4*dn**3-2953993.0_cp*a**4*dn**2-6458116.0_cp* &
      &              a**4*dn+25083492.0_cp*a**4+224.0_cp*a**2*b**2*dm**4*dn**2-3360.0_cp &
      &              *a**2*b**2*dm**4*dn+5824.0_cp*a**2*b**2*dm**4-4800.0_cp*a**2*b**2*dm**2 &
      &              *dn**4-49008.0_cp*a**2*b**2*dm**2*dn**3+78784.0_cp*a**2*b**2*dm**2*dn**2 &
      &              +1051152.0_cp*a**2*b**2*dm**2*dn-1948576.0_cp*a**2*b**2*dm**2+7920.0_cp &
      &              *a**2*b**2*dn**6+231660.0_cp*a**2*b**2*dn**5+2242500.0_cp*a**2*b**2 &
      &              *dn**4+5777124.0_cp*a**2*b**2*dn**3-23298180.0_cp*a**2*b**2*dn**2-77497296.0_cp &
      &              *a**2*b**2*dn+158170320.0_cp*a**2*b**2-96.0_cp*a**2*dm**4*dn**2*r** &
      &              2+1440.0_cp*a**2*dm**4*dn*r**2-2496.0_cp*a**2*dm**4*r**2+1280.0_cp* &
      &              a**2*dm**2*dn**4*r**2+13220.0_cp*a**2*dm**2*dn**3*r**2-20816.0_cp*a**2 &
      &              *dm**2*dn**2*r**2-289100.0_cp*a**2*dm**2*dn*r**2+535224.0_cp*a**2*dm**2 &
      &              *r**2-1440.0_cp*a**2*dn**6*r**2-41868.0_cp*a**2*dn**5*r**2-403364.0_cp &
      &              *a**2*dn**4*r**2-1035344.0_cp*a**2*dn**3*r**2+4183700.0_cp*a**2*dn**2* &
      &              r**2+13942244.0_cp*a**2*dn*r**2-28450776.0_cp*a**2*r**2+896.0_cp*b**4* &
      &              dm**4*dn**2-2688.0_cp*b**4*dm**4*dn+1792.0_cp*b**4*dm**4-8064.0_cp* &
      &              b**4*dm**2*dn**4-98784.0_cp*b**4*dm**2*dn**3-116480.0_cp*b**4*dm**2 &
      &              *dn**2+1161888.0_cp*b**4*dm**2*dn-938560.0_cp*b**4*dm**2+12672.0_cp &
      &              *b**4*dn**6+370656.0_cp*b**4*dn**5+3720672.0_cp*b**4*dn**4+12282912.0_cp &
      &              *b**4*dn**3-16635744.0_cp*b**4*dn**2-104171904.0_cp*b**4*dn+104420736.0_cp &
      &              *b**4-1280.0_cp*b**2*dm**4*dn**2*r**2+3840.0_cp*b**2*dm**4*dn*r**2-2560.0_cp &
      &              *b**2*dm**4*r**2+7168.0_cp*b**2*dm**2*dn**4*r**2+88816.0_cp*b**2*dm**2 &
      &              *dn**3*r**2+110144.0_cp*b**2*dm**2*dn**2*r**2-1059664.0_cp*b**2*dm**2* &
      &              dn*r**2+853536.0_cp*b**2*dm**2*r**2-7680.0_cp*b**2*dn**6*r**2-223296.0_cp &
      &              *b**2*dn**5*r**2-2230720.0_cp*b**2*dn**4*r**2-7340304.0_cp*b**2*dn**3* &
      &              r**2+9948544.0_cp*b**2*dn**2*r**2+62403312.0_cp*b**2*dn*r**2-62549856.0_cp &
      &              *b**2*r**2+384.0_cp*dm**4*dn**2*r**4-1152.0_cp*dm**4*dn*r**4+768.0_cp &
      &              *dm**4*r**4-1152.0_cp*dm**2*dn**4*r**4-14544.0_cp*dm**2*dn**3*r**4-19648.0_cp &
      &              *dm**2*dn**2*r**4+178032.0_cp*dm**2*dn*r**4-142688.0_cp*dm**2*r**4+768.0_cp &
      &              *dn**6*r**4+22128.0_cp*dn**5*r**4+219216.0_cp*dn**4*r**4+714208.0_cp &
      &              *dn**3*r**4-987600.0_cp*dn**2*r**4-6059248.0_cp*dn*r**4+6090528.0_cp &
      &              *r**4)/(2048.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp &
      &              )*(dn+3.0_cp))
      stencil(ku-5)=-a**6*(6.0_cp*a**6*dm**4*dn**3+6.0_cp*a**6*dm**4*dn** &
      &              2-324.0_cp*a**6*dm**4*dn+816.0_cp*a**6*dm**4+40.0_cp*a**6*dm**2*dn**5- &
      &              69.0_cp*a**6*dm**2*dn**4-4138.0_cp*a**6*dm**2*dn**3+6627.0_cp*a**6*dm**2 &
      &              *dn**2+86112.0_cp*a**6*dm**2*dn-221964.0_cp*a**6*dm**2-110.0_cp*a** &
      &              6*dn**7-2475.0_cp*a**6*dn**6-12280.0_cp*a**6*dn**5+68407.0_cp*a**6*dn**4 &
      &              +481422.0_cp*a**6*dn**3-1221508.0_cp*a**6*dn**2-5610168.0_cp*a**6*dn &
      &              +14836320.0_cp*a**6+168.0_cp*a**4*b**2*dm**4*dn**3+1008.0_cp*a**4*b**2 &
      &              *dm**4*dn**2-11592.0_cp*a**4*b**2*dm**4*dn+21168.0_cp*a**4*b**2*dm**4+ &
      &              2340.0_cp*a**4*b**2*dm**2*dn**5+3996.0_cp*a**4*b**2*dm**2*dn**4-177972.0_cp &
      &              *a**4*b**2*dm**2*dn**3+19188.0_cp*a**4*b**2*dm**2*dn**2+3357648.0_cp &
      &              *a**4*b**2*dm**2*dn-6332688.0_cp*a**4*b**2*dm**2-5940.0_cp*a**4*b** &
      &              2*dn**7-133650.0_cp*a**4*b**2*dn**6-725310.0_cp*a**4*b**2*dn**5+2632662.0_cp &
      &              *a**4*b**2*dn**4+23611050.0_cp*a**4*b**2*dn**3-35194500.0_cp*a**4*b**2 &
      &              *dn**2-237236040.0_cp*a**4*b**2*dn+464386608.0_cp*a**4*b**2-24.0_cp &
      &              *a**4*dm**4*dn**3*r**2-144.0_cp*a**4*dm**4*dn**2*r**2+1656.0_cp*a** &
      &              4*dm**4*dn*r**2-3024.0_cp*a**4*dm**4*r**2-208.0_cp*a**4*dm**2*dn**5 &
      &              *r**2-366.0_cp*a**4*dm**2*dn**4*r**2+16084.0_cp*a**4*dm**2*dn**3*r**2- &
      &              1146.0_cp*a**4*dm**2*dn**2*r**2-308268.0_cp*a**4*dm**2*dn*r**2+581040.0_cp &
      &              *a**4*dm**2*r**2+360.0_cp*a**4*dn**7*r**2+8046.0_cp*a**4*dn**6*r**2 &
      &              +43372.0_cp*a**4*dn**5*r**2-158334.0_cp*a**4*dn**4*r**2-1413052.0_cp &
      &              *a**4*dn**3*r**2+2104872.0_cp*a**4*dn**2*r**2+14232816.0_cp*a**4*dn &
      &              *r**2-27857088.0_cp*a**4*r**2+6720.0_cp*a**2*b**4*dm**4*dn**2-33600.0_cp &
      &              *a**2*b**4*dm**4*dn+40320.0_cp*a**2*b**4*dm**4+13440.0_cp*a**2*b**4 &
      &              *dm**2*dn**5+62496.0_cp*a**2*b**4*dm**2*dn**4-670656.0_cp*a**2*b**4 &
      &              *dm**2*dn**3-1117536.0_cp*a**2*b**4*dm**2*dn**2+11432064.0_cp*a**2* &
      &              b**4*dm**2*dn-14458752.0_cp*a**2*b**4*dm**2-31680.0_cp*a**2*b**4*dn**7 &
      &              -712800.0_cp*a**2*b**4*dn**6-4200000.0_cp*a**2*b**4*dn**5+8380512.0_cp &
      &              *a**2*b**4*dn**4+110197248.0_cp*a**2*b**4*dn**3-58367232.0_cp*a**2* &
      &              b**4*dn**2-921203712.0_cp*a**2*b**4*dn+1244284416.0_cp*a**2*b**4-5760.0_cp &
      &              *a**2*b**2*dm**4*dn**2*r**2+28800.0_cp*a**2*b**2*dm**4*dn*r**2-34560.0_cp &
      &              *a**2*b**2*dm**4*r**2-7168.0_cp*a**2*b**2*dm**2*dn**5*r**2-33936.0_cp &
      &              *a**2*b**2*dm**2*dn**4*r**2+361312.0_cp*a**2*b**2*dm**2*dn**3*r**2+626352.0_cp &
      &              *a**2*b**2*dm**2*dn**2*r**2-6278400.0_cp*a**2*b**2*dm**2*dn*r**2+7933248.0_cp &
      &              *a**2*b**2*dm**2*r**2+11520.0_cp*a**2*b**2*dn**7*r**2+257472.0_cp*a**2 &
      &              *b**2*dn**6*r**2+1507072.0_cp*a**2*b**2*dn**5*r**2-3035088.0_cp*a** &
      &              2*b**2*dn**4*r**2-39556384.0_cp*a**2*b**2*dn**3*r**2+20913264.0_cp* &
      &              a**2*b**2*dn**2*r**2+331354944.0_cp*a**2*b**2*dn*r**2-447529536.0_cp &
      &              *a**2*b**2*r**2+576.0_cp*a**2*dm**4*dn**2*r**4-2880.0_cp*a**2*dm**4 &
      &              *dn*r**4+3456.0_cp*a**2*dm**4*r**4+384.0_cp*a**2*dm**2*dn**5*r**4+1872.0_cp &
      &              *a**2*dm**2*dn**4*r**4-19680.0_cp*a**2*dm**2*dn**3*r**4-36528.0_cp* &
      &              a**2*dm**2*dn**2*r**4+353856.0_cp*a**2*dm**2*dn*r**4-446400.0_cp*a**2* &
      &              dm**2*r**4-384.0_cp*a**2*dn**7*r**4-8496.0_cp*a**2*dn**6*r**4-49152.0_cp &
      &              *a**2*dn**5*r**4+101200.0_cp*a**2*dn**4*r**4+1287264.0_cp*a**2*dn** &
      &              3*r**4-702400.0_cp*a**2*dn**2*r**4-10712448.0_cp*a**2*dn*r**4+14482944.0_cp &
      &              *a**2*r**4-896.0_cp*b**6*dm**4*dn**3+5376.0_cp*b**6*dm**4*dn**2-9856.0_cp &
      &              *b**6*dm**4*dn+5376.0_cp*b**6*dm**4+13440.0_cp*b**6*dm**2*dn**5+96768.0_cp &
      &              *b**6*dm**2*dn**4-332416.0_cp*b**6*dm**2*dn**3-1634304.0_cp*b**6*dm**2 &
      &              *dn**2+5361664.0_cp*b**6*dm**2*dn-3505152.0_cp*b**6*dm**2-29568.0_cp &
      &              *b**6*dn**7-665280.0_cp*b**6*dn**6-4229568.0_cp*b**6*dn**5+2538816.0_cp &
      &              *b**6*dn**4+85366848.0_cp*b**6*dn**3+33801600.0_cp*b**6*dn**2-560270592.0_cp &
      &              *b**6*dn+443487744.0_cp*b**6+1920.0_cp*b**4*dm**4*dn**3*r**2-11520.0_cp &
      &              *b**4*dm**4*dn**2*r**2+21120.0_cp*b**4*dm**4*dn*r**2-11520.0_cp*b** &
      &              4*dm**4*r**2-17920.0_cp*b**4*dm**2*dn**5*r**2-131040.0_cp*b**4*dm** &
      &              2*dn**4*r**2+439360.0_cp*b**4*dm**2*dn**3*r**2+2252640.0_cp*b**4*dm**2 &
      &              *dn**2*r**2-7312320.0_cp*b**4*dm**2*dn*r**2+4769280.0_cp*b**4*dm**2 &
      &              *r**2+26880.0_cp*b**4*dn**7*r**2+600768.0_cp*b**4*dn**6*r**2+3794560.0_cp &
      &              *b**4*dn**5*r**2-2341472.0_cp*b**4*dn**4*r**2-76570240.0_cp*b**4*dn**3 &
      &              *r**2-30349408.0_cp*b**4*dn**2*r**2+503268480.0_cp*b**4*dn*r**2-398429568.0_cp &
      &              *b**4*r**2-1152.0_cp*b**2*dm**4*dn**3*r**4+6912.0_cp*b**2*dm**4*dn**2* &
      &              r**4-12672.0_cp*b**2*dm**4*dn*r**4+6912.0_cp*b**2*dm**4*r**4+5760.0_cp &
      &              *b**2*dm**2*dn**5*r**4+43200.0_cp*b**2*dm**2*dn**4*r**4-138624.0_cp &
      &              *b**2*dm**2*dn**3*r**4-766656.0_cp*b**2*dm**2*dn**2*r**4+2443776.0_cp &
      &              *b**2*dm**2*dn*r**4-1587456.0_cp*b**2*dm**2*r**4-5376.0_cp*b**2*dn**7* &
      &              r**4-118944.0_cp*b**2*dn**6*r**4-742656.0_cp*b**2*dn**5*r**4+492416.0_cp &
      &              *b**2*dn**4*r**4+14986368.0_cp*b**2*dn**3*r**4+5605792.0_cp*b**2*dn**2 &
      &              *r**4-97909632.0_cp*b**2*dn*r**4+77692032.0_cp*b**2*r**4+128.0_cp*dm**4 &
      &              *dn**3*r**6-768.0_cp*dm**4*dn**2*r**6+1408.0_cp*dm**4*dn*r**6-768.0_cp &
      &              *dm**4*r**6-256.0_cp*dm**2*dn**5*r**6-2016.0_cp*dm**2*dn**4*r**6+5824.0_cp &
      &              *dm**2*dn**3*r**6+38496.0_cp*dm**2*dn**2*r**6-118080.0_cp*dm**2*dn* &
      &              r**6+76032.0_cp*dm**2*r**6+128.0_cp*dn**7*r**6+2784.0_cp*dn**6*r**6 &
      &              +16960.0_cp*dn**5*r**6-13536.0_cp*dn**4*r**6-341568.0_cp*dn**3*r**6 &
      &              -100224.0_cp*dn**2*r**6+2177280.0_cp*dn*r**6-1741824.0_cp*r**6)/(2048.0_cp &
      &              *dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn &
      &              +3.0_cp))
      stencil(ku-4)=-a**5*b*(88.0_cp*a**6*dm**4*dn**3-240.0_cp*a**6*dm**4 &
      &              *dn**2-2872.0_cp*a**6*dm**4*dn+9360.0_cp*a**6*dm**4+420.0_cp*a**6*dm**2 &
      &              *dn**5-3345.0_cp*a**6*dm**2*dn**4-45691.0_cp*a**6*dm**2*dn**3+167295.0_cp &
      &              *a**6*dm**2*dn**2+778639.0_cp*a**6*dm**2*dn-2624310.0_cp*a**6*dm**2 &
      &              -1980.0_cp*a**6*dn**7-37125.0_cp*a**6*dn**6-117390.0_cp*a**6*dn**5+1229160.0_cp &
      &              *a**6*dn**4+4964682.0_cp*a**6*dn**3-21161055.0_cp*a**6*dn**2-51851088.0_cp &
      &              *a**6*dn+180678060.0_cp*a**6+896.0_cp*a**4*b**2*dm**4*dn**3-30464.0_cp &
      &              *a**4*b**2*dm**4*dn+67200.0_cp*a**4*b**2*dm**4+7680.0_cp*a**4*b**2*dm**2 &
      &              *dn**5-19200.0_cp*a**4*b**2*dm**2*dn**4-603392.0_cp*a**4*b**2*dm**2 &
      &              *dn**3+1076160.0_cp*a**4*b**2*dm**2*dn**2+9131648.0_cp*a**4*b**2*dm**2 &
      &              *dn-21099840.0_cp*a**4*b**2*dm**2-31680.0_cp*a**4*b**2*dn**7-594000.0_cp &
      &              *a**4*b**2*dn**6-2209920.0_cp*a**4*b**2*dn**5+14949600.0_cp*a**4*b**2* &
      &              dn**4+75004224.0_cp*a**4*b**2*dn**3-200982480.0_cp*a**4*b**2*dn**2-668088576.0_cp &
      &              *a**4*b**2*dn+1616397120.0_cp*a**4*b**2-384.0_cp*a**4*dm**4*dn**3*r**2 &
      &              +13056.0_cp*a**4*dm**4*dn*r**2-28800.0_cp*a**4*dm**4*r**2-2048.0_cp &
      &              *a**4*dm**2*dn**5*r**2+5120.0_cp*a**4*dm**2*dn**4*r**2+164096.0_cp* &
      &              a**4*dm**2*dn**3*r**2-290000.0_cp*a**4*dm**2*dn**2*r**2-2513376.0_cp &
      &              *a**4*dm**2*dn*r**2+5802480.0_cp*a**4*dm**2*r**2+5760.0_cp*a**4*dn**7* &
      &              r**2+107280.0_cp*a**4*dn**6*r**2+396032.0_cp*a**4*dn**5*r**2-2693600.0_cp &
      &              *a**4*dn**4*r**2-13469504.0_cp*a**4*dn**3*r**2+36108800.0_cp*a**4*dn**2 &
      &              *r**2+120220704.0_cp*a**4*dn*r**2-290823120.0_cp*a**4*r**2+896.0_cp &
      &              *a**2*b**4*dm**4*dn**3+5376.0_cp*a**2*b**4*dm**4*dn**2-43904.0_cp*a**2 &
      &              *b**4*dm**4*dn+59136.0_cp*a**2*b**4*dm**4+24192.0_cp*a**2*b**4*dm** &
      &              2*dn**5+42336.0_cp*a**2*b**4*dm**2*dn**4-1231328.0_cp*a**2*b**4*dm**2* &
      &              dn**3-100128.0_cp*a**2*b**4*dm**2*dn**2+16039520.0_cp*a**2*b**4*dm**2* &
      &              dn-23279424.0_cp*a**2*b**4*dm**2-88704.0_cp*a**2*b**4*dn**7-1663200.0_cp &
      &              *a**2*b**4*dn**6-7116480.0_cp*a**2*b**4*dn**5+28651392.0_cp*a**2*b**4* &
      &              dn**4+189193536.0_cp*a**2*b**4*dn**3-258578208.0_cp*a**2*b**4*dn**2 &
      &              -1373879808.0_cp*a**2*b**4*dn+2155628160.0_cp*a**2*b**4-1280.0_cp*a**2 &
      &              *b**2*dm**4*dn**3*r**2-7680.0_cp*a**2*b**2*dm**4*dn**2*r**2+62720.0_cp &
      &              *a**2*b**2*dm**4*dn*r**2-84480.0_cp*a**2*b**2*dm**4*r**2-21504.0_cp &
      &              *a**2*b**2*dm**2*dn**5*r**2-38640.0_cp*a**2*b**2*dm**2*dn**4*r**2+1110192.0_cp &
      &              *a**2*b**2*dm**2*dn**3*r**2+123600.0_cp*a**2*b**2*dm**2*dn**2*r**2-14662704.0_cp &
      &              *a**2*b**2*dm**2*dn*r**2+21255840.0_cp*a**2*b**2*dm**2*r**2+53760.0_cp &
      &              *a**2*b**2*dn**7*r**2+1001280.0_cp*a**2*b**2*dn**6*r**2+4252416.0_cp &
      &              *a**2*b**2*dn**5*r**2-17239600.0_cp*a**2*b**2*dn**4*r**2-113215248.0_cp &
      &              *a**2*b**2*dn**3*r**2+154887760.0_cp*a**2*b**2*dn**2*r**2+823422096.0_cp &
      &              *a**2*b**2*dn*r**2-1291880160.0_cp*a**2*b**2*r**2+384.0_cp*a**2*dm**4* &
      &              dn**3*r**4+2304.0_cp*a**2*dm**4*dn**2*r**4-18816.0_cp*a**2*dm**4*dn &
      &              *r**4+25344.0_cp*a**2*dm**4*r**4+3456.0_cp*a**2*dm**2*dn**5*r**4+6480.0_cp &
      &              *a**2*dm**2*dn**4*r**4-182800.0_cp*a**2*dm**2*dn**3*r**4-30192.0_cp &
      &              *a**2*dm**2*dn**2*r**4+2473744.0_cp*a**2*dm**2*dn*r**4-3578592.0_cp &
      &              *a**2*dm**2*r**4-5376.0_cp*a**2*dn**7*r**4-99120.0_cp*a**2*dn**6*r**4- &
      &              415488.0_cp*a**2*dn**5*r**4+1707920.0_cp*a**2*dn**4*r**4+11047248.0_cp &
      &              *a**2*dn**3*r**4-15267584.0_cp*a**2*dn**2*r**4-79922640.0_cp*a**2*dn &
      &              *r**4+125538336.0_cp*a**2*r**4-512.0_cp*b**6*dm**4*dn**3+3072.0_cp* &
      &              b**6*dm**4*dn**2-5632.0_cp*b**6*dm**4*dn+3072.0_cp*b**6*dm**4+15360.0_cp &
      &              *b**6*dm**2*dn**5+79104.0_cp*b**6*dm**2*dn**4-384256.0_cp*b**6*dm** &
      &              2*dn**3-1054464.0_cp*b**6*dm**2*dn**2+4190464.0_cp*b**6*dm**2*dn-2846208.0_cp &
      &              *b**6*dm**2-50688.0_cp*b**6*dn**7-950400.0_cp*b**6*dn**6-4597248.0_cp &
      &              *b**6*dn**5+8825088.0_cp*b**6*dn**4+91407360.0_cp*b**6*dn**3-20284032.0_cp &
      &              *b**6*dn**2-502548480.0_cp*b**6*dn+428198400.0_cp*b**6+1536.0_cp*b**4* &
      &              dm**4*dn**3*r**2-9216.0_cp*b**4*dm**4*dn**2*r**2+16896.0_cp*b**4*dm**4 &
      &              *dn*r**2-9216.0_cp*b**4*dm**4*r**2-28672.0_cp*b**4*dm**2*dn**5*r**2 &
      &              -150080.0_cp*b**4*dm**2*dn**4*r**2+716608.0_cp*b**4*dm**2*dn**3*r** &
      &              2+2032832.0_cp*b**4*dm**2*dn**2*r**2-7974720.0_cp*b**4*dm**2*dn*r** &
      &              2+5404032.0_cp*b**4*dm**2*r**2+64512.0_cp*b**4*dn**7*r**2+1201536.0_cp &
      &              *b**4*dn**6*r**2+5770240.0_cp*b**4*dn**5*r**2-11206720.0_cp*b**4*dn**4 &
      &              *r**2-114814144.0_cp*b**4*dn**3*r**2+25750336.0_cp*b**4*dn**2*r**2+631781568.0_cp &
      &              *b**4*dn*r**2-538547328.0_cp*b**4*r**2-1536.0_cp*b**2*dm**4*dn**3*r**4 &
      &              +9216.0_cp*b**2*dm**4*dn**2*r**4-16896.0_cp*b**2*dm**4*dn*r**4+9216.0_cp &
      &              *b**2*dm**4*r**4+15360.0_cp*b**2*dm**2*dn**5*r**4+82560.0_cp*b**2*dm**2 &
      &              *dn**4*r**4-382592.0_cp*b**2*dm**2*dn**3*r**4-1150848.0_cp*b**2*dm**2* &
      &              dn**2*r**4+4416128.0_cp*b**2*dm**2*dn*r**4-2980608.0_cp*b**2*dm**2* &
      &              r**4-21504.0_cp*b**2*dn**7*r**4-396480.0_cp*b**2*dn**6*r**4-1880064.0_cp &
      &              *b**2*dn**5*r**4+3750400.0_cp*b**2*dn**4*r**4+37438080.0_cp*b**2*dn**3 &
      &              *r**4-9069376.0_cp*b**2*dn**2*r**4-205123200.0_cp*b**2*dn*r**4+175302144.0_cp &
      &              *b**2*r**4+512.0_cp*dm**4*dn**3*r**6-3072.0_cp*dm**4*dn**2*r**6+5632.0_cp &
      &              *dm**4*dn*r**6-3072.0_cp*dm**4*r**6-2048.0_cp*dm**2*dn**5*r**6-11584.0_cp &
      &              *dm**2*dn**4*r**6+50240.0_cp*dm**2*dn**3*r**6+172480.0_cp*dm**2*dn**2* &
      &              r**6-631872.0_cp*dm**2*dn*r**6+422784.0_cp*dm**2*r**6+1536.0_cp*dn**7* &
      &              r**6+27840.0_cp*dn**6*r**6+128512.0_cp*dn**5*r**6-274624.0_cp*dn**4 &
      &              *r**6-2556352.0_cp*dn**3*r**6+774016.0_cp*dn**2*r**6+13758912.0_cp*dn &
      &              *r**6-11859840.0_cp*r**6)/(2048.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn &
      &              -1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-3)=a**4*(-17.0_cp*a**8*dm**4*dn**3+162.0_cp*a**8*dm**4*dn**2 &
      &              +293.0_cp*a**8*dm**4*dn-2718.0_cp*a**8*dm**4-30.0_cp*a**8*dm**2*dn**5+ &
      &              1437.0_cp*a**8*dm**2*dn**4+5516.0_cp*a**8*dm**2*dn**3-66921.0_cp*a**8* &
      &              dm**2*dn**2-74234.0_cp*a**8*dm**2*dn+729072.0_cp*a**8*dm**2+495.0_cp &
      &              *a**8*dn**7+7425.0_cp*a**8*dn**6+3435.0_cp*a**8*dn**5-332661.0_cp*a**8 &
      &              *dn**4-449754.0_cp*a**8*dn**3+6507084.0_cp*a**8*dn**2+4628376.0_cp* &
      &              a**8*dn-48008160.0_cp*a**8-896.0_cp*a**6*b**2*dm**4*dn**3+5376.0_cp &
      &              *a**6*b**2*dm**4*dn**2+17024.0_cp*a**6*b**2*dm**4*dn-88704.0_cp*a** &
      &              6*b**2*dm**4-2880.0_cp*a**6*b**2*dm**2*dn**5+64512.0_cp*a**6*b**2*dm**2 &
      &              *dn**4+352064.0_cp*a**6*b**2*dm**2*dn**3-2584704.0_cp*a**6*b**2*dm**2* &
      &              dn**2-4664576.0_cp*a**6*b**2*dm**2*dn+25818624.0_cp*a**6*b**2*dm**2 &
      &              +31680.0_cp*a**6*b**2*dn**7+475200.0_cp*a**6*b**2*dn**6+551520.0_cp &
      &              *a**6*b**2*dn**5-17516736.0_cp*a**6*b**2*dn**4-33870240.0_cp*a**6*b**2 &
      &              *dn**3+284757120.0_cp*a**6*b**2*dn**2+315123840.0_cp*a**6*b**2*dn-1841647104.0_cp &
      &              *a**6*b**2+128.0_cp*a**6*dm**4*dn**3*r**2-768.0_cp*a**6*dm**4*dn**2 &
      &              *r**2-2432.0_cp*a**6*dm**4*dn*r**2+12672.0_cp*a**6*dm**4*r**2+256.0_cp &
      &              *a**6*dm**2*dn**5*r**2-5792.0_cp*a**6*dm**2*dn**4*r**2-32128.0_cp*a**6 &
      &              *dm**2*dn**3*r**2+235328.0_cp*a**6*dm**2*dn**2*r**2+428736.0_cp*a** &
      &              6*dm**2*dn*r**2-2370240.0_cp*a**6*dm**2*r**2-1920.0_cp*a**6*dn**7*r**2 &
      &              -28608.0_cp*a**6*dn**6*r**2-32704.0_cp*a**6*dn**5*r**2+1050592.0_cp &
      &              *a**6*dn**4*r**2+2028544.0_cp*a**6*dn**3*r**2-17063296.0_cp*a**6*dn**2 &
      &              *r**2-18911232.0_cp*a**6*dn*r**2+110490624.0_cp*a**6*r**2-4480.0_cp &
      &              *a**4*b**4*dm**4*dn**3+13440.0_cp*a**4*b**4*dm**4*dn**2+85120.0_cp* &
      &              a**4*b**4*dm**4*dn-255360.0_cp*a**4*b**4*dm**4-26880.0_cp*a**4*b**4 &
      &              *dm**2*dn**5+259392.0_cp*a**4*b**4*dm**2*dn**4+2205952.0_cp*a**4*b**4* &
      &              dm**2*dn**3-8941632.0_cp*a**4*b**4*dm**2*dn**2-26322688.0_cp*a**4*b**4 &
      &              *dm**2*dn+85403136.0_cp*a**4*b**4*dm**2+221760.0_cp*a**4*b**4*dn**7 &
      &              +3326400.0_cp*a**4*b**4*dn**6+6182400.0_cp*a**4*b**4*dn**5-96202176.0_cp &
      &              *a**4*b**4*dn**4-251662656.0_cp*a**4*b**4*dn**3+1233609216.0_cp*a** &
      &              4*b**4*dn**2+1992824064.0_cp*a**4*b**4*dn-6905945088.0_cp*a**4*b**4 &
      &              +3840.0_cp*a**4*b**2*dm**4*dn**3*r**2-11520.0_cp*a**4*b**2*dm**4*dn**2 &
      &              *r**2-72960.0_cp*a**4*b**2*dm**4*dn*r**2+218880.0_cp*a**4*b**2*dm** &
      &              4*r**2+14336.0_cp*a**4*b**2*dm**2*dn**5*r**2-139552.0_cp*a**4*b**2*dm**2 &
      &              *dn**4*r**2-1203584.0_cp*a**4*b**2*dm**2*dn**3*r**2+4865824.0_cp*a**4* &
      &              b**2*dm**2*dn**2*r**2+14480640.0_cp*a**4*b**2*dm**2*dn*r**2-46917504.0_cp &
      &              *a**4*b**2*dm**2*r**2-80640.0_cp*a**4*b**2*dn**7*r**2-1201536.0_cp* &
      &              a**4*b**2*dn**6*r**2-2207744.0_cp*a**4*b**2*dn**5*r**2+34644064.0_cp &
      &              *a**4*b**2*dn**4*r**2+90399488.0_cp*a**4*b**2*dn**3*r**2-443489632.0_cp &
      &              *a**4*b**2*dn**2*r**2-717072768.0_cp*a**4*b**2*dn*r**2+2484430848.0_cp &
      &              *a**4*b**2*r**2-384.0_cp*a**4*dm**4*dn**3*r**4+1152.0_cp*a**4*dm**4 &
      &              *dn**2*r**4+7296.0_cp*a**4*dm**4*dn*r**4-21888.0_cp*a**4*dm**4*r**4 &
      &              -768.0_cp*a**4*dm**2*dn**5*r**4+7584.0_cp*a**4*dm**2*dn**4*r**4+67072.0_cp &
      &              *a**4*dm**2*dn**3*r**4-269856.0_cp*a**4*dm**2*dn**2*r**4-818560.0_cp &
      &              *a**4*dm**2*dn*r**4+2645760.0_cp*a**4*dm**2*r**4+2688.0_cp*a**4*dn**7* &
      &              r**4+39648.0_cp*a**4*dn**6*r**4+71424.0_cp*a**4*dn**5*r**4-1135840.0_cp &
      &              *a**4*dn**4*r**4-2932224.0_cp*a**4*dn**3*r**4+14431360.0_cp*a**4*dn**2 &
      &              *r**4+23162880.0_cp*a**4*dn*r**4-80335872.0_cp*a**4*r**4-3584.0_cp* &
      &              a**2*b**6*dm**4*dn**3+68096.0_cp*a**2*b**6*dm**4*dn-107520.0_cp*a** &
      &              2*b**6*dm**4-53760.0_cp*a**2*b**6*dm**2*dn**5+107520.0_cp*a**2*b**6 &
      &              *dm**2*dn**4+2713088.0_cp*a**2*b**6*dm**2*dn**3-4042752.0_cp*a**2*b**6 &
      &              *dm**2*dn**2-26980352.0_cp*a**2*b**6*dm**2*dn+48427008.0_cp*a**2*b**6* &
      &              dm**2+354816.0_cp*a**2*b**6*dn**7+5322240.0_cp*a**2*b**6*dn**6+13606656.0_cp &
      &              *a**2*b**6*dn**5-111659520.0_cp*a**2*b**6*dn**4-392324352.0_cp*a**2 &
      &              *b**6*dn**3+1017762816.0_cp*a**2*b**6*dn**2+2482357248.0_cp*a**2*b**6* &
      &              dn-4932071424.0_cp*a**2*b**6+7680.0_cp*a**2*b**4*dm**4*dn**3*r**2-145920.0_cp &
      &              *a**2*b**4*dm**4*dn*r**2+230400.0_cp*a**2*b**4*dm**4*r**2+71680.0_cp &
      &              *a**2*b**4*dm**2*dn**5*r**2-143360.0_cp*a**2*b**4*dm**2*dn**4*r**2-3681280.0_cp &
      &              *a**2*b**4*dm**2*dn**3*r**2+5438720.0_cp*a**2*b**4*dm**2*dn**2*r**2 &
      &              +36944640.0_cp*a**2*b**4*dm**2*dn*r**2-66193920.0_cp*a**2*b**4*dm** &
      &              2*r**2-322560.0_cp*a**2*b**4*dn**7*r**2-4806144.0_cp*a**2*b**4*dn** &
      &              6*r**2-12167680.0_cp*a**2*b**4*dn**5*r**2+100653056.0_cp*a**2*b**4*dn**4 &
      &              *r**2+352180480.0_cp*a**2*b**4*dn**3*r**2-915036416.0_cp*a**2*b**4*dn**2 &
      &              *r**2-2231078400.0_cp*a**2*b**4*dn*r**2+4432656384.0_cp*a**2*b**4*r**2 &
      &              -4608.0_cp*a**2*b**2*dm**4*dn**3*r**4+87552.0_cp*a**2*b**2*dm**4*dn &
      &              *r**4-138240.0_cp*a**2*b**2*dm**4*r**4-23040.0_cp*a**2*b**2*dm**2*dn**5 &
      &              *r**4+46080.0_cp*a**2*b**2*dm**2*dn**4*r**4+1219584.0_cp*a**2*b**2*dm**2 &
      &              *dn**3*r**4-1774080.0_cp*a**2*b**2*dm**2*dn**2*r**4-12435456.0_cp*a**2 &
      &              *b**2*dm**2*dn*r**4+22210560.0_cp*a**2*b**2*dm**2*r**4+64512.0_cp*a**2 &
      &              *b**2*dn**7*r**4+951552.0_cp*a**2*b**2*dn**6*r**4+2368512.0_cp*a**2 &
      &              *b**2*dn**5*r**4-19865088.0_cp*a**2*b**2*dn**4*r**4-68680704.0_cp*a**2 &
      &              *b**2*dn**3*r**4+179467008.0_cp*a**2*b**2*dn**2*r**4+433460736.0_cp &
      &              *a**2*b**2*dn*r**4-862451712.0_cp*a**2*b**2*r**4+512.0_cp*a**2*dm** &
      &              4*dn**3*r**6-9728.0_cp*a**2*dm**4*dn*r**6+15360.0_cp*a**2*dm**4*r** &
      &              6+1024.0_cp*a**2*dm**2*dn**5*r**6-2048.0_cp*a**2*dm**2*dn**4*r**6-57856.0_cp &
      &              *a**2*dm**2*dn**3*r**6+81152.0_cp*a**2*dm**2*dn**2*r**6+610560.0_cp &
      &              *a**2*dm**2*dn*r**6-1082880.0_cp*a**2*dm**2*r**6-1536.0_cp*a**2*dn**7* &
      &              r**6-22272.0_cp*a**2*dn**6*r**6-53504.0_cp*a**2*dn**5*r**6+461824.0_cp &
      &              *a**2*dn**4*r**6+1544192.0_cp*a**2*dn**3*r**6-4102144.0_cp*a**2*dn**2* &
      &              r**6-9584640.0_cp*a**2*dn*r**6+19169280.0_cp*a**2*r**6+256.0_cp*b** &
      &              8*dm**4*dn**3-1536.0_cp*b**8*dm**4*dn**2+2816.0_cp*b**8*dm**4*dn-1536.0_cp &
      &              *b**8*dm**4-23040.0_cp*b**8*dm**2*dn**5-71424.0_cp*b**8*dm**2*dn**4 &
      &              +533504.0_cp*b**8*dm**2*dn**3+658176.0_cp*b**8*dm**2*dn**2-3923456.0_cp &
      &              *b**8*dm**2*dn+2826240.0_cp*b**8*dm**2+126720.0_cp*b**8*dn**7+1900800.0_cp &
      &              *b**8*dn**6+6186240.0_cp*b**8*dn**5-24784128.0_cp*b**8*dn**4-124406784.0_cp &
      &              *b**8*dn**3+114729984.0_cp*b**8*dn**2+570802176.0_cp*b**8*dn-544555008.0_cp &
      &              *b**8-1024.0_cp*b**6*dm**4*dn**3*r**2+6144.0_cp*b**6*dm**4*dn**2*r**2- &
      &              11264.0_cp*b**6*dm**4*dn*r**2+6144.0_cp*b**6*dm**4*r**2+57344.0_cp* &
      &              b**6*dm**2*dn**5*r**2+180992.0_cp*b**6*dm**2*dn**4*r**2-1332224.0_cp &
      &              *b**6*dm**2*dn**3*r**2-1692416.0_cp*b**6*dm**2*dn**2*r**2+9910272.0_cp &
      &              *b**6*dm**2*dn*r**2-7123968.0_cp*b**6*dm**2*r**2-215040.0_cp*b**6*dn**7 &
      &              *r**2-3204096.0_cp*b**6*dn**6*r**2-10336256.0_cp*b**6*dn**5*r**2+41819904.0_cp &
      &              *b**6*dn**4*r**2+208361984.0_cp*b**6*dn**3*r**2-193244928.0_cp*b**6 &
      &              *dn**2*r**2-956505600.0_cp*b**6*dn*r**2+913324032.0_cp*b**6*r**2+1536.0_cp &
      &              *b**4*dm**4*dn**3*r**4-9216.0_cp*b**4*dm**4*dn**2*r**4+16896.0_cp*b**4 &
      &              *dm**4*dn*r**4-9216.0_cp*b**4*dm**4*r**4-46080.0_cp*b**4*dm**2*dn** &
      &              5*r**4-149760.0_cp*b**4*dm**2*dn**4*r**4+1075712.0_cp*b**4*dm**2*dn**3 &
      &              *r**4+1436928.0_cp*b**4*dm**2*dn**2*r**4-8165888.0_cp*b**4*dm**2*dn &
      &              *r**4+5849088.0_cp*b**4*dm**2*r**4+107520.0_cp*b**4*dn**7*r**4+1585920.0_cp &
      &              *b**4*dn**6*r**4+5038080.0_cp*b**4*dn**5*r**4-20783360.0_cp*b**4*dn**4 &
      &              *r**4-101876736.0_cp*b**4*dn**3*r**4+96190976.0_cp*b**4*dn**2*r**4+466566144.0_cp &
      &              *b**4*dn*r**4-446828544.0_cp*b**4*r**4-1024.0_cp*b**2*dm**4*dn**3*r**6 &
      &              +6144.0_cp*b**2*dm**4*dn**2*r**6-11264.0_cp*b**2*dm**4*dn*r**6+6144.0_cp &
      &              *b**2*dm**4*r**6+12288.0_cp*b**2*dm**2*dn**5*r**6+42240.0_cp*b**2*dm**2 &
      &              *dn**4*r**6-288768.0_cp*b**2*dm**2*dn**3*r**6-429312.0_cp*b**2*dm** &
      &              2*dn**2*r**6+2294784.0_cp*b**2*dm**2*dn*r**6-1631232.0_cp*b**2*dm** &
      &              2*r**6-15360.0_cp*b**2*dn**7*r**6-222720.0_cp*b**2*dn**6*r**6-685056.0_cp &
      &              *b**2*dn**5*r**6+2957568.0_cp*b**2*dn**4*r**6+13902336.0_cp*b**2*dn**3 &
      &              *r**6-13821696.0_cp*b**2*dn**2*r**6-62995968.0_cp*b**2*dn*r**6+60880896.0_cp &
      &              *b**2*r**6+256.0_cp*dm**4*dn**3*r**8-1536.0_cp*dm**4*dn**2*r**8+2816.0_cp &
      &              *dm**4*dn*r**8-1536.0_cp*dm**4*r**8-512.0_cp*dm**2*dn**5*r**8-2048.0_cp &
      &              *dm**2*dn**4*r**8+11776.0_cp*dm**2*dn**3*r**8+26624.0_cp*dm**2*dn** &
      &              2*r**8-115712.0_cp*dm**2*dn*r**8+79872.0_cp*dm**2*r**8+256.0_cp*dn**7* &
      &              r**8+3584.0_cp*dn**6*r**8+9984.0_cp*dn**5*r**8-49664.0_cp*dn**4*r** &
      &              8-203776.0_cp*dn**3*r**8+239616.0_cp*dn**2*r**8+884736.0_cp*dn*r**8 &
      &              -884736.0_cp*r**8)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-2)=-a**3*b*(24.0_cp*a**8*dm**4*dn**2-600.0_cp*a**8*dm**4 &
      &              *dn+2064.0_cp*a**8*dm**4-60.0_cp*a**8*dm**2*dn**4-4779.0_cp*a**8*dm**2 &
      &              *dn**3+13806.0_cp*a**8*dm**2*dn**2+162591.0_cp*a**8*dm**2*dn-570774.0_cp &
      &              *a**8*dm**2-1980.0_cp*a**8*dn**6-16335.0_cp*a**8*dn**5+97455.0_cp*a**8 &
      &              *dn**4+750987.0_cp*a**8*dn**3-2754471.0_cp*a**8*dn**2-10837188.0_cp &
      &              *a**8*dn+38725884.0_cp*a**8+448.0_cp*a**6*b**2*dm**4*dn**2-6720.0_cp &
      &              *a**6*b**2*dm**4*dn+18368.0_cp*a**6*b**2*dm**4-68544.0_cp*a**6*b**2 &
      &              *dm**2*dn**3+119168.0_cp*a**6*b**2*dm**2*dn**2+1984416.0_cp*a**6*b**2* &
      &              dm**2*dn-5607392.0_cp*a**6*b**2*dm**2-36960.0_cp*a**6*b**2*dn**6-304920.0_cp &
      &              *a**6*b**2*dn**5+1432200.0_cp*a**6*b**2*dn**4+11877432.0_cp*a**6*b**2* &
      &              dn**3-33540360.0_cp*a**6*b**2*dn**2-143995488.0_cp*a**6*b**2*dn+417705120.0_cp &
      &              *a**6*b**2-192.0_cp*a**6*dm**4*dn**2*r**2+2880.0_cp*a**6*dm**4*dn*r**2 &
      &              -7872.0_cp*a**6*dm**4*r**2+18480.0_cp*a**6*dm**2*dn**3*r**2-31392.0_cp &
      &              *a**6*dm**2*dn**2*r**2-546360.0_cp*a**6*dm**2*dn*r**2+1542888.0_cp* &
      &              a**6*dm**2*r**2+6720.0_cp*a**6*dn**6*r**2+54936.0_cp*a**6*dn**5*r** &
      &              2-258888.0_cp*a**6*dn**4*r**2-2132872.0_cp*a**6*dn**3*r**2+6022440.0_cp &
      &              *a**6*dn**2*r**2+25914472.0_cp*a**6*dn*r**2-75163512.0_cp*a**6*r**2 &
      &              +1344.0_cp*a**4*b**4*dm**4*dn**2-12096.0_cp*a**4*b**4*dm**4*dn+24192.0_cp &
      &              *a**4*b**4*dm**4+4032.0_cp*a**4*b**4*dm**2*dn**4-172368.0_cp*a**4*b**4 &
      &              *dm**2*dn**3+81312.0_cp*a**4*b**4*dm**2*dn**2+4139856.0_cp*a**4*b** &
      &              4*dm**2*dn-8824032.0_cp*a**4*b**4*dm**2-133056.0_cp*a**4*b**4*dn**6 &
      &              -1097712.0_cp*a**4*b**4*dn**5+3762864.0_cp*a**4*b**4*dn**4+35051184.0_cp &
      &              *a**4*b**4*dn**3-69008688.0_cp*a**4*b**4*dn**2-343635264.0_cp*a**4* &
      &              b**4*dn+764551872.0_cp*a**4*b**4-1920.0_cp*a**4*b**2*dm**4*dn**2*r**2+ &
      &              17280.0_cp*a**4*b**2*dm**4*dn*r**2-34560.0_cp*a**4*b**2*dm**4*r**2-3584.0_cp &
      &              *a**4*b**2*dm**2*dn**4*r**2+154728.0_cp*a**4*b**2*dm**2*dn**3*r**2-65392.0_cp &
      &              *a**4*b**2*dm**2*dn**2*r**2-3788712.0_cp*a**4*b**2*dm**2*dn*r**2+8067312.0_cp &
      &              *a**4*b**2*dm**2*r**2+80640.0_cp*a**4*b**2*dn**6*r**2+659232.0_cp*a**4 &
      &              *b**2*dn**5*r**2-2272480.0_cp*a**4*b**2*dn**4*r**2-20986392.0_cp*a**4* &
      &              b**2*dn**3*r**2+41331568.0_cp*a**4*b**2*dn**2*r**2+206009976.0_cp*a**4 &
      &              *b**2*dn*r**2-458290512.0_cp*a**4*b**2*r**2+576.0_cp*a**4*dm**4*dn**2* &
      &              r**4-5184.0_cp*a**4*dm**4*dn*r**4+10368.0_cp*a**4*dm**4*r**4+576.0_cp &
      &              *a**4*dm**2*dn**4*r**4-25272.0_cp*a**4*dm**2*dn**3*r**4+8400.0_cp*a**4 &
      &              *dm**2*dn**2*r**4+640440.0_cp*a**4*dm**2*dn*r**4-1361232.0_cp*a**4*dm**2 &
      &              *r**4-8064.0_cp*a**4*dn**6*r**4-65016.0_cp*a**4*dn**5*r**4+226152.0_cp &
      &              *a**4*dn**4*r**4+2053104.0_cp*a**4*dn**3*r**4-4066968.0_cp*a**4*dn**2* &
      &              r**4-19989816.0_cp*a**4*dn*r**4+44497584.0_cp*a**4*r**4+768.0_cp*a**2* &
      &              b**6*dm**4*dn**2-3840.0_cp*a**2*b**6*dm**4*dn+4608.0_cp*a**2*b**6*dm**4 &
      &              +7680.0_cp*a**2*b**6*dm**2*dn**4-93312.0_cp*a**2*b**6*dm**2*dn**3-109824.0_cp &
      &              *a**2*b**6*dm**2*dn**2+1822848.0_cp*a**2*b**6*dm**2*dn-2582784.0_cp &
      &              *a**2*b**6*dm**2-126720.0_cp*a**2*b**6*dn**6-1045440.0_cp*a**2*b**6 &
      &              *dn**5+2256960.0_cp*a**2*b**6*dn**4+26041536.0_cp*a**2*b**6*dn**3-28467264.0_cp &
      &              *a**2*b**6*dn**2-194298624.0_cp*a**2*b**6*dn+299586816.0_cp*a**2*b**6- &
      &              2304.0_cp*a**2*b**4*dm**4*dn**2*r**2+11520.0_cp*a**2*b**4*dm**4*dn* &
      &              r**2-13824.0_cp*a**2*b**4*dm**4*r**2-14336.0_cp*a**2*b**4*dm**2*dn**4* &
      &              r**2+175392.0_cp*a**2*b**4*dm**2*dn**3*r**2+216896.0_cp*a**2*b**4*dm**2 &
      &              *dn**2*r**2-3485088.0_cp*a**2*b**4*dm**2*dn*r**2+4928832.0_cp*a**2* &
      &              b**4*dm**2*r**2+161280.0_cp*a**2*b**4*dn**6*r**2+1318464.0_cp*a**2* &
      &              b**4*dn**5*r**2-2876608.0_cp*a**2*b**4*dn**4*r**2-32756640.0_cp*a** &
      &              2*b**4*dn**3*r**2+35897728.0_cp*a**2*b**4*dn**2*r**2+244419936.0_cp &
      &              *a**2*b**4*dn*r**2-376864704.0_cp*a**2*b**4*r**2+2304.0_cp*a**2*b** &
      &              2*dm**4*dn**2*r**4-11520.0_cp*a**2*b**2*dm**4*dn*r**4+13824.0_cp*a**2* &
      &              b**2*dm**4*r**4+7680.0_cp*a**2*b**2*dm**2*dn**4*r**4-95040.0_cp*a** &
      &              2*b**2*dm**2*dn**3*r**4-127872.0_cp*a**2*b**2*dm**2*dn**2*r**4+1945920.0_cp &
      &              *a**2*b**2*dm**2*dn*r**4-2742912.0_cp*a**2*b**2*dm**2*r**4-53760.0_cp &
      &              *a**2*b**2*dn**6*r**4-433440.0_cp*a**2*b**2*dn**5*r**4+962400.0_cp* &
      &              a**2*b**2*dn**4*r**4+10701280.0_cp*a**2*b**2*dn**3*r**4-11877984.0_cp &
      &              *a**2*b**2*dn**2*r**4-79284160.0_cp*a**2*b**2*dn*r**4+122382336.0_cp &
      &              *a**2*b**2*r**4-768.0_cp*a**2*dm**4*dn**2*r**6+3840.0_cp*a**2*dm**4 &
      &              *dn*r**6-4608.0_cp*a**2*dm**4*r**6-1024.0_cp*a**2*dm**2*dn**4*r**6+12960.0_cp &
      &              *a**2*dm**2*dn**3*r**6+20800.0_cp*a**2*dm**2*dn**2*r**6-283680.0_cp &
      &              *a**2*dm**2*dn*r**6+396864.0_cp*a**2*dm**2*r**6+3840.0_cp*a**2*dn** &
      &              6*r**6+30240.0_cp*a**2*dn**5*r**6-69472.0_cp*a**2*dn**4*r**6-733824.0_cp &
      &              *a**2*dn**3*r**6+845728.0_cp*a**2*dn**2*r**6+5295456.0_cp*a**2*dn*r**6 &
      &              -8205120.0_cp*a**2*r**6+2560.0_cp*b**8*dm**2*dn**4-4992.0_cp*b**8*dm**2 &
      &              *dn**3-34048.0_cp*b**8*dm**2*dn**2+98688.0_cp*b**8*dm**2*dn-62208.0_cp &
      &              *b**8*dm**2-28160.0_cp*b**8*dn**6-232320.0_cp*b**8*dn**5+206720.0_cp &
      &              *b**8*dn**4+4155776.0_cp*b**8*dn**3-717696.0_cp*b**8*dn**2-21062144.0_cp &
      &              *b**8*dn+17677824.0_cp*b**8-8192.0_cp*b**6*dm**2*dn**4*r**2+15744.0_cp &
      &              *b**6*dm**2*dn**3*r**2+110336.0_cp*b**6*dm**2*dn**2*r**2-318336.0_cp &
      &              *b**6*dm**2*dn*r**2+200448.0_cp*b**6*dm**2*r**2+61440.0_cp*b**6*dn**6* &
      &              r**2+502272.0_cp*b**6*dn**5*r**2-460288.0_cp*b**6*dn**4*r**2-8967808.0_cp &
      &              *b**6*dn**3*r**2+1616128.0_cp*b**6*dn**2*r**2+45389440.0_cp*b**6*dn &
      &              *r**2-38141184.0_cp*b**6*r**2+9216.0_cp*b**4*dm**2*dn**4*r**4-17280.0_cp &
      &              *b**4*dm**2*dn**3*r**4-126720.0_cp*b**4*dm**2*dn**2*r**4+362880.0_cp &
      &              *b**4*dm**2*dn*r**4-228096.0_cp*b**4*dm**2*r**4-43008.0_cp*b**4*dn**6* &
      &              r**4-346752.0_cp*b**4*dn**5*r**4+333696.0_cp*b**4*dn**4*r**4+6172160.0_cp &
      &              *b**4*dn**3*r**4-1222272.0_cp*b**4*dn**2*r**4-31084928.0_cp*b**4*dn &
      &              *r**4+26191104.0_cp*b**4*r**4-4096.0_cp*b**2*dm**2*dn**4*r**6+7296.0_cp &
      &              *b**2*dm**2*dn**3*r**6+58624.0_cp*b**2*dm**2*dn**2*r**6-165504.0_cp &
      &              *b**2*dm**2*dn*r**6+103680.0_cp*b**2*dm**2*r**6+10240.0_cp*b**2*dn**6* &
      &              r**6+80640.0_cp*b**2*dn**5*r**6-85248.0_cp*b**2*dn**4*r**6-1426560.0_cp &
      &              *b**2*dn**3*r**6+350720.0_cp*b**2*dn**2*r**6+7067520.0_cp*b**2*dn*r**6 &
      &              -5997312.0_cp*b**2*r**6+512.0_cp*dm**2*dn**4*r**8-768.0_cp*dm**2*dn**3 &
      &              *r**8-8192.0_cp*dm**2*dn**2*r**8+22272.0_cp*dm**2*dn*r**8-13824.0_cp &
      &              *dm**2*r**8-512.0_cp*dn**6*r**8-3840.0_cp*dn**5*r**8+5120.0_cp*dn** &
      &              4*r**8+66432.0_cp*dn**3*r**8-26880.0_cp*dn**2*r**8-309888.0_cp*dn*r**8 &
      &              +269568.0_cp*r**8)/(1024.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku-1)=-a**2*(-2.0_cp*a**10*dm**4*dn**2-50.0_cp*a**10*dm**4*dn &
      &              +288.0_cp*a**10*dm**4-24.0_cp*a**10*dm**2*dn**4-387.0_cp*a**10*dm** &
      &              2*dn**3+2948.0_cp*a**10*dm**2*dn**2+13133.0_cp*a**10*dm**2*dn-76698.0_cp &
      &              *a**10*dm**2-198.0_cp*a**10*dn**6-1089.0_cp*a**10*dn**5+13242.0_cp* &
      &              a**10*dn**4+55341.0_cp*a**10*dn**3-398532.0_cp*a**10*dn**2-845460.0_cp &
      &              *a**10*dn+5010768.0_cp*a**10-56.0_cp*a**8*b**2*dm**4*dn**2-2240.0_cp &
      &              *a**8*b**2*dm**4*dn+10584.0_cp*a**8*b**2*dm**4-1260.0_cp*a**8*b**2*dm**2 &
      &              *dn**4-21924.0_cp*a**8*b**2*dm**2*dn**3+133028.0_cp*a**8*b**2*dm**2 &
      &              *dn**2+629636.0_cp*a**8*b**2*dm**2*dn-3040632.0_cp*a**8*b**2*dm**2-13860.0_cp &
      &              *a**8*b**2*dn**6-76230.0_cp*a**8*b**2*dn**5+781830.0_cp*a**8*b**2*dn**4 &
      &              +3338622.0_cp*a**8*b**2*dn**3-19926522.0_cp*a**8*b**2*dn**2-43387848.0_cp &
      &              *a**8*b**2*dn+213931368.0_cp*a**8*b**2+8.0_cp*a**8*dm**4*dn**2*r**2 &
      &              +320.0_cp*a**8*dm**4*dn*r**2-1512.0_cp*a**8*dm**4*r**2+112.0_cp*a** &
      &              8*dm**2*dn**4*r**2+1974.0_cp*a**8*dm**2*dn**3*r**2-11992.0_cp*a**8*dm**2 &
      &              *dn**2*r**2-57846.0_cp*a**8*dm**2*dn*r**2+279216.0_cp*a**8*dm**2*r**2+ &
      &              840.0_cp*a**8*dn**6*r**2+4578.0_cp*a**8*dn**5*r**2-46984.0_cp*a**8*dn**4 &
      &              *r**2-199850.0_cp*a**8*dn**3*r**2+1193144.0_cp*a**8*dn**2*r**2+2603592.0_cp &
      &              *a**8*dn*r**2-12835872.0_cp*a**8*r**2-10080.0_cp*a**6*b**4*dm**4*dn &
      &              +36960.0_cp*a**6*b**4*dm**4-6720.0_cp*a**6*b**4*dm**2*dn**4-134064.0_cp &
      &              *a**6*b**4*dm**2*dn**3+603456.0_cp*a**6*b**4*dm**2*dn**2+3152016.0_cp &
      &              *a**6*b**4*dm**2*dn-12020064.0_cp*a**6*b**4*dm**2-110880.0_cp*a**6* &
      &              b**4*dn**6-609840.0_cp*a**6*b**4*dn**5+5093760.0_cp*a**6*b**4*dn**4 &
      &              +22426992.0_cp*a**6*b**4*dn**3-106161888.0_cp*a**6*b**4*dn**2-240262848.0_cp &
      &              *a**6*b**4*dn+945689472.0_cp*a**6*b**4+8640.0_cp*a**6*b**2*dm**4*dn &
      &              *r**2-31680.0_cp*a**6*b**2*dm**4*r**2+3584.0_cp*a**6*b**2*dm**2*dn**4* &
      &              r**2+72408.0_cp*a**6*b**2*dm**2*dn**3*r**2-325472.0_cp*a**6*b**2*dm**2 &
      &              *dn**2*r**2-1733832.0_cp*a**6*b**2*dm**2*dn*r**2+6606864.0_cp*a**6* &
      &              b**2*dm**2*r**2+40320.0_cp*a**6*b**2*dn**6*r**2+219744.0_cp*a**6*b**2* &
      &              dn**5*r**2-1838144.0_cp*a**6*b**2*dn**4*r**2-8056776.0_cp*a**6*b**2 &
      &              *dn**3*r**2+38149664.0_cp*a**6*b**2*dn**2*r**2+86457240.0_cp*a**6*b**2 &
      &              *dn*r**2-340250832.0_cp*a**6*b**2*r**2-864.0_cp*a**6*dm**4*dn*r**4+3168.0_cp &
      &              *a**6*dm**4*r**4-192.0_cp*a**6*dm**2*dn**4*r**4-3960.0_cp*a**6*dm** &
      &              2*dn**3*r**4+17760.0_cp*a**6*dm**2*dn**2*r**4+97992.0_cp*a**6*dm**2 &
      &              *dn*r**4-372912.0_cp*a**6*dm**2*r**4-1344.0_cp*a**6*dn**6*r**4-7224.0_cp &
      &              *a**6*dn**5*r**4+60528.0_cp*a**6*dn**4*r**4+262440.0_cp*a**6*dn**3* &
      &              r**4-1243776.0_cp*a**6*dn**2*r**4-2793312.0_cp*a**6*dn*r**4+10998144.0_cp &
      &              *a**6*r**4+448.0_cp*a**4*b**6*dm**4*dn**2-8960.0_cp*a**4*b**6*dm**4 &
      &              *dn+22848.0_cp*a**4*b**6*dm**4-6720.0_cp*a**4*b**6*dm**2*dn**4-185472.0_cp &
      &              *a**4*b**6*dm**2*dn**3+537152.0_cp*a**4*b**6*dm**2*dn**2+3398528.0_cp &
      &              *a**4*b**6*dm**2*dn-9477888.0_cp*a**4*b**6*dm**2-221760.0_cp*a**4*b**6 &
      &              *dn**6-1219680.0_cp*a**4*b**6*dn**5+7865760.0_cp*a**4*b**6*dn**4+36290016.0_cp &
      &              *a**4*b**6*dn**3-126854112.0_cp*a**4*b**6*dn**2-305875584.0_cp*a**4 &
      &              *b**6*dn+900402048.0_cp*a**4*b**6-960.0_cp*a**4*b**4*dm**4*dn**2*r**2+ &
      &              19200.0_cp*a**4*b**4*dm**4*dn*r**2-48960.0_cp*a**4*b**4*dm**4*r**2+8960.0_cp &
      &              *a**4*b**4*dm**2*dn**4*r**2+250320.0_cp*a**4*b**4*dm**2*dn**3*r**2-720320.0_cp &
      &              *a**4*b**4*dm**2*dn**2*r**2-4657680.0_cp*a**4*b**4*dm**2*dn*r**2+12971520.0_cp &
      &              *a**4*b**4*dm**2*r**2+201600.0_cp*a**4*b**4*dn**6*r**2+1098720.0_cp &
      &              *a**4*b**4*dn**5*r**2-7105280.0_cp*a**4*b**4*dn**4*r**2-32603760.0_cp &
      &              *a**4*b**4*dn**3*r**2+114029600.0_cp*a**4*b**4*dn**2*r**2+274983120.0_cp &
      &              *a**4*b**4*dn*r**2-809341920.0_cp*a**4*b**4*r**2+576.0_cp*a**4*b**2 &
      &              *dm**4*dn**2*r**4-11520.0_cp*a**4*b**2*dm**4*dn*r**4+29376.0_cp*a** &
      &              4*b**2*dm**4*r**4-2880.0_cp*a**4*b**2*dm**2*dn**4*r**4-82080.0_cp*a**4 &
      &              *b**2*dm**2*dn**3*r**4+233472.0_cp*a**4*b**2*dm**2*dn**2*r**4+1570080.0_cp &
      &              *a**4*b**2*dm**2*dn*r**4-4362048.0_cp*a**4*b**2*dm**2*r**4-40320.0_cp &
      &              *a**4*b**2*dn**6*r**4-216720.0_cp*a**4*b**2*dn**5*r**4+1406880.0_cp &
      &              *a**4*b**2*dn**4*r**4+6380160.0_cp*a**4*b**2*dn**3*r**4-22360896.0_cp &
      &              *a**4*b**2*dn**2*r**4-53420400.0_cp*a**4*b**2*dn*r**4+157343904.0_cp &
      &              *a**4*b**2*r**4-64.0_cp*a**4*dm**4*dn**2*r**6+1280.0_cp*a**4*dm**4*dn &
      &              *r**6-3264.0_cp*a**4*dm**4*r**6+128.0_cp*a**4*dm**2*dn**4*r**6+3792.0_cp &
      &              *a**4*dm**2*dn**3*r**6-10496.0_cp*a**4*dm**2*dn**2*r**6-77328.0_cp* &
      &              a**4*dm**2*dn*r**6+213696.0_cp*a**4*dm**2*r**6+960.0_cp*a**4*dn**6* &
      &              r**6+5040.0_cp*a**4*dn**5*r**6-32896.0_cp*a**4*dn**4*r**6-144816.0_cp &
      &              *a**4*dn**3*r**6+510592.0_cp*a**4*dn**2*r**6+1180224.0_cp*a**4*dn*r**6 &
      &              -3485952.0_cp*a**4*r**6+256.0_cp*a**2*b**8*dm**4*dn**2-1280.0_cp*a**2* &
      &              b**8*dm**4*dn+1536.0_cp*a**2*b**8*dm**4-58752.0_cp*a**2*b**8*dm**2*dn**3 &
      &              +68096.0_cp*a**2*b**8*dm**2*dn**2+775808.0_cp*a**2*b**8*dm**2*dn-1353984.0_cp &
      &              *a**2*b**8*dm**2-126720.0_cp*a**2*b**8*dn**6-696960.0_cp*a**2*b**8*dn**5 &
      &              +3168000.0_cp*a**2*b**8*dn**4+15843456.0_cp*a**2*b**8*dn**3-35665920.0_cp &
      &              *a**2*b**8*dn**2-97288704.0_cp*a**2*b**8*dn+190218240.0_cp*a**2*b** &
      &              8-1024.0_cp*a**2*b**6*dm**4*dn**2*r**2+5120.0_cp*a**2*b**6*dm**4*dn &
      &              *r**2-6144.0_cp*a**2*b**6*dm**4*r**2+147840.0_cp*a**2*b**6*dm**2*dn**3 &
      &              *r**2-167424.0_cp*a**2*b**6*dm**2*dn**2*r**2-1971840.0_cp*a**2*b**6 &
      &              *dm**2*dn*r**2+3430656.0_cp*a**2*b**6*dm**2*r**2+215040.0_cp*a**2*b**6 &
      &              *dn**6*r**2+1171968.0_cp*a**2*b**6*dn**5*r**2-5354496.0_cp*a**2*b** &
      &              6*dn**4*r**2-26585216.0_cp*a**2*b**6*dn**3*r**2+59944960.0_cp*a**2* &
      &              b**6*dn**2*r**2+163149696.0_cp*a**2*b**6*dn*r**2-318991104.0_cp*a** &
      &              2*b**6*r**2+1536.0_cp*a**2*b**4*dm**4*dn**2*r**4-7680.0_cp*a**2*b** &
      &              4*dm**4*dn*r**4+9216.0_cp*a**2*b**4*dm**4*r**4-120960.0_cp*a**2*b** &
      &              4*dm**2*dn**3*r**4+131072.0_cp*a**2*b**4*dm**2*dn**2*r**4+1642880.0_cp &
      &              *a**2*b**4*dm**2*dn*r**4-2842368.0_cp*a**2*b**4*dm**2*r**4-107520.0_cp &
      &              *a**2*b**4*dn**6*r**4-577920.0_cp*a**2*b**4*dn**5*r**4+2661120.0_cp &
      &              *a**2*b**4*dn**4*r**4+13032320.0_cp*a**2*b**4*dn**3*r**4-29526016.0_cp &
      &              *a**2*b**4*dn**2*r**4-79514880.0_cp*a**2*b**4*dn*r**4+155672064.0_cp &
      &              *a**2*b**4*r**4-1024.0_cp*a**2*b**2*dm**4*dn**2*r**6+5120.0_cp*a**2 &
      &              *b**2*dm**4*dn*r**6-6144.0_cp*a**2*b**2*dm**4*r**6+33408.0_cp*a**2* &
      &              b**2*dm**2*dn**3*r**6-32256.0_cp*a**2*b**2*dm**2*dn**2*r**6-473472.0_cp &
      &              *a**2*b**2*dm**2*dn*r**6+808704.0_cp*a**2*b**2*dm**2*r**6+15360.0_cp &
      &              *a**2*b**2*dn**6*r**6+80640.0_cp*a**2*b**2*dn**5*r**6-376320.0_cp*a**2 &
      &              *b**2*dn**4*r**6-1786752.0_cp*a**2*b**2*dn**3*r**6+4098048.0_cp*a** &
      &              2*b**2*dn**2*r**6+10689408.0_cp*a**2*b**2*dn*r**6-21019392.0_cp*a** &
      &              2*b**2*r**6+256.0_cp*a**2*dm**4*dn**2*r**8-1280.0_cp*a**2*dm**4*dn* &
      &              r**8+1536.0_cp*a**2*dm**4*r**8-1536.0_cp*a**2*dm**2*dn**3*r**8+512.0_cp &
      &              *a**2*dm**2*dn**2*r**8+26624.0_cp*a**2*dm**2*dn*r**8-43008.0_cp*a** &
      &              2*dm**2*r**8-256.0_cp*a**2*dn**6*r**8-1280.0_cp*a**2*dn**5*r**8+6144.0_cp &
      &              *a**2*dn**4*r**8+26624.0_cp*a**2*dn**3*r**8-63488.0_cp*a**2*dn**2*r**8 &
      &              -147456.0_cp*a**2*dn*r**8+294912.0_cp*a**2*r**8+512.0_cp*b**10*dm** &
      &              2*dn**4-1536.0_cp*b**10*dm**2*dn**3-3584.0_cp*b**10*dm**2*dn**2+13824.0_cp &
      &              *b**10*dm**2*dn-9216.0_cp*b**10*dm**2-16896.0_cp*b**10*dn**6-92928.0_cp &
      &              *b**10*dn**5+245504.0_cp*b**10*dn**4+1459968.0_cp*b**10*dn**3-1448192.0_cp &
      &              *b**10*dn**2-5612544.0_cp*b**10*dn+5465088.0_cp*b**10-2048.0_cp*b** &
      &              8*dm**2*dn**4*r**2+6144.0_cp*b**8*dm**2*dn**3*r**2+14336.0_cp*b**8*dm**2 &
      &              *dn**2*r**2-55296.0_cp*b**8*dm**2*dn*r**2+36864.0_cp*b**8*dm**2*r** &
      &              2+46080.0_cp*b**8*dn**6*r**2+251136.0_cp*b**8*dn**5*r**2-670720.0_cp &
      &              *b**8*dn**4*r**2-3941376.0_cp*b**8*dn**3*r**2+3943936.0_cp*b**8*dn**2* &
      &              r**2+15130368.0_cp*b**8*dn*r**2-14759424.0_cp*b**8*r**2+3072.0_cp*b**6 &
      &              *dm**2*dn**4*r**4-9216.0_cp*b**6*dm**2*dn**3*r**4-21504.0_cp*b**6*dm**2 &
      &              *dn**2*r**4+82944.0_cp*b**6*dm**2*dn*r**4-55296.0_cp*b**6*dm**2*r** &
      &              4-43008.0_cp*b**6*dn**6*r**4-231168.0_cp*b**6*dn**5*r**4+628224.0_cp &
      &              *b**6*dn**4*r**4+3620352.0_cp*b**6*dn**3*r**4-3677184.0_cp*b**6*dn**2* &
      &              r**4-13858560.0_cp*b**6*dn*r**4+13561344.0_cp*b**6*r**4-2048.0_cp*b**4 &
      &              *dm**2*dn**4*r**6+6144.0_cp*b**4*dm**2*dn**3*r**6+14336.0_cp*b**4*dm**2 &
      &              *dn**2*r**6-55296.0_cp*b**4*dm**2*dn*r**6+36864.0_cp*b**4*dm**2*r** &
      &              6+15360.0_cp*b**4*dn**6*r**6+80640.0_cp*b**4*dn**5*r**6-226304.0_cp &
      &              *b**4*dn**4*r**6-1256448.0_cp*b**4*dn**3*r**6+1315328.0_cp*b**4*dn**2* &
      &              r**6+4776192.0_cp*b**4*dn*r**6-4704768.0_cp*b**4*r**6+512.0_cp*b**2 &
      &              *dm**2*dn**4*r**8-1536.0_cp*b**2*dm**2*dn**3*r**8-3584.0_cp*b**2*dm**2 &
      &              *dn**2*r**8+13824.0_cp*b**2*dm**2*dn*r**8-9216.0_cp*b**2*dm**2*r**8 &
      &              -1536.0_cp*b**2*dn**6*r**8-7680.0_cp*b**2*dn**5*r**8+23296.0_cp*b** &
      &              2*dn**4*r**8+117504.0_cp*b**2*dn**3*r**8-133888.0_cp*b**2*dn**2*r** &
      &              8-435456.0_cp*b**2*dn*r**8+437760.0_cp*b**2*r**8)/(1024.0_cp*dn*(dn &
      &              -3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+3.0_cp))
      stencil(ku)=a*b*(56.0_cp*a**10*dm**4*dn**2+280.0_cp*a**10*dm**4*dn-3024.0_cp &
      &              *a**10*dm**4+420.0_cp*a**10*dm**2*dn**4+2583.0_cp*a**10*dm**2*dn**3 &
      &              -42098.0_cp*a**10*dm**2*dn**2-77287.0_cp*a**10*dm**2*dn+831222.0_cp &
      &              *a**10*dm**2+2772.0_cp*a**10*dn**6+7623.0_cp*a**10*dn**5-191541.0_cp &
      &              *a**10*dn**4-366723.0_cp*a**10*dn**3+5285469.0_cp*a**10*dn**2+5122404.0_cp &
      &              *a**10*dn-56029428.0_cp*a**10+672.0_cp*a**8*b**2*dm**4*dn**2+3360.0_cp &
      &              *a**8*b**2*dm**4*dn-29568.0_cp*a**8*b**2*dm**4+6720.0_cp*a**8*b**2*dm**2 &
      &              *dn**4+41328.0_cp*a**8*b**2*dm**2*dn**3-547680.0_cp*a**8*b**2*dm**2 &
      &              *dn**2-1020432.0_cp*a**8*b**2*dm**2*dn+8924160.0_cp*a**8*b**2*dm**2 &
      &              +55440.0_cp*a**8*b**2*dn**6+152460.0_cp*a**8*b**2*dn**5-3250380.0_cp &
      &              *a**8*b**2*dn**4-6263964.0_cp*a**8*b**2*dn**3+74991756.0_cp*a**8*b**2* &
      &              dn**2+73298736.0_cp*a**8*b**2*dn-656967024.0_cp*a**8*b**2-288.0_cp* &
      &              a**8*dm**4*dn**2*r**2-1440.0_cp*a**8*dm**4*dn*r**2+12672.0_cp*a**8*dm**4 &
      &              *r**2-1792.0_cp*a**8*dm**2*dn**4*r**2-11172.0_cp*a**8*dm**2*dn**3*r**2 &
      &              +148744.0_cp*a**8*dm**2*dn**2*r**2+280908.0_cp*a**8*dm**2*dn*r**2-2456064.0_cp &
      &              *a**8*dm**2*r**2-10080.0_cp*a**8*dn**6*r**2-27468.0_cp*a**8*dn**5*r**2 &
      &              +585676.0_cp*a**8*dn**4*r**2+1125096.0_cp*a**8*dn**3*r**2-13475428.0_cp &
      &              *a**8*dn**2*r**2-13191876.0_cp*a**8*dn*r**2+118223856.0_cp*a**8*r** &
      &              2+1344.0_cp*a**6*b**4*dm**4*dn**2+6720.0_cp*a**6*b**4*dm**4*dn-45696.0_cp &
      &              *a**6*b**4*dm**4+20160.0_cp*a**6*b**4*dm**2*dn**4+123984.0_cp*a**6* &
      &              b**4*dm**2*dn**3-1265376.0_cp*a**6*b**4*dm**2*dn**2-2412816.0_cp*a**6* &
      &              b**4*dm**2*dn+16240224.0_cp*a**6*b**4*dm**2+221760.0_cp*a**6*b**4*dn**6 &
      &              +609840.0_cp*a**6*b**4*dn**5-10679760.0_cp*a**6*b**4*dn**4-20773872.0_cp &
      &              *a**6*b**4*dn**3+198127440.0_cp*a**6*b**4*dn**2+196112448.0_cp*a**6 &
      &              *b**4*dn-1372472640.0_cp*a**6*b**4-1920.0_cp*a**6*b**2*dm**4*dn**2* &
      &              r**2-9600.0_cp*a**6*b**2*dm**4*dn*r**2+65280.0_cp*a**6*b**2*dm**4*r**2 &
      &              -17920.0_cp*a**6*b**2*dm**2*dn**4*r**2-111720.0_cp*a**6*b**2*dm**2*dn**3 &
      &              *r**2+1143760.0_cp*a**6*b**2*dm**2*dn**2*r**2+2207880.0_cp*a**6*b** &
      &              2*dm**2*dn*r**2-14854320.0_cp*a**6*b**2*dm**2*r**2-134400.0_cp*a**6 &
      &              *b**2*dn**6*r**2-366240.0_cp*a**6*b**2*dn**5*r**2+6418720.0_cp*a**6 &
      &              *b**2*dn**4*r**2+12441240.0_cp*a**6*b**2*dn**3*r**2-118703440.0_cp* &
      &              a**6*b**2*dn**2*r**2-117578520.0_cp*a**6*b**2*dn*r**2+822756240.0_cp &
      &              *a**6*b**2*r**2+576.0_cp*a**6*dm**4*dn**2*r**4+2880.0_cp*a**6*dm**4 &
      &              *dn*r**4-19584.0_cp*a**6*dm**4*r**4+2880.0_cp*a**6*dm**2*dn**4*r**4 &
      &              +18360.0_cp*a**6*dm**2*dn**3*r**4-189168.0_cp*a**6*dm**2*dn**2*r**4 &
      &              -373080.0_cp*a**6*dm**2*dn*r**4+2508432.0_cp*a**6*dm**2*r**4+13440.0_cp &
      &              *a**6*dn**6*r**4+36120.0_cp*a**6*dn**5*r**4-633240.0_cp*a**6*dn**4* &
      &              r**4-1216080.0_cp*a**6*dn**3*r**4+11601384.0_cp*a**6*dn**2*r**4+11407320.0_cp &
      &              *a**6*dn*r**4-79859376.0_cp*a**6*r**4+512.0_cp*a**4*b**6*dm**4*dn** &
      &              2+2560.0_cp*a**4*b**6*dm**4*dn-12288.0_cp*a**4*b**6*dm**4+15360.0_cp &
      &              *a**4*b**6*dm**2*dn**4+94464.0_cp*a**4*b**6*dm**2*dn**3-676352.0_cp &
      &              *a**4*b**6*dm**2*dn**2-1344256.0_cp*a**4*b**6*dm**2*dn+6325248.0_cp &
      &              *a**4*b**6*dm**2+253440.0_cp*a**4*b**6*dn**6+696960.0_cp*a**4*b**6*dn**5 &
      &              -9552000.0_cp*a**4*b**6*dn**4-18847872.0_cp*a**4*b**6*dn**3+134078592.0_cp &
      &              *a**4*b**6*dn**2+135479808.0_cp*a**4*b**6*dn-684661248.0_cp*a**4*b**6- &
      &              1536.0_cp*a**4*b**4*dm**4*dn**2*r**2-7680.0_cp*a**4*b**4*dm**4*dn*r**2 &
      &              +36864.0_cp*a**4*b**4*dm**4*r**2-28672.0_cp*a**4*b**4*dm**2*dn**4*r**2 &
      &              -178752.0_cp*a**4*b**4*dm**2*dn**3*r**2+1280128.0_cp*a**4*b**4*dm** &
      &              2*dn**2*r**2+2570688.0_cp*a**4*b**4*dm**2*dn*r**2-12084480.0_cp*a** &
      &              4*b**4*dm**2*r**2-322560.0_cp*a**4*b**4*dn**6*r**2-878976.0_cp*a**4 &
      &              *b**4*dn**5*r**2+12068224.0_cp*a**4*b**4*dn**4*r**2+23714880.0_cp*a**4 &
      &              *b**4*dn**3*r**2-168782848.0_cp*a**4*b**4*dn**2*r**2-170455104.0_cp &
      &              *a**4*b**4*dn*r**2+861320448.0_cp*a**4*b**4*r**2+1536.0_cp*a**4*b** &
      &              2*dm**4*dn**2*r**4+7680.0_cp*a**4*b**2*dm**4*dn*r**4-36864.0_cp*a** &
      &              4*b**2*dm**4*r**4+15360.0_cp*a**4*b**2*dm**2*dn**4*r**4+97920.0_cp* &
      &              a**4*b**2*dm**2*dn**3*r**4-702208.0_cp*a**4*b**2*dm**2*dn**2*r**4-1435520.0_cp &
      &              *a**4*b**2*dm**2*dn*r**4+6738432.0_cp*a**4*b**2*dm**2*r**4+107520.0_cp &
      &              *a**4*b**2*dn**6*r**4+288960.0_cp*a**4*b**2*dn**5*r**4-3975360.0_cp &
      &              *a**4*b**2*dn**4*r**4-7737920.0_cp*a**4*b**2*dn**3*r**4+55112384.0_cp &
      &              *a**4*b**2*dn**2*r**4+55271040.0_cp*a**4*b**2*dn*r**4-279495936.0_cp &
      &              *a**4*b**2*r**4-512.0_cp*a**4*dm**4*dn**2*r**6-2560.0_cp*a**4*dm**4 &
      &              *dn*r**6+12288.0_cp*a**4*dm**4*r**6-2048.0_cp*a**4*dm**2*dn**4*r**6 &
      &              -13632.0_cp*a**4*dm**2*dn**3*r**6+98432.0_cp*a**4*dm**2*dn**2*r**6+209088.0_cp &
      &              *a**4*dm**2*dn*r**6-979200.0_cp*a**4*dm**2*r**6-7680.0_cp*a**4*dn** &
      &              6*r**6-20160.0_cp*a**4*dn**5*r**6+277696.0_cp*a**4*dn**4*r**6+527232.0_cp &
      &              *a**4*dn**3*r**6-3763264.0_cp*a**4*dn**2*r**6-3682368.0_cp*a**4*dn* &
      &              r**6+18685440.0_cp*a**4*r**6+2560.0_cp*a**2*b**8*dm**2*dn**4+15744.0_cp &
      &              *a**2*b**8*dm**2*dn**3-64768.0_cp*a**2*b**8*dm**2*dn**2-141696.0_cp &
      &              *a**2*b**8*dm**2*dn+375552.0_cp*a**2*b**8*dm**2+84480.0_cp*a**2*b** &
      &              8*dn**6+232320.0_cp*a**2*b**8*dn**5-2299520.0_cp*a**2*b**8*dn**4-4651392.0_cp &
      &              *a**2*b**8*dn**3+21920384.0_cp*a**2*b**8*dn**2+23044608.0_cp*a**2*b**8 &
      &              *dn-72608256.0_cp*a**2*b**8-8192.0_cp*a**2*b**6*dm**2*dn**4*r**2-51072.0_cp &
      &              *a**2*b**6*dm**2*dn**3*r**2+208640.0_cp*a**2*b**6*dm**2*dn**2*r**2+459648.0_cp &
      &              *a**2*b**6*dm**2*dn*r**2-1214208.0_cp*a**2*b**6*dm**2*r**2-184320.0_cp &
      &              *a**2*b**6*dn**6*r**2-502272.0_cp*a**2*b**6*dn**5*r**2+4989440.0_cp &
      &              *a**2*b**6*dn**4*r**2+10040448.0_cp*a**2*b**6*dn**3*r**2-47369984.0_cp &
      &              *a**2*b**6*dn**2*r**2-49680000.0_cp*a**2*b**6*dn*r**2+156554496.0_cp &
      &              *a**2*b**6*r**2+9216.0_cp*a**2*b**4*dm**2*dn**4*r**4+58752.0_cp*a** &
      &              2*b**4*dm**2*dn**3*r**4-237312.0_cp*a**2*b**4*dm**2*dn**2*r**4-528768.0_cp &
      &              *a**2*b**4*dm**2*dn*r**4+1389312.0_cp*a**2*b**4*dm**2*r**4+129024.0_cp &
      &              *a**2*b**4*dn**6*r**4+346752.0_cp*a**2*b**4*dn**5*r**4-3461760.0_cp &
      &              *a**2*b**4*dn**4*r**4-6896640.0_cp*a**2*b**4*dn**3*r**4+32620416.0_cp &
      &              *a**2*b**4*dn**2*r**4+33982848.0_cp*a**2*b**4*dn*r**4-107239680.0_cp &
      &              *a**2*b**4*r**4-4096.0_cp*a**2*b**2*dm**2*dn**4*r**6-27264.0_cp*a** &
      &              2*b**2*dm**2*dn**3*r**6+107776.0_cp*a**2*b**2*dm**2*dn**2*r**6+245376.0_cp &
      &              *a**2*b**2*dm**2*dn*r**6-638208.0_cp*a**2*b**2*dm**2*r**6-30720.0_cp &
      &              *a**2*b**2*dn**6*r**6-80640.0_cp*a**2*b**2*dn**5*r**6+810752.0_cp*a**2 &
      &              *b**2*dn**4*r**6+1578624.0_cp*a**2*b**2*dn**3*r**6-7514624.0_cp*a** &
      &              2*b**2*dn**2*r**6-7675776.0_cp*a**2*b**2*dn*r**6+24355584.0_cp*a**2 &
      &              *b**2*r**6+512.0_cp*a**2*dm**2*dn**4*r**8+3840.0_cp*a**2*dm**2*dn** &
      &              3*r**8-14336.0_cp*a**2*dm**2*dn**2*r**8-34560.0_cp*a**2*dm**2*dn*r**8+ &
      &              87552.0_cp*a**2*dm**2*r**8+1536.0_cp*a**2*dn**6*r**8+3840.0_cp*a**2 &
      &              *dn**5*r**8-38912.0_cp*a**2*dn**4*r**8-71040.0_cp*a**2*dn**3*r**8+343808.0_cp &
      &              *a**2*dn**2*r**8+328320.0_cp*a**2*dn*r**8-1062144.0_cp*a**2*r**8+6144.0_cp &
      &              *b**10*dn**6+16896.0_cp*b**10*dn**5-102912.0_cp*b**10*dn**4-219648.0_cp &
      &              *b**10*dn**3+520704.0_cp*b**10*dn**2+608256.0_cp*b**10*dn-829440.0_cp &
      &              *b**10-20480.0_cp*b**8*dn**6*r**2-55808.0_cp*b**8*dn**5*r**2+342528.0_cp &
      &              *b**8*dn**4*r**2+725504.0_cp*b**8*dn**3*r**2-1729024.0_cp*b**8*dn** &
      &              2*r**2-2009088.0_cp*b**8*dn*r**2+2746368.0_cp*b**8*r**2+24576.0_cp* &
      &              b**6*dn**6*r**4+66048.0_cp*b**6*dn**5*r**4-410112.0_cp*b**6*dn**4*r**4 &
      &              -858624.0_cp*b**6*dn**3*r**4+2062848.0_cp*b**6*dn**2*r**4+2377728.0_cp &
      &              *b**6*dn*r**4-3262464.0_cp*b**6*r**4-12288.0_cp*b**4*dn**6*r**6-32256.0_cp &
      &              *b**4*dn**5*r**6+204288.0_cp*b**4*dn**4*r**6+419328.0_cp*b**4*dn**3 &
      &              *r**6-1021440.0_cp*b**4*dn**2*r**6-1161216.0_cp*b**4*dn*r**6+1603584.0_cp &
      &              *b**4*r**6+2048.0_cp*b**2*dn**6*r**8+5120.0_cp*b**2*dn**5*r**8-33792.0_cp &
      &              *b**2*dn**4*r**8-66560.0_cp*b**2*dn**3*r**8+166912.0_cp*b**2*dn**2* &
      &              r**8+184320.0_cp*b**2*dn*r**8-258048.0_cp*b**2*r**8)/(1024.0_cp*dn* &
      &              (dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+1)=(7.0_cp*a**12*dm**4*dn**2-343.0_cp*a**12*dm**4+42.0_cp &
      &              *a**12*dm**2*dn**4-4480.0_cp*a**12*dm**2*dn**2+90958.0_cp*a**12*dm**2+ &
      &              231.0_cp*a**12*dn**6-17745.0_cp*a**12*dn**4+528570.0_cp*a**12*dn**2 &
      &              -5927544.0_cp*a**12+336.0_cp*a**10*b**2*dm**4*dn**2-13104.0_cp*a**10 &
      &              *b**2*dm**4+2520.0_cp*a**10*b**2*dm**2*dn**4-221592.0_cp*a**10*b**2 &
      &              *dm**2*dn**2+3735648.0_cp*a**10*b**2*dm**2+16632.0_cp*a**10*b**2*dn**6 &
      &              -1103508.0_cp*a**10*b**2*dn**4+27906732.0_cp*a**10*b**2*dn**2-261718128.0_cp &
      &              *a**10*b**2-48.0_cp*a**10*dm**4*dn**2*r**2+1872.0_cp*a**10*dm**4*r**2- &
      &              224.0_cp*a**10*dm**2*dn**4*r**2+20096.0_cp*a**10*dm**2*dn**2*r**2-343080.0_cp &
      &              *a**10*dm**2*r**2-1008.0_cp*a**10*dn**6*r**2+66248.0_cp*a**10*dn**4 &
      &              *r**2-1671536.0_cp*a**10*dn**2*r**2+15703488.0_cp*a**10*r**2+1680.0_cp &
      &              *a**8*b**4*dm**4*dn**2-48720.0_cp*a**8*b**4*dm**4+16800.0_cp*a**8*b**4 &
      &              *dm**2*dn**4-1162560.0_cp*a**8*b**4*dm**2*dn**2+15587040.0_cp*a**8* &
      &              b**4*dm**2+138600.0_cp*a**8*b**4*dn**6-7744800.0_cp*a**8*b**4*dn**4 &
      &              +161114520.0_cp*a**8*b**4*dn**2-1216464480.0_cp*a**8*b**4-1440.0_cp &
      &              *a**8*b**2*dm**4*dn**2*r**2+41760.0_cp*a**8*b**2*dm**4*r**2-8960.0_cp &
      &              *a**8*b**2*dm**2*dn**4*r**2+632000.0_cp*a**8*b**2*dm**2*dn**2*r**2-8569440.0_cp &
      &              *a**8*b**2*dm**2*r**2-50400.0_cp*a**8*b**2*dn**6*r**2+2791040.0_cp* &
      &              a**8*b**2*dn**4*r**2-57911600.0_cp*a**8*b**2*dn**2*r**2+437690160.0_cp &
      &              *a**8*b**2*r**2+144.0_cp*a**8*dm**4*dn**2*r**4-4176.0_cp*a**8*dm**4 &
      &              *r**4+480.0_cp*a**8*dm**2*dn**4*r**4-34992.0_cp*a**8*dm**2*dn**2*r**4+ &
      &              483888.0_cp*a**8*dm**2*r**4+1680.0_cp*a**8*dn**6*r**4-91680.0_cp*a**8* &
      &              dn**4*r**4+1884336.0_cp*a**8*dn**2*r**4-14145984.0_cp*a**8*r**4+1792.0_cp &
      &              *a**6*b**6*dm**4*dn**2-34048.0_cp*a**6*b**6*dm**4+26880.0_cp*a**6*b**6 &
      &              *dm**2*dn**4-1356544.0_cp*a**6*b**6*dm**2*dn**2+13490176.0_cp*a**6* &
      &              b**6*dm**2+295680.0_cp*a**6*b**6*dn**6-13426560.0_cp*a**6*b**6*dn** &
      &              4+219343488.0_cp*a**6*b**6*dn**2-1258209792.0_cp*a**6*b**6-3840.0_cp &
      &              *a**6*b**4*dm**4*dn**2*r**2+72960.0_cp*a**6*b**4*dm**4*r**2-35840.0_cp &
      &              *a**6*b**4*dm**2*dn**4*r**2+1840640.0_cp*a**6*b**4*dm**2*dn**2*r**2 &
      &              -18472320.0_cp*a**6*b**4*dm**2*r**2-268800.0_cp*a**6*b**4*dn**6*r** &
      &              2+12104960.0_cp*a**6*b**4*dn**4*r**2-197164160.0_cp*a**6*b**4*dn**2 &
      &              *r**2+1131022080.0_cp*a**6*b**4*r**2+2304.0_cp*a**6*b**2*dm**4*dn** &
      &              2*r**4-43776.0_cp*a**6*b**2*dm**4*r**4+11520.0_cp*a**6*b**2*dm**2*dn**4 &
      &              *r**4-609792.0_cp*a**6*b**2*dm**2*dn**2*r**4+6217728.0_cp*a**6*b**2 &
      &              *dm**2*r**4+53760.0_cp*a**6*b**2*dn**6*r**4-2388480.0_cp*a**6*b**2*dn**4 &
      &              *r**4+38555136.0_cp*a**6*b**2*dn**2*r**4-219826944.0_cp*a**6*b**2*r**4 &
      &              -256.0_cp*a**6*dm**4*dn**2*r**6+4864.0_cp*a**6*dm**4*r**6-512.0_cp* &
      &              a**6*dm**2*dn**4*r**6+28928.0_cp*a**6*dm**2*dn**2*r**6-305280.0_cp* &
      &              a**6*dm**2*r**6-1280.0_cp*a**6*dn**6*r**6+55424.0_cp*a**6*dn**4*r** &
      &              6-872448.0_cp*a**6*dn**2*r**6+4866048.0_cp*a**6*r**6+384.0_cp*a**4* &
      &              b**8*dm**4*dn**2-3456.0_cp*a**4*b**8*dm**4+11520.0_cp*a**4*b**8*dm**2* &
      &              dn**4-365568.0_cp*a**4*b**8*dm**2*dn**2+2356992.0_cp*a**4*b**8*dm** &
      &              2+190080.0_cp*a**4*b**8*dn**6-6641280.0_cp*a**4*b**8*dn**4+79082496.0_cp &
      &              *a**4*b**8*dn**2-312367104.0_cp*a**4*b**8-1536.0_cp*a**4*b**6*dm**4 &
      &              *dn**2*r**2+13824.0_cp*a**4*b**6*dm**4*r**2-28672.0_cp*a**4*b**6*dm**2 &
      &              *dn**4*r**2+922624.0_cp*a**4*b**6*dm**2*dn**2*r**2-5981184.0_cp*a** &
      &              4*b**6*dm**2*r**2-322560.0_cp*a**4*b**6*dn**6*r**2+11189248.0_cp*a**4* &
      &              b**6*dn**4*r**2-132779776.0_cp*a**4*b**6*dn**2*r**2+523835136.0_cp* &
      &              a**4*b**6*r**2+2304.0_cp*a**4*b**4*dm**4*dn**2*r**4-20736.0_cp*a**4 &
      &              *b**4*dm**4*r**4+23040.0_cp*a**4*b**4*dm**2*dn**4*r**4-759552.0_cp* &
      &              a**4*b**4*dm**2*dn**2*r**4+4969728.0_cp*a**4*b**4*dm**2*r**4+161280.0_cp &
      &              *a**4*b**4*dn**6*r**4-5529600.0_cp*a**4*b**4*dn**4*r**4+65089536.0_cp &
      &              *a**4*b**4*dn**2*r**4-255481344.0_cp*a**4*b**4*r**4-1536.0_cp*a**4* &
      &              b**2*dm**4*dn**2*r**6+13824.0_cp*a**4*b**2*dm**4*r**6-6144.0_cp*a** &
      &              4*b**2*dm**2*dn**4*r**6+213504.0_cp*a**4*b**2*dm**2*dn**2*r**6-1423872.0_cp &
      &              *a**4*b**2*dm**2*r**6-23040.0_cp*a**4*b**2*dn**6*r**6+772608.0_cp*a**4 &
      &              *b**2*dn**4*r**6-8912640.0_cp*a**4*b**2*dn**2*r**6+34428672.0_cp*a**4* &
      &              b**2*r**6+384.0_cp*a**4*dm**4*dn**2*r**8-3456.0_cp*a**4*dm**4*r**8+256.0_cp &
      &              *a**4*dm**2*dn**4*r**8-11008.0_cp*a**4*dm**2*dn**2*r**8+78336.0_cp* &
      &              a**4*dm**2*r**8+384.0_cp*a**4*dn**6*r**8-12160.0_cp*a**4*dn**4*r**8 &
      &              +131584.0_cp*a**4*dn**2*r**8-479232.0_cp*a**4*r**8+1024.0_cp*a**2*b**10 &
      &              *dm**2*dn**4-13312.0_cp*a**2*b**10*dm**2*dn**2+36864.0_cp*a**2*b**10 &
      &              *dm**2+33792.0_cp*a**2*b**10*dn**6-826880.0_cp*a**2*b**10*dn**4+6255104.0_cp &
      &              *a**2*b**10*dn**2-13953024.0_cp*a**2*b**10-4096.0_cp*a**2*b**8*dm** &
      &              2*dn**4*r**2+53248.0_cp*a**2*b**8*dm**2*dn**2*r**2-147456.0_cp*a**2 &
      &              *b**8*dm**2*r**2-92160.0_cp*a**2*b**8*dn**6*r**2+2243584.0_cp*a**2* &
      &              b**8*dn**4*r**2-16909312.0_cp*a**2*b**8*dn**2*r**2+37638144.0_cp*a**2* &
      &              b**8*r**2+6144.0_cp*a**2*b**6*dm**2*dn**4*r**4-79872.0_cp*a**2*b**6 &
      &              *dm**2*dn**2*r**4+221184.0_cp*a**2*b**6*dm**2*r**4+86016.0_cp*a**2* &
      &              b**6*dn**6*r**4-2076672.0_cp*a**2*b**6*dn**4*r**4+15556608.0_cp*a** &
      &              2*b**6*dn**2*r**4-34504704.0_cp*a**2*b**6*r**4-4096.0_cp*a**2*b**4*dm**2 &
      &              *dn**4*r**6+53248.0_cp*a**2*b**4*dm**2*dn**2*r**6-147456.0_cp*a**2* &
      &              b**4*dm**2*r**6-30720.0_cp*a**2*b**4*dn**6*r**6+730112.0_cp*a**2*b**4* &
      &              dn**4*r**6-5405696.0_cp*a**2*b**4*dn**2*r**6+11907072.0_cp*a**2*b** &
      &              4*r**6+1024.0_cp*a**2*b**2*dm**2*dn**4*r**8-13312.0_cp*a**2*b**2*dm**2 &
      &              *dn**2*r**8+36864.0_cp*a**2*b**2*dm**2*r**8+3072.0_cp*a**2*b**2*dn**6* &
      &              r**8-70144.0_cp*a**2*b**2*dn**4*r**8+503296.0_cp*a**2*b**2*dn**2*r**8- &
      &              1087488.0_cp*a**2*b**2*r**8+1024.0_cp*b**12*dn**6-14336.0_cp*b**12*dn**4 &
      &              +50176.0_cp*b**12*dn**2-36864.0_cp*b**12-4096.0_cp*b**10*dn**6*r**2 &
      &              +57344.0_cp*b**10*dn**4*r**2-200704.0_cp*b**10*dn**2*r**2+147456.0_cp &
      &              *b**10*r**2+6144.0_cp*b**8*dn**6*r**4-86016.0_cp*b**8*dn**4*r**4+301056.0_cp &
      &              *b**8*dn**2*r**4-221184.0_cp*b**8*r**4-4096.0_cp*b**6*dn**6*r**6+57344.0_cp &
      &              *b**6*dn**4*r**6-200704.0_cp*b**6*dn**2*r**6+147456.0_cp*b**6*r**6+1024.0_cp &
      &              *b**4*dn**6*r**8-14336.0_cp*b**4*dn**4*r**8+50176.0_cp*b**4*dn**2*r**8 &
      &              -36864.0_cp*b**4*r**8)/(1024.0_cp*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+2)=a*b*(56.0_cp*a**10*dm**4*dn**2-280.0_cp*a**10*dm**4*dn &
      &              -3024.0_cp*a**10*dm**4+420.0_cp*a**10*dm**2*dn**4-2583.0_cp*a**10*dm**2 &
      &              *dn**3-42098.0_cp*a**10*dm**2*dn**2+77287.0_cp*a**10*dm**2*dn+831222.0_cp &
      &              *a**10*dm**2+2772.0_cp*a**10*dn**6-7623.0_cp*a**10*dn**5-191541.0_cp &
      &              *a**10*dn**4+366723.0_cp*a**10*dn**3+5285469.0_cp*a**10*dn**2-5122404.0_cp &
      &              *a**10*dn-56029428.0_cp*a**10+672.0_cp*a**8*b**2*dm**4*dn**2-3360.0_cp &
      &              *a**8*b**2*dm**4*dn-29568.0_cp*a**8*b**2*dm**4+6720.0_cp*a**8*b**2*dm**2 &
      &              *dn**4-41328.0_cp*a**8*b**2*dm**2*dn**3-547680.0_cp*a**8*b**2*dm**2 &
      &              *dn**2+1020432.0_cp*a**8*b**2*dm**2*dn+8924160.0_cp*a**8*b**2*dm**2 &
      &              +55440.0_cp*a**8*b**2*dn**6-152460.0_cp*a**8*b**2*dn**5-3250380.0_cp &
      &              *a**8*b**2*dn**4+6263964.0_cp*a**8*b**2*dn**3+74991756.0_cp*a**8*b**2* &
      &              dn**2-73298736.0_cp*a**8*b**2*dn-656967024.0_cp*a**8*b**2-288.0_cp* &
      &              a**8*dm**4*dn**2*r**2+1440.0_cp*a**8*dm**4*dn*r**2+12672.0_cp*a**8*dm**4 &
      &              *r**2-1792.0_cp*a**8*dm**2*dn**4*r**2+11172.0_cp*a**8*dm**2*dn**3*r**2 &
      &              +148744.0_cp*a**8*dm**2*dn**2*r**2-280908.0_cp*a**8*dm**2*dn*r**2-2456064.0_cp &
      &              *a**8*dm**2*r**2-10080.0_cp*a**8*dn**6*r**2+27468.0_cp*a**8*dn**5*r**2 &
      &              +585676.0_cp*a**8*dn**4*r**2-1125096.0_cp*a**8*dn**3*r**2-13475428.0_cp &
      &              *a**8*dn**2*r**2+13191876.0_cp*a**8*dn*r**2+118223856.0_cp*a**8*r** &
      &              2+1344.0_cp*a**6*b**4*dm**4*dn**2-6720.0_cp*a**6*b**4*dm**4*dn-45696.0_cp &
      &              *a**6*b**4*dm**4+20160.0_cp*a**6*b**4*dm**2*dn**4-123984.0_cp*a**6* &
      &              b**4*dm**2*dn**3-1265376.0_cp*a**6*b**4*dm**2*dn**2+2412816.0_cp*a**6* &
      &              b**4*dm**2*dn+16240224.0_cp*a**6*b**4*dm**2+221760.0_cp*a**6*b**4*dn**6 &
      &              -609840.0_cp*a**6*b**4*dn**5-10679760.0_cp*a**6*b**4*dn**4+20773872.0_cp &
      &              *a**6*b**4*dn**3+198127440.0_cp*a**6*b**4*dn**2-196112448.0_cp*a**6 &
      &              *b**4*dn-1372472640.0_cp*a**6*b**4-1920.0_cp*a**6*b**2*dm**4*dn**2* &
      &              r**2+9600.0_cp*a**6*b**2*dm**4*dn*r**2+65280.0_cp*a**6*b**2*dm**4*r**2 &
      &              -17920.0_cp*a**6*b**2*dm**2*dn**4*r**2+111720.0_cp*a**6*b**2*dm**2*dn**3 &
      &              *r**2+1143760.0_cp*a**6*b**2*dm**2*dn**2*r**2-2207880.0_cp*a**6*b** &
      &              2*dm**2*dn*r**2-14854320.0_cp*a**6*b**2*dm**2*r**2-134400.0_cp*a**6 &
      &              *b**2*dn**6*r**2+366240.0_cp*a**6*b**2*dn**5*r**2+6418720.0_cp*a**6 &
      &              *b**2*dn**4*r**2-12441240.0_cp*a**6*b**2*dn**3*r**2-118703440.0_cp* &
      &              a**6*b**2*dn**2*r**2+117578520.0_cp*a**6*b**2*dn*r**2+822756240.0_cp &
      &              *a**6*b**2*r**2+576.0_cp*a**6*dm**4*dn**2*r**4-2880.0_cp*a**6*dm**4 &
      &              *dn*r**4-19584.0_cp*a**6*dm**4*r**4+2880.0_cp*a**6*dm**2*dn**4*r**4 &
      &              -18360.0_cp*a**6*dm**2*dn**3*r**4-189168.0_cp*a**6*dm**2*dn**2*r**4 &
      &              +373080.0_cp*a**6*dm**2*dn*r**4+2508432.0_cp*a**6*dm**2*r**4+13440.0_cp &
      &              *a**6*dn**6*r**4-36120.0_cp*a**6*dn**5*r**4-633240.0_cp*a**6*dn**4* &
      &              r**4+1216080.0_cp*a**6*dn**3*r**4+11601384.0_cp*a**6*dn**2*r**4-11407320.0_cp &
      &              *a**6*dn*r**4-79859376.0_cp*a**6*r**4+512.0_cp*a**4*b**6*dm**4*dn** &
      &              2-2560.0_cp*a**4*b**6*dm**4*dn-12288.0_cp*a**4*b**6*dm**4+15360.0_cp &
      &              *a**4*b**6*dm**2*dn**4-94464.0_cp*a**4*b**6*dm**2*dn**3-676352.0_cp &
      &              *a**4*b**6*dm**2*dn**2+1344256.0_cp*a**4*b**6*dm**2*dn+6325248.0_cp &
      &              *a**4*b**6*dm**2+253440.0_cp*a**4*b**6*dn**6-696960.0_cp*a**4*b**6*dn**5 &
      &              -9552000.0_cp*a**4*b**6*dn**4+18847872.0_cp*a**4*b**6*dn**3+134078592.0_cp &
      &              *a**4*b**6*dn**2-135479808.0_cp*a**4*b**6*dn-684661248.0_cp*a**4*b**6- &
      &              1536.0_cp*a**4*b**4*dm**4*dn**2*r**2+7680.0_cp*a**4*b**4*dm**4*dn*r**2 &
      &              +36864.0_cp*a**4*b**4*dm**4*r**2-28672.0_cp*a**4*b**4*dm**2*dn**4*r**2 &
      &              +178752.0_cp*a**4*b**4*dm**2*dn**3*r**2+1280128.0_cp*a**4*b**4*dm** &
      &              2*dn**2*r**2-2570688.0_cp*a**4*b**4*dm**2*dn*r**2-12084480.0_cp*a** &
      &              4*b**4*dm**2*r**2-322560.0_cp*a**4*b**4*dn**6*r**2+878976.0_cp*a**4 &
      &              *b**4*dn**5*r**2+12068224.0_cp*a**4*b**4*dn**4*r**2-23714880.0_cp*a**4 &
      &              *b**4*dn**3*r**2-168782848.0_cp*a**4*b**4*dn**2*r**2+170455104.0_cp &
      &              *a**4*b**4*dn*r**2+861320448.0_cp*a**4*b**4*r**2+1536.0_cp*a**4*b** &
      &              2*dm**4*dn**2*r**4-7680.0_cp*a**4*b**2*dm**4*dn*r**4-36864.0_cp*a** &
      &              4*b**2*dm**4*r**4+15360.0_cp*a**4*b**2*dm**2*dn**4*r**4-97920.0_cp* &
      &              a**4*b**2*dm**2*dn**3*r**4-702208.0_cp*a**4*b**2*dm**2*dn**2*r**4+1435520.0_cp &
      &              *a**4*b**2*dm**2*dn*r**4+6738432.0_cp*a**4*b**2*dm**2*r**4+107520.0_cp &
      &              *a**4*b**2*dn**6*r**4-288960.0_cp*a**4*b**2*dn**5*r**4-3975360.0_cp &
      &              *a**4*b**2*dn**4*r**4+7737920.0_cp*a**4*b**2*dn**3*r**4+55112384.0_cp &
      &              *a**4*b**2*dn**2*r**4-55271040.0_cp*a**4*b**2*dn*r**4-279495936.0_cp &
      &              *a**4*b**2*r**4-512.0_cp*a**4*dm**4*dn**2*r**6+2560.0_cp*a**4*dm**4 &
      &              *dn*r**6+12288.0_cp*a**4*dm**4*r**6-2048.0_cp*a**4*dm**2*dn**4*r**6 &
      &              +13632.0_cp*a**4*dm**2*dn**3*r**6+98432.0_cp*a**4*dm**2*dn**2*r**6-209088.0_cp &
      &              *a**4*dm**2*dn*r**6-979200.0_cp*a**4*dm**2*r**6-7680.0_cp*a**4*dn** &
      &              6*r**6+20160.0_cp*a**4*dn**5*r**6+277696.0_cp*a**4*dn**4*r**6-527232.0_cp &
      &              *a**4*dn**3*r**6-3763264.0_cp*a**4*dn**2*r**6+3682368.0_cp*a**4*dn* &
      &              r**6+18685440.0_cp*a**4*r**6+2560.0_cp*a**2*b**8*dm**2*dn**4-15744.0_cp &
      &              *a**2*b**8*dm**2*dn**3-64768.0_cp*a**2*b**8*dm**2*dn**2+141696.0_cp &
      &              *a**2*b**8*dm**2*dn+375552.0_cp*a**2*b**8*dm**2+84480.0_cp*a**2*b** &
      &              8*dn**6-232320.0_cp*a**2*b**8*dn**5-2299520.0_cp*a**2*b**8*dn**4+4651392.0_cp &
      &              *a**2*b**8*dn**3+21920384.0_cp*a**2*b**8*dn**2-23044608.0_cp*a**2*b**8 &
      &              *dn-72608256.0_cp*a**2*b**8-8192.0_cp*a**2*b**6*dm**2*dn**4*r**2+51072.0_cp &
      &              *a**2*b**6*dm**2*dn**3*r**2+208640.0_cp*a**2*b**6*dm**2*dn**2*r**2-459648.0_cp &
      &              *a**2*b**6*dm**2*dn*r**2-1214208.0_cp*a**2*b**6*dm**2*r**2-184320.0_cp &
      &              *a**2*b**6*dn**6*r**2+502272.0_cp*a**2*b**6*dn**5*r**2+4989440.0_cp &
      &              *a**2*b**6*dn**4*r**2-10040448.0_cp*a**2*b**6*dn**3*r**2-47369984.0_cp &
      &              *a**2*b**6*dn**2*r**2+49680000.0_cp*a**2*b**6*dn*r**2+156554496.0_cp &
      &              *a**2*b**6*r**2+9216.0_cp*a**2*b**4*dm**2*dn**4*r**4-58752.0_cp*a** &
      &              2*b**4*dm**2*dn**3*r**4-237312.0_cp*a**2*b**4*dm**2*dn**2*r**4+528768.0_cp &
      &              *a**2*b**4*dm**2*dn*r**4+1389312.0_cp*a**2*b**4*dm**2*r**4+129024.0_cp &
      &              *a**2*b**4*dn**6*r**4-346752.0_cp*a**2*b**4*dn**5*r**4-3461760.0_cp &
      &              *a**2*b**4*dn**4*r**4+6896640.0_cp*a**2*b**4*dn**3*r**4+32620416.0_cp &
      &              *a**2*b**4*dn**2*r**4-33982848.0_cp*a**2*b**4*dn*r**4-107239680.0_cp &
      &              *a**2*b**4*r**4-4096.0_cp*a**2*b**2*dm**2*dn**4*r**6+27264.0_cp*a** &
      &              2*b**2*dm**2*dn**3*r**6+107776.0_cp*a**2*b**2*dm**2*dn**2*r**6-245376.0_cp &
      &              *a**2*b**2*dm**2*dn*r**6-638208.0_cp*a**2*b**2*dm**2*r**6-30720.0_cp &
      &              *a**2*b**2*dn**6*r**6+80640.0_cp*a**2*b**2*dn**5*r**6+810752.0_cp*a**2 &
      &              *b**2*dn**4*r**6-1578624.0_cp*a**2*b**2*dn**3*r**6-7514624.0_cp*a** &
      &              2*b**2*dn**2*r**6+7675776.0_cp*a**2*b**2*dn*r**6+24355584.0_cp*a**2 &
      &              *b**2*r**6+512.0_cp*a**2*dm**2*dn**4*r**8-3840.0_cp*a**2*dm**2*dn** &
      &              3*r**8-14336.0_cp*a**2*dm**2*dn**2*r**8+34560.0_cp*a**2*dm**2*dn*r**8+ &
      &              87552.0_cp*a**2*dm**2*r**8+1536.0_cp*a**2*dn**6*r**8-3840.0_cp*a**2 &
      &              *dn**5*r**8-38912.0_cp*a**2*dn**4*r**8+71040.0_cp*a**2*dn**3*r**8+343808.0_cp &
      &              *a**2*dn**2*r**8-328320.0_cp*a**2*dn*r**8-1062144.0_cp*a**2*r**8+6144.0_cp &
      &              *b**10*dn**6-16896.0_cp*b**10*dn**5-102912.0_cp*b**10*dn**4+219648.0_cp &
      &              *b**10*dn**3+520704.0_cp*b**10*dn**2-608256.0_cp*b**10*dn-829440.0_cp &
      &              *b**10-20480.0_cp*b**8*dn**6*r**2+55808.0_cp*b**8*dn**5*r**2+342528.0_cp &
      &              *b**8*dn**4*r**2-725504.0_cp*b**8*dn**3*r**2-1729024.0_cp*b**8*dn** &
      &              2*r**2+2009088.0_cp*b**8*dn*r**2+2746368.0_cp*b**8*r**2+24576.0_cp* &
      &              b**6*dn**6*r**4-66048.0_cp*b**6*dn**5*r**4-410112.0_cp*b**6*dn**4*r**4 &
      &              +858624.0_cp*b**6*dn**3*r**4+2062848.0_cp*b**6*dn**2*r**4-2377728.0_cp &
      &              *b**6*dn*r**4-3262464.0_cp*b**6*r**4-12288.0_cp*b**4*dn**6*r**6+32256.0_cp &
      &              *b**4*dn**5*r**6+204288.0_cp*b**4*dn**4*r**6-419328.0_cp*b**4*dn**3 &
      &              *r**6-1021440.0_cp*b**4*dn**2*r**6+1161216.0_cp*b**4*dn*r**6+1603584.0_cp &
      &              *b**4*r**6+2048.0_cp*b**2*dn**6*r**8-5120.0_cp*b**2*dn**5*r**8-33792.0_cp &
      &              *b**2*dn**4*r**8+66560.0_cp*b**2*dn**3*r**8+166912.0_cp*b**2*dn**2* &
      &              r**8-184320.0_cp*b**2*dn*r**8-258048.0_cp*b**2*r**8)/(1024.0_cp*dn* &
      &              (dn-3.0_cp)*(dn-2.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+3)=-a**2*(-2.0_cp*a**10*dm**4*dn**2+50.0_cp*a**10*dm**4*dn &
      &              +288.0_cp*a**10*dm**4-24.0_cp*a**10*dm**2*dn**4+387.0_cp*a**10*dm** &
      &              2*dn**3+2948.0_cp*a**10*dm**2*dn**2-13133.0_cp*a**10*dm**2*dn-76698.0_cp &
      &              *a**10*dm**2-198.0_cp*a**10*dn**6+1089.0_cp*a**10*dn**5+13242.0_cp* &
      &              a**10*dn**4-55341.0_cp*a**10*dn**3-398532.0_cp*a**10*dn**2+845460.0_cp &
      &              *a**10*dn+5010768.0_cp*a**10-56.0_cp*a**8*b**2*dm**4*dn**2+2240.0_cp &
      &              *a**8*b**2*dm**4*dn+10584.0_cp*a**8*b**2*dm**4-1260.0_cp*a**8*b**2*dm**2 &
      &              *dn**4+21924.0_cp*a**8*b**2*dm**2*dn**3+133028.0_cp*a**8*b**2*dm**2 &
      &              *dn**2-629636.0_cp*a**8*b**2*dm**2*dn-3040632.0_cp*a**8*b**2*dm**2-13860.0_cp &
      &              *a**8*b**2*dn**6+76230.0_cp*a**8*b**2*dn**5+781830.0_cp*a**8*b**2*dn**4 &
      &              -3338622.0_cp*a**8*b**2*dn**3-19926522.0_cp*a**8*b**2*dn**2+43387848.0_cp &
      &              *a**8*b**2*dn+213931368.0_cp*a**8*b**2+8.0_cp*a**8*dm**4*dn**2*r**2 &
      &              -320.0_cp*a**8*dm**4*dn*r**2-1512.0_cp*a**8*dm**4*r**2+112.0_cp*a** &
      &              8*dm**2*dn**4*r**2-1974.0_cp*a**8*dm**2*dn**3*r**2-11992.0_cp*a**8*dm**2 &
      &              *dn**2*r**2+57846.0_cp*a**8*dm**2*dn*r**2+279216.0_cp*a**8*dm**2*r**2+ &
      &              840.0_cp*a**8*dn**6*r**2-4578.0_cp*a**8*dn**5*r**2-46984.0_cp*a**8*dn**4 &
      &              *r**2+199850.0_cp*a**8*dn**3*r**2+1193144.0_cp*a**8*dn**2*r**2-2603592.0_cp &
      &              *a**8*dn*r**2-12835872.0_cp*a**8*r**2+10080.0_cp*a**6*b**4*dm**4*dn &
      &              +36960.0_cp*a**6*b**4*dm**4-6720.0_cp*a**6*b**4*dm**2*dn**4+134064.0_cp &
      &              *a**6*b**4*dm**2*dn**3+603456.0_cp*a**6*b**4*dm**2*dn**2-3152016.0_cp &
      &              *a**6*b**4*dm**2*dn-12020064.0_cp*a**6*b**4*dm**2-110880.0_cp*a**6* &
      &              b**4*dn**6+609840.0_cp*a**6*b**4*dn**5+5093760.0_cp*a**6*b**4*dn**4 &
      &              -22426992.0_cp*a**6*b**4*dn**3-106161888.0_cp*a**6*b**4*dn**2+240262848.0_cp &
      &              *a**6*b**4*dn+945689472.0_cp*a**6*b**4-8640.0_cp*a**6*b**2*dm**4*dn &
      &              *r**2-31680.0_cp*a**6*b**2*dm**4*r**2+3584.0_cp*a**6*b**2*dm**2*dn**4* &
      &              r**2-72408.0_cp*a**6*b**2*dm**2*dn**3*r**2-325472.0_cp*a**6*b**2*dm**2 &
      &              *dn**2*r**2+1733832.0_cp*a**6*b**2*dm**2*dn*r**2+6606864.0_cp*a**6* &
      &              b**2*dm**2*r**2+40320.0_cp*a**6*b**2*dn**6*r**2-219744.0_cp*a**6*b**2* &
      &              dn**5*r**2-1838144.0_cp*a**6*b**2*dn**4*r**2+8056776.0_cp*a**6*b**2 &
      &              *dn**3*r**2+38149664.0_cp*a**6*b**2*dn**2*r**2-86457240.0_cp*a**6*b**2 &
      &              *dn*r**2-340250832.0_cp*a**6*b**2*r**2+864.0_cp*a**6*dm**4*dn*r**4+3168.0_cp &
      &              *a**6*dm**4*r**4-192.0_cp*a**6*dm**2*dn**4*r**4+3960.0_cp*a**6*dm** &
      &              2*dn**3*r**4+17760.0_cp*a**6*dm**2*dn**2*r**4-97992.0_cp*a**6*dm**2 &
      &              *dn*r**4-372912.0_cp*a**6*dm**2*r**4-1344.0_cp*a**6*dn**6*r**4+7224.0_cp &
      &              *a**6*dn**5*r**4+60528.0_cp*a**6*dn**4*r**4-262440.0_cp*a**6*dn**3* &
      &              r**4-1243776.0_cp*a**6*dn**2*r**4+2793312.0_cp*a**6*dn*r**4+10998144.0_cp &
      &              *a**6*r**4+448.0_cp*a**4*b**6*dm**4*dn**2+8960.0_cp*a**4*b**6*dm**4 &
      &              *dn+22848.0_cp*a**4*b**6*dm**4-6720.0_cp*a**4*b**6*dm**2*dn**4+185472.0_cp &
      &              *a**4*b**6*dm**2*dn**3+537152.0_cp*a**4*b**6*dm**2*dn**2-3398528.0_cp &
      &              *a**4*b**6*dm**2*dn-9477888.0_cp*a**4*b**6*dm**2-221760.0_cp*a**4*b**6 &
      &              *dn**6+1219680.0_cp*a**4*b**6*dn**5+7865760.0_cp*a**4*b**6*dn**4-36290016.0_cp &
      &              *a**4*b**6*dn**3-126854112.0_cp*a**4*b**6*dn**2+305875584.0_cp*a**4 &
      &              *b**6*dn+900402048.0_cp*a**4*b**6-960.0_cp*a**4*b**4*dm**4*dn**2*r**2- &
      &              19200.0_cp*a**4*b**4*dm**4*dn*r**2-48960.0_cp*a**4*b**4*dm**4*r**2+8960.0_cp &
      &              *a**4*b**4*dm**2*dn**4*r**2-250320.0_cp*a**4*b**4*dm**2*dn**3*r**2-720320.0_cp &
      &              *a**4*b**4*dm**2*dn**2*r**2+4657680.0_cp*a**4*b**4*dm**2*dn*r**2+12971520.0_cp &
      &              *a**4*b**4*dm**2*r**2+201600.0_cp*a**4*b**4*dn**6*r**2-1098720.0_cp &
      &              *a**4*b**4*dn**5*r**2-7105280.0_cp*a**4*b**4*dn**4*r**2+32603760.0_cp &
      &              *a**4*b**4*dn**3*r**2+114029600.0_cp*a**4*b**4*dn**2*r**2-274983120.0_cp &
      &              *a**4*b**4*dn*r**2-809341920.0_cp*a**4*b**4*r**2+576.0_cp*a**4*b**2 &
      &              *dm**4*dn**2*r**4+11520.0_cp*a**4*b**2*dm**4*dn*r**4+29376.0_cp*a** &
      &              4*b**2*dm**4*r**4-2880.0_cp*a**4*b**2*dm**2*dn**4*r**4+82080.0_cp*a**4 &
      &              *b**2*dm**2*dn**3*r**4+233472.0_cp*a**4*b**2*dm**2*dn**2*r**4-1570080.0_cp &
      &              *a**4*b**2*dm**2*dn*r**4-4362048.0_cp*a**4*b**2*dm**2*r**4-40320.0_cp &
      &              *a**4*b**2*dn**6*r**4+216720.0_cp*a**4*b**2*dn**5*r**4+1406880.0_cp &
      &              *a**4*b**2*dn**4*r**4-6380160.0_cp*a**4*b**2*dn**3*r**4-22360896.0_cp &
      &              *a**4*b**2*dn**2*r**4+53420400.0_cp*a**4*b**2*dn*r**4+157343904.0_cp &
      &              *a**4*b**2*r**4-64.0_cp*a**4*dm**4*dn**2*r**6-1280.0_cp*a**4*dm**4*dn &
      &              *r**6-3264.0_cp*a**4*dm**4*r**6+128.0_cp*a**4*dm**2*dn**4*r**6-3792.0_cp &
      &              *a**4*dm**2*dn**3*r**6-10496.0_cp*a**4*dm**2*dn**2*r**6+77328.0_cp* &
      &              a**4*dm**2*dn*r**6+213696.0_cp*a**4*dm**2*r**6+960.0_cp*a**4*dn**6* &
      &              r**6-5040.0_cp*a**4*dn**5*r**6-32896.0_cp*a**4*dn**4*r**6+144816.0_cp &
      &              *a**4*dn**3*r**6+510592.0_cp*a**4*dn**2*r**6-1180224.0_cp*a**4*dn*r**6 &
      &              -3485952.0_cp*a**4*r**6+256.0_cp*a**2*b**8*dm**4*dn**2+1280.0_cp*a**2* &
      &              b**8*dm**4*dn+1536.0_cp*a**2*b**8*dm**4+58752.0_cp*a**2*b**8*dm**2*dn**3 &
      &              +68096.0_cp*a**2*b**8*dm**2*dn**2-775808.0_cp*a**2*b**8*dm**2*dn-1353984.0_cp &
      &              *a**2*b**8*dm**2-126720.0_cp*a**2*b**8*dn**6+696960.0_cp*a**2*b**8*dn**5 &
      &              +3168000.0_cp*a**2*b**8*dn**4-15843456.0_cp*a**2*b**8*dn**3-35665920.0_cp &
      &              *a**2*b**8*dn**2+97288704.0_cp*a**2*b**8*dn+190218240.0_cp*a**2*b** &
      &              8-1024.0_cp*a**2*b**6*dm**4*dn**2*r**2-5120.0_cp*a**2*b**6*dm**4*dn &
      &              *r**2-6144.0_cp*a**2*b**6*dm**4*r**2-147840.0_cp*a**2*b**6*dm**2*dn**3 &
      &              *r**2-167424.0_cp*a**2*b**6*dm**2*dn**2*r**2+1971840.0_cp*a**2*b**6 &
      &              *dm**2*dn*r**2+3430656.0_cp*a**2*b**6*dm**2*r**2+215040.0_cp*a**2*b**6 &
      &              *dn**6*r**2-1171968.0_cp*a**2*b**6*dn**5*r**2-5354496.0_cp*a**2*b** &
      &              6*dn**4*r**2+26585216.0_cp*a**2*b**6*dn**3*r**2+59944960.0_cp*a**2* &
      &              b**6*dn**2*r**2-163149696.0_cp*a**2*b**6*dn*r**2-318991104.0_cp*a** &
      &              2*b**6*r**2+1536.0_cp*a**2*b**4*dm**4*dn**2*r**4+7680.0_cp*a**2*b** &
      &              4*dm**4*dn*r**4+9216.0_cp*a**2*b**4*dm**4*r**4+120960.0_cp*a**2*b** &
      &              4*dm**2*dn**3*r**4+131072.0_cp*a**2*b**4*dm**2*dn**2*r**4-1642880.0_cp &
      &              *a**2*b**4*dm**2*dn*r**4-2842368.0_cp*a**2*b**4*dm**2*r**4-107520.0_cp &
      &              *a**2*b**4*dn**6*r**4+577920.0_cp*a**2*b**4*dn**5*r**4+2661120.0_cp &
      &              *a**2*b**4*dn**4*r**4-13032320.0_cp*a**2*b**4*dn**3*r**4-29526016.0_cp &
      &              *a**2*b**4*dn**2*r**4+79514880.0_cp*a**2*b**4*dn*r**4+155672064.0_cp &
      &              *a**2*b**4*r**4-1024.0_cp*a**2*b**2*dm**4*dn**2*r**6-5120.0_cp*a**2 &
      &              *b**2*dm**4*dn*r**6-6144.0_cp*a**2*b**2*dm**4*r**6-33408.0_cp*a**2* &
      &              b**2*dm**2*dn**3*r**6-32256.0_cp*a**2*b**2*dm**2*dn**2*r**6+473472.0_cp &
      &              *a**2*b**2*dm**2*dn*r**6+808704.0_cp*a**2*b**2*dm**2*r**6+15360.0_cp &
      &              *a**2*b**2*dn**6*r**6-80640.0_cp*a**2*b**2*dn**5*r**6-376320.0_cp*a**2 &
      &              *b**2*dn**4*r**6+1786752.0_cp*a**2*b**2*dn**3*r**6+4098048.0_cp*a** &
      &              2*b**2*dn**2*r**6-10689408.0_cp*a**2*b**2*dn*r**6-21019392.0_cp*a** &
      &              2*b**2*r**6+256.0_cp*a**2*dm**4*dn**2*r**8+1280.0_cp*a**2*dm**4*dn* &
      &              r**8+1536.0_cp*a**2*dm**4*r**8+1536.0_cp*a**2*dm**2*dn**3*r**8+512.0_cp &
      &              *a**2*dm**2*dn**2*r**8-26624.0_cp*a**2*dm**2*dn*r**8-43008.0_cp*a** &
      &              2*dm**2*r**8-256.0_cp*a**2*dn**6*r**8+1280.0_cp*a**2*dn**5*r**8+6144.0_cp &
      &              *a**2*dn**4*r**8-26624.0_cp*a**2*dn**3*r**8-63488.0_cp*a**2*dn**2*r**8 &
      &              +147456.0_cp*a**2*dn*r**8+294912.0_cp*a**2*r**8+512.0_cp*b**10*dm** &
      &              2*dn**4+1536.0_cp*b**10*dm**2*dn**3-3584.0_cp*b**10*dm**2*dn**2-13824.0_cp &
      &              *b**10*dm**2*dn-9216.0_cp*b**10*dm**2-16896.0_cp*b**10*dn**6+92928.0_cp &
      &              *b**10*dn**5+245504.0_cp*b**10*dn**4-1459968.0_cp*b**10*dn**3-1448192.0_cp &
      &              *b**10*dn**2+5612544.0_cp*b**10*dn+5465088.0_cp*b**10-2048.0_cp*b** &
      &              8*dm**2*dn**4*r**2-6144.0_cp*b**8*dm**2*dn**3*r**2+14336.0_cp*b**8*dm**2 &
      &              *dn**2*r**2+55296.0_cp*b**8*dm**2*dn*r**2+36864.0_cp*b**8*dm**2*r** &
      &              2+46080.0_cp*b**8*dn**6*r**2-251136.0_cp*b**8*dn**5*r**2-670720.0_cp &
      &              *b**8*dn**4*r**2+3941376.0_cp*b**8*dn**3*r**2+3943936.0_cp*b**8*dn**2* &
      &              r**2-15130368.0_cp*b**8*dn*r**2-14759424.0_cp*b**8*r**2+3072.0_cp*b**6 &
      &              *dm**2*dn**4*r**4+9216.0_cp*b**6*dm**2*dn**3*r**4-21504.0_cp*b**6*dm**2 &
      &              *dn**2*r**4-82944.0_cp*b**6*dm**2*dn*r**4-55296.0_cp*b**6*dm**2*r** &
      &              4-43008.0_cp*b**6*dn**6*r**4+231168.0_cp*b**6*dn**5*r**4+628224.0_cp &
      &              *b**6*dn**4*r**4-3620352.0_cp*b**6*dn**3*r**4-3677184.0_cp*b**6*dn**2* &
      &              r**4+13858560.0_cp*b**6*dn*r**4+13561344.0_cp*b**6*r**4-2048.0_cp*b**4 &
      &              *dm**2*dn**4*r**6-6144.0_cp*b**4*dm**2*dn**3*r**6+14336.0_cp*b**4*dm**2 &
      &              *dn**2*r**6+55296.0_cp*b**4*dm**2*dn*r**6+36864.0_cp*b**4*dm**2*r** &
      &              6+15360.0_cp*b**4*dn**6*r**6-80640.0_cp*b**4*dn**5*r**6-226304.0_cp &
      &              *b**4*dn**4*r**6+1256448.0_cp*b**4*dn**3*r**6+1315328.0_cp*b**4*dn**2* &
      &              r**6-4776192.0_cp*b**4*dn*r**6-4704768.0_cp*b**4*r**6+512.0_cp*b**2 &
      &              *dm**2*dn**4*r**8+1536.0_cp*b**2*dm**2*dn**3*r**8-3584.0_cp*b**2*dm**2 &
      &              *dn**2*r**8-13824.0_cp*b**2*dm**2*dn*r**8-9216.0_cp*b**2*dm**2*r**8 &
      &              -1536.0_cp*b**2*dn**6*r**8+7680.0_cp*b**2*dn**5*r**8+23296.0_cp*b** &
      &              2*dn**4*r**8-117504.0_cp*b**2*dn**3*r**8-133888.0_cp*b**2*dn**2*r** &
      &              8+435456.0_cp*b**2*dn*r**8+437760.0_cp*b**2*r**8)/(1024.0_cp*dn*(dn &
      &              -3.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+4)=-a**3*b*(24.0_cp*a**8*dm**4*dn**2+600.0_cp*a**8*dm**4 &
      &              *dn+2064.0_cp*a**8*dm**4-60.0_cp*a**8*dm**2*dn**4+4779.0_cp*a**8*dm**2 &
      &              *dn**3+13806.0_cp*a**8*dm**2*dn**2-162591.0_cp*a**8*dm**2*dn-570774.0_cp &
      &              *a**8*dm**2-1980.0_cp*a**8*dn**6+16335.0_cp*a**8*dn**5+97455.0_cp*a**8 &
      &              *dn**4-750987.0_cp*a**8*dn**3-2754471.0_cp*a**8*dn**2+10837188.0_cp &
      &              *a**8*dn+38725884.0_cp*a**8+448.0_cp*a**6*b**2*dm**4*dn**2+6720.0_cp &
      &              *a**6*b**2*dm**4*dn+18368.0_cp*a**6*b**2*dm**4+68544.0_cp*a**6*b**2 &
      &              *dm**2*dn**3+119168.0_cp*a**6*b**2*dm**2*dn**2-1984416.0_cp*a**6*b**2* &
      &              dm**2*dn-5607392.0_cp*a**6*b**2*dm**2-36960.0_cp*a**6*b**2*dn**6+304920.0_cp &
      &              *a**6*b**2*dn**5+1432200.0_cp*a**6*b**2*dn**4-11877432.0_cp*a**6*b**2* &
      &              dn**3-33540360.0_cp*a**6*b**2*dn**2+143995488.0_cp*a**6*b**2*dn+417705120.0_cp &
      &              *a**6*b**2-192.0_cp*a**6*dm**4*dn**2*r**2-2880.0_cp*a**6*dm**4*dn*r**2 &
      &              -7872.0_cp*a**6*dm**4*r**2-18480.0_cp*a**6*dm**2*dn**3*r**2-31392.0_cp &
      &              *a**6*dm**2*dn**2*r**2+546360.0_cp*a**6*dm**2*dn*r**2+1542888.0_cp* &
      &              a**6*dm**2*r**2+6720.0_cp*a**6*dn**6*r**2-54936.0_cp*a**6*dn**5*r** &
      &              2-258888.0_cp*a**6*dn**4*r**2+2132872.0_cp*a**6*dn**3*r**2+6022440.0_cp &
      &              *a**6*dn**2*r**2-25914472.0_cp*a**6*dn*r**2-75163512.0_cp*a**6*r**2 &
      &              +1344.0_cp*a**4*b**4*dm**4*dn**2+12096.0_cp*a**4*b**4*dm**4*dn+24192.0_cp &
      &              *a**4*b**4*dm**4+4032.0_cp*a**4*b**4*dm**2*dn**4+172368.0_cp*a**4*b**4 &
      &              *dm**2*dn**3+81312.0_cp*a**4*b**4*dm**2*dn**2-4139856.0_cp*a**4*b** &
      &              4*dm**2*dn-8824032.0_cp*a**4*b**4*dm**2-133056.0_cp*a**4*b**4*dn**6 &
      &              +1097712.0_cp*a**4*b**4*dn**5+3762864.0_cp*a**4*b**4*dn**4-35051184.0_cp &
      &              *a**4*b**4*dn**3-69008688.0_cp*a**4*b**4*dn**2+343635264.0_cp*a**4* &
      &              b**4*dn+764551872.0_cp*a**4*b**4-1920.0_cp*a**4*b**2*dm**4*dn**2*r**2- &
      &              17280.0_cp*a**4*b**2*dm**4*dn*r**2-34560.0_cp*a**4*b**2*dm**4*r**2-3584.0_cp &
      &              *a**4*b**2*dm**2*dn**4*r**2-154728.0_cp*a**4*b**2*dm**2*dn**3*r**2-65392.0_cp &
      &              *a**4*b**2*dm**2*dn**2*r**2+3788712.0_cp*a**4*b**2*dm**2*dn*r**2+8067312.0_cp &
      &              *a**4*b**2*dm**2*r**2+80640.0_cp*a**4*b**2*dn**6*r**2-659232.0_cp*a**4 &
      &              *b**2*dn**5*r**2-2272480.0_cp*a**4*b**2*dn**4*r**2+20986392.0_cp*a**4* &
      &              b**2*dn**3*r**2+41331568.0_cp*a**4*b**2*dn**2*r**2-206009976.0_cp*a**4 &
      &              *b**2*dn*r**2-458290512.0_cp*a**4*b**2*r**2+576.0_cp*a**4*dm**4*dn**2* &
      &              r**4+5184.0_cp*a**4*dm**4*dn*r**4+10368.0_cp*a**4*dm**4*r**4+576.0_cp &
      &              *a**4*dm**2*dn**4*r**4+25272.0_cp*a**4*dm**2*dn**3*r**4+8400.0_cp*a**4 &
      &              *dm**2*dn**2*r**4-640440.0_cp*a**4*dm**2*dn*r**4-1361232.0_cp*a**4*dm**2 &
      &              *r**4-8064.0_cp*a**4*dn**6*r**4+65016.0_cp*a**4*dn**5*r**4+226152.0_cp &
      &              *a**4*dn**4*r**4-2053104.0_cp*a**4*dn**3*r**4-4066968.0_cp*a**4*dn**2* &
      &              r**4+19989816.0_cp*a**4*dn*r**4+44497584.0_cp*a**4*r**4+768.0_cp*a**2* &
      &              b**6*dm**4*dn**2+3840.0_cp*a**2*b**6*dm**4*dn+4608.0_cp*a**2*b**6*dm**4 &
      &              +7680.0_cp*a**2*b**6*dm**2*dn**4+93312.0_cp*a**2*b**6*dm**2*dn**3-109824.0_cp &
      &              *a**2*b**6*dm**2*dn**2-1822848.0_cp*a**2*b**6*dm**2*dn-2582784.0_cp &
      &              *a**2*b**6*dm**2-126720.0_cp*a**2*b**6*dn**6+1045440.0_cp*a**2*b**6 &
      &              *dn**5+2256960.0_cp*a**2*b**6*dn**4-26041536.0_cp*a**2*b**6*dn**3-28467264.0_cp &
      &              *a**2*b**6*dn**2+194298624.0_cp*a**2*b**6*dn+299586816.0_cp*a**2*b**6- &
      &              2304.0_cp*a**2*b**4*dm**4*dn**2*r**2-11520.0_cp*a**2*b**4*dm**4*dn* &
      &              r**2-13824.0_cp*a**2*b**4*dm**4*r**2-14336.0_cp*a**2*b**4*dm**2*dn**4* &
      &              r**2-175392.0_cp*a**2*b**4*dm**2*dn**3*r**2+216896.0_cp*a**2*b**4*dm**2 &
      &              *dn**2*r**2+3485088.0_cp*a**2*b**4*dm**2*dn*r**2+4928832.0_cp*a**2* &
      &              b**4*dm**2*r**2+161280.0_cp*a**2*b**4*dn**6*r**2-1318464.0_cp*a**2* &
      &              b**4*dn**5*r**2-2876608.0_cp*a**2*b**4*dn**4*r**2+32756640.0_cp*a** &
      &              2*b**4*dn**3*r**2+35897728.0_cp*a**2*b**4*dn**2*r**2-244419936.0_cp &
      &              *a**2*b**4*dn*r**2-376864704.0_cp*a**2*b**4*r**2+2304.0_cp*a**2*b** &
      &              2*dm**4*dn**2*r**4+11520.0_cp*a**2*b**2*dm**4*dn*r**4+13824.0_cp*a**2* &
      &              b**2*dm**4*r**4+7680.0_cp*a**2*b**2*dm**2*dn**4*r**4+95040.0_cp*a** &
      &              2*b**2*dm**2*dn**3*r**4-127872.0_cp*a**2*b**2*dm**2*dn**2*r**4-1945920.0_cp &
      &              *a**2*b**2*dm**2*dn*r**4-2742912.0_cp*a**2*b**2*dm**2*r**4-53760.0_cp &
      &              *a**2*b**2*dn**6*r**4+433440.0_cp*a**2*b**2*dn**5*r**4+962400.0_cp* &
      &              a**2*b**2*dn**4*r**4-10701280.0_cp*a**2*b**2*dn**3*r**4-11877984.0_cp &
      &              *a**2*b**2*dn**2*r**4+79284160.0_cp*a**2*b**2*dn*r**4+122382336.0_cp &
      &              *a**2*b**2*r**4-768.0_cp*a**2*dm**4*dn**2*r**6-3840.0_cp*a**2*dm**4 &
      &              *dn*r**6-4608.0_cp*a**2*dm**4*r**6-1024.0_cp*a**2*dm**2*dn**4*r**6-12960.0_cp &
      &              *a**2*dm**2*dn**3*r**6+20800.0_cp*a**2*dm**2*dn**2*r**6+283680.0_cp &
      &              *a**2*dm**2*dn*r**6+396864.0_cp*a**2*dm**2*r**6+3840.0_cp*a**2*dn** &
      &              6*r**6-30240.0_cp*a**2*dn**5*r**6-69472.0_cp*a**2*dn**4*r**6+733824.0_cp &
      &              *a**2*dn**3*r**6+845728.0_cp*a**2*dn**2*r**6-5295456.0_cp*a**2*dn*r**6 &
      &              -8205120.0_cp*a**2*r**6+2560.0_cp*b**8*dm**2*dn**4+4992.0_cp*b**8*dm**2 &
      &              *dn**3-34048.0_cp*b**8*dm**2*dn**2-98688.0_cp*b**8*dm**2*dn-62208.0_cp &
      &              *b**8*dm**2-28160.0_cp*b**8*dn**6+232320.0_cp*b**8*dn**5+206720.0_cp &
      &              *b**8*dn**4-4155776.0_cp*b**8*dn**3-717696.0_cp*b**8*dn**2+21062144.0_cp &
      &              *b**8*dn+17677824.0_cp*b**8-8192.0_cp*b**6*dm**2*dn**4*r**2-15744.0_cp &
      &              *b**6*dm**2*dn**3*r**2+110336.0_cp*b**6*dm**2*dn**2*r**2+318336.0_cp &
      &              *b**6*dm**2*dn*r**2+200448.0_cp*b**6*dm**2*r**2+61440.0_cp*b**6*dn**6* &
      &              r**2-502272.0_cp*b**6*dn**5*r**2-460288.0_cp*b**6*dn**4*r**2+8967808.0_cp &
      &              *b**6*dn**3*r**2+1616128.0_cp*b**6*dn**2*r**2-45389440.0_cp*b**6*dn &
      &              *r**2-38141184.0_cp*b**6*r**2+9216.0_cp*b**4*dm**2*dn**4*r**4+17280.0_cp &
      &              *b**4*dm**2*dn**3*r**4-126720.0_cp*b**4*dm**2*dn**2*r**4-362880.0_cp &
      &              *b**4*dm**2*dn*r**4-228096.0_cp*b**4*dm**2*r**4-43008.0_cp*b**4*dn**6* &
      &              r**4+346752.0_cp*b**4*dn**5*r**4+333696.0_cp*b**4*dn**4*r**4-6172160.0_cp &
      &              *b**4*dn**3*r**4-1222272.0_cp*b**4*dn**2*r**4+31084928.0_cp*b**4*dn &
      &              *r**4+26191104.0_cp*b**4*r**4-4096.0_cp*b**2*dm**2*dn**4*r**6-7296.0_cp &
      &              *b**2*dm**2*dn**3*r**6+58624.0_cp*b**2*dm**2*dn**2*r**6+165504.0_cp &
      &              *b**2*dm**2*dn*r**6+103680.0_cp*b**2*dm**2*r**6+10240.0_cp*b**2*dn**6* &
      &              r**6-80640.0_cp*b**2*dn**5*r**6-85248.0_cp*b**2*dn**4*r**6+1426560.0_cp &
      &              *b**2*dn**3*r**6+350720.0_cp*b**2*dn**2*r**6-7067520.0_cp*b**2*dn*r**6 &
      &              -5997312.0_cp*b**2*r**6+512.0_cp*dm**2*dn**4*r**8+768.0_cp*dm**2*dn**3 &
      &              *r**8-8192.0_cp*dm**2*dn**2*r**8-22272.0_cp*dm**2*dn*r**8-13824.0_cp &
      &              *dm**2*r**8-512.0_cp*dn**6*r**8+3840.0_cp*dn**5*r**8+5120.0_cp*dn** &
      &              4*r**8-66432.0_cp*dn**3*r**8-26880.0_cp*dn**2*r**8+309888.0_cp*dn*r**8 &
      &              +269568.0_cp*r**8)/(1024.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+5)=a**4*(-17.0_cp*a**8*dm**4*dn**3-162.0_cp*a**8*dm**4*dn**2 &
      &              +293.0_cp*a**8*dm**4*dn+2718.0_cp*a**8*dm**4-30.0_cp*a**8*dm**2*dn**5- &
      &              1437.0_cp*a**8*dm**2*dn**4+5516.0_cp*a**8*dm**2*dn**3+66921.0_cp*a**8* &
      &              dm**2*dn**2-74234.0_cp*a**8*dm**2*dn-729072.0_cp*a**8*dm**2+495.0_cp &
      &              *a**8*dn**7-7425.0_cp*a**8*dn**6+3435.0_cp*a**8*dn**5+332661.0_cp*a**8 &
      &              *dn**4-449754.0_cp*a**8*dn**3-6507084.0_cp*a**8*dn**2+4628376.0_cp* &
      &              a**8*dn+48008160.0_cp*a**8-896.0_cp*a**6*b**2*dm**4*dn**3-5376.0_cp &
      &              *a**6*b**2*dm**4*dn**2+17024.0_cp*a**6*b**2*dm**4*dn+88704.0_cp*a** &
      &              6*b**2*dm**4-2880.0_cp*a**6*b**2*dm**2*dn**5-64512.0_cp*a**6*b**2*dm**2 &
      &              *dn**4+352064.0_cp*a**6*b**2*dm**2*dn**3+2584704.0_cp*a**6*b**2*dm**2* &
      &              dn**2-4664576.0_cp*a**6*b**2*dm**2*dn-25818624.0_cp*a**6*b**2*dm**2 &
      &              +31680.0_cp*a**6*b**2*dn**7-475200.0_cp*a**6*b**2*dn**6+551520.0_cp &
      &              *a**6*b**2*dn**5+17516736.0_cp*a**6*b**2*dn**4-33870240.0_cp*a**6*b**2 &
      &              *dn**3-284757120.0_cp*a**6*b**2*dn**2+315123840.0_cp*a**6*b**2*dn+1841647104.0_cp &
      &              *a**6*b**2+128.0_cp*a**6*dm**4*dn**3*r**2+768.0_cp*a**6*dm**4*dn**2 &
      &              *r**2-2432.0_cp*a**6*dm**4*dn*r**2-12672.0_cp*a**6*dm**4*r**2+256.0_cp &
      &              *a**6*dm**2*dn**5*r**2+5792.0_cp*a**6*dm**2*dn**4*r**2-32128.0_cp*a**6 &
      &              *dm**2*dn**3*r**2-235328.0_cp*a**6*dm**2*dn**2*r**2+428736.0_cp*a** &
      &              6*dm**2*dn*r**2+2370240.0_cp*a**6*dm**2*r**2-1920.0_cp*a**6*dn**7*r**2 &
      &              +28608.0_cp*a**6*dn**6*r**2-32704.0_cp*a**6*dn**5*r**2-1050592.0_cp &
      &              *a**6*dn**4*r**2+2028544.0_cp*a**6*dn**3*r**2+17063296.0_cp*a**6*dn**2 &
      &              *r**2-18911232.0_cp*a**6*dn*r**2-110490624.0_cp*a**6*r**2-4480.0_cp &
      &              *a**4*b**4*dm**4*dn**3-13440.0_cp*a**4*b**4*dm**4*dn**2+85120.0_cp* &
      &              a**4*b**4*dm**4*dn+255360.0_cp*a**4*b**4*dm**4-26880.0_cp*a**4*b**4 &
      &              *dm**2*dn**5-259392.0_cp*a**4*b**4*dm**2*dn**4+2205952.0_cp*a**4*b**4* &
      &              dm**2*dn**3+8941632.0_cp*a**4*b**4*dm**2*dn**2-26322688.0_cp*a**4*b**4 &
      &              *dm**2*dn-85403136.0_cp*a**4*b**4*dm**2+221760.0_cp*a**4*b**4*dn**7 &
      &              -3326400.0_cp*a**4*b**4*dn**6+6182400.0_cp*a**4*b**4*dn**5+96202176.0_cp &
      &              *a**4*b**4*dn**4-251662656.0_cp*a**4*b**4*dn**3-1233609216.0_cp*a** &
      &              4*b**4*dn**2+1992824064.0_cp*a**4*b**4*dn+6905945088.0_cp*a**4*b**4 &
      &              +3840.0_cp*a**4*b**2*dm**4*dn**3*r**2+11520.0_cp*a**4*b**2*dm**4*dn**2 &
      &              *r**2-72960.0_cp*a**4*b**2*dm**4*dn*r**2-218880.0_cp*a**4*b**2*dm** &
      &              4*r**2+14336.0_cp*a**4*b**2*dm**2*dn**5*r**2+139552.0_cp*a**4*b**2*dm**2 &
      &              *dn**4*r**2-1203584.0_cp*a**4*b**2*dm**2*dn**3*r**2-4865824.0_cp*a**4* &
      &              b**2*dm**2*dn**2*r**2+14480640.0_cp*a**4*b**2*dm**2*dn*r**2+46917504.0_cp &
      &              *a**4*b**2*dm**2*r**2-80640.0_cp*a**4*b**2*dn**7*r**2+1201536.0_cp* &
      &              a**4*b**2*dn**6*r**2-2207744.0_cp*a**4*b**2*dn**5*r**2-34644064.0_cp &
      &              *a**4*b**2*dn**4*r**2+90399488.0_cp*a**4*b**2*dn**3*r**2+443489632.0_cp &
      &              *a**4*b**2*dn**2*r**2-717072768.0_cp*a**4*b**2*dn*r**2-2484430848.0_cp &
      &              *a**4*b**2*r**2-384.0_cp*a**4*dm**4*dn**3*r**4-1152.0_cp*a**4*dm**4 &
      &              *dn**2*r**4+7296.0_cp*a**4*dm**4*dn*r**4+21888.0_cp*a**4*dm**4*r**4 &
      &              -768.0_cp*a**4*dm**2*dn**5*r**4-7584.0_cp*a**4*dm**2*dn**4*r**4+67072.0_cp &
      &              *a**4*dm**2*dn**3*r**4+269856.0_cp*a**4*dm**2*dn**2*r**4-818560.0_cp &
      &              *a**4*dm**2*dn*r**4-2645760.0_cp*a**4*dm**2*r**4+2688.0_cp*a**4*dn**7* &
      &              r**4-39648.0_cp*a**4*dn**6*r**4+71424.0_cp*a**4*dn**5*r**4+1135840.0_cp &
      &              *a**4*dn**4*r**4-2932224.0_cp*a**4*dn**3*r**4-14431360.0_cp*a**4*dn**2 &
      &              *r**4+23162880.0_cp*a**4*dn*r**4+80335872.0_cp*a**4*r**4-3584.0_cp* &
      &              a**2*b**6*dm**4*dn**3+68096.0_cp*a**2*b**6*dm**4*dn+107520.0_cp*a** &
      &              2*b**6*dm**4-53760.0_cp*a**2*b**6*dm**2*dn**5-107520.0_cp*a**2*b**6 &
      &              *dm**2*dn**4+2713088.0_cp*a**2*b**6*dm**2*dn**3+4042752.0_cp*a**2*b**6 &
      &              *dm**2*dn**2-26980352.0_cp*a**2*b**6*dm**2*dn-48427008.0_cp*a**2*b**6* &
      &              dm**2+354816.0_cp*a**2*b**6*dn**7-5322240.0_cp*a**2*b**6*dn**6+13606656.0_cp &
      &              *a**2*b**6*dn**5+111659520.0_cp*a**2*b**6*dn**4-392324352.0_cp*a**2 &
      &              *b**6*dn**3-1017762816.0_cp*a**2*b**6*dn**2+2482357248.0_cp*a**2*b**6* &
      &              dn+4932071424.0_cp*a**2*b**6+7680.0_cp*a**2*b**4*dm**4*dn**3*r**2-145920.0_cp &
      &              *a**2*b**4*dm**4*dn*r**2-230400.0_cp*a**2*b**4*dm**4*r**2+71680.0_cp &
      &              *a**2*b**4*dm**2*dn**5*r**2+143360.0_cp*a**2*b**4*dm**2*dn**4*r**2-3681280.0_cp &
      &              *a**2*b**4*dm**2*dn**3*r**2-5438720.0_cp*a**2*b**4*dm**2*dn**2*r**2 &
      &              +36944640.0_cp*a**2*b**4*dm**2*dn*r**2+66193920.0_cp*a**2*b**4*dm** &
      &              2*r**2-322560.0_cp*a**2*b**4*dn**7*r**2+4806144.0_cp*a**2*b**4*dn** &
      &              6*r**2-12167680.0_cp*a**2*b**4*dn**5*r**2-100653056.0_cp*a**2*b**4*dn**4 &
      &              *r**2+352180480.0_cp*a**2*b**4*dn**3*r**2+915036416.0_cp*a**2*b**4*dn**2 &
      &              *r**2-2231078400.0_cp*a**2*b**4*dn*r**2-4432656384.0_cp*a**2*b**4*r**2 &
      &              -4608.0_cp*a**2*b**2*dm**4*dn**3*r**4+87552.0_cp*a**2*b**2*dm**4*dn &
      &              *r**4+138240.0_cp*a**2*b**2*dm**4*r**4-23040.0_cp*a**2*b**2*dm**2*dn**5 &
      &              *r**4-46080.0_cp*a**2*b**2*dm**2*dn**4*r**4+1219584.0_cp*a**2*b**2*dm**2 &
      &              *dn**3*r**4+1774080.0_cp*a**2*b**2*dm**2*dn**2*r**4-12435456.0_cp*a**2 &
      &              *b**2*dm**2*dn*r**4-22210560.0_cp*a**2*b**2*dm**2*r**4+64512.0_cp*a**2 &
      &              *b**2*dn**7*r**4-951552.0_cp*a**2*b**2*dn**6*r**4+2368512.0_cp*a**2 &
      &              *b**2*dn**5*r**4+19865088.0_cp*a**2*b**2*dn**4*r**4-68680704.0_cp*a**2 &
      &              *b**2*dn**3*r**4-179467008.0_cp*a**2*b**2*dn**2*r**4+433460736.0_cp &
      &              *a**2*b**2*dn*r**4+862451712.0_cp*a**2*b**2*r**4+512.0_cp*a**2*dm** &
      &              4*dn**3*r**6-9728.0_cp*a**2*dm**4*dn*r**6-15360.0_cp*a**2*dm**4*r** &
      &              6+1024.0_cp*a**2*dm**2*dn**5*r**6+2048.0_cp*a**2*dm**2*dn**4*r**6-57856.0_cp &
      &              *a**2*dm**2*dn**3*r**6-81152.0_cp*a**2*dm**2*dn**2*r**6+610560.0_cp &
      &              *a**2*dm**2*dn*r**6+1082880.0_cp*a**2*dm**2*r**6-1536.0_cp*a**2*dn**7* &
      &              r**6+22272.0_cp*a**2*dn**6*r**6-53504.0_cp*a**2*dn**5*r**6-461824.0_cp &
      &              *a**2*dn**4*r**6+1544192.0_cp*a**2*dn**3*r**6+4102144.0_cp*a**2*dn**2* &
      &              r**6-9584640.0_cp*a**2*dn*r**6-19169280.0_cp*a**2*r**6+256.0_cp*b** &
      &              8*dm**4*dn**3+1536.0_cp*b**8*dm**4*dn**2+2816.0_cp*b**8*dm**4*dn+1536.0_cp &
      &              *b**8*dm**4-23040.0_cp*b**8*dm**2*dn**5+71424.0_cp*b**8*dm**2*dn**4 &
      &              +533504.0_cp*b**8*dm**2*dn**3-658176.0_cp*b**8*dm**2*dn**2-3923456.0_cp &
      &              *b**8*dm**2*dn-2826240.0_cp*b**8*dm**2+126720.0_cp*b**8*dn**7-1900800.0_cp &
      &              *b**8*dn**6+6186240.0_cp*b**8*dn**5+24784128.0_cp*b**8*dn**4-124406784.0_cp &
      &              *b**8*dn**3-114729984.0_cp*b**8*dn**2+570802176.0_cp*b**8*dn+544555008.0_cp &
      &              *b**8-1024.0_cp*b**6*dm**4*dn**3*r**2-6144.0_cp*b**6*dm**4*dn**2*r**2- &
      &              11264.0_cp*b**6*dm**4*dn*r**2-6144.0_cp*b**6*dm**4*r**2+57344.0_cp* &
      &              b**6*dm**2*dn**5*r**2-180992.0_cp*b**6*dm**2*dn**4*r**2-1332224.0_cp &
      &              *b**6*dm**2*dn**3*r**2+1692416.0_cp*b**6*dm**2*dn**2*r**2+9910272.0_cp &
      &              *b**6*dm**2*dn*r**2+7123968.0_cp*b**6*dm**2*r**2-215040.0_cp*b**6*dn**7 &
      &              *r**2+3204096.0_cp*b**6*dn**6*r**2-10336256.0_cp*b**6*dn**5*r**2-41819904.0_cp &
      &              *b**6*dn**4*r**2+208361984.0_cp*b**6*dn**3*r**2+193244928.0_cp*b**6 &
      &              *dn**2*r**2-956505600.0_cp*b**6*dn*r**2-913324032.0_cp*b**6*r**2+1536.0_cp &
      &              *b**4*dm**4*dn**3*r**4+9216.0_cp*b**4*dm**4*dn**2*r**4+16896.0_cp*b**4 &
      &              *dm**4*dn*r**4+9216.0_cp*b**4*dm**4*r**4-46080.0_cp*b**4*dm**2*dn** &
      &              5*r**4+149760.0_cp*b**4*dm**2*dn**4*r**4+1075712.0_cp*b**4*dm**2*dn**3 &
      &              *r**4-1436928.0_cp*b**4*dm**2*dn**2*r**4-8165888.0_cp*b**4*dm**2*dn &
      &              *r**4-5849088.0_cp*b**4*dm**2*r**4+107520.0_cp*b**4*dn**7*r**4-1585920.0_cp &
      &              *b**4*dn**6*r**4+5038080.0_cp*b**4*dn**5*r**4+20783360.0_cp*b**4*dn**4 &
      &              *r**4-101876736.0_cp*b**4*dn**3*r**4-96190976.0_cp*b**4*dn**2*r**4+466566144.0_cp &
      &              *b**4*dn*r**4+446828544.0_cp*b**4*r**4-1024.0_cp*b**2*dm**4*dn**3*r**6 &
      &              -6144.0_cp*b**2*dm**4*dn**2*r**6-11264.0_cp*b**2*dm**4*dn*r**6-6144.0_cp &
      &              *b**2*dm**4*r**6+12288.0_cp*b**2*dm**2*dn**5*r**6-42240.0_cp*b**2*dm**2 &
      &              *dn**4*r**6-288768.0_cp*b**2*dm**2*dn**3*r**6+429312.0_cp*b**2*dm** &
      &              2*dn**2*r**6+2294784.0_cp*b**2*dm**2*dn*r**6+1631232.0_cp*b**2*dm** &
      &              2*r**6-15360.0_cp*b**2*dn**7*r**6+222720.0_cp*b**2*dn**6*r**6-685056.0_cp &
      &              *b**2*dn**5*r**6-2957568.0_cp*b**2*dn**4*r**6+13902336.0_cp*b**2*dn**3 &
      &              *r**6+13821696.0_cp*b**2*dn**2*r**6-62995968.0_cp*b**2*dn*r**6-60880896.0_cp &
      &              *b**2*r**6+256.0_cp*dm**4*dn**3*r**8+1536.0_cp*dm**4*dn**2*r**8+2816.0_cp &
      &              *dm**4*dn*r**8+1536.0_cp*dm**4*r**8-512.0_cp*dm**2*dn**5*r**8+2048.0_cp &
      &              *dm**2*dn**4*r**8+11776.0_cp*dm**2*dn**3*r**8-26624.0_cp*dm**2*dn** &
      &              2*r**8-115712.0_cp*dm**2*dn*r**8-79872.0_cp*dm**2*r**8+256.0_cp*dn**7* &
      &              r**8-3584.0_cp*dn**6*r**8+9984.0_cp*dn**5*r**8+49664.0_cp*dn**4*r** &
      &              8-203776.0_cp*dn**3*r**8-239616.0_cp*dn**2*r**8+884736.0_cp*dn*r**8 &
      &              +884736.0_cp*r**8)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+6)=-a**5*b*(88.0_cp*a**6*dm**4*dn**3+240.0_cp*a**6*dm**4 &
      &              *dn**2-2872.0_cp*a**6*dm**4*dn-9360.0_cp*a**6*dm**4+420.0_cp*a**6*dm**2 &
      &              *dn**5+3345.0_cp*a**6*dm**2*dn**4-45691.0_cp*a**6*dm**2*dn**3-167295.0_cp &
      &              *a**6*dm**2*dn**2+778639.0_cp*a**6*dm**2*dn+2624310.0_cp*a**6*dm**2 &
      &              -1980.0_cp*a**6*dn**7+37125.0_cp*a**6*dn**6-117390.0_cp*a**6*dn**5-1229160.0_cp &
      &              *a**6*dn**4+4964682.0_cp*a**6*dn**3+21161055.0_cp*a**6*dn**2-51851088.0_cp &
      &              *a**6*dn-180678060.0_cp*a**6+896.0_cp*a**4*b**2*dm**4*dn**3-30464.0_cp &
      &              *a**4*b**2*dm**4*dn-67200.0_cp*a**4*b**2*dm**4+7680.0_cp*a**4*b**2*dm**2 &
      &              *dn**5+19200.0_cp*a**4*b**2*dm**2*dn**4-603392.0_cp*a**4*b**2*dm**2 &
      &              *dn**3-1076160.0_cp*a**4*b**2*dm**2*dn**2+9131648.0_cp*a**4*b**2*dm**2 &
      &              *dn+21099840.0_cp*a**4*b**2*dm**2-31680.0_cp*a**4*b**2*dn**7+594000.0_cp &
      &              *a**4*b**2*dn**6-2209920.0_cp*a**4*b**2*dn**5-14949600.0_cp*a**4*b**2* &
      &              dn**4+75004224.0_cp*a**4*b**2*dn**3+200982480.0_cp*a**4*b**2*dn**2-668088576.0_cp &
      &              *a**4*b**2*dn-1616397120.0_cp*a**4*b**2-384.0_cp*a**4*dm**4*dn**3*r**2 &
      &              +13056.0_cp*a**4*dm**4*dn*r**2+28800.0_cp*a**4*dm**4*r**2-2048.0_cp &
      &              *a**4*dm**2*dn**5*r**2-5120.0_cp*a**4*dm**2*dn**4*r**2+164096.0_cp* &
      &              a**4*dm**2*dn**3*r**2+290000.0_cp*a**4*dm**2*dn**2*r**2-2513376.0_cp &
      &              *a**4*dm**2*dn*r**2-5802480.0_cp*a**4*dm**2*r**2+5760.0_cp*a**4*dn**7* &
      &              r**2-107280.0_cp*a**4*dn**6*r**2+396032.0_cp*a**4*dn**5*r**2+2693600.0_cp &
      &              *a**4*dn**4*r**2-13469504.0_cp*a**4*dn**3*r**2-36108800.0_cp*a**4*dn**2 &
      &              *r**2+120220704.0_cp*a**4*dn*r**2+290823120.0_cp*a**4*r**2+896.0_cp &
      &              *a**2*b**4*dm**4*dn**3-5376.0_cp*a**2*b**4*dm**4*dn**2-43904.0_cp*a**2 &
      &              *b**4*dm**4*dn-59136.0_cp*a**2*b**4*dm**4+24192.0_cp*a**2*b**4*dm** &
      &              2*dn**5-42336.0_cp*a**2*b**4*dm**2*dn**4-1231328.0_cp*a**2*b**4*dm**2* &
      &              dn**3+100128.0_cp*a**2*b**4*dm**2*dn**2+16039520.0_cp*a**2*b**4*dm**2* &
      &              dn+23279424.0_cp*a**2*b**4*dm**2-88704.0_cp*a**2*b**4*dn**7+1663200.0_cp &
      &              *a**2*b**4*dn**6-7116480.0_cp*a**2*b**4*dn**5-28651392.0_cp*a**2*b**4* &
      &              dn**4+189193536.0_cp*a**2*b**4*dn**3+258578208.0_cp*a**2*b**4*dn**2 &
      &              -1373879808.0_cp*a**2*b**4*dn-2155628160.0_cp*a**2*b**4-1280.0_cp*a**2 &
      &              *b**2*dm**4*dn**3*r**2+7680.0_cp*a**2*b**2*dm**4*dn**2*r**2+62720.0_cp &
      &              *a**2*b**2*dm**4*dn*r**2+84480.0_cp*a**2*b**2*dm**4*r**2-21504.0_cp &
      &              *a**2*b**2*dm**2*dn**5*r**2+38640.0_cp*a**2*b**2*dm**2*dn**4*r**2+1110192.0_cp &
      &              *a**2*b**2*dm**2*dn**3*r**2-123600.0_cp*a**2*b**2*dm**2*dn**2*r**2-14662704.0_cp &
      &              *a**2*b**2*dm**2*dn*r**2-21255840.0_cp*a**2*b**2*dm**2*r**2+53760.0_cp &
      &              *a**2*b**2*dn**7*r**2-1001280.0_cp*a**2*b**2*dn**6*r**2+4252416.0_cp &
      &              *a**2*b**2*dn**5*r**2+17239600.0_cp*a**2*b**2*dn**4*r**2-113215248.0_cp &
      &              *a**2*b**2*dn**3*r**2-154887760.0_cp*a**2*b**2*dn**2*r**2+823422096.0_cp &
      &              *a**2*b**2*dn*r**2+1291880160.0_cp*a**2*b**2*r**2+384.0_cp*a**2*dm**4* &
      &              dn**3*r**4-2304.0_cp*a**2*dm**4*dn**2*r**4-18816.0_cp*a**2*dm**4*dn &
      &              *r**4-25344.0_cp*a**2*dm**4*r**4+3456.0_cp*a**2*dm**2*dn**5*r**4-6480.0_cp &
      &              *a**2*dm**2*dn**4*r**4-182800.0_cp*a**2*dm**2*dn**3*r**4+30192.0_cp &
      &              *a**2*dm**2*dn**2*r**4+2473744.0_cp*a**2*dm**2*dn*r**4+3578592.0_cp &
      &              *a**2*dm**2*r**4-5376.0_cp*a**2*dn**7*r**4+99120.0_cp*a**2*dn**6*r**4- &
      &              415488.0_cp*a**2*dn**5*r**4-1707920.0_cp*a**2*dn**4*r**4+11047248.0_cp &
      &              *a**2*dn**3*r**4+15267584.0_cp*a**2*dn**2*r**4-79922640.0_cp*a**2*dn &
      &              *r**4-125538336.0_cp*a**2*r**4-512.0_cp*b**6*dm**4*dn**3-3072.0_cp* &
      &              b**6*dm**4*dn**2-5632.0_cp*b**6*dm**4*dn-3072.0_cp*b**6*dm**4+15360.0_cp &
      &              *b**6*dm**2*dn**5-79104.0_cp*b**6*dm**2*dn**4-384256.0_cp*b**6*dm** &
      &              2*dn**3+1054464.0_cp*b**6*dm**2*dn**2+4190464.0_cp*b**6*dm**2*dn+2846208.0_cp &
      &              *b**6*dm**2-50688.0_cp*b**6*dn**7+950400.0_cp*b**6*dn**6-4597248.0_cp &
      &              *b**6*dn**5-8825088.0_cp*b**6*dn**4+91407360.0_cp*b**6*dn**3+20284032.0_cp &
      &              *b**6*dn**2-502548480.0_cp*b**6*dn-428198400.0_cp*b**6+1536.0_cp*b**4* &
      &              dm**4*dn**3*r**2+9216.0_cp*b**4*dm**4*dn**2*r**2+16896.0_cp*b**4*dm**4 &
      &              *dn*r**2+9216.0_cp*b**4*dm**4*r**2-28672.0_cp*b**4*dm**2*dn**5*r**2 &
      &              +150080.0_cp*b**4*dm**2*dn**4*r**2+716608.0_cp*b**4*dm**2*dn**3*r** &
      &              2-2032832.0_cp*b**4*dm**2*dn**2*r**2-7974720.0_cp*b**4*dm**2*dn*r** &
      &              2-5404032.0_cp*b**4*dm**2*r**2+64512.0_cp*b**4*dn**7*r**2-1201536.0_cp &
      &              *b**4*dn**6*r**2+5770240.0_cp*b**4*dn**5*r**2+11206720.0_cp*b**4*dn**4 &
      &              *r**2-114814144.0_cp*b**4*dn**3*r**2-25750336.0_cp*b**4*dn**2*r**2+631781568.0_cp &
      &              *b**4*dn*r**2+538547328.0_cp*b**4*r**2-1536.0_cp*b**2*dm**4*dn**3*r**4 &
      &              -9216.0_cp*b**2*dm**4*dn**2*r**4-16896.0_cp*b**2*dm**4*dn*r**4-9216.0_cp &
      &              *b**2*dm**4*r**4+15360.0_cp*b**2*dm**2*dn**5*r**4-82560.0_cp*b**2*dm**2 &
      &              *dn**4*r**4-382592.0_cp*b**2*dm**2*dn**3*r**4+1150848.0_cp*b**2*dm**2* &
      &              dn**2*r**4+4416128.0_cp*b**2*dm**2*dn*r**4+2980608.0_cp*b**2*dm**2* &
      &              r**4-21504.0_cp*b**2*dn**7*r**4+396480.0_cp*b**2*dn**6*r**4-1880064.0_cp &
      &              *b**2*dn**5*r**4-3750400.0_cp*b**2*dn**4*r**4+37438080.0_cp*b**2*dn**3 &
      &              *r**4+9069376.0_cp*b**2*dn**2*r**4-205123200.0_cp*b**2*dn*r**4-175302144.0_cp &
      &              *b**2*r**4+512.0_cp*dm**4*dn**3*r**6+3072.0_cp*dm**4*dn**2*r**6+5632.0_cp &
      &              *dm**4*dn*r**6+3072.0_cp*dm**4*r**6-2048.0_cp*dm**2*dn**5*r**6+11584.0_cp &
      &              *dm**2*dn**4*r**6+50240.0_cp*dm**2*dn**3*r**6-172480.0_cp*dm**2*dn**2* &
      &              r**6-631872.0_cp*dm**2*dn*r**6-422784.0_cp*dm**2*r**6+1536.0_cp*dn**7* &
      &              r**6-27840.0_cp*dn**6*r**6+128512.0_cp*dn**5*r**6+274624.0_cp*dn**4 &
      &              *r**6-2556352.0_cp*dn**3*r**6-774016.0_cp*dn**2*r**6+13758912.0_cp*dn &
      &              *r**6+11859840.0_cp*r**6)/(2048.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn &
      &              -1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+7)=-a**6*(6.0_cp*a**6*dm**4*dn**3-6.0_cp*a**6*dm**4*dn** &
      &              2-324.0_cp*a**6*dm**4*dn-816.0_cp*a**6*dm**4+40.0_cp*a**6*dm**2*dn**5+ &
      &              69.0_cp*a**6*dm**2*dn**4-4138.0_cp*a**6*dm**2*dn**3-6627.0_cp*a**6*dm**2 &
      &              *dn**2+86112.0_cp*a**6*dm**2*dn+221964.0_cp*a**6*dm**2-110.0_cp*a** &
      &              6*dn**7+2475.0_cp*a**6*dn**6-12280.0_cp*a**6*dn**5-68407.0_cp*a**6*dn**4 &
      &              +481422.0_cp*a**6*dn**3+1221508.0_cp*a**6*dn**2-5610168.0_cp*a**6*dn &
      &              -14836320.0_cp*a**6+168.0_cp*a**4*b**2*dm**4*dn**3-1008.0_cp*a**4*b**2 &
      &              *dm**4*dn**2-11592.0_cp*a**4*b**2*dm**4*dn-21168.0_cp*a**4*b**2*dm**4+ &
      &              2340.0_cp*a**4*b**2*dm**2*dn**5-3996.0_cp*a**4*b**2*dm**2*dn**4-177972.0_cp &
      &              *a**4*b**2*dm**2*dn**3-19188.0_cp*a**4*b**2*dm**2*dn**2+3357648.0_cp &
      &              *a**4*b**2*dm**2*dn+6332688.0_cp*a**4*b**2*dm**2-5940.0_cp*a**4*b** &
      &              2*dn**7+133650.0_cp*a**4*b**2*dn**6-725310.0_cp*a**4*b**2*dn**5-2632662.0_cp &
      &              *a**4*b**2*dn**4+23611050.0_cp*a**4*b**2*dn**3+35194500.0_cp*a**4*b**2 &
      &              *dn**2-237236040.0_cp*a**4*b**2*dn-464386608.0_cp*a**4*b**2-24.0_cp &
      &              *a**4*dm**4*dn**3*r**2+144.0_cp*a**4*dm**4*dn**2*r**2+1656.0_cp*a** &
      &              4*dm**4*dn*r**2+3024.0_cp*a**4*dm**4*r**2-208.0_cp*a**4*dm**2*dn**5 &
      &              *r**2+366.0_cp*a**4*dm**2*dn**4*r**2+16084.0_cp*a**4*dm**2*dn**3*r**2+ &
      &              1146.0_cp*a**4*dm**2*dn**2*r**2-308268.0_cp*a**4*dm**2*dn*r**2-581040.0_cp &
      &              *a**4*dm**2*r**2+360.0_cp*a**4*dn**7*r**2-8046.0_cp*a**4*dn**6*r**2 &
      &              +43372.0_cp*a**4*dn**5*r**2+158334.0_cp*a**4*dn**4*r**2-1413052.0_cp &
      &              *a**4*dn**3*r**2-2104872.0_cp*a**4*dn**2*r**2+14232816.0_cp*a**4*dn &
      &              *r**2+27857088.0_cp*a**4*r**2-6720.0_cp*a**2*b**4*dm**4*dn**2-33600.0_cp &
      &              *a**2*b**4*dm**4*dn-40320.0_cp*a**2*b**4*dm**4+13440.0_cp*a**2*b**4 &
      &              *dm**2*dn**5-62496.0_cp*a**2*b**4*dm**2*dn**4-670656.0_cp*a**2*b**4 &
      &              *dm**2*dn**3+1117536.0_cp*a**2*b**4*dm**2*dn**2+11432064.0_cp*a**2* &
      &              b**4*dm**2*dn+14458752.0_cp*a**2*b**4*dm**2-31680.0_cp*a**2*b**4*dn**7 &
      &              +712800.0_cp*a**2*b**4*dn**6-4200000.0_cp*a**2*b**4*dn**5-8380512.0_cp &
      &              *a**2*b**4*dn**4+110197248.0_cp*a**2*b**4*dn**3+58367232.0_cp*a**2* &
      &              b**4*dn**2-921203712.0_cp*a**2*b**4*dn-1244284416.0_cp*a**2*b**4+5760.0_cp &
      &              *a**2*b**2*dm**4*dn**2*r**2+28800.0_cp*a**2*b**2*dm**4*dn*r**2+34560.0_cp &
      &              *a**2*b**2*dm**4*r**2-7168.0_cp*a**2*b**2*dm**2*dn**5*r**2+33936.0_cp &
      &              *a**2*b**2*dm**2*dn**4*r**2+361312.0_cp*a**2*b**2*dm**2*dn**3*r**2-626352.0_cp &
      &              *a**2*b**2*dm**2*dn**2*r**2-6278400.0_cp*a**2*b**2*dm**2*dn*r**2-7933248.0_cp &
      &              *a**2*b**2*dm**2*r**2+11520.0_cp*a**2*b**2*dn**7*r**2-257472.0_cp*a**2 &
      &              *b**2*dn**6*r**2+1507072.0_cp*a**2*b**2*dn**5*r**2+3035088.0_cp*a** &
      &              2*b**2*dn**4*r**2-39556384.0_cp*a**2*b**2*dn**3*r**2-20913264.0_cp* &
      &              a**2*b**2*dn**2*r**2+331354944.0_cp*a**2*b**2*dn*r**2+447529536.0_cp &
      &              *a**2*b**2*r**2-576.0_cp*a**2*dm**4*dn**2*r**4-2880.0_cp*a**2*dm**4 &
      &              *dn*r**4-3456.0_cp*a**2*dm**4*r**4+384.0_cp*a**2*dm**2*dn**5*r**4-1872.0_cp &
      &              *a**2*dm**2*dn**4*r**4-19680.0_cp*a**2*dm**2*dn**3*r**4+36528.0_cp* &
      &              a**2*dm**2*dn**2*r**4+353856.0_cp*a**2*dm**2*dn*r**4+446400.0_cp*a**2* &
      &              dm**2*r**4-384.0_cp*a**2*dn**7*r**4+8496.0_cp*a**2*dn**6*r**4-49152.0_cp &
      &              *a**2*dn**5*r**4-101200.0_cp*a**2*dn**4*r**4+1287264.0_cp*a**2*dn** &
      &              3*r**4+702400.0_cp*a**2*dn**2*r**4-10712448.0_cp*a**2*dn*r**4-14482944.0_cp &
      &              *a**2*r**4-896.0_cp*b**6*dm**4*dn**3-5376.0_cp*b**6*dm**4*dn**2-9856.0_cp &
      &              *b**6*dm**4*dn-5376.0_cp*b**6*dm**4+13440.0_cp*b**6*dm**2*dn**5-96768.0_cp &
      &              *b**6*dm**2*dn**4-332416.0_cp*b**6*dm**2*dn**3+1634304.0_cp*b**6*dm**2 &
      &              *dn**2+5361664.0_cp*b**6*dm**2*dn+3505152.0_cp*b**6*dm**2-29568.0_cp &
      &              *b**6*dn**7+665280.0_cp*b**6*dn**6-4229568.0_cp*b**6*dn**5-2538816.0_cp &
      &              *b**6*dn**4+85366848.0_cp*b**6*dn**3-33801600.0_cp*b**6*dn**2-560270592.0_cp &
      &              *b**6*dn-443487744.0_cp*b**6+1920.0_cp*b**4*dm**4*dn**3*r**2+11520.0_cp &
      &              *b**4*dm**4*dn**2*r**2+21120.0_cp*b**4*dm**4*dn*r**2+11520.0_cp*b** &
      &              4*dm**4*r**2-17920.0_cp*b**4*dm**2*dn**5*r**2+131040.0_cp*b**4*dm** &
      &              2*dn**4*r**2+439360.0_cp*b**4*dm**2*dn**3*r**2-2252640.0_cp*b**4*dm**2 &
      &              *dn**2*r**2-7312320.0_cp*b**4*dm**2*dn*r**2-4769280.0_cp*b**4*dm**2 &
      &              *r**2+26880.0_cp*b**4*dn**7*r**2-600768.0_cp*b**4*dn**6*r**2+3794560.0_cp &
      &              *b**4*dn**5*r**2+2341472.0_cp*b**4*dn**4*r**2-76570240.0_cp*b**4*dn**3 &
      &              *r**2+30349408.0_cp*b**4*dn**2*r**2+503268480.0_cp*b**4*dn*r**2+398429568.0_cp &
      &              *b**4*r**2-1152.0_cp*b**2*dm**4*dn**3*r**4-6912.0_cp*b**2*dm**4*dn**2* &
      &              r**4-12672.0_cp*b**2*dm**4*dn*r**4-6912.0_cp*b**2*dm**4*r**4+5760.0_cp &
      &              *b**2*dm**2*dn**5*r**4-43200.0_cp*b**2*dm**2*dn**4*r**4-138624.0_cp &
      &              *b**2*dm**2*dn**3*r**4+766656.0_cp*b**2*dm**2*dn**2*r**4+2443776.0_cp &
      &              *b**2*dm**2*dn*r**4+1587456.0_cp*b**2*dm**2*r**4-5376.0_cp*b**2*dn**7* &
      &              r**4+118944.0_cp*b**2*dn**6*r**4-742656.0_cp*b**2*dn**5*r**4-492416.0_cp &
      &              *b**2*dn**4*r**4+14986368.0_cp*b**2*dn**3*r**4-5605792.0_cp*b**2*dn**2 &
      &              *r**4-97909632.0_cp*b**2*dn*r**4-77692032.0_cp*b**2*r**4+128.0_cp*dm**4 &
      &              *dn**3*r**6+768.0_cp*dm**4*dn**2*r**6+1408.0_cp*dm**4*dn*r**6+768.0_cp &
      &              *dm**4*r**6-256.0_cp*dm**2*dn**5*r**6+2016.0_cp*dm**2*dn**4*r**6+5824.0_cp &
      &              *dm**2*dn**3*r**6-38496.0_cp*dm**2*dn**2*r**6-118080.0_cp*dm**2*dn* &
      &              r**6-76032.0_cp*dm**2*r**6+128.0_cp*dn**7*r**6-2784.0_cp*dn**6*r**6 &
      &              +16960.0_cp*dn**5*r**6+13536.0_cp*dn**4*r**6-341568.0_cp*dn**3*r**6 &
      &              +100224.0_cp*dn**2*r**6+2177280.0_cp*dn*r**6+1741824.0_cp*r**6)/(2048.0_cp &
      &              *dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn &
      &              +3.0_cp))
      stencil(ku+8)=a**7*b*(-8.0_cp*a**4*dm**4*dn**2+360.0_cp*a**4*dm**4*dn &
      &              +1232.0_cp*a**4*dm**4-380.0_cp*a**4*dm**2*dn**4+3023.0_cp*a**4*dm** &
      &              2*dn**3+18448.0_cp*a**4*dm**2*dn**2-98177.0_cp*a**4*dm**2*dn-354662.0_cp &
      &              *a**4*dm**2+660.0_cp*a**4*dn**6-19305.0_cp*a**4*dn**5+179965.0_cp*a**4 &
      &              *dn**4-323119.0_cp*a**4*dn**3-2953993.0_cp*a**4*dn**2+6458116.0_cp* &
      &              a**4*dn+25083492.0_cp*a**4+224.0_cp*a**2*b**2*dm**4*dn**2+3360.0_cp &
      &              *a**2*b**2*dm**4*dn+5824.0_cp*a**2*b**2*dm**4-4800.0_cp*a**2*b**2*dm**2 &
      &              *dn**4+49008.0_cp*a**2*b**2*dm**2*dn**3+78784.0_cp*a**2*b**2*dm**2*dn**2 &
      &              -1051152.0_cp*a**2*b**2*dm**2*dn-1948576.0_cp*a**2*b**2*dm**2+7920.0_cp &
      &              *a**2*b**2*dn**6-231660.0_cp*a**2*b**2*dn**5+2242500.0_cp*a**2*b**2 &
      &              *dn**4-5777124.0_cp*a**2*b**2*dn**3-23298180.0_cp*a**2*b**2*dn**2+77497296.0_cp &
      &              *a**2*b**2*dn+158170320.0_cp*a**2*b**2-96.0_cp*a**2*dm**4*dn**2*r** &
      &              2-1440.0_cp*a**2*dm**4*dn*r**2-2496.0_cp*a**2*dm**4*r**2+1280.0_cp* &
      &              a**2*dm**2*dn**4*r**2-13220.0_cp*a**2*dm**2*dn**3*r**2-20816.0_cp*a**2 &
      &              *dm**2*dn**2*r**2+289100.0_cp*a**2*dm**2*dn*r**2+535224.0_cp*a**2*dm**2 &
      &              *r**2-1440.0_cp*a**2*dn**6*r**2+41868.0_cp*a**2*dn**5*r**2-403364.0_cp &
      &              *a**2*dn**4*r**2+1035344.0_cp*a**2*dn**3*r**2+4183700.0_cp*a**2*dn**2* &
      &              r**2-13942244.0_cp*a**2*dn*r**2-28450776.0_cp*a**2*r**2+896.0_cp*b**4* &
      &              dm**4*dn**2+2688.0_cp*b**4*dm**4*dn+1792.0_cp*b**4*dm**4-8064.0_cp* &
      &              b**4*dm**2*dn**4+98784.0_cp*b**4*dm**2*dn**3-116480.0_cp*b**4*dm**2 &
      &              *dn**2-1161888.0_cp*b**4*dm**2*dn-938560.0_cp*b**4*dm**2+12672.0_cp &
      &              *b**4*dn**6-370656.0_cp*b**4*dn**5+3720672.0_cp*b**4*dn**4-12282912.0_cp &
      &              *b**4*dn**3-16635744.0_cp*b**4*dn**2+104171904.0_cp*b**4*dn+104420736.0_cp &
      &              *b**4-1280.0_cp*b**2*dm**4*dn**2*r**2-3840.0_cp*b**2*dm**4*dn*r**2-2560.0_cp &
      &              *b**2*dm**4*r**2+7168.0_cp*b**2*dm**2*dn**4*r**2-88816.0_cp*b**2*dm**2 &
      &              *dn**3*r**2+110144.0_cp*b**2*dm**2*dn**2*r**2+1059664.0_cp*b**2*dm**2* &
      &              dn*r**2+853536.0_cp*b**2*dm**2*r**2-7680.0_cp*b**2*dn**6*r**2+223296.0_cp &
      &              *b**2*dn**5*r**2-2230720.0_cp*b**2*dn**4*r**2+7340304.0_cp*b**2*dn**3* &
      &              r**2+9948544.0_cp*b**2*dn**2*r**2-62403312.0_cp*b**2*dn*r**2-62549856.0_cp &
      &              *b**2*r**2+384.0_cp*dm**4*dn**2*r**4+1152.0_cp*dm**4*dn*r**4+768.0_cp &
      &              *dm**4*r**4-1152.0_cp*dm**2*dn**4*r**4+14544.0_cp*dm**2*dn**3*r**4-19648.0_cp &
      &              *dm**2*dn**2*r**4-178032.0_cp*dm**2*dn*r**4-142688.0_cp*dm**2*r**4+768.0_cp &
      &              *dn**6*r**4-22128.0_cp*dn**5*r**4+219216.0_cp*dn**4*r**4-714208.0_cp &
      &              *dn**3*r**4-987600.0_cp*dn**2*r**4+6059248.0_cp*dn*r**4+6090528.0_cp &
      &              *r**4)/(2048.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp))
      stencil(ku+9)=a**8*(a**4*dm**4*dn**2+33.0_cp*a**4*dm**4*dn+92.0_cp* &
      &              a**4*dm**4-26.0_cp*a**4*dm**2*dn**4+280.0_cp*a**4*dm**2*dn**3+904.0_cp &
      &              *a**4*dm**2*dn**2-8842.0_cp*a**4*dm**2*dn-25664.0_cp*a**4*dm**2+33.0_cp &
      &              *a**4*dn**6-1089.0_cp*a**4*dn**5+11788.0_cp*a**4*dn**4-31790.0_cp*a**4 &
      &              *dn**3-176752.0_cp*a**4*dn**2+572264.0_cp*a**4*dn+1760448.0_cp*a**4 &
      &              +112.0_cp*a**2*b**2*dm**4*dn**2+1008.0_cp*a**2*b**2*dm**4*dn+1568.0_cp &
      &              *a**2*b**2*dm**4-1080.0_cp*a**2*b**2*dm**2*dn**4+13608.0_cp*a**2*b**2* &
      &              dm**2*dn**3+4736.0_cp*a**2*b**2*dm**2*dn**2-302112.0_cp*a**2*b**2*dm**2 &
      &              *dn-497024.0_cp*a**2*b**2*dm**2+1320.0_cp*a**2*b**2*dn**6-43560.0_cp &
      &              *a**2*b**2*dn**5+485340.0_cp*a**2*b**2*dn**4-1627524.0_cp*a**2*b**2 &
      &              *dn**3-4433232.0_cp*a**2*b**2*dn**2+21494160.0_cp*a**2*b**2*dn+38457216.0_cp &
      &              *a**2*b**2-16.0_cp*a**2*dm**4*dn**2*r**2-144.0_cp*a**2*dm**4*dn*r** &
      &              2-224.0_cp*a**2*dm**4*r**2+96.0_cp*a**2*dm**2*dn**4*r**2-1224.0_cp* &
      &              a**2*dm**2*dn**3*r**2-360.0_cp*a**2*dm**2*dn**2*r**2+27720.0_cp*a** &
      &              2*dm**2*dn*r**2+45552.0_cp*a**2*dm**2*r**2-80.0_cp*a**2*dn**6*r**2+2624.0_cp &
      &              *a**2*dn**5*r**2-29096.0_cp*a**2*dn**4*r**2+97264.0_cp*a**2*dn**3*r**2 &
      &              +265120.0_cp*a**2*dn**2*r**2-1289280.0_cp*a**2*dn*r**2-2306304.0_cp &
      &              *a**2*r**2+560.0_cp*b**4*dm**4*dn**2+1680.0_cp*b**4*dm**4*dn+1120.0_cp &
      &              *b**4*dm**4-3360.0_cp*b**4*dm**2*dn**4+48048.0_cp*b**4*dm**2*dn**3-84560.0_cp &
      &              *b**4*dm**2*dn**2-640416.0_cp*b**4*dm**2*dn-504448.0_cp*b**4*dm**2+3960.0_cp &
      &              *b**4*dn**6-130680.0_cp*b**4*dn**5+1497480.0_cp*b**4*dn**4-5950344.0_cp &
      &              *b**4*dn**3-5013600.0_cp*b**4*dn**2+53374368.0_cp*b**4*dn+50805504.0_cp &
      &              *b**4-480.0_cp*b**2*dm**4*dn**2*r**2-1440.0_cp*b**2*dm**4*dn*r**2-960.0_cp &
      &              *b**2*dm**4*r**2+1792.0_cp*b**2*dm**2*dn**4*r**2-25928.0_cp*b**2*dm**2 &
      &              *dn**3*r**2+47576.0_cp*b**2*dm**2*dn**2*r**2+351104.0_cp*b**2*dm**2 &
      &              *dn*r**2+275808.0_cp*b**2*dm**2*r**2-1440.0_cp*b**2*dn**6*r**2+47232.0_cp &
      &              *b**2*dn**5*r**2-538624.0_cp*b**2*dn**4*r**2+2134040.0_cp*b**2*dn** &
      &              3*r**2+1793896.0_cp*b**2*dn**2*r**2-19190864.0_cp*b**2*dn*r**2-18263424.0_cp &
      &              *b**2*r**2+48.0_cp*dm**4*dn**2*r**4+144.0_cp*dm**4*dn*r**4+96.0_cp*dm**4 &
      &              *r**4-96.0_cp*dm**2*dn**4*r**4+1416.0_cp*dm**2*dn**3*r**4-2792.0_cp &
      &              *dm**2*dn**2*r**4-19728.0_cp*dm**2*dn*r**4-15424.0_cp*dm**2*r**4+48.0_cp &
      &              *dn**6*r**4-1560.0_cp*dn**5*r**4+17640.0_cp*dn**4*r**4-69232.0_cp*dn**3 &
      &              *r**4-60000.0_cp*dn**2*r**4+620608.0_cp*dn*r**4+592128.0_cp*r**4)/(2048.0_cp &
      &              *dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+10)=-a**9*b*(-24.0_cp*a**2*dm**4*dn-120.0_cp*a**2*dm**4+140.0_cp &
      &              *a**2*dm**2*dn**3-2365.0_cp*a**2*dm**2*dn**2+6048.0_cp*a**2*dm**2*dn &
      &              +36417.0_cp*a**2*dm**2-132.0_cp*a**2*dn**5+5115.0_cp*a**2*dn**4-71305.0_cp &
      &              *a**2*dn**3+392435.0_cp*a**2*dn**2-318219.0_cp*a**2*dn-2706966.0_cp &
      &              *a**2-224.0_cp*b**2*dm**4*dn-224.0_cp*b**2*dm**4+960.0_cp*b**2*dm** &
      &              2*dn**3-17616.0_cp*b**2*dm**2*dn**2+71744.0_cp*b**2*dm**2*dn+90320.0_cp &
      &              *b**2*dm**2-880.0_cp*b**2*dn**5+34100.0_cp*b**2*dn**4-484580.0_cp*b**2 &
      &              *dn**3+2898148.0_cp*b**2*dn**2-4984044.0_cp*b**2*dn-8401752.0_cp*b**2+ &
      &              96.0_cp*dm**4*dn*r**2+96.0_cp*dm**4*r**2-256.0_cp*dm**2*dn**3*r**2+4748.0_cp &
      &              *dm**2*dn**2*r**2-19728.0_cp*dm**2*dn*r**2-24732.0_cp*dm**2*r**2+160.0_cp &
      &              *dn**5*r**2-6164.0_cp*dn**4*r**2+87204.0_cp*dn**3*r**2-520360.0_cp*dn**2 &
      &              *r**2+896540.0_cp*dn*r**2+1510428.0_cp*r**2)/(2048.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+11)=-a**10*(-2.0_cp*a**2*dm**4*dn-8.0_cp*a**2*dm**4+8.0_cp &
      &              *a**2*dm**2*dn**3-153.0_cp*a**2*dm**2*dn**2+507.0_cp*a**2*dm**2*dn+2342.0_cp &
      &              *a**2*dm**2-6.0_cp*a**2*dn**5+255.0_cp*a**2*dn**4-3919.0_cp*a**2*dn**3 &
      &              +24204.0_cp*a**2*dn**2-28316.0_cp*a**2*dn-168240.0_cp*a**2-56.0_cp* &
      &              b**2*dm**4*dn-56.0_cp*b**2*dm**4+180.0_cp*b**2*dm**2*dn**3-3672.0_cp &
      &              *b**2*dm**2*dn**2+16892.0_cp*b**2*dm**2*dn+20744.0_cp*b**2*dm**2-132.0_cp &
      &              *b**2*dn**5+5610.0_cp*b**2*dn**4-87600.0_cp*b**2*dn**3+578706.0_cp* &
      &              b**2*dn**2-1136232.0_cp*b**2*dn-1808280.0_cp*b**2+8.0_cp*dm**4*dn*r**2 &
      &              +8.0_cp*dm**4*r**2-16.0_cp*dm**2*dn**3*r**2+330.0_cp*dm**2*dn**2*r**2- &
      &              1550.0_cp*dm**2*dn*r**2-1896.0_cp*dm**2*r**2+8.0_cp*dn**5*r**2-338.0_cp &
      &              *dn**4*r**2+5254.0_cp*dn**3*r**2-34632.0_cp*dn**2*r**2+68152.0_cp*dn &
      &              *r**2+108384.0_cp*r**2)/(2048.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp))
      stencil(ku+12)=a**11*b*(8.0_cp*dm**4-20.0_cp*dm**2*dn**2+469.0_cp*dm**2 &
      &              *dn-2771.0_cp*dm**2+12.0_cp*dn**4-567.0_cp*dn**3+10012.0_cp*dn**2-78299.0_cp &
      &              *dn+228822.0_cp)/(2048.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp))
      stencil(ku+13)=a**12*(dm-dn+14.0_cp)*(dm+dn-14.0_cp)*(dm**2-dn**2+23.0_cp &
      &              *dn-132.0_cp)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp))
      stencil(ku+14:len_stencil) = 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4rmult4hmult8laplrot2
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
