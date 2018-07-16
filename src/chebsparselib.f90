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

   public :: rmult1, intcheb1, intcheb2, intcheb4, intcheb2rmult1, eye,     &
   &         rmult2, intcheb2rmult2, intcheb1rmult1, intcheb2rmult2lapl,    &
   &         intcheb4rmult4lapl2, intcheb4rmult4lapl, intcheb4rmult4,       &
   &         intcheb4rmult4laplrot, intcheb4rmult4laplrot2,                 &
   &         intcheb2rmult2laplrot, intcheb4rmult3

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

      call mirror_stencil(idx, stencil, len_stencil)

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

      call mirror_stencil(idx, stencil, len_stencil)

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

      call mirror_stencil(idx, stencil, len_stencil)

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

      call mirror_stencil(idx, stencil, len_stencil)

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

      call mirror_stencil(idx, stencil, len_stencil)

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

      call mirror_stencil(idx, stencil, len_stencil)

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

      call mirror_stencil(idx, stencil, len_stencil)

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

      call mirror_stencil(idx, stencil, len_stencil)

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

      call mirror_stencil(idx, stencil, len_stencil)

   end function intcheb2rmult2lapl
!------------------------------------------------------------------------------
   function intcheb2rmult2laplrot(a, b, m, n, len_stencil) result(stencil)

      !-- Input variables
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
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

      stencil(1:ku-4) = 0.0_cp
      stencil(ku-3)=a**4*(dm**2-dn**2-3.0_cp*dn-2.0_cp)/(16.0_cp*dn*(dn+1.0_cp &
      &             ))
      stencil(ku-2)=a**3*b*(dm**2-2.0_cp*dn**2-5.0_cp*dn-3.0_cp)/(4.0_cp*dn &
      &             *(dn+1.0_cp))
      stencil(ku-1)=-a**2*(2.0_cp*a**2*dm**2*dn-a**2*dm**2-4.0_cp*a**2*dn**2 &
      &             -2.0_cp*a**2*dn+6.0_cp*a**2+4.0_cp*a*b*dm**2*dn-4.0_cp*a*b*dm**2-4.0_cp &
      &             *a*b*dn**3-12.0_cp*a*b*dn**2+16.0_cp*a*b+10.0_cp*b**2*dn**3+6.0_cp* &
      &             b**2*dn**2-12.0_cp*b**2*dn-4.0_cp*b**2)/(8.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp &
      &             ))
      stencil(ku)=-a*b*(a**2*dm**2+2.0_cp*a**2*dn**2-5.0_cp*a**2*dn+3.0_cp &
      &             *a**2-8.0_cp*a*b*dn**2-4.0_cp*a*b*dn+12.0_cp*a*b+4.0_cp*b**2*dn**2-6.0_cp &
      &             *b**2*dn+2.0_cp*b**2)/(4.0_cp*dn*(dn-1.0_cp))
      stencil(ku+1)=a*(3.0_cp*a**3*dm**2+a**3*dn**2-5.0_cp*a**3+8.0_cp*a**2* &
      &             b*dm**2+8.0_cp*a**2*b*dn**2-16.0_cp*a**2*b-12.0_cp*a*b**2*dn**2+8.0_cp &
      &             *a*b**2+16.0_cp*b**3*dn**2-16.0_cp*b**3)/(8.0_cp*(dn-1.0_cp)*(dn+1.0_cp &
      &             ))
      stencil(ku+2)=-a*b*(a**2*dm**2+2.0_cp*a**2*dn**2+5.0_cp*a**2*dn+3.0_cp &
      &             *a**2-8.0_cp*a*b*dn**2+4.0_cp*a*b*dn+12.0_cp*a*b+4.0_cp*b**2*dn**2+6.0_cp &
      &             *b**2*dn+2.0_cp*b**2)/(4.0_cp*dn*(dn+1.0_cp))
      stencil(ku+3)=-a**2*(2.0_cp*a**2*dm**2*dn+a**2*dm**2+4.0_cp*a**2*dn**2 &
      &             -2.0_cp*a**2*dn-6.0_cp*a**2+4.0_cp*a*b*dm**2*dn+4.0_cp*a*b*dm**2-4.0_cp &
      &             *a*b*dn**3+12.0_cp*a*b*dn**2-16.0_cp*a*b+10.0_cp*b**2*dn**3-6.0_cp* &
      &             b**2*dn**2-12.0_cp*b**2*dn+4.0_cp*b**2)/(8.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp &
      &             ))
      stencil(ku+4)=a**3*b*(dm**2-2.0_cp*dn**2+5.0_cp*dn-3.0_cp)/(4.0_cp*dn &
      &             *(dn-1.0_cp))
      stencil(ku+5)=a**4*(dm**2-dn**2+3.0_cp*dn-2.0_cp)/(16.0_cp*dn*(dn-1.0_cp &
      &             ))
      stencil(ku+6:len_stencil) = 0.0_cp

      call mirror_stencil(n, stencil, len_stencil)

   end function intcheb2rmult2laplrot
!------------------------------------------------------------------------------
   function intcheb4rmult3(a, b, n, len_stencil) result(stencil)

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

      dn = real(n, cp)
      ku = (len_stencil-1)/2

      stencil(1:ku-7) = 0.0_cp
      stencil(ku-6)=a**7/(128.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-5)=3.0_cp*a**6*b/(64.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-4)=-a**5*(a**2*dn+11.0_cp*a**2-12.0_cp*b**2*dn+12.0_cp*b**2)/  &
      &             (128.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-3)=-a**4*b*(3.0_cp*a**2*dn+15.0_cp*a**2-2.0_cp*b**2*dn+2.0_cp &
      &             *b**2)/(32.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*    &
      &             (dn+3.0_cp))
      stencil(ku-2)=-3.0_cp*a**5*(a**2*dn-6.0_cp*a**2+12.0_cp*b**2*dn-24.0_cp &
      &             *b**2)/(128.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*  &
      &             (dn+2.0_cp))
      stencil(ku-1)=-a**4*b*(3.0_cp*a**2*dn-51.0_cp*a**2+16.0_cp*b**2*dn-32.0_cp &
      &             *b**2)/(64.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*      &
      &             (dn+3.0_cp))
      stencil(ku)  =3.0_cp*a**5*(a**2*dn**2+5.0_cp*a**2*dn-34.0_cp*a**2+8.0_cp &
      &             *b**2*dn**2+40.0_cp*b**2*dn-192.0_cp*b**2)/(128.0_cp*dn* &
      &             (dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+1)=3.0_cp*a**4*b*(a**2*dn**2-19.0_cp*a**2+2.0_cp*b**2*dn**2    &
      &             -18.0_cp*b**2)/(16.0_cp*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*&
      &             (dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+2)=3.0_cp*a**5*(a**2*dn**2-5.0_cp*a**2*dn-34.0_cp*a**2+8.0_cp &
      &             *b**2*dn**2-40.0_cp*b**2*dn-192.0_cp*b**2)/(128.0_cp*dn*   &
      &             (dn-3.0_cp)*(dn-2.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+3)=-a**4*b*(3.0_cp*a**2*dn+51.0_cp*a**2+16.0_cp*b**2*dn+32.0_cp &
      &             *b**2)/(64.0_cp*dn*(dn-3.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*      &
      &             (dn+2.0_cp))
      stencil(ku+4)=-3.0_cp*a**5*(a**2*dn+6.0_cp*a**2+12.0_cp*b**2*dn+24.0_cp &
      &             *b**2)/(128.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*  &
      &             (dn+2.0_cp))
      stencil(ku+5)=-a**4*b*(3.0_cp*a**2*dn-15.0_cp*a**2-2.0_cp*b**2*dn-2.0_cp &
      &             *b**2)/(32.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*    &
      &             (dn+1.0_cp))
      stencil(ku+6)=-a**5*(a**2*dn-11.0_cp*a**2-12.0_cp*b**2*dn-12.0_cp*b**2)/ &
      &             (128.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+7)=3.0_cp*a**6*b/(64.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp))
      stencil(ku+8)=a**7/(128.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp))
      stencil(ku+9:len_stencil) = 0.0_cp

      call mirror_stencil(n, stencil, len_stencil)

   end function intcheb4rmult3
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

      call mirror_stencil(n, stencil, len_stencil)

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

      call mirror_stencil(n, stencil, len_stencil)

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

      call mirror_stencil(n, stencil, len_stencil)

   end function intcheb4rmult4lapl2
!------------------------------------------------------------------------------
   function intcheb4rmult4laplrot(a, b, m, n, len_stencil) result(stencil)

      !-- Input variables
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
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

      stencil(1:ku-8) = 0.0_cp
      stencil(ku-7)=-a**8*(dm**2-dn**2-11.0_cp*dn-30.0_cp)/(256.0_cp*dn*(dn &
      &             +1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-6)=-a**7*b*(2.0_cp*dm**2-3.0_cp*dn**2-30.0_cp*dn-75.0_cp &
      &             )/(64.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-5)=a**6*(2.0_cp*a**2*dm**2*dn+4.0_cp*a**2*dm**2-9.0_cp*a**2 &
      &             *dn**2-57.0_cp*a**2*dn-54.0_cp*a**2+4.0_cp*a*b*dm**2*dn-4.0_cp*a*b*dm**2 &
      &             -4.0_cp*a*b*dn**3-44.0_cp*a*b*dn**2-96.0_cp*a*b*dn+144.0_cp*a*b-10.0_cp &
      &             *b**2*dm**2*dn+10.0_cp*b**2*dm**2+28.0_cp*b**2*dn**3+218.0_cp*b**2*dn**2 &
      &             +294.0_cp*b**2*dn-540.0_cp*b**2)/(128.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp &
      &             )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-4)=a**5*b*(6.0_cp*a**2*dm**2*dn+18.0_cp*a**2*dm**2+a**2*dn**3 &
      &             -29.0_cp*a**2*dn**2-221.0_cp*a**2*dn-327.0_cp*a**2+8.0_cp*a*b*dm**2 &
      &             *dn-8.0_cp*a*b*dm**2-16.0_cp*a*b*dn**3-148.0_cp*a*b*dn**2-256.0_cp* &
      &             a*b*dn+420.0_cp*a*b-4.0_cp*b**2*dm**2*dn+4.0_cp*b**2*dm**2+32.0_cp* &
      &             b**2*dn**3+206.0_cp*b**2*dn**2+200.0_cp*b**2*dn-438.0_cp*b**2)/(64.0_cp &
      &             *dn*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-3)=-a**4*(2.0_cp*a**4*dm**2*dn**2+12.0_cp*a**4*dm**2*dn-2.0_cp &
      &             *a**4*dm**2+2.0_cp*a**4*dn**4+11.0_cp*a**4*dn**3-65.0_cp*a**4*dn**2 &
      &             -236.0_cp*a**4*dn+252.0_cp*a**4+8.0_cp*a**3*b*dm**2*dn**2+24.0_cp*a**3 &
      &             *b*dm**2*dn-80.0_cp*a**3*b*dm**2+8.0_cp*a**3*b*dn**4+8.0_cp*a**3*b*dn**3 &
      &             -320.0_cp*a**3*b*dn**2-416.0_cp*a**3*b*dn+1920.0_cp*a**3*b-28.0_cp* &
      &             a**2*b**2*dm**2*dn**2-36.0_cp*a**2*b**2*dm**2*dn+184.0_cp*a**2*b**2 &
      &             *dm**2-8.0_cp*a**2*b**2*dn**4+208.0_cp*a**2*b**2*dn**3+1040.0_cp*a**2* &
      &             b**2*dn**2-592.0_cp*a**2*b**2*dn-4512.0_cp*a**2*b**2-16.0_cp*a*b**3 &
      &             *dm**2*dn**2+48.0_cp*a*b**3*dm**2*dn-32.0_cp*a*b**3*dm**2+96.0_cp*a &
      &             *b**3*dn**4+528.0_cp*a*b**3*dn**3-528.0_cp*a*b**3*dn**2-3552.0_cp*a &
      &             *b**3*dn+3456.0_cp*a*b**3-72.0_cp*b**4*dn**4-216.0_cp*b**4*dn**3+528.0_cp &
      &             *b**4*dn**2+1008.0_cp*b**4*dn-1248.0_cp*b**4)/(128.0_cp*dn*(dn-2.0_cp &
      &             )*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-2)=-a**3*b*(6.0_cp*a**4*dm**2*dn+12.0_cp*a**4*dm**2+5.0_cp &
      &             *a**4*dn**3+6.0_cp*a**4*dn**2-173.0_cp*a**4*dn-42.0_cp*a**4+24.0_cp &
      &             *a**3*b*dm**2*dn-48.0_cp*a**3*b*dm**2+16.0_cp*a**3*b*dn**3-132.0_cp &
      &             *a**3*b*dn**2-532.0_cp*a**3*b*dn+1464.0_cp*a**3*b-12.0_cp*a**2*b**2 &
      &             *dm**2*dn+24.0_cp*a**2*b**2*dm**2+198.0_cp*a**2*b**2*dn**2+102.0_cp &
      &             *a**2*b**2*dn-996.0_cp*a**2*b**2+64.0_cp*a*b**3*dn**3+48.0_cp*a*b** &
      &             3*dn**2-592.0_cp*a*b**3*dn+480.0_cp*a*b**3-16.0_cp*b**4*dn**3+24.0_cp &
      &             *b**4*dn**2+40.0_cp*b**4*dn-48.0_cp*b**4)/(64.0_cp*dn*(dn-2.0_cp)*(dn &
      &             -1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku-1)=-a**3*(2.0_cp*a**5*dm**2*dn**2-22.0_cp*a**5*dm**2*dn+36.0_cp &
      &             *a**5*dm**2-27.0_cp*a**5*dn**3+42.0_cp*a**5*dn**2+495.0_cp*a**5*dn-1062.0_cp &
      &             *a**5+4.0_cp*a**4*b*dm**2*dn**2-80.0_cp*a**4*b*dm**2*dn+204.0_cp*a**4* &
      &             b*dm**2-4.0_cp*a**4*b*dn**4-96.0_cp*a**4*b*dn**3+292.0_cp*a**4*b*dn**2 &
      &             +1584.0_cp*a**4*b*dn-4464.0_cp*a**4*b+22.0_cp*a**3*b**2*dm**2*dn**2 &
      &             +40.0_cp*a**3*b**2*dm**2*dn-318.0_cp*a**3*b**2*dm**2+28.0_cp*a**3*b**2 &
      &             *dn**4-30.0_cp*a**3*b**2*dn**3-856.0_cp*a**3*b**2*dn**2+630.0_cp*a**3* &
      &             b**2*dn+4356.0_cp*a**3*b**2+64.0_cp*a**2*b**3*dm**2*dn**2-320.0_cp* &
      &             a**2*b**3*dm**2*dn+384.0_cp*a**2*b**3*dm**2-864.0_cp*a**2*b**3*dn** &
      &             3+1152.0_cp*a**2*b**3*dn**2+10656.0_cp*a**2*b**3*dn-19008.0_cp*a**2 &
      &             *b**3+32.0_cp*a*b**4*dn**4+192.0_cp*a*b**4*dn**3-896.0_cp*a*b**4*dn**2 &
      &             -1248.0_cp*a*b**4*dn+4032.0_cp*a*b**4+64.0_cp*b**5*dn**4-192.0_cp*b**5 &
      &             *dn**3-448.0_cp*b**5*dn**2+1728.0_cp*b**5*dn-1152.0_cp*b**5)/(128.0_cp &
      &             *dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+3.0_cp))
      stencil(ku)=a**3*b*(2.0_cp*a**4*dm**2*dn**2+10.0_cp*a**4*dm**2*dn+12.0_cp &
      &             *a**4*dm**2+a**4*dn**4+39.0_cp*a**4*dn**3-115.0_cp*a**4*dn**2-531.0_cp &
      &             *a**4*dn+1134.0_cp*a**4+16.0_cp*a**3*b*dm**2*dn**2+80.0_cp*a**3*b*dm**2 &
      &             *dn-384.0_cp*a**3*b*dm**2+32.0_cp*a**3*b*dn**4+168.0_cp*a**3*b*dn** &
      &             3-1232.0_cp*a**3*b*dn**2-2232.0_cp*a**3*b*dn+10656.0_cp*a**3*b-8.0_cp &
      &             *a**2*b**2*dm**2*dn**2-40.0_cp*a**2*b**2*dm**2*dn+192.0_cp*a**2*b** &
      &             2*dm**2-32.0_cp*a**2*b**2*dn**4+12.0_cp*a**2*b**2*dn**3+728.0_cp*a**2* &
      &             b**2*dn**2+12.0_cp*a**2*b**2*dn-4320.0_cp*a**2*b**2+64.0_cp*a*b**3*dn**4 &
      &             +336.0_cp*a*b**3*dn**3-1504.0_cp*a*b**3*dn**2-3024.0_cp*a*b**3*dn+8352.0_cp &
      &             *a*b**3-16.0_cp*b**4*dn**4+24.0_cp*b**4*dn**3+160.0_cp*b**4*dn**2-216.0_cp &
      &             *b**4*dn-144.0_cp*b**4)/(64.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &             )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+1)=a**3*(5.0_cp*a**5*dm**2*dn**2-65.0_cp*a**5*dm**2+3.0_cp &
      &             *a**5*dn**4-171.0_cp*a**5*dn**2+1656.0_cp*a**5+16.0_cp*a**4*b*dm**2 &
      &             *dn**2-304.0_cp*a**4*b*dm**2+16.0_cp*a**4*b*dn**4-688.0_cp*a**4*b*dn**2 &
      &             +6336.0_cp*a**4*b+8.0_cp*a**3*b**2*dm**2*dn**2+328.0_cp*a**3*b**2*dm**2 &
      &             -16.0_cp*a**3*b**2*dn**4+256.0_cp*a**3*b**2*dn**2-3168.0_cp*a**3*b**2+ &
      &             96.0_cp*a**2*b**3*dm**2*dn**2-864.0_cp*a**2*b**3*dm**2+192.0_cp*a** &
      &             2*b**3*dn**4-5376.0_cp*a**2*b**3*dn**2+32832.0_cp*a**2*b**3-80.0_cp &
      &             *a*b**4*dn**4+1328.0_cp*a*b**4*dn**2-5472.0_cp*a*b**4+128.0_cp*b**5 &
      &             *dn**4-1664.0_cp*b**5*dn**2+4608.0_cp*b**5)/(128.0_cp*(dn-3.0_cp)*(dn &
      &             -2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+2)=a**3*b*(2.0_cp*a**4*dm**2*dn**2-10.0_cp*a**4*dm**2*dn &
      &             +12.0_cp*a**4*dm**2+a**4*dn**4-39.0_cp*a**4*dn**3-115.0_cp*a**4*dn**2+ &
      &             531.0_cp*a**4*dn+1134.0_cp*a**4+16.0_cp*a**3*b*dm**2*dn**2-80.0_cp* &
      &             a**3*b*dm**2*dn-384.0_cp*a**3*b*dm**2+32.0_cp*a**3*b*dn**4-168.0_cp &
      &             *a**3*b*dn**3-1232.0_cp*a**3*b*dn**2+2232.0_cp*a**3*b*dn+10656.0_cp &
      &             *a**3*b-8.0_cp*a**2*b**2*dm**2*dn**2+40.0_cp*a**2*b**2*dm**2*dn+192.0_cp &
      &             *a**2*b**2*dm**2-32.0_cp*a**2*b**2*dn**4-12.0_cp*a**2*b**2*dn**3+728.0_cp &
      &             *a**2*b**2*dn**2-12.0_cp*a**2*b**2*dn-4320.0_cp*a**2*b**2+64.0_cp*a &
      &             *b**3*dn**4-336.0_cp*a*b**3*dn**3-1504.0_cp*a*b**3*dn**2+3024.0_cp* &
      &             a*b**3*dn+8352.0_cp*a*b**3-16.0_cp*b**4*dn**4-24.0_cp*b**4*dn**3+160.0_cp &
      &             *b**4*dn**2+216.0_cp*b**4*dn-144.0_cp*b**4)/(64.0_cp*dn*(dn-3.0_cp) &
      &             *(dn-2.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+3)=-a**3*(2.0_cp*a**5*dm**2*dn**2+22.0_cp*a**5*dm**2*dn+36.0_cp &
      &             *a**5*dm**2+27.0_cp*a**5*dn**3+42.0_cp*a**5*dn**2-495.0_cp*a**5*dn-1062.0_cp &
      &             *a**5+4.0_cp*a**4*b*dm**2*dn**2+80.0_cp*a**4*b*dm**2*dn+204.0_cp*a**4* &
      &             b*dm**2-4.0_cp*a**4*b*dn**4+96.0_cp*a**4*b*dn**3+292.0_cp*a**4*b*dn**2 &
      &             -1584.0_cp*a**4*b*dn-4464.0_cp*a**4*b+22.0_cp*a**3*b**2*dm**2*dn**2 &
      &             -40.0_cp*a**3*b**2*dm**2*dn-318.0_cp*a**3*b**2*dm**2+28.0_cp*a**3*b**2 &
      &             *dn**4+30.0_cp*a**3*b**2*dn**3-856.0_cp*a**3*b**2*dn**2-630.0_cp*a**3* &
      &             b**2*dn+4356.0_cp*a**3*b**2+64.0_cp*a**2*b**3*dm**2*dn**2+320.0_cp* &
      &             a**2*b**3*dm**2*dn+384.0_cp*a**2*b**3*dm**2+864.0_cp*a**2*b**3*dn** &
      &             3+1152.0_cp*a**2*b**3*dn**2-10656.0_cp*a**2*b**3*dn-19008.0_cp*a**2 &
      &             *b**3+32.0_cp*a*b**4*dn**4-192.0_cp*a*b**4*dn**3-896.0_cp*a*b**4*dn**2 &
      &             +1248.0_cp*a*b**4*dn+4032.0_cp*a*b**4+64.0_cp*b**5*dn**4+192.0_cp*b**5 &
      &             *dn**3-448.0_cp*b**5*dn**2-1728.0_cp*b**5*dn-1152.0_cp*b**5)/(128.0_cp &
      &             *dn*(dn-3.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+4)=-a**3*b*(6.0_cp*a**4*dm**2*dn-12.0_cp*a**4*dm**2+5.0_cp &
      &             *a**4*dn**3-6.0_cp*a**4*dn**2-173.0_cp*a**4*dn+42.0_cp*a**4+24.0_cp &
      &             *a**3*b*dm**2*dn+48.0_cp*a**3*b*dm**2+16.0_cp*a**3*b*dn**3+132.0_cp &
      &             *a**3*b*dn**2-532.0_cp*a**3*b*dn-1464.0_cp*a**3*b-12.0_cp*a**2*b**2 &
      &             *dm**2*dn-24.0_cp*a**2*b**2*dm**2-198.0_cp*a**2*b**2*dn**2+102.0_cp &
      &             *a**2*b**2*dn+996.0_cp*a**2*b**2+64.0_cp*a*b**3*dn**3-48.0_cp*a*b** &
      &             3*dn**2-592.0_cp*a*b**3*dn-480.0_cp*a*b**3-16.0_cp*b**4*dn**3-24.0_cp &
      &             *b**4*dn**2+40.0_cp*b**4*dn+48.0_cp*b**4)/(64.0_cp*dn*(dn-2.0_cp)*(dn &
      &             -1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+5)=-a**4*(2.0_cp*a**4*dm**2*dn**2-12.0_cp*a**4*dm**2*dn-2.0_cp &
      &             *a**4*dm**2+2.0_cp*a**4*dn**4-11.0_cp*a**4*dn**3-65.0_cp*a**4*dn**2 &
      &             +236.0_cp*a**4*dn+252.0_cp*a**4+8.0_cp*a**3*b*dm**2*dn**2-24.0_cp*a**3 &
      &             *b*dm**2*dn-80.0_cp*a**3*b*dm**2+8.0_cp*a**3*b*dn**4-8.0_cp*a**3*b*dn**3 &
      &             -320.0_cp*a**3*b*dn**2+416.0_cp*a**3*b*dn+1920.0_cp*a**3*b-28.0_cp* &
      &             a**2*b**2*dm**2*dn**2+36.0_cp*a**2*b**2*dm**2*dn+184.0_cp*a**2*b**2 &
      &             *dm**2-8.0_cp*a**2*b**2*dn**4-208.0_cp*a**2*b**2*dn**3+1040.0_cp*a**2* &
      &             b**2*dn**2+592.0_cp*a**2*b**2*dn-4512.0_cp*a**2*b**2-16.0_cp*a*b**3 &
      &             *dm**2*dn**2-48.0_cp*a*b**3*dm**2*dn-32.0_cp*a*b**3*dm**2+96.0_cp*a &
      &             *b**3*dn**4-528.0_cp*a*b**3*dn**3-528.0_cp*a*b**3*dn**2+3552.0_cp*a &
      &             *b**3*dn+3456.0_cp*a*b**3-72.0_cp*b**4*dn**4+216.0_cp*b**4*dn**3+528.0_cp &
      &             *b**4*dn**2-1008.0_cp*b**4*dn-1248.0_cp*b**4)/(128.0_cp*dn*(dn-3.0_cp &
      &             )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+6)=a**5*b*(6.0_cp*a**2*dm**2*dn-18.0_cp*a**2*dm**2+a**2*dn**3 &
      &             +29.0_cp*a**2*dn**2-221.0_cp*a**2*dn+327.0_cp*a**2+8.0_cp*a*b*dm**2 &
      &             *dn+8.0_cp*a*b*dm**2-16.0_cp*a*b*dn**3+148.0_cp*a*b*dn**2-256.0_cp* &
      &             a*b*dn-420.0_cp*a*b-4.0_cp*b**2*dm**2*dn-4.0_cp*b**2*dm**2+32.0_cp* &
      &             b**2*dn**3-206.0_cp*b**2*dn**2+200.0_cp*b**2*dn+438.0_cp*b**2)/(64.0_cp &
      &             *dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+7)=a**6*(2.0_cp*a**2*dm**2*dn-4.0_cp*a**2*dm**2+9.0_cp*a**2 &
      &             *dn**2-57.0_cp*a**2*dn+54.0_cp*a**2+4.0_cp*a*b*dm**2*dn+4.0_cp*a*b*dm**2 &
      &             -4.0_cp*a*b*dn**3+44.0_cp*a*b*dn**2-96.0_cp*a*b*dn-144.0_cp*a*b-10.0_cp &
      &             *b**2*dm**2*dn-10.0_cp*b**2*dm**2+28.0_cp*b**2*dn**3-218.0_cp*b**2*dn**2 &
      &             +294.0_cp*b**2*dn+540.0_cp*b**2)/(128.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp &
      &             )*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+8)=-a**7*b*(2.0_cp*dm**2-3.0_cp*dn**2+30.0_cp*dn-75.0_cp &
      &             )/(64.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp))
      stencil(ku+9)=-a**8*(dm**2-dn**2+11.0_cp*dn-30.0_cp)/(256.0_cp*dn*(dn &
      &             -3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp))
      stencil(ku+10:len_stencil) = 0.0_cp

      call mirror_stencil(n, stencil, len_stencil)

   end function intcheb4rmult4laplrot
!------------------------------------------------------------------------------
   function intcheb4rmult4laplrot2(a, b, m, n, len_stencil) result(stencil)

      !-- Input variables
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
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

      stencil(1:ku-6) = 0.0_cp
      stencil(ku-5)=a**6*(dm-dn-6.0_cp)*(dm+dn+6.0_cp)*(dm**2-dn**2-7.0_cp &
      &             *dn-12.0_cp)/(64.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-4)=a**5*b*(2.0_cp*dm**4-8.0_cp*dm**2*dn**2-67.0_cp*dm**2 &
      &             *dn-147.0_cp*dm**2+6.0_cp*dn**4+101.0_cp*dn**3+628.0_cp*dn**2+1707.0_cp &
      &             *dn+1710.0_cp)/(32.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-3)=-a**4*(3.0_cp*a**2*dm**4*dn+3.0_cp*a**2*dm**4-2.0_cp* &
      &             a**2*dm**2*dn**3-30.0_cp*a**2*dm**2*dn**2-108.0_cp*a**2*dm**2*dn-46.0_cp &
      &             *a**2*dm**2-a**2*dn**5-3.0_cp*a**2*dn**4+68.0_cp*a**2*dn**3+420.0_cp &
      &             *a**2*dn**2+608.0_cp*a**2*dn-192.0_cp*a**2+4.0_cp*a*b*dm**4*dn-4.0_cp &
      &             *a*b*dm**4-8.0_cp*a*b*dm**2*dn**3-72.0_cp*a*b*dm**2*dn**2-128.0_cp* &
      &             a*b*dm**2*dn+208.0_cp*a*b*dm**2+4.0_cp*a*b*dn**5+76.0_cp*a*b*dn**4+512.0_cp &
      &             *a*b*dn**3+1328.0_cp*a*b*dn**2+384.0_cp*a*b*dn-2304.0_cp*a*b+20.0_cp &
      &             *b**2*dm**2*dn**3+114.0_cp*b**2*dm**2*dn**2+86.0_cp*b**2*dm**2*dn-220.0_cp &
      &             *b**2*dm**2-28.0_cp*b**2*dn**5-372.0_cp*b**2*dn**4-1710.0_cp*b**2*dn**3 &
      &             -2742.0_cp*b**2*dn**2+772.0_cp*b**2*dn+4080.0_cp*b**2)/(32.0_cp*dn* &
      &             (dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-2)=-a**3*b*(6.0_cp*a**2*dm**4-8.0_cp*a**2*dm**2*dn**2-91.0_cp &
      &             *a**2*dm**2*dn-165.0_cp*a**2*dm**2-14.0_cp*a**2*dn**4-55.0_cp*a**2*dn**3 &
      &             +212.0_cp*a**2*dn**2+967.0_cp*a**2*dn+498.0_cp*a**2-32.0_cp*a*b*dm**2* &
      &             dn**2-112.0_cp*a*b*dm**2*dn+144.0_cp*a*b*dm**2+32.0_cp*a*b*dn**4+400.0_cp &
      &             *a*b*dn**3+1488.0_cp*a*b*dn**2+888.0_cp*a*b*dn-2808.0_cp*a*b+16.0_cp &
      &             *b**2*dm**2*dn**2+12.0_cp*b**2*dm**2*dn-28.0_cp*b**2*dm**2-64.0_cp* &
      &             b**2*dn**4-480.0_cp*b**2*dn**3-904.0_cp*b**2*dn**2+308.0_cp*b**2*dn &
      &             +1140.0_cp*b**2)/(32.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku-1)=a**2*(15.0_cp*a**4*dm**4*dn-15.0_cp*a**4*dm**4+2.0_cp &
      &             *a**4*dm**2*dn**3-57.0_cp*a**4*dm**2*dn**2-249.0_cp*a**4*dm**2*dn+440.0_cp &
      &             *a**4*dm**2-a**4*dn**5-48.0_cp*a**4*dn**4-83.0_cp*a**4*dn**3+840.0_cp &
      &             *a**4*dn**2+1092.0_cp*a**4*dn-3600.0_cp*a**4+32.0_cp*a**3*b*dm**4*dn &
      &             -64.0_cp*a**3*b*dm**4-192.0_cp*a**3*b*dm**2*dn**2-512.0_cp*a**3*b*dm**2 &
      &             *dn+1792.0_cp*a**3*b*dm**2-32.0_cp*a**3*b*dn**5-256.0_cp*a**3*b*dn**4+ &
      &             3328.0_cp*a**3*b*dn**2+2048.0_cp*a**3*b*dn-12288.0_cp*a**3*b+32.0_cp &
      &             *a**2*b**2*dm**2*dn**3+216.0_cp*a**2*b**2*dm**2*dn**2-24.0_cp*a**2* &
      &             b**2*dm**2*dn-1072.0_cp*a**2*b**2*dm**2+128.0_cp*a**2*b**2*dn**5+384.0_cp &
      &             *a**2*b**2*dn**4-1488.0_cp*a**2*b**2*dn**3-3888.0_cp*a**2*b**2*dn** &
      &             2+4624.0_cp*a**2*b**2*dn+7968.0_cp*a**2*b**2+64.0_cp*a*b**3*dm**2*dn**3 &
      &             -448.0_cp*a*b**3*dm**2*dn+384.0_cp*a*b**3*dm**2-192.0_cp*a*b**3*dn**5- &
      &             1536.0_cp*a*b**3*dn**4-1696.0_cp*a*b**3*dn**3+9600.0_cp*a*b**3*dn** &
      &             2+12064.0_cp*a*b**3*dn-18240.0_cp*a*b**3+144.0_cp*b**4*dn**5+512.0_cp &
      &             *b**4*dn**4-800.0_cp*b**4*dn**3-2720.0_cp*b**4*dn**2+1616.0_cp*b**4 &
      &             *dn+1248.0_cp*b**4)/(64.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp) &
      &             *(dn+3.0_cp))
      stencil(ku)=a*b*(2.0_cp*a**4*dm**4*dn+16.0_cp*a**4*dm**4-33.0_cp*a**4* &
      &             dm**2*dn**2-105.0_cp*a**4*dm**2*dn-108.0_cp*a**4*dm**2+6.0_cp*a**4*dn**5 &
      &             -7.0_cp*a**4*dn**4-40.0_cp*a**4*dn**3+439.0_cp*a**4*dn**2+490.0_cp* &
      &             a**4*dn-1536.0_cp*a**4-16.0_cp*a**3*b*dm**2*dn**3-168.0_cp*a**3*b*dm**2 &
      &             *dn**2-56.0_cp*a**3*b*dm**2*dn+912.0_cp*a**3*b*dm**2-48.0_cp*a**3*b &
      &             *dn**5-264.0_cp*a**3*b*dn**4+424.0_cp*a**3*b*dn**3+3492.0_cp*a**3*b &
      &             *dn**2-268.0_cp*a**3*b*dn-11064.0_cp*a**3*b+8.0_cp*a**2*b**2*dm**2*dn**3 &
      &             +18.0_cp*a**2*b**2*dm**2*dn**2-38.0_cp*a**2*b**2*dm**2*dn-60.0_cp*a**2 &
      &             *b**2*dm**2+64.0_cp*a**2*b**2*dn**5+192.0_cp*a**2*b**2*dn**4-452.0_cp &
      &             *a**2*b**2*dn**3-1218.0_cp*a**2*b**2*dn**2+922.0_cp*a**2*b**2*dn+1524.0_cp &
      &             *a**2*b**2-64.0_cp*a*b**3*dn**5-352.0_cp*a*b**3*dn**4+2080.0_cp*a*b**3 &
      &             *dn**2+1024.0_cp*a*b**3*dn-2688.0_cp*a*b**3+16.0_cp*b**4*dn**5+24.0_cp &
      &             *b**4*dn**4-128.0_cp*b**4*dn**3-72.0_cp*b**4*dn**2+256.0_cp*b**4*dn &
      &             -96.0_cp*b**4)/(16.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &             ))
      stencil(ku+1)=-a*(5.0_cp*a**5*dm**4*dn**2-35.0_cp*a**5*dm**4+2.0_cp &
      &             *a**5*dm**2*dn**4-110.0_cp*a**5*dm**2*dn**2+738.0_cp*a**5*dm**2+a** &
      &             5*dn**6-77.0_cp*a**5*dn**4+1156.0_cp*a**5*dn**2-4896.0_cp*a**5+12.0_cp &
      &             *a**4*b*dm**4*dn**2-108.0_cp*a**4*b*dm**4+8.0_cp*a**4*b*dm**2*dn**4 &
      &             -344.0_cp*a**4*b*dm**2*dn**2+2448.0_cp*a**4*b*dm**2+12.0_cp*a**4*b*dn**6 &
      &             -380.0_cp*a**4*b*dn**4+4112.0_cp*a**4*b*dn**2-14976.0_cp*a**4*b-4.0_cp &
      &             *a**3*b**2*dm**2*dn**4+148.0_cp*a**3*b**2*dm**2*dn**2-1008.0_cp*a** &
      &             3*b**2*dm**2-36.0_cp*a**3*b**2*dn**6+490.0_cp*a**3*b**2*dn**4-1774.0_cp &
      &             *a**3*b**2*dn**2+2520.0_cp*a**3*b**2+32.0_cp*a**2*b**3*dm**2*dn**4-416.0_cp &
      &             *a**2*b**3*dm**2*dn**2+1152.0_cp*a**2*b**3*dm**2+96.0_cp*a**2*b**3*dn**6 &
      &             -2192.0_cp*a**2*b**3*dn**4+15728.0_cp*a**2*b**3*dn**2-33984.0_cp*a**2* &
      &             b**3-56.0_cp*a*b**4*dn**6+704.0_cp*a*b**4*dn**4-1704.0_cp*a*b**4*dn**2 &
      &             -864.0_cp*a*b**4+32.0_cp*b**5*dn**6-448.0_cp*b**5*dn**4+1568.0_cp*b**5 &
      &             *dn**2-1152.0_cp*b**5)/(16.0_cp*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp) &
      &             *(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+2)=a*b*(2.0_cp*a**4*dm**4*dn-16.0_cp*a**4*dm**4+33.0_cp* &
      &             a**4*dm**2*dn**2-105.0_cp*a**4*dm**2*dn+108.0_cp*a**4*dm**2+6.0_cp* &
      &             a**4*dn**5+7.0_cp*a**4*dn**4-40.0_cp*a**4*dn**3-439.0_cp*a**4*dn**2 &
      &             +490.0_cp*a**4*dn+1536.0_cp*a**4-16.0_cp*a**3*b*dm**2*dn**3+168.0_cp &
      &             *a**3*b*dm**2*dn**2-56.0_cp*a**3*b*dm**2*dn-912.0_cp*a**3*b*dm**2-48.0_cp &
      &             *a**3*b*dn**5+264.0_cp*a**3*b*dn**4+424.0_cp*a**3*b*dn**3-3492.0_cp &
      &             *a**3*b*dn**2-268.0_cp*a**3*b*dn+11064.0_cp*a**3*b+8.0_cp*a**2*b**2 &
      &             *dm**2*dn**3-18.0_cp*a**2*b**2*dm**2*dn**2-38.0_cp*a**2*b**2*dm**2*dn &
      &             +60.0_cp*a**2*b**2*dm**2+64.0_cp*a**2*b**2*dn**5-192.0_cp*a**2*b**2 &
      &             *dn**4-452.0_cp*a**2*b**2*dn**3+1218.0_cp*a**2*b**2*dn**2+922.0_cp* &
      &             a**2*b**2*dn-1524.0_cp*a**2*b**2-64.0_cp*a*b**3*dn**5+352.0_cp*a*b**3* &
      &             dn**4-2080.0_cp*a*b**3*dn**2+1024.0_cp*a*b**3*dn+2688.0_cp*a*b**3+16.0_cp &
      &             *b**4*dn**5-24.0_cp*b**4*dn**4-128.0_cp*b**4*dn**3+72.0_cp*b**4*dn**2+ &
      &             256.0_cp*b**4*dn+96.0_cp*b**4)/(16.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)* &
      &             (dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+3)=a**2*(15.0_cp*a**4*dm**4*dn+15.0_cp*a**4*dm**4+2.0_cp &
      &             *a**4*dm**2*dn**3+57.0_cp*a**4*dm**2*dn**2-249.0_cp*a**4*dm**2*dn-440.0_cp &
      &             *a**4*dm**2-a**4*dn**5+48.0_cp*a**4*dn**4-83.0_cp*a**4*dn**3-840.0_cp &
      &             *a**4*dn**2+1092.0_cp*a**4*dn+3600.0_cp*a**4+32.0_cp*a**3*b*dm**4*dn &
      &             +64.0_cp*a**3*b*dm**4+192.0_cp*a**3*b*dm**2*dn**2-512.0_cp*a**3*b*dm**2 &
      &             *dn-1792.0_cp*a**3*b*dm**2-32.0_cp*a**3*b*dn**5+256.0_cp*a**3*b*dn**4- &
      &             3328.0_cp*a**3*b*dn**2+2048.0_cp*a**3*b*dn+12288.0_cp*a**3*b+32.0_cp &
      &             *a**2*b**2*dm**2*dn**3-216.0_cp*a**2*b**2*dm**2*dn**2-24.0_cp*a**2* &
      &             b**2*dm**2*dn+1072.0_cp*a**2*b**2*dm**2+128.0_cp*a**2*b**2*dn**5-384.0_cp &
      &             *a**2*b**2*dn**4-1488.0_cp*a**2*b**2*dn**3+3888.0_cp*a**2*b**2*dn** &
      &             2+4624.0_cp*a**2*b**2*dn-7968.0_cp*a**2*b**2+64.0_cp*a*b**3*dm**2*dn**3 &
      &             -448.0_cp*a*b**3*dm**2*dn-384.0_cp*a*b**3*dm**2-192.0_cp*a*b**3*dn**5+ &
      &             1536.0_cp*a*b**3*dn**4-1696.0_cp*a*b**3*dn**3-9600.0_cp*a*b**3*dn** &
      &             2+12064.0_cp*a*b**3*dn+18240.0_cp*a*b**3+144.0_cp*b**4*dn**5-512.0_cp &
      &             *b**4*dn**4-800.0_cp*b**4*dn**3+2720.0_cp*b**4*dn**2+1616.0_cp*b**4 &
      &             *dn-1248.0_cp*b**4)/(64.0_cp*dn*(dn-3.0_cp)*(dn-1.0_cp)*(dn+1.0_cp) &
      &             *(dn+2.0_cp))
      stencil(ku+4)=-a**3*b*(6.0_cp*a**2*dm**4-8.0_cp*a**2*dm**2*dn**2+91.0_cp &
      &             *a**2*dm**2*dn-165.0_cp*a**2*dm**2-14.0_cp*a**2*dn**4+55.0_cp*a**2*dn**3 &
      &             +212.0_cp*a**2*dn**2-967.0_cp*a**2*dn+498.0_cp*a**2-32.0_cp*a*b*dm**2* &
      &             dn**2+112.0_cp*a*b*dm**2*dn+144.0_cp*a*b*dm**2+32.0_cp*a*b*dn**4-400.0_cp &
      &             *a*b*dn**3+1488.0_cp*a*b*dn**2-888.0_cp*a*b*dn-2808.0_cp*a*b+16.0_cp &
      &             *b**2*dm**2*dn**2-12.0_cp*b**2*dm**2*dn-28.0_cp*b**2*dm**2-64.0_cp* &
      &             b**2*dn**4+480.0_cp*b**2*dn**3-904.0_cp*b**2*dn**2-308.0_cp*b**2*dn &
      &             +1140.0_cp*b**2)/(32.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+5)=-a**4*(3.0_cp*a**2*dm**4*dn-3.0_cp*a**2*dm**4-2.0_cp* &
      &             a**2*dm**2*dn**3+30.0_cp*a**2*dm**2*dn**2-108.0_cp*a**2*dm**2*dn+46.0_cp &
      &             *a**2*dm**2-a**2*dn**5+3.0_cp*a**2*dn**4+68.0_cp*a**2*dn**3-420.0_cp &
      &             *a**2*dn**2+608.0_cp*a**2*dn+192.0_cp*a**2+4.0_cp*a*b*dm**4*dn+4.0_cp &
      &             *a*b*dm**4-8.0_cp*a*b*dm**2*dn**3+72.0_cp*a*b*dm**2*dn**2-128.0_cp* &
      &             a*b*dm**2*dn-208.0_cp*a*b*dm**2+4.0_cp*a*b*dn**5-76.0_cp*a*b*dn**4+512.0_cp &
      &             *a*b*dn**3-1328.0_cp*a*b*dn**2+384.0_cp*a*b*dn+2304.0_cp*a*b+20.0_cp &
      &             *b**2*dm**2*dn**3-114.0_cp*b**2*dm**2*dn**2+86.0_cp*b**2*dm**2*dn+220.0_cp &
      &             *b**2*dm**2-28.0_cp*b**2*dn**5+372.0_cp*b**2*dn**4-1710.0_cp*b**2*dn**3 &
      &             +2742.0_cp*b**2*dn**2+772.0_cp*b**2*dn-4080.0_cp*b**2)/(32.0_cp*dn* &
      &             (dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+6)=a**5*b*(2.0_cp*dm**4-8.0_cp*dm**2*dn**2+67.0_cp*dm**2 &
      &             *dn-147.0_cp*dm**2+6.0_cp*dn**4-101.0_cp*dn**3+628.0_cp*dn**2-1707.0_cp &
      &             *dn+1710.0_cp)/(32.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp))
      stencil(ku+7)=a**6*(dm-dn+6.0_cp)*(dm+dn-6.0_cp)*(dm**2-dn**2+7.0_cp &
      &             *dn-12.0_cp)/(64.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp))
      stencil(ku+8:len_stencil) = 0.0_cp

      call mirror_stencil(n, stencil, len_stencil)

   end function intcheb4rmult4laplrot2
!------------------------------------------------------------------------------
   subroutine mirror_stencil(idx, stencil, len_stencil)

      !-- Input variables
      integer,  intent(in) :: idx
      integer,  intent(in) :: len_stencil

      !-- Output variable
      real(cp), intent(inout) :: stencil(len_stencil)

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

   end subroutine mirror_stencil
!------------------------------------------------------------------------------
end module chebsparselib
