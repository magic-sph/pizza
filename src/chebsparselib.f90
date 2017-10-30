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
   &         intcheb4rmult4hmult6, intcheb4rmult4hmult8laplrot,          &
   &         intcheb4rmult4hmult8laplrot2, intcheb4hmult2,               &
   &         intcheb2rmult2hmult2, intcheb2rmult2hmult2lapl,             &
   &         intcheb2rmult2hmult4, intcheb2rmult2hmult4laplaxi,          &
   &         intcheb4hmult2der1

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
   function intcheb2rmult2hmult2(a, b, n, len_stencil) result(stencil)

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

      stencil(1:ku-6) = 0.0_cp
      stencil(ku-5)=-a**6/(64.0_cp*dn*(dn+1.0_cp))
      stencil(ku-4)=-a**5*b/(8.0_cp*dn*(dn+1.0_cp))
      stencil(ku-3)=a**4*(a**2*dn+4.0_cp*a*b*dn-4.0_cp*a*b-10.0_cp*b**2*dn &
      &             +10.0_cp*b**2)/(32.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku-2)=a**3*b*(a**2*dn+a**2+4.0_cp*a*b*dn-4.0_cp*a*b-2.0_cp* &
      &             b**2*dn+2.0_cp*b**2)/(8.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku-1)=a**3*(a**3*dn-3.0_cp*a**3-16.0_cp*a**2*b+16.0_cp*a*b**2*  &
      &             dn+24.0_cp*a*b**2+32.0_cp*b**3*dn-32.0_cp*b**3)/(64.0_cp* &
      &             dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku)  =-a**3*b**2*(2.0_cp*a-b)/(4.0_cp*dn*(dn-1.0_cp))
      stencil(ku+1)=-a**3*(a**3+4.0_cp*a**2*b-2.0_cp*a*b**2+16.0_cp*b**3) &
      &             /(16.0_cp*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+2)=-a**3*b**2*(2.0_cp*a-b)/(4.0_cp*dn*(dn+1.0_cp))
      stencil(ku+3)=a**3*(a**3*dn+3.0_cp*a**3+16.0_cp*a**2*b+16.0_cp*a*b**2* &
      &             dn-24.0_cp*a*b**2+32.0_cp*b**3*dn+32.0_cp*b**3)/(64.0_cp*&
      &             dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+4)=a**3*b*(a**2*dn-a**2+4.0_cp*a*b*dn+4.0_cp*a*b-2.0_cp* &
      &             b**2*dn-2.0_cp*b**2)/(8.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+5)=a**4*(a**2*dn+4.0_cp*a*b*dn+4.0_cp*a*b-10.0_cp*b**2*dn &
      &             -10.0_cp*b**2)/(32.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+6)=-a**5*b/(8.0_cp*dn*(dn-1.0_cp))
      stencil(ku+7)=-a**6/(64.0_cp*dn*(dn-1.0_cp))
      stencil(ku+8:len_stencil) = 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb2rmult2hmult2
!------------------------------------------------------------------------------
   function intcheb2rmult2hmult4(a, b, n, len_stencil) result(stencil)

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

      stencil(1:ku-8) = 0.0_cp
      stencil(ku-7)=a**8/(256.0_cp*dn*(dn+1.0_cp))
      stencil(ku-6)=3.0_cp*a**7*b/(64.0_cp*dn*(dn+1.0_cp))
      stencil(ku-5)=-a**6*(2.0_cp*a**2*dn-a**2+8.0_cp*a*b*dn-8.0_cp*a*b-26.0_cp &
      &             *b**2*dn+26.0_cp*b**2)/(128.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku-4)=-a**5*b*(7.0_cp*a**2*dn-a**2+32.0_cp*a*b*dn-32.0_cp*a &
      &             *b-24.0_cp*b**2*dn+24.0_cp*b**2)/(64.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp &
      &             ))
      stencil(ku-3)=a**4*(2.0_cp*a**4*dn+a**4+16.0_cp*a**3*b*dn+4.0_cp*a**2* &
      &             b**2*dn-56.0_cp*a**2*b**2-160.0_cp*a*b**3*dn+160.0_cp*a*b**3+32.0_cp &
      &             *b**4*dn-32.0_cp*b**4)/(128.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku-2)=a**4*b*(3.0_cp*a**3*dn+5.0_cp*a**3+32.0_cp*a**2*b*dn+32.0_cp &
      &             *a**2*b+56.0_cp*a*b**2*dn-104.0_cp*a*b**2-64.0_cp*b**3*dn+64.0_cp*b**3 &
      &             )/(64.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku-1)=a**4*(2.0_cp*a**4*dn-3.0_cp*a**4+8.0_cp*a**3*b*dn-24.0_cp &
      &             *a**3*b+6.0_cp*a**2*b**2*dn-66.0_cp*a**2*b**2+128.0_cp*a*b**3*dn+192.0_cp &
      &             *a*b**3+128.0_cp*b**4*dn-192.0_cp*b**4)/(128.0_cp*dn*(dn-1.0_cp)*(dn &
      &             +1.0_cp))
      stencil(ku)=a**4*b*(a**3-80.0_cp*a*b**2+64.0_cp*b**3)/(64.0_cp*dn*(dn &
      &             -1.0_cp))
      stencil(ku+1)=-a**4*(5.0_cp*a**4+32.0_cp*a**3*b+72.0_cp*a**2*b**2-64.0_cp &
      &             *a*b**3+320.0_cp*b**4)/(128.0_cp*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+2)=a**4*b*(a**3-80.0_cp*a*b**2+64.0_cp*b**3)/(64.0_cp*dn &
      &             *(dn+1.0_cp))
      stencil(ku+3)=a**4*(2.0_cp*a**4*dn+3.0_cp*a**4+8.0_cp*a**3*b*dn+24.0_cp &
      &             *a**3*b+6.0_cp*a**2*b**2*dn+66.0_cp*a**2*b**2+128.0_cp*a*b**3*dn-192.0_cp &
      &             *a*b**3+128.0_cp*b**4*dn+192.0_cp*b**4)/(128.0_cp*dn*(dn-1.0_cp)*(dn &
      &             +1.0_cp))
      stencil(ku+4)=a**4*b*(3.0_cp*a**3*dn-5.0_cp*a**3+32.0_cp*a**2*b*dn-32.0_cp &
      &             *a**2*b+56.0_cp*a*b**2*dn+104.0_cp*a*b**2-64.0_cp*b**3*dn-64.0_cp*b**3 &
      &             )/(64.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+5)=a**4*(2.0_cp*a**4*dn-a**4+16.0_cp*a**3*b*dn+4.0_cp*a**2* &
      &             b**2*dn+56.0_cp*a**2*b**2-160.0_cp*a*b**3*dn-160.0_cp*a*b**3+32.0_cp &
      &             *b**4*dn+32.0_cp*b**4)/(128.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+6)=-a**5*b*(7.0_cp*a**2*dn+a**2+32.0_cp*a*b*dn+32.0_cp*a &
      &             *b-24.0_cp*b**2*dn-24.0_cp*b**2)/(64.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp &
      &             ))
      stencil(ku+7)=-a**6*(2.0_cp*a**2*dn+a**2+8.0_cp*a*b*dn+8.0_cp*a*b-26.0_cp &
      &             *b**2*dn-26.0_cp*b**2)/(128.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+8)=3.0_cp*a**7*b/(64.0_cp*dn*(dn-1.0_cp))
      stencil(ku+9)=a**8/(256.0_cp*dn*(dn-1.0_cp))
      stencil(ku+10:len_stencil) = 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb2rmult2hmult4
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
   function intcheb2rmult2hmult2lapl(a, b, m, n, len_stencil) result(stencil)

      !-- Input variables
      real(cp), intent(in) :: a
      real(cp), intent(in) :: b
      integer,  intent(in) :: m
      integer,  intent(in) :: n
      integer,  intent(in) :: len_stencil

      !-- Output variable 
      real(cp) :: stencil(len_stencil)

      !-- Local variables
      real(cp) :: dn, dm
      integer :: ku

      dn = real(n, cp)
      dm = real(m, cp)
      ku = (len_stencil-1)/2

      stencil(1:ku-4) = 0.0_cp
      stencil(ku-3)=a**4*(dm-dn-4.0_cp)*(dm+dn+4.0_cp)/(16.0_cp*dn*(dn+1.0_cp))
      stencil(ku-2)=a**3*b*(2.0_cp*dm**2-4.0_cp*dn**2-25.0_cp*dn-39.0_cp) &
      &             /(8.0_cp*dn*(dn+1.0_cp))
      stencil(ku-1)=-a**2*(2.0_cp*a**2*dm**2*dn-a**2*dm**2+a**2*dn**2-2.0_cp   &
      &             *a**2*dn-8.0_cp*a**2+4.0_cp*a*b*dm**2*dn-4.0_cp*a*b*dm**2- &
      &             4.0_cp*a*b*dn**3-12.0_cp*a*b*dn**2+16.0_cp*a*b+10.0_cp*b**2&
      &             *dn**3+36.0_cp*b**2*dn**2+6.0_cp*b**2*dn-52.0_cp*b**2)/    &
      &             (8.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku)  =-a*b*(2.0_cp*a**2*dm**2+4.0_cp*a**2*dn**2+5.0_cp*a**2*dn   &
      &             -27.0_cp*a**2-16.0_cp*a*b*dn**2-8.0_cp*a*b*dn+24.0_cp*a*b+ &
      &             8.0_cp*b**2*dn**2+8.0_cp*b**2*dn-16.0_cp*b**2)/(8.0_cp*dn* &
      &             (dn-1.0_cp))
      stencil(ku+1)=a*(3.0_cp*a**3*dm**2+a**3*dn**2+4.0_cp*a**3+8.0_cp*a**2*   &
      &             b*dm**2+8.0_cp*a**2*b*dn**2-16.0_cp*a**2*b-12.0_cp*a*b**2* &
      &             dn**2+44.0_cp*a*b**2+16.0_cp*b**3*dn**2-16.0_cp*b**3)/     &
      &             (8.0_cp*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+2)=-a*b*(2.0_cp*a**2*dm**2+4.0_cp*a**2*dn**2-5.0_cp*a**2*dn-  &
      &             27.0_cp*a**2-16.0_cp*a*b*dn**2+8.0_cp*a*b*dn+24.0_cp*a*b+  &
      &             8.0_cp*b**2*dn**2-8.0_cp*b**2*dn-16.0_cp*b**2)/(8.0_cp*dn* &
      &             (dn+1.0_cp))
      stencil(ku+3)=-a**2*(2.0_cp*a**2*dm**2*dn+a**2*dm**2-a**2*dn**2-2.0_cp   &
      &             *a**2*dn+8.0_cp*a**2+4.0_cp*a*b*dm**2*dn+4.0_cp*a*b*dm**2- &
      &             4.0_cp*a*b*dn**3+12.0_cp*a*b*dn**2-16.0_cp*a*b+10.0_cp*b**2&
      &             *dn**3-36.0_cp*b**2*dn**2+6.0_cp*b**2*dn+52.0_cp*b**2)/    &
      &             (8.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+4)=a**3*b*(2.0_cp*dm**2-4.0_cp*dn**2+25.0_cp*dn-39.0_cp) &
      &             /(8.0_cp*dn*(dn-1.0_cp))
      stencil(ku+5)=a**4*(dm-dn+4.0_cp)*(dm+dn-4.0_cp)/(16.0_cp*dn*(dn-1.0_cp))
      stencil(ku+6:len_stencil) = 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb2rmult2hmult2lapl
!------------------------------------------------------------------------------
   function intcheb2rmult2hmult4laplaxi(a, b, n, len_stencil) result(stencil)

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

      stencil(1:ku-6) = 0.0_cp
      stencil(ku-5)=a**6*(dn+5.0_cp)*(dn+7.0_cp)/(64.0_cp*dn*(dn+1.0_cp))
      stencil(ku-4)=a**5*b*(6.0_cp*dn**2+61.0_cp*dn+151.0_cp)/(32.0_cp*dn &
      &             *(dn+1.0_cp))
      stencil(ku-3)=-a**4*(a**2*dn**3+6.0_cp*a**2*dn**2+7.0_cp*a**2*dn+10.0_cp &
      &             *a**2+8.0_cp*a*b*dn**3+56.0_cp*a*b*dn**2+56.0_cp*a*b*dn-120.0_cp*a* &
      &             b-26.0_cp*b**2*dn**3-192.0_cp*b**2*dn**2-230.0_cp*b**2*dn+448.0_cp* &
      &             b**2)/(32.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku-2)=-a**3*b*(2.0_cp*a**2*dn**3+3.0_cp*a**2*dn**2+30.0_cp* &
      &             a**2*dn+157.0_cp*a**2+64.0_cp*a*b*dn**3+336.0_cp*a*b*dn**2+192.0_cp &
      &             *a*b*dn-592.0_cp*a*b-48.0_cp*b**2*dn**3-272.0_cp*b**2*dn**2-208.0_cp &
      &             *b**2*dn+528.0_cp*b**2)/(32.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku-1)=-a**2*(a**4*dn**3+9.0_cp*a**4*dn**2+23.0_cp*a**4*dn-17.0_cp &
      &             *a**4+32.0_cp*a**3*b*dn**2-288.0_cp*a**3*b-112.0_cp*a**2*b**2*dn**3 &
      &             -384.0_cp*a**2*b**2*dn**2+368.0_cp*a**2*b**2*dn+1152.0_cp*a**2*b**2 &
      &             +320.0_cp*a*b**3*dn**3+1152.0_cp*a*b**3*dn**2+192.0_cp*a*b**3*dn-1664.0_cp &
      &             *a*b**3-64.0_cp*b**4*dn**3-256.0_cp*b**4*dn**2-64.0_cp*b**4*dn+384.0_cp &
      &             *b**4)/(64.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku)=-a**2*b*(2.0_cp*a**3*dn**2+3.0_cp*a**3*dn+27.0_cp*a**3+32.0_cp &
      &             *a**2*b*dn**2+40.0_cp*a**2*b*dn-200.0_cp*a**2*b-104.0_cp*a*b**2*dn**2- &
      &             64.0_cp*a*b**2*dn+296.0_cp*a*b**2+64.0_cp*b**3*dn**2+64.0_cp*b**3*dn &
      &             -128.0_cp*b**3)/(16.0_cp*dn*(dn-1.0_cp))
      stencil(ku+1)=a**2*(a**4*dn**2+7.0_cp*a**4+8.0_cp*a**3*b*dn**2+56.0_cp &
      &             *a**3*b+46.0_cp*a**2*b**2*dn**2-174.0_cp*a**2*b**2-96.0_cp*a*b**3*dn**2 &
      &             +352.0_cp*a*b**3+96.0_cp*b**4*dn**2-160.0_cp*b**4)/(16.0_cp*(dn-1.0_cp &
      &             )*(dn+1.0_cp))
      stencil(ku+2)=-a**2*b*(2.0_cp*a**3*dn**2-3.0_cp*a**3*dn+27.0_cp*a** &
      &             3+32.0_cp*a**2*b*dn**2-40.0_cp*a**2*b*dn-200.0_cp*a**2*b-104.0_cp*a &
      &             *b**2*dn**2+64.0_cp*a*b**2*dn+296.0_cp*a*b**2+64.0_cp*b**3*dn**2-64.0_cp &
      &             *b**3*dn-128.0_cp*b**3)/(16.0_cp*dn*(dn+1.0_cp))
      stencil(ku+3)=-a**2*(a**4*dn**3-9.0_cp*a**4*dn**2+23.0_cp*a**4*dn+17.0_cp &
      &             *a**4-32.0_cp*a**3*b*dn**2+288.0_cp*a**3*b-112.0_cp*a**2*b**2*dn**3 &
      &             +384.0_cp*a**2*b**2*dn**2+368.0_cp*a**2*b**2*dn-1152.0_cp*a**2*b**2 &
      &             +320.0_cp*a*b**3*dn**3-1152.0_cp*a*b**3*dn**2+192.0_cp*a*b**3*dn+1664.0_cp &
      &             *a*b**3-64.0_cp*b**4*dn**3+256.0_cp*b**4*dn**2-64.0_cp*b**4*dn-384.0_cp &
      &             *b**4)/(64.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+4)=-a**3*b*(2.0_cp*a**2*dn**3-3.0_cp*a**2*dn**2+30.0_cp* &
      &             a**2*dn-157.0_cp*a**2+64.0_cp*a*b*dn**3-336.0_cp*a*b*dn**2+192.0_cp &
      &             *a*b*dn+592.0_cp*a*b-48.0_cp*b**2*dn**3+272.0_cp*b**2*dn**2-208.0_cp &
      &             *b**2*dn-528.0_cp*b**2)/(32.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+5)=-a**4*(a**2*dn**3-6.0_cp*a**2*dn**2+7.0_cp*a**2*dn-10.0_cp &
      &             *a**2+8.0_cp*a*b*dn**3-56.0_cp*a*b*dn**2+56.0_cp*a*b*dn+120.0_cp*a* &
      &             b-26.0_cp*b**2*dn**3+192.0_cp*b**2*dn**2-230.0_cp*b**2*dn-448.0_cp* &
      &             b**2)/(32.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+6)=a**5*b*(6.0_cp*dn**2-61.0_cp*dn+151.0_cp)/(32.0_cp*dn &
      &             *(dn-1.0_cp))
      stencil(ku+7)=a**6*(dn-7.0_cp)*(dn-5.0_cp)/(64.0_cp*dn*(dn-1.0_cp))
      stencil(ku+8:len_stencil) = 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb2rmult2hmult4laplaxi
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
   function intcheb4hmult2(a, b, n, len_stencil) result(stencil)

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

      stencil(1:ku-6) = 0.0_cp
      stencil(ku-5)=-a**6/(64.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-4)=-a**5*b/(16.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-3)=a**5*(3.0_cp*a*dn+3.0_cp*a+4.0_cp*b*dn-4.0_cp*b)/(32.0_cp &
      &             *dn*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-2)=3.0_cp*a**5*b/(16.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku-1)=-a**5*(15.0_cp*a*dn-15.0_cp*a+32.0_cp*b*dn-64.0_cp*b) &
      &             /(64.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+3.0_cp))
      stencil(ku)  =-a**5*b*(dn+8.0_cp)/(8.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)* &
      &             (dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+1)=a**5*(5.0_cp*a*dn**2-35.0_cp*a+12.0_cp*b*dn**2-108.0_cp &
      &             *b)/(16.0_cp*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*       &
      &             (dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+2)=-a**5*b*(dn-8.0_cp)/(8.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp &
      &             )*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+3)=-a**5*(15.0_cp*a*dn+15.0_cp*a+32.0_cp*b*dn+64.0_cp*b) &
      &             /(64.0_cp*dn*(dn-3.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+4)=3.0_cp*a**5*b/(16.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+5)=a**5*(3.0_cp*a*dn-3.0_cp*a+4.0_cp*b*dn+4.0_cp*b)/(32.0_cp &
      &             *dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+6)=-a**5*b/(16.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp))
      stencil(ku+7)=-a**6/(64.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp))
      stencil(ku+8:len_stencil) = 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4hmult2
!------------------------------------------------------------------------------
   function intcheb4hmult2der1(a, b, n, len_stencil) result(stencil)

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

      stencil(1:ku-5) = 0.0_cp
      stencil(ku-4)=a**5*(dn+5.0_cp)/(32.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)* &
      &             (dn+3.0_cp))
      stencil(ku-3)=a**4*b*(dn+4.0_cp)/(8.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp) &
      &             *(dn+3.0_cp))
      stencil(ku-2)=-a**4*(5.0_cp*a*dn+7.0_cp*a+8.0_cp*b*dn-8.0_cp*b)/(32.0_cp &
      &             *dn*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku-1)=-a**4*b*(dn+5.0_cp)/(4.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp &
      &             )*(dn+3.0_cp))
      stencil(ku)  =a**4*(5.0_cp*a*dn**2+9.0_cp*a*dn-8.0_cp*a+12.0_cp*b*dn**2+   &
      &             12.0_cp*b*dn-72.0_cp*b)/(16.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)* &
      &             (dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+1)=3.0_cp*a**4*b/(2.0_cp*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)* &
      &             (dn+2.0_cp))
      stencil(ku+2)=-a**4*(5.0_cp*a*dn**2-9.0_cp*a*dn-8.0_cp*a+12.0_cp*b*dn**2   &
      &             -12.0_cp*b*dn-72.0_cp*b)/(16.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*&
      &             (dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+3)=a**4*b*(dn-5.0_cp)/(4.0_cp*dn*(dn-3.0_cp)*(dn-1.0_cp) &
      &             *(dn+1.0_cp))
      stencil(ku+4)=a**4*(5.0_cp*a*dn-7.0_cp*a+8.0_cp*b*dn+8.0_cp*b)/(32.0_cp &
      &             *dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+5)=-a**4*b*(dn-4.0_cp)/(8.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp) &
      &             *(dn-1.0_cp))
      stencil(ku+6)=-a**5*(dn-5.0_cp)/(32.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)   &
      &             *(dn-1.0_cp))
      stencil(ku+7:len_stencil) = 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4hmult2der1
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
   function intcheb4rmult4hmult6(a, b, n, len_stencil) result(stencil)
      !
      ! Fourth integral of r^4*(ro^2-r^2)^3 operator
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

      dn = real(n, cp)
      ku = (len_stencil-1)/2

      stencil(1:ku-14) = 0.0_cp
      stencil(ku-13)=-a**14/(16384.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku-12)=-5.0_cp*a**13*b/(4096.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp &
      &              )*(dn+3.0_cp))
      stencil(ku-11)=3.0_cp*a**12*(a**2*dn+a**2+4.0_cp*a*b*dn-4.0_cp*a*b-28.0_cp &
      &              *b**2*dn+28.0_cp*b**2)/(8192.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp &
      &              )*(dn+3.0_cp))
      stencil(ku-10)=a**11*b*(23.0_cp*a**2*dn+37.0_cp*a**2+96.0_cp*a*b*dn &
      &              -96.0_cp*a*b-192.0_cp*b**2*dn+192.0_cp*b**2)/(4096.0_cp*dn*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-9)=-a**10*(11.0_cp*a**4*dn**2+27.0_cp*a**4*dn-38.0_cp*a**4+ &
      &              96.0_cp*a**3*b*dn**2-384.0_cp*a**3*b-384.0_cp*a**2*b**2*dn**2-864.0_cp &
      &              *a**2*b**2*dn+3264.0_cp*a**2*b**2-2496.0_cp*a*b**3*dn**2+7488.0_cp* &
      &              a*b**3*dn-4992.0_cp*a*b**3+2064.0_cp*b**4*dn**2-6192.0_cp*b**4*dn+4128.0_cp &
      &              *b**4)/(16384.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp &
      &              )*(dn+3.0_cp))
      stencil(ku-8)=-3.0_cp*a**9*b*(5.0_cp*a**4*dn**2+21.0_cp*a**4*dn-12.0_cp &
      &              *a**4+48.0_cp*a**3*b*dn**2+48.0_cp*a**3*b*dn-288.0_cp*a**3*b+16.0_cp &
      &              *a**2*b**2*dn**2-432.0_cp*a**2*b**2*dn+800.0_cp*a**2*b**2-352.0_cp* &
      &              a*b**3*dn**2+1056.0_cp*a*b**3*dn-704.0_cp*a*b**3+136.0_cp*b**4*dn** &
      &              2-408.0_cp*b**4*dn+272.0_cp*b**4)/(2048.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-7)=-a**8*(a**6*dn**3-24.0_cp*a**6*dn**2+41.0_cp*a**6*dn+36.0_cp &
      &              *a**6-12.0_cp*a**5*b*dn**3-144.0_cp*a**5*b*dn**2+588.0_cp*a**5*b*dn &
      &              -144.0_cp*a**5*b+36.0_cp*a**4*b**2*dn**3+432.0_cp*a**4*b**2*dn**2-324.0_cp &
      &              *a**4*b**2*dn-3888.0_cp*a**4*b**2+1312.0_cp*a**3*b**3*dn**3-384.0_cp &
      &              *a**3*b**3*dn**2-23008.0_cp*a**3*b**3*dn+37056.0_cp*a**3*b**3+1752.0_cp &
      &              *a**2*b**4*dn**3-16704.0_cp*a**2*b**4*dn**2+50232.0_cp*a**2*b**4*dn &
      &              -47664.0_cp*a**2*b**4-3936.0_cp*a*b**5*dn**3+23616.0_cp*a*b**5*dn** &
      &              2-43296.0_cp*a*b**5*dn+23616.0_cp*a*b**5+704.0_cp*b**6*dn**3-4224.0_cp &
      &              *b**6*dn**2+7744.0_cp*b**6*dn-4224.0_cp*b**6)/(4096.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-6)=-a**7*b*(7.0_cp*a**6*dn**3-114.0_cp*a**6*dn**2+47.0_cp &
      &              *a**6*dn+396.0_cp*a**6-48.0_cp*a**5*b*dn**3-864.0_cp*a**5*b*dn**2+2352.0_cp &
      &              *a**5*b*dn+2016.0_cp*a**5*b-336.0_cp*a**4*b**2*dn**3+288.0_cp*a**4* &
      &              b**2*dn**2+10704.0_cp*a**4*b**2*dn-25632.0_cp*a**4*b**2+1504.0_cp*a**3 &
      &              *b**3*dn**3+3648.0_cp*a**3*b**3*dn**2-46816.0_cp*a**3*b**3*dn+67008.0_cp &
      &              *a**3*b**3+2520.0_cp*a**2*b**4*dn**3-20016.0_cp*a**2*b**4*dn**2+52200.0_cp &
      &              *a**2*b**4*dn-44496.0_cp*a**2*b**4-1920.0_cp*a*b**5*dn**3+11520.0_cp &
      &              *a*b**5*dn**2-21120.0_cp*a*b**5*dn+11520.0_cp*a*b**5+128.0_cp*b**6*dn**3 &
      &              -768.0_cp*b**6*dn**2+1408.0_cp*b**6*dn-768.0_cp*b**6)/(2048.0_cp*dn &
      &              *(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku-5)=3.0_cp*a**7*(13.0_cp*a**7*dn**3-38.0_cp*a**7*dn**2-97.0_cp &
      &              *a**7*dn+162.0_cp*a**7+96.0_cp*a**6*b*dn**3-480.0_cp*a**6*b*dn**2-384.0_cp &
      &              *a**6*b*dn+1920.0_cp*a**6*b+128.0_cp*a**5*b**2*dn**3-1056.0_cp*a**5 &
      &              *b**2*dn**2+3808.0_cp*a**5*b**2*dn+1344.0_cp*a**5*b**2+576.0_cp*a** &
      &              4*b**3*dn**3+7552.0_cp*a**4*b**3*dn**2+1216.0_cp*a**4*b**3*dn-87168.0_cp &
      &              *a**4*b**3+5904.0_cp*a**3*b**4*dn**3+864.0_cp*a**3*b**4*dn**2-157776.0_cp &
      &              *a**3*b**4*dn+306144.0_cp*a**3*b**4-5120.0_cp*a**2*b**5*dn**3-32256.0_cp &
      &              *a**2*b**5*dn**2+258560.0_cp*a**2*b**5*dn-347136.0_cp*a**2*b**5-8192.0_cp &
      &              *a*b**6*dn**3+60416.0_cp*a*b**6*dn**2-146432.0_cp*a*b**6*dn+116736.0_cp &
      &              *a*b**6+2048.0_cp*b**7*dn**3-12288.0_cp*b**7*dn**2+22528.0_cp*b**7*dn &
      &              -12288.0_cp*b**7)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-4)=a**7*b*(65.0_cp*a**6*dn**3-78.0_cp*a**6*dn**2-785.0_cp &
      &              *a**6*dn+342.0_cp*a**6+480.0_cp*a**5*b*dn**3-1728.0_cp*a**5*b*dn**2 &
      &              -6240.0_cp*a**5*b*dn+9792.0_cp*a**5*b+704.0_cp*a**4*b**2*dn**3-8832.0_cp &
      &              *a**4*b**2*dn**2+1984.0_cp*a**4*b**2*dn+77568.0_cp*a**4*b**2-256.0_cp &
      &              *a**3*b**3*dn**3+12288.0_cp*a**3*b**3*dn**2+70144.0_cp*a**3*b**3*dn &
      &              -314112.0_cp*a**3*b**3+12864.0_cp*a**2*b**4*dn**3-6912.0_cp*a**2*b**4* &
      &              dn**2-258816.0_cp*a**2*b**4*dn+491328.0_cp*a**2*b**4-2816.0_cp*a*b**5* &
      &              dn**3-29184.0_cp*a*b**5*dn**2+199424.0_cp*a*b**5*dn-259584.0_cp*a*b**5 &
      &              -2816.0_cp*b**6*dn**3+19968.0_cp*b**6*dn**2-46336.0_cp*b**6*dn+35328.0_cp &
      &              *b**6)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-3)=-a**7*(19.0_cp*a**7*dn**3+60.0_cp*a**7*dn**2-301.0_cp &
      &              *a**7*dn-90.0_cp*a**7+204.0_cp*a**6*b*dn**3+360.0_cp*a**6*b*dn**2-3516.0_cp &
      &              *a**6*b*dn+360.0_cp*a**6*b+492.0_cp*a**5*b**2*dn**3-216.0_cp*a**5*b**2 &
      &              *dn**2-11868.0_cp*a**5*b**2*dn+17064.0_cp*a**5*b**2-1280.0_cp*a**4* &
      &              b**3*dn**3+1536.0_cp*a**4*b**3*dn**2+24320.0_cp*a**4*b**3*dn+97536.0_cp &
      &              *a**4*b**3-3008.0_cp*a**3*b**4*dn**3+69888.0_cp*a**3*b**4*dn**2+41792.0_cp &
      &              *a**3*b**4*dn-796992.0_cp*a**3*b**4+13056.0_cp*a**2*b**5*dn**3-76032.0_cp &
      &              *a**2*b**5*dn**2-340224.0_cp*a**2*b**5*dn+1352448.0_cp*a**2*b**5-36352.0_cp &
      &              *a*b**6*dn**3+53760.0_cp*a*b**6*dn**2+506368.0_cp*a*b**6*dn-1021440.0_cp &
      &              *a*b**6+2048.0_cp*b**7*dn**3+24576.0_cp*b**7*dn**2-161792.0_cp*b**7 &
      &              *dn+208896.0_cp*b**7)/(8192.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-2)=-3.0_cp*a**7*b*(17.0_cp*a**6*dn**2+3.0_cp*a**6*dn-102.0_cp &
      &              *a**6+192.0_cp*a**5*b*dn**2-192.0_cp*a**5*b*dn-1152.0_cp*a**5*b+768.0_cp &
      &              *a**4*b**2*dn**2-2560.0_cp*a**4*b**2*dn-3712.0_cp*a**4*b**2+768.0_cp &
      &              *a**3*b**3*dn**2-11520.0_cp*a**3*b**3*dn+41728.0_cp*a**3*b**3+64.0_cp &
      &              *a**2*b**4*dn**2+27456.0_cp*a**2*b**4*dn-88384.0_cp*a**2*b**4+4864.0_cp &
      &              *a*b**5*dn**2-39680.0_cp*a*b**5*dn+75264.0_cp*a*b**5-3328.0_cp*b**6 &
      &              *dn**2+17664.0_cp*b**6*dn-23040.0_cp*b**6)/(4096.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku-1)=-a**7*(27.0_cp*a**7*dn**2-324.0_cp*a**7*dn+465.0_cp*a**7 &
      &              +192.0_cp*a**6*b*dn**2-3264.0_cp*a**6*b*dn+5760.0_cp*a**6*b+768.0_cp &
      &              *a**5*b**2*dn**2-12480.0_cp*a**5*b**2*dn+30528.0_cp*a**5*b**2+4224.0_cp &
      &              *a**4*b**3*dn**2-15360.0_cp*a**4*b**3*dn+61824.0_cp*a**4*b**3+19744.0_cp &
      &              *a**3*b**4*dn**2+17920.0_cp*a**3*b**4*dn-548256.0_cp*a**3*b**4+9216.0_cp &
      &              *a**2*b**5*dn**2-382464.0_cp*a**2*b**5*dn+1442304.0_cp*a**2*b**5+24576.0_cp &
      &              *a*b**6*dn**2+347136.0_cp*a*b**6*dn-1330176.0_cp*a*b**6+38912.0_cp* &
      &              b**7*dn**2-286720.0_cp*b**7*dn+509952.0_cp*b**7)/(16384.0_cp*dn*(dn &
      &              -3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+3.0_cp))
      stencil(ku)=a**7*b*(3.0_cp*a**6*dn**2+15.0_cp*a**6*dn+48.0_cp*a**6+48.0_cp &
      &              *a**5*b*dn**2+240.0_cp*a**5*b*dn+288.0_cp*a**5*b+304.0_cp*a**4*b**2 &
      &              *dn**2+1520.0_cp*a**4*b**2*dn-1536.0_cp*a**4*b**2+864.0_cp*a**3*b** &
      &              3*dn**2+4320.0_cp*a**3*b**3*dn-45696.0_cp*a**3*b**3-1704.0_cp*a**2* &
      &              b**4*dn**2-8520.0_cp*a**2*b**4*dn+97056.0_cp*a**2*b**4+3392.0_cp*a* &
      &              b**5*dn**2+16960.0_cp*a*b**5*dn-110208.0_cp*a*b**5-1728.0_cp*b**6*dn**2 &
      &              -8640.0_cp*b**6*dn+43392.0_cp*b**6)/(1024.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp &
      &              )*(dn-1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+1)=3.0_cp*a**7*(3.0_cp*a**7*dn**2-37.0_cp*a**7+28.0_cp*a**6 &
      &              *b*dn**2-412.0_cp*a**6*b+108.0_cp*a**5*b**2*dn**2-1932.0_cp*a**5*b**2+ &
      &              224.0_cp*a**4*b**3*dn**2-3616.0_cp*a**4*b**3+424.0_cp*a**3*b**4*dn**2+ &
      &              22344.0_cp*a**3*b**4+2912.0_cp*a**2*b**5*dn**2-79328.0_cp*a**2*b**5 &
      &              -1728.0_cp*a*b**6*dn**2+70592.0_cp*a*b**6+3072.0_cp*b**7*dn**2-37888.0_cp &
      &              *b**7)/(2048.0_cp*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn &
      &              +2.0_cp)*(dn+3.0_cp))
      stencil(ku+2)=a**7*b*(3.0_cp*a**6*dn**2-15.0_cp*a**6*dn+48.0_cp*a** &
      &              6+48.0_cp*a**5*b*dn**2-240.0_cp*a**5*b*dn+288.0_cp*a**5*b+304.0_cp* &
      &              a**4*b**2*dn**2-1520.0_cp*a**4*b**2*dn-1536.0_cp*a**4*b**2+864.0_cp &
      &              *a**3*b**3*dn**2-4320.0_cp*a**3*b**3*dn-45696.0_cp*a**3*b**3-1704.0_cp &
      &              *a**2*b**4*dn**2+8520.0_cp*a**2*b**4*dn+97056.0_cp*a**2*b**4+3392.0_cp &
      &              *a*b**5*dn**2-16960.0_cp*a*b**5*dn-110208.0_cp*a*b**5-1728.0_cp*b** &
      &              6*dn**2+8640.0_cp*b**6*dn+43392.0_cp*b**6)/(1024.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+3)=-a**7*(27.0_cp*a**7*dn**2+324.0_cp*a**7*dn+465.0_cp*a**7 &
      &              +192.0_cp*a**6*b*dn**2+3264.0_cp*a**6*b*dn+5760.0_cp*a**6*b+768.0_cp &
      &              *a**5*b**2*dn**2+12480.0_cp*a**5*b**2*dn+30528.0_cp*a**5*b**2+4224.0_cp &
      &              *a**4*b**3*dn**2+15360.0_cp*a**4*b**3*dn+61824.0_cp*a**4*b**3+19744.0_cp &
      &              *a**3*b**4*dn**2-17920.0_cp*a**3*b**4*dn-548256.0_cp*a**3*b**4+9216.0_cp &
      &              *a**2*b**5*dn**2+382464.0_cp*a**2*b**5*dn+1442304.0_cp*a**2*b**5+24576.0_cp &
      &              *a*b**6*dn**2-347136.0_cp*a*b**6*dn-1330176.0_cp*a*b**6+38912.0_cp* &
      &              b**7*dn**2+286720.0_cp*b**7*dn+509952.0_cp*b**7)/(16384.0_cp*dn*(dn &
      &              -3.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+4)=-3.0_cp*a**7*b*(17.0_cp*a**6*dn**2-3.0_cp*a**6*dn-102.0_cp &
      &              *a**6+192.0_cp*a**5*b*dn**2+192.0_cp*a**5*b*dn-1152.0_cp*a**5*b+768.0_cp &
      &              *a**4*b**2*dn**2+2560.0_cp*a**4*b**2*dn-3712.0_cp*a**4*b**2+768.0_cp &
      &              *a**3*b**3*dn**2+11520.0_cp*a**3*b**3*dn+41728.0_cp*a**3*b**3+64.0_cp &
      &              *a**2*b**4*dn**2-27456.0_cp*a**2*b**4*dn-88384.0_cp*a**2*b**4+4864.0_cp &
      &              *a*b**5*dn**2+39680.0_cp*a*b**5*dn+75264.0_cp*a*b**5-3328.0_cp*b**6 &
      &              *dn**2-17664.0_cp*b**6*dn-23040.0_cp*b**6)/(4096.0_cp*dn*(dn-2.0_cp &
      &              )*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+5)=-a**7*(19.0_cp*a**7*dn**3-60.0_cp*a**7*dn**2-301.0_cp &
      &              *a**7*dn+90.0_cp*a**7+204.0_cp*a**6*b*dn**3-360.0_cp*a**6*b*dn**2-3516.0_cp &
      &              *a**6*b*dn-360.0_cp*a**6*b+492.0_cp*a**5*b**2*dn**3+216.0_cp*a**5*b**2 &
      &              *dn**2-11868.0_cp*a**5*b**2*dn-17064.0_cp*a**5*b**2-1280.0_cp*a**4* &
      &              b**3*dn**3-1536.0_cp*a**4*b**3*dn**2+24320.0_cp*a**4*b**3*dn-97536.0_cp &
      &              *a**4*b**3-3008.0_cp*a**3*b**4*dn**3-69888.0_cp*a**3*b**4*dn**2+41792.0_cp &
      &              *a**3*b**4*dn+796992.0_cp*a**3*b**4+13056.0_cp*a**2*b**5*dn**3+76032.0_cp &
      &              *a**2*b**5*dn**2-340224.0_cp*a**2*b**5*dn-1352448.0_cp*a**2*b**5-36352.0_cp &
      &              *a*b**6*dn**3-53760.0_cp*a*b**6*dn**2+506368.0_cp*a*b**6*dn+1021440.0_cp &
      &              *a*b**6+2048.0_cp*b**7*dn**3-24576.0_cp*b**7*dn**2-161792.0_cp*b**7 &
      &              *dn-208896.0_cp*b**7)/(8192.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+6)=a**7*b*(65.0_cp*a**6*dn**3+78.0_cp*a**6*dn**2-785.0_cp &
      &              *a**6*dn-342.0_cp*a**6+480.0_cp*a**5*b*dn**3+1728.0_cp*a**5*b*dn**2 &
      &              -6240.0_cp*a**5*b*dn-9792.0_cp*a**5*b+704.0_cp*a**4*b**2*dn**3+8832.0_cp &
      &              *a**4*b**2*dn**2+1984.0_cp*a**4*b**2*dn-77568.0_cp*a**4*b**2-256.0_cp &
      &              *a**3*b**3*dn**3-12288.0_cp*a**3*b**3*dn**2+70144.0_cp*a**3*b**3*dn &
      &              +314112.0_cp*a**3*b**3+12864.0_cp*a**2*b**4*dn**3+6912.0_cp*a**2*b**4* &
      &              dn**2-258816.0_cp*a**2*b**4*dn-491328.0_cp*a**2*b**4-2816.0_cp*a*b**5* &
      &              dn**3+29184.0_cp*a*b**5*dn**2+199424.0_cp*a*b**5*dn+259584.0_cp*a*b**5 &
      &              -2816.0_cp*b**6*dn**3-19968.0_cp*b**6*dn**2-46336.0_cp*b**6*dn-35328.0_cp &
      &              *b**6)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+7)=3.0_cp*a**7*(13.0_cp*a**7*dn**3+38.0_cp*a**7*dn**2-97.0_cp &
      &              *a**7*dn-162.0_cp*a**7+96.0_cp*a**6*b*dn**3+480.0_cp*a**6*b*dn**2-384.0_cp &
      &              *a**6*b*dn-1920.0_cp*a**6*b+128.0_cp*a**5*b**2*dn**3+1056.0_cp*a**5 &
      &              *b**2*dn**2+3808.0_cp*a**5*b**2*dn-1344.0_cp*a**5*b**2+576.0_cp*a** &
      &              4*b**3*dn**3-7552.0_cp*a**4*b**3*dn**2+1216.0_cp*a**4*b**3*dn+87168.0_cp &
      &              *a**4*b**3+5904.0_cp*a**3*b**4*dn**3-864.0_cp*a**3*b**4*dn**2-157776.0_cp &
      &              *a**3*b**4*dn-306144.0_cp*a**3*b**4-5120.0_cp*a**2*b**5*dn**3+32256.0_cp &
      &              *a**2*b**5*dn**2+258560.0_cp*a**2*b**5*dn+347136.0_cp*a**2*b**5-8192.0_cp &
      &              *a*b**6*dn**3-60416.0_cp*a*b**6*dn**2-146432.0_cp*a*b**6*dn-116736.0_cp &
      &              *a*b**6+2048.0_cp*b**7*dn**3+12288.0_cp*b**7*dn**2+22528.0_cp*b**7*dn &
      &              +12288.0_cp*b**7)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+8)=-a**7*b*(7.0_cp*a**6*dn**3+114.0_cp*a**6*dn**2+47.0_cp &
      &              *a**6*dn-396.0_cp*a**6-48.0_cp*a**5*b*dn**3+864.0_cp*a**5*b*dn**2+2352.0_cp &
      &              *a**5*b*dn-2016.0_cp*a**5*b-336.0_cp*a**4*b**2*dn**3-288.0_cp*a**4* &
      &              b**2*dn**2+10704.0_cp*a**4*b**2*dn+25632.0_cp*a**4*b**2+1504.0_cp*a**3 &
      &              *b**3*dn**3-3648.0_cp*a**3*b**3*dn**2-46816.0_cp*a**3*b**3*dn-67008.0_cp &
      &              *a**3*b**3+2520.0_cp*a**2*b**4*dn**3+20016.0_cp*a**2*b**4*dn**2+52200.0_cp &
      &              *a**2*b**4*dn+44496.0_cp*a**2*b**4-1920.0_cp*a*b**5*dn**3-11520.0_cp &
      &              *a*b**5*dn**2-21120.0_cp*a*b**5*dn-11520.0_cp*a*b**5+128.0_cp*b**6*dn**3 &
      &              +768.0_cp*b**6*dn**2+1408.0_cp*b**6*dn+768.0_cp*b**6)/(2048.0_cp*dn &
      &              *(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku+9)=-a**8*(a**6*dn**3+24.0_cp*a**6*dn**2+41.0_cp*a**6*dn-36.0_cp &
      &              *a**6-12.0_cp*a**5*b*dn**3+144.0_cp*a**5*b*dn**2+588.0_cp*a**5*b*dn &
      &              +144.0_cp*a**5*b+36.0_cp*a**4*b**2*dn**3-432.0_cp*a**4*b**2*dn**2-324.0_cp &
      &              *a**4*b**2*dn+3888.0_cp*a**4*b**2+1312.0_cp*a**3*b**3*dn**3+384.0_cp &
      &              *a**3*b**3*dn**2-23008.0_cp*a**3*b**3*dn-37056.0_cp*a**3*b**3+1752.0_cp &
      &              *a**2*b**4*dn**3+16704.0_cp*a**2*b**4*dn**2+50232.0_cp*a**2*b**4*dn &
      &              +47664.0_cp*a**2*b**4-3936.0_cp*a*b**5*dn**3-23616.0_cp*a*b**5*dn** &
      &              2-43296.0_cp*a*b**5*dn-23616.0_cp*a*b**5+704.0_cp*b**6*dn**3+4224.0_cp &
      &              *b**6*dn**2+7744.0_cp*b**6*dn+4224.0_cp*b**6)/(4096.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+10)=-3.0_cp*a**9*b*(5.0_cp*a**4*dn**2-21.0_cp*a**4*dn-12.0_cp &
      &              *a**4+48.0_cp*a**3*b*dn**2-48.0_cp*a**3*b*dn-288.0_cp*a**3*b+16.0_cp &
      &              *a**2*b**2*dn**2+432.0_cp*a**2*b**2*dn+800.0_cp*a**2*b**2-352.0_cp* &
      &              a*b**3*dn**2-1056.0_cp*a*b**3*dn-704.0_cp*a*b**3+136.0_cp*b**4*dn** &
      &              2+408.0_cp*b**4*dn+272.0_cp*b**4)/(2048.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp &
      &              )*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+11)=-a**10*(11.0_cp*a**4*dn**2-27.0_cp*a**4*dn-38.0_cp*a**4 &
      &              +96.0_cp*a**3*b*dn**2-384.0_cp*a**3*b-384.0_cp*a**2*b**2*dn**2+864.0_cp &
      &              *a**2*b**2*dn+3264.0_cp*a**2*b**2-2496.0_cp*a*b**3*dn**2-7488.0_cp* &
      &              a*b**3*dn-4992.0_cp*a*b**3+2064.0_cp*b**4*dn**2+6192.0_cp*b**4*dn+4128.0_cp &
      &              *b**4)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp))
      stencil(ku+12)=a**11*b*(23.0_cp*a**2*dn-37.0_cp*a**2+96.0_cp*a*b*dn &
      &              +96.0_cp*a*b-192.0_cp*b**2*dn-192.0_cp*b**2)/(4096.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+13)=3.0_cp*a**12*(a**2*dn-a**2+4.0_cp*a*b*dn+4.0_cp*a*b-28.0_cp &
      &              *b**2*dn-28.0_cp*b**2)/(8192.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp))
      stencil(ku+14)=-5.0_cp*a**13*b/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp &
      &              )*(dn-1.0_cp))
      stencil(ku+15)=-a**14/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              ))
      stencil(ku+16:len_stencil) = 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4rmult4hmult6
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

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4rmult4laplrot
!------------------------------------------------------------------------------
   function intcheb4rmult4hmult8laplrot(a, b, m, n, len_stencil) result(stencil)

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

      stencil(1:ku-14) = 0.0_cp
      stencil(ku-13)=-a**14*(dm**2-dn**2-27.0_cp*dn-182.0_cp)/(16384.0_cp &
      &              *dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-12)=-a**13*b*(5.0_cp*dm**2-6.0_cp*dn**2-151.0_cp*dn-949.0_cp &
      &              )/(4096.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-11)=a**12*(5.0_cp*a**2*dm**2*dn+a**2*dm**2-3.0_cp*a**2*dn**3 &
      &              -71.0_cp*a**2*dn**2-458.0_cp*a**2*dn-404.0_cp*a**2+16.0_cp*a*b*dm** &
      &              2*dn-16.0_cp*a*b*dm**2-16.0_cp*a*b*dn**3-356.0_cp*a*b*dn**2-1796.0_cp &
      &              *a*b*dn+2168.0_cp*a*b-82.0_cp*b**2*dm**2*dn+82.0_cp*b**2*dm**2+124.0_cp &
      &              *b**2*dn**3+2770.0_cp*b**2*dn**2+13974.0_cp*b**2*dn-16868.0_cp*b**2 &
      &              )/(8192.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-10)=a**11*b*(39.0_cp*a**2*dm**2*dn+21.0_cp*a**2*dm**2-26.0_cp &
      &              *a**2*dn**3-583.0_cp*a**2*dn**2-3780.0_cp*a**2*dn-5259.0_cp*a**2+128.0_cp &
      &              *a*b*dm**2*dn-128.0_cp*a*b*dm**2-160.0_cp*a*b*dn**3-3268.0_cp*a*b*dn**2 &
      &              -14984.0_cp*a*b*dn+18412.0_cp*a*b-176.0_cp*b**2*dm**2*dn+176.0_cp*b**2 &
      &              *dm**2+360.0_cp*b**2*dn**3+7386.0_cp*b**2*dn**2+33868.0_cp*b**2*dn-41614.0_cp &
      &              *b**2)/(4096.0_cp*dn*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku-9)=-a**10*(43.0_cp*a**4*dm**2*dn**2-21.0_cp*a**4*dm**2*dn &
      &              -70.0_cp*a**4*dm**2-11.0_cp*a**4*dn**4-222.0_cp*a**4*dn**3-1375.0_cp &
      &              *a**4*dn**2-1908.0_cp*a**4*dn+3348.0_cp*a**4+256.0_cp*a**3*b*dm**2*dn**2 &
      &              -384.0_cp*a**3*b*dm**2*dn-256.0_cp*a**3*b*dm**2-128.0_cp*a**3*b*dn**4- &
      &              2384.0_cp*a**3*b*dn**3-10960.0_cp*a**3*b*dn**2+8000.0_cp*a**3*b*dn+48960.0_cp &
      &              *a**3*b-560.0_cp*a**2*b**2*dm**2*dn**2-288.0_cp*a**2*b**2*dm**2*dn+2816.0_cp &
      &              *a**2*b**2*dm**2+320.0_cp*a**2*b**2*dn**4+7544.0_cp*a**2*b**2*dn**3 &
      &              +59480.0_cp*a**2*b**2*dn**2+102304.0_cp*a**2*b**2*dn-508000.0_cp*a**2* &
      &              b**2-3200.0_cp*a*b**3*dm**2*dn**2+9600.0_cp*a*b**3*dm**2*dn-6400.0_cp &
      &              *a*b**3*dm**2+5376.0_cp*a*b**3*dn**4+89280.0_cp*a*b**3*dn**3+212544.0_cp &
      &              *a*b**3*dn**2-1343232.0_cp*a*b**3*dn+1036032.0_cp*a*b**3+1664.0_cp* &
      &              b**4*dm**2*dn**2-4992.0_cp*b**4*dm**2*dn+3328.0_cp*b**4*dm**2-5136.0_cp &
      &              *b**4*dn**4-85776.0_cp*b**4*dn**3-204128.0_cp*b**4*dn**2+1289856.0_cp &
      &              *b**4*dn-994816.0_cp*b**4)/(16384.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn &
      &              +1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-8)=-a**9*b*(63.0_cp*a**4*dm**2*dn**2+15.0_cp*a**4*dm**2*dn &
      &              -132.0_cp*a**4*dm**2-14.0_cp*a**4*dn**4-288.0_cp*a**4*dn**3-2122.0_cp &
      &              *a**4*dn**2-5430.0_cp*a**4*dn+1746.0_cp*a**4+384.0_cp*a**3*b*dm**2*dn**2 &
      &              -384.0_cp*a**3*b*dm**2*dn-768.0_cp*a**3*b*dm**2-208.0_cp*a**3*b*dn**4- &
      &              3750.0_cp*a**3*b*dn**3-19148.0_cp*a**3*b*dn**2-4050.0_cp*a**3*b*dn+118020.0_cp &
      &              *a**3*b+232.0_cp*a**2*b**2*dm**2*dn**2-1752.0_cp*a**2*b**2*dm**2*dn &
      &              +2576.0_cp*a**2*b**2*dm**2-492.0_cp*a**2*b**2*dn**4-5577.0_cp*a**2* &
      &              b**2*dn**3+17026.0_cp*a**2*b**2*dn**2+200709.0_cp*a**2*b**2*dn-417034.0_cp &
      &              *a**2*b**2-1216.0_cp*a*b**3*dm**2*dn**2+3648.0_cp*a*b**3*dm**2*dn-2432.0_cp &
      &              *a*b**3*dm**2+3072.0_cp*a*b**3*dn**4+45456.0_cp*a*b**3*dn**3+85856.0_cp &
      &              *a*b**3*dn**2-621840.0_cp*a*b**3*dn+487456.0_cp*a*b**3+256.0_cp*b** &
      &              4*dm**2*dn**2-768.0_cp*b**4*dm**2*dn+512.0_cp*b**4*dm**2-1440.0_cp* &
      &              b**4*dn**4-21444.0_cp*b**4*dn**3-40472.0_cp*b**4*dn**2+293124.0_cp* &
      &              b**4*dn-229768.0_cp*b**4)/(2048.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn &
      &              +1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-7)=a**8*(25.0_cp*a**6*dm**2*dn**3-48.0_cp*a**6*dm**2*dn**2- &
      &              115.0_cp*a**6*dm**2*dn+132.0_cp*a**6*dm**2+a**6*dn**5-213.0_cp*a**6 &
      &              *dn**3-936.0_cp*a**6*dn**2+2696.0_cp*a**6*dn+2544.0_cp*a**6+208.0_cp &
      &              *a**5*b*dm**2*dn**3-576.0_cp*a**5*b*dm**2*dn**2-592.0_cp*a**5*b*dm**2* &
      &              dn+1344.0_cp*a**5*b*dm**2-16.0_cp*a**5*b*dn**5-344.0_cp*a**5*b*dn** &
      &              4-2720.0_cp*a**5*b*dn**3-2152.0_cp*a**5*b*dn**2+43536.0_cp*a**5*b*dn &
      &              -6048.0_cp*a**5*b+46.0_cp*a**4*b**2*dm**2*dn**3-1464.0_cp*a**4*b**2 &
      &              *dm**2*dn**2+3986.0_cp*a**4*b**2*dm**2*dn-24.0_cp*a**4*b**2*dm**2-116.0_cp &
      &              *a**4*b**2*dn**5-988.0_cp*a**4*b**2*dn**4+6540.0_cp*a**4*b**2*dn**3 &
      &              +64156.0_cp*a**4*b**2*dn**2-9064.0_cp*a**4*b**2*dn-618576.0_cp*a**4 &
      &              *b**2-2880.0_cp*a**3*b**3*dm**2*dn**3+7680.0_cp*a**3*b**3*dm**2*dn**2+ &
      &              16320.0_cp*a**3*b**3*dm**2*dn-40320.0_cp*a**3*b**3*dm**2+1408.0_cp* &
      &              a**3*b**3*dn**5+25920.0_cp*a**3*b**3*dn**4+135360.0_cp*a**3*b**3*dn**3 &
      &              -332928.0_cp*a**3*b**3*dn**2-2889088.0_cp*a**3*b**3*dn+5567232.0_cp &
      &              *a**3*b**3-4032.0_cp*a**2*b**4*dm**2*dn**3+29184.0_cp*a**2*b**4*dm**2* &
      &              dn**2-69312.0_cp*a**2*b**4*dm**2*dn+54144.0_cp*a**2*b**4*dm**2+9432.0_cp &
      &              *a**2*b**4*dn**5+81408.0_cp*a**2*b**4*dn**4-355336.0_cp*a**2*b**4*dn**3 &
      &              -1956144.0_cp*a**2*b**4*dn**2+9250432.0_cp*a**2*b**4*dn-9437952.0_cp &
      &              *a**2*b**4+3584.0_cp*a*b**5*dm**2*dn**3-21504.0_cp*a*b**5*dm**2*dn**2+ &
      &              39424.0_cp*a*b**5*dm**2*dn-21504.0_cp*a*b**5*dm**2-16512.0_cp*a*b** &
      &              5*dn**5-165024.0_cp*a*b**5*dn**4+346112.0_cp*a*b**5*dn**3+3535008.0_cp &
      &              *a*b**5*dn**2-10040576.0_cp*a*b**5*dn+6340992.0_cp*a*b**5-256.0_cp* &
      &              b**6*dm**2*dn**3+1536.0_cp*b**6*dm**2*dn**2-2816.0_cp*b**6*dm**2*dn &
      &              +1536.0_cp*b**6*dm**2+3968.0_cp*b**6*dn**5+40032.0_cp*b**6*dn**4-83936.0_cp &
      &              *b**6*dn**3-854304.0_cp*b**6*dn**2+2426976.0_cp*b**6*dn-1532736.0_cp &
      &              *b**6)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-6)=a**7*b*(105.0_cp*a**6*dm**2*dn**3-78.0_cp*a**6*dm**2*dn**2 &
      &              -735.0_cp*a**6*dm**2*dn+372.0_cp*a**6*dm**2+14.0_cp*a**6*dn**5+94.0_cp &
      &              *a**6*dn**4-850.0_cp*a**6*dn**3-6124.0_cp*a**6*dn**2+2936.0_cp*a**6 &
      &              *dn+25242.0_cp*a**6+896.0_cp*a**5*b*dm**2*dn**3-1536.0_cp*a**5*b*dm**2 &
      &              *dn**2-5504.0_cp*a**5*b*dm**2*dn+6144.0_cp*a**5*b*dm**2+16.0_cp*a** &
      &              5*b*dn**5-878.0_cp*a**5*b*dn**4-15422.0_cp*a**5*b*dn**3-39910.0_cp* &
      &              a**5*b*dn**2+180382.0_cp*a**5*b*dn+301668.0_cp*a**5*b+2040.0_cp*a** &
      &              4*b**2*dm**2*dn**3-8400.0_cp*a**4*b**2*dm**2*dn**2-2040.0_cp*a**4*b**2 &
      &              *dm**2*dn+26640.0_cp*a**4*b**2*dm**2-804.0_cp*a**4*b**2*dn**5-12345.0_cp &
      &              *a**4*b**2*dn**4-44065.0_cp*a**4*b**2*dn**3+280971.0_cp*a**4*b**2*dn**2 &
      &              +1231309.0_cp*a**4*b**2*dn-3837594.0_cp*a**4*b**2-1344.0_cp*a**3*b**3* &
      &              dm**2*dn**3-6528.0_cp*a**3*b**3*dm**2*dn**2+58176.0_cp*a**3*b**3*dm**2 &
      &              *dn-79488.0_cp*a**3*b**3*dm**2-1536.0_cp*a**3*b**3*dn**5+15824.0_cp &
      &              *a**3*b**3*dn**4+354448.0_cp*a**3*b**3*dn**3-78512.0_cp*a**3*b**3*dn**2 &
      &              -7046416.0_cp*a**3*b**3*dn+11367264.0_cp*a**3*b**3-4352.0_cp*a**2*b**4 &
      &              *dm**2*dn**3+29184.0_cp*a**2*b**4*dm**2*dn**2-63232.0_cp*a**2*b**4*dm**2 &
      &              *dn+44544.0_cp*a**2*b**4*dm**2+17568.0_cp*a**2*b**4*dn**5+128796.0_cp &
      &              *a**2*b**4*dn**4-566132.0_cp*a**2*b**4*dn**3-2490948.0_cp*a**2*b**4 &
      &              *dn**2+11132660.0_cp*a**2*b**4*dn-10395384.0_cp*a**2*b**4+1024.0_cp &
      &              *a*b**5*dm**2*dn**3-6144.0_cp*a*b**5*dm**2*dn**2+11264.0_cp*a*b**5*dm**2 &
      &              *dn-6144.0_cp*a*b**5*dm**2-13056.0_cp*a*b**5*dn**5-107040.0_cp*a*b**5* &
      &              dn**4+311520.0_cp*a*b**5*dn**3+1981920.0_cp*a*b**5*dn**2-6116064.0_cp &
      &              *a*b**5*dn+3942720.0_cp*a*b**5+1536.0_cp*b**6*dn**5+12736.0_cp*b**6 &
      &              *dn**4-37056.0_cp*b**6*dn**3-234304.0_cp*b**6*dn**2+723648.0_cp*b** &
      &              6*dn-466560.0_cp*b**6)/(2048.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-5)=-a**6*(121.0_cp*a**8*dm**2*dn**3+66.0_cp*a**8*dm**2*dn**2 &
      &              -1069.0_cp*a**8*dm**2*dn+186.0_cp*a**8*dm**2+39.0_cp*a**8*dn**5+343.0_cp &
      &              *a**8*dn**4+9.0_cp*a**8*dn**3-6427.0_cp*a**8*dn**2-6960.0_cp*a**8*dn &
      &              +16500.0_cp*a**8+1280.0_cp*a**7*b*dm**2*dn**3-384.0_cp*a**7*b*dm**2 &
      &              *dn**2-10880.0_cp*a**7*b*dm**2*dn+5376.0_cp*a**7*b*dm**2+384.0_cp*a**7 &
      &              *b*dn**5+2832.0_cp*a**7*b*dn**4-7008.0_cp*a**7*b*dn**3-77424.0_cp*a**7 &
      &              *b*dn**2-1152.0_cp*a**7*b*dn+279744.0_cp*a**7*b+2704.0_cp*a**6*b**2 &
      &              *dm**2*dn**3-9264.0_cp*a**6*b**2*dm**2*dn**2-18976.0_cp*a**6*b**2*dm**2 &
      &              *dn+47616.0_cp*a**6*b**2*dm**2+576.0_cp*a**6*b**2*dn**5-1048.0_cp*a**6 &
      &              *b**2*dn**4-62160.0_cp*a**6*b**2*dn**3-80504.0_cp*a**6*b**2*dn**2+1113024.0_cp &
      &              *a**6*b**2*dn+1241184.0_cp*a**6*b**2-12928.0_cp*a**5*b**3*dm**2*dn**3- &
      &              22272.0_cp*a**5*b**3*dm**2*dn**2+164992.0_cp*a**5*b**3*dm**2*dn+54528.0_cp &
      &              *a**5*b**3*dm**2-3840.0_cp*a**5*b**3*dn**5-25280.0_cp*a**5*b**3*dn**4+ &
      &              273408.0_cp*a**5*b**3*dn**3+3034304.0_cp*a**5*b**3*dn**2+1809408.0_cp &
      &              *a**5*b**3*dn-37138176.0_cp*a**5*b**3-56192.0_cp*a**4*b**4*dm**2*dn**3 &
      &              +123648.0_cp*a**4*b**4*dm**2*dn**2+549248.0_cp*a**4*b**4*dm**2*dn-1243392.0_cp &
      &              *a**4*b**4*dm**2+10544.0_cp*a**4*b**4*dn**5+451680.0_cp*a**4*b**4*dn**4 &
      &              +3158416.0_cp*a**4*b**4*dn**3-9692064.0_cp*a**4*b**4*dn**2-62068608.0_cp &
      &              *a**4*b**4*dn+149008896.0_cp*a**4*b**4-16384.0_cp*a**3*b**5*dm**2*dn**3 &
      &              +270336.0_cp*a**3*b**5*dm**2*dn**2-1040384.0_cp*a**3*b**5*dm**2*dn+1130496.0_cp &
      &              *a**3*b**5*dm**2+126976.0_cp*a**3*b**5*dn**5+176384.0_cp*a**3*b**5*dn**4 &
      &              -8717824.0_cp*a**3*b**5*dn**3-2374400.0_cp*a**3*b**5*dn**2+139253760.0_cp &
      &              *a**3*b**5*dn-206152704.0_cp*a**3*b**5+24576.0_cp*a**2*b**6*dm**2*dn**3 &
      &              -159744.0_cp*a**2*b**6*dm**2*dn**2+331776.0_cp*a**2*b**6*dm**2*dn-221184.0_cp &
      &              *a**2*b**6*dm**2-270336.0_cp*a**2*b**6*dn**5-1556224.0_cp*a**2*b**6 &
      &              *dn**4+8076800.0_cp*a**2*b**6*dn**3+24163072.0_cp*a**2*b**6*dn**2-116079104.0_cp &
      &              *a**2*b**6*dn+104441856.0_cp*a**2*b**6+90112.0_cp*a*b**7*dn**5+577536.0_cp &
      &              *a*b**7*dn**4-2267136.0_cp*a*b**7*dn**3-8945664.0_cp*a*b**7*dn**2+31250432.0_cp &
      &              *a*b**7*dn-20705280.0_cp*a*b**7-4096.0_cp*b**8*dn**5-26624.0_cp*b** &
      &              8*dn**4+104448.0_cp*b**8*dn**3+407552.0_cp*b**8*dn**2-1427456.0_cp* &
      &              b**8*dn+946176.0_cp*b**8)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn &
      &              -1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-4)=-a**6*b*(175.0_cp*a**7*dm**2*dn**3+366.0_cp*a**7*dm** &
      &              2*dn**2-1375.0_cp*a**7*dm**2*dn-1014.0_cp*a**7*dm**2+70.0_cp*a**7*dn**5 &
      &              +649.0_cp*a**7*dn**4+1135.0_cp*a**7*dn**3-7507.0_cp*a**7*dn**2-25085.0_cp &
      &              *a**7*dn-3486.0_cp*a**7+1920.0_cp*a**6*b*dm**2*dn**3+2304.0_cp*a**6 &
      &              *b*dm**2*dn**2-17280.0_cp*a**6*b*dm**2*dn-5376.0_cp*a**6*b*dm**2+800.0_cp &
      &              *a**6*b*dn**5+6412.0_cp*a**6*b*dn**4-700.0_cp*a**6*b*dn**3-136084.0_cp &
      &              *a**6*b*dn**2-278020.0_cp*a**6*b*dn+215784.0_cp*a**6*b+7280.0_cp*a**5* &
      &              b**2*dm**2*dn**3-2400.0_cp*a**5*b**2*dm**2*dn**2-77360.0_cp*a**5*b**2* &
      &              dm**2*dn+36000.0_cp*a**5*b**2*dm**2+2936.0_cp*a**5*b**2*dn**5+11802.0_cp &
      &              *a**5*b**2*dn**4-133250.0_cp*a**5*b**2*dn**3-843462.0_cp*a**5*b**2*dn**2 &
      &              +820314.0_cp*a**5*b**2*dn+9652476.0_cp*a**5*b**2+6656.0_cp*a**4*b** &
      &              3*dm**2*dn**3-43008.0_cp*a**4*b**3*dm**2*dn**2-57344.0_cp*a**4*b**3 &
      &              *dm**2*dn+379392.0_cp*a**4*b**3*dm**2-1024.0_cp*a**4*b**3*dn**5-89856.0_cp &
      &              *a**4*b**3*dn**4-372480.0_cp*a**4*b**3*dn**3+4176000.0_cp*a**4*b**3 &
      &              *dn**2+9854464.0_cp*a**4*b**3*dn-49563264.0_cp*a**4*b**3-18432.0_cp &
      &              *a**3*b**4*dm**2*dn**3+380928.0_cp*a**3*b**4*dm**2*dn-645120.0_cp*a**3 &
      &              *b**4*dm**2-10496.0_cp*a**3*b**4*dn**5+303008.0_cp*a**3*b**4*dn**4+2617696.0_cp &
      &              *a**3*b**4*dn**3-7289024.0_cp*a**3*b**4*dn**2-39065888.0_cp*a**3*b**4* &
      &              dn+90127968.0_cp*a**3*b**4-6144.0_cp*a**2*b**5*dm**2*dn**3+61440.0_cp &
      &              *a**2*b**5*dm**2*dn**2-190464.0_cp*a**2*b**5*dm**2*dn+184320.0_cp*a**2 &
      &              *b**5*dm**2+82432.0_cp*a**2*b**5*dn**5+120384.0_cp*a**2*b**5*dn**4-3701568.0_cp &
      &              *a**2*b**5*dn**3+145728.0_cp*a**2*b**5*dn**2+44383040.0_cp*a**2*b** &
      &              5*dn-64300416.0_cp*a**2*b**5-66560.0_cp*a*b**6*dn**5-272256.0_cp*a* &
      &              b**6*dn**4+1833088.0_cp*a*b**6*dn**3+3174528.0_cp*a*b**6*dn**2-18960512.0_cp &
      &              *a*b**6*dn+17044224.0_cp*a*b**6+8192.0_cp*b**7*dn**5+37888.0_cp*b** &
      &              7*dn**4-203776.0_cp*b**7*dn**3-461824.0_cp*b**7*dn**2+1989632.0_cp* &
      &              b**7*dn-1370112.0_cp*b**7)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn &
      &              -1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-3)=a**6*(11.0_cp*a**8*dm**2*dn**3+264.0_cp*a**8*dm**2*dn**2 &
      &              -329.0_cp*a**8*dm**2*dn-546.0_cp*a**8*dm**2+19.0_cp*a**8*dn**5+168.0_cp &
      &              *a**8*dn**4+489.0_cp*a**8*dn**3-636.0_cp*a**8*dn**2-5800.0_cp*a**8*dn &
      &              -2904.0_cp*a**8+240.0_cp*a**7*b*dm**2*dn**3+2592.0_cp*a**7*b*dm**2*dn**2 &
      &              -5040.0_cp*a**7*b*dm**2*dn-5088.0_cp*a**7*b*dm**2+272.0_cp*a**7*b*dn**5 &
      &              +2204.0_cp*a**7*b*dn**4+4336.0_cp*a**7*b*dn**3-18428.0_cp*a**7*b*dn**2 &
      &              -80016.0_cp*a**7*b*dn-14832.0_cp*a**7*b+882.0_cp*a**6*b**2*dm**2*dn**3 &
      &              +7452.0_cp*a**6*b**2*dm**2*dn**2-29658.0_cp*a**6*b**2*dm**2*dn-5988.0_cp &
      &              *a**6*b**2*dm**2+1348.0_cp*a**6*b**2*dn**5+8254.0_cp*a**6*b**2*dn** &
      &              4-7140.0_cp*a**6*b**2*dn**3-181438.0_cp*a**6*b**2*dn**2-333688.0_cp &
      &              *a**6*b**2*dn+458568.0_cp*a**6*b**2-3584.0_cp*a**5*b**3*dm**2*dn**3 &
      &              -6144.0_cp*a**5*b**3*dm**2*dn**2-54784.0_cp*a**5*b**3*dm**2*dn+124416.0_cp &
      &              *a**5*b**3*dm**2+1024.0_cp*a**5*b**3*dn**5-17664.0_cp*a**5*b**3*dn**4- &
      &              201216.0_cp*a**5*b**3*dn**3-486912.0_cp*a**5*b**3*dn**2+3599360.0_cp &
      &              *a**5*b**3*dn+15209472.0_cp*a**5*b**3-31232.0_cp*a**4*b**4*dm**2*dn**3 &
      &              -43008.0_cp*a**4*b**4*dm**2*dn**2+224768.0_cp*a**4*b**4*dm**2*dn+655872.0_cp &
      &              *a**4*b**4*dm**2-15296.0_cp*a**4*b**4*dn**5-199680.0_cp*a**4*b**4*dn**4 &
      &              +332096.0_cp*a**4*b**4*dn**3+11346624.0_cp*a**4*b**4*dn**2+2598912.0_cp &
      &              *a**4*b**4*dn-113396736.0_cp*a**4*b**4-61440.0_cp*a**3*b**5*dm**2*dn**3 &
      &              +184320.0_cp*a**3*b**5*dm**2*dn**2+675840.0_cp*a**3*b**5*dm**2*dn-2027520.0_cp &
      &              *a**3*b**5*dm**2-7168.0_cp*a**3*b**5*dn**5+867968.0_cp*a**3*b**5*dn**4 &
      &              +3436032.0_cp*a**3*b**5*dn**3-28279424.0_cp*a**3*b**5*dn**2-50155520.0_cp &
      &              *a**3*b**5*dn+243644928.0_cp*a**3*b**5+18432.0_cp*a**2*b**6*dm**2*dn**3 &
      &              +43008.0_cp*a**2*b**6*dm**2*dn**2-595968.0_cp*a**2*b**6*dm**2*dn+903168.0_cp &
      &              *a**2*b**6*dm**2+78848.0_cp*a**2*b**6*dn**5-1127808.0_cp*a**2*b**6*dn**4 &
      &              -7985920.0_cp*a**2*b**6*dn**3+27421824.0_cp*a**2*b**6*dn**2+88531712.0_cp &
      &              *a**2*b**6*dn-224579328.0_cp*a**2*b**6-188416.0_cp*a*b**7*dn**5-73728.0_cp &
      &              *a*b**7*dn**4+6248448.0_cp*a*b**7*dn**3-4153344.0_cp*a*b**7*dn**2-53405696.0_cp &
      &              *a*b**7*dn+80646144.0_cp*a*b**7+53248.0_cp*b**8*dn**5+126976.0_cp*b**8 &
      &              *dn**4-1296384.0_cp*b**8*dn**3-851968.0_cp*b**8*dn**2+9402368.0_cp* &
      &              b**8*dn-8761344.0_cp*b**8)/(8192.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn &
      &              -1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-2)=a**6*b*(45.0_cp*a**7*dm**2*dn**2+279.0_cp*a**7*dm**2*dn &
      &              -462.0_cp*a**7*dm**2+26.0_cp*a**7*dn**4+225.0_cp*a**7*dn**3+796.0_cp &
      &              *a**7*dn**2+123.0_cp*a**7*dn-4554.0_cp*a**7+768.0_cp*a**6*b*dm**2*dn**2 &
      &              +2304.0_cp*a**6*b*dm**2*dn-6144.0_cp*a**6*b*dm**2+448.0_cp*a**6*b*dn**4 &
      &              +3480.0_cp*a**6*b*dn**3+8624.0_cp*a**6*b*dn**2-9912.0_cp*a**6*b*dn-93648.0_cp &
      &              *a**6*b+4640.0_cp*a**5*b**2*dm**2*dn**2+4320.0_cp*a**5*b**2*dm**2*dn &
      &              -34880.0_cp*a**5*b**2*dm**2+3216.0_cp*a**5*b**2*dn**4+21956.0_cp*a**5* &
      &              b**2*dn**3+14552.0_cp*a**5*b**2*dn**2-211124.0_cp*a**5*b**2*dn-1861304.0_cp &
      &              *a**5*b**2+11776.0_cp*a**4*b**3*dm**2*dn**2-23040.0_cp*a**4*b**3*dm**2 &
      &              *dn-85504.0_cp*a**4*b**3*dm**2+11264.0_cp*a**4*b**3*dn**4+52224.0_cp &
      &              *a**4*b**3*dn**3-623104.0_cp*a**4*b**3*dn**2-3359616.0_cp*a**4*b**3 &
      &              *dn+17658752.0_cp*a**4*b**3+2048.0_cp*a**3*b**4*dm**2*dn**2-129024.0_cp &
      &              *a**3*b**4*dm**2*dn+378880.0_cp*a**3*b**4*dm**2-8960.0_cp*a**3*b**4 &
      &              *dn**4-437600.0_cp*a**3*b**4*dn**3+1527232.0_cp*a**3*b**4*dn**2+13877312.0_cp &
      &              *a**3*b**4*dn-44858144.0_cp*a**3*b**4-30720.0_cp*a**2*b**5*dm**2*dn**2 &
      &              +178176.0_cp*a**2*b**5*dm**2*dn-258048.0_cp*a**2*b**5*dm**2+10752.0_cp &
      &              *a**2*b**5*dn**4+978112.0_cp*a**2*b**5*dn**3-1163904.0_cp*a**2*b**5 &
      &              *dn**2-21589696.0_cp*a**2*b**5*dn+47964288.0_cp*a**2*b**5+44032.0_cp &
      &              *a*b**6*dn**4-752256.0_cp*a*b**6*dn**3-269824.0_cp*a*b**6*dn**2+13739904.0_cp &
      &              *a*b**6*dn-22046976.0_cp*a*b**6-40960.0_cp*b**7*dn**4+164864.0_cp*b**7 &
      &              *dn**3+483328.0_cp*b**7*dn**2-2886656.0_cp*b**7*dn+3176448.0_cp*b** &
      &              7)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn &
      &              +2.0_cp))
      stencil(ku-1)=a**6*(165.0_cp*a**8*dm**2*dn**2-924.0_cp*a**8*dm**2*dn &
      &              +783.0_cp*a**8*dm**2+27.0_cp*a**8*dn**4-117.0_cp*a**8*dn**3-497.0_cp &
      &              *a**8*dn**2-1595.0_cp*a**8*dn+3462.0_cp*a**8+1536.0_cp*a**7*b*dm**2 &
      &              *dn**2-9984.0_cp*a**7*b*dm**2*dn+9984.0_cp*a**7*b*dm**2+256.0_cp*a**7* &
      &              b*dn**4-2016.0_cp*a**7*b*dn**3-6784.0_cp*a**7*b*dn**2-14304.0_cp*a**7* &
      &              b*dn+52416.0_cp*a**7*b+6240.0_cp*a**6*b**2*dm**2*dn**2-45600.0_cp*a**6 &
      &              *b**2*dm**2*dn+59520.0_cp*a**6*b**2*dm**2+896.0_cp*a**6*b**2*dn**4-15024.0_cp &
      &              *a**6*b**2*dn**3-32128.0_cp*a**6*b**2*dn**2-17264.0_cp*a**6*b**2*dn &
      &              +480096.0_cp*a**6*b**2+16640.0_cp*a**5*b**3*dm**2*dn**2-102400.0_cp &
      &              *a**5*b**3*dm**2*dn+203520.0_cp*a**5*b**3*dm**2+1536.0_cp*a**5*b**3 &
      &              *dn**4-57216.0_cp*a**5*b**3*dn**3+37632.0_cp*a**5*b**3*dn**2+499584.0_cp &
      &              *a**5*b**3*dn+7055616.0_cp*a**5*b**3+43776.0_cp*a**4*b**4*dm**2*dn**2- &
      &              61440.0_cp*a**4*b**4*dm**2*dn+297216.0_cp*a**4*b**4*dm**2+9504.0_cp &
      &              *a**4*b**4*dn**4-43904.0_cp*a**4*b**4*dn**3+2611296.0_cp*a**4*b**4*dn**2 &
      &              +8441216.0_cp*a**4*b**4*dn-81128448.0_cp*a**4*b**4+114688.0_cp*a**3 &
      &              *b**5*dm**2*dn**2+188416.0_cp*a**3*b**5*dm**2*dn-1941504.0_cp*a**3* &
      &              b**5*dm**2+159744.0_cp*a**3*b**5*dn**4+1475328.0_cp*a**3*b**5*dn**3 &
      &              -11774976.0_cp*a**3*b**5*dn**2-46291200.0_cp*a**3*b**5*dn+231737856.0_cp &
      &              *a**3*b**5+90112.0_cp*a**2*b**6*dm**2*dn**2-856064.0_cp*a**2*b**6*dm**2 &
      &              *dn+1781760.0_cp*a**2*b**6*dm**2-172032.0_cp*a**2*b**6*dn**4-3599616.0_cp &
      &              *a**2*b**6*dn**3+16242688.0_cp*a**2*b**6*dn**2+76892416.0_cp*a**2*b**6 &
      &              *dn-275303424.0_cp*a**2*b**6+221184.0_cp*a*b**7*dn**4+3682304.0_cp* &
      &              a*b**7*dn**3-9992192.0_cp*a*b**7*dn**2-58003456.0_cp*a*b**7*dn+146601984.0_cp &
      &              *a*b**7+61440.0_cp*b**8*dn**4-1234944.0_cp*b**8*dn**3+1067008.0_cp* &
      &              b**8*dn**2+15210496.0_cp*b**8*dn-26867712.0_cp*b**8)/(16384.0_cp*dn &
      &              *(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+3.0_cp))
      stencil(ku)=a**6*b*(3.0_cp*a**7*dm**2*dn**2+15.0_cp*a**7*dm**2*dn-192.0_cp &
      &              *a**7*dm**2+2.0_cp*a**7*dn**4+18.0_cp*a**7*dn**3-78.0_cp*a**7*dn**2 &
      &              -212.0_cp*a**7*dn-750.0_cp*a**7-1920.0_cp*a**6*b*dm**2+16.0_cp*a**6 &
      &              *b*dn**4+174.0_cp*a**6*b*dn**3-1452.0_cp*a**6*b*dn**2-3406.0_cp*a** &
      &              6*b*dn-14868.0_cp*a**6*b-200.0_cp*a**5*b**2*dm**2*dn**2-1000.0_cp*a**5 &
      &              *b**2*dm**2*dn-7440.0_cp*a**5*b**2*dm**2-4.0_cp*a**5*b**2*dn**4+445.0_cp &
      &              *a**5*b**2*dn**3-14726.0_cp*a**5*b**2*dn**2-28525.0_cp*a**5*b**2*dn &
      &              -370782.0_cp*a**5*b**2-1216.0_cp*a**4*b**3*dm**2*dn**2-6080.0_cp*a**4* &
      &              b**3*dm**2*dn-5376.0_cp*a**4*b**3*dm**2-768.0_cp*a**4*b**3*dn**4-2832.0_cp &
      &              *a**4*b**3*dn**3-176160.0_cp*a**4*b**3*dn**2-279312.0_cp*a**4*b**3*dn &
      &              +5547648.0_cp*a**4*b**3-2816.0_cp*a**3*b**4*dm**2*dn**2-14080.0_cp* &
      &              a**3*b**4*dm**2*dn+144384.0_cp*a**3*b**4*dm**2-9888.0_cp*a**3*b**4*dn**4 &
      &              -52980.0_cp*a**3*b**4*dn**3+948792.0_cp*a**3*b**4*dn**2+1588860.0_cp &
      &              *a**3*b**4*dn-16669680.0_cp*a**3*b**4+5632.0_cp*a**2*b**5*dm**2*dn**2+ &
      &              28160.0_cp*a**2*b**5*dm**2*dn-150528.0_cp*a**2*b**5*dm**2+24448.0_cp &
      &              *a**2*b**5*dn**4+136112.0_cp*a**2*b**5*dn**3-1659872.0_cp*a**2*b**5 &
      &              *dn**2-2897008.0_cp*a**2*b**5*dn+21252000.0_cp*a**2*b**5-28416.0_cp &
      &              *a*b**6*dn**4-153504.0_cp*a*b**6*dn**3+1261952.0_cp*a*b**6*dn**2+2298976.0_cp &
      &              *a*b**6*dn-12194112.0_cp*a*b**6+12288.0_cp*b**7*dn**4+69120.0_cp*b**7* &
      &              dn**3-364544.0_cp*b**7*dn**2-724480.0_cp*b**7*dn+2592768.0_cp*b**7) &
      &              /(1024.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku+1)=-a**6*(33.0_cp*a**8*dm**2*dn**2-267.0_cp*a**8*dm**2+9.0_cp &
      &              *a**8*dn**4-29.0_cp*a**8*dn**2-568.0_cp*a**8+336.0_cp*a**7*b*dm**2*dn**2 &
      &              -3024.0_cp*a**7*b*dm**2+112.0_cp*a**7*b*dn**4-608.0_cp*a**7*b*dn**2 &
      &              -7280.0_cp*a**7*b+1470.0_cp*a**6*b**2*dm**2*dn**2-15390.0_cp*a**6*b**2 &
      &              *dm**2+620.0_cp*a**6*b**2*dn**4-5812.0_cp*a**6*b**2*dn**2-54952.0_cp &
      &              *a**6*b**2+3520.0_cp*a**5*b**3*dm**2*dn**2-43840.0_cp*a**5*b**3*dm**2+ &
      &              1920.0_cp*a**5*b**3*dn**4-38976.0_cp*a**5*b**3*dn**2-794496.0_cp*a**5* &
      &              b**3+4928.0_cp*a**4*b**4*dm**2*dn**2-55232.0_cp*a**4*b**4*dm**2+2808.0_cp &
      &              *a**4*b**4*dn**4-361320.0_cp*a**4*b**4*dn**2+11221632.0_cp*a**4*b** &
      &              4+5632.0_cp*a**3*b**5*dm**2*dn**2+246272.0_cp*a**3*b**5*dm**2-11904.0_cp &
      &              *a**3*b**5*dn**4+2037504.0_cp*a**3*b**5*dn**2-35142912.0_cp*a**3*b**5+ &
      &              25344.0_cp*a**2*b**6*dm**2*dn**2-366336.0_cp*a**2*b**6*dm**2+67968.0_cp &
      &              *a**2*b**6*dn**4-3835232.0_cp*a**2*b**6*dn**2+45332960.0_cp*a**2*b**6- &
      &              61440.0_cp*a*b**7*dn**4+3034112.0_cp*a*b**7*dn**2-26866688.0_cp*a*b**7 &
      &              +43008.0_cp*b**8*dn**4-1025024.0_cp*b**8*dn**2+5946368.0_cp*b**8)/(2048.0_cp &
      &              *(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku+2)=a**6*b*(3.0_cp*a**7*dm**2*dn**2-15.0_cp*a**7*dm**2*dn &
      &              -192.0_cp*a**7*dm**2+2.0_cp*a**7*dn**4-18.0_cp*a**7*dn**3-78.0_cp*a**7 &
      &              *dn**2+212.0_cp*a**7*dn-750.0_cp*a**7-1920.0_cp*a**6*b*dm**2+16.0_cp &
      &              *a**6*b*dn**4-174.0_cp*a**6*b*dn**3-1452.0_cp*a**6*b*dn**2+3406.0_cp &
      &              *a**6*b*dn-14868.0_cp*a**6*b-200.0_cp*a**5*b**2*dm**2*dn**2+1000.0_cp &
      &              *a**5*b**2*dm**2*dn-7440.0_cp*a**5*b**2*dm**2-4.0_cp*a**5*b**2*dn** &
      &              4-445.0_cp*a**5*b**2*dn**3-14726.0_cp*a**5*b**2*dn**2+28525.0_cp*a**5* &
      &              b**2*dn-370782.0_cp*a**5*b**2-1216.0_cp*a**4*b**3*dm**2*dn**2+6080.0_cp &
      &              *a**4*b**3*dm**2*dn-5376.0_cp*a**4*b**3*dm**2-768.0_cp*a**4*b**3*dn**4 &
      &              +2832.0_cp*a**4*b**3*dn**3-176160.0_cp*a**4*b**3*dn**2+279312.0_cp* &
      &              a**4*b**3*dn+5547648.0_cp*a**4*b**3-2816.0_cp*a**3*b**4*dm**2*dn**2 &
      &              +14080.0_cp*a**3*b**4*dm**2*dn+144384.0_cp*a**3*b**4*dm**2-9888.0_cp &
      &              *a**3*b**4*dn**4+52980.0_cp*a**3*b**4*dn**3+948792.0_cp*a**3*b**4*dn**2 &
      &              -1588860.0_cp*a**3*b**4*dn-16669680.0_cp*a**3*b**4+5632.0_cp*a**2*b**5 &
      &              *dm**2*dn**2-28160.0_cp*a**2*b**5*dm**2*dn-150528.0_cp*a**2*b**5*dm**2 &
      &              +24448.0_cp*a**2*b**5*dn**4-136112.0_cp*a**2*b**5*dn**3-1659872.0_cp &
      &              *a**2*b**5*dn**2+2897008.0_cp*a**2*b**5*dn+21252000.0_cp*a**2*b**5-28416.0_cp &
      &              *a*b**6*dn**4+153504.0_cp*a*b**6*dn**3+1261952.0_cp*a*b**6*dn**2-2298976.0_cp &
      &              *a*b**6*dn-12194112.0_cp*a*b**6+12288.0_cp*b**7*dn**4-69120.0_cp*b**7* &
      &              dn**3-364544.0_cp*b**7*dn**2+724480.0_cp*b**7*dn+2592768.0_cp*b**7) &
      &              /(1024.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku+3)=a**6*(165.0_cp*a**8*dm**2*dn**2+924.0_cp*a**8*dm**2*dn &
      &              +783.0_cp*a**8*dm**2+27.0_cp*a**8*dn**4+117.0_cp*a**8*dn**3-497.0_cp &
      &              *a**8*dn**2+1595.0_cp*a**8*dn+3462.0_cp*a**8+1536.0_cp*a**7*b*dm**2 &
      &              *dn**2+9984.0_cp*a**7*b*dm**2*dn+9984.0_cp*a**7*b*dm**2+256.0_cp*a**7* &
      &              b*dn**4+2016.0_cp*a**7*b*dn**3-6784.0_cp*a**7*b*dn**2+14304.0_cp*a**7* &
      &              b*dn+52416.0_cp*a**7*b+6240.0_cp*a**6*b**2*dm**2*dn**2+45600.0_cp*a**6 &
      &              *b**2*dm**2*dn+59520.0_cp*a**6*b**2*dm**2+896.0_cp*a**6*b**2*dn**4+15024.0_cp &
      &              *a**6*b**2*dn**3-32128.0_cp*a**6*b**2*dn**2+17264.0_cp*a**6*b**2*dn &
      &              +480096.0_cp*a**6*b**2+16640.0_cp*a**5*b**3*dm**2*dn**2+102400.0_cp &
      &              *a**5*b**3*dm**2*dn+203520.0_cp*a**5*b**3*dm**2+1536.0_cp*a**5*b**3 &
      &              *dn**4+57216.0_cp*a**5*b**3*dn**3+37632.0_cp*a**5*b**3*dn**2-499584.0_cp &
      &              *a**5*b**3*dn+7055616.0_cp*a**5*b**3+43776.0_cp*a**4*b**4*dm**2*dn**2+ &
      &              61440.0_cp*a**4*b**4*dm**2*dn+297216.0_cp*a**4*b**4*dm**2+9504.0_cp &
      &              *a**4*b**4*dn**4+43904.0_cp*a**4*b**4*dn**3+2611296.0_cp*a**4*b**4*dn**2 &
      &              -8441216.0_cp*a**4*b**4*dn-81128448.0_cp*a**4*b**4+114688.0_cp*a**3 &
      &              *b**5*dm**2*dn**2-188416.0_cp*a**3*b**5*dm**2*dn-1941504.0_cp*a**3* &
      &              b**5*dm**2+159744.0_cp*a**3*b**5*dn**4-1475328.0_cp*a**3*b**5*dn**3 &
      &              -11774976.0_cp*a**3*b**5*dn**2+46291200.0_cp*a**3*b**5*dn+231737856.0_cp &
      &              *a**3*b**5+90112.0_cp*a**2*b**6*dm**2*dn**2+856064.0_cp*a**2*b**6*dm**2 &
      &              *dn+1781760.0_cp*a**2*b**6*dm**2-172032.0_cp*a**2*b**6*dn**4+3599616.0_cp &
      &              *a**2*b**6*dn**3+16242688.0_cp*a**2*b**6*dn**2-76892416.0_cp*a**2*b**6 &
      &              *dn-275303424.0_cp*a**2*b**6+221184.0_cp*a*b**7*dn**4-3682304.0_cp* &
      &              a*b**7*dn**3-9992192.0_cp*a*b**7*dn**2+58003456.0_cp*a*b**7*dn+146601984.0_cp &
      &              *a*b**7+61440.0_cp*b**8*dn**4+1234944.0_cp*b**8*dn**3+1067008.0_cp* &
      &              b**8*dn**2-15210496.0_cp*b**8*dn-26867712.0_cp*b**8)/(16384.0_cp*dn &
      &              *(dn-3.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+4)=a**6*b*(45.0_cp*a**7*dm**2*dn**2-279.0_cp*a**7*dm**2*dn &
      &              -462.0_cp*a**7*dm**2+26.0_cp*a**7*dn**4-225.0_cp*a**7*dn**3+796.0_cp &
      &              *a**7*dn**2-123.0_cp*a**7*dn-4554.0_cp*a**7+768.0_cp*a**6*b*dm**2*dn**2 &
      &              -2304.0_cp*a**6*b*dm**2*dn-6144.0_cp*a**6*b*dm**2+448.0_cp*a**6*b*dn**4 &
      &              -3480.0_cp*a**6*b*dn**3+8624.0_cp*a**6*b*dn**2+9912.0_cp*a**6*b*dn-93648.0_cp &
      &              *a**6*b+4640.0_cp*a**5*b**2*dm**2*dn**2-4320.0_cp*a**5*b**2*dm**2*dn &
      &              -34880.0_cp*a**5*b**2*dm**2+3216.0_cp*a**5*b**2*dn**4-21956.0_cp*a**5* &
      &              b**2*dn**3+14552.0_cp*a**5*b**2*dn**2+211124.0_cp*a**5*b**2*dn-1861304.0_cp &
      &              *a**5*b**2+11776.0_cp*a**4*b**3*dm**2*dn**2+23040.0_cp*a**4*b**3*dm**2 &
      &              *dn-85504.0_cp*a**4*b**3*dm**2+11264.0_cp*a**4*b**3*dn**4-52224.0_cp &
      &              *a**4*b**3*dn**3-623104.0_cp*a**4*b**3*dn**2+3359616.0_cp*a**4*b**3 &
      &              *dn+17658752.0_cp*a**4*b**3+2048.0_cp*a**3*b**4*dm**2*dn**2+129024.0_cp &
      &              *a**3*b**4*dm**2*dn+378880.0_cp*a**3*b**4*dm**2-8960.0_cp*a**3*b**4 &
      &              *dn**4+437600.0_cp*a**3*b**4*dn**3+1527232.0_cp*a**3*b**4*dn**2-13877312.0_cp &
      &              *a**3*b**4*dn-44858144.0_cp*a**3*b**4-30720.0_cp*a**2*b**5*dm**2*dn**2 &
      &              -178176.0_cp*a**2*b**5*dm**2*dn-258048.0_cp*a**2*b**5*dm**2+10752.0_cp &
      &              *a**2*b**5*dn**4-978112.0_cp*a**2*b**5*dn**3-1163904.0_cp*a**2*b**5 &
      &              *dn**2+21589696.0_cp*a**2*b**5*dn+47964288.0_cp*a**2*b**5+44032.0_cp &
      &              *a*b**6*dn**4+752256.0_cp*a*b**6*dn**3-269824.0_cp*a*b**6*dn**2-13739904.0_cp &
      &              *a*b**6*dn-22046976.0_cp*a*b**6-40960.0_cp*b**7*dn**4-164864.0_cp*b**7 &
      &              *dn**3+483328.0_cp*b**7*dn**2+2886656.0_cp*b**7*dn+3176448.0_cp*b** &
      &              7)/(4096.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn &
      &              +3.0_cp))
      stencil(ku+5)=a**6*(11.0_cp*a**8*dm**2*dn**3-264.0_cp*a**8*dm**2*dn**2 &
      &              -329.0_cp*a**8*dm**2*dn+546.0_cp*a**8*dm**2+19.0_cp*a**8*dn**5-168.0_cp &
      &              *a**8*dn**4+489.0_cp*a**8*dn**3+636.0_cp*a**8*dn**2-5800.0_cp*a**8*dn &
      &              +2904.0_cp*a**8+240.0_cp*a**7*b*dm**2*dn**3-2592.0_cp*a**7*b*dm**2*dn**2 &
      &              -5040.0_cp*a**7*b*dm**2*dn+5088.0_cp*a**7*b*dm**2+272.0_cp*a**7*b*dn**5 &
      &              -2204.0_cp*a**7*b*dn**4+4336.0_cp*a**7*b*dn**3+18428.0_cp*a**7*b*dn**2 &
      &              -80016.0_cp*a**7*b*dn+14832.0_cp*a**7*b+882.0_cp*a**6*b**2*dm**2*dn**3 &
      &              -7452.0_cp*a**6*b**2*dm**2*dn**2-29658.0_cp*a**6*b**2*dm**2*dn+5988.0_cp &
      &              *a**6*b**2*dm**2+1348.0_cp*a**6*b**2*dn**5-8254.0_cp*a**6*b**2*dn** &
      &              4-7140.0_cp*a**6*b**2*dn**3+181438.0_cp*a**6*b**2*dn**2-333688.0_cp &
      &              *a**6*b**2*dn-458568.0_cp*a**6*b**2-3584.0_cp*a**5*b**3*dm**2*dn**3 &
      &              +6144.0_cp*a**5*b**3*dm**2*dn**2-54784.0_cp*a**5*b**3*dm**2*dn-124416.0_cp &
      &              *a**5*b**3*dm**2+1024.0_cp*a**5*b**3*dn**5+17664.0_cp*a**5*b**3*dn**4- &
      &              201216.0_cp*a**5*b**3*dn**3+486912.0_cp*a**5*b**3*dn**2+3599360.0_cp &
      &              *a**5*b**3*dn-15209472.0_cp*a**5*b**3-31232.0_cp*a**4*b**4*dm**2*dn**3 &
      &              +43008.0_cp*a**4*b**4*dm**2*dn**2+224768.0_cp*a**4*b**4*dm**2*dn-655872.0_cp &
      &              *a**4*b**4*dm**2-15296.0_cp*a**4*b**4*dn**5+199680.0_cp*a**4*b**4*dn**4 &
      &              +332096.0_cp*a**4*b**4*dn**3-11346624.0_cp*a**4*b**4*dn**2+2598912.0_cp &
      &              *a**4*b**4*dn+113396736.0_cp*a**4*b**4-61440.0_cp*a**3*b**5*dm**2*dn**3 &
      &              -184320.0_cp*a**3*b**5*dm**2*dn**2+675840.0_cp*a**3*b**5*dm**2*dn+2027520.0_cp &
      &              *a**3*b**5*dm**2-7168.0_cp*a**3*b**5*dn**5-867968.0_cp*a**3*b**5*dn**4 &
      &              +3436032.0_cp*a**3*b**5*dn**3+28279424.0_cp*a**3*b**5*dn**2-50155520.0_cp &
      &              *a**3*b**5*dn-243644928.0_cp*a**3*b**5+18432.0_cp*a**2*b**6*dm**2*dn**3 &
      &              -43008.0_cp*a**2*b**6*dm**2*dn**2-595968.0_cp*a**2*b**6*dm**2*dn-903168.0_cp &
      &              *a**2*b**6*dm**2+78848.0_cp*a**2*b**6*dn**5+1127808.0_cp*a**2*b**6*dn**4 &
      &              -7985920.0_cp*a**2*b**6*dn**3-27421824.0_cp*a**2*b**6*dn**2+88531712.0_cp &
      &              *a**2*b**6*dn+224579328.0_cp*a**2*b**6-188416.0_cp*a*b**7*dn**5+73728.0_cp &
      &              *a*b**7*dn**4+6248448.0_cp*a*b**7*dn**3+4153344.0_cp*a*b**7*dn**2-53405696.0_cp &
      &              *a*b**7*dn-80646144.0_cp*a*b**7+53248.0_cp*b**8*dn**5-126976.0_cp*b**8 &
      &              *dn**4-1296384.0_cp*b**8*dn**3+851968.0_cp*b**8*dn**2+9402368.0_cp* &
      &              b**8*dn+8761344.0_cp*b**8)/(8192.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn &
      &              -1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+6)=-a**6*b*(175.0_cp*a**7*dm**2*dn**3-366.0_cp*a**7*dm** &
      &              2*dn**2-1375.0_cp*a**7*dm**2*dn+1014.0_cp*a**7*dm**2+70.0_cp*a**7*dn**5 &
      &              -649.0_cp*a**7*dn**4+1135.0_cp*a**7*dn**3+7507.0_cp*a**7*dn**2-25085.0_cp &
      &              *a**7*dn+3486.0_cp*a**7+1920.0_cp*a**6*b*dm**2*dn**3-2304.0_cp*a**6 &
      &              *b*dm**2*dn**2-17280.0_cp*a**6*b*dm**2*dn+5376.0_cp*a**6*b*dm**2+800.0_cp &
      &              *a**6*b*dn**5-6412.0_cp*a**6*b*dn**4-700.0_cp*a**6*b*dn**3+136084.0_cp &
      &              *a**6*b*dn**2-278020.0_cp*a**6*b*dn-215784.0_cp*a**6*b+7280.0_cp*a**5* &
      &              b**2*dm**2*dn**3+2400.0_cp*a**5*b**2*dm**2*dn**2-77360.0_cp*a**5*b**2* &
      &              dm**2*dn-36000.0_cp*a**5*b**2*dm**2+2936.0_cp*a**5*b**2*dn**5-11802.0_cp &
      &              *a**5*b**2*dn**4-133250.0_cp*a**5*b**2*dn**3+843462.0_cp*a**5*b**2*dn**2 &
      &              +820314.0_cp*a**5*b**2*dn-9652476.0_cp*a**5*b**2+6656.0_cp*a**4*b** &
      &              3*dm**2*dn**3+43008.0_cp*a**4*b**3*dm**2*dn**2-57344.0_cp*a**4*b**3 &
      &              *dm**2*dn-379392.0_cp*a**4*b**3*dm**2-1024.0_cp*a**4*b**3*dn**5+89856.0_cp &
      &              *a**4*b**3*dn**4-372480.0_cp*a**4*b**3*dn**3-4176000.0_cp*a**4*b**3 &
      &              *dn**2+9854464.0_cp*a**4*b**3*dn+49563264.0_cp*a**4*b**3-18432.0_cp &
      &              *a**3*b**4*dm**2*dn**3+380928.0_cp*a**3*b**4*dm**2*dn+645120.0_cp*a**3 &
      &              *b**4*dm**2-10496.0_cp*a**3*b**4*dn**5-303008.0_cp*a**3*b**4*dn**4+2617696.0_cp &
      &              *a**3*b**4*dn**3+7289024.0_cp*a**3*b**4*dn**2-39065888.0_cp*a**3*b**4* &
      &              dn-90127968.0_cp*a**3*b**4-6144.0_cp*a**2*b**5*dm**2*dn**3-61440.0_cp &
      &              *a**2*b**5*dm**2*dn**2-190464.0_cp*a**2*b**5*dm**2*dn-184320.0_cp*a**2 &
      &              *b**5*dm**2+82432.0_cp*a**2*b**5*dn**5-120384.0_cp*a**2*b**5*dn**4-3701568.0_cp &
      &              *a**2*b**5*dn**3-145728.0_cp*a**2*b**5*dn**2+44383040.0_cp*a**2*b** &
      &              5*dn+64300416.0_cp*a**2*b**5-66560.0_cp*a*b**6*dn**5+272256.0_cp*a* &
      &              b**6*dn**4+1833088.0_cp*a*b**6*dn**3-3174528.0_cp*a*b**6*dn**2-18960512.0_cp &
      &              *a*b**6*dn-17044224.0_cp*a*b**6+8192.0_cp*b**7*dn**5-37888.0_cp*b** &
      &              7*dn**4-203776.0_cp*b**7*dn**3+461824.0_cp*b**7*dn**2+1989632.0_cp* &
      &              b**7*dn+1370112.0_cp*b**7)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn &
      &              -1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+7)=-a**6*(121.0_cp*a**8*dm**2*dn**3-66.0_cp*a**8*dm**2*dn**2 &
      &              -1069.0_cp*a**8*dm**2*dn-186.0_cp*a**8*dm**2+39.0_cp*a**8*dn**5-343.0_cp &
      &              *a**8*dn**4+9.0_cp*a**8*dn**3+6427.0_cp*a**8*dn**2-6960.0_cp*a**8*dn &
      &              -16500.0_cp*a**8+1280.0_cp*a**7*b*dm**2*dn**3+384.0_cp*a**7*b*dm**2 &
      &              *dn**2-10880.0_cp*a**7*b*dm**2*dn-5376.0_cp*a**7*b*dm**2+384.0_cp*a**7 &
      &              *b*dn**5-2832.0_cp*a**7*b*dn**4-7008.0_cp*a**7*b*dn**3+77424.0_cp*a**7 &
      &              *b*dn**2-1152.0_cp*a**7*b*dn-279744.0_cp*a**7*b+2704.0_cp*a**6*b**2 &
      &              *dm**2*dn**3+9264.0_cp*a**6*b**2*dm**2*dn**2-18976.0_cp*a**6*b**2*dm**2 &
      &              *dn-47616.0_cp*a**6*b**2*dm**2+576.0_cp*a**6*b**2*dn**5+1048.0_cp*a**6 &
      &              *b**2*dn**4-62160.0_cp*a**6*b**2*dn**3+80504.0_cp*a**6*b**2*dn**2+1113024.0_cp &
      &              *a**6*b**2*dn-1241184.0_cp*a**6*b**2-12928.0_cp*a**5*b**3*dm**2*dn**3+ &
      &              22272.0_cp*a**5*b**3*dm**2*dn**2+164992.0_cp*a**5*b**3*dm**2*dn-54528.0_cp &
      &              *a**5*b**3*dm**2-3840.0_cp*a**5*b**3*dn**5+25280.0_cp*a**5*b**3*dn**4+ &
      &              273408.0_cp*a**5*b**3*dn**3-3034304.0_cp*a**5*b**3*dn**2+1809408.0_cp &
      &              *a**5*b**3*dn+37138176.0_cp*a**5*b**3-56192.0_cp*a**4*b**4*dm**2*dn**3 &
      &              -123648.0_cp*a**4*b**4*dm**2*dn**2+549248.0_cp*a**4*b**4*dm**2*dn+1243392.0_cp &
      &              *a**4*b**4*dm**2+10544.0_cp*a**4*b**4*dn**5-451680.0_cp*a**4*b**4*dn**4 &
      &              +3158416.0_cp*a**4*b**4*dn**3+9692064.0_cp*a**4*b**4*dn**2-62068608.0_cp &
      &              *a**4*b**4*dn-149008896.0_cp*a**4*b**4-16384.0_cp*a**3*b**5*dm**2*dn**3 &
      &              -270336.0_cp*a**3*b**5*dm**2*dn**2-1040384.0_cp*a**3*b**5*dm**2*dn-1130496.0_cp &
      &              *a**3*b**5*dm**2+126976.0_cp*a**3*b**5*dn**5-176384.0_cp*a**3*b**5*dn**4 &
      &              -8717824.0_cp*a**3*b**5*dn**3+2374400.0_cp*a**3*b**5*dn**2+139253760.0_cp &
      &              *a**3*b**5*dn+206152704.0_cp*a**3*b**5+24576.0_cp*a**2*b**6*dm**2*dn**3 &
      &              +159744.0_cp*a**2*b**6*dm**2*dn**2+331776.0_cp*a**2*b**6*dm**2*dn+221184.0_cp &
      &              *a**2*b**6*dm**2-270336.0_cp*a**2*b**6*dn**5+1556224.0_cp*a**2*b**6 &
      &              *dn**4+8076800.0_cp*a**2*b**6*dn**3-24163072.0_cp*a**2*b**6*dn**2-116079104.0_cp &
      &              *a**2*b**6*dn-104441856.0_cp*a**2*b**6+90112.0_cp*a*b**7*dn**5-577536.0_cp &
      &              *a*b**7*dn**4-2267136.0_cp*a*b**7*dn**3+8945664.0_cp*a*b**7*dn**2+31250432.0_cp &
      &              *a*b**7*dn+20705280.0_cp*a*b**7-4096.0_cp*b**8*dn**5+26624.0_cp*b** &
      &              8*dn**4+104448.0_cp*b**8*dn**3-407552.0_cp*b**8*dn**2-1427456.0_cp* &
      &              b**8*dn-946176.0_cp*b**8)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn &
      &              -1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+8)=a**7*b*(105.0_cp*a**6*dm**2*dn**3+78.0_cp*a**6*dm**2*dn**2 &
      &              -735.0_cp*a**6*dm**2*dn-372.0_cp*a**6*dm**2+14.0_cp*a**6*dn**5-94.0_cp &
      &              *a**6*dn**4-850.0_cp*a**6*dn**3+6124.0_cp*a**6*dn**2+2936.0_cp*a**6 &
      &              *dn-25242.0_cp*a**6+896.0_cp*a**5*b*dm**2*dn**3+1536.0_cp*a**5*b*dm**2 &
      &              *dn**2-5504.0_cp*a**5*b*dm**2*dn-6144.0_cp*a**5*b*dm**2+16.0_cp*a** &
      &              5*b*dn**5+878.0_cp*a**5*b*dn**4-15422.0_cp*a**5*b*dn**3+39910.0_cp* &
      &              a**5*b*dn**2+180382.0_cp*a**5*b*dn-301668.0_cp*a**5*b+2040.0_cp*a** &
      &              4*b**2*dm**2*dn**3+8400.0_cp*a**4*b**2*dm**2*dn**2-2040.0_cp*a**4*b**2 &
      &              *dm**2*dn-26640.0_cp*a**4*b**2*dm**2-804.0_cp*a**4*b**2*dn**5+12345.0_cp &
      &              *a**4*b**2*dn**4-44065.0_cp*a**4*b**2*dn**3-280971.0_cp*a**4*b**2*dn**2 &
      &              +1231309.0_cp*a**4*b**2*dn+3837594.0_cp*a**4*b**2-1344.0_cp*a**3*b**3* &
      &              dm**2*dn**3+6528.0_cp*a**3*b**3*dm**2*dn**2+58176.0_cp*a**3*b**3*dm**2 &
      &              *dn+79488.0_cp*a**3*b**3*dm**2-1536.0_cp*a**3*b**3*dn**5-15824.0_cp &
      &              *a**3*b**3*dn**4+354448.0_cp*a**3*b**3*dn**3+78512.0_cp*a**3*b**3*dn**2 &
      &              -7046416.0_cp*a**3*b**3*dn-11367264.0_cp*a**3*b**3-4352.0_cp*a**2*b**4 &
      &              *dm**2*dn**3-29184.0_cp*a**2*b**4*dm**2*dn**2-63232.0_cp*a**2*b**4*dm**2 &
      &              *dn-44544.0_cp*a**2*b**4*dm**2+17568.0_cp*a**2*b**4*dn**5-128796.0_cp &
      &              *a**2*b**4*dn**4-566132.0_cp*a**2*b**4*dn**3+2490948.0_cp*a**2*b**4 &
      &              *dn**2+11132660.0_cp*a**2*b**4*dn+10395384.0_cp*a**2*b**4+1024.0_cp &
      &              *a*b**5*dm**2*dn**3+6144.0_cp*a*b**5*dm**2*dn**2+11264.0_cp*a*b**5*dm**2 &
      &              *dn+6144.0_cp*a*b**5*dm**2-13056.0_cp*a*b**5*dn**5+107040.0_cp*a*b**5* &
      &              dn**4+311520.0_cp*a*b**5*dn**3-1981920.0_cp*a*b**5*dn**2-6116064.0_cp &
      &              *a*b**5*dn-3942720.0_cp*a*b**5+1536.0_cp*b**6*dn**5-12736.0_cp*b**6 &
      &              *dn**4-37056.0_cp*b**6*dn**3+234304.0_cp*b**6*dn**2+723648.0_cp*b** &
      &              6*dn+466560.0_cp*b**6)/(2048.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+9)=a**8*(25.0_cp*a**6*dm**2*dn**3+48.0_cp*a**6*dm**2*dn**2- &
      &              115.0_cp*a**6*dm**2*dn-132.0_cp*a**6*dm**2+a**6*dn**5-213.0_cp*a**6 &
      &              *dn**3+936.0_cp*a**6*dn**2+2696.0_cp*a**6*dn-2544.0_cp*a**6+208.0_cp &
      &              *a**5*b*dm**2*dn**3+576.0_cp*a**5*b*dm**2*dn**2-592.0_cp*a**5*b*dm**2* &
      &              dn-1344.0_cp*a**5*b*dm**2-16.0_cp*a**5*b*dn**5+344.0_cp*a**5*b*dn** &
      &              4-2720.0_cp*a**5*b*dn**3+2152.0_cp*a**5*b*dn**2+43536.0_cp*a**5*b*dn &
      &              +6048.0_cp*a**5*b+46.0_cp*a**4*b**2*dm**2*dn**3+1464.0_cp*a**4*b**2 &
      &              *dm**2*dn**2+3986.0_cp*a**4*b**2*dm**2*dn+24.0_cp*a**4*b**2*dm**2-116.0_cp &
      &              *a**4*b**2*dn**5+988.0_cp*a**4*b**2*dn**4+6540.0_cp*a**4*b**2*dn**3 &
      &              -64156.0_cp*a**4*b**2*dn**2-9064.0_cp*a**4*b**2*dn+618576.0_cp*a**4 &
      &              *b**2-2880.0_cp*a**3*b**3*dm**2*dn**3-7680.0_cp*a**3*b**3*dm**2*dn**2+ &
      &              16320.0_cp*a**3*b**3*dm**2*dn+40320.0_cp*a**3*b**3*dm**2+1408.0_cp* &
      &              a**3*b**3*dn**5-25920.0_cp*a**3*b**3*dn**4+135360.0_cp*a**3*b**3*dn**3 &
      &              +332928.0_cp*a**3*b**3*dn**2-2889088.0_cp*a**3*b**3*dn-5567232.0_cp &
      &              *a**3*b**3-4032.0_cp*a**2*b**4*dm**2*dn**3-29184.0_cp*a**2*b**4*dm**2* &
      &              dn**2-69312.0_cp*a**2*b**4*dm**2*dn-54144.0_cp*a**2*b**4*dm**2+9432.0_cp &
      &              *a**2*b**4*dn**5-81408.0_cp*a**2*b**4*dn**4-355336.0_cp*a**2*b**4*dn**3 &
      &              +1956144.0_cp*a**2*b**4*dn**2+9250432.0_cp*a**2*b**4*dn+9437952.0_cp &
      &              *a**2*b**4+3584.0_cp*a*b**5*dm**2*dn**3+21504.0_cp*a*b**5*dm**2*dn**2+ &
      &              39424.0_cp*a*b**5*dm**2*dn+21504.0_cp*a*b**5*dm**2-16512.0_cp*a*b** &
      &              5*dn**5+165024.0_cp*a*b**5*dn**4+346112.0_cp*a*b**5*dn**3-3535008.0_cp &
      &              *a*b**5*dn**2-10040576.0_cp*a*b**5*dn-6340992.0_cp*a*b**5-256.0_cp* &
      &              b**6*dm**2*dn**3-1536.0_cp*b**6*dm**2*dn**2-2816.0_cp*b**6*dm**2*dn &
      &              -1536.0_cp*b**6*dm**2+3968.0_cp*b**6*dn**5-40032.0_cp*b**6*dn**4-83936.0_cp &
      &              *b**6*dn**3+854304.0_cp*b**6*dn**2+2426976.0_cp*b**6*dn+1532736.0_cp &
      &              *b**6)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+10)=-a**9*b*(63.0_cp*a**4*dm**2*dn**2-15.0_cp*a**4*dm**2 &
      &              *dn-132.0_cp*a**4*dm**2-14.0_cp*a**4*dn**4+288.0_cp*a**4*dn**3-2122.0_cp &
      &              *a**4*dn**2+5430.0_cp*a**4*dn+1746.0_cp*a**4+384.0_cp*a**3*b*dm**2*dn**2 &
      &              +384.0_cp*a**3*b*dm**2*dn-768.0_cp*a**3*b*dm**2-208.0_cp*a**3*b*dn**4+ &
      &              3750.0_cp*a**3*b*dn**3-19148.0_cp*a**3*b*dn**2+4050.0_cp*a**3*b*dn+118020.0_cp &
      &              *a**3*b+232.0_cp*a**2*b**2*dm**2*dn**2+1752.0_cp*a**2*b**2*dm**2*dn &
      &              +2576.0_cp*a**2*b**2*dm**2-492.0_cp*a**2*b**2*dn**4+5577.0_cp*a**2* &
      &              b**2*dn**3+17026.0_cp*a**2*b**2*dn**2-200709.0_cp*a**2*b**2*dn-417034.0_cp &
      &              *a**2*b**2-1216.0_cp*a*b**3*dm**2*dn**2-3648.0_cp*a*b**3*dm**2*dn-2432.0_cp &
      &              *a*b**3*dm**2+3072.0_cp*a*b**3*dn**4-45456.0_cp*a*b**3*dn**3+85856.0_cp &
      &              *a*b**3*dn**2+621840.0_cp*a*b**3*dn+487456.0_cp*a*b**3+256.0_cp*b** &
      &              4*dm**2*dn**2+768.0_cp*b**4*dm**2*dn+512.0_cp*b**4*dm**2-1440.0_cp* &
      &              b**4*dn**4+21444.0_cp*b**4*dn**3-40472.0_cp*b**4*dn**2-293124.0_cp* &
      &              b**4*dn-229768.0_cp*b**4)/(2048.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn &
      &              -1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+11)=-a**10*(43.0_cp*a**4*dm**2*dn**2+21.0_cp*a**4*dm**2*dn &
      &              -70.0_cp*a**4*dm**2-11.0_cp*a**4*dn**4+222.0_cp*a**4*dn**3-1375.0_cp &
      &              *a**4*dn**2+1908.0_cp*a**4*dn+3348.0_cp*a**4+256.0_cp*a**3*b*dm**2*dn**2 &
      &              +384.0_cp*a**3*b*dm**2*dn-256.0_cp*a**3*b*dm**2-128.0_cp*a**3*b*dn**4+ &
      &              2384.0_cp*a**3*b*dn**3-10960.0_cp*a**3*b*dn**2-8000.0_cp*a**3*b*dn+48960.0_cp &
      &              *a**3*b-560.0_cp*a**2*b**2*dm**2*dn**2+288.0_cp*a**2*b**2*dm**2*dn+2816.0_cp &
      &              *a**2*b**2*dm**2+320.0_cp*a**2*b**2*dn**4-7544.0_cp*a**2*b**2*dn**3 &
      &              +59480.0_cp*a**2*b**2*dn**2-102304.0_cp*a**2*b**2*dn-508000.0_cp*a**2* &
      &              b**2-3200.0_cp*a*b**3*dm**2*dn**2-9600.0_cp*a*b**3*dm**2*dn-6400.0_cp &
      &              *a*b**3*dm**2+5376.0_cp*a*b**3*dn**4-89280.0_cp*a*b**3*dn**3+212544.0_cp &
      &              *a*b**3*dn**2+1343232.0_cp*a*b**3*dn+1036032.0_cp*a*b**3+1664.0_cp* &
      &              b**4*dm**2*dn**2+4992.0_cp*b**4*dm**2*dn+3328.0_cp*b**4*dm**2-5136.0_cp &
      &              *b**4*dn**4+85776.0_cp*b**4*dn**3-204128.0_cp*b**4*dn**2-1289856.0_cp &
      &              *b**4*dn-994816.0_cp*b**4)/(16384.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn &
      &              -1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+12)=a**11*b*(39.0_cp*a**2*dm**2*dn-21.0_cp*a**2*dm**2-26.0_cp &
      &              *a**2*dn**3+583.0_cp*a**2*dn**2-3780.0_cp*a**2*dn+5259.0_cp*a**2+128.0_cp &
      &              *a*b*dm**2*dn+128.0_cp*a*b*dm**2-160.0_cp*a*b*dn**3+3268.0_cp*a*b*dn**2 &
      &              -14984.0_cp*a*b*dn-18412.0_cp*a*b-176.0_cp*b**2*dm**2*dn-176.0_cp*b**2 &
      &              *dm**2+360.0_cp*b**2*dn**3-7386.0_cp*b**2*dn**2+33868.0_cp*b**2*dn+41614.0_cp &
      &              *b**2)/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              ))
      stencil(ku+13)=a**12*(5.0_cp*a**2*dm**2*dn-a**2*dm**2-3.0_cp*a**2*dn**3 &
      &              +71.0_cp*a**2*dn**2-458.0_cp*a**2*dn+404.0_cp*a**2+16.0_cp*a*b*dm** &
      &              2*dn+16.0_cp*a*b*dm**2-16.0_cp*a*b*dn**3+356.0_cp*a*b*dn**2-1796.0_cp &
      &              *a*b*dn-2168.0_cp*a*b-82.0_cp*b**2*dm**2*dn-82.0_cp*b**2*dm**2+124.0_cp &
      &              *b**2*dn**3-2770.0_cp*b**2*dn**2+13974.0_cp*b**2*dn+16868.0_cp*b**2 &
      &              )/(8192.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+14)=-a**13*b*(5.0_cp*dm**2-6.0_cp*dn**2+151.0_cp*dn-949.0_cp &
      &              )/(4096.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp))
      stencil(ku+15)=-a**14*(dm**2-dn**2+27.0_cp*dn-182.0_cp)/(16384.0_cp &
      &              *dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp))
      stencil(ku+16:len_stencil) = 0.0_cp

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4rmult4hmult8laplrot
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

      stencil = mirror_stencil(n, len_stencil)

   end function intcheb4rmult4laplrot2
!------------------------------------------------------------------------------
   function intcheb4rmult4hmult8laplrot2(a, b, m, n, len_stencil) result(stencil)

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

      stencil(1:ku-12) = 0.0_cp
      stencil(ku-11)=a**12*(dm-dn-14.0_cp)*(dm+dn+14.0_cp)*(dm**2-dn**2-23.0_cp &
      &              *dn-132.0_cp)/(4096.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-10)=a**11*b*(8.0_cp*dm**4-20.0_cp*dm**2*dn**2-469.0_cp*dm**2 &
      &              *dn-2771.0_cp*dm**2+12.0_cp*dn**4+567.0_cp*dn**3+10012.0_cp*dn**2+78299.0_cp &
      &              *dn+228822.0_cp)/(2048.0_cp*dn*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-9)=-a**10*(6.0_cp*a**2*dm**4*dn-8.0_cp*a**2*dm**2*dn**3-177.0_cp &
      &              *a**2*dm**2*dn**2-1043.0_cp*a**2*dm**2*dn-446.0_cp*a**2*dm**2+2.0_cp &
      &              *a**2*dn**5+83.0_cp*a**2*dn**4+1335.0_cp*a**2*dn**3+10428.0_cp*a**2 &
      &              *dn**2+39836.0_cp*a**2*dn+59856.0_cp*a**2+16.0_cp*a*b*dm**4*dn-16.0_cp &
      &              *a*b*dm**4-32.0_cp*a*b*dm**2*dn**3-660.0_cp*a*b*dm**2*dn**2-3100.0_cp &
      &              *a*b*dm**2*dn+3792.0_cp*a*b*dm**2+16.0_cp*a*b*dn**5+676.0_cp*a*b*dn**4 &
      &              +10508.0_cp*a*b*dn**3+69264.0_cp*a*b*dn**2+136304.0_cp*a*b*dn-216768.0_cp &
      &              *a*b-48.0_cp*b**2*dm**4*dn+48.0_cp*b**2*dm**4+164.0_cp*b**2*dm**2*dn**3 &
      &              +3342.0_cp*b**2*dm**2*dn**2+15342.0_cp*b**2*dm**2*dn-18848.0_cp*b** &
      &              2*dm**2-124.0_cp*b**2*dn**5-5272.0_cp*b**2*dn**4-82346.0_cp*b**2*dn**3 &
      &              -544074.0_cp*b**2*dn**2-1068080.0_cp*b**2*dn+1699896.0_cp*b**2)/(2048.0_cp &
      &              *dn*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-8)=-a**9*b*(72.0_cp*a**2*dm**4*dn+24.0_cp*a**2*dm**4-116.0_cp &
      &              *a**2*dm**2*dn**3-2383.0_cp*a**2*dm**2*dn**2-13680.0_cp*a**2*dm**2*dn &
      &              -11685.0_cp*a**2*dm**2+28.0_cp*a**2*dn**5+1049.0_cp*a**2*dn**4+15899.0_cp &
      &              *a**2*dn**3+127925.0_cp*a**2*dn**2+578321.0_cp*a**2*dn+1196538.0_cp &
      &              *a**2+192.0_cp*a*b*dm**4*dn-192.0_cp*a*b*dm**4-512.0_cp*a*b*dm**2*dn**3 &
      &              -9496.0_cp*a*b*dm**2*dn**2-39456.0_cp*a*b*dm**2*dn+49464.0_cp*a*b*dm**2 &
      &              +320.0_cp*a*b*dn**5+12328.0_cp*a*b*dn**4+174408.0_cp*a*b*dn**3+1040720.0_cp &
      &              *a*b*dn**2+1793080.0_cp*a*b*dn-3020856.0_cp*a*b-128.0_cp*b**2*dm**4 &
      &              *dn+128.0_cp*b**2*dm**4+704.0_cp*b**2*dm**2*dn**3+12868.0_cp*b**2*dm**2 &
      &              *dn**2+52016.0_cp*b**2*dm**2*dn-65588.0_cp*b**2*dm**2-720.0_cp*b**2 &
      &              *dn**5-27936.0_cp*b**2*dn**4-397376.0_cp*b**2*dn**3-2377788.0_cp*b**2* &
      &              dn**2-4087504.0_cp*b**2*dn+6891324.0_cp*b**2)/(2048.0_cp*dn*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-7)=a**8*(33.0_cp*a**4*dm**4*dn**2-33.0_cp*a**4*dm**4*dn-36.0_cp &
      &              *a**4*dm**4-26.0_cp*a**4*dm**2*dn**4-472.0_cp*a**4*dm**2*dn**3-2248.0_cp &
      &              *a**4*dm**2*dn**2+850.0_cp*a**4*dm**2*dn+4464.0_cp*a**4*dm**2+a**4*dn**6 &
      &              +25.0_cp*a**4*dn**5+332.0_cp*a**4*dn**4+3758.0_cp*a**4*dn**3+28368.0_cp &
      &              *a**4*dn**2+96408.0_cp*a**4*dn+46272.0_cp*a**4+160.0_cp*a**3*b*dm** &
      &              4*dn**2-288.0_cp*a**3*b*dm**4*dn-64.0_cp*a**3*b*dm**4-192.0_cp*a**3 &
      &              *b*dm**2*dn**4-3216.0_cp*a**3*b*dm**2*dn**3-11888.0_cp*a**3*b*dm**2 &
      &              *dn**2+23472.0_cp*a**3*b*dm**2*dn+29408.0_cp*a**3*b*dm**2+32.0_cp*a**3 &
      &              *b*dn**6+992.0_cp*a**3*b*dn**5+12368.0_cp*a**3*b*dn**4+82400.0_cp*a**3 &
      &              *b*dn**3+290240.0_cp*a**3*b*dn**2+96128.0_cp*a**3*b*dn-2244096.0_cp &
      &              *a**3*b-96.0_cp*a**2*b**2*dm**4*dn**2-288.0_cp*a**2*b**2*dm**4*dn+960.0_cp &
      &              *a**2*b**2*dm**4+232.0_cp*a**2*b**2*dm**2*dn**4+5048.0_cp*a**2*b**2 &
      &              *dm**2*dn**3+35200.0_cp*a**2*b**2*dm**2*dn**2+41656.0_cp*a**2*b**2*dm**2 &
      &              *dn-268208.0_cp*a**2*b**2*dm**2+88.0_cp*a**2*b**2*dn**6+3064.0_cp*a**2 &
      &              *b**2*dn**5+23460.0_cp*a**2*b**2*dn**4-188388.0_cp*a**2*b**2*dn**3-2734216.0_cp &
      &              *a**2*b**2*dn**2-4737664.0_cp*a**2*b**2*dn+21440256.0_cp*a**2*b**2-768.0_cp &
      &              *a*b**3*dm**4*dn**2+2304.0_cp*a*b**3*dm**4*dn-1536.0_cp*a*b**3*dm** &
      &              4+3200.0_cp*a*b**3*dm**2*dn**4+46192.0_cp*a*b**3*dm**2*dn**3+83984.0_cp &
      &              *a*b**3*dm**2*dn**2-623296.0_cp*a*b**3*dm**2*dn+489920.0_cp*a*b**3*dm**2 &
      &              -2688.0_cp*a*b**3*dn**6-88224.0_cp*a*b**3*dn**5-1006688.0_cp*a*b**3 &
      &              *dn**4-3991152.0_cp*a*b**3*dn**3+3347792.0_cp*a*b**3*dn**2+35899296.0_cp &
      &              *a*b**3*dn-34158336.0_cp*a*b**3+128.0_cp*b**4*dm**4*dn**2-384.0_cp* &
      &              b**4*dm**4*dn+256.0_cp*b**4*dm**4-1664.0_cp*b**4*dm**2*dn**4-23536.0_cp &
      &              *b**4*dm**2*dn**3-39776.0_cp*b**4*dm**2*dn**2+309040.0_cp*b**4*dm** &
      &              2*dn-244064.0_cp*b**4*dm**2+2568.0_cp*b**4*dn**6+85008.0_cp*b**4*dn**5 &
      &              +976496.0_cp*b**4*dn**4+3885536.0_cp*b**4*dn**3-3279704.0_cp*b**4*dn**2 &
      &              -34804112.0_cp*b**4*dn+33134208.0_cp*b**4)/(2048.0_cp*dn*(dn-2.0_cp &
      &              )*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-6)=a**7*b*(280.0_cp*a**4*dm**4*dn**2-72.0_cp*a**4*dm**4*dn &
      &              -496.0_cp*a**4*dm**4-252.0_cp*a**4*dm**2*dn**4-4347.0_cp*a**4*dm**2 &
      &              *dn**3-22016.0_cp*a**4*dm**2*dn**2-12891.0_cp*a**4*dm**2*dn+37874.0_cp &
      &              *a**4*dm**2-12.0_cp*a**4*dn**6-435.0_cp*a**4*dn**5-4183.0_cp*a**4*dn**4 &
      &              +1983.0_cp*a**4*dn**3+242107.0_cp*a**4*dn**2+1424880.0_cp*a**4*dn+2723244.0_cp &
      &              *a**4+1344.0_cp*a**3*b*dm**4*dn**2-1728.0_cp*a**3*b*dm**4*dn-1920.0_cp &
      &              *a**3*b*dm**4-2048.0_cp*a**3*b*dm**2*dn**4-31736.0_cp*a**3*b*dm**2*dn**3 &
      &              -120224.0_cp*a**3*b*dm**2*dn**2+133928.0_cp*a**3*b*dm**2*dn+499696.0_cp &
      &              *a**3*b*dm**2+192.0_cp*a**3*b*dn**6+4776.0_cp*a**3*b*dn**5+70136.0_cp &
      &              *a**3*b*dn**4+786144.0_cp*a**3*b*dn**3+4417000.0_cp*a**3*b*dn**2+3647496.0_cp &
      &              *a**3*b*dn-32539440.0_cp*a**3*b+1152.0_cp*a**2*b**2*dm**4*dn**2-4992.0_cp &
      &              *a**2*b**2*dm**4*dn+5376.0_cp*a**2*b**2*dm**4-3264.0_cp*a**2*b**2*dm**2 &
      &              *dn**4-34236.0_cp*a**2*b**2*dm**2*dn**3+50224.0_cp*a**2*b**2*dm**2*dn**2 &
      &              +770580.0_cp*a**2*b**2*dm**2*dn-1415944.0_cp*a**2*b**2*dm**2+3408.0_cp &
      &              *a**2*b**2*dn**6+99264.0_cp*a**2*b**2*dn**5+923712.0_cp*a**2*b**2*dn**4 &
      &              +1686724.0_cp*a**2*b**2*dn**3-15091536.0_cp*a**2*b**2*dn**2-37507228.0_cp &
      &              *a**2*b**2*dn+103712856.0_cp*a**2*b**2-1024.0_cp*a*b**3*dm**4*dn**2 &
      &              +3072.0_cp*a*b**3*dm**4*dn-2048.0_cp*a*b**3*dm**4+9728.0_cp*a*b**3*dm**2 &
      &              *dn**4+119456.0_cp*a*b**3*dm**2*dn**3+141696.0_cp*a*b**3*dm**2*dn** &
      &              2-1407200.0_cp*a*b**3*dm**2*dn+1136320.0_cp*a*b**3*dm**2-12288.0_cp &
      &              *a*b**3*dn**6-358080.0_cp*a*b**3*dn**5-3584576.0_cp*a*b**3*dn**4-11823776.0_cp &
      &              *a*b**3*dn**3+15946688.0_cp*a*b**3*dn**2+100569632.0_cp*a*b**3*dn-100737600.0_cp &
      &              *a*b**3-2048.0_cp*b**4*dm**2*dn**4-24512.0_cp*b**4*dm**2*dn**3-25984.0_cp &
      &              *b**4*dm**2*dn**2+280256.0_cp*b**4*dm**2*dn-227712.0_cp*b**4*dm**2+5760.0_cp &
      &              *b**4*dn**6+169488.0_cp*b**4*dn**5+1709168.0_cp*b**4*dn**4+5656816.0_cp &
      &              *b**4*dn**3-7674800.0_cp*b**4*dn**2-47827840.0_cp*b**4*dn+47961408.0_cp &
      &              *b**4)/(2048.0_cp*dn*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp &
      &              )*(dn+3.0_cp))
      stencil(ku-5)=-a**6*(110.0_cp*a**6*dm**4*dn**3-330.0_cp*a**6*dm**4*dn**2 &
      &              -140.0_cp*a**6*dm**4*dn+480.0_cp*a**6*dm**4-40.0_cp*a**6*dm**2*dn** &
      &              5-579.0_cp*a**6*dm**2*dn**4-1910.0_cp*a**6*dm**2*dn**3+7449.0_cp*a**6* &
      &              dm**2*dn**2+13620.0_cp*a**6*dm**2*dn-11292.0_cp*a**6*dm**2-6.0_cp*a**6 &
      &              *dn**7-141.0_cp*a**6*dn**6-1100.0_cp*a**6*dn**5-2263.0_cp*a**6*dn** &
      &              4+14066.0_cp*a**6*dn**3+80740.0_cp*a**6*dn**2+87480.0_cp*a**6*dn-279648.0_cp &
      &              *a**6+720.0_cp*a**5*b*dm**4*dn**3-2592.0_cp*a**5*b*dm**4*dn**2+240.0_cp &
      &              *a**5*b*dm**4*dn+3168.0_cp*a**5*b*dm**4-416.0_cp*a**5*b*dm**2*dn**5 &
      &              -5340.0_cp*a**5*b*dm**2*dn**4-11608.0_cp*a**5*b*dm**2*dn**3+82572.0_cp &
      &              *a**5*b*dm**2*dn**2+90408.0_cp*a**5*b*dm**2*dn-167328.0_cp*a**5*b*dm**2 &
      &              -48.0_cp*a**5*b*dn**7-1188.0_cp*a**5*b*dn**6-8104.0_cp*a**5*b*dn**5 &
      &              +6916.0_cp*a**5*b*dn**4+273544.0_cp*a**5*b*dn**3+798800.0_cp*a**5*b &
      &              *dn**2-1320480.0_cp*a**5*b*dn-8233344.0_cp*a**5*b+912.0_cp*a**4*b** &
      &              2*dm**4*dn**3-6048.0_cp*a**4*b**2*dm**4*dn**2+10032.0_cp*a**4*b**2*dm**4 &
      &              *dn-288.0_cp*a**4*b**2*dm**4-812.0_cp*a**4*b**2*dm**2*dn**5-6114.0_cp &
      &              *a**4*b**2*dm**2*dn**4+30080.0_cp*a**4*b**2*dm**2*dn**3+236010.0_cp &
      &              *a**4*b**2*dm**2*dn**2-433308.0_cp*a**4*b**2*dm**2*dn-943776.0_cp*a**4 &
      &              *b**2*dm**2+180.0_cp*a**4*b**2*dn**7+3708.0_cp*a**4*b**2*dn**6+41966.0_cp &
      &              *a**4*b**2*dn**5+335816.0_cp*a**4*b**2*dn**4+228046.0_cp*a**4*b**2*dn**3 &
      &              -12288332.0_cp*a**4*b**2*dn**2-21173400.0_cp*a**4*b**2*dn+127462320.0_cp &
      &              *a**4*b**2-2048.0_cp*a**3*b**3*dm**4*dn**3+3072.0_cp*a**3*b**3*dm** &
      &              4*dn**2+23552.0_cp*a**3*b**3*dm**4*dn-43008.0_cp*a**3*b**3*dm**4+5120.0_cp &
      &              *a**3*b**3*dm**2*dn**5+72096.0_cp*a**3*b**3*dm**2*dn**4+205888.0_cp &
      &              *a**3*b**3*dm**2*dn**3-1190112.0_cp*a**3*b**3*dm**2*dn**2-3727872.0_cp &
      &              *a**3*b**3*dm**2*dn+9251712.0_cp*a**3*b**3*dm**2+2560.0_cp*a**3*b** &
      &              3*dn**7+60864.0_cp*a**3*b**3*dn**6+186112.0_cp*a**3*b**3*dn**5-3966432.0_cp &
      &              *a**3*b**3*dn**4-20849600.0_cp*a**3*b**3*dn**3+59435616.0_cp*a**3*b**3 &
      &              *dn**2+271767168.0_cp*a**3*b**3*dn-561195648.0_cp*a**3*b**3-3072.0_cp &
      &              *a**2*b**4*dm**4*dn**3+19968.0_cp*a**2*b**4*dm**4*dn**2-41472.0_cp* &
      &              a**2*b**4*dm**4*dn+27648.0_cp*a**2*b**4*dm**4+19456.0_cp*a**2*b**4*dm**2 &
      &              *dn**5+128352.0_cp*a**2*b**4*dm**2*dn**4-634048.0_cp*a**2*b**4*dm** &
      &              2*dn**3-2297568.0_cp*a**2*b**4*dm**2*dn**2+11086656.0_cp*a**2*b**4*dm**2 &
      &              *dn-10586880.0_cp*a**2*b**4*dm**2-24000.0_cp*a**2*b**4*dn**7-534960.0_cp &
      &              *a**2*b**4*dn**6-3149056.0_cp*a**2*b**4*dn**5+5856608.0_cp*a**2*b** &
      &              4*dn**4+80152576.0_cp*a**2*b**4*dn**3-36374384.0_cp*a**2*b**4*dn**2 &
      &              -652091328.0_cp*a**2*b**4*dn+852833088.0_cp*a**2*b**4-14336.0_cp*a* &
      &              b**5*dm**2*dn**5-101376.0_cp*a*b**5*dm**2*dn**4+359168.0_cp*a*b**5*dm**2 &
      &              *dn**3+1669632.0_cp*a*b**5*dm**2*dn**2-5558016.0_cp*a*b**5*dm**2*dn &
      &              +3644928.0_cp*a*b**5*dm**2+33024.0_cp*a*b**5*dn**7+742464.0_cp*a*b**5* &
      &              dn**6+4720256.0_cp*a*b**5*dn**5-2794496.0_cp*a*b**5*dn**4-95244416.0_cp &
      &              *a*b**5*dn**3-38876992.0_cp*a*b**5*dn**2+627962112.0_cp*a*b**5*dn-496541952.0_cp &
      &              *a*b**5+1024.0_cp*b**6*dm**2*dn**5+6912.0_cp*b**6*dm**2*dn**4-25856.0_cp &
      &              *b**6*dm**2*dn**3-109824.0_cp*b**6*dm**2*dn**2+375040.0_cp*b**6*dm**2* &
      &              dn-247296.0_cp*b**6*dm**2-7936.0_cp*b**6*dn**7-180672.0_cp*b**6*dn**6- &
      &              1160704.0_cp*b**6*dn**5+676224.0_cp*b**6*dn**4+23441408.0_cp*b**6*dn**3 &
      &              +8957760.0_cp*b**6*dn**2-152734464.0_cp*b**6*dn+121008384.0_cp*b**6 &
      &              )/(2048.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn &
      &              +2.0_cp)*(dn+3.0_cp))
      stencil(ku-4)=-a**5*b*(600.0_cp*a**6*dm**4*dn**3-1008.0_cp*a**6*dm**4* &
      &              dn**2-3000.0_cp*a**6*dm**4*dn+2832.0_cp*a**6*dm**4-220.0_cp*a**6*dm**2 &
      &              *dn**5-3329.0_cp*a**6*dm**2*dn**4-14155.0_cp*a**6*dm**2*dn**3+19583.0_cp &
      &              *a**6*dm**2*dn**2+107135.0_cp*a**6*dm**2*dn+22362.0_cp*a**6*dm**2-60.0_cp &
      &              *a**6*dn**7-1125.0_cp*a**6*dn**6-8334.0_cp*a**6*dn**5-31144.0_cp*a**6* &
      &              dn**4-13926.0_cp*a**6*dn**3+454177.0_cp*a**6*dn**2+2205888.0_cp*a** &
      &              6*dn+3533436.0_cp*a**6+3840.0_cp*a**5*b*dm**4*dn**3-9216.0_cp*a**5* &
      &              b*dm**4*dn**2-15360.0_cp*a**5*b*dm**4*dn+25344.0_cp*a**5*b*dm**4-2560.0_cp &
      &              *a**5*b*dm**2*dn**5-33344.0_cp*a**5*b*dm**2*dn**4-101568.0_cp*a**5* &
      &              b*dm**2*dn**3+334112.0_cp*a**5*b*dm**2*dn**2+1076992.0_cp*a**5*b*dm**2 &
      &              *dn-172704.0_cp*a**5*b*dm**2-768.0_cp*a**5*b*dn**7-14880.0_cp*a**5* &
      &              b*dn**6-98816.0_cp*a**5*b*dn**5-203264.0_cp*a**5*b*dn**4+1911872.0_cp &
      &              *a**5*b*dn**3+15791360.0_cp*a**5*b*dn**2+3304320.0_cp*a**5*b*dn-150651936.0_cp &
      &              *a**5*b+7680.0_cp*a**4*b**2*dm**4*dn**3-30720.0_cp*a**4*b**2*dm**4*dn**2 &
      &              +69120.0_cp*a**4*b**2*dm**4-10496.0_cp*a**4*b**2*dm**2*dn**5-105040.0_cp &
      &              *a**4*b**2*dm**2*dn**4-54896.0_cp*a**4*b**2*dm**2*dn**3+2164960.0_cp &
      &              *a**4*b**2*dm**2*dn**2+1736080.0_cp*a**4*b**2*dm**2*dn-12151920.0_cp &
      &              *a**4*b**2*dm**2-2880.0_cp*a**4*b**2*dn**7-59040.0_cp*a**4*b**2*dn**6- &
      &              6784.0_cp*a**4*b**2*dn**5+4894960.0_cp*a**4*b**2*dn**4+13695760.0_cp &
      &              *a**4*b**2*dn**3-99050560.0_cp*a**4*b**2*dn**2-202721136.0_cp*a**4* &
      &              b**2*dn+784328400.0_cp*a**4*b**2+3072.0_cp*a**3*b**3*dm**4*dn**3-30720.0_cp &
      &              *a**3*b**3*dm**4*dn**2+95232.0_cp*a**3*b**3*dm**4*dn-92160.0_cp*a** &
      &              3*b**3*dm**4-8704.0_cp*a**3*b**3*dm**2*dn**5+47200.0_cp*a**3*b**3*dm**2 &
      &              *dn**4+963616.0_cp*a**3*b**3*dm**2*dn**3-1027360.0_cp*a**3*b**3*dm**2* &
      &              dn**2-14403360.0_cp*a**3*b**3*dm**2*dn+24730560.0_cp*a**3*b**3*dm** &
      &              2+30720.0_cp*a**3*b**3*dn**7+576960.0_cp*a**3*b**3*dn**6+1892864.0_cp &
      &              *a**3*b**3*dn**5-18138400.0_cp*a**3*b**3*dn**4-83616224.0_cp*a**3*b**3 &
      &              *dn**3+227908000.0_cp*a**3*b**3*dn**2+781839072.0_cp*a**3*b**3*dn-1617595200.0_cp &
      &              *a**3*b**3+38912.0_cp*a**2*b**4*dm**2*dn**5+181696.0_cp*a**2*b**4*dm**2 &
      &              *dn**4-1129280.0_cp*a**2*b**4*dm**2*dn**3-2291776.0_cp*a**2*b**4*dm**2 &
      &              *dn**2+12894528.0_cp*a**2*b**4*dm**2*dn-11740032.0_cp*a**2*b**4*dm**2- &
      &              81792.0_cp*a**2*b**4*dn**7-1520784.0_cp*a**2*b**4*dn**6-6862016.0_cp &
      &              *a**2*b**4*dn**5+20296032.0_cp*a**2*b**4*dn**4+158494592.0_cp*a**2* &
      &              b**4*dn**3-136013712.0_cp*a**2*b**4*dn**2-1022954304.0_cp*a**2*b**4 &
      &              *dn+1324654272.0_cp*a**2*b**4-8192.0_cp*a*b**5*dm**2*dn**5-39424.0_cp &
      &              *a*b**5*dm**2*dn**4+204288.0_cp*a*b**5*dm**2*dn**3+497152.0_cp*a*b**5* &
      &              dm**2*dn**2-2076160.0_cp*a*b**5*dm**2*dn+1422336.0_cp*a*b**5*dm**2+52224.0_cp &
      &              *a*b**5*dn**7+984192.0_cp*a*b**5*dn**6+4791296.0_cp*a*b**5*dn**5-9059584.0_cp &
      &              *a*b**5*dn**4-95214080.0_cp*a*b**5*dn**3+19867264.0_cp*a*b**5*dn**2 &
      &              +525623808.0_cp*a*b**5*dn-447045120.0_cp*a*b**5-6144.0_cp*b**6*dn** &
      &              7-117504.0_cp*b**6*dn**6-578560.0_cp*b**6*dn**5+1094144.0_cp*b**6*dn**4 &
      &              +11474944.0_cp*b**6*dn**3-2829056.0_cp*b**6*dn**2-62131200.0_cp*b** &
      &              6*dn+53093376.0_cp*b**6)/(2048.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku-3)=a**4*(495.0_cp*a**8*dm**4*dn**3-990.0_cp*a**8*dm**4*dn**2 &
      &              -1755.0_cp*a**8*dm**4*dn+1890.0_cp*a**8*dm**4-30.0_cp*a**8*dm**2*dn**5 &
      &              -867.0_cp*a**8*dm**2*dn**4-5620.0_cp*a**8*dm**2*dn**3+6327.0_cp*a** &
      &              8*dm**2*dn**2+30790.0_cp*a**8*dm**2*dn+1584.0_cp*a**8*dm**2-17.0_cp &
      &              *a**8*dn**7-223.0_cp*a**8*dn**6-1365.0_cp*a**8*dn**5-5749.0_cp*a**8 &
      &              *dn**4-13018.0_cp*a**8*dn**3+12620.0_cp*a**8*dn**2+180120.0_cp*a**8 &
      &              *dn+431136.0_cp*a**8+3840.0_cp*a**7*b*dm**4*dn**3-9216.0_cp*a**7*b*dm**4 &
      &              *dn**2-11520.0_cp*a**7*b*dm**4*dn+17664.0_cp*a**7*b*dm**4-512.0_cp* &
      &              a**7*b*dm**2*dn**5-9920.0_cp*a**7*b*dm**2*dn**4-48896.0_cp*a**7*b*dm**2 &
      &              *dn**3+91136.0_cp*a**7*b*dm**2*dn**2+320896.0_cp*a**7*b*dm**2*dn-15744.0_cp &
      &              *a**7*b*dm**2-256.0_cp*a**7*b*dn**7-3584.0_cp*a**7*b*dn**6-20864.0_cp &
      &              *a**7*b*dn**5-68544.0_cp*a**7*b*dn**4-36864.0_cp*a**7*b*dn**3+902912.0_cp &
      &              *a**7*b*dn**2+4399104.0_cp*a**7*b*dn+7575552.0_cp*a**7*b+9984.0_cp* &
      &              a**6*b**2*dm**4*dn**3-36864.0_cp*a**6*b**2*dm**4*dn**2-5376.0_cp*a**6* &
      &              b**2*dm**4*dn+66816.0_cp*a**6*b**2*dm**4-2624.0_cp*a**6*b**2*dm**2*dn**5 &
      &              -35072.0_cp*a**6*b**2*dm**2*dn**4-88512.0_cp*a**6*b**2*dm**2*dn**3+656672.0_cp &
      &              *a**6*b**2*dm**2*dn**2+1111232.0_cp*a**6*b**2*dm**2*dn-1022016.0_cp &
      &              *a**6*b**2*dm**2-1472.0_cp*a**6*b**2*dn**7-21952.0_cp*a**6*b**2*dn**6- &
      &              99936.0_cp*a**6*b**2*dn**5-7872.0_cp*a**6*b**2*dn**4+3643232.0_cp*a**6 &
      &              *b**2*dn**3+21614752.0_cp*a**6*b**2*dn**2-30415104.0_cp*a**6*b**2*dn &
      &              -267545088.0_cp*a**6*b**2+6144.0_cp*a**5*b**3*dm**4*dn**3-67584.0_cp &
      &              *a**5*b**3*dm**4*dn**2+129024.0_cp*a**5*b**3*dm**4*dn+55296.0_cp*a**5* &
      &              b**3*dm**4-1024.0_cp*a**5*b**3*dm**2*dn**5+33344.0_cp*a**5*b**3*dm**2* &
      &              dn**4+509184.0_cp*a**5*b**3*dm**2*dn**3+2094016.0_cp*a**5*b**3*dm** &
      &              2*dn**2-4554752.0_cp*a**5*b**3*dm**2*dn-21381888.0_cp*a**5*b**3*dm**2- &
      &              1024.0_cp*a**5*b**3*dn**7-19328.0_cp*a**5*b**3*dn**6+722944.0_cp*a**5* &
      &              b**3*dn**5+9485120.0_cp*a**5*b**3*dn**4-2766336.0_cp*a**5*b**3*dn** &
      &              3-262940352.0_cp*a**5*b**3*dn**2-127774464.0_cp*a**5*b**3*dn+1896837120.0_cp &
      &              *a**5*b**3-9216.0_cp*a**4*b**4*dm**4*dn**3-21504.0_cp*a**4*b**4*dm**4* &
      &              dn**2+297984.0_cp*a**4*b**4*dm**4*dn-451584.0_cp*a**4*b**4*dm**4+37888.0_cp &
      &              *a**4*b**4*dm**2*dn**5+570304.0_cp*a**4*b**4*dm**2*dn**4+1406336.0_cp &
      &              *a**4*b**4*dm**2*dn**3-11473216.0_cp*a**4*b**4*dm**2*dn**2-23014272.0_cp &
      &              *a**4*b**4*dm**2*dn+78929280.0_cp*a**4*b**4*dm**2+80320.0_cp*a**4*b**4 &
      &              *dn**7+1229600.0_cp*a**4*b**4*dn**6+748032.0_cp*a**4*b**4*dn**5-54200384.0_cp &
      &              *a**4*b**4*dn**4-108542272.0_cp*a**4*b**4*dn**3+810423072.0_cp*a**4 &
      &              *b**4*dn**2+1108388736.0_cp*a**4*b**4*dn-4551911424.0_cp*a**4*b**4+90112.0_cp &
      &              *a**3*b**5*dm**2*dn**5+16384.0_cp*a**3*b**5*dm**2*dn**4-3644416.0_cp &
      &              *a**3*b**5*dm**2*dn**3+2920448.0_cp*a**3*b**5*dm**2*dn**2+34563072.0_cp &
      &              *a**3*b**5*dm**2*dn-54798336.0_cp*a**3*b**5*dm**2-259072.0_cp*a**3* &
      &              b**5*dn**7-3849728.0_cp*a**3*b**5*dn**6-8172032.0_cp*a**3*b**5*dn** &
      &              5+97853440.0_cp*a**3*b**5*dn**4+298031616.0_cp*a**3*b**5*dn**3-1015069184.0_cp &
      &              *a**3*b**5*dn**2-2129931264.0_cp*a**3*b**5*dn+4911280128.0_cp*a**3* &
      &              b**5-53248.0_cp*a**2*b**6*dm**2*dn**5-133120.0_cp*a**2*b**6*dm**2*dn**4 &
      &              +1313792.0_cp*a**2*b**6*dm**2*dn**3+937984.0_cp*a**2*b**6*dm**2*dn**2- &
      &              9763840.0_cp*a**2*b**6*dm**2*dn+9099264.0_cp*a**2*b**6*dm**2+302080.0_cp &
      &              *a**2*b**6*dn**7+4516352.0_cp*a**2*b**6*dn**6+13649920.0_cp*a**2*b**6* &
      &              dn**5-70317056.0_cp*a**2*b**6*dn**4-307349504.0_cp*a**2*b**6*dn**3+461376000.0_cp &
      &              *a**2*b**6*dn**2+1597879296.0_cp*a**2*b**6*dn-2221903872.0_cp*a**2* &
      &              b**6-90112.0_cp*a*b**7*dn**7-1372160.0_cp*a*b**7*dn**6-4550656.0_cp &
      &              *a*b**7*dn**5+17854464.0_cp*a*b**7*dn**4+91000832.0_cp*a*b**7*dn**3 &
      &              -82739200.0_cp*a*b**7*dn**2-417644544.0_cp*a*b**7*dn+397541376.0_cp &
      &              *a*b**7+4096.0_cp*b**8*dn**7+63488.0_cp*b**8*dn**6+212992.0_cp*b**8 &
      &              *dn**5-839680.0_cp*b**8*dn**4-4222976.0_cp*b**8*dn**3+4093952.0_cp* &
      &              b**8*dn**2+18751488.0_cp*b**8*dn-18063360.0_cp*b**8)/(4096.0_cp*dn* &
      &              (dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku-2)=a**4*b*(360.0_cp*a**7*dm**4*dn**2-936.0_cp*a**7*dm**4 &
      &              *dn+48.0_cp*a**7*dm**4-4.0_cp*a**7*dm**2*dn**4-621.0_cp*a**7*dm**2*dn**3 &
      &              -3422.0_cp*a**7*dm**2*dn**2+4737.0_cp*a**7*dm**2*dn+6078.0_cp*a**7*dm**2 &
      &              -4.0_cp*a**7*dn**6+15.0_cp*a**7*dn**5-367.0_cp*a**7*dn**4-3827.0_cp &
      &              *a**7*dn**3-19849.0_cp*a**7*dn**2-73036.0_cp*a**7*dn-124404.0_cp*a**7+ &
      &              2688.0_cp*a**6*b*dm**4*dn**2-8064.0_cp*a**6*b*dm**4*dn+1920.0_cp*a**6* &
      &              b*dm**4-256.0_cp*a**6*b*dm**2*dn**4-7488.0_cp*a**6*b*dm**2*dn**3-30080.0_cp &
      &              *a**6*b*dm**2*dn**2+54864.0_cp*a**6*b*dm**2*dn+88560.0_cp*a**6*b*dm**2 &
      &              -128.0_cp*a**6*b*dn**6-528.0_cp*a**6*b*dn**5-10960.0_cp*a**6*b*dn** &
      &              4-75184.0_cp*a**6*b*dn**3-636336.0_cp*a**6*b*dn**2-1163312.0_cp*a** &
      &              6*b*dn+19410864.0_cp*a**6*b+7424.0_cp*a**5*b**2*dm**4*dn**2-28416.0_cp &
      &              *a**5*b**2*dm**4*dn+17152.0_cp*a**5*b**2*dm**4-2432.0_cp*a**5*b**2*dm**2 &
      &              *dn**4-38184.0_cp*a**5*b**2*dm**2*dn**3-86160.0_cp*a**5*b**2*dm**2*dn**2 &
      &              +359184.0_cp*a**5*b**2*dm**2*dn+1237928.0_cp*a**5*b**2*dm**2-1760.0_cp &
      &              *a**5*b**2*dn**6-12432.0_cp*a**5*b**2*dn**5-236176.0_cp*a**5*b**2*dn**4 &
      &              -1204248.0_cp*a**5*b**2*dn**3+10182144.0_cp*a**5*b**2*dn**2+33471600.0_cp &
      &              *a**5*b**2*dn-152092728.0_cp*a**5*b**2+7680.0_cp*a**4*b**3*dm**4*dn**2 &
      &              -44544.0_cp*a**4*b**3*dm**4*dn+64512.0_cp*a**4*b**3*dm**4-9472.0_cp &
      &              *a**4*b**3*dm**2*dn**4-88176.0_cp*a**4*b**3*dm**2*dn**3+299680.0_cp &
      &              *a**4*b**3*dm**2*dn**2+2651376.0_cp*a**4*b**3*dm**2*dn-7503264.0_cp &
      &              *a**4*b**3*dm**2-23552.0_cp*a**4*b**3*dn**6-198240.0_cp*a**4*b**3*dn**5 &
      &              +1404960.0_cp*a**4*b**3*dn**4+10470896.0_cp*a**4*b**3*dn**3-36396928.0_cp &
      &              *a**4*b**3*dn**2-145884560.0_cp*a**4*b**3*dn+434051808.0_cp*a**4*b**3- &
      &              1024.0_cp*a**3*b**4*dm**2*dn**4+204960.0_cp*a**3*b**4*dm**2*dn**3+34816.0_cp &
      &              *a**3*b**4*dm**2*dn**2-4366176.0_cp*a**3*b**4*dm**2*dn+7334208.0_cp &
      &              *a**3*b**4*dm**2+89408.0_cp*a**3*b**4*dn**6+738024.0_cp*a**3*b**4*dn**5 &
      &              -2985624.0_cp*a**3*b**4*dn**4-25985576.0_cp*a**3*b**4*dn**3+52271416.0_cp &
      &              *a**3*b**4*dn**2+256232576.0_cp*a**3*b**4*dn-540212640.0_cp*a**3*b**4+ &
      &              20480.0_cp*a**2*b**5*dm**2*dn**4-82176.0_cp*a**2*b**5*dm**2*dn**3-253952.0_cp &
      &              *a**2*b**5*dm**2*dn**2+1499904.0_cp*a**2*b**5*dm**2*dn-1654272.0_cp &
      &              *a**2*b**5*dm**2-134656.0_cp*a**2*b**5*dn**6-1095360.0_cp*a**2*b**5 &
      &              *dn**5+2403904.0_cp*a**2*b**5*dn**4+27233472.0_cp*a**2*b**5*dn**3-29977920.0_cp &
      &              *a**2*b**5*dn**2-203132928.0_cp*a**2*b**5*dn+313516800.0_cp*a**2*b**5+ &
      &              72704.0_cp*a*b**6*dn**6+606336.0_cp*a*b**6*dn**5-679808.0_cp*a*b**6 &
      &              *dn**4-11779200.0_cp*a*b**6*dn**3+4811136.0_cp*a*b**6*dn**2+67651584.0_cp &
      &              *a*b**6*dn-73492992.0_cp*a*b**6-8192.0_cp*b**7*dn**6-70656.0_cp*b** &
      &              7*dn**5+56320.0_cp*b**7*dn**4+1274880.0_cp*b**7*dn**3-232448.0_cp*b**7 &
      &              *dn**2-6365184.0_cp*b**7*dn+5345280.0_cp*b**7)/(1024.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku-1)=-a**4*(198.0_cp*a**8*dm**4*dn**2-594.0_cp*a**8*dm**4*dn &
      &              +216.0_cp*a**8*dm**4+24.0_cp*a**8*dm**2*dn**4-117.0_cp*a**8*dm**2*dn**3 &
      &              -1268.0_cp*a**8*dm**2*dn**2+2575.0_cp*a**8*dm**2*dn+294.0_cp*a**8*dm**2 &
      &              +2.0_cp*a**8*dn**6+25.0_cp*a**8*dn**5+34.0_cp*a**8*dn**4-261.0_cp*a**8 &
      &              *dn**3-2060.0_cp*a**8*dn**2-2412.0_cp*a**8*dn-18000.0_cp*a**8+1680.0_cp &
      &              *a**7*b*dm**4*dn**2-5376.0_cp*a**7*b*dm**4*dn+2352.0_cp*a**7*b*dm** &
      &              4+224.0_cp*a**7*b*dm**2*dn**4-1428.0_cp*a**7*b*dm**2*dn**3-11824.0_cp &
      &              *a**7*b*dm**2*dn**2+25300.0_cp*a**7*b*dm**2*dn+4896.0_cp*a**7*b*dm**2+ &
      &              16.0_cp*a**7*b*dn**6+260.0_cp*a**7*b*dn**5-80.0_cp*a**7*b*dn**4-5844.0_cp &
      &              *a**7*b*dn**3-33168.0_cp*a**7*b*dn**2-64368.0_cp*a**7*b*dn-235584.0_cp &
      &              *a**7*b+5712.0_cp*a**6*b**2*dm**4*dn**2-21504.0_cp*a**6*b**2*dm**4*dn &
      &              +13680.0_cp*a**6*b**2*dm**4+836.0_cp*a**6*b**2*dm**2*dn**4-7638.0_cp &
      &              *a**6*b**2*dm**2*dn**3-43348.0_cp*a**6*b**2*dm**2*dn**2+121894.0_cp &
      &              *a**6*b**2*dm**2*dn+46632.0_cp*a**6*b**2*dm**2-28.0_cp*a**6*b**2*dn**6 &
      &              +748.0_cp*a**6*b**2*dn**5-7682.0_cp*a**6*b**2*dn**4-59220.0_cp*a**6 &
      &              *b**2*dn**3-561890.0_cp*a**6*b**2*dn**2-678744.0_cp*a**6*b**2*dn+19564056.0_cp &
      &              *a**6*b**2+9216.0_cp*a**5*b**3*dm**4*dn**2-47616.0_cp*a**5*b**3*dm**4* &
      &              dn+50688.0_cp*a**5*b**3*dm**4+1536.0_cp*a**5*b**3*dm**2*dn**4-21360.0_cp &
      &              *a**5*b**3*dm**2*dn**3-49472.0_cp*a**5*b**3*dm**2*dn**2+418768.0_cp &
      &              *a**5*b**3*dm**2*dn+917856.0_cp*a**5*b**3*dm**2-1280.0_cp*a**5*b**3 &
      &              *dn**6-4768.0_cp*a**5*b**3*dn**5-192064.0_cp*a**5*b**3*dn**4-729008.0_cp &
      &              *a**5*b**3*dn**3+12054336.0_cp*a**5*b**3*dn**2+24059376.0_cp*a**5*b**3 &
      &              *dn-182951712.0_cp*a**5*b**3+5632.0_cp*a**4*b**4*dm**4*dn**2-53504.0_cp &
      &              *a**4*b**4*dm**4*dn+111360.0_cp*a**4*b**4*dm**4+2560.0_cp*a**4*b**4 &
      &              *dm**2*dn**4-15120.0_cp*a**4*b**4*dm**2*dn**3+415872.0_cp*a**4*b**4 &
      &              *dm**2*dn**2+1855312.0_cp*a**4*b**4*dm**2*dn-9724992.0_cp*a**4*b**4 &
      &              *dm**2-20512.0_cp*a**4*b**4*dn**6-115640.0_cp*a**4*b**4*dn**5+2031088.0_cp &
      &              *a**4*b**4*dn**4+8266000.0_cp*a**4*b**4*dn**3-56665888.0_cp*a**4*b**4* &
      &              dn**2-121347096.0_cp*a**4*b**4*dn+567446544.0_cp*a**4*b**4+23552.0_cp &
      &              *a**3*b**5*dm**2*dn**4+244224.0_cp*a**3*b**5*dm**2*dn**3-776576.0_cp &
      &              *a**3*b**5*dm**2*dn**2-4463616.0_cp*a**3*b**5*dm**2*dn+11878272.0_cp &
      &              *a**3*b**5*dm**2+116608.0_cp*a**3*b**5*dn**6+644000.0_cp*a**3*b**5*dn**5 &
      &              -5371520.0_cp*a**3*b**5*dn**4-23629056.0_cp*a**3*b**5*dn**3+102374976.0_cp &
      &              *a**3*b**5*dn**2+235108704.0_cp*a**3*b**5*dn-795122496.0_cp*a**3*b**5+ &
      &              7680.0_cp*a**2*b**6*dm**2*dn**4-157824.0_cp*a**2*b**6*dm**2*dn**3+131968.0_cp &
      &              *a**2*b**6*dm**2*dn**2+1983616.0_cp*a**2*b**6*dm**2*dn-3499392.0_cp &
      &              *a**2*b**6*dm**2-186496.0_cp*a**2*b**6*dn**6-1013600.0_cp*a**2*b**6 &
      &              *dn**5+6138304.0_cp*a**2*b**6*dn**4+28473408.0_cp*a**2*b**6*dn**3-85512320.0_cp &
      &              *a**2*b**6*dn**2-213557472.0_cp*a**2*b**6*dn+526559040.0_cp*a**2*b**6+ &
      &              139264.0_cp*a*b**7*dn**6+763904.0_cp*a*b**7*dn**5-2981888.0_cp*a*b**7* &
      &              dn**4-15616000.0_cp*a*b**7*dn**3+29966336.0_cp*a*b**7*dn**2+86900736.0_cp &
      &              *a*b**7*dn-154386432.0_cp*a*b**7-28672.0_cp*b**8*dn**6-164864.0_cp* &
      &              b**8*dn**5+460800.0_cp*b**8*dn**4+2754560.0_cp*b**8*dn**3-3350528.0_cp &
      &              *b**8*dn**2-11682816.0_cp*b**8*dn+14469120.0_cp*b**8)/(1024.0_cp*dn &
      &              *(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+3.0_cp))
      stencil(ku)=-a**4*b*(168.0_cp*a**7*dm**4*dn**2+840.0_cp*a**7*dm**4*dn &
      &              -2352.0_cp*a**7*dm**4+28.0_cp*a**7*dm**2*dn**4+21.0_cp*a**7*dm**2*dn**3 &
      &              -1574.0_cp*a**7*dm**2*dn**2-5069.0_cp*a**7*dm**2*dn+8058.0_cp*a**7*dm**2 &
      &              +12.0_cp*a**7*dn**6+45.0_cp*a**7*dn**5+321.0_cp*a**7*dn**4+1515.0_cp &
      &              *a**7*dn**3+8031.0_cp*a**7*dn**2+16200.0_cp*a**7*dn+41652.0_cp*a**7 &
      &              +1344.0_cp*a**6*b*dm**4*dn**2+6720.0_cp*a**6*b*dm**4*dn-20736.0_cp* &
      &              a**6*b*dm**4+256.0_cp*a**6*b*dm**2*dn**4-24.0_cp*a**6*b*dm**2*dn**3 &
      &              -16720.0_cp*a**6*b*dm**2*dn**2-47544.0_cp*a**6*b*dm**2*dn+53184.0_cp &
      &              *a**6*b*dm**2+192.0_cp*a**6*b*dn**6+696.0_cp*a**6*b*dn**5+6728.0_cp &
      &              *a**6*b*dn**4+19056.0_cp*a**6*b*dn**3+374440.0_cp*a**6*b*dn**2+222120.0_cp &
      &              *a**6*b*dn-20625696.0_cp*a**6*b+4224.0_cp*a**5*b**2*dm**4*dn**2+21120.0_cp &
      &              *a**5*b**2*dm**4*dn-78336.0_cp*a**5*b**2*dm**4+832.0_cp*a**5*b**2*dm**2 &
      &              *dn**4-2292.0_cp*a**5*b**2*dm**2*dn**3-90456.0_cp*a**5*b**2*dm**2*dn**2 &
      &              -208372.0_cp*a**5*b**2*dm**2*dn-528048.0_cp*a**5*b**2*dm**2+1744.0_cp &
      &              *a**5*b**2*dn**6+5968.0_cp*a**5*b**2*dn**5+167920.0_cp*a**5*b**2*dn**4 &
      &              +300604.0_cp*a**5*b**2*dn**3-13363528.0_cp*a**5*b**2*dn**2-12709284.0_cp &
      &              *a**5*b**2*dn+200000016.0_cp*a**5*b**2+5632.0_cp*a**4*b**3*dm**4*dn**2 &
      &              +28160.0_cp*a**4*b**3*dm**4*dn-150528.0_cp*a**4*b**3*dm**4-256.0_cp &
      &              *a**4*b**3*dm**2*dn**4-20496.0_cp*a**4*b**3*dm**2*dn**3-534496.0_cp &
      &              *a**4*b**3*dm**2*dn**2-900016.0_cp*a**4*b**3*dm**2*dn+11231520.0_cp &
      &              *a**4*b**3*dm**2+20480.0_cp*a**4*b**3*dn**6+63200.0_cp*a**4*b**3*dn**5 &
      &              -2372064.0_cp*a**4*b**3*dn**4-4572144.0_cp*a**4*b**3*dn**3+70316288.0_cp &
      &              *a**4*b**3*dn**2+68285136.0_cp*a**4*b**3*dn-666389088.0_cp*a**4*b** &
      &              3-21504.0_cp*a**3*b**4*dm**2*dn**4-153504.0_cp*a**3*b**4*dm**2*dn** &
      &              3+1391616.0_cp*a**3*b**4*dm**2*dn**2+2751456.0_cp*a**3*b**4*dm**2*dn &
      &              -15497280.0_cp*a**3*b**4*dm**2-128960.0_cp*a**3*b**4*dn**6-349160.0_cp &
      &              *a**3*b**4*dn**5+7278984.0_cp*a**3*b**4*dn**4+13987304.0_cp*a**3*b**4* &
      &              dn**3-144087592.0_cp*a**3*b**4*dn**2-141705216.0_cp*a**3*b**4*dn+986827104.0_cp &
      &              *a**3*b**4+24576.0_cp*a**2*b**5*dm**2*dn**4+142848.0_cp*a**2*b**5*dm**2 &
      &              *dn**3-745472.0_cp*a**2*b**5*dm**2*dn**2-1510912.0_cp*a**2*b**5*dm**2* &
      &              dn+5394432.0_cp*a**2*b**5*dm**2+232448.0_cp*a**2*b**5*dn**6+640640.0_cp &
      &              *a**2*b**5*dn**5-9423488.0_cp*a**2*b**5*dn**4-18437760.0_cp*a**2*b**5* &
      &              dn**3+137034880.0_cp*a**2*b**5*dn**2+137763840.0_cp*a**2*b**5*dn-710613504.0_cp &
      &              *a**2*b**5-186368.0_cp*a*b**6*dn**6-503552.0_cp*a*b**6*dn**5+5602048.0_cp &
      &              *a*b**6*dn**4+11220736.0_cp*a*b**6*dn**3-60329728.0_cp*a*b**6*dn**2 &
      &              -62748672.0_cp*a*b**6*dn+232713216.0_cp*a*b**6+57344.0_cp*b**7*dn** &
      &              6+164864.0_cp*b**7*dn**5-1238016.0_cp*b**7*dn**4-2622464.0_cp*b**7*dn**3 &
      &              +9425920.0_cp*b**7*dn**2+10248192.0_cp*b**7*dn-26357760.0_cp*b**7)/ &
      &              (1024.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku+1)=a**4*(231.0_cp*a**8*dm**4*dn**2-1239.0_cp*a**8*dm**4+42.0_cp &
      &              *a**8*dm**2*dn**4-1456.0_cp*a**8*dm**2*dn**2+4822.0_cp*a**8*dm**2+7.0_cp &
      &              *a**8*dn**6+87.0_cp*a**8*dn**4+506.0_cp*a**8*dn**2+16776.0_cp*a**8+2016.0_cp &
      &              *a**7*b*dm**4*dn**2-11424.0_cp*a**7*b*dm**4+448.0_cp*a**7*b*dm**2*dn**4 &
      &              -14272.0_cp*a**7*b*dm**2*dn**2+44400.0_cp*a**7*b*dm**2+96.0_cp*a**7 &
      &              *b*dn**6+1040.0_cp*a**7*b*dn**4+12256.0_cp*a**7*b*dn**2+185472.0_cp &
      &              *a**7*b+7392.0_cp*a**6*b**2*dm**4*dn**2-48288.0_cp*a**6*b**2*dm**4+2104.0_cp &
      &              *a**6*b**2*dm**2*dn**4-63352.0_cp*a**6*b**2*dm**2*dn**2+171384.0_cp &
      &              *a**6*b**2*dm**2+648.0_cp*a**6*b**2*dn**6+8564.0_cp*a**6*b**2*dn**4 &
      &              +373036.0_cp*a**6*b**2*dn**2-20113920.0_cp*a**6*b**2+14080.0_cp*a** &
      &              5*b**3*dm**4*dn**2-121600.0_cp*a**5*b**3*dm**4+5504.0_cp*a**5*b**3*dm**2 &
      &              *dn**4-178496.0_cp*a**5*b**3*dm**2*dn**2-299520.0_cp*a**5*b**3*dm** &
      &              2+3200.0_cp*a**5*b**3*dn**6+163456.0_cp*a**5*b**3*dn**4-13595040.0_cp &
      &              *a**5*b**3*dn**2+207844704.0_cp*a**5*b**3+12672.0_cp*a**4*b**4*dm** &
      &              4*dn**2-183168.0_cp*a**4*b**4*dm**4+7296.0_cp*a**4*b**4*dm**2*dn**4 &
      &              -596784.0_cp*a**4*b**4*dm**2*dn**2+11736240.0_cp*a**4*b**4*dm**2+23320.0_cp &
      &              *a**4*b**4*dn**6-2479936.0_cp*a**4*b**4*dn**4+75514776.0_cp*a**4*b**4* &
      &              dn**2-699046128.0_cp*a**4*b**4-12288.0_cp*a**3*b**5*dm**2*dn**4+1605120.0_cp &
      &              *a**3*b**5*dm**2*dn**2-16937472.0_cp*a**3*b**5*dm**2-128512.0_cp*a**3* &
      &              b**5*dn**6+7979264.0_cp*a**3*b**5*dn**4-159716608.0_cp*a**3*b**5*dn**2 &
      &              +1061991936.0_cp*a**3*b**5+43008.0_cp*a**2*b**6*dm**2*dn**4-1050112.0_cp &
      &              *a**2*b**6*dm**2*dn**2+6192640.0_cp*a**2*b**6*dm**2+254464.0_cp*a** &
      &              2*b**6*dn**6-10773504.0_cp*a**2*b**6*dn**4+158520320.0_cp*a**2*b**6 &
      &              *dn**2-781097472.0_cp*a**2*b**6-200704.0_cp*a*b**7*dn**6+6707200.0_cp &
      &              *a*b**7*dn**4-73721856.0_cp*a*b**7*dn**2+266526720.0_cp*a*b**7+71680.0_cp &
      &              *b**8*dn**6-1679360.0_cp*b**8*dn**4+12851200.0_cp*b**8*dn**2-31887360.0_cp &
      &              *b**8)/(1024.0_cp*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn &
      &              +2.0_cp)*(dn+3.0_cp))
      stencil(ku+2)=-a**4*b*(168.0_cp*a**7*dm**4*dn**2-840.0_cp*a**7*dm** &
      &              4*dn-2352.0_cp*a**7*dm**4+28.0_cp*a**7*dm**2*dn**4-21.0_cp*a**7*dm**2* &
      &              dn**3-1574.0_cp*a**7*dm**2*dn**2+5069.0_cp*a**7*dm**2*dn+8058.0_cp* &
      &              a**7*dm**2+12.0_cp*a**7*dn**6-45.0_cp*a**7*dn**5+321.0_cp*a**7*dn** &
      &              4-1515.0_cp*a**7*dn**3+8031.0_cp*a**7*dn**2-16200.0_cp*a**7*dn+41652.0_cp &
      &              *a**7+1344.0_cp*a**6*b*dm**4*dn**2-6720.0_cp*a**6*b*dm**4*dn-20736.0_cp &
      &              *a**6*b*dm**4+256.0_cp*a**6*b*dm**2*dn**4+24.0_cp*a**6*b*dm**2*dn** &
      &              3-16720.0_cp*a**6*b*dm**2*dn**2+47544.0_cp*a**6*b*dm**2*dn+53184.0_cp &
      &              *a**6*b*dm**2+192.0_cp*a**6*b*dn**6-696.0_cp*a**6*b*dn**5+6728.0_cp &
      &              *a**6*b*dn**4-19056.0_cp*a**6*b*dn**3+374440.0_cp*a**6*b*dn**2-222120.0_cp &
      &              *a**6*b*dn-20625696.0_cp*a**6*b+4224.0_cp*a**5*b**2*dm**4*dn**2-21120.0_cp &
      &              *a**5*b**2*dm**4*dn-78336.0_cp*a**5*b**2*dm**4+832.0_cp*a**5*b**2*dm**2 &
      &              *dn**4+2292.0_cp*a**5*b**2*dm**2*dn**3-90456.0_cp*a**5*b**2*dm**2*dn**2 &
      &              +208372.0_cp*a**5*b**2*dm**2*dn-528048.0_cp*a**5*b**2*dm**2+1744.0_cp &
      &              *a**5*b**2*dn**6-5968.0_cp*a**5*b**2*dn**5+167920.0_cp*a**5*b**2*dn**4 &
      &              -300604.0_cp*a**5*b**2*dn**3-13363528.0_cp*a**5*b**2*dn**2+12709284.0_cp &
      &              *a**5*b**2*dn+200000016.0_cp*a**5*b**2+5632.0_cp*a**4*b**3*dm**4*dn**2 &
      &              -28160.0_cp*a**4*b**3*dm**4*dn-150528.0_cp*a**4*b**3*dm**4-256.0_cp &
      &              *a**4*b**3*dm**2*dn**4+20496.0_cp*a**4*b**3*dm**2*dn**3-534496.0_cp &
      &              *a**4*b**3*dm**2*dn**2+900016.0_cp*a**4*b**3*dm**2*dn+11231520.0_cp &
      &              *a**4*b**3*dm**2+20480.0_cp*a**4*b**3*dn**6-63200.0_cp*a**4*b**3*dn**5 &
      &              -2372064.0_cp*a**4*b**3*dn**4+4572144.0_cp*a**4*b**3*dn**3+70316288.0_cp &
      &              *a**4*b**3*dn**2-68285136.0_cp*a**4*b**3*dn-666389088.0_cp*a**4*b** &
      &              3-21504.0_cp*a**3*b**4*dm**2*dn**4+153504.0_cp*a**3*b**4*dm**2*dn** &
      &              3+1391616.0_cp*a**3*b**4*dm**2*dn**2-2751456.0_cp*a**3*b**4*dm**2*dn &
      &              -15497280.0_cp*a**3*b**4*dm**2-128960.0_cp*a**3*b**4*dn**6+349160.0_cp &
      &              *a**3*b**4*dn**5+7278984.0_cp*a**3*b**4*dn**4-13987304.0_cp*a**3*b**4* &
      &              dn**3-144087592.0_cp*a**3*b**4*dn**2+141705216.0_cp*a**3*b**4*dn+986827104.0_cp &
      &              *a**3*b**4+24576.0_cp*a**2*b**5*dm**2*dn**4-142848.0_cp*a**2*b**5*dm**2 &
      &              *dn**3-745472.0_cp*a**2*b**5*dm**2*dn**2+1510912.0_cp*a**2*b**5*dm**2* &
      &              dn+5394432.0_cp*a**2*b**5*dm**2+232448.0_cp*a**2*b**5*dn**6-640640.0_cp &
      &              *a**2*b**5*dn**5-9423488.0_cp*a**2*b**5*dn**4+18437760.0_cp*a**2*b**5* &
      &              dn**3+137034880.0_cp*a**2*b**5*dn**2-137763840.0_cp*a**2*b**5*dn-710613504.0_cp &
      &              *a**2*b**5-186368.0_cp*a*b**6*dn**6+503552.0_cp*a*b**6*dn**5+5602048.0_cp &
      &              *a*b**6*dn**4-11220736.0_cp*a*b**6*dn**3-60329728.0_cp*a*b**6*dn**2 &
      &              +62748672.0_cp*a*b**6*dn+232713216.0_cp*a*b**6+57344.0_cp*b**7*dn** &
      &              6-164864.0_cp*b**7*dn**5-1238016.0_cp*b**7*dn**4+2622464.0_cp*b**7*dn**3 &
      &              +9425920.0_cp*b**7*dn**2-10248192.0_cp*b**7*dn-26357760.0_cp*b**7)/ &
      &              (1024.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku+3)=-a**4*(198.0_cp*a**8*dm**4*dn**2+594.0_cp*a**8*dm**4*dn &
      &              +216.0_cp*a**8*dm**4+24.0_cp*a**8*dm**2*dn**4+117.0_cp*a**8*dm**2*dn**3 &
      &              -1268.0_cp*a**8*dm**2*dn**2-2575.0_cp*a**8*dm**2*dn+294.0_cp*a**8*dm**2 &
      &              +2.0_cp*a**8*dn**6-25.0_cp*a**8*dn**5+34.0_cp*a**8*dn**4+261.0_cp*a**8 &
      &              *dn**3-2060.0_cp*a**8*dn**2+2412.0_cp*a**8*dn-18000.0_cp*a**8+1680.0_cp &
      &              *a**7*b*dm**4*dn**2+5376.0_cp*a**7*b*dm**4*dn+2352.0_cp*a**7*b*dm** &
      &              4+224.0_cp*a**7*b*dm**2*dn**4+1428.0_cp*a**7*b*dm**2*dn**3-11824.0_cp &
      &              *a**7*b*dm**2*dn**2-25300.0_cp*a**7*b*dm**2*dn+4896.0_cp*a**7*b*dm**2+ &
      &              16.0_cp*a**7*b*dn**6-260.0_cp*a**7*b*dn**5-80.0_cp*a**7*b*dn**4+5844.0_cp &
      &              *a**7*b*dn**3-33168.0_cp*a**7*b*dn**2+64368.0_cp*a**7*b*dn-235584.0_cp &
      &              *a**7*b+5712.0_cp*a**6*b**2*dm**4*dn**2+21504.0_cp*a**6*b**2*dm**4*dn &
      &              +13680.0_cp*a**6*b**2*dm**4+836.0_cp*a**6*b**2*dm**2*dn**4+7638.0_cp &
      &              *a**6*b**2*dm**2*dn**3-43348.0_cp*a**6*b**2*dm**2*dn**2-121894.0_cp &
      &              *a**6*b**2*dm**2*dn+46632.0_cp*a**6*b**2*dm**2-28.0_cp*a**6*b**2*dn**6 &
      &              -748.0_cp*a**6*b**2*dn**5-7682.0_cp*a**6*b**2*dn**4+59220.0_cp*a**6 &
      &              *b**2*dn**3-561890.0_cp*a**6*b**2*dn**2+678744.0_cp*a**6*b**2*dn+19564056.0_cp &
      &              *a**6*b**2+9216.0_cp*a**5*b**3*dm**4*dn**2+47616.0_cp*a**5*b**3*dm**4* &
      &              dn+50688.0_cp*a**5*b**3*dm**4+1536.0_cp*a**5*b**3*dm**2*dn**4+21360.0_cp &
      &              *a**5*b**3*dm**2*dn**3-49472.0_cp*a**5*b**3*dm**2*dn**2-418768.0_cp &
      &              *a**5*b**3*dm**2*dn+917856.0_cp*a**5*b**3*dm**2-1280.0_cp*a**5*b**3 &
      &              *dn**6+4768.0_cp*a**5*b**3*dn**5-192064.0_cp*a**5*b**3*dn**4+729008.0_cp &
      &              *a**5*b**3*dn**3+12054336.0_cp*a**5*b**3*dn**2-24059376.0_cp*a**5*b**3 &
      &              *dn-182951712.0_cp*a**5*b**3+5632.0_cp*a**4*b**4*dm**4*dn**2+53504.0_cp &
      &              *a**4*b**4*dm**4*dn+111360.0_cp*a**4*b**4*dm**4+2560.0_cp*a**4*b**4 &
      &              *dm**2*dn**4+15120.0_cp*a**4*b**4*dm**2*dn**3+415872.0_cp*a**4*b**4 &
      &              *dm**2*dn**2-1855312.0_cp*a**4*b**4*dm**2*dn-9724992.0_cp*a**4*b**4 &
      &              *dm**2-20512.0_cp*a**4*b**4*dn**6+115640.0_cp*a**4*b**4*dn**5+2031088.0_cp &
      &              *a**4*b**4*dn**4-8266000.0_cp*a**4*b**4*dn**3-56665888.0_cp*a**4*b**4* &
      &              dn**2+121347096.0_cp*a**4*b**4*dn+567446544.0_cp*a**4*b**4+23552.0_cp &
      &              *a**3*b**5*dm**2*dn**4-244224.0_cp*a**3*b**5*dm**2*dn**3-776576.0_cp &
      &              *a**3*b**5*dm**2*dn**2+4463616.0_cp*a**3*b**5*dm**2*dn+11878272.0_cp &
      &              *a**3*b**5*dm**2+116608.0_cp*a**3*b**5*dn**6-644000.0_cp*a**3*b**5*dn**5 &
      &              -5371520.0_cp*a**3*b**5*dn**4+23629056.0_cp*a**3*b**5*dn**3+102374976.0_cp &
      &              *a**3*b**5*dn**2-235108704.0_cp*a**3*b**5*dn-795122496.0_cp*a**3*b**5+ &
      &              7680.0_cp*a**2*b**6*dm**2*dn**4+157824.0_cp*a**2*b**6*dm**2*dn**3+131968.0_cp &
      &              *a**2*b**6*dm**2*dn**2-1983616.0_cp*a**2*b**6*dm**2*dn-3499392.0_cp &
      &              *a**2*b**6*dm**2-186496.0_cp*a**2*b**6*dn**6+1013600.0_cp*a**2*b**6 &
      &              *dn**5+6138304.0_cp*a**2*b**6*dn**4-28473408.0_cp*a**2*b**6*dn**3-85512320.0_cp &
      &              *a**2*b**6*dn**2+213557472.0_cp*a**2*b**6*dn+526559040.0_cp*a**2*b**6+ &
      &              139264.0_cp*a*b**7*dn**6-763904.0_cp*a*b**7*dn**5-2981888.0_cp*a*b**7* &
      &              dn**4+15616000.0_cp*a*b**7*dn**3+29966336.0_cp*a*b**7*dn**2-86900736.0_cp &
      &              *a*b**7*dn-154386432.0_cp*a*b**7-28672.0_cp*b**8*dn**6+164864.0_cp* &
      &              b**8*dn**5+460800.0_cp*b**8*dn**4-2754560.0_cp*b**8*dn**3-3350528.0_cp &
      &              *b**8*dn**2+11682816.0_cp*b**8*dn+14469120.0_cp*b**8)/(1024.0_cp*dn &
      &              *(dn-3.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+4)=a**4*b*(360.0_cp*a**7*dm**4*dn**2+936.0_cp*a**7*dm**4 &
      &              *dn+48.0_cp*a**7*dm**4-4.0_cp*a**7*dm**2*dn**4+621.0_cp*a**7*dm**2*dn**3 &
      &              -3422.0_cp*a**7*dm**2*dn**2-4737.0_cp*a**7*dm**2*dn+6078.0_cp*a**7*dm**2 &
      &              -4.0_cp*a**7*dn**6-15.0_cp*a**7*dn**5-367.0_cp*a**7*dn**4+3827.0_cp &
      &              *a**7*dn**3-19849.0_cp*a**7*dn**2+73036.0_cp*a**7*dn-124404.0_cp*a**7+ &
      &              2688.0_cp*a**6*b*dm**4*dn**2+8064.0_cp*a**6*b*dm**4*dn+1920.0_cp*a**6* &
      &              b*dm**4-256.0_cp*a**6*b*dm**2*dn**4+7488.0_cp*a**6*b*dm**2*dn**3-30080.0_cp &
      &              *a**6*b*dm**2*dn**2-54864.0_cp*a**6*b*dm**2*dn+88560.0_cp*a**6*b*dm**2 &
      &              -128.0_cp*a**6*b*dn**6+528.0_cp*a**6*b*dn**5-10960.0_cp*a**6*b*dn** &
      &              4+75184.0_cp*a**6*b*dn**3-636336.0_cp*a**6*b*dn**2+1163312.0_cp*a** &
      &              6*b*dn+19410864.0_cp*a**6*b+7424.0_cp*a**5*b**2*dm**4*dn**2+28416.0_cp &
      &              *a**5*b**2*dm**4*dn+17152.0_cp*a**5*b**2*dm**4-2432.0_cp*a**5*b**2*dm**2 &
      &              *dn**4+38184.0_cp*a**5*b**2*dm**2*dn**3-86160.0_cp*a**5*b**2*dm**2*dn**2 &
      &              -359184.0_cp*a**5*b**2*dm**2*dn+1237928.0_cp*a**5*b**2*dm**2-1760.0_cp &
      &              *a**5*b**2*dn**6+12432.0_cp*a**5*b**2*dn**5-236176.0_cp*a**5*b**2*dn**4 &
      &              +1204248.0_cp*a**5*b**2*dn**3+10182144.0_cp*a**5*b**2*dn**2-33471600.0_cp &
      &              *a**5*b**2*dn-152092728.0_cp*a**5*b**2+7680.0_cp*a**4*b**3*dm**4*dn**2 &
      &              +44544.0_cp*a**4*b**3*dm**4*dn+64512.0_cp*a**4*b**3*dm**4-9472.0_cp &
      &              *a**4*b**3*dm**2*dn**4+88176.0_cp*a**4*b**3*dm**2*dn**3+299680.0_cp &
      &              *a**4*b**3*dm**2*dn**2-2651376.0_cp*a**4*b**3*dm**2*dn-7503264.0_cp &
      &              *a**4*b**3*dm**2-23552.0_cp*a**4*b**3*dn**6+198240.0_cp*a**4*b**3*dn**5 &
      &              +1404960.0_cp*a**4*b**3*dn**4-10470896.0_cp*a**4*b**3*dn**3-36396928.0_cp &
      &              *a**4*b**3*dn**2+145884560.0_cp*a**4*b**3*dn+434051808.0_cp*a**4*b**3- &
      &              1024.0_cp*a**3*b**4*dm**2*dn**4-204960.0_cp*a**3*b**4*dm**2*dn**3+34816.0_cp &
      &              *a**3*b**4*dm**2*dn**2+4366176.0_cp*a**3*b**4*dm**2*dn+7334208.0_cp &
      &              *a**3*b**4*dm**2+89408.0_cp*a**3*b**4*dn**6-738024.0_cp*a**3*b**4*dn**5 &
      &              -2985624.0_cp*a**3*b**4*dn**4+25985576.0_cp*a**3*b**4*dn**3+52271416.0_cp &
      &              *a**3*b**4*dn**2-256232576.0_cp*a**3*b**4*dn-540212640.0_cp*a**3*b**4+ &
      &              20480.0_cp*a**2*b**5*dm**2*dn**4+82176.0_cp*a**2*b**5*dm**2*dn**3-253952.0_cp &
      &              *a**2*b**5*dm**2*dn**2-1499904.0_cp*a**2*b**5*dm**2*dn-1654272.0_cp &
      &              *a**2*b**5*dm**2-134656.0_cp*a**2*b**5*dn**6+1095360.0_cp*a**2*b**5 &
      &              *dn**5+2403904.0_cp*a**2*b**5*dn**4-27233472.0_cp*a**2*b**5*dn**3-29977920.0_cp &
      &              *a**2*b**5*dn**2+203132928.0_cp*a**2*b**5*dn+313516800.0_cp*a**2*b**5+ &
      &              72704.0_cp*a*b**6*dn**6-606336.0_cp*a*b**6*dn**5-679808.0_cp*a*b**6 &
      &              *dn**4+11779200.0_cp*a*b**6*dn**3+4811136.0_cp*a*b**6*dn**2-67651584.0_cp &
      &              *a*b**6*dn-73492992.0_cp*a*b**6-8192.0_cp*b**7*dn**6+70656.0_cp*b** &
      &              7*dn**5+56320.0_cp*b**7*dn**4-1274880.0_cp*b**7*dn**3-232448.0_cp*b**7 &
      &              *dn**2+6365184.0_cp*b**7*dn+5345280.0_cp*b**7)/(1024.0_cp*dn*(dn-2.0_cp &
      &              )*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+5)=a**4*(495.0_cp*a**8*dm**4*dn**3+990.0_cp*a**8*dm**4*dn**2 &
      &              -1755.0_cp*a**8*dm**4*dn-1890.0_cp*a**8*dm**4-30.0_cp*a**8*dm**2*dn**5 &
      &              +867.0_cp*a**8*dm**2*dn**4-5620.0_cp*a**8*dm**2*dn**3-6327.0_cp*a** &
      &              8*dm**2*dn**2+30790.0_cp*a**8*dm**2*dn-1584.0_cp*a**8*dm**2-17.0_cp &
      &              *a**8*dn**7+223.0_cp*a**8*dn**6-1365.0_cp*a**8*dn**5+5749.0_cp*a**8 &
      &              *dn**4-13018.0_cp*a**8*dn**3-12620.0_cp*a**8*dn**2+180120.0_cp*a**8 &
      &              *dn-431136.0_cp*a**8+3840.0_cp*a**7*b*dm**4*dn**3+9216.0_cp*a**7*b*dm**4 &
      &              *dn**2-11520.0_cp*a**7*b*dm**4*dn-17664.0_cp*a**7*b*dm**4-512.0_cp* &
      &              a**7*b*dm**2*dn**5+9920.0_cp*a**7*b*dm**2*dn**4-48896.0_cp*a**7*b*dm**2 &
      &              *dn**3-91136.0_cp*a**7*b*dm**2*dn**2+320896.0_cp*a**7*b*dm**2*dn+15744.0_cp &
      &              *a**7*b*dm**2-256.0_cp*a**7*b*dn**7+3584.0_cp*a**7*b*dn**6-20864.0_cp &
      &              *a**7*b*dn**5+68544.0_cp*a**7*b*dn**4-36864.0_cp*a**7*b*dn**3-902912.0_cp &
      &              *a**7*b*dn**2+4399104.0_cp*a**7*b*dn-7575552.0_cp*a**7*b+9984.0_cp* &
      &              a**6*b**2*dm**4*dn**3+36864.0_cp*a**6*b**2*dm**4*dn**2-5376.0_cp*a**6* &
      &              b**2*dm**4*dn-66816.0_cp*a**6*b**2*dm**4-2624.0_cp*a**6*b**2*dm**2*dn**5 &
      &              +35072.0_cp*a**6*b**2*dm**2*dn**4-88512.0_cp*a**6*b**2*dm**2*dn**3-656672.0_cp &
      &              *a**6*b**2*dm**2*dn**2+1111232.0_cp*a**6*b**2*dm**2*dn+1022016.0_cp &
      &              *a**6*b**2*dm**2-1472.0_cp*a**6*b**2*dn**7+21952.0_cp*a**6*b**2*dn**6- &
      &              99936.0_cp*a**6*b**2*dn**5+7872.0_cp*a**6*b**2*dn**4+3643232.0_cp*a**6 &
      &              *b**2*dn**3-21614752.0_cp*a**6*b**2*dn**2-30415104.0_cp*a**6*b**2*dn &
      &              +267545088.0_cp*a**6*b**2+6144.0_cp*a**5*b**3*dm**4*dn**3+67584.0_cp &
      &              *a**5*b**3*dm**4*dn**2+129024.0_cp*a**5*b**3*dm**4*dn-55296.0_cp*a**5* &
      &              b**3*dm**4-1024.0_cp*a**5*b**3*dm**2*dn**5-33344.0_cp*a**5*b**3*dm**2* &
      &              dn**4+509184.0_cp*a**5*b**3*dm**2*dn**3-2094016.0_cp*a**5*b**3*dm** &
      &              2*dn**2-4554752.0_cp*a**5*b**3*dm**2*dn+21381888.0_cp*a**5*b**3*dm**2- &
      &              1024.0_cp*a**5*b**3*dn**7+19328.0_cp*a**5*b**3*dn**6+722944.0_cp*a**5* &
      &              b**3*dn**5-9485120.0_cp*a**5*b**3*dn**4-2766336.0_cp*a**5*b**3*dn** &
      &              3+262940352.0_cp*a**5*b**3*dn**2-127774464.0_cp*a**5*b**3*dn-1896837120.0_cp &
      &              *a**5*b**3-9216.0_cp*a**4*b**4*dm**4*dn**3+21504.0_cp*a**4*b**4*dm**4* &
      &              dn**2+297984.0_cp*a**4*b**4*dm**4*dn+451584.0_cp*a**4*b**4*dm**4+37888.0_cp &
      &              *a**4*b**4*dm**2*dn**5-570304.0_cp*a**4*b**4*dm**2*dn**4+1406336.0_cp &
      &              *a**4*b**4*dm**2*dn**3+11473216.0_cp*a**4*b**4*dm**2*dn**2-23014272.0_cp &
      &              *a**4*b**4*dm**2*dn-78929280.0_cp*a**4*b**4*dm**2+80320.0_cp*a**4*b**4 &
      &              *dn**7-1229600.0_cp*a**4*b**4*dn**6+748032.0_cp*a**4*b**4*dn**5+54200384.0_cp &
      &              *a**4*b**4*dn**4-108542272.0_cp*a**4*b**4*dn**3-810423072.0_cp*a**4 &
      &              *b**4*dn**2+1108388736.0_cp*a**4*b**4*dn+4551911424.0_cp*a**4*b**4+90112.0_cp &
      &              *a**3*b**5*dm**2*dn**5-16384.0_cp*a**3*b**5*dm**2*dn**4-3644416.0_cp &
      &              *a**3*b**5*dm**2*dn**3-2920448.0_cp*a**3*b**5*dm**2*dn**2+34563072.0_cp &
      &              *a**3*b**5*dm**2*dn+54798336.0_cp*a**3*b**5*dm**2-259072.0_cp*a**3* &
      &              b**5*dn**7+3849728.0_cp*a**3*b**5*dn**6-8172032.0_cp*a**3*b**5*dn** &
      &              5-97853440.0_cp*a**3*b**5*dn**4+298031616.0_cp*a**3*b**5*dn**3+1015069184.0_cp &
      &              *a**3*b**5*dn**2-2129931264.0_cp*a**3*b**5*dn-4911280128.0_cp*a**3* &
      &              b**5-53248.0_cp*a**2*b**6*dm**2*dn**5+133120.0_cp*a**2*b**6*dm**2*dn**4 &
      &              +1313792.0_cp*a**2*b**6*dm**2*dn**3-937984.0_cp*a**2*b**6*dm**2*dn**2- &
      &              9763840.0_cp*a**2*b**6*dm**2*dn-9099264.0_cp*a**2*b**6*dm**2+302080.0_cp &
      &              *a**2*b**6*dn**7-4516352.0_cp*a**2*b**6*dn**6+13649920.0_cp*a**2*b**6* &
      &              dn**5+70317056.0_cp*a**2*b**6*dn**4-307349504.0_cp*a**2*b**6*dn**3-461376000.0_cp &
      &              *a**2*b**6*dn**2+1597879296.0_cp*a**2*b**6*dn+2221903872.0_cp*a**2* &
      &              b**6-90112.0_cp*a*b**7*dn**7+1372160.0_cp*a*b**7*dn**6-4550656.0_cp &
      &              *a*b**7*dn**5-17854464.0_cp*a*b**7*dn**4+91000832.0_cp*a*b**7*dn**3 &
      &              +82739200.0_cp*a*b**7*dn**2-417644544.0_cp*a*b**7*dn-397541376.0_cp &
      &              *a*b**7+4096.0_cp*b**8*dn**7-63488.0_cp*b**8*dn**6+212992.0_cp*b**8 &
      &              *dn**5+839680.0_cp*b**8*dn**4-4222976.0_cp*b**8*dn**3-4093952.0_cp* &
      &              b**8*dn**2+18751488.0_cp*b**8*dn+18063360.0_cp*b**8)/(4096.0_cp*dn* &
      &              (dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp &
      &              ))
      stencil(ku+6)=-a**5*b*(600.0_cp*a**6*dm**4*dn**3+1008.0_cp*a**6*dm**4* &
      &              dn**2-3000.0_cp*a**6*dm**4*dn-2832.0_cp*a**6*dm**4-220.0_cp*a**6*dm**2 &
      &              *dn**5+3329.0_cp*a**6*dm**2*dn**4-14155.0_cp*a**6*dm**2*dn**3-19583.0_cp &
      &              *a**6*dm**2*dn**2+107135.0_cp*a**6*dm**2*dn-22362.0_cp*a**6*dm**2-60.0_cp &
      &              *a**6*dn**7+1125.0_cp*a**6*dn**6-8334.0_cp*a**6*dn**5+31144.0_cp*a**6* &
      &              dn**4-13926.0_cp*a**6*dn**3-454177.0_cp*a**6*dn**2+2205888.0_cp*a** &
      &              6*dn-3533436.0_cp*a**6+3840.0_cp*a**5*b*dm**4*dn**3+9216.0_cp*a**5* &
      &              b*dm**4*dn**2-15360.0_cp*a**5*b*dm**4*dn-25344.0_cp*a**5*b*dm**4-2560.0_cp &
      &              *a**5*b*dm**2*dn**5+33344.0_cp*a**5*b*dm**2*dn**4-101568.0_cp*a**5* &
      &              b*dm**2*dn**3-334112.0_cp*a**5*b*dm**2*dn**2+1076992.0_cp*a**5*b*dm**2 &
      &              *dn+172704.0_cp*a**5*b*dm**2-768.0_cp*a**5*b*dn**7+14880.0_cp*a**5* &
      &              b*dn**6-98816.0_cp*a**5*b*dn**5+203264.0_cp*a**5*b*dn**4+1911872.0_cp &
      &              *a**5*b*dn**3-15791360.0_cp*a**5*b*dn**2+3304320.0_cp*a**5*b*dn+150651936.0_cp &
      &              *a**5*b+7680.0_cp*a**4*b**2*dm**4*dn**3+30720.0_cp*a**4*b**2*dm**4*dn**2 &
      &              -69120.0_cp*a**4*b**2*dm**4-10496.0_cp*a**4*b**2*dm**2*dn**5+105040.0_cp &
      &              *a**4*b**2*dm**2*dn**4-54896.0_cp*a**4*b**2*dm**2*dn**3-2164960.0_cp &
      &              *a**4*b**2*dm**2*dn**2+1736080.0_cp*a**4*b**2*dm**2*dn+12151920.0_cp &
      &              *a**4*b**2*dm**2-2880.0_cp*a**4*b**2*dn**7+59040.0_cp*a**4*b**2*dn**6- &
      &              6784.0_cp*a**4*b**2*dn**5-4894960.0_cp*a**4*b**2*dn**4+13695760.0_cp &
      &              *a**4*b**2*dn**3+99050560.0_cp*a**4*b**2*dn**2-202721136.0_cp*a**4* &
      &              b**2*dn-784328400.0_cp*a**4*b**2+3072.0_cp*a**3*b**3*dm**4*dn**3+30720.0_cp &
      &              *a**3*b**3*dm**4*dn**2+95232.0_cp*a**3*b**3*dm**4*dn+92160.0_cp*a** &
      &              3*b**3*dm**4-8704.0_cp*a**3*b**3*dm**2*dn**5-47200.0_cp*a**3*b**3*dm**2 &
      &              *dn**4+963616.0_cp*a**3*b**3*dm**2*dn**3+1027360.0_cp*a**3*b**3*dm**2* &
      &              dn**2-14403360.0_cp*a**3*b**3*dm**2*dn-24730560.0_cp*a**3*b**3*dm** &
      &              2+30720.0_cp*a**3*b**3*dn**7-576960.0_cp*a**3*b**3*dn**6+1892864.0_cp &
      &              *a**3*b**3*dn**5+18138400.0_cp*a**3*b**3*dn**4-83616224.0_cp*a**3*b**3 &
      &              *dn**3-227908000.0_cp*a**3*b**3*dn**2+781839072.0_cp*a**3*b**3*dn+1617595200.0_cp &
      &              *a**3*b**3+38912.0_cp*a**2*b**4*dm**2*dn**5-181696.0_cp*a**2*b**4*dm**2 &
      &              *dn**4-1129280.0_cp*a**2*b**4*dm**2*dn**3+2291776.0_cp*a**2*b**4*dm**2 &
      &              *dn**2+12894528.0_cp*a**2*b**4*dm**2*dn+11740032.0_cp*a**2*b**4*dm**2- &
      &              81792.0_cp*a**2*b**4*dn**7+1520784.0_cp*a**2*b**4*dn**6-6862016.0_cp &
      &              *a**2*b**4*dn**5-20296032.0_cp*a**2*b**4*dn**4+158494592.0_cp*a**2* &
      &              b**4*dn**3+136013712.0_cp*a**2*b**4*dn**2-1022954304.0_cp*a**2*b**4 &
      &              *dn-1324654272.0_cp*a**2*b**4-8192.0_cp*a*b**5*dm**2*dn**5+39424.0_cp &
      &              *a*b**5*dm**2*dn**4+204288.0_cp*a*b**5*dm**2*dn**3-497152.0_cp*a*b**5* &
      &              dm**2*dn**2-2076160.0_cp*a*b**5*dm**2*dn-1422336.0_cp*a*b**5*dm**2+52224.0_cp &
      &              *a*b**5*dn**7-984192.0_cp*a*b**5*dn**6+4791296.0_cp*a*b**5*dn**5+9059584.0_cp &
      &              *a*b**5*dn**4-95214080.0_cp*a*b**5*dn**3-19867264.0_cp*a*b**5*dn**2 &
      &              +525623808.0_cp*a*b**5*dn+447045120.0_cp*a*b**5-6144.0_cp*b**6*dn** &
      &              7+117504.0_cp*b**6*dn**6-578560.0_cp*b**6*dn**5-1094144.0_cp*b**6*dn**4 &
      &              +11474944.0_cp*b**6*dn**3+2829056.0_cp*b**6*dn**2-62131200.0_cp*b** &
      &              6*dn-53093376.0_cp*b**6)/(2048.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp &
      &              )*(dn+1.0_cp)*(dn+2.0_cp)*(dn+3.0_cp))
      stencil(ku+7)=-a**6*(110.0_cp*a**6*dm**4*dn**3+330.0_cp*a**6*dm**4*dn**2 &
      &              -140.0_cp*a**6*dm**4*dn-480.0_cp*a**6*dm**4-40.0_cp*a**6*dm**2*dn** &
      &              5+579.0_cp*a**6*dm**2*dn**4-1910.0_cp*a**6*dm**2*dn**3-7449.0_cp*a**6* &
      &              dm**2*dn**2+13620.0_cp*a**6*dm**2*dn+11292.0_cp*a**6*dm**2-6.0_cp*a**6 &
      &              *dn**7+141.0_cp*a**6*dn**6-1100.0_cp*a**6*dn**5+2263.0_cp*a**6*dn** &
      &              4+14066.0_cp*a**6*dn**3-80740.0_cp*a**6*dn**2+87480.0_cp*a**6*dn+279648.0_cp &
      &              *a**6+720.0_cp*a**5*b*dm**4*dn**3+2592.0_cp*a**5*b*dm**4*dn**2+240.0_cp &
      &              *a**5*b*dm**4*dn-3168.0_cp*a**5*b*dm**4-416.0_cp*a**5*b*dm**2*dn**5 &
      &              +5340.0_cp*a**5*b*dm**2*dn**4-11608.0_cp*a**5*b*dm**2*dn**3-82572.0_cp &
      &              *a**5*b*dm**2*dn**2+90408.0_cp*a**5*b*dm**2*dn+167328.0_cp*a**5*b*dm**2 &
      &              -48.0_cp*a**5*b*dn**7+1188.0_cp*a**5*b*dn**6-8104.0_cp*a**5*b*dn**5 &
      &              -6916.0_cp*a**5*b*dn**4+273544.0_cp*a**5*b*dn**3-798800.0_cp*a**5*b &
      &              *dn**2-1320480.0_cp*a**5*b*dn+8233344.0_cp*a**5*b+912.0_cp*a**4*b** &
      &              2*dm**4*dn**3+6048.0_cp*a**4*b**2*dm**4*dn**2+10032.0_cp*a**4*b**2*dm**4 &
      &              *dn+288.0_cp*a**4*b**2*dm**4-812.0_cp*a**4*b**2*dm**2*dn**5+6114.0_cp &
      &              *a**4*b**2*dm**2*dn**4+30080.0_cp*a**4*b**2*dm**2*dn**3-236010.0_cp &
      &              *a**4*b**2*dm**2*dn**2-433308.0_cp*a**4*b**2*dm**2*dn+943776.0_cp*a**4 &
      &              *b**2*dm**2+180.0_cp*a**4*b**2*dn**7-3708.0_cp*a**4*b**2*dn**6+41966.0_cp &
      &              *a**4*b**2*dn**5-335816.0_cp*a**4*b**2*dn**4+228046.0_cp*a**4*b**2*dn**3 &
      &              +12288332.0_cp*a**4*b**2*dn**2-21173400.0_cp*a**4*b**2*dn-127462320.0_cp &
      &              *a**4*b**2-2048.0_cp*a**3*b**3*dm**4*dn**3-3072.0_cp*a**3*b**3*dm** &
      &              4*dn**2+23552.0_cp*a**3*b**3*dm**4*dn+43008.0_cp*a**3*b**3*dm**4+5120.0_cp &
      &              *a**3*b**3*dm**2*dn**5-72096.0_cp*a**3*b**3*dm**2*dn**4+205888.0_cp &
      &              *a**3*b**3*dm**2*dn**3+1190112.0_cp*a**3*b**3*dm**2*dn**2-3727872.0_cp &
      &              *a**3*b**3*dm**2*dn-9251712.0_cp*a**3*b**3*dm**2+2560.0_cp*a**3*b** &
      &              3*dn**7-60864.0_cp*a**3*b**3*dn**6+186112.0_cp*a**3*b**3*dn**5+3966432.0_cp &
      &              *a**3*b**3*dn**4-20849600.0_cp*a**3*b**3*dn**3-59435616.0_cp*a**3*b**3 &
      &              *dn**2+271767168.0_cp*a**3*b**3*dn+561195648.0_cp*a**3*b**3-3072.0_cp &
      &              *a**2*b**4*dm**4*dn**3-19968.0_cp*a**2*b**4*dm**4*dn**2-41472.0_cp* &
      &              a**2*b**4*dm**4*dn-27648.0_cp*a**2*b**4*dm**4+19456.0_cp*a**2*b**4*dm**2 &
      &              *dn**5-128352.0_cp*a**2*b**4*dm**2*dn**4-634048.0_cp*a**2*b**4*dm** &
      &              2*dn**3+2297568.0_cp*a**2*b**4*dm**2*dn**2+11086656.0_cp*a**2*b**4*dm**2 &
      &              *dn+10586880.0_cp*a**2*b**4*dm**2-24000.0_cp*a**2*b**4*dn**7+534960.0_cp &
      &              *a**2*b**4*dn**6-3149056.0_cp*a**2*b**4*dn**5-5856608.0_cp*a**2*b** &
      &              4*dn**4+80152576.0_cp*a**2*b**4*dn**3+36374384.0_cp*a**2*b**4*dn**2 &
      &              -652091328.0_cp*a**2*b**4*dn-852833088.0_cp*a**2*b**4-14336.0_cp*a* &
      &              b**5*dm**2*dn**5+101376.0_cp*a*b**5*dm**2*dn**4+359168.0_cp*a*b**5*dm**2 &
      &              *dn**3-1669632.0_cp*a*b**5*dm**2*dn**2-5558016.0_cp*a*b**5*dm**2*dn &
      &              -3644928.0_cp*a*b**5*dm**2+33024.0_cp*a*b**5*dn**7-742464.0_cp*a*b**5* &
      &              dn**6+4720256.0_cp*a*b**5*dn**5+2794496.0_cp*a*b**5*dn**4-95244416.0_cp &
      &              *a*b**5*dn**3+38876992.0_cp*a*b**5*dn**2+627962112.0_cp*a*b**5*dn+496541952.0_cp &
      &              *a*b**5+1024.0_cp*b**6*dm**2*dn**5-6912.0_cp*b**6*dm**2*dn**4-25856.0_cp &
      &              *b**6*dm**2*dn**3+109824.0_cp*b**6*dm**2*dn**2+375040.0_cp*b**6*dm**2* &
      &              dn+247296.0_cp*b**6*dm**2-7936.0_cp*b**6*dn**7+180672.0_cp*b**6*dn**6- &
      &              1160704.0_cp*b**6*dn**5-676224.0_cp*b**6*dn**4+23441408.0_cp*b**6*dn**3 &
      &              -8957760.0_cp*b**6*dn**2-152734464.0_cp*b**6*dn-121008384.0_cp*b**6 &
      &              )/(2048.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn &
      &              +2.0_cp)*(dn+3.0_cp))
      stencil(ku+8)=a**7*b*(280.0_cp*a**4*dm**4*dn**2+72.0_cp*a**4*dm**4*dn &
      &              -496.0_cp*a**4*dm**4-252.0_cp*a**4*dm**2*dn**4+4347.0_cp*a**4*dm**2 &
      &              *dn**3-22016.0_cp*a**4*dm**2*dn**2+12891.0_cp*a**4*dm**2*dn+37874.0_cp &
      &              *a**4*dm**2-12.0_cp*a**4*dn**6+435.0_cp*a**4*dn**5-4183.0_cp*a**4*dn**4 &
      &              -1983.0_cp*a**4*dn**3+242107.0_cp*a**4*dn**2-1424880.0_cp*a**4*dn+2723244.0_cp &
      &              *a**4+1344.0_cp*a**3*b*dm**4*dn**2+1728.0_cp*a**3*b*dm**4*dn-1920.0_cp &
      &              *a**3*b*dm**4-2048.0_cp*a**3*b*dm**2*dn**4+31736.0_cp*a**3*b*dm**2*dn**3 &
      &              -120224.0_cp*a**3*b*dm**2*dn**2-133928.0_cp*a**3*b*dm**2*dn+499696.0_cp &
      &              *a**3*b*dm**2+192.0_cp*a**3*b*dn**6-4776.0_cp*a**3*b*dn**5+70136.0_cp &
      &              *a**3*b*dn**4-786144.0_cp*a**3*b*dn**3+4417000.0_cp*a**3*b*dn**2-3647496.0_cp &
      &              *a**3*b*dn-32539440.0_cp*a**3*b+1152.0_cp*a**2*b**2*dm**4*dn**2+4992.0_cp &
      &              *a**2*b**2*dm**4*dn+5376.0_cp*a**2*b**2*dm**4-3264.0_cp*a**2*b**2*dm**2 &
      &              *dn**4+34236.0_cp*a**2*b**2*dm**2*dn**3+50224.0_cp*a**2*b**2*dm**2*dn**2 &
      &              -770580.0_cp*a**2*b**2*dm**2*dn-1415944.0_cp*a**2*b**2*dm**2+3408.0_cp &
      &              *a**2*b**2*dn**6-99264.0_cp*a**2*b**2*dn**5+923712.0_cp*a**2*b**2*dn**4 &
      &              -1686724.0_cp*a**2*b**2*dn**3-15091536.0_cp*a**2*b**2*dn**2+37507228.0_cp &
      &              *a**2*b**2*dn+103712856.0_cp*a**2*b**2-1024.0_cp*a*b**3*dm**4*dn**2 &
      &              -3072.0_cp*a*b**3*dm**4*dn-2048.0_cp*a*b**3*dm**4+9728.0_cp*a*b**3*dm**2 &
      &              *dn**4-119456.0_cp*a*b**3*dm**2*dn**3+141696.0_cp*a*b**3*dm**2*dn** &
      &              2+1407200.0_cp*a*b**3*dm**2*dn+1136320.0_cp*a*b**3*dm**2-12288.0_cp &
      &              *a*b**3*dn**6+358080.0_cp*a*b**3*dn**5-3584576.0_cp*a*b**3*dn**4+11823776.0_cp &
      &              *a*b**3*dn**3+15946688.0_cp*a*b**3*dn**2-100569632.0_cp*a*b**3*dn-100737600.0_cp &
      &              *a*b**3-2048.0_cp*b**4*dm**2*dn**4+24512.0_cp*b**4*dm**2*dn**3-25984.0_cp &
      &              *b**4*dm**2*dn**2-280256.0_cp*b**4*dm**2*dn-227712.0_cp*b**4*dm**2+5760.0_cp &
      &              *b**4*dn**6-169488.0_cp*b**4*dn**5+1709168.0_cp*b**4*dn**4-5656816.0_cp &
      &              *b**4*dn**3-7674800.0_cp*b**4*dn**2+47827840.0_cp*b**4*dn+47961408.0_cp &
      &              *b**4)/(2048.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              )*(dn+2.0_cp))
      stencil(ku+9)=a**8*(33.0_cp*a**4*dm**4*dn**2+33.0_cp*a**4*dm**4*dn-36.0_cp &
      &              *a**4*dm**4-26.0_cp*a**4*dm**2*dn**4+472.0_cp*a**4*dm**2*dn**3-2248.0_cp &
      &              *a**4*dm**2*dn**2-850.0_cp*a**4*dm**2*dn+4464.0_cp*a**4*dm**2+a**4*dn**6 &
      &              -25.0_cp*a**4*dn**5+332.0_cp*a**4*dn**4-3758.0_cp*a**4*dn**3+28368.0_cp &
      &              *a**4*dn**2-96408.0_cp*a**4*dn+46272.0_cp*a**4+160.0_cp*a**3*b*dm** &
      &              4*dn**2+288.0_cp*a**3*b*dm**4*dn-64.0_cp*a**3*b*dm**4-192.0_cp*a**3 &
      &              *b*dm**2*dn**4+3216.0_cp*a**3*b*dm**2*dn**3-11888.0_cp*a**3*b*dm**2 &
      &              *dn**2-23472.0_cp*a**3*b*dm**2*dn+29408.0_cp*a**3*b*dm**2+32.0_cp*a**3 &
      &              *b*dn**6-992.0_cp*a**3*b*dn**5+12368.0_cp*a**3*b*dn**4-82400.0_cp*a**3 &
      &              *b*dn**3+290240.0_cp*a**3*b*dn**2-96128.0_cp*a**3*b*dn-2244096.0_cp &
      &              *a**3*b-96.0_cp*a**2*b**2*dm**4*dn**2+288.0_cp*a**2*b**2*dm**4*dn+960.0_cp &
      &              *a**2*b**2*dm**4+232.0_cp*a**2*b**2*dm**2*dn**4-5048.0_cp*a**2*b**2 &
      &              *dm**2*dn**3+35200.0_cp*a**2*b**2*dm**2*dn**2-41656.0_cp*a**2*b**2*dm**2 &
      &              *dn-268208.0_cp*a**2*b**2*dm**2+88.0_cp*a**2*b**2*dn**6-3064.0_cp*a**2 &
      &              *b**2*dn**5+23460.0_cp*a**2*b**2*dn**4+188388.0_cp*a**2*b**2*dn**3-2734216.0_cp &
      &              *a**2*b**2*dn**2+4737664.0_cp*a**2*b**2*dn+21440256.0_cp*a**2*b**2-768.0_cp &
      &              *a*b**3*dm**4*dn**2-2304.0_cp*a*b**3*dm**4*dn-1536.0_cp*a*b**3*dm** &
      &              4+3200.0_cp*a*b**3*dm**2*dn**4-46192.0_cp*a*b**3*dm**2*dn**3+83984.0_cp &
      &              *a*b**3*dm**2*dn**2+623296.0_cp*a*b**3*dm**2*dn+489920.0_cp*a*b**3*dm**2 &
      &              -2688.0_cp*a*b**3*dn**6+88224.0_cp*a*b**3*dn**5-1006688.0_cp*a*b**3 &
      &              *dn**4+3991152.0_cp*a*b**3*dn**3+3347792.0_cp*a*b**3*dn**2-35899296.0_cp &
      &              *a*b**3*dn-34158336.0_cp*a*b**3+128.0_cp*b**4*dm**4*dn**2+384.0_cp* &
      &              b**4*dm**4*dn+256.0_cp*b**4*dm**4-1664.0_cp*b**4*dm**2*dn**4+23536.0_cp &
      &              *b**4*dm**2*dn**3-39776.0_cp*b**4*dm**2*dn**2-309040.0_cp*b**4*dm** &
      &              2*dn-244064.0_cp*b**4*dm**2+2568.0_cp*b**4*dn**6-85008.0_cp*b**4*dn**5 &
      &              +976496.0_cp*b**4*dn**4-3885536.0_cp*b**4*dn**3-3279704.0_cp*b**4*dn**2 &
      &              +34804112.0_cp*b**4*dn+33134208.0_cp*b**4)/(2048.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp)*(dn+2.0_cp))
      stencil(ku+10)=-a**9*b*(72.0_cp*a**2*dm**4*dn-24.0_cp*a**2*dm**4-116.0_cp &
      &              *a**2*dm**2*dn**3+2383.0_cp*a**2*dm**2*dn**2-13680.0_cp*a**2*dm**2*dn &
      &              +11685.0_cp*a**2*dm**2+28.0_cp*a**2*dn**5-1049.0_cp*a**2*dn**4+15899.0_cp &
      &              *a**2*dn**3-127925.0_cp*a**2*dn**2+578321.0_cp*a**2*dn-1196538.0_cp &
      &              *a**2+192.0_cp*a*b*dm**4*dn+192.0_cp*a*b*dm**4-512.0_cp*a*b*dm**2*dn**3 &
      &              +9496.0_cp*a*b*dm**2*dn**2-39456.0_cp*a*b*dm**2*dn-49464.0_cp*a*b*dm**2 &
      &              +320.0_cp*a*b*dn**5-12328.0_cp*a*b*dn**4+174408.0_cp*a*b*dn**3-1040720.0_cp &
      &              *a*b*dn**2+1793080.0_cp*a*b*dn+3020856.0_cp*a*b-128.0_cp*b**2*dm**4 &
      &              *dn-128.0_cp*b**2*dm**4+704.0_cp*b**2*dm**2*dn**3-12868.0_cp*b**2*dm**2 &
      &              *dn**2+52016.0_cp*b**2*dm**2*dn+65588.0_cp*b**2*dm**2-720.0_cp*b**2 &
      &              *dn**5+27936.0_cp*b**2*dn**4-397376.0_cp*b**2*dn**3+2377788.0_cp*b**2* &
      &              dn**2-4087504.0_cp*b**2*dn-6891324.0_cp*b**2)/(2048.0_cp*dn*(dn-3.0_cp &
      &              )*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp))
      stencil(ku+11)=-a**10*(6.0_cp*a**2*dm**4*dn-8.0_cp*a**2*dm**2*dn**3 &
      &              +177.0_cp*a**2*dm**2*dn**2-1043.0_cp*a**2*dm**2*dn+446.0_cp*a**2*dm**2 &
      &              +2.0_cp*a**2*dn**5-83.0_cp*a**2*dn**4+1335.0_cp*a**2*dn**3-10428.0_cp &
      &              *a**2*dn**2+39836.0_cp*a**2*dn-59856.0_cp*a**2+16.0_cp*a*b*dm**4*dn &
      &              +16.0_cp*a*b*dm**4-32.0_cp*a*b*dm**2*dn**3+660.0_cp*a*b*dm**2*dn**2 &
      &              -3100.0_cp*a*b*dm**2*dn-3792.0_cp*a*b*dm**2+16.0_cp*a*b*dn**5-676.0_cp &
      &              *a*b*dn**4+10508.0_cp*a*b*dn**3-69264.0_cp*a*b*dn**2+136304.0_cp*a* &
      &              b*dn+216768.0_cp*a*b-48.0_cp*b**2*dm**4*dn-48.0_cp*b**2*dm**4+164.0_cp &
      &              *b**2*dm**2*dn**3-3342.0_cp*b**2*dm**2*dn**2+15342.0_cp*b**2*dm**2*dn &
      &              +18848.0_cp*b**2*dm**2-124.0_cp*b**2*dn**5+5272.0_cp*b**2*dn**4-82346.0_cp &
      &              *b**2*dn**3+544074.0_cp*b**2*dn**2-1068080.0_cp*b**2*dn-1699896.0_cp &
      &              *b**2)/(2048.0_cp*dn*(dn-3.0_cp)*(dn-2.0_cp)*(dn-1.0_cp)*(dn+1.0_cp &
      &              ))
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
