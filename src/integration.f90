module integration
   !
   ! Radial integration functions
   !

   use precision_mod
   use constants, only: half, one, two
   use radial_scheme, only: type_rscheme

   implicit none

   private

   public :: rInt_R, simps

contains

   real(cp) function rInt_R(f,r,r_scheme) result(rInt)
      !
      !   Same as function rInt but for a radial dependent mapping function
      !   dr_fac2.
      !

      !-- Input variables:
      real(cp),            intent(in) :: f(:)    ! Input function
      real(cp),            intent(in) :: r(:)    ! Radius
      class(type_rscheme), intent(in) :: r_scheme! Radial scheme (FD or Cheb)
              
      !-- Local variables
      real(cp), allocatable :: f2(:)
      real(cp) :: h1, h2
      integer :: n_r, nCheb, n_r_max
                 
      n_r_max=size(f)

      !--- Integrals:
      if ( r_scheme%version == 'cheb' ) then

         allocate( f2(n_r_max) )

         do n_r=1,n_r_max
            f2(n_r)=f(n_r)/r_scheme%drx(n_r)
         end do

         !-- Transform to cheb space:
         call r_scheme%costf1(f2,n_r_max)
         f2(1)      =half*f2(1)
         f2(n_r_max)=half*f2(n_r_max)

         !-- Sum contribution:
         rInt=f2(1)            ! This is zero order contribution
         do nCheb=3,n_r_max,2 ! Only even chebs contribute
            rInt=rInt-one/real(nCheb*(nCheb-2),cp)*f2(nCheb)
         end do

         !-- Remaining renormalisation:
         rInt=two*sqrt(two/real(n_r_max-1,cp))*rInt

         deallocate( f2 )

      else

         if ( mod(n_r_max,2)==1 ) then ! Odd number (Simpson ok)

            rInt = 0.0_cp
            do n_r=2,n_r_max-1,2
               h2=r(n_r+1)-r(n_r)
               h1=r(n_r)-r(n_r-1)
               rInt=rInt+(h1+h2)/6.0_cp*( f(n_r-1)*(two*h1-h2)/h1     +&
               &                      f(n_r)  *(h1+h2)*(h1+h2)/(h1*h2)+&
               &                      f(n_r+1)*(two*h2-h1)/h2 )
            end do

            rInt = -rInt

         else ! Even number (twice simpson + trapz on the first and last points)

            rInt = half*(r(2)-r(1))*(f(2)+f(1))
            do n_r=3,n_r_max,2
               h2=r(n_r+1)-r(n_r)
               h1=r(n_r)-r(n_r-1)
               rInt=rInt+(h1+h2)/6.0_cp*( f(n_r-1)*(two*h1-h2)/h1     +&
               &                      f(n_r)  *(h1+h2)*(h1+h2)/(h1*h2)+&
               &                      f(n_r+1)*(two*h2-h1)/h2 )
            end do
            rInt = rInt+half*(r(n_r_max)-r(n_r_max-1))*(f(n_r_max)+f(n_r_max-1))
            do n_r=2,n_r_max-1,2
               h2=r(n_r+1)-r(n_r)
               h1=r(n_r)-r(n_r-1)
               rInt=rInt+(h1+h2)/6.0_cp*( f(n_r-1)*(two*h1-h2)/h1     +&
               &                      f(n_r)  *(h1+h2)*(h1+h2)/(h1*h2)+&
               &                      f(n_r+1)*(two*h2-h1)/h2 )
            end do
            rInt = -half*rInt

         end if

      end if

   end function rInt_R
!------------------------------------------------------------------------------
   real(cp) function simps(f,r) result(rInt)
      !
      ! Simpson's method to integrate a function
      !

      !-- Input variables:
      real(cp), intent(in) :: f(:)    ! Input function
      real(cp), intent(in) :: r(:)    ! Radius
              
      !-- Local variables
      real(cp) :: h1, h2
      integer :: n_r, n_r_max
                 
      n_r_max=size(f)

      if ( mod(n_r_max,2)==1 ) then ! Odd number (Simpson ok)

         rInt = 0.0_cp
         do n_r=2,n_r_max-1,2
            h2=r(n_r+1)-r(n_r)
            h1=r(n_r)-r(n_r-1)
            rInt=rInt+(h1+h2)/6.0_cp*( f(n_r-1)*(two*h1-h2)/h1     +&
            &                      f(n_r)  *(h1+h2)*(h1+h2)/(h1*h2)+&
            &                      f(n_r+1)*(two*h2-h1)/h2 )
         end do

         rInt = -rInt

      else ! Even number (twice simpson + trapz on the first and last points)

         rInt = half*(r(2)-r(1))*(f(2)+f(1))
         do n_r=3,n_r_max,2
            h2=r(n_r+1)-r(n_r)
            h1=r(n_r)-r(n_r-1)
            rInt=rInt+(h1+h2)/6.0_cp*( f(n_r-1)*(two*h1-h2)/h1     +&
            &                      f(n_r)  *(h1+h2)*(h1+h2)/(h1*h2)+&
            &                      f(n_r+1)*(two*h2-h1)/h2 )
         end do
         rInt = rInt+half*(r(n_r_max)-r(n_r_max-1))*(f(n_r_max)+f(n_r_max-1))
         do n_r=2,n_r_max-1,2
            h2=r(n_r+1)-r(n_r)
            h1=r(n_r)-r(n_r-1)
            rInt=rInt+(h1+h2)/6.0_cp*( f(n_r-1)*(two*h1-h2)/h1     +&
            &                      f(n_r)  *(h1+h2)*(h1+h2)/(h1*h2)+&
            &                      f(n_r+1)*(two*h2-h1)/h2 )
         end do
         rInt = -half*rInt

      end if

   end function simps
!------------------------------------------------------------------------------
end module integration
