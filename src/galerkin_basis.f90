module galerkin

   use precision_mod
   use constants, only: one, two
   use matrix_types, only: type_bandmat_real

   implicit none

   private

   interface galerkin2cheb
      module procedure galerkin2cheb_real
      module procedure galerkin2cheb_complex
   end interface galerkin2cheb

   public :: galerkin2cheb, get_galerkin_stencil, destroy_galerkin_stencil

contains

   subroutine get_galerkin_stencil(sten, n_r_max, i_bound_type)

      !-- Input variables
      integer, intent(in) :: i_bound_type
      integer, intent(in) :: n_r_max

      !-- Output variables
      type(type_bandmat_real), intent(out) :: sten

      !-- Local variables
      integer :: kl, ku, i, m


      select case (i_bound_type)

         case(1) ! u(ri)=u(ro)=0 ! Dirichlet on both boundaries

            kl = 2
            ku = 0
            call sten%initialize(kl, ku, n_r_max)
            do i=1,n_r_max
               sten%dat(1,i)=-one
               if ( i == 1 .or. i==n_r_max ) sten%dat(1,i)=two*sten%dat(1,i)
               if ( i < n_r_max-1 ) then
                  sten%dat(3,i) = one
                  if ( i == n_r_max-2 ) sten%dat(3,i)=two*sten%dat(3,i)
               end if
            end do

         case(2) ! u(ri)=u'(ri)=u(ro)=u'(ro)=0
            
            kl = 4
            ku = 0
            call sten%initialize(kl, ku, n_r_max)
            do i=1,n_r_max
               m = i-1
               sten%dat(1,i) = one
               if ( i == 1 .or. i==n_r_max ) sten%dat(1,i)=two*sten%dat(1,i)
               if ( i < n_r_max-1 ) then
                  sten%dat(3,i) = -two*real(m+2,cp)/real(m+3,cp)
                  if ( i == n_r_max-2 ) sten%dat(3,i)=two*sten%dat(3,i)
               end if
               if ( i < n_r_max-3 ) then
                  sten%dat(5,i) = real(m+1,cp)/real(m+3,cp)
                  if ( i == n_r_max-4 ) sten%dat(3,i)=two*sten%dat(3,i)
               end if
            end do

      end select

   end subroutine get_galerkin_stencil
!-------------------------------------------------------------------------------
   subroutine galerkin2cheb_real(sten, arr)

      !-- Input variable
      type(type_bandmat_real), intent(in) :: sten

      !-- Output variable
      real(cp), intent(inout) :: arr(sten%nlines)

      call sten%mat_vec_mul(arr)

   end subroutine galerkin2cheb_real
!-------------------------------------------------------------------------------
   subroutine galerkin2cheb_complex(sten, arr)

      !-- Input variable
      type(type_bandmat_real), intent(in) :: sten

      !-- Output variable
      complex(cp), intent(inout) :: arr(sten%nlines)

      call sten%mat_vec_mul(arr)

   end subroutine galerkin2cheb_complex
!-------------------------------------------------------------------------------
   subroutine destroy_galerkin_stencil(sten)
      !
      ! Memory deallocation of the stencils
      !
      type(type_bandmat_real), intent(in) :: sten

      call sten%finalize()

   end subroutine destroy_galerkin_stencil
!-------------------------------------------------------------------------------
end module galerkin
