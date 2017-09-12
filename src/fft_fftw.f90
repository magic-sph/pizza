module fourier
   
   use iso_c_binding
   use namelists, only: fftw_plan_flag
   use constants, only: zero, half, two
   use truncation, only: n_m_max, n_phi_max
   use precision_mod

   implicit none

   include 'fftw3.f03'

   private

   type(c_ptr) :: plan_forward, plan_backward

   public :: initialize_fourier, finalize_fourier, fft, ifft

contains

   subroutine initialize_fourier(n_phi_max)

      !-- Input variable
      integer, intent(in) :: n_phi_max

      !-- Local variables
      real(cp) :: real_array(n_phi_max)
      complex(cp) :: complex_array(n_phi_max/2+1)

      plan_forward = fftw_plan_dft_r2c_1d(n_phi_max,real_array,complex_array, &
                     &                    fftw_plan_flag)
      plan_backward = fftw_plan_dft_c2r_1d(n_phi_max,complex_array,real_array, &
                      &                    fftw_plan_flag)

   end subroutine initialize_fourier
!------------------------------------------------------------------------------
   subroutine finalize_fourier

      call fftw_destroy_plan(plan_forward)
      call fftw_destroy_plan(plan_backward)

   end subroutine finalize_fourier
!------------------------------------------------------------------------------
   subroutine fft(array_in, array_out)

      real(cp), intent(inout) :: array_in(*)
      complex(cp), intent(out) :: array_out(*)

      complex(cp) :: tmp(n_phi_max/2+1)
      integer :: n_m
      
      call fftw_execute_dft_r2c(plan_forward, array_in, tmp)

      do n_m=1,n_m_max
         array_out(n_m) = tmp(n_m)/real(n_phi_max,cp)
      end do

   end subroutine fft
!------------------------------------------------------------------------------
   subroutine ifft(array_in, array_out)

      complex(cp), intent(in) :: array_in(*)

      real(cp), intent(out) :: array_out(*)

      complex(cp) :: tmp(n_phi_max/2+1)
      integer :: n_m

      do n_m=1,n_m_max
         tmp(n_m)= array_in(n_m)
      end do
      do n_m=n_m_max+1,n_phi_max/2+1
         tmp(n_m)=zero
      end do

      call fftw_execute_dft_c2r(plan_backward, tmp, array_out)

   end subroutine ifft
!------------------------------------------------------------------------------
end module fourier
