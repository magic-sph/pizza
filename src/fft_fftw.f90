module fourier
   
   use iso_c_binding
   use namelists, only: fftw_plan_flag
   use constants, only: zero, half, two
   use precision_mod

   implicit none

   include 'fftw3.f03'

   private

   type(c_ptr) :: plan_forward, plan_backward
   type(c_ptr) :: plan_backward_3D, plan_forward_3D

   public :: initialize_fourier, finalize_fourier, fft, ifft

contains

   subroutine initialize_fourier(n_phi_max, n_phi_max_3D, l_3D)

      !-- Input variable
      integer, intent(in) :: n_phi_max
      integer, intent(in) :: n_phi_max_3D
      logical, intent(in) :: l_3D

      !-- Local variables
      real(cp), allocatable :: real_array(:)
      complex(cp), allocatable :: complex_array(:)

      allocate( real_array(n_phi_max), complex_array(n_phi_max/2+1) )
      plan_forward = fftw_plan_dft_r2c_1d(n_phi_max,real_array,complex_array, &
                     &                    fftw_plan_flag)
      plan_backward = fftw_plan_dft_c2r_1d(n_phi_max,complex_array,real_array, &
                      &                    fftw_plan_flag)
      deallocate( real_array, complex_array )

      if ( l_3D ) then
         allocate( real_array(n_phi_max_3D), complex_array(n_phi_max_3D/2+1) )
         plan_forward_3D = fftw_plan_dft_r2c_1d(n_phi_max_3D,real_array, &
                           &                    complex_array, fftw_plan_flag)
         plan_backward_3D = fftw_plan_dft_c2r_1d(n_phi_max_3D,complex_array, &
                            &                    real_array, fftw_plan_flag)
         deallocate( real_array, complex_array )
      end if

   end subroutine initialize_fourier
!------------------------------------------------------------------------------
   subroutine finalize_fourier(l_3D)

      logical, intent(in) :: l_3D

      if ( l_3D ) then
         call fftw_destroy_plan(plan_backward_3D)
         call fftw_destroy_plan(plan_forward_3D)
      end if
      call fftw_destroy_plan(plan_forward)
      call fftw_destroy_plan(plan_backward)

   end subroutine finalize_fourier
!------------------------------------------------------------------------------
   subroutine fft(array_in, array_out, l_3D)

      !-- Input variables
      real(cp),          intent(inout) :: array_in(:)
      logical, optional, intent(in) :: l_3D

      !-- Output variable
      complex(cp), intent(out) :: array_out(:)

      !-- Local variables
      complex(cp) :: tmp(size(array_in)/2+1)
      integer :: n_m, n_phi_max, n_m_max
      logical :: l_3D_loc

      l_3D_loc = .false.
      if ( present(l_3D) ) l_3D_loc=l_3D

      n_phi_max = size(array_in)
      n_m_max = size(array_out)
      
      if ( l_3D_loc ) then
         call fftw_execute_dft_r2c(plan_forward_3D, array_in, tmp)
      else
         call fftw_execute_dft_r2c(plan_forward, array_in, tmp)
      end if

      do n_m=1,n_m_max
         array_out(n_m) = tmp(n_m)/real(n_phi_max,cp)
      end do

   end subroutine fft
!------------------------------------------------------------------------------
   subroutine ifft(array_in, array_out, l_3D)

      !-- Input variables
      complex(cp),       intent(in) :: array_in(:)
      logical, optional, intent(in) :: l_3D

      !-- Output variable
      real(cp), intent(out) :: array_out(:)

      !-- Local variables
      complex(cp) :: tmp(size(array_out)/2+1)
      integer :: n_m, n_m_max, n_phi_max
      logical :: l_3D_loc

      l_3D_loc = .false.
      if ( present(l_3D) ) l_3D_loc=l_3D

      n_m_max = size(array_in)
      n_phi_max = size(array_out)

      do n_m=1,n_m_max
         tmp(n_m)= array_in(n_m)
      end do
      do n_m=n_m_max+1,n_phi_max/2+1
         tmp(n_m)=zero
      end do

      if ( l_3D_loc ) then
         call fftw_execute_dft_c2r(plan_backward_3D, tmp, array_out)
      else
         call fftw_execute_dft_c2r(plan_backward, tmp, array_out)
      end if

   end subroutine ifft
!------------------------------------------------------------------------------
end module fourier
