#define dct_many 0
#define dft_loop 1
#define dft_many 2
#define DCT_VERSION dft_loop
module dct_fftw

   use iso_c_binding
   use precision_mod
   use namelists, only: fftw_plan_flag
   use mem_alloc, only: bytes_allocated
   use blocking, only: nm_per_rank
   use constants, only: half

   implicit none

   include 'fftw3.f03'

   private

   complex(cp), target, allocatable :: work(:,:)
   real(cp), pointer :: work_r(:,:)

   type, public :: costf_t
      real(cp) :: cheb_fac
#if   (DCT_VERSION==dct_many)
      type(c_ptr) :: plan
#elif (DCT_VERSION==dft_loop)
      type(c_ptr) :: plan_fft_back, plan_fft_forw   ! FFTW single plan for FFT
#elif (DCT_VERSION==dft_many)
      type(c_ptr) :: plan_fft_forw
#endif
      type(c_ptr) :: plan_1d
   contains
      procedure :: initialize
      procedure :: finalize
      procedure, private :: costf_complex_2d
      procedure, private :: costf_real_1d
      generic :: costf => costf_real_1d, costf_complex_2d
   end type costf_t

contains

   subroutine initialize(this, nMstart, nMstop, n_r_max, no_work_array)

      class(costf_t) :: this

      integer, intent(in) :: nMstart
      integer, intent(in) :: nMstop
      integer, intent(in) :: n_r_max
      logical, optional, intent(in) :: no_work_array

      !--Local variables
      integer :: plan_size(1)
      integer(C_INT) :: plan_type(1)
#if   (DCT_VERSION==dft_loop)
      complex(cp) :: array_cplx_1d(2*n_r_max-2), array_cplx_out_1d(2*n_r_max-2)
#elif (DCT_VERSION==dct_many)
      integer :: inembed(1), istride, idist
      integer :: onembed(1), ostride, odist
      real(cp) :: array_in(2*(nMstop-nMstart+1), n_r_max)
      real(cp) :: array_out(2*(nMstop-nMstart+1), n_r_max)
#elif (DCT_VERSION==dft_many)
      integer :: inembed(1), istride, idist
      integer :: onembed(1), ostride, odist
      complex(cp) :: array_in(nMstop-nMstart+1, 2*n_r_max-2)
      complex(cp) :: array_out(nMstop-nMstart+1, 2*n_r_max-2)
#endif
      real(cp) :: array_in_1d(n_r_max), array_out_1d(n_r_max)
      logical :: l_work_array

#if   (DCT_VERSION==dft_loop)
      plan_size(1) = 2*n_r_max-2
      this%plan_fft_back = fftw_plan_dft(1, plan_size, array_cplx_1d,      &
                           &             array_cplx_out_1d, FFTW_BACKWARD, &
                           &             fftw_plan_flag)
      this%plan_fft_forw = fftw_plan_dft(1, plan_size, array_cplx_1d,      &
                           &             array_cplx_out_1d, FFTW_FORWARD,  &
                           &             fftw_plan_flag)
#elif (DCT_VERSION==dft_many)
      plan_size(1) = 2*n_r_max-2
      inembed(1) = 0
      onembed(1) = 0
      istride = nMstop-nMstart+1
      ostride = nMstop-nMstart+1
      idist   = 1
      odist   = 1
      this%plan_fft_forw = fftw_plan_many_dft(1, plan_size, nM_per_rank,          &
                           &                  array_in, inembed, istride, idist,  &
                           &                  array_out, onembed, ostride, odist, &
                           &                  FFTW_FORWARD,fftw_plan_flag)
#elif (DCT_VERSION==dct_many)
      inembed(1) = 0
      onembed(1) = 0
      plan_size(1) = n_r_max
      plan_type(1) = FFTW_REDFT00
      istride = 2*(nMstop-nMstart+1)
      ostride = 2*(nMstop-nMstart+1)
      idist   = 1
      odist   = 1

      this%plan = fftw_plan_many_r2r(1, plan_size, 2*nM_per_rank, array_in, &
                  &                  inembed, istride, idist, array_out,    &
                  &                  onembed, ostride, odist,               &
                  &                  plan_type, fftw_plan_flag)
#endif

      plan_size(1) = n_r_max
      plan_type(1) = FFTW_REDFT00
      this%plan_1d = fftw_plan_r2r(1, plan_size, array_in_1d, array_out_1d, &
                     &             plan_type, fftw_plan_flag)

      this%cheb_fac = sqrt(half/(n_r_max-1))

      if ( present(no_work_array) ) then
         l_work_array= .not. no_work_array
      else
         l_work_array=.true.
      end if

      if ( l_work_array ) then
         allocate( work(nMstart:nMstop,n_r_max) )
         call c_f_pointer(c_loc(work), work_r, [2*nm_per_rank, n_r_max])

         bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*n_r_max* &
         &                 SIZEOF_DEF_COMPLEX
      end if

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this, no_work_array)

      class(costf_t) :: this

      !-- Input variable:
      logical, optional, intent(in) :: no_work_array

      !--Local variable:
      logical :: l_work_array

      if ( present(no_work_array) ) then
         l_work_array=.not. no_work_array
      else
         l_work_array=.true.
      end if

      if ( l_work_array ) then
         deallocate( work )
      end if

      call fftw_destroy_plan(this%plan_1d)
#if   (DCT_VERSION==dct_many)
      call fftw_destroy_plan(this%plan)
#elif (DCT_VERSION==dft_loop)
      call fftw_destroy_plan(this%plan_fft_back)
      call fftw_destroy_plan(this%plan_fft_forw)
#elif (DCT_VERSION==dft_many)
      call fftw_destroy_plan(this%plan_fft_forw)
#endif

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine costf_complex_2d(this,array_in,nMstart,nMstop,n_r_max)

      class(costf_t), intent(in) :: this

      !-- Input variables
      integer, intent(in) :: nMstart
      integer, intent(in) :: nMstop
      integer, intent(in) :: n_r_max

#if (DCT_VERSION==dft_loop)
      !-- Output variables:
      complex(cp), intent(inout) :: array_in(nMstart:nMstop,n_r_max)

      !-- Local variables:
      integer :: n_m
      complex(cp) :: work_1d(2*n_r_max-2), work_1d_out(2*n_r_max-2)

      do n_m=nMstart,nMstop
         work_1d(1:n_r_max) = array_in(n_m,1:n_r_max)
         work_1d(n_r_max+1:)= array_in(n_m,n_r_max-1:2:-1)
         call fftw_execute_dft(this%plan_fft_forw, work_1d, work_1d_out)
         array_in(n_m,:)=this%cheb_fac*work_1d_out(1:n_r_max)
      end do

#elif (DCT_VERSION==dft_many)
      !-- Output variables:
      complex(cp), intent(inout) :: array_in(nMstart:nMstop,n_r_max)

      !-- Local variables:
      integer :: n_r
      complex(cp) :: tmp_in(nMstart:nMstop,2*n_r_max-2)
      complex(cp) :: tmp_out(nMstart:nMstop,2*n_r_max-2)

      do n_r=1,n_r_max
         tmp_in(:,n_r)=array_in(:,n_r)
      end do
      do n_r=n_r_max+1,2*n_r_max-2
         tmp_in(:,n_r)=array_in(:,2*n_r_max-n_r)
      end do

      call fftw_execute_dft(this%plan_fft_forw, tmp_in, tmp_out)

      !-- Copy output onto a
      do n_r=1,n_r_max
         array_in(:,n_r)=this%cheb_fac*tmp_out(:,n_r)
      end do

#elif (DCT_VERSION==dct_many)
      !-- Output variables:
      complex(cp), target, intent(inout) :: array_in(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r
      real(cp), pointer :: r_input(:,:)

      call c_f_pointer(c_loc(array_in), r_input, [2*nm_per_rank, n_r_max])
      call fftw_execute_r2r(this%plan, r_input, work_r)

      do n_r=1,n_r_max
         array_in(nMstart:nMstop,n_r)=this%cheb_fac*work(nMstart:nMstop,n_r)
      end do
#endif

   end subroutine costf_complex_2d
!------------------------------------------------------------------------------
   subroutine costf_real_1d(this, array_in, n_r_max)

      class(costf_t), intent(in) :: this

      !-- Input variable:
      integer, intent(in) :: n_r_max

      !-- Output variables:
      real(cp), intent(inout) :: array_in(n_r_max)

      !-- Local variables
      real(cp) :: work_1d(n_r_max)

      call fftw_execute_r2r(this%plan_1d, array_in, work_1d)
      array_in(:)=this%cheb_fac*work_1d(:)

   end subroutine costf_real_1d
!------------------------------------------------------------------------------
end module dct_fftw
