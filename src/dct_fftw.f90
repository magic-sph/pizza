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
   use constants, only: half, zero, ci, two, one, third

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
      integer, allocatable :: der(:), der2(:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure, private :: costf_complex_2d
      procedure, private :: costf_real_1d
      generic :: costf => costf_real_1d, costf_complex_2d
      procedure :: get_dr_fft
      procedure :: get_ddr_fft
   end type costf_t

contains

   subroutine initialize(this, nMstart, nMstop, n_r_max, n_cheb_max, no_work_array)

      class(costf_t) :: this

      integer, intent(in) :: nMstart
      integer, intent(in) :: nMstop
      integer, intent(in) :: n_r_max
      integer, intent(in) :: n_cheb_max
      logical, optional, intent(in) :: no_work_array

      !--Local variables
      integer :: plan_size(1)
      integer(C_INT) :: plan_type(1)
#if   (DCT_VERSION==dft_loop)
      integer :: k
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
      allocate ( this%der(2*n_r_max-2), this%der2(2*n_r_max-2) )
      bytes_allocated=bytes_allocated+(2*n_r_max-2)*SIZEOF_INTEGER
      this%der(:)=0
      this%der(2*n_r_max-2:2*n_r_max-n_cheb_max:-1)=[(-k,k=1,n_cheb_max-1)]
      this%der(1:n_cheb_max)=[(k-1,k=1,n_cheb_max)]
      this%der2(:)=this%der(:)*this%der(:)
      this%der(n_r_max)=0
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
      deallocate(this%der,this%der2)
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
   subroutine get_dr_fft(this,f,df,xcheb,nMstart,nMstop,n_r_max,n_cheb_max,l_dct_in)

      class(costf_t), intent(in) :: this

      !-- Input variables
      integer,     intent(in) :: nMstart ! Starting index (OMP)
      integer,     intent(in) :: nMstop  ! Stopping index (OMP)
      integer,     intent(in) :: n_r_max  ! Max cheb
      integer,     intent(in) :: n_cheb_max  ! Max cheb
      real(cp),    intent(in) :: xcheb(n_r_max) ! Gauss-Lobatto grid
      complex(cp), intent(in) :: f(nMstart:nMstop,n_r_max) ! Array to be transformed
      logical,     intent(in) :: l_dct_in ! Do we need a DCT for the input array?

      !-- Output variables:
      complex(cp), intent(out) :: df(nMstart:nMstop,n_r_max)  ! Radial derivative

#if (DCT_VERSION==dft_loop)
      !-- Local variables:
      integer :: n_m, k
      complex(cp) :: tot
      complex(cp) :: work_1d(2*n_r_max-2), work_1d_out(2*n_r_max-2)

      do n_m=nMstart,nMstop
         if ( l_dct_in ) then
            work_1d(1:n_r_max) =f(n_m,1:n_r_max)
            work_1d(n_r_max+1:)=f(n_m,n_r_max-1:2:-1)
            call fftw_execute_dft(this%plan_fft_forw, work_1d, work_1d_out)
         else
            work_1d_out(1:n_r_max) =f(n_m,1:n_r_max) / this%cheb_fac
            work_1d_out(n_r_max+1:)=f(n_m,n_r_max-1:2:-1) / this%cheb_fac
         end if

         work_1d_out(n_r_max)=half*work_1d_out(n_r_max)

         !--  Boundary points = tau lines
         tot=zero
         do k=1,n_cheb_max-1
            tot=tot+k**2 * work_1d_out(k+1)
         end do
         df(n_m,1)=tot/(n_r_max-1)
         tot=zero
         do k=1,n_cheb_max-1
            tot=tot+(-1)**(k+1)*k**2*work_1d_out(k+1)
         end do
         df(n_m,n_r_max)=tot/(n_r_max-1)

         work_1d_out(n_r_max)=two*work_1d_out(n_r_max)

         !-- Derivatives in FFT space
         work_1d_out(:)=ci*this%der(:)*work_1d_out(:)
         call fftw_execute_dft(this%plan_fft_back, work_1d_out, work_1d)

         !-- Bring back to Gauss-Lobatto grid for bulk points
         df(n_m,2:n_r_max-1)=-work_1d(2:n_r_max-1) / sqrt(one-xcheb(2:n_r_max-1)**2) /&
         &                   (2*n_r_max-2)
      end do
#endif

   end subroutine get_dr_fft
!------------------------------------------------------------------------------
   subroutine get_ddr_fft(this,f,df,ddf,xcheb,nMstart,nMstop,n_r_max,n_cheb_max,l_dct_in)

      class(costf_t), intent(in) :: this

      !-- Input variables
      integer,     intent(in) :: nMstart ! Starting index (OMP)
      integer,     intent(in) :: nMstop  ! Stopping index (OMP)
      integer,     intent(in) :: n_r_max   ! Number of radial grid points
      integer,     intent(in) :: n_cheb_max  ! Max cheb
      real(cp),    intent(in) :: xcheb(n_r_max) ! Gauss-Lobatto grid
      complex(cp), intent(in) :: f(nMstart:nMstop,n_r_max) ! Array to be transformed
      logical,     intent(in) :: l_dct_in ! Do we need a DCT for the input array?

      !-- Output variables:
      complex(cp), intent(out) :: df(nMstart:nMstop,n_r_max)  ! Radial derivative
      complex(cp), intent(out) :: ddf(nMstart:nMstop,n_r_max)  ! 2nd radial derivative

#if (DCT_VERSION==dft_loop)
      !-- Local variables:
      integer :: n_m, k
      complex(cp) :: tot
      complex(cp) :: work_1(2*n_r_max-2), tmp(2*n_r_max-2)
      complex(cp) :: work_2(2*n_r_max-2)

      do n_m=nMstart,nMstop
         if ( l_dct_in ) then
            work_1(1:n_r_max) =f(n_m,1:n_r_max)
            work_1(n_r_max+1:)=f(n_m,n_r_max-1:2:-1)
            call fftw_execute_dft(this%plan_fft_forw, work_1, tmp)
         else
            tmp(1:n_r_max) =f(n_m,1:n_r_max) / this%cheb_fac
            tmp(n_r_max+1:)=f(n_m,n_r_max-1:2:-1) / this%cheb_fac
         end if

         tmp(n_r_max)=half*tmp(n_r_max)

         !--  Boundary points = tau lines
         tot=zero
         do k=1,n_cheb_max-1
            tot=tot+k**2 * tmp(k+1)
         end do
         df(n_m,1)=tot/(n_r_max-1)
         tot=zero
         do k=2,n_cheb_max-1
            tot=tot+k**2*(k**2-1) * tmp(k+1)
         end do
         ddf(n_m,1)=third * tot/(n_r_max-1)
         tot=zero
         do k=1,n_cheb_max-1
            tot=tot+(-1)**(k+1)*k**2*tmp(k+1)
         end do
         df(n_m,n_r_max)=tot/(n_r_max-1)
         tot=zero
         do k=2,n_cheb_max-1
            tot=tot+(-1)**k*k**2*(k**2-1)*tmp(k+1)
         end do
         ddf(n_m,n_r_max)=third*tot/(n_r_max-1)

         tmp(n_r_max)=two*tmp(n_r_max)

         !-- Derivatives in Fourier space
         work_2(:)=ci*this%der(:)*tmp(:)
         call fftw_execute_dft(this%plan_fft_back, work_2, work_1)
         tmp(:)   =-this%der2(:)*tmp(:)
         call fftw_execute_dft(this%plan_fft_back, tmp, work_2)

         !-- Bring back to Gauss-Lobatto grid for bulk points
         df(n_m,2:n_r_max-1)=-work_1(2:n_r_max-1) /                 &
         &                         sqrt(one-xcheb(2:n_r_max-1)**2) /&
         &                         (2*n_r_max-2)
         ddf(n_m,2:n_r_max-1)=-work_1(2:n_r_max-1)*                       &
         &                          xcheb(2:n_r_max-1) /                  &
         &                         (one-xcheb(2:n_r_max-1)**2)**1.5_cp /  &
         &                         (2*n_r_max-2)+work_2(2:n_r_max-1)/     &
         &                         (one-xcheb(2:n_r_max-1)**2)/           &
         &                         (2*n_r_max-2)
      end do
#endif

   end subroutine get_ddr_fft
!------------------------------------------------------------------------------
end module dct_fftw
