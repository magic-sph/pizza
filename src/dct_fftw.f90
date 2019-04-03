module dct_fftw

   use iso_c_binding
   use precision_mod
   use namelists, only: fftw_plan_flag
   use mem_alloc, only: bytes_allocated
   use constants, only: half

   implicit none

   include 'fftw3.f03'

   private

   !complex(cp), target, allocatable :: work(:,:)
   !real(cp), pointer :: work_r(:,:)

   type, public :: costf_t
      real(cp) :: cheb_fac
      type(c_ptr) :: plan
      type(c_ptr) :: plan_1d
      complex(cp), pointer :: work(:,:)
      real(cp), pointer :: work_r(:,:)
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
      integer :: inembed(1), istride, idist, plan_size(1)
      integer :: onembed(1), ostride, odist, isize
      integer(C_INT) :: plan_type(1)
      real(cp) :: array_in(1:2*(nMstop-nMstart+1), n_r_max)
      real(cp) :: array_out(1:2*(nMstop-nMstart+1), n_r_max)
      real(cp) :: array_in_1d(n_r_max)
      real(cp) :: array_out_1d(n_r_max)
      logical :: l_work_array

      inembed(1) = 0
      onembed(1) = 0
      plan_size(1) = n_r_max
      plan_type(1) = FFTW_REDFT00
      istride = 2*(nMstop-nMstart+1)
      ostride = 2*(nMstop-nMstart+1)
      isize   = 2*(nMstop-nMstart+1)
      idist   = 1
      odist   = 1

      this%plan = fftw_plan_many_r2r(1, plan_size, isize, array_in,         &
                  &                  inembed, istride, idist, array_out,    &
                  &                  onembed, ostride, odist,               &
                  &                  plan_type, fftw_plan_flag)

      this%plan_1d = fftw_plan_r2r(1, plan_size, array_in_1d, array_out_1d, &
                     &             plan_type, fftw_plan_flag)

      this%cheb_fac = sqrt(half/(n_r_max-1))

      if ( present(no_work_array) ) then
         l_work_array= .not. no_work_array
      else
         l_work_array=.true.
      end if

      if ( l_work_array ) then
         allocate( this%work(nMstart:nMstop,n_r_max) )
         call c_f_pointer(c_loc(this%work), this%work_r, [isize, n_r_max])

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
         deallocate( this%work )
      end if

      call fftw_destroy_plan(this%plan_1d)
      call fftw_destroy_plan(this%plan)

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine costf_complex_2d(this,array_in,nMstart,nMstop,n_r_max)

      class(costf_t), intent(in) :: this

      !-- Input variables
      integer, intent(in) :: nMstart
      integer, intent(in) :: nMstop
      integer, intent(in) :: n_r_max

      !-- Output variables:
      complex(cp), target, intent(inout) :: array_in(nMstart:nMstop,n_r_max)

      !-- Local variables
      real(cp), pointer :: r_input(:,:)
      integer :: n_r, n_m, isize

      isize = 2*(nMstop-nMstart+1)

      call c_f_pointer(c_loc(array_in), r_input, [isize, n_r_max])
      call fftw_execute_r2r(this%plan, r_input, this%work_r)

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            array_in(n_m,n_r)=this%cheb_fac*this%work(n_m,n_r)
         end do
      end do

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
      integer :: n_r

      call fftw_execute_r2r(this%plan_1d, array_in, work_1d)

      do n_r=1,n_r_max
         array_in(n_r)=this%cheb_fac*work_1d(n_r)
      end do

   end subroutine costf_real_1d
!------------------------------------------------------------------------------
end module dct_fftw
