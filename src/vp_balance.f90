module vp_balance
   !
   ! This module is used to compute the balance of the axisymmetric
   ! uphi equation.
   !

   use precision_mod
   use time_schemes, only: type_tscheme
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, m2idx
   use blocking, only: nMstart, nMstop, l_rank_has_m0
   use namelists, only: ra, pr, ek, radratio, sc, raxi, tag
   use radial_functions, only: r

   implicit none

   private

   !-- This object is used to compute and store the azimuthal force balance
   type, public :: vp_bal_type
      real(cp), allocatable :: rey_stress(:)
      real(cp), allocatable :: dvpdt(:)
      real(cp), allocatable :: visc(:)
      real(cp), allocatable :: pump(:)
      logical :: l_calc
      integer :: n_calls
      integer :: n_vphi_bal_file
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: initialize_dvpdt
      procedure :: finalize_dvpdt
      procedure :: write_outputs
   end type vp_bal_type

contains

   subroutine initialize(this)
      !
      ! Memory allocation and opening of the files
      !

      class(vp_bal_type) :: this

      if ( l_rank_has_m0 ) then
         open(newunit=this%n_vphi_bal_file, file='vphi_bal.'//tag, &
         &    form='unformatted', status='new', access='stream')

         allocate( this%rey_stress(n_r_max) )
         allocate( this%dvpdt(n_r_max) )
         allocate( this%visc(n_r_max) )
         allocate( this%pump(n_r_max) )
         this%n_calls = 0
         this%l_calc = .false.
         bytes_allocated = bytes_allocated+4*n_r_max*SIZEOF_DEF_REAL
      end if

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation and closing of the files
      !

      class(vp_bal_type) :: this

      if ( l_rank_has_m0 ) then
         close(this%n_vphi_bal_file)
         deallocate( this%pump, this%visc )
         deallocate( this%dvpdt, this%rey_stress )
      end if
      
   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine initialize_dvpdt(this, up_Mloc, tscheme)

      class(vp_bal_type) :: this

      !-- Input variables
      complex(cp),         intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      class(type_tscheme), intent(in) :: tscheme

      !-- Local variables:
      integer :: n_r, m0

      if ( l_rank_has_m0 ) then
         m0 = m2idx(0)
         do n_r=1,n_r_max
            this%dvpdt(n_r)=real(up_Mloc(m0,n_r))/tscheme%dt(1)
         end do
      end if

   end subroutine initialize_dvpdt
!------------------------------------------------------------------------------
   subroutine finalize_dvpdt(this, up_Mloc, tscheme)

      class(vp_bal_type) :: this

      !-- Input variables
      complex(cp),         intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      class(type_tscheme), intent(in) :: tscheme

      !-- Local variables:
      integer :: n_r, m0

      if ( l_rank_has_m0 ) then
         m0 = m2idx(0)
         do n_r=1,n_r_max
            this%dvpdt(n_r)=real(up_Mloc(m0,n_r))/tscheme%dt(1)-this%dvpdt(n_r)
         end do
      end if

   end subroutine finalize_dvpdt
!------------------------------------------------------------------------------
   subroutine write_outputs(this, time, up_Mloc)
      !
      ! This subroutine writes the vphi_bal.TAG file. This is a Fortran
      ! unformatted file.
      !

      class(vp_bal_type) :: this

      !-- Input variables
      real(cp),    intent(in) :: time
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)

      !-- Local variable
      integer :: idx_m0

      this%n_calls = this%n_calls+1

      if ( l_rank_has_m0 ) then
         idx_m0 = m2idx(0)
         if ( this%n_calls == 1 ) then
            write(this%n_vphi_bal_file) ra, ek, pr, radratio, raxi, sc
            write(this%n_vphi_bal_file) r
         end if
         write(this%n_vphi_bal_file) time
         write(this%n_vphi_bal_file) real(up_Mloc(idx_m0,:))
         write(this%n_vphi_bal_file) this%dvpdt
         write(this%n_vphi_bal_file) this%rey_stress
         write(this%n_vphi_bal_file) this%pump
         write(this%n_vphi_bal_file) this%visc
      end if

   end subroutine write_outputs
!------------------------------------------------------------------------------
end module vp_balance
