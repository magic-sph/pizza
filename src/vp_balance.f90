module vp_balance
   !
   ! This module is used to compute the balance of the axisymmetric
   ! uphi equation.
   !

   use precision_mod
   use time_schemes, only: type_tscheme
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, m2idx, n_m_max
   use blocking, only: nMstart, nMstop, l_rank_has_m0, nRstart, nRstop
   use namelists, only: ra, pr, ek, radratio, sc, raxi, tag, l_finite_diff
   use radial_functions, only: r
   use communications, only: gather_from_Rloc

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
      integer :: nRmin, nRmax, mMin, mMax
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

      !-- Local variables:

      if ( l_finite_diff ) then
         this%nRmin = nRstart
         this%nRmax = nRstop
         this%mMin = 1
         this%mMax = n_m_max
      else
         this%nRmin = 1
         this%nRmax = n_r_max
         this%mMin = nMstart
         this%mMax = nMstop
      end if

      if ( l_finite_diff .or. l_rank_has_m0 ) then
         allocate( this%rey_stress(this%nRmin:this%nRmax) )
         allocate( this%dvpdt(this%nRmin:this%nRmax) )
         allocate( this%visc(this%nRmin:this%nRmax) )
         allocate( this%pump(this%nRmin:this%nRmax) )
         bytes_allocated = bytes_allocated+4*(this%nRmax-this%nRmin+1)*SIZEOF_DEF_REAL
      end if

      if ( l_rank_has_m0 ) then
         open(newunit=this%n_vphi_bal_file, file='vphi_bal.'//tag, &
         &    form='unformatted', status='new', access='stream')

         this%n_calls = 0
         this%l_calc = .false.
      end if

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation and closing of the files
      !

      class(vp_bal_type) :: this

      if ( l_finite_diff .or. l_rank_has_m0 ) then
         deallocate( this%pump, this%visc, this%dvpdt, this%rey_stress )
      end if
      if ( l_rank_has_m0 ) then
         close(this%n_vphi_bal_file)
      end if
      
   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine initialize_dvpdt(this, up, tscheme)

      class(vp_bal_type) :: this

      !-- Input variables
      complex(cp),         intent(in) :: up(this%mMin:this%mMax,this%nRmin:this%nRmax)
      class(type_tscheme), intent(in) :: tscheme

      !-- Local variables:
      integer :: n_r, m0

      if ( l_finite_diff .or. l_rank_has_m0 ) then
         m0 = m2idx(0)
         do n_r=this%nRmin,this%nRmax
            this%dvpdt(n_r)=real(up(m0,n_r))/tscheme%dt(1)
         end do
      end if

   end subroutine initialize_dvpdt
!------------------------------------------------------------------------------
   subroutine finalize_dvpdt(this, up, tscheme)

      class(vp_bal_type) :: this

      !-- Input variables
      complex(cp),         intent(in) :: up(this%mMin:this%mMax,this%nRmin:this%nRmax)
      class(type_tscheme), intent(in) :: tscheme

      !-- Local variables:
      integer :: n_r, m0

      if ( l_finite_diff .or. l_rank_has_m0 ) then
         m0 = m2idx(0)
         do n_r=this%nRmin,this%nRmax
            this%dvpdt(n_r)=real(up(m0,n_r))/tscheme%dt(1)-this%dvpdt(n_r)
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
      real(cp) :: dvpdt_glob(n_r_max), rey_glob(n_r_max), pump_glob(n_r_max)
      real(cp) :: visc_glob(n_r_max)
      integer :: idx_m0

      this%n_calls = this%n_calls+1

      if ( l_finite_diff ) then ! mpi gather on rank 0
         call gather_from_Rloc(this%dvpdt, dvpdt_glob, 0)
         call gather_from_Rloc(this%rey_stress, rey_glob, 0)
         call gather_from_Rloc(this%pump, pump_glob, 0)
         call gather_from_Rloc(this%visc, visc_glob, 0)
      else
         dvpdt_glob = this%dvpdt
         rey_glob = this%rey_stress
         pump_glob = this%pump
         visc_glob = this%visc
      end if

      if ( l_rank_has_m0 ) then
         idx_m0 = m2idx(0)
         if ( this%n_calls == 1 ) then
            write(this%n_vphi_bal_file) ra, ek, pr, radratio, raxi, sc
            write(this%n_vphi_bal_file) r
         end if
         write(this%n_vphi_bal_file) time
         write(this%n_vphi_bal_file) real(up_Mloc(idx_m0,:))
         write(this%n_vphi_bal_file) dvpdt_glob
         write(this%n_vphi_bal_file) rey_glob
         write(this%n_vphi_bal_file) pump_glob
         write(this%n_vphi_bal_file) visc_glob
      end if

   end subroutine write_outputs
!------------------------------------------------------------------------------
end module vp_balance
