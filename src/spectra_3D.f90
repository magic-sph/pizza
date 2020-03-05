module spectra_3D
   !
   ! This module handles the calculation and the writing of spectra for spherical quantities.
   !   -spec_3D_#.TAG files for snapshot spectra_3D
   !   -spec_3D_avg.TAG files for time-averaged spectra (and standard deviation)
   !   -2D_spec_avg.TAG files which correspond to 2-D spectra in a (r,lm) plane

   use precision_mod
   use parallel_mod
   use constants, only: half
   use mem_alloc, only: bytes_allocated
   use communications, only: reduce_radial_on_rank
   use blocking, only: lmStart, lmStop
   use truncation_3D, only: n_r_max_3D, l_max
   use blocking_lm, only: lo_map
   use namelists, only: tag, l_mag_3D, DyMagFac
   use useful, only: cc2real, getMSD2, round_off
   use radial_functions, only: rscheme_3D, r_3D, or2_3D
   use integration, only: rInt_R
   use mean_sd, only: mean_sd_type

   implicit none

   private

   type, public :: spectra_3D_type
      type(mean_sd_type) :: bpol2L
      type(mean_sd_type) :: btor2L
      integer :: n_calls
      integer :: ispec_counter
      real(cp) :: dt
      real(cp) :: timeLast
      logical :: l_calc
   contains 
      procedure :: initialize
      procedure :: finalize
      procedure :: compute_spectra_3D
   end type spectra_3D_type

contains

   subroutine initialize(this)
      !
      ! Memory allocation and initial values
      !
      class(spectra_3D_type) :: this

      if ( l_mag_3D ) then
         call this%bpol2L%initialize(1,l_max)
         call this%btor2L%initialize(1,l_max)
      end if

      this%ispec_counter = 1
      this%n_calls = 0
      this%dt = 0.0_cp
      this%timeLast = 0.0_cp
      this%l_calc = .false.

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation
      !
      class(spectra_3D_type) :: this

      if ( l_mag_3D ) then
         call this%bpol2L%finalize()
         call this%btor2L%finalize()
      end if

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine compute_spectra_3D(this, time, l_stop_time, b_3D, db_3D, &
              &                    aj_3D)

      class(spectra_3D_type) :: this

      !-- Input variables:
      real(cp),    intent(in) :: time
      logical,     intent(in) :: l_stop_time
      complex(cp), intent(in) :: b_3D(lmStart:lmStop,n_r_max_3D)
      complex(cp), intent(in) :: db_3D(lmStart:lmStop,n_r_max_3D)
      complex(cp), intent(in) :: aj_3D(lmStart:lmStop,n_r_max_3D)

      !-- Local variables:
      real(cp) :: E_pol_l(l_max), E_tor_l(l_max)
      real(cp) :: e_pol_r_l(n_r_max_3D,l_max), e_tor_r_l(n_r_max_3D,l_max)
      real(cp) :: bpol2_SD_global(l_max), btor2_SD_global(l_max)
      real(cp) :: e_pol_tmp, e_tor_tmp, dLh, Efac
      integer :: n_r, lm, l, m
      integer :: file_handle

      this%n_calls = this%n_calls+1
      this%dt = time-this%timeLast

      !-- This is not cache-friendly but hopefully it's happening only
      !-- once in a while (otherwise we need (n_r, lm) arrays
      Efac        = 0.0_cp
      E_pol_l(:)   = 0.0_cp
      E_tor_l(:)   = 0.0_cp
      do n_r=1,n_r_max_3D
         e_pol_r_l(n_r,:)    = 0.0_cp
         e_tor_r_l(n_r,:)    = 0.0_cp
         do lm=max(lmStart,2),lmStop
            l  =lo_map%lm2l(lm)
            m  =lo_map%lm2m(lm)
            dLh= real(l*(l+1),cp)

            e_pol_tmp= dLh*( dLh*or2_3D(n_r)* cc2real(b_3D(lm,n_r),m)     &
            &                                +cc2real(db_3D(lm,n_r),m) )
            e_tor_tmp= dLh * cc2real(aj_3D(lm,n_r),m)

            !----- l-spectra:
            e_pol_r_l(n_r,l) = e_pol_r_l(n_r,l) + e_pol_tmp
            e_tor_r_l(n_r,l) = e_tor_r_l(n_r,l) + e_tor_tmp
            !!----- m-spectra:
            !e_pol_r_m(n_r,mc) = e_pol_r_m(n_r,l) + e_pol_tmp
            !e_tor_r_m(n_r,mc) = e_tor_r_m(n_r,l) + e_pol_tmp
         end do
      end do
      call reduce_radial_on_rank(e_pol_r_l, n_r_max_3D, 0)
      call reduce_radial_on_rank(e_tor_r_l, n_r_max_3D, 0)
      !call reduce_radial_on_rank(e_pol_r_m, n_r_max_3D, 0)
      !call reduce_radial_on_rank(e_pol_r_m, n_r_max_3D, 0)

      if ( rank == 0 ) then
         Efac=half*DyMagFac!*eScale
         do l=1,l_max
         !--Radial integration
            E_pol_l(l) =Efac*rInt_R(e_pol_r_l(:,l),r_3D,rscheme_3D)
            E_tor_l(l) =Efac*rInt_R(e_tor_r_l(:,l),r_3D,rscheme_3D)
            !-- Mean and SD of l-spectra
            call getMSD2(this%bpol2L%mean(l), this%bpol2L%SD(l), E_pol_l(l), &
                 &       this%n_calls, this%dt, time)
            call getMSD2(this%btor2L%mean(l), this%btor2L%SD(l), E_tor_l(l), &
                 &       this%n_calls, this%dt, time)
         end do
         !do m=1,l_max+1 ! Note: counter m is actual order+1
         !   E_pol_m(l) =Efac*rInt_R(e_pol_r_l(:,m),r_3D,rscheme_3D)
         !   E_tor_m(l) =Efac*rInt_R(e_tor_r_l(:,m),r_3D,rscheme_3D)
         !end do
      end if

      this%timeLast = time

      !if ( l_stop_time ) call this%write_spectra_avg()

      !----------------
      !- First write the time-average l spectra
      !----------------

      !-- Only rank==0 writes the spec_avg.TAG file
      bpol2_SD_global(:) = 0.0_cp
      btor2_SD_global(:) = 0.0_cp
      if ( rank == 0 ) then
         open(newunit=file_handle, file='spec_3D_avg.'//tag)
         do l=1,l_max
            bpol2_SD_global(l) =sqrt(this%bpol2L%SD(l)/this%timeLast)
            btor2_SD_global(l) =sqrt(this%btor2L%SD(l)/this%timeLast)
            write(file_handle, '(I5, 4es16.8)') l,                                &
            &     round_off(this%bpol2L%mean(l)), round_off(bpol2_SD_global(l)), &
            &     round_off(this%btor2L%mean(l)), round_off(btor2_SD_global(l))
         end do
         close(file_handle)
      end if

   end subroutine compute_spectra_3D
!----------------------------------------------------------------------
end module spectra_3D
