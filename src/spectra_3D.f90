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
   use namelists, only: tag, l_mag_3D, DyMagFac, &
       &                l_2D_spectra, l_2D_SD
   use useful, only: cc2real, getMSD2, round_off
   use radial_functions, only: rscheme_3D, r_3D, or2_3D
   use integration, only: rInt_R
   use mean_sd, only: mean_sd_type, mean_sd_2D_type

   implicit none

   private

   type, public :: spectra_3D_type
      type(mean_sd_2D_type) :: bpol2LR
      type(mean_sd_2D_type) :: btor2LR
      type(mean_sd_2D_type) :: bpol2MR
      type(mean_sd_2D_type) :: btor2MR
      type(mean_sd_type) :: bpol2L
      type(mean_sd_type) :: btor2L
      type(mean_sd_type) :: bpol2M
      type(mean_sd_type) :: btor2M
      type(mean_sd_type) :: bpol2R
      type(mean_sd_type) :: btor2R
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
         if ( l_2D_spectra ) then
            call this%bpol2LR%initialize(1,l_max,n_r_max_3D,l_2D_SD)
            call this%btor2LR%initialize(1,l_max,n_r_max_3D,l_2D_SD)
            call this%bpol2MR%initialize(1,l_max+1,n_r_max_3D,l_2D_SD)
            call this%btor2MR%initialize(1,l_max+1,n_r_max_3D,l_2D_SD)
         end if

         call this%bpol2L%initialize(1,l_max)
         call this%btor2L%initialize(1,l_max)
         call this%bpol2M%initialize(1,l_max+1)
         call this%btor2M%initialize(1,l_max+1)
         call this%bpol2R%initialize(1,n_r_max_3D)
         call this%btor2R%initialize(1,n_r_max_3D)
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
         call this%btor2R%finalize()
         call this%bpol2R%finalize()
         call this%btor2M%finalize()
         call this%bpol2M%finalize()
         call this%btor2L%finalize()
         call this%bpol2L%finalize()
         if ( l_2D_spectra ) then
            call this%btor2MR%finalize()
            call this%bpol2MR%finalize()
            call this%btor2LR%finalize()
            call this%bpol2LR%finalize()
         end if
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
      real(cp) :: E_pol_m(l_max+1), E_tor_m(l_max+1)
      real(cp) :: E_pol_r(n_r_max_3D), E_tor_r(n_r_max_3D)
      real(cp) :: e_pol_r_l(l_max,n_r_max_3D), e_tor_r_l(l_max,n_r_max_3D)
      real(cp) :: e_pol_r_m(l_max+1,n_r_max_3D), e_tor_r_m(l_max+1,n_r_max_3D)
      real(cp) :: e_pol_tmp, e_tor_tmp, dLh, Efac
      integer :: n_r, lm, l, m
      integer :: file_handle

      this%n_calls = this%n_calls+1
      this%dt = time-this%timeLast

      !-- This is not cache-friendly but hopefully it's happening only
      !-- once in a while (otherwise we need (lm, n_r) arrays
      Efac       = 0.0_cp
      E_pol_l(:) = 0.0_cp
      E_tor_l(:) = 0.0_cp
      E_pol_m(:) = 0.0_cp
      E_tor_m(:) = 0.0_cp
      E_pol_r(:) = 0.0_cp
      E_tor_r(:) = 0.0_cp
      do n_r=1,n_r_max_3D
         e_pol_r_l(:,n_r) = 0.0_cp
         e_tor_r_l(:,n_r) = 0.0_cp
         do lm=max(lmStart,2),lmStop
            l  =lo_map%lm2l(lm)
            m  =lo_map%lm2m(lm)
            dLh= real(l*(l+1),cp)

            e_pol_tmp= dLh*( dLh*or2_3D(n_r)* cc2real(b_3D(lm,n_r),m)     &
            &                                +cc2real(db_3D(lm,n_r),m) )
            e_tor_tmp= dLh * cc2real(aj_3D(lm,n_r),m)

            !----- l-spectra:
            e_pol_r_l(l,n_r) = e_pol_r_l(l,n_r) + e_pol_tmp
            e_tor_r_l(l,n_r) = e_tor_r_l(l,n_r) + e_tor_tmp
            !----- m-spectra:
            !e_pol_r_m(m+1,n_r) = e_pol_r_m(m+1,n_r) + e_pol_tmp
            !e_tor_r_m(m+1,n_r) = e_tor_r_m(m+1,n_r) + e_pol_tmp
            !!----- r-spectra:
            E_pol_r(n_r) = E_pol_r(n_r) + e_pol_tmp
            E_tor_r(n_r) = E_tor_r(n_r) + e_pol_tmp
         end do
         !-- A reduction is needed here
         call reduce_radial_on_rank(e_pol_r_l(:,n_r), l_max, 0)
         call reduce_radial_on_rank(e_tor_r_l(:,n_r), l_max, 0)
         !call reduce_radial_on_rank(e_pol_r_m(:,n_r), l_max+1, 0)
         !call reduce_radial_on_rank(e_pol_r_m(:,n_r), l_max+1, 0)
       end do
       call reduce_radial_on_rank(E_pol_r, n_r_max_3D, 0)
       call reduce_radial_on_rank(E_tor_r, n_r_max_3D, 0)

      if ( rank == 0 ) then ! switch to rank 0 for the postprocess
         Efac=half*DyMagFac!*eScale
         do l=1,l_max
         !--Radial integration
            E_pol_l(l) =Efac*rInt_R(e_pol_r_l(l,:),r_3D,rscheme_3D)
            E_tor_l(l) =Efac*rInt_R(e_tor_r_l(l,:),r_3D,rscheme_3D)
            !-- Mean and SD of l-spectra
            call getMSD2(this%bpol2L%mean(l), this%bpol2L%SD(l), E_pol_l(l), &
                 &       this%n_calls, this%dt, time)
            call getMSD2(this%btor2L%mean(l), this%btor2L%SD(l), E_tor_l(l), &
                 &       this%n_calls, this%dt, time)
         end do
         !-- m-spectra
         !do m=1,l_max+1 !- Note: counter m is actual order+1
         !   E_pol_m(m) =Efac*rInt_R(e_pol_r_l(m,:),r_3D,rscheme_3D)
         !   E_tor_m(m) =Efac*rInt_R(e_tor_r_l(m,:),r_3D,rscheme_3D)
         !   call getMSD2(this%bpol2M%mean(m), this%bpol2M%SD(m), E_pol_m(m), &
         !        &       this%n_calls, this%dt, time)
         !   call getMSD2(this%btor2M%mean(m), this%btor2M%SD(m), E_tor_m(m), &
         !        &       this%n_calls, this%dt, time)
         !end do
         !-- r-spectra
         do n_r=1,n_r_max_3D
            call getMSD2(this%bpol2R%mean(n_r), this%bpol2R%SD(n_r), E_pol_r(n_r), &
                 &       this%n_calls, this%dt, time)
            call getMSD2(this%btor2R%mean(n_r), this%btor2R%SD(n_r), E_tor_r(n_r), &
                 &       this%n_calls, this%dt, time)
         end do

         !-- Averaging of 2D (l/m-R) spectra
         if ( l_2D_spectra ) then
            do n_r=1,n_r_max_3D
               do l=1,l_max
                  call getMSD2(this%bpol2LR%mean(l,n_r), this%bpol2LR%SD(l,n_r), &
                       &       e_pol_r_l(l,n_r),this%n_calls, this%dt, time)
                  call getMSD2(this%btor2LR%mean(l,n_r), this%btor2LR%SD(l,n_r), &
                       &       e_tor_r_l(l,n_r),this%n_calls, this%dt, time)
                  !call getMSD2(this%bpol2MR%mean(l+1,n_r), this%bpol2MR%SD(l+1,n_r), &
                  !     &       e_pol_r_m(l+1,n_r),this%n_calls, this%dt, time)
                  !call getMSD2(this%btor2MR%mean(l+1,n_r), this%btor2MR%SD(l+1,n_r), &
                  !     &       e_tor_r_m(l+1,n_r),this%n_calls, this%dt, time)
               end do
            end do
         end if
      end if

      this%timeLast = time

      if ( l_stop_time ) then !call this%write_spectra_avg()

         !----------------
         !- First write the time-average l spectra
         !----------------

         !-- Only rank==0 writes the spec_avg.TAG file
         if ( rank == 0 ) then
            open(newunit=file_handle, file='spec_3D_avg.'//tag)
            do l=1,l_max
               this%bpol2L%SD(l) =sqrt(this%bpol2L%SD(l)/this%timeLast)
               this%btor2L%SD(l) =sqrt(this%btor2L%SD(l)/this%timeLast)
               !this%bpol2M%SD(l+1) =sqrt(this%bpol2M%SD(l+1)/this%timeLast)
               !this%btor2M%SD(l+1) =sqrt(this%btor2M%SD(l+1)/this%timeLast)
               write(file_handle, '(I5, 8es16.8)') l,                              &
               &     round_off(this%bpol2L%mean(l)), round_off(this%bpol2L%SD(l)), &
               &     round_off(this%btor2L%mean(l)), round_off(this%btor2L%SD(l))!, &
               !&     round_off(this%bpol2M%mean(l+1)), round_off(this%bpol2M%SD(l+1)), &
               !&     round_off(this%btor2M%mean(l+1)), round_off(this%btor2M%SD(l+1))
            end do
            close(file_handle)

            open(newunit=file_handle, file='spec_r_3D_avg.'//tag)
            do n_r=1,n_r_max_3D
               this%bpol2R%SD(n_r) =sqrt(this%bpol2R%SD(n_r)/this%timeLast)
               this%btor2R%SD(n_r) =sqrt(this%btor2R%SD(n_r)/this%timeLast)
               write(file_handle, '(I5, 4es16.8)') n_r,                              &
               &     round_off(this%bpol2R%mean(n_r)), round_off(this%bpol2R%SD(n_r)), &
               &     round_off(this%btor2R%mean(n_r)), round_off(this%btor2R%SD(n_r))
            end do
            close(file_handle)

            if ( l_2D_spectra ) then
               open(newunit=file_handle, file='2D_spec_3D_avg.'//tag)
               do n_r=1,n_r_max_3D
                  do l=1,l_max
                     this%bpol2LR%SD(l,n_r) =sqrt(this%bpol2LR%SD(l,n_r)/this%timeLast)
                     this%btor2LR%SD(l,n_r) =sqrt(this%btor2LR%SD(l,n_r)/this%timeLast)
                     this%bpol2MR%SD(l+1,n_r) =sqrt(this%bpol2MR%SD(l+1,n_r)/this%timeLast)
                     this%btor2MR%SD(l+1,n_r) =sqrt(this%btor2MR%SD(l+1,n_r)/this%timeLast)
                     write(file_handle, '(2I5, 8es16.8)') l, n_r,                                  &
                     &     round_off(this%bpol2LR%mean(l,n_r)), round_off(this%bpol2LR%SD(l,n_r)), &
                     &     round_off(this%btor2LR%mean(l,n_r)), round_off(this%btor2LR%SD(l,n_r)), &
                     &     round_off(this%bpol2MR%mean(l+1,n_r)), round_off(this%bpol2MR%SD(l+1,n_r)), &
                     &     round_off(this%btor2MR%mean(l+1,n_r)), round_off(this%btor2MR%SD(l+1,n_r))
                  end do
               end do
            close(file_handle)
            end if
         end if

      end if

   end subroutine compute_spectra_3D
!----------------------------------------------------------------------
end module spectra_3D
