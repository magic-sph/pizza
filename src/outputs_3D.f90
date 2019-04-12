module outputs_3D
   !
   ! This module computes the outputs that are related to fiels expressed in 3-D
   !

   use parallel_mod
   use precision_mod
   use constants, only: osq4pi
   use truncation_3D, only: n_r_max_3D, n_theta_max, n_phi_max_3D
   use blocking, only: lmStart, lmStop, nRstart3D, nRstop3D
   use namelists, only: l_heat_3D, tag
   use radial_functions, only: rscheme_3D, dtcond_3D, tcond_3D, r_3D
   use radial_der, only: get_dr
   use mean_sd, only: mean_sd_type
   use time_schemes, only: type_tscheme
   use useful, only: round_off, getMSD2
   use output_frames, only: write_snapshot_3D

   implicit none

   private

   type(mean_sd_type) :: tempR
   integer :: n_heat_file, n_calls
   real(cp) :: timeLast_rad, timeAvg_rad

   public :: initialize_outputs_3D, finalize_outputs_3D, write_outputs_3D

contains

   subroutine initialize_outputs_3D

      character(len=144) :: file_name

      timeAvg_rad = 0.0_cp

      if ( rank == 0 ) then
         if ( l_heat_3D ) then
            file_name = 'heat_3D.'//tag
            open(newunit=n_heat_file, file=file_name, status='new')
         end if
         call tempR%initialize(1,n_r_max_3D)
         n_calls     =0
         timeLast_rad=0.0_cp
      end if

   end subroutine initialize_outputs_3D
!---------------------------------------------------------------------------------
   subroutine finalize_outputs_3D

      if ( rank == 0 ) then
         call tempR%finalize()
         if ( l_heat_3D ) close(n_heat_file)
      end if

   end subroutine finalize_outputs_3D
!---------------------------------------------------------------------------------
   subroutine write_outputs_3D(time, tscheme, l_log, l_stop_time, temp_3D)

      !-- Input variables
      real(cp),            intent(in) :: time
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: l_log
      logical,             intent(in) :: l_stop_time
      complex(cp),         intent(in) :: temp_3D(lmStart:lmStop,n_r_max_3D)

      timeAvg_rad  = timeAvg_rad  + tscheme%dt(1)

      if ( l_log ) then
         call write_time_series(time, temp_3D)
         call get_radial_averages(timeAvg_rad, l_stop_time, temp_3D)
      end if

   end subroutine write_outputs_3D
!---------------------------------------------------------------------------------
   subroutine write_time_series(time, temp_3D)

      !-- Input variables
      real(cp),    intent(in) :: time
      complex(cp), intent(in) :: temp_3D(lmStart:lmStop,n_r_max_3D)

      !-- Local variables
      real(cp) :: temp0(n_r_max_3D), dtemp0(n_r_max_3D)
      real(cp) :: NuTop, NuBot, tTop, tBot, NuDelta, beta_t

      if ( rank == 0 ) then
         temp0(:) = osq4pi*real(temp_3D(1, :))
         call get_dr(temp0, dtemp0, n_r_max_3D, rscheme_3D)

         !-- Nusselt at top and bottom
         NuTop = dtemp0(1)/dtcond_3D(1)
         NuBot = dtemp0(n_r_max_3D)/dtcond_3D(n_r_max_3D)

         !-- Nusset based on the ratio of temperature contrast
         NuDelta = (tcond_3D(n_r_max_3D)-tcond_3D(1)) / &
         &         (temp0(n_r_max_3D)-temp0(1))

         !-- Temperature gradient at mid-shell
         beta_t = dtemp0(n_r_max_3D/2)

         !-- Top and bottom temperature
         tTop = temp0(1)
         tBot = temp0(n_r_max_3D)

         write(n_heat_file, '(1P, ES20.12, 7ES16.8)') time, NuTop, NuBot,   &
         &                                            NuDelta, tTop, tBot,  &
         &                                            beta_t

      end if

   end subroutine write_time_series
!---------------------------------------------------------------------------------
   subroutine get_radial_averages(timeAvg_rad, l_stop_time, temp_3D)

      !-- Input variables
      real(cp),    intent(in) :: timeAvg_rad
      logical,     intent(in) :: l_stop_time
      complex(cp), intent(in) :: temp_3D(lmStart:lmStop,n_r_max_3D)

      !-- Local variables
      real(cp) :: dtAvg, dat
      integer :: n_r, file_handle

      if ( rank == 0 ) then

         n_calls = n_calls+1
         dtAvg = timeAvg_rad-timeLast_rad

         do n_r=1,n_r_max_3D
            if ( l_heat_3D ) then
               dat = real(temp_3D(1,n_r))*osq4pi
               call getMSD2(tempR%mean(n_r), tempR%SD(n_r), dat, &
                    &       n_calls, dtAvg, timeAvg_rad)
            end if
         end do
         timeLast_rad = timeAvg_rad

         if ( l_stop_time ) then
            open(newunit=file_handle, file='radial_profiles_3D.'//tag)
            do n_r=1,n_r_max_3D
               if ( l_heat_3D ) then
                  tempR%SD(n_r)=sqrt(tempR%SD(n_r)/timeAvg_rad)
               end if
               write(file_handle, '(es20.12, 2es16.8)') r_3D(n_r),        &
               &     round_off(tempR%mean(n_r)), round_off(tempR%SD(n_r))
            end do
            close(file_handle)
         end if

      end if

   end subroutine get_radial_averages
!---------------------------------------------------------------------------------
end module outputs_3D
