module outputs_3D
   !
   ! This module computes the outputs that are related to fiels expressed in 3-D
   !

   use parallel_mod
   use precision_mod
   use constants, only: osq4pi
   use truncation_3D, only: n_r_max_3D
   use blocking, only: lmStart, lmStop
   use namelists, only: l_heat_3D, tag
   use radial_functions, only: rscheme_3D, dtcond_3D, tcond_3D
   use time_schemes, only: type_tscheme
   use radial_der, only: get_dr

   implicit none

   private

   integer :: n_heat_file

   public :: initialize_outputs_3D, finalize_outputs_3D, write_outputs_3D


contains

   subroutine initialize_outputs_3D

      character(len=144) :: file_name

      if ( rank == 0 ) then
         if ( l_heat_3D ) then
            file_name = 'heat_3D.'//tag
            open(newunit=n_heat_file, file=file_name, status='new')
         end if
      end if

   end subroutine initialize_outputs_3D
!---------------------------------------------------------------------------------
   subroutine finalize_outputs_3D

      if ( rank == 0 ) then
         if ( l_heat_3D ) close(n_heat_file)
      end if

   end subroutine finalize_outputs_3D
!---------------------------------------------------------------------------------
   subroutine write_outputs_3D(time, tscheme, n_time_step, l_log, temp_3D)

      !-- Input variables
      real(cp),            intent(in) :: time
      class(type_tscheme), intent(in) :: tscheme
      integer,             intent(in) :: n_time_step
      logical,             intent(in) :: l_log

      complex(cp),         intent(in) :: temp_3D(lmStart:lmStop,n_r_max_3D)

      if ( l_log ) then
         call write_time_series(time, temp_3D)
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
end module outputs_3D
