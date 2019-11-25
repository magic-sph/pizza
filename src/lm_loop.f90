module LMLoop_mod

   use precision_mod
   use parallel_mod, only: rank
   use truncation_3D, only: n_r_max_3D
   use namelists, only: l_heat_3D, l_mag_3D
   use blocking, only: lmStart, lmStop
   use update_temp_3D_mod, only: update_temp_3D, finish_exp_temp_3D
   use update_mag_3D_mod, only: update_mag_3D, finish_exp_mag_3D
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray


   implicit none

   private

   public :: LMLoop, finish_explicit_assembly_3D

contains

   subroutine LMLoop(temp_3D_LMloc, dtemp_3D_LMloc, dTdt_3D,         &
              &      b_3D_LMloc, db_3D_LMloc, ddb_3D_LMloc, dBdt_3D, &
              &      aj_3D_LMloc, dj_3D_LMloc, djdt_3D,              &
              &      tscheme, lMat)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lMat

      !-- Output variables
      complex(cp),       intent(inout) :: temp_3D_LMloc(lmStart:lmStop,n_r_max_3D)
      complex(cp),       intent(out) :: dtemp_3D_LMloc(lmStart:lmStop,n_r_max_3D)
      complex(cp),       intent(inout) :: b_3D_LMloc(lmStart:lmStop,n_r_max_3D)
      complex(cp),       intent(out) :: db_3D_LMloc(lmStart:lmStop,n_r_max_3D)
      complex(cp),       intent(out) :: ddb_3D_LMloc(lmStart:lmStop,n_r_max_3D)
      complex(cp),       intent(inout) :: aj_3D_LMloc(lmStart:lmStop,n_r_max_3D)
      complex(cp),       intent(out) :: dj_3D_LMloc(lmStart:lmStop,n_r_max_3D)
      type(type_tarray), intent(inout) :: dTdt_3D
      type(type_tarray), intent(inout) :: dBdt_3D
      type(type_tarray), intent(inout) :: djdt_3D

      !--- Local counter
      integer :: nLMB

      nLMB=1+rank

      if ( l_heat_3D ) then
         call update_temp_3D( temp_3D_LMloc, dtemp_3D_LMloc, dTdt_3D, &
              &               tscheme, lMat, nLMB )
      end if

      if ( l_mag_3D ) then
         call update_mag_3D( b_3D_LMloc, db_3D_LMloc, ddb_3D_LMloc, dBdt_3D, &
              &              aj_3D_LMloc, dj_3D_LMloc, djdt_3D,              &
              &              tscheme, lMat, nLMB )
      end if

   end subroutine LMLoop
!--------------------------------------------------------------------------------
   subroutine finish_explicit_assembly_3D(dVrT_LMloc, dTdt_3D, dVxBh_LMloc, &
              &                           djdt_3D, tscheme)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      complex(cp),         intent(inout) :: dVrT_LMloc(lmStart:lmStop,n_r_max_3D)
      complex(cp),         intent(inout) :: dVxBh_LMloc(lmStart:lmStop,n_r_max_3D)

      !-- Output variables
      type(type_tarray),   intent(inout) :: dTdt_3D
      type(type_tarray),   intent(inout) :: djdt_3D

      if ( l_heat_3D ) then
         call finish_exp_temp_3D(dVrT_LMloc, dTdt_3D%expl(:,:,tscheme%istage))
      end if

      if ( l_mag_3D ) then
         call finish_exp_mag_3D(dVxBh_LMloc, djdt_3D%expl(:,:,tscheme%istage))
      end if

   end subroutine finish_explicit_assembly_3D
!--------------------------------------------------------------------------------
end module LMLoop_mod
