module LMLoop_mod

   use precision_mod
   use parallel_mod, only: rank
   use truncation_3D, only: n_r_max_3D
   use namelists, only: l_heat_3D
   use blocking_lm, only: llm, ulm
   use update_temp_3D_mod, only: update_temp_3D, finish_exp_temp_3D
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray


   implicit none

   private

   public :: LMLoop, finish_explicit_assembly_3D

contains

   subroutine LMLoop(temp_3D_LMloc, dtemp_3D_LMloc, dTdt_3D, tscheme, lMat)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lMat

      !-- Output variables
      complex(cp),         intent(out) :: temp_3D_LMloc(llm:ulm,n_r_max_3D)
      complex(cp),         intent(out) :: dtemp_3D_LMloc(llm:ulm,n_r_max_3D)
      type(type_tarray),   intent(inout) :: dTdt_3D

      !--- Local counter
      integer :: nLMB

      nLMB=1+rank

      if ( l_heat_3D ) then 
         call update_temp_3D( temp_3D_LMloc, dtemp_3D_LMloc, dTdt_3D, &
              &               tscheme, lMat, nLMB )
      end if

   end subroutine LMLoop
!--------------------------------------------------------------------------------
   subroutine finish_explicit_assembly_3D(dVrT_LMloc, dTdt_3D, tscheme)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      complex(cp),         intent(inout) :: dVrT_LMloc(llm:ulm,n_r_max_3D)

      !-- Output variables
      type(type_tarray),   intent(inout) :: dTdt_3D

      if ( l_heat_3D ) then
         call finish_exp_temp_3D(dVrT_LMloc, dTdt_3D%expl(:,:,tscheme%istage))
      end if

   end subroutine finish_explicit_assembly_3D
!--------------------------------------------------------------------------------
end module LMLoop_mod
