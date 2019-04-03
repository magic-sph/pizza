module rloop

   use precision_mod
   use parallel_mod
   use constants, only: ci, one, half
   use mem_alloc, only: bytes_allocated
   use namelists, only: ek, tadvz_fac, CorFac
   use radial_functions, only: or1, r, beta, oheight, dtcond, ekpump
   use nonlinear_lm_mod, only: get_td
   use general_arrays_mod, only: get_nl_shtns
   use blocking, only: nRstart, nRstop, nRstart3D, nRstop3D
   use truncation, only: n_m_max, n_phi_max, idx2m, m2idx
   use truncation_3D, only: n_m_max, n_phi_max, idx2m, m2idx
   use courant_mod, only: courant
   use fourier, only: fft, ifft
   use useful, only: cc22real
   use time_schemes, only: type_tscheme
   use timers_mod, only: timers_type

   implicit none

   private

   public :: radial_loop, initialize_radial_loop, finalize_radial_loop

contains

!------------------------------------------------------------------------------
   subroutine radial_loop_3D( temp_Rloc, dtempdt_Rloc, dVrT_Rloc)

      !-- Input variables
      complex(cp), intent(in) :: temp_Rloc(lm_max, nRstart3D:nRstop3D)

      !-- Output variables
      complex(cp), intent(out) :: dtempdt_Rloc(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(out) :: dVrT_Rloc(lm_max, nRstart3D:nRstop3D)

      !-- Local variables
      real(cp) :: usom, runStart, runStop
      complex(cp) :: us_fluct
      integer :: n_r, n_phi, n_m, m, idx_m0

      do n_r=nRstart3D,nRstop3D

         call transform_to_grid_space_shtns(this%gsa, time)

         call gsa%get_nl(ur_3D_Rloc(:,:,n_r), ut_3D_Rloc(:,:,n_r), &
              &          up_3D_Rloc(:,:,n_r), this%n_r)

         call transform_to_lm_space_shtns(this%gsa, this%nl_lm)

         !-- Partial calculation of time derivatives (horizontal parts):
         !   input flm...  is in (l,m) space at radial grid points n_r !
         !   get_td finally calculates the d*dt terms needed for the
         !   time step performed in 's_LMLoop.f' . This should be distributed
         !   over the different models that 's_LMLoop.f' parallelizes over.
         call get_td(this%n_r, dVTrLM, dtempdt_Rloc)

      end do

   end subroutine radial_loop_3D
!-------------------------------------------------------------------------------
   subroutine transform_to_grid_space_shtns(this, temp, gsa)

      complex(cp), intent(in) :: temp(lm_max)
      type(grid_space_arrays_t) :: gsa

      call scal_to_spat(temp, gsa%Tc)

   end subroutine transform_to_grid_space_shtns
!-------------------------------------------------------------------------------
   subroutine transform_to_lm_space_shtns(this, gsa, nl_lm)

      class(rIterThetaBlocking_shtns_t) :: this
      type(grid_space_arrays_t) :: gsa
      type(nonlinear_lm_t) :: nl_lm

      call shtns_load_cfg(1)

      call spat_to_SH(gsa%VTr, nl_lm%VTrLM)
      call spat_to_SH(gsa%VTt, nl_lm%VTtLM)
      call spat_to_SH(gsa%VTp, nl_lm%VTpLM)

      call shtns_load_cfg(0)

   end subroutine transform_to_lm_space_shtns
!------------------------------------------------------------------------------
end module rloop
