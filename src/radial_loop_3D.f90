module rloop_3D

   use precision_mod
   use nonlinear_lm_mod, only: nonlinear_lm_t
   use grid_space_arrays_mod, only: grid_space_arrays_t
   use blocking, only: nRstart3D, nRstop3D
   use truncation_3D, only: lm_max, lmP_max, n_phi_max_3D, n_theta_max
   use shtns, only: spat_to_SH, scal_to_spat

   implicit none

   private

   type(grid_space_arrays_t) :: gsa
   type(nonlinear_lm_t) :: nl_lm

   public :: radial_loop_3D, initialize_radial_loop_3D, finalize_radial_loop_3D

contains

   subroutine initialize_radial_loop_3D(lmP_max)

      integer, intent(in) :: lmP_max

      call gsa%initialize()
      call nl_lm%initialize(lmP_max)

   end subroutine initialize_radial_loop_3D
!------------------------------------------------------------------------------
   subroutine finalize_radial_loop_3D

      call nl_lm%finalize()
      call gsa%finalize()

   end subroutine finalize_radial_loop_3D
!------------------------------------------------------------------------------
   subroutine radial_loop_3D( ur, ut, up, temp, dtempdt, dVrTLM)

      !-- Input variables
      complex(cp), intent(in) :: temp(lm_max, nRstart3D:nRstop3D)
      real(cp),    intent(in) :: ur(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp),    intent(in) :: ut(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp),    intent(in) :: up(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)

      !-- Output variables
      complex(cp), intent(out) :: dtempdt(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(out) :: dVrTLM(lm_max, nRstart3D:nRstop3D)

      !-- Local variables
      integer :: n_r

      do n_r=nRstart3D,nRstop3D

         call transform_to_grid_space_shtns(temp(:,n_r), gsa)

         call gsa%get_nl(ur(:,:,n_r), ut(:,:,n_r), up(:,:,n_r), n_r)

         call transform_to_lm_space_shtns(gsa, nl_lm)

         !-- Partial calculation of time derivatives (horizontal parts):
         !   input flm...  is in (l,m) space at radial grid points n_r !
         !   get_td finally calculates the d*dt terms needed for the
         !   time step performed in 's_LMLoop.f' . This should be distributed
         !   over the different models that 's_LMLoop.f' parallelizes over.
         call nl_lm%get_td(dVrTLM(:,n_r), dtempdt(:,n_r))

      end do

   end subroutine radial_loop_3D
!-------------------------------------------------------------------------------
   subroutine transform_to_grid_space_shtns(temp, gsa)

      complex(cp), intent(in) :: temp(lm_max)
      type(grid_space_arrays_t) :: gsa

      call scal_to_spat(temp, gsa%Tc)

   end subroutine transform_to_grid_space_shtns
!-------------------------------------------------------------------------------
   subroutine transform_to_lm_space_shtns(gsa, nl_lm)

      type(grid_space_arrays_t) :: gsa
      type(nonlinear_lm_t) :: nl_lm

      call shtns_load_cfg(1)

      call spat_to_SH(gsa%VTr, nl_lm%VTrLM)
      call spat_to_SH(gsa%VTt, nl_lm%VTtLM)
      call spat_to_SH(gsa%VTp, nl_lm%VTpLM)

      call shtns_load_cfg(0)

   end subroutine transform_to_lm_space_shtns
!------------------------------------------------------------------------------
end module rloop_3D
