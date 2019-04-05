module rloop_3D

   use precision_mod
   use constants, only: ci
   use nonlinear_lm_mod, only: nonlinear_lm_t
   use grid_space_arrays_mod, only: grid_space_arrays_t
   use blocking, only: nRstart3D, nRstop3D, nRstart, nRstop
   use namelists, only: BuoFac
   use truncation_3D, only: lm_max, lmP_max, n_phi_max_3D, n_theta_max
   use truncation, only: idx2m, n_m_max
   use shtns, only: spat_to_SH, scal_to_spat
   use z_functions, only: zfunc_type

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
   subroutine radial_loop_3D( ur, ut, up, temp, dtempdt, dVrTLM, dpsidt_Rloc, &
              &               zinterp)

      !-- Input variables
      complex(cp), intent(in) :: temp(lm_max, nRstart3D:nRstop3D)
      real(cp),    intent(in) :: ur(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp),    intent(in) :: ut(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp),    intent(in) :: up(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      type(zfunc_type), intent(in) :: zinterp

      !-- Output variables
      complex(cp), intent(out) :: dtempdt(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(out) :: dVrTLM(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(inout) :: dpsidt_Rloc(n_m_max, nRstart:nRstop)

      !-- Local variables
      real(cp) :: buo_tmp(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      complex(cp) :: buo_tmp_Rloc(n_m_max,nRstart:nRstop)
      integer :: n_r, n_m, m

      do n_r=nRstart3D,nRstop3D

         call transform_to_grid_space_shtns(temp(:,n_r), gsa)

         call gsa%get_nl(ur(:,:,n_r), ut(:,:,n_r), up(:,:,n_r), n_r, &
              &          buo_tmp(:,:,n_r))

         call transform_to_lm_space_shtns(gsa, nl_lm)

         call nl_lm%get_td(dVrTLM(:,n_r), dtempdt(:,n_r))

      end do

      !-- Compute z-averaging of buoyancy
      call zinterp%compute_avg(buo_tmp, buo_tmp_Rloc)

#ifdef DEBUG
      block

         use truncation, only: n_r_max
         use radial_functions, only: r

         integer :: n_r, file_handle

         open(newunit=file_handle, file='buo', status='new')
         do n_r=1,n_r_max
            write(file_handle, '(3ES20.12)') r(n_r), real(buo_tmp_Rloc(6,n_r)), real(buo_tmp_Rloc(1,n_r))
         end do
         close(file_handle)

         open(newunit=file_handle, file='buo_3D', status='new', form='unformatted',&
         &    access='stream')
         write(file_handle) buo_tmp
         close(file_handle)

      end block
#endif

      !-- Finish assembling buoyancy and sum it with dpsidt
      do n_m=1,n_m_max
         m = idx2m(n_m)
         do n_r=nRstart,nRstop
            dpsidt_Rloc(n_m,n_r)=dpsidt_Rloc(n_m,n_r)-BuoFac*ci*m* &
            &                    buo_tmp_Rloc(n_m,n_r)
         end do
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
