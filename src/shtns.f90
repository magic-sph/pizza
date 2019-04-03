module shtns

   use precision_mod, only: cp
   use constants, only: ci
   use truncation_3D, only: m_max_3D, l_max, n_theta_max, n_phi_max_3D, &
       &                    minc_3D, lm_max, lmP_max
   use parallel_mod

   implicit none

   include "shtns.f"

   private

   public :: init_shtns, scal_to_spat, spat_to_SH

contains

   subroutine init_shtns()

      integer :: norm

      if ( rank == 0 ) then
         write(*,*) ''
         call shtns_verbose(1)
      end if

      call shtns_use_threads(0)

      norm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE

      call shtns_set_size(l_max, m_max_3D/minc_3D, minc_3D, norm)
      !call shtns_precompute(SHT_GAUSS, SHT_PHI_CONTIGUOUS, &
      call shtns_precompute(SHT_QUICK_INIT, SHT_PHI_CONTIGUOUS, &
           &                1.e-10_cp, n_theta_max, n_phi_max_3D)
      call shtns_save_cfg(0)

      if ( rank == 0 ) then
         call shtns_verbose(0)
      end if

      call shtns_set_size(l_max+1, m_max_3D/minc_3D, minc_3D, norm)
      call shtns_precompute(SHT_QUICK_INIT, SHT_PHI_CONTIGUOUS, &
           &                1.e-10_cp, n_theta_max, n_phi_max_3D)
      call shtns_save_cfg(1)

      call shtns_load_cfg(0)

   end subroutine
!------------------------------------------------------------------------------
   subroutine scal_to_spat(Slm, fieldc)
      ! transform a spherical harmonic field into grid space
      complex(cp), intent(in) :: Slm(lm_max)
      real(cp), intent(out) :: fieldc(n_phi_max_3D, n_theta_max)

      call shtns_SH_to_spat(Slm, fieldc)

   end subroutine scal_to_spat
!------------------------------------------------------------------------------
   subroutine spat_to_SH(f, fLM)

      real(cp), intent(in) :: f(n_phi_max_3D, n_theta_max)
      complex(cp), intent(out) :: fLM(lmP_max)

      call shtns_load_cfg(1)
      call shtns_spat_to_sh(f, fLM)
      call shtns_load_cfg(0)

   end subroutine spat_to_SH
!------------------------------------------------------------------------------
end module shtns
