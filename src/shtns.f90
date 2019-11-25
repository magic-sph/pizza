module shtns

   use precision_mod, only: cp
   use iso_fortran_env, only: output_unit
   use constants, only: ci
   use truncation_3D, only: m_max_3D, l_max, n_theta_max, n_phi_max_3D, &
       &                    minc_3D, lm_max, lmP_max
   use horizontal, only: dLh
   use radial_functions, only: or2_3D, or1_3D
   use parallel_mod

   implicit none

   include "shtns.f"

   private

   public :: init_shtns, scal_to_spat, spat_to_SH, torpol_to_spat, &
   &         torpol_to_curl_spat, scal_axi_to_grad_spat

contains

   subroutine init_shtns()

      integer :: norm

      if ( rank == 0 ) then
         write(output_unit,*) ''
         call shtns_verbose(1)
      end if

      call shtns_use_threads(0)

      norm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE

      call shtns_set_size(l_max, m_max_3D/minc_3D, minc_3D, norm)
!      call shtns_precompute(SHT_GAUSS, SHT_PHI_CONTIGUOUS, &
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
   subroutine torpol_to_spat(Wlm, dWlm, Zlm, n_r, fieldrc, fieldtc, fieldpc)
      complex(cp), intent(in) :: Wlm(lm_max), dWlm(lm_max), Zlm(lm_max)
      integer, intent(in) :: n_r
      real(cp), intent(out) :: fieldrc(n_phi_max_3D, n_theta_max)
      real(cp), intent(out) :: fieldtc(n_phi_max_3D, n_theta_max)
      real(cp), intent(out) :: fieldpc(n_phi_max_3D, n_theta_max)

      !-- Local variable
      complex(cp) :: Qlm(lm_max)
      integer :: lm

      do lm = 1, lm_max
         Qlm(lm) = dLh(lm) * or2_3D(n_r) * Wlm(lm)
         !Qlm(lm) = dLh(lm) * Wlm(lm)
      end do

      call shtns_qst_to_spat(Qlm, dWlm, Zlm, fieldrc, fieldtc, fieldpc)

      fieldtc(:,:) = or1_3D(n_r) * fieldtc(:,:)
      fieldpc(:,:) = or1_3D(n_r) * fieldpc(:,:)

   end subroutine torpol_to_spat
!------------------------------------------------------------------------------
   subroutine torpol_to_curl_spat(Blm, ddBlm, Jlm, dJlm, n_r, &
              &                   curlfieldrc, curlfieldtc, curlfieldpc)
      complex(cp), intent(in) :: Blm(lm_max), ddBlm(lm_max)
      complex(cp), intent(in) :: Jlm(lm_max), dJlm(lm_max)
      integer, intent(in) :: n_r
      real(cp), intent(out) :: curlfieldrc(n_phi_max_3D, n_theta_max)
      real(cp), intent(out) :: curlfieldtc(n_phi_max_3D, n_theta_max)
      real(cp), intent(out) :: curlfieldpc(n_phi_max_3D, n_theta_max)

      !-- Local variable
      complex(cp) :: Qlm(lm_max), Tlm(lm_max)
      integer :: lm

      do lm = 1, lm_max
         !Qlm(lm) = dLh(lm) * or2_3D(n_r) * Jlm(lm)
         Qlm(lm) = dLh(lm) * Jlm(lm)
         Tlm(lm) = or2_3D(n_r) * dLh(lm) * Blm(lm) - ddBlm(lm)
      end do

      call shtns_qst_to_spat(Qlm, dJlm, Tlm, curlfieldrc, curlfieldtc, curlfieldpc)

      !curlfieldtc(:,:) = or1_3D(n_r) * curlfieldtc(:,:)
      !curlfieldpc(:,:) = or1_3D(n_r) * curlfieldpc(:,:)

   end subroutine torpol_to_curl_spat
!------------------------------------------------------------------------------
   subroutine spat_to_SH(f, fLM)

      !-- Input variable
      real(cp), intent(in) :: f(n_phi_max_3D, n_theta_max)

      !-- Output variable
      complex(cp), intent(out) :: fLM(:)

      !-- Local variable
      integer :: nlm

      nlm = size(fLM)

      if ( nlm == lmP_max ) then
         call shtns_load_cfg(1)
      else
         call shtns_load_cfg(0)
      end if
      call shtns_spat_to_sh(f, fLM)
      call shtns_load_cfg(0)

   end subroutine spat_to_SH
!------------------------------------------------------------------------------
   subroutine scal_axi_to_grad_spat(Saxi_l, gradtc)
      ! transform a scalar spherical harmonic field into it's gradient
      ! on the grid

      !-- Input variable
      complex(cp), intent(in) :: Saxi_l(:)

      !-- Output variable
      real(cp), intent(out) :: gradtc(n_theta_max)

      !-- Local variable
      complex(cp) :: tmpt(n_theta_max),tmpp(n_theta_max)

      call shtns_load_cfg(0)
      call shtns_sph_to_spat_ml(0, Saxi_l, tmpt, tmpp, l_max)
      gradtc(:)=real(tmpt)

   end subroutine scal_axi_to_grad_spat
!------------------------------------------------------------------------------
end module shtns
