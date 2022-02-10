module shtns
  !
   ! This module contains is a wrapper of the SHTns routines used in pizza
   !

   use iso_c_binding
   use iso_fortran_env, only: output_unit
   use precision_mod, only: cp
   use constants, only: ci
   use truncation_3D, only: m_max_3D, l_max, n_theta_max, n_phi_max_3D, &
       &                    minc_3D, lm_max, lmP_max
   use horizontal, only: dLh
   use radial_functions, only: or2_3D, or1_3D
   use parallel_mod

   implicit none

   include "shtns.f03"

   private

   public :: init_shtns, scal_to_spat, scal_to_SH, spat_to_qst, &
   &         torpol_to_spat,torpol_to_curl_spat,                &
   &         scal_axi_to_grad_spat

   type(c_ptr) :: sht_l, sht_lP

contains

   subroutine init_shtns()

      integer :: norm, nthreads, layout
      real(cp) :: eps_polar
      type(shtns_info), pointer :: sht_info

      if ( rank == 0 ) then
         write(output_unit,*) ''
         call shtns_verbose(1)
      end if

      nthreads = shtns_use_threads(0)

      norm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE
      !layout = SHT_QUICK_INIT+SHT_PHI_CONTIGUOUS
      layout = SHT_GAUSS+SHT_PHI_CONTIGUOUS
      eps_polar = 1.e-10_cp

      sht_l = shtns_create(l_max, m_max_3D/minc_3D, minc_3D, norm)
      call shtns_set_grid(sht_l, layout, eps_polar, n_theta_max, n_phi_max_3D)

      call c_f_pointer(cptr=sht_l, fptr=sht_info)

      if ( rank == 0 ) then
         call shtns_verbose(0)
      end if

      sht_lP = shtns_create(l_max+1, m_max_3D/minc_3D, minc_3D, norm)
      call shtns_set_grid(sht_lP, layout, eps_polar, n_theta_max, n_phi_max_3D)

   end subroutine
!------------------------------------------------------------------------------
   subroutine scal_to_spat(Slm, fieldc)

      ! transform a spherical harmonic field into grid space
      complex(cp), intent(inout) :: Slm(lm_max)
      real(cp), intent(out) :: fieldc(n_phi_max_3D, n_theta_max)

      call SH_to_spat(sht_l, Slm, fieldc)

   end subroutine scal_to_spat
!------------------------------------------------------------------------------
   subroutine torpol_to_spat(Wlm, dWlm, Zlm, n_r, fieldrc, fieldtc, fieldpc)

      !-- Input variables
      complex(cp), intent(inout) :: Wlm(lm_max), dWlm(lm_max), Zlm(lm_max)
      integer,     intent(in) :: n_r

      !-- Output variables
      real(cp), intent(out) :: fieldrc(n_phi_max_3D, n_theta_max)
      real(cp), intent(out) :: fieldtc(n_phi_max_3D, n_theta_max)
      real(cp), intent(out) :: fieldpc(n_phi_max_3D, n_theta_max)

      !-- Local variables
      complex(cp) :: Qlm(lm_max)
      integer :: lm

      do lm = 1, lm_max
         Qlm(lm) = dLh(lm) * or2_3D(n_r) * Wlm(lm)
      end do

      call SHqst_to_spat(sht_l, Qlm, dWlm, Zlm, fieldrc, fieldtc, fieldpc)

      fieldtc(:,:) = or1_3D(n_r) * fieldtc(:,:)
      fieldpc(:,:) = or1_3D(n_r) * fieldpc(:,:)

   end subroutine torpol_to_spat
!------------------------------------------------------------------------------
   subroutine torpol_to_curl_spat(Blm, ddBlm, Jlm, dJlm, n_r, &
              &                   curlfieldrc, curlfieldtc, curlfieldpc)

      !-- Input variables
      complex(cp), intent(inout) :: Blm(lm_max), ddBlm(lm_max)
      complex(cp), intent(inout) :: Jlm(lm_max), dJlm(lm_max)
      integer,     intent(in) :: n_r

      !-- Output variables
      real(cp), intent(out) :: curlfieldrc(n_phi_max_3D, n_theta_max)
      real(cp), intent(out) :: curlfieldtc(n_phi_max_3D, n_theta_max)
      real(cp), intent(out) :: curlfieldpc(n_phi_max_3D, n_theta_max)

      !-- Local variables
      complex(cp) :: Qlm(lm_max), Tlm(lm_max)
      integer :: lm

      do lm = 1, lm_max
         Qlm(lm) = dLh(lm) * or2_3D(n_r) * Jlm(lm)
         Tlm(lm) = or2_3D(n_r) * dLh(lm) * Blm(lm) - ddBlm(lm)
      end do

      call SHqst_to_spat(sht_l, Qlm, dJlm, Tlm, curlfieldrc, curlfieldtc, curlfieldpc)

      curlfieldtc(:,:) = or1_3D(n_r) * curlfieldtc(:,:)
      curlfieldpc(:,:) = or1_3D(n_r) * curlfieldpc(:,:)

   end subroutine torpol_to_curl_spat
!------------------------------------------------------------------------------
   subroutine spat_to_qst(f, g, h, fLM, gLM, hLM)

      !-- Input variables
      real(cp), intent(inout) :: f(:,:)
      real(cp), intent(inout) :: g(:,:)
      real(cp), intent(inout) :: h(:,:)

      !-- Output variables
      complex(cp), intent(out) :: fLM(:)
      complex(cp), intent(out) :: gLM(:)
      complex(cp), intent(out) :: hLM(:)

      call spat_to_SHqst(sht_l, f, g, h, fLM, gLM, hLM)

   end subroutine spat_to_qst
!------------------------------------------------------------------------------
   subroutine scal_to_SH(f, fLM)

      !-- Input variable
      real(cp), intent(inout) :: f(:,:)

      !-- Output variable
      complex(cp), intent(out) :: fLM(:)

      !-- Local variable
      integer :: nlm

      nlm = size(fLM)

      if ( nlm == lmP_max ) then
         call spat_to_SH(sht_lP, f, fLM)
      else
         call spat_to_SH(sht_l, f, fLM)
      end if

   end subroutine scal_to_SH
!------------------------------------------------------------------------------
   subroutine scal_axi_to_grad_spat(Saxi_l, gradtc)
      ! transform a scalar spherical harmonic field into it's gradient
      ! on the grid

      !-- Input variable
      complex(cp), intent(inout) :: Saxi_l(:)

      !-- Output variable
      real(cp), intent(out) :: gradtc(n_theta_max)

      !-- Local variable
      complex(cp) :: tmpt(n_theta_max),tmpp(n_theta_max)

      call SHsph_to_spat_ml(sht_l, 0, Saxi_l, tmpt, tmpp, l_max)
      gradtc(:)=real(tmpt)

   end subroutine scal_axi_to_grad_spat
!------------------------------------------------------------------------------
end module shtns
