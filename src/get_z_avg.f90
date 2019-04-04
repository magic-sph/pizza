module z_functions
   !
   ! This module implements the calculation of the z-averaging
   ! Include reconstruction-interpolation of 3D vel from 2D
   !

   use precision_mod
   use parallel_mod
   use mem_alloc
   use parallel_mod, only: n_procs
   use communications, only: help_transp,my_transp_r2all
   use constants, only: zero, half, one, two, pi
   use blocking, only: nRstart, nRstop, nRstart3D, nRstop3D
   use blocking_lm, only: llm, ulm
   use truncation, only: n_r_max, n_m_max, minc, idx2m, m2idx
   use truncation_3D, only: n_r_max_3D, n_z_max, n_m_max_3D, n_theta_max, minc_3D, idx2m3D
   use namelists, only: r_cmb, tadvz_fac, l_ek_pump
   use horizontal, only: cost, sint
   use radial_functions, only: r, r_3D, beta, height

   implicit none

   private

   integer, public, allocatable :: interp_zr_mat(:,:)
   integer, public, allocatable :: interp_zt_mat(:,:)
   real(cp), public, allocatable :: interp_wt_mat(:,:)
   real(cp), public, allocatable :: usr_Rloc(:,:)
   real(cp), public, allocatable :: upp_Rloc(:,:)
   real(cp), public, allocatable :: dVzT_Rloc(:,:)
   real(cp), public, allocatable :: urr_Rloc(:,:,:)
   real(cp), public, allocatable :: uth_Rloc(:,:,:)
   real(cp), public, allocatable :: uph_Rloc(:,:,:)

   public :: initialize_zfunctions, finalize_zfunctions, compute_z_avg, &
   &         get_zavg_mat, extrapolate_z_fields

contains

   subroutine initialize_zfunctions

      allocate( interp_zr_mat(4*n_z_max,n_r_max) )
      allocate( interp_zt_mat(4*n_z_max,n_r_max) )
      allocate( interp_wt_mat(4*n_z_max,n_r_max) )

      interp_zr_mat(:,:)=1
      interp_zt_mat(:,:)=1
      interp_wt_mat(:,:)=0.0_cp

      bytes_allocated = bytes_allocated+8*(n_z_max*n_r_max)*SIZEOF_INTEGER
      bytes_allocated = bytes_allocated+4*(n_z_max*n_r_max)*SIZEOF_DEF_REAL

      allocate( usr_Rloc(n_m_max_3D,nRstart:nRstop) )
      allocate( upp_Rloc(n_m_max_3D,nRstart:nRstop) )
      allocate( dVzT_Rloc(n_m_max_3D,nRstart:nRstop) )
      allocate( urr_Rloc(n_m_max_3D,n_theta_max,nRstart3D:nRstop3D) )
      allocate( uth_Rloc(n_m_max_3D,n_theta_max,nRstart3D:nRstop3D) )
      allocate( uph_Rloc(n_m_max_3D,n_theta_max,nRstart3D:nRstop3D) )

      usr_Rloc(:,:)=0.0_cp
      upp_Rloc(:,:)=0.0_cp
      dVzT_Rloc(:,:)=0.0_cp
      urr_Rloc(:,:,:)=0.0_cp
      uth_Rloc(:,:,:)=0.0_cp
      uph_Rloc(:,:,:)=0.0_cp

      bytes_allocated = bytes_allocated+3*(n_m_max_3D*      &
      &                 (nRstop-nRstart+1))*SIZEOF_DEF_REAL
      bytes_allocated = bytes_allocated+3*(n_m_max_3D*n_theta_max* &
      &                 (nRstop3D-nRstart3D+1))*SIZEOF_DEF_REAL

   end subroutine initialize_zfunctions
!--------------------------------------------------------------------------------
   subroutine finalize_zfunctions

      deallocate( interp_zr_mat, interp_zt_mat, interp_wt_mat )

      deallocate( urr_Rloc, uth_Rloc, uph_Rloc )

   end subroutine finalize_zfunctions
!--------------------------------------------------------------------------------
   subroutine prepare_z_extension(us_Rloc,up_Rloc,usr_Rloc,upp_Rloc,dVzT_Rloc)

      !-- Input variables
      complex(cp), intent(in) :: us_Rloc(n_m_max,nRstart:nRstop)
      complex(cp), intent(in) :: up_Rloc(n_m_max,nRstart:nRstop)

      !-- Output variables
      real(cp), intent(out) :: usr_Rloc(n_m_max_3D,nRstart:nRstop)
      real(cp), intent(out) :: upp_Rloc(n_m_max_3D,nRstart:nRstop)
      real(cp), intent(out) :: dVzT_Rloc(n_m_max_3D,nRstart:nRstop)

      !-- Local variables
      complex(cp) :: usm3D_Rloc(n_m_max_3D,nRstart:nRstop)
      complex(cp) :: upm3D_Rloc(n_m_max_3D,nRstart:nRstop)
      integer :: n_m_3D, n_m, n_r, m3D
      real(cp) :: z, dsZ

      do n_m_3D=1,n_m_max_3D
         m3D = idx2m3D(n_m_3D)
         n_m = m2idx(m3D)
         if( n_m /= -1 ) then
            usm3D_Rloc(n_m_3D,nRstart:nRstop) = us_Rloc(n_m,nRstart:nRstop)
            upm3D_Rloc(n_m_3D,nRstart:nRstop) = us_Rloc(n_m,nRstart:nRstop)
         else
            usm3D_Rloc(n_m_3D,nRstart:nRstop) = zero
            upm3D_Rloc(n_m_3D,nRstart:nRstop) = zero
         end if
      end do

      do n_r=nRstart,nRstop
         call ifft(usm3D_Rloc(:,n_r), usr_Rloc(:,n_r))
         call ifft(upm3D_Rloc(:,n_r), upp_Rloc(:,n_r))
         !---- Get resolution in s and z for z-integral:
         dsZ  =r_cmb/real(n_r_max,cp)  ! Step in s controlled by n_r_max
         z = (n_r-half)*dsZ
         dVzT_Rloc(:,n_r) = beta(n_r)*z*usr_Rloc(:,n_r)
      end do

   end subroutine prepare_z_extension
!--------------------------------------------------------------------------------
   subroutine extrapolate_z_fields(usr_Rloc,upp_Rloc,dVzT_Rloc,urr_Rloc,uth_Rloc,uph_Rloc)

      !-- Input variables
      real(cp), intent(in) :: usr_Rloc(n_m_max_3D,nRstart:nRstop)
      real(cp), intent(in) :: upp_Rloc(n_m_max_3D,nRstart:nRstop)
      real(cp), intent(in) :: dVzT_Rloc(n_m_max_3D,nRstart:nRstop)

      !-- Output variables
      real(cp), intent(out) :: urr_Rloc(n_m_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp), intent(out) :: uth_Rloc(n_m_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp), intent(out) :: uph_Rloc(n_m_max_3D,n_theta_max,nRstart3D:nRstop3D)

      !-- Local arrays
      real(cp) :: usr(n_m_max_3D,n_r_max)
      real(cp) :: upp(n_m_max_3D,n_r_max)
      real(cp) :: dVzT(n_m_max_3D,n_r_max)
      !-- Local variables
      type(help_transp) :: r2all_fields
      integer :: n_r, n_r_r, n_m, n_t_t, n_t_ct
      real(cp) :: s_r, z_r, z_eta
      real(cp) :: vs, vz, vrr, vth, vph
      real(cp) :: alpha_r1, alpha_r2
      logical :: l_non_lin

      l_non_lin=.true.

      call my_transp_r2all(r2all_fields,usr_Rloc,usr)
      call my_transp_r2all(r2all_fields,upp_Rloc,upp)
      call my_transp_r2all(r2all_fields,dVzT_Rloc,dVzT)

      !-- Compute thermal wind by a linear interpolation
      do n_r_r=nRstart3D,nRstop3D
         do n_t_t=1,n_theta_max/2
            n_t_ct=n_theta_max+1-n_t_t
            s_r = r_3D(n_r_r)*sint(n_t_t)
            z_r = r_3D(n_r_r)*cost(n_t_t)
            n_r = 0
            do while ( r(n_r+1) < s_r  )
               n_r = n_r+1
            end do
            alpha_r2 = (s_r-r(n_r))/(r(n_r+1)-r(n_r))
            alpha_r1 = one - alpha_r2
            z_eta = -s_r/(one-s_r**2)*z_r
            do n_m=1,n_m_max_3D
               vs = alpha_r1*usr(n_m,n_r) + alpha_r2*usr(n_m,n_r+1)
               vz = z_eta*vs
               if ( l_ek_pump ) then
                  vz = vz + z_r*(alpha_r1*dVzT(n_m,n_r) +  & 
                  &              alpha_r2*dVzT(n_m,n_r+1))
               end if
               vrr= vz*cost(n_t_t) + vz*sint(n_t_t)
               urr_Rloc(n_m,n_t_t,n_r_r) = vrr
               urr_Rloc(n_m,n_t_ct,n_r_r)= vrr
               if ( l_non_lin ) then
                  vth= vz*cost(n_t_t) - vz*sint(n_t_t)
                  vph= alpha_r1*upp(n_m,n_r) + alpha_r2*upp(n_m,n_r+1)
                  uth_Rloc(n_m,n_t_t,n_r_r) = vth
                  uth_Rloc(n_m,n_t_ct,n_r_r)=-vth
                  uph_Rloc(n_m,n_t_t,n_r_r) = vph
                  uph_Rloc(n_m,n_t_ct,n_r_r)= vph
               end if
            end do
         end do
      end do

      !-- No slip on the sphere surface
      if ( rank==n_procs-1 ) then
         urr_Rloc(:,:,n_r_max_3D) = 0.0_cp
         if ( l_non_lin ) then
            uth_Rloc(:,:,n_r_max_3D) = 0.0_cp
            uph_Rloc(:,:,n_r_max_3D) = 0.0_cp
         end if
      end if

   end subroutine extrapolate_z_fields
!--------------------------------------------------------------------------------
   subroutine compute_z_avg(work_Rloc,zavg_Rloc,adam_Rloc)

      !-- Input variables
      complex(cp), intent(in) :: work_Rloc(n_m_max_3D,n_theta_max,nRstart3D:nRstop3D)

      !-- Output variables
      complex(cp), intent(out) :: zavg_Rloc(n_m_max_3D,nRstart:nRstop)
      complex(cp), intent(out) :: adam_Rloc(n_m_max,nRstart:nRstop)

      !-- Local variables
      integer :: m, n_r, n_m, n_z
      integer :: n_z_m, n_r_r, n_t_t, n_t_ct
      real(cp) :: czavg
      complex(cp) :: tmp(n_m_max,nRstart:nRstop)

      !!-- Transforms to grid for computing the z-avg
      !!-- ==> better to give directly on grid as input, I think!
      !do n_r=nRstart,nRstop
      !    call Spec_Spat_Lorentz_shell(Lorentzp(:,ir),Lorentzt(:,ir),   &
      !                                 Lorentzs(:,ir),Lorentz(:,:,ir))
      !   work_Rloc(:,:,n_r) = Lorentz(:,:,n_r)
      !   call fft_thetab(work_Rloc,+1)
      !enddo

      !!-- need an exchange of points at the boundary domain of each proc 
      !call exch1_2dr3(work_Rloc,n_m_max_3D,n_theta_max,n_r_max_3D,nRStart,nRstop)

      !-- z-average in spacial space
      do n_r=nRstart,nRstop
         do n_z=1,4*n_z_max
            n_t_t = interp_zt_mat(n_z,n_r)
            n_t_ct= n_theta_max+1-n_t_t
            n_r_r = interp_zr_mat(n_z,n_r)
            czavg = interp_wt_mat(n_z,n_r)
            do n_m=1,n_m_max_3D,minc_3D
               !if ( n_r == 1 .or. n_r == n_r_max ) then
               !   zavg_Rloc(n_m,n_r) = 0.0_cp
               !else
                  zavg_Rloc(n_m,n_r) = zavg_Rloc(n_m,n_r) + czavg*  &
                  &                   (work_Rloc(n_m,n_t_t ,n_r_r)  &
                  &                  + work_Rloc(n_m,n_t_ct,n_r_r))
               !endif
            end do
         end do
      end do

      !-- Transforms back to spectral for computing an 'adam' for the velocity field
      do n_r=nRstart,nRstop
         tmp(:,n_r) = zavg_Rloc(:,n_r)
         !call fft_thetab(tmp,-1)
      end do

      !-- Finally get the z-average back in spectral space
      do n_r=nRstart,nRstop
         do n_m=1,n_m_max,minc
            m = idx2m(n_m)
            if ( m /= 0 ) then
                n_z_m = 2*m+1
                adam_Rloc(n_z_m,n_r) = adam_Rloc(n_z_m,n_r) + tmp(n_z_m,n_r)
                adam_Rloc(n_z_m+1,n_r)=adam_Rloc(n_z_m+1,n_r)+tmp(n_z_m+1,n_r)
            else
                adam_Rloc(n_m,n_r) = adam_Rloc(n_m,n_r) + tmp(n_m,n_r)
            end if
         end do
      end do

   end subroutine compute_z_avg
!--------------------------------------------------------------------------------
   subroutine get_zavg_mat

      !-- Local variables
      integer :: n_r, n_z, n_r_r, n_t_t
      real(cp) :: r_r, z_r, c_t, h
      real(cp) :: norm, alpha_r, alpha_t

      !-- Get z grid
      !-- for interpolation: n_r_max points on z-axis
      do n_r=1,n_r_max-1
         h = half*height(n_r)
         n_r_r = n_r_max_3D-1
         n_t_t = 2
         do n_z=n_z_max,0,-1
            z_r = h*n_z/n_z_max
            r_r = sqrt(r(n_r)**2 + z_r**2)
            c_t = z_r/r_r
            do while ( r_3D(n_r_r) > r_r  )
               n_r_r = n_r_r-1
            end do
            do while ( cost(n_t_t) > c_t  )
               n_t_t = n_t_t+1
            end do
            !-- Compute coeffs for bilinear interpolation
            !if ( n_z /= 0 ) then
            !   !-- Multiply by 2 because of the symmetries
            !   norm = two/(two*n_z_max+1)
            !else
               norm = one/(two*n_z_max+1)
            !endif
            alpha_r = (r_r-r_3D(n_r_r))/(r_3D(n_r_r+1)-r_3D(n_r_r))
            alpha_t = (c_t-cost(n_t_t))/(cost(n_t_t-1)-cost(n_t_t))
            !-- Coordinates neighbourg 1: n_r_r, n_t_t
            interp_zr_mat(4*n_z-3,n_r)=n_r_r
            interp_zt_mat(4*n_z-3,n_r)=n_t_t
            interp_wt_mat(4*n_z-3,n_r)=norm*(one-alpha_r)*(one-alpha_t)
            !-- Coordinates neighbourg 2: n_r_r, n_t_t-1
            interp_zr_mat(4*n_z-2,n_r)=n_r_r
            interp_zt_mat(4*n_z-2,n_r)=n_t_t-1
            interp_wt_mat(4*n_z-2,n_r)=norm*(one-alpha_r)*alpha_t
            !-- Coordinates neighbourg 3: n_r_r+1, n_t_t
            interp_zr_mat(4*n_z-1,n_r)=n_r_r+1
            interp_zt_mat(4*n_z-1,n_r)=n_t_t
            interp_wt_mat(4*n_z-1,n_r)=norm*alpha_r*(one-alpha_t)
            !-- Coordinates neighbourg 4: n_r_r+1, n_t_t-1
            interp_zr_mat(4*n_z,n_r)  =n_r_r+1
            interp_zt_mat(4*n_z,n_r)  =n_t_t-1
            interp_wt_mat(4*n_z,n_r)  =norm*alpha_r*alpha_t
         end do
      end do

   end subroutine get_zavg_mat
!--------------------------------------------------------------------------------
end module z_functions
