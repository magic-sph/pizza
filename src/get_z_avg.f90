module z_functions
   !
   ! This module implements the calculation of the z-averaging
   ! Include reconstruction-interpolation of 3D vel from 2D
   !

   use precision_mod
   use parallel_mod
   use mem_alloc
   use fourier, only: ifft, fft
   use communications, only: allgather_from_rloc
   use constants, only: zero, half, one, two
   use blocking, only: nRstart, nRstop, nRstart3D, nRstop3D
   use truncation, only: n_r_max, n_m_max, minc, idx2m, m2idx
   use truncation_3D, only: n_r_max_3D, n_z_max, n_m_max_3D, n_theta_max, &
       &                    minc_3D, idx2m3D, n_phi_max_3D
   use namelists, only: r_cmb, l_ek_pump
   use horizontal, only: cost, sint
   use radial_functions, only: r, r_3D, beta, height

   implicit none

   private

   type, public :: zfunc_type

      integer, allocatable :: interp_zr_mat(:,:)
      integer, allocatable :: interp_zt_mat(:,:)
      real(cp), allocatable :: interp_wt_mat(:,:)
      real(cp), allocatable :: us_phys_Rloc(:,:)
      real(cp), allocatable :: up_phys_Rloc(:,:)
      real(cp), allocatable :: dVzT_Rloc(:,:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: compute_avg
      procedure :: prepare_extension
      procedure :: extrapolate
      procedure :: fill_mat
   end type zfunc_type

contains

   subroutine initialize(this)

      class(zfunc_type) :: this

      allocate( this%interp_zr_mat(4*n_z_max,n_r_max) )
      allocate( this%interp_zt_mat(4*n_z_max,n_r_max) )
      allocate( this%interp_wt_mat(4*n_z_max,n_r_max) )

      this%interp_zr_mat(:,:)=1
      this%interp_zt_mat(:,:)=1
      this%interp_wt_mat(:,:)=0.0_cp

      bytes_allocated = bytes_allocated+8*(n_z_max*n_r_max)*SIZEOF_INTEGER
      bytes_allocated = bytes_allocated+4*(n_z_max*n_r_max)*SIZEOF_DEF_REAL

      allocate( this%us_phys_Rloc(n_phi_max_3D,nRstart:nRstop) )
      allocate( this%up_phys_Rloc(n_phi_max_3D,nRstart:nRstop) )
      allocate( this%dVzT_Rloc(n_phi_max_3D,nRstart:nRstop) )

      this%us_phys_Rloc(:,:)=0.0_cp
      this%up_phys_Rloc(:,:)=0.0_cp
      this%dVzT_Rloc(:,:)   =0.0_cp

      bytes_allocated = bytes_allocated+3*(n_phi_max_3D*      &
      &                 (nRstop-nRstart+1))*SIZEOF_DEF_REAL

   end subroutine initialize
!--------------------------------------------------------------------------------
   subroutine finalize(this)

      class(zfunc_type) :: this

      deallocate( this%interp_zr_mat, this%interp_zt_mat, this%interp_wt_mat )
      deallocate( this%us_phys_Rloc, this%up_phys_Rloc, this%dVzT_Rloc )

   end subroutine finalize
!--------------------------------------------------------------------------------
   subroutine prepare_extension(this, us_Rloc,up_Rloc)

      class(zfunc_type) :: this

      !-- Input variables
      complex(cp), intent(in) :: us_Rloc(n_m_max,nRstart:nRstop)
      complex(cp), intent(in) :: up_Rloc(n_m_max,nRstart:nRstop)

      !-- Local variables
      complex(cp) :: usm3D_Rloc(n_m_max_3D,nRstart:nRstop)
      complex(cp) :: upm3D_Rloc(n_m_max_3D,nRstart:nRstop)
      integer :: n_m_3D, n_m, n_r, m3D
      real(cp) :: z, dsZ

      do n_m_3D=1,n_m_max_3D
         m3D = idx2m3D(n_m_3D)
         if ( m3D < size(m2idx) ) then ! a little bit weird
            n_m = m2idx(m3D)
         else
            n_m = -1
         end if
         if( n_m /= -1 ) then
            usm3D_Rloc(n_m_3D,nRstart:nRstop) = us_Rloc(n_m,nRstart:nRstop)
            upm3D_Rloc(n_m_3D,nRstart:nRstop) = up_Rloc(n_m,nRstart:nRstop)
         else
            usm3D_Rloc(n_m_3D,nRstart:nRstop) = zero
            upm3D_Rloc(n_m_3D,nRstart:nRstop) = zero
         end if
      end do

      do n_r=nRstart,nRstop
         call ifft(usm3D_Rloc(:,n_r), this%us_phys_Rloc(:,n_r), l_3D=.true.)
         call ifft(upm3D_Rloc(:,n_r), this%up_phys_Rloc(:,n_r), l_3D=.true.)
         !---- Get resolution in s and z for z-integral:
         dsZ  =r_cmb/real(n_r_max,cp)  ! Step in s controlled by n_r_max
         z = (n_r-half)*dsZ
         this%dVzT_Rloc(:,n_r) = beta(n_r)*z*this%us_phys_Rloc(:,n_r)
      end do

   end subroutine prepare_extension
!--------------------------------------------------------------------------------
   subroutine extrapolate(this, ur_Rloc, ut_Rloc, up_Rloc)

      class(zfunc_type) :: this

      !-- Output variables
      real(cp), intent(out) :: ur_Rloc(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp), intent(out) :: ut_Rloc(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp), intent(out) :: up_Rloc(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)

      !-- Local arrays
      real(cp) :: usr(n_phi_max_3D,n_r_max)
      real(cp) :: upp(n_phi_max_3D,n_r_max)
      real(cp) :: dVzT(n_phi_max_3D,n_r_max)
      !-- Local variables
      integer :: n_r, n_r_r, n_phi, n_t_t, n_t_ct
      real(cp) :: s_r, z_r, z_eta
      real(cp) :: vs, vz, vrr, vth, vph
      real(cp) :: alpha_r1, alpha_r2

      call allgather_from_rloc(this%us_phys_Rloc,usr,n_phi_max_3D)
      call allgather_from_rloc(this%up_phys_Rloc,upp,n_phi_max_3D)
      call allgather_from_rloc(this%dVzT_Rloc,dVzT,n_phi_max_3D)

      !-- Compute thermal wind by a linear interpolation
      do n_r_r=nRstart3D,nRstop3D
         do n_t_t=1,n_theta_max/2
            n_t_ct=n_theta_max+1-n_t_t
            s_r = r_3D(n_r_r)*sint(n_t_t)
            z_r = r_3D(n_r_r)*cost(n_t_t)
            n_r = n_r_max_3D+1
            do while ( r(n_r-1) < s_r  )
               n_r = n_r-1
            end do
            alpha_r2 = (s_r-r(n_r))/(r(n_r-1)-r(n_r))
            alpha_r1 = one - alpha_r2
            z_eta = -s_r/(one-s_r**2)*z_r
            do n_phi=1,n_phi_max_3D
               vs = alpha_r1*usr(n_phi,n_r) + alpha_r2*usr(n_phi,n_r-1)
               vz = z_eta*vs
               if ( l_ek_pump ) then
                  vz = vz + z_r*(alpha_r1*dVzT(n_phi,n_r) +  & 
                  &              alpha_r2*dVzT(n_phi,n_r-1))
               end if
               vrr= vz*cost(n_t_t) + vz*sint(n_t_t)
               ur_Rloc(n_phi,n_t_t,n_r_r) = vrr
               ur_Rloc(n_phi,n_t_ct,n_r_r)= vrr
               vth= vz*cost(n_t_t) - vz*sint(n_t_t)
               vph= alpha_r1*upp(n_phi,n_r) + alpha_r2*upp(n_phi,n_r-1)
               ut_Rloc(n_phi,n_t_t,n_r_r) = vth
               ut_Rloc(n_phi,n_t_ct,n_r_r)=-vth
               up_Rloc(n_phi,n_t_t,n_r_r) = vph
               up_Rloc(n_phi,n_t_ct,n_r_r)= vph
            end do
         end do
      end do

      !-- No slip on the spherical surface
      if ( nRstart3D == 1 ) then
         ur_Rloc(:,:,1) = 0.0_cp
         ut_Rloc(:,:,1) = 0.0_cp
         up_Rloc(:,:,1) = 0.0_cp
      end if

   end subroutine extrapolate
!--------------------------------------------------------------------------------
   subroutine compute_avg(this,work_Rloc,zavg_Rloc)

      class(zfunc_type) :: this

      !-- Input variables
      real(cp), intent(in) :: work_Rloc(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)

      !-- Output variables
      complex(cp), intent(out) :: zavg_Rloc(n_m_max,nRstart:nRstop)

      !-- Local variables
      integer :: n_r, n_phi, n_z, n_m_QG, n_m, m3D
      integer :: n_r_r, n_t_t, n_t_ct
      complex(cp) :: tmp_hat(n_m_max_3D)
      real(cp) :: czavg
      real(cp) :: tmp(n_phi_max_3D,n_r_max)

      tmp(:,:)=0.0_cp
      zavg_Rloc(:,:)=zero

      !-- z-average in spatial space
      do n_r=1,n_r_max
         do n_z=1,4*n_z_max
            n_t_t = this%interp_zt_mat(n_z,n_r)
            n_t_ct= n_theta_max+1-n_t_t
            n_r_r = this%interp_zr_mat(n_z,n_r)
            czavg = this%interp_wt_mat(n_z,n_r)
            if ( n_r_r >= nRstart3D .and. n_r_r <= nRstop3D ) then
               do n_phi=1,n_phi_max_3D
                  tmp(n_phi,n_r) = tmp(n_phi,n_r) + czavg*             &
                  &                   (work_Rloc(n_phi,n_t_t ,n_r_r)   &
                  &                  + work_Rloc(n_phi,n_t_ct,n_r_r))
               end do
            end if
         end do
      end do

      call MPI_Allreduce(MPI_IN_PLACE, tmp, n_phi_max_3D, MPI_DEF_REAL, &
           &             MPI_SUM, MPI_COMM_WORLD, ierr)


      !-- Transforms back to spectral space
      do n_r=nRstart,nRstop
         call fft(tmp(:,n_r), tmp_hat(:), l_3D=.true.)
         do n_m=1,n_m_max_3D
            m3D = idx2m3D(n_m)
            if ( m3D < size(m2idx) ) then ! a little bit weird
               n_m_QG = m2idx(m3D)
            else
               n_m_QG = -1
            end if
            if ( n_m_QG > 0 ) then
               zavg_Rloc(n_m_QG,n_r)=tmp_hat(n_m)
            end if
         end do
      end do

   end subroutine compute_avg
!--------------------------------------------------------------------------------
   subroutine fill_mat(this)

      class(zfunc_type) :: this

      !-- Local variables
      integer :: n_r, n_z, n_r_r, n_t_t
      real(cp) :: r_r, z_r, c_t, h
      real(cp) :: norm, alpha_r, alpha_t

      !-- Get z grid
      !-- for interpolation: n_r_max points on z-axis
      do n_r=1,n_r_max-1
         h = half*height(n_r)
         !n_r_r = n_r_max_3D-1
         n_r_r = 2
         n_t_t = 2
         do n_z=n_z_max,4,-1
            z_r = h*n_z/n_z_max
            r_r = sqrt(r(n_r)**2 + z_r**2)
            c_t = z_r/r_r
            do while ( r_3D(n_r_r) > r_r  )
               n_r_r = n_r_r+1
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
            alpha_r = (r_r-r_3D(n_r_r))/(r_3D(n_r_r-1)-r_3D(n_r_r))
            alpha_t = (c_t-cost(n_t_t))/(cost(n_t_t-1)-cost(n_t_t))
            !-- Coordinates neighbourg 1: n_r_r, n_t_t
            this%interp_zr_mat(4*n_z-3,n_r)=n_r_r
            this%interp_zt_mat(4*n_z-3,n_r)=n_t_t
            this%interp_wt_mat(4*n_z-3,n_r)=norm*(one-alpha_r)*(one-alpha_t)
            !-- Coordinates neighbourg 2: n_r_r, n_t_t-1
            this%interp_zr_mat(4*n_z-2,n_r)=n_r_r
            this%interp_zt_mat(4*n_z-2,n_r)=n_t_t-1
            this%interp_wt_mat(4*n_z-2,n_r)=norm*(one-alpha_r)*alpha_t
            !-- Coordinates neighbourg 3: n_r_r+1, n_t_t
            this%interp_zr_mat(4*n_z-1,n_r)=n_r_r+1
            this%interp_zt_mat(4*n_z-1,n_r)=n_t_t
            this%interp_wt_mat(4*n_z-1,n_r)=norm*alpha_r*(one-alpha_t)
            !-- Coordinates neighbourg 4: n_r_r+1, n_t_t-1
            this%interp_zr_mat(4*n_z,n_r)  =n_r_r+1
            this%interp_zt_mat(4*n_z,n_r)  =n_t_t-1
            this%interp_wt_mat(4*n_z,n_r)  =norm*alpha_r*alpha_t
         end do
      end do

   end subroutine fill_mat
!--------------------------------------------------------------------------------
end module z_functions
