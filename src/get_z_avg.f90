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
   use constants, only: zero, half, one, two, ci, pi
   use blocking, only: nRstart, nRstop, nRstart3D, nRstop3D
   use truncation, only: n_r_max, n_m_max, minc, idx2m, m2idx
   use truncation_3D, only: n_r_max_3D, n_z_max, n_m_max_3D, n_theta_max, &
       &                    minc_3D, idx2m3D, n_phi_max_3D
   use namelists, only: r_icb, r_cmb, l_ek_pump, ktopv, CorFac, ek, ra, &
       &                BuoFac
   use horizontal, only: theta, cost, sint
   use radial_functions, only: r, r_3D, beta, height, oheight, ekpump, or1_3D, &
       &                       rgrav_3D

   implicit none

   private

   type, public :: zfunc_type

      integer, allocatable :: nzp_thw(:,:)
      integer, allocatable :: interp_zr_mat(:,:)
      integer, allocatable :: interp_zt_mat(:,:)
      integer, allocatable :: interp_zp_thw(:,:,:,:)
      integer, allocatable :: interp_zpb_thw(:,:,:)
      real(cp), allocatable :: interp_wt_mat(:,:)
      real(cp), allocatable :: interp_wt_thw(:,:,:,:)
      real(cp), allocatable :: interp_wtb_thw(:,:,:,:)
      real(cp), allocatable :: us_phys_Rloc(:,:)
      real(cp), allocatable :: up_phys_Rloc(:,:)
      real(cp), allocatable :: ek_phys_Rloc(:,:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: compute_avg
      procedure :: prepare_extension
      procedure :: extrapolate
      procedure :: compute_thermal_wind
      procedure :: fill_mat
   end type zfunc_type

contains

   subroutine initialize(this)

      class(zfunc_type) :: this

      !-- Local variable
      integer :: n_size

      allocate( this%interp_zr_mat(4*n_z_max,n_r_max) )
      allocate( this%interp_zt_mat(4*n_z_max,n_r_max) )
      allocate( this%interp_wt_mat(4*n_z_max,n_r_max) )

      this%interp_zr_mat(:,:)=1
      this%interp_zt_mat(:,:)=1
      this%interp_wt_mat(:,:)=0.0_cp

      bytes_allocated = bytes_allocated+8*(n_z_max*n_r_max)*SIZEOF_INTEGER
      bytes_allocated = bytes_allocated+4*(n_z_max*n_r_max)*SIZEOF_DEF_REAL

      n_size=nRstop3D-nRstart3D+1
      allocate( this%nzp_thw(n_theta_max/2,nRstart3D:nRstop3D) )
      allocate( this%interp_zp_thw(n_size,n_theta_max/2,nRstart3D:nRstop3D,2) )
      allocate( this%interp_zpb_thw(n_procs,n_theta_max/2,nRstart3D:nRstop3D) )
      allocate( this%interp_wt_thw(n_size,n_theta_max/2,nRstart3D:nRstop3D,2) )
      allocate( this%interp_wtb_thw(n_procs,n_theta_max/2,nRstart3D:nRstop3D,2) )

      this%nzp_thw(:,:)=1
      this%interp_zp_thw(:,:,:,:)=1
      this%interp_zpb_thw(:,:,:)=1
      this%interp_wt_thw(:,:,:,:)=0.0_cp
      this%interp_wtb_thw(:,:,:,:)=0.0_cp

      bytes_allocated = bytes_allocated+((2*n_size+n_procs+1)*n_theta_max/2* &
      &                 (nRstop3D-nRstart3D+1))*SIZEOF_INTEGER
      bytes_allocated = bytes_allocated+2*((n_size+n_procs)*n_theta_max/2* &
      &                 (nRstop3D-nRstart3D+1))*SIZEOF_DEF_REAL

      allocate( this%us_phys_Rloc(n_phi_max_3D,nRstart:nRstop) )
      allocate( this%up_phys_Rloc(n_phi_max_3D,nRstart:nRstop) )
      allocate( this%ek_phys_Rloc(n_phi_max_3D,nRstart:nRstop) )

      this%us_phys_Rloc(:,:)=0.0_cp
      this%up_phys_Rloc(:,:)=0.0_cp
      this%ek_phys_Rloc(:,:)=0.0_cp

      bytes_allocated = bytes_allocated+3*(n_phi_max_3D*      &
      &                 (nRstop-nRstart+1))*SIZEOF_DEF_REAL

   end subroutine initialize
!--------------------------------------------------------------------------------
   subroutine finalize(this)

      class(zfunc_type) :: this

      deallocate( this%interp_zr_mat, this%interp_zt_mat, this%interp_wt_mat )
      deallocate( this%nzp_thw, this%interp_zp_thw, this%interp_wt_thw )
      deallocate( this%interp_zpb_thw, this%interp_wtb_thw )
      deallocate( this%us_phys_Rloc, this%up_phys_Rloc, this%ek_phys_Rloc )

   end subroutine finalize
!--------------------------------------------------------------------------------
   subroutine prepare_extension(this, us_Rloc, up_Rloc, om_Rloc)

      class(zfunc_type) :: this

      !-- Input variables
      complex(cp), intent(in) :: us_Rloc(n_m_max,nRstart:nRstop)
      complex(cp), intent(in) :: up_Rloc(n_m_max,nRstart:nRstop)
      complex(cp), intent(in) :: om_Rloc(n_m_max,nRstart:nRstop)

      !-- Local variables
      complex(cp) :: usm3D_Rloc(n_m_max_3D,nRstart:nRstop)
      complex(cp) :: upm3D_Rloc(n_m_max_3D,nRstart:nRstop)
      complex(cp) :: ekpump_m3D(n_m_max_3D,nRstart:nRstop)
      integer :: n_m_3D, n_m, n_r, m3D

      do n_r=nRstart,nRstop
         do n_m_3D=1,n_m_max_3D
            m3D = idx2m3D(n_m_3D)
            if ( m3D < size(m2idx) ) then ! a little bit weird
               n_m = m2idx(m3D)
            else
               n_m = -1
            end if
            if ( n_m /= -1 ) then
               usm3D_Rloc(n_m_3D,n_r) = us_Rloc(n_m,n_r)
               upm3D_Rloc(n_m_3D,n_r) = up_Rloc(n_m,n_r)
               ekpump_m3D(n_m_3D,n_r) = zero
               if ( l_ek_pump ) then
                  if( m3D /= 0 ) then
                     ekpump_m3D(n_m_3D,n_r)=ekpump(n_r)*(-om_Rloc(n_m,n_r) &
                     &              +half*beta(n_r)*      up_Rloc(n_m,n_r) &
                     &   +beta(n_r)*( -ci*real(m3D,cp)+                    &
                     &    5.0_cp*r_cmb*oheight(n_r) )*    us_Rloc(n_m,n_r) )
                  else
                     ekpump_m3D(n_m_3D,n_r)=0.0_cp
                     !ekpump_m3D(n_m_3D,n_r)=-CorFac*                &
                     !&                                 ekpump(n_r)* &
                     !&                                 up_Rloc(n_m,n_r)
                  end if
               end if
            else
               usm3D_Rloc(n_m_3D,n_r) = zero
               upm3D_Rloc(n_m_3D,n_r) = zero
               ekpump_m3D(n_m_3D,n_r) = zero
            end if
         end do
      end do

      do n_r=nRstart,nRstop
         call ifft(usm3D_Rloc(:,n_r), this%us_phys_Rloc(:,n_r), l_3D=.true.)
         call ifft(upm3D_Rloc(:,n_r), this%up_phys_Rloc(:,n_r), l_3D=.true.)
         if( l_ek_pump ) &
         &   call ifft(ekpump_m3D(:,n_r), this%ek_phys_Rloc(:,n_r), l_3D=.true.)
      end do

      !-- Boundary point: fix Ek-pumping to zero
      if ( rank == 0 ) this%ek_phys_Rloc(:,1)=0.0_cp

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
      real(cp) :: ekp(n_phi_max_3D,n_r_max)
      !-- Local variables
      integer :: n_r, n_r_r, n_phi, n_t_t, n_t_ct
      real(cp) :: s_r, z_r, z_eta
      real(cp) :: vs, vz, vrr, vth, vph
      real(cp) :: alpha_r1, alpha_r2

      call allgather_from_rloc(this%us_phys_Rloc,usr,n_phi_max_3D)
      call allgather_from_rloc(this%up_phys_Rloc,upp,n_phi_max_3D)
      if( l_ek_pump ) call allgather_from_rloc(this%ek_phys_Rloc,ekp,n_phi_max_3D)

      !-- Compute 3D velocity fields by a linear interpolation
      do n_r_r=nRstart3D,nRstop3D
         do n_t_t=1,n_theta_max/2
            n_t_ct=n_theta_max+1-n_t_t
            s_r = r_3D(n_r_r)*sint(n_t_t)
            z_r = r_3D(n_r_r)*cost(n_t_t)
            if ( s_r >= r_icb ) then !-- Outside TC

               n_r = 1
               do while ( r(n_r) > s_r .and. n_r < n_r_max )
                  n_r = n_r+1
               end do
               alpha_r2 = (s_r-r(n_r))/(r(n_r-1)-r(n_r))
               alpha_r1 = one - alpha_r2
               z_eta = -s_r/(r_cmb*r_cmb-s_r*s_r)*z_r ! \beta * z
               do n_phi=1,n_phi_max_3D
                  vs = alpha_r1*usr(n_phi,n_r) + alpha_r2*usr(n_phi,n_r-1)
                  !-- vz = beta*vs
                  vz = z_eta*vs
                  if ( l_ek_pump ) then
                     !-- vz = beta*vz+ekpump
                     vz = vz + z_r*(alpha_r1*ekp(n_phi,n_r) +  & 
                     &              alpha_r2*ekp(n_phi,n_r-1))
                  end if
                  vrr= vz*cost(n_t_t) + vs*sint(n_t_t)
                  ur_Rloc(n_phi,n_t_t,n_r_r) = vrr
                  ur_Rloc(n_phi,n_t_ct,n_r_r)= vrr
                  vth= vs*cost(n_t_t) - vz*sint(n_t_t)
                  vph= alpha_r1*upp(n_phi,n_r) + alpha_r2*upp(n_phi,n_r-1)
                  ut_Rloc(n_phi,n_t_t,n_r_r) = vth
                  ut_Rloc(n_phi,n_t_ct,n_r_r)=-vth
                  up_Rloc(n_phi,n_t_t,n_r_r) = vph
                  up_Rloc(n_phi,n_t_ct,n_r_r)= vph
               end do

            else !-- Inside the tangent cylinder

               do n_phi=1,n_phi_max_3D
                  ur_Rloc(n_phi,n_t_t,n_r_r) =0.0_cp
                  ur_Rloc(n_phi,n_t_ct,n_r_r)=0.0_cp
                  ut_Rloc(n_phi,n_t_t,n_r_r) =0.0_cp
                  ut_Rloc(n_phi,n_t_ct,n_r_r)=0.0_cp
                  up_Rloc(n_phi,n_t_t,n_r_r) =0.0_cp
                  up_Rloc(n_phi,n_t_ct,n_r_r)=0.0_cp
               end do

            end if ! Inside/outside TC
         end do
      end do

      !-- CMB values
      if ( nRstart3D == 1 ) then
         ur_Rloc(:,:,1) = 0.0_cp
         if ( ktopv == 2 ) then !-- Rigid boundaries
            ut_Rloc(:,:,1) = 0.0_cp
            up_Rloc(:,:,1) = 0.0_cp
         end if
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

      !-- TG: dirty fix: needs to be improved
      do n_r=1,n_r_max
         call MPI_Allreduce(MPI_IN_PLACE, tmp(:,n_r), n_phi_max_3D, MPI_DEF_REAL, &
              &             MPI_SUM, MPI_COMM_WORLD, ierr)
      !call MPI_reduce(MPI_IN_PLACE, tmp(:,n_r), n_phi_max_3D, MPI_DEF_REAL, &
      !     &             MPI_SUM, MPI_COMM_WORLD, ierr)
      end do

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
   subroutine compute_thermal_wind(this, dTdth, up_Rloc)

      !-- Input variables
      class(zfunc_type) :: this
      real(cp), intent(in) :: dTdth(n_theta_max,nRstart3D:nRstop3D)

      !-- Output variables - modified (inout)
      real(cp), intent(inout) :: up_Rloc(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)

      !-- Local arrays
      real(cp) :: thw_Rloc(n_theta_max/2,nRstart3D:nRstop3D)
      real(cp) :: dTzdt(n_theta_max/2,nRstart3D:nRstop3D)
      real(cp) :: tmp(n_theta_max/2,n_procs)
      !-- Local variables
      integer :: n_theta, n_z, n_z_r, n_z_t, n_t_cth
      integer :: n_r, n_rank
      real(cp) :: thwFac

      thwFac=BuoFac/CorFac

      !-- Remaining term for the temperature gradient
      dTzdt(:,:)=0.0_cp
      do n_r=max(2,nRstart3D),nRstop3D
         do n_theta=1,n_theta_max/2
            !dTzdt(n_theta,n_r)=thwFac*r_3D(n_r)* dTdth(n_theta,n_r)
            !-- TG I don't understand the r factor in the above equation
            !-- Th wind should be
            !-- duphi/dz = Ra/Pr  * g / r * dT/dtheta
            dTzdt(n_theta,n_r)=thwFac*rgrav_3D(n_r)*or1_3D(n_r)* dTdth(n_theta,n_r)
         end do
      end do

      !-- Compute thermal wind
      thw_Rloc(:,:)=0.0_cp
      do n_r=nRstart3D,nRstop3D
         n_theta=1
         do n_z=1,this%nzp_thw(n_theta,n_r)
            n_z_r = this%interp_zp_thw(n_z,n_theta,n_r,1)
            n_z_t = this%interp_zp_thw(n_z,n_theta,n_r,2)
            thw_Rloc(n_theta,n_r)=thw_Rloc(n_theta,n_r) -                &
            &                     this%interp_wt_thw(n_z,n_theta,n_r,1)* &
            &                     dTzdt(n_z_t,n_z_r)
         end do
         do n_theta=2,n_theta_max/2
            do n_z=1,this%nzp_thw(n_theta,n_r)
               n_z_r = this%interp_zp_thw(n_z,n_theta,n_r,1)
               n_z_t = this%interp_zp_thw(n_z,n_theta,n_r,2)
               if ( n_z_t > 1 ) then
                  thw_Rloc(n_theta,n_r)= thw_Rloc(n_theta,n_r) -                &
                  &                     (this%interp_wt_thw(n_z,n_theta,n_r,1)* &
                  & dTzdt(n_z_t,n_z_r) + this%interp_wt_thw(n_z,n_theta,n_r,2)* &
                  & dTzdt(n_z_t-1,n_z_r))
               end if
            end do
         end do
      end do
      if(n_procs>1) then
         tmp(:,:)=0.0_cp
         do n_rank=1,n_procs-1
            !print*, 'before BCast, n_rank, tmp, target', n_rank, tmp(1:2,n_rank), thw_Rloc(1:2,nRstart3D)
            if( n_rank == rank ) tmp(:,n_rank) = thw_Rloc(:,nRstart3D)
            call MPI_Bcast(tmp(:,n_rank), n_theta_max/2, MPI_DEF_REAL, &
            &              n_rank, MPI_COMM_WORLD, ierr)
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            !print*, 'after BCast, n_rank, tmp', n_rank, tmp(1:2,n_rank)
         end do
         do n_rank=1,n_procs
            if( rank < n_rank-1 ) then
               do n_r=nRstart3D,nRstop3D
                  n_theta=1
                  n_z_t  =this%interp_zpb_thw(n_rank,n_theta,n_r)
                  thw_Rloc(n_theta,n_r)=thw_Rloc(n_theta,n_r) +                    &
                  &                     this%interp_wtb_thw(n_rank,n_theta,n_r,1)* &
                  &                     tmp(n_z_t,n_rank)
                  do n_theta=2,n_theta_max/2
                     n_z_t  =this%interp_zpb_thw(n_rank,n_theta,n_r)
                     if ( n_z_t > 1 ) then
                        thw_Rloc(n_theta,n_r)= thw_Rloc(n_theta,n_r) +                   &
                        &                    (this%interp_wtb_thw(n_rank,n_theta,n_r,1)* &
                        & tmp(n_z_t,n_rank) + this%interp_wtb_thw(n_rank,n_theta,n_r,2)* &
                        & tmp(n_z_t-1,n_rank))
                     end if
                  end do
               end do
            end if
         end do
      end if

      !-- Add thermal wind to u_phi
      do n_r=nRstart3D,nRstop3D
         do n_theta=1,n_theta_max/2
            n_t_cth=n_theta_max+1-n_theta
            up_Rloc(:,n_theta,n_r)=up_Rloc(:,n_theta,n_r) + thw_Rloc(n_theta,n_r)
            up_Rloc(:,n_t_cth,n_r)=up_Rloc(:,n_t_cth,n_r) + thw_Rloc(n_theta,n_r)
         end do
      end do

   end subroutine compute_thermal_wind
!--------------------------------------------------------------------------------
   subroutine fill_mat(this)

      class(zfunc_type) :: this

      !-- Local arrays
      real(cp) :: zz(0:(nRstop3D-nRstart3D)+1)

      !-- Local variables
      integer :: n_r, n_t, n_z, n_r_r, n_t_t
      integer :: n_start, n_rank
      real(cp) :: r_r, z_r, c_t, h
      real(cp) :: norm, alpha_r, alpha_t
      real(cp) :: s_r, dz, th

      !-- Get z grid
      !-- for interpolation: n_r_max points on z-axis
      do n_r=1,n_r_max!-1 loop on all s
         h = half*height(n_r)
         !n_r_r = n_r_max_3D-1
         n_r_r = 2
         n_t_t = 2
         do n_z=n_z_max,1,-1
            z_r = h*n_z/n_z_max
            r_r = sqrt(r(n_r)**2 + z_r**2)
            c_t = z_r/r_r
            !n_r_r=minloc(abs(r_3D(:)-r_r),dim=1)
            !if (n_r_r==1) n_r_r=2
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
            !-- Coordinates neighbourg 3: n_r_r-1, n_t_t
            this%interp_zr_mat(4*n_z-1,n_r)=n_r_r-1
            this%interp_zt_mat(4*n_z-1,n_r)=n_t_t
            this%interp_wt_mat(4*n_z-1,n_r)=norm*alpha_r*(one-alpha_t)
            !-- Coordinates neighbourg 4: n_r_r-1, n_t_t-1
            this%interp_zr_mat(4*n_z,n_r)  =n_r_r-1
            this%interp_zt_mat(4*n_z,n_r)  =n_t_t-1
            this%interp_wt_mat(4*n_z,n_r)  =norm*alpha_r*alpha_t
         end do
      end do

      !-- Get theta weights for thermal wind calculation
      zz(:) = 0.0_cp
      !do n_r=n_r_max_3D,2,-1
      do n_r=2,n_r_max_3D !-- Flip the loop since the radii are decreasing
      if( n_r>=nRstart3D .and. n_r<=nRstop3D ) then
         do n_t=1,n_theta_max/2
            s_r = r_3D(n_r)*sint(n_t)
            n_z = 0
            !n_start = min(nRstop,n_r_max_3D-1)
            !do n_r_r=n_start,n_r,-1 !-- Normal loop
            n_start = max(2, nRstart3D)
            do n_r_r=n_start,n_r !-- Flip the loop since our radii are decreasing.
               n_z = n_z+1
               th  = asin(s_r*or1_3D(n_r_r))
               zz(n_z)=sqrt(r_3D(n_r_r)**2-s_r**2)
               if( th < theta(1) ) then
                  n_t_t = 1
                  this%interp_zp_thw(n_z,n_t,n_r,1)=n_r_r
                  this%interp_zp_thw(n_z,n_t,n_r,2)=n_t_t
                  this%interp_wt_thw(n_z,n_t,n_r,1)=one
                  this%interp_wt_thw(n_z,n_t,n_r,2)=0.0_cp
               else
                  n_t_t = 2
                  do while( .not.(th>=theta(n_t_t-1) .and. & 
                  &         th<=theta(n_t_t))  )
                     n_t_t=n_t_t+1
                  end do
                  if ( n_r_r==n_r ) then
                     this%interp_zp_thw(n_z,n_t,n_r,1)=n_r
                     this%interp_zp_thw(n_z,n_t,n_r,2)=n_t
                     this%interp_wt_thw(n_z,n_t,n_r,1)=one
                     this%interp_wt_thw(n_z,n_t,n_r,2)=0.0_cp
                  else
                     this%interp_zp_thw(n_z,n_t,n_r,1)=n_r_r
                     this%interp_zp_thw(n_z,n_t,n_r,2)=n_t_t
                     this%interp_wt_thw(n_z,n_t,n_r,1)=(th-theta(n_t_t-1))/ &
                     &                          (theta(n_t_t)-theta(n_t_t-1))
                     this%interp_wt_thw(n_z,n_t,n_r,2)=(theta(n_t_t)-th)/   &
                     &                          (theta(n_t_t)-theta(n_t_t-1))
                  end if
               end if
            end do
            this%nzp_thw(n_t,n_r)=n_z
            zz(0)=sqrt(r_3D(n_start-1)**2-s_r**2)
            do n_z=1,this%nzp_thw(n_t,n_r)
               dz=zz(n_z-1)-zz(n_z)
               this%interp_wt_thw(n_z,n_t,n_r,1)=dz* &
                    this%interp_wt_thw(n_z,n_t,n_r,1)
               this%interp_wt_thw(n_z,n_t,n_r,2)=dz* &
                    this%interp_wt_thw(n_z,n_t,n_r,2)
            end do
            if( rank /= n_procs-1 ) then
               do n_rank=1,n_procs-1
                  n_r_r=nRstart3D
                  th  = asin(s_r*or1_3D(n_r_r))
                  if( th < theta(1) ) then
                     n_t_t = 1
                     this%interp_zpb_thw(n_rank,n_t,n_r)  =n_t_t
                     this%interp_wtb_thw(n_rank,n_t,n_r,1)=one
                     this%interp_wtb_thw(n_rank,n_t,n_r,2)=0.0_cp
                  else
                     n_t_t = 2
                     do while( .not.(th>=theta(n_t_t-1) .and. &
                     &         th<=theta(n_t_t))  )
                        n_t_t=n_t_t+1
                  end do
                     this%interp_zpb_thw(n_rank,n_t,n_r)  =n_t_t
                     this%interp_wtb_thw(n_rank,n_t,n_r,1)=(th-theta(n_t_t-1))/ &
                     &                          (theta(n_t_t)-theta(n_t_t-1))
                     this%interp_wtb_thw(n_rank,n_t,n_r,2)=(theta(n_t_t)-th)/   &
                     &                          (theta(n_t_t)-theta(n_t_t-1))
                  end if
               end do
            end if
         end do
      end if
      end do

#ifdef DEBUG
      do n_r=nRstart3D,nRstop3D
         do n_t=1,n_theta_max/2
            do n_z=1,n_procs-1
               if( this%interp_zpb_thw(n_z,n_t,n_r)<1 .or. &
               &   this%interp_zpb_thw(n_z,n_t,n_r)>n_theta_max/2 ) then
                  print*, 'WARNING!!:: segFault in fill mat interp_zpb_thw!!!'
               end if
               if( this%interp_zpb_thw(n_z,n_t,n_r) == 1 .and. &
               !&   (this%interp_wtb_thw(n_z,n_t,n_r,1) /= 0.0_cp .or. &
               &   this%interp_wtb_thw(n_z,n_t,n_r,2) /= 0.0_cp) then!)  then
                  print*, 'zpb! n_z_r, n_z_t, wt1, wt2', this%interp_zpb_thw(n_z,n_t,n_r),            &
                  &        this%interp_zpb_thw(n_z,n_t,n_r), this%interp_wtb_thw(n_z,n_t,n_r,1), &
                  &        this%interp_wtb_thw(n_z,n_t,n_r,2)
               end if
            end do
            do n_z=1,this%nzp_thw(n_t,n_r)
               if( this%interp_zp_thw(n_z,n_t,n_r,1)<nRstart3D .or. &
               &   this%interp_zp_thw(n_z,n_t,n_r,1)>nRstop3D .or. &
               &   this%interp_zp_thw(n_z,n_t,n_r,2)<1 .or. &
               &   this%interp_zp_thw(n_z,n_t,n_r,2)>n_theta_max/2 ) then
                  print*, 'WARNING!!:: segFault in fill mat interp_zp_thw!!!'
               end if
               if( this%interp_zp_thw(n_z,n_t,n_r,2) == 1 .and. &
               !&   (this%interp_wt_thw(n_z,n_t,n_r,1) /= 0.0_cp .or. &
               &   this%interp_wt_thw(n_z,n_t,n_r,2) /= 0.0_cp) then!)  then
                  print*, '   n_r_r, n_t_t, n_z, rank', n_r, n_t, n_z, rank
                  print*, 'zp n_z_r, n_z_t,  wt1, wt2', this%interp_zp_thw(n_z,n_t,n_r,1),           &
                  &        this%interp_zp_thw(n_z,n_t,n_r,2), this%interp_wt_thw(n_z,n_t,n_r,1), &
                  &        this%interp_wt_thw(n_z,n_t,n_r,2)
               end if
            end do
         end do
      end do
#endif

   end subroutine fill_mat
!--------------------------------------------------------------------------------
end module z_functions
