module z_functions
   !
   ! This module implements the calculation of the z-averaging
   ! Include reconstruction-interpolation of 3D vel from 2D
   !

   use precision_mod
   use parallel_mod
   use mem_alloc
   use shtns, only: scal_to_spat
   use fourier, only: ifft, fft
   use communications, only: allgather_from_rloc, allgather_from_rloc_3D, &
   &                         transp_lm2r,lm2r_fields, transp_r2lm,r2lm_fields!, transp_r2m3D, r2m3D_fields
   use constants, only: zero, half, one, two, ci, pi
   use blocking, only: nRstart, nRstop, nRstart3D, nRstop3D, lmStart, lmStop!, nphiStart, nphiStop!, radial_balance_3D
   use blocking_lm, only: lm2l, lm2m, lm2lmP
   use truncation, only: n_r_max, n_m_max, minc, idx2m, m2idx
   use truncation_3D, only: n_r_max_3D, n_z_max, n_m_max_3D, n_theta_max, &
       &                    minc_3D, idx2m3D, n_phi_max_3D, lm_max
   use namelists, only: r_icb, r_cmb, l_ek_pump, ktopv, CorFac, ek, ra, &
       &                BuoFac
   use horizontal, only: theta, cost, sint
   use radial_functions, only: r, r_3D, beta, height, oheight, ekpump, or1_3D, &
       &                       rgrav_3D, rscheme_3D
   use radial_der, only: get_dr

   implicit none

   private

   type, public :: zfunc_type

      integer, allocatable :: nzp_thw(:,:)
      integer, allocatable :: interp_zr_mat(:,:)
      integer, allocatable :: interp_zt_mat(:,:)
      integer, allocatable :: interp_zp_thw(:,:,:,:)
      !integer, allocatable :: interp_zpb_thw(:,:,:)
      real(cp), allocatable :: interp_wt_mat(:,:)
      real(cp), allocatable :: interp_wt_thw(:,:,:,:)
      !real(cp), allocatable :: interp_wtb_thw(:,:,:,:)
      real(cp), allocatable :: us_phys_Rloc(:,:)
      real(cp), allocatable :: up_phys_Rloc(:,:)
      real(cp), allocatable :: ek_phys_Rloc(:,:)
      !real(cp), allocatable :: us_phys_Mloc(:,:)
      !real(cp), allocatable :: up_phys_Mloc(:,:)
      !real(cp), allocatable :: ek_phys_Mloc(:,:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: compute_avg
      procedure :: prepare_extension
      procedure :: extrapolate
      procedure :: compute_thermal_wind
      procedure :: compute_lorentz_force
      procedure :: fill_mat
   end type zfunc_type

contains

   subroutine initialize(this)

      class(zfunc_type) :: this

      !-- Local variable
      integer :: n_size

      allocate( this%interp_zr_mat(4*n_z_max,n_r_max) )!nRstart:nRstop) )!
      allocate( this%interp_zt_mat(4*n_z_max,n_r_max) )!nRstart:nRstop) )!
      allocate( this%interp_wt_mat(4*n_z_max,n_r_max) )!nRstart:nRstop) )!

      this%interp_zr_mat(:,:)=1
      this%interp_zt_mat(:,:)=1
      this%interp_wt_mat(:,:)=0.0_cp

      bytes_allocated = bytes_allocated+8*(n_z_max*n_r_max)*SIZEOF_INTEGER
      bytes_allocated = bytes_allocated+4*(n_z_max*n_r_max)*SIZEOF_DEF_REAL

      n_size=n_z_max/2!nRstop3D-nRstart3D+1
      allocate( this%nzp_thw(n_theta_max/2,nRstart3D:nRstop3D) )!n_r_max_3D) )!
      allocate( this%interp_zp_thw(2,n_size,n_theta_max/2,nRstart3D:nRstop3D) )!n_r_max_3D) )!
      !allocate( this%interp_zpb_thw(n_theta_max/2,nRstart3D:nRstop3D,0:n_procs-1) )
      allocate( this%interp_wt_thw(2,n_size,n_theta_max/2,nRstart3D:nRstop3D) )!n_r_max_3D) )!
      !allocate( this%interp_wtb_thw(2,n_theta_max/2,nRstart3D:nRstop3D,0:n_procs-1))

      this%nzp_thw(:,:)=1
      this%interp_zp_thw(:,:,:,:)=1
      !this%interp_zpb_thw(:,:,:)=1
      this%interp_wt_thw(:,:,:,:)=0.0_cp
      !this%interp_wtb_thw(:,:,:,:)=0.0_cp

      bytes_allocated = bytes_allocated+(2*n_size*n_theta_max/2* &
      &                 (nRstop-nRstart+1))*SIZEOF_INTEGER
      !&                 n_r_max_3D)*SIZEOF_INTEGER
      bytes_allocated = bytes_allocated+2*(n_size*n_theta_max/2* &
      &                 (nRstop-nRstart+1))*SIZEOF_DEF_REAL
      !&                 n_r_max_3D)*SIZEOF_DEF_REAL

      allocate( this%us_phys_Rloc(n_phi_max_3D,nRstart:nRstop) )!n_r_max) )!
      allocate( this%up_phys_Rloc(n_phi_max_3D,nRstart:nRstop) )!n_r_max) )!
      allocate( this%ek_phys_Rloc(n_phi_max_3D,nRstart:nRstop) )!n_r_max) )!

      !allocate( this%us_phys_Mloc(nphiStart:nphiStop,n_r_max) )!n_r_max) )!
      !allocate( this%up_phys_Mloc(nphiStart:nphiStop,n_r_max) )!n_r_max) )!
      !allocate( this%ek_phys_Mloc(nphiStart:nphiStop,n_r_max) )!n_r_max) )!

      this%us_phys_Rloc(:,:)=0.0_cp
      this%up_phys_Rloc(:,:)=0.0_cp
      this%ek_phys_Rloc(:,:)=0.0_cp

      !this%us_phys_Mloc(:,:)=0.0_cp
      !this%up_phys_Mloc(:,:)=0.0_cp
      !this%ek_phys_Mloc(:,:)=0.0_cp

      bytes_allocated = bytes_allocated+3*(n_phi_max_3D*      &
      &                 (nRstop-nRstart+1))*SIZEOF_DEF_REAL
      !&                 n_r_max)*SIZEOF_DEF_REAL

   end subroutine initialize
!--------------------------------------------------------------------------------
   subroutine finalize(this)

      class(zfunc_type) :: this

      deallocate( this%interp_zr_mat, this%interp_zt_mat, this%interp_wt_mat )
      deallocate( this%nzp_thw, this%interp_zp_thw, this%interp_wt_thw )
      !deallocate( this%interp_zpb_thw, this%interp_wtb_thw )
      deallocate( this%us_phys_Rloc, this%up_phys_Rloc, this%ek_phys_Rloc )
      !deallocate( this%us_phys_Mloc, this%up_phys_Mloc, this%ek_phys_Mloc )

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

      !real(cp) :: usm3D_Mloc(n_phi_max_3D,nRstart:nRstop)
      !real(cp) :: upm3D_Mloc(n_phi_max_3D,nRstart:nRstop)
      !real(cp) :: ekpump_p3D(n_phi_max_3D,nRstart:nRstop)
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
         !call ifft(usm3D_Rloc(:,n_r), usm3D_Mloc(:,n_r), l_3D=.true.)
         !call ifft(upm3D_Rloc(:,n_r), upm3D_Mloc(:,n_r), l_3D=.true.)
         !if( l_ek_pump ) &
         !&   call ifft(ekpump_m3D(:,n_r), ekpump_p3D(:,n_r), l_3D=.true.)
      end do

      !call transp_r2m3D(r2m3D_fields, usm3D_Mloc, this%us_phys_Mloc)
      !call transp_r2m3D(r2m3D_fields, upm3D_Mloc, this%up_phys_Mloc)
      !if( l_ek_pump ) call transp_r2m3D(r2m3D_fields, ekpump_p3D, this%ek_phys_Mloc)

      !-- Boundary point: fix Ek-pumping to zero
      if ( rank == 0 ) this%ek_phys_Rloc(:,1)=0.0_cp
      !if ( rank == 0 ) this%ek_phys_Mloc(:,1)=0.0_cp

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
      integer :: n_r, n_r_r, n_phi, n_th_NHS, n_th_SHS
      real(cp) :: s_r, z_r, z_eta
      real(cp) :: vs, vz, vrr, vth, vph
      real(cp) :: alpha_r1, alpha_r2

      call allgather_from_rloc(this%us_phys_Rloc,usr,n_phi_max_3D)
      call allgather_from_rloc(this%up_phys_Rloc,upp,n_phi_max_3D)
      if( l_ek_pump ) call allgather_from_rloc(this%ek_phys_Rloc,ekp,n_phi_max_3D)

      !-- Compute 3D velocity fields by a linear interpolation
      do n_r_r=nRstart3D,nRstop3D
         do n_th_NHS=1,n_theta_max/2
            n_th_SHS=n_theta_max+1-n_th_NHS
            s_r = r_3D(n_r_r)*sint(n_th_NHS)
            z_r = r_3D(n_r_r)*cost(n_th_NHS)
            if ( s_r >= r_icb ) then !-- Outside TC
               n_r = 1
               do while ( r(n_r) > s_r .and. n_r < n_r_max )
                  n_r = n_r+1
               end do
               alpha_r2 = (s_r-r(n_r))/(r(n_r-1)-r(n_r))
               alpha_r1 = one - alpha_r2
               z_eta = -s_r/(r_cmb*r_cmb-s_r*s_r)*z_r ! \beta * z
               do n_phi=1,n_phi_max_3D
               !do n_phi=nphiStart,nphiStop
                  vs = alpha_r1*usr(n_phi,n_r) + alpha_r2*usr(n_phi,n_r-1)
                  !vs = alpha_r1*this%us_phys_Mloc(n_phi,n_r) + alpha_r2*this%us_phys_Mloc(n_phi,n_r-1)
                  !-- vz = beta*vs
                  vz = z_eta*vs
                  if ( l_ek_pump ) then
                     !-- vz = beta*vz+ekpump
                     vz = vz + z_r*(alpha_r1*ekp(n_phi,n_r) +  & 
                     &              alpha_r2*ekp(n_phi,n_r-1))
                     !vz = vz + z_r*(alpha_r1*this%ek_phys_Mloc(n_phi,n_r) +  &
                     !&              alpha_r2*this%ek_phys_Mloc(n_phi,n_r-1))
                  end if
                  vrr= vz*cost(n_th_NHS) + vs*sint(n_th_NHS)
                  ur_Rloc(n_phi,n_th_NHS,n_r_r)= vrr
                  ur_Rloc(n_phi,n_th_SHS,n_r_r)= vrr
                  vth= vs*cost(n_th_NHS) - vz*sint(n_th_NHS)
                  vph= alpha_r1*upp(n_phi,n_r) + alpha_r2*upp(n_phi,n_r-1)
                  !vph= alpha_r1*this%up_phys_Mloc(n_phi,n_r) + alpha_r2*this%up_phys_Mloc(n_phi,n_r-1)
                  ut_Rloc(n_phi,n_th_NHS,n_r_r)= vth
                  ut_Rloc(n_phi,n_th_SHS,n_r_r)=-vth
                  up_Rloc(n_phi,n_th_NHS,n_r_r)= vph
                  up_Rloc(n_phi,n_th_SHS,n_r_r)= vph
               end do

            else !-- Inside the tangent cylinder

               do n_phi=1,n_phi_max_3D
               !do n_phi=nphiStart,nphiStop
                  ur_Rloc(n_phi,n_th_NHS,n_r_r)=0.0_cp
                  ur_Rloc(n_phi,n_th_SHS,n_r_r)=0.0_cp
                  ut_Rloc(n_phi,n_th_NHS,n_r_r)=0.0_cp
                  ut_Rloc(n_phi,n_th_SHS,n_r_r)=0.0_cp
                  up_Rloc(n_phi,n_th_NHS,n_r_r)=0.0_cp
                  up_Rloc(n_phi,n_th_SHS,n_r_r)=0.0_cp
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
      integer :: n_r_r, n_th_NHS, n_th_SHS
      complex(cp) :: tmp_hat(n_m_max_3D)
      real(cp) :: czavg
      real(cp) :: tmp(n_phi_max_3D,n_r_max)!nRstart:nRstop)!

      tmp(:,:)      =0.0_cp
      zavg_Rloc(:,:)=zero

      !-- z-average in spatial space
      do n_r=1,n_r_max
      !do n_r=nRstart,nRstop
         do n_z=1,4*n_z_max
            n_th_NHS= this%interp_zt_mat(n_z,n_r)
            n_th_SHS= n_theta_max+1-n_th_NHS
            n_r_r = this%interp_zr_mat(n_z,n_r)
            czavg = this%interp_wt_mat(n_z,n_r)
            if ( n_r_r >= nRstart3D .and. n_r_r <= nRstop3D ) then
               do n_phi=1,n_phi_max_3D
                  tmp(n_phi,n_r) = tmp(n_phi,n_r) + czavg*              &
                  &                   (work_Rloc(n_phi,n_th_NHS,n_r_r)  &
                  &                  + work_Rloc(n_phi,n_th_SHS,n_r_r))
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
   subroutine compute_thermal_wind(this, dTdth_Rloc, up_Rloc)

      !-- Input variables
      class(zfunc_type) :: this
      real(cp), intent(in) :: dTdth_Rloc(n_theta_max,nRstart3D:nRstop3D)

      !-- Output variables - modified (inout)
      real(cp), intent(inout) :: up_Rloc(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)

      !-- Local arrays
      real(cp) :: dTdth(n_theta_max,n_r_max_3D)!nRstart3D:nRstop3D)!
      real(cp) :: thw_Rloc(n_theta_max/2,n_r_max_3D)!nRstart3D:nRstop3D)!
      real(cp) :: dTzdt(n_theta_max/2,n_r_max_3D)!nRstart3D:nRstop3D)!
      !real(cp) :: tmp(n_theta_max/2)
      !real(cp) :: Zwb(n_theta_max/2,0:n_procs-1)
      !-- Local variables
      integer :: n_th_NHS, n_z, n_z_r, n_z_t, n_th_SHS
      integer :: n_r, n_p
      real(cp) :: thwFac

      thwFac=BuoFac/CorFac

      call allgather_from_rloc_3D(dTdth_Rloc,dTdth,n_theta_max)
      !dTdth(:,:)=dTdth_Rloc(:,:)

      !-- Remaining term for the temperature gradient
      dTzdt(:,:)=0.0_cp
      do n_r=2,n_r_max_3D!max(2,nRstart3D),nRstop3D
      !do n_r=max(2,nRstart3D),nRstop3D
         do n_th_NHS=1,n_theta_max/2
            !dTzdt(n_th_NHS,n_r)=thwFac*r_3D(n_r)* dTdth(n_th_NHS,n_r)
            !-- TG I don't understand the r factor in the above equation
            !-- Th wind should be
            !-- duphi/dz = Ra/Pr  * g / r * dT/dtheta
            dTzdt(n_th_NHS,n_r)=thwFac*rgrav_3D(n_r)*or1_3D(n_r)* &
            !&                   cos(theta(n_th_NHS))
            &                   dTdth(n_th_NHS,n_r)
            !&                   dTdth_Rloc(n_th_NHS,n_r)
         end do
      end do

      !-- Compute thermal wind
      thw_Rloc(:,:)=0.0_cp
      do n_r=nRstart3D,nRstop3D
         n_th_NHS=1
         do n_z=1,this%nzp_thw(n_th_NHS,n_r)
            n_z_r = this%interp_zp_thw(1,n_z,n_th_NHS,n_r)
            n_z_t = this%interp_zp_thw(2,n_z,n_th_NHS,n_r)
            thw_Rloc(n_th_NHS,n_r)=thw_Rloc(n_th_NHS,n_r) -               &
            &                     this%interp_wt_thw(1,n_z,n_th_NHS,n_r)* &
            &                     dTzdt(n_z_t,n_z_r)
         end do
         do n_th_NHS=2,n_theta_max/2
            do n_z=1,this%nzp_thw(n_th_NHS,n_r)
               n_z_r = this%interp_zp_thw(1,n_z,n_th_NHS,n_r)
               n_z_t = this%interp_zp_thw(2,n_z,n_th_NHS,n_r)
               if ( n_z_t > 1 ) then
                  thw_Rloc(n_th_NHS,n_r)= thw_Rloc(n_th_NHS,n_r) -               &
                  &                     (this%interp_wt_thw(1,n_z,n_th_NHS,n_r)* &
                  & dTzdt(n_z_t,n_z_r) + this%interp_wt_thw(2,n_z,n_th_NHS,n_r)* &
                  & dTzdt(n_z_t-1,n_z_r))
               end if
            end do
         end do
      end do
      !if( n_procs>1 ) then
      !   Zwb(:,:)=0.0_cp
      !   do n_p=0,n_procs-1
      !      if( n_p == rank ) tmp(:) = thw_Rloc(:,nRstop3D)
      !      call MPI_Bcast(tmp, n_theta_max/2, MPI_DEF_REAL, &
      !      &              n_p, MPI_COMM_WORLD, ierr)
      !      ZWb(:,n_p)=tmp(:)
      !      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      !   end do
      !   do n_p=0,n_procs-1
      !      if( rank > n_p ) then
      !         do n_r=nRstart3D,nRstop3D
      !            n_th_NHS=1
      !            n_z_t  =this%interp_zpb_thw(n_th_NHS,n_r,n_p)
      !            thw_Rloc(n_th_NHS,n_r)=thw_Rloc(n_th_NHS,n_r) +               &
      !            &                     this%interp_wtb_thw(1,n_th_NHS,n_r,n_p)*&
      !            &                     Zwb(n_z_t,n_p)
      !            do n_th_NHS=2,n_theta_max/2
      !               n_z_t  =this%interp_zpb_thw(n_th_NHS,n_r,n_p)
      !               if ( n_z_t > 1 ) then
      !                  thw_Rloc(n_th_NHS,n_r)= thw_Rloc(n_th_NHS,n_r) +        &
      !                  &             (this%interp_wtb_thw(1,n_th_NHS,n_r,n_p)* &
      !                  &             Zwb(n_z_t,n_p) +                          &
      !                  &             this%interp_wtb_thw(2,n_th_NHS,n_r,n_p)*  &
      !                  &             Zwb(n_z_t-1,n_p))
      !               end if
      !            end do
      !         end do
      !      end if
      !   end do
      !end if

      !-- Add thermal wind to u_phi
      do n_r=nRstart3D,nRstop3D
         do n_th_NHS=1,n_theta_max/2
            n_th_SHS=n_theta_max+1-n_th_NHS
            up_Rloc(:,n_th_NHS,n_r)=up_Rloc(:,n_th_NHS,n_r) + thw_Rloc(n_th_NHS,n_r)
            up_Rloc(:,n_th_SHS,n_r)=up_Rloc(:,n_th_SHS,n_r) + thw_Rloc(n_th_NHS,n_r)
         end do
      end do

   end subroutine compute_thermal_wind
!--------------------------------------------------------------------------------
   subroutine compute_lorentz_force(this, lrf0, lrfLM, lrf_Rloc, lrf0_Rloc)

      !-- Input variables
      class(zfunc_type) :: this
      real(cp), intent(in) :: lrf0(n_theta_max,nRstart3D:nRstop3D)
      complex(cp), intent(inout) :: lrfLM(lm_max,nRstart3D:nRstop3D,2)

      !-- Output variables - modified (inout)
      real(cp), intent(inout) :: lrf_Rloc(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      complex(cp), intent(inout) :: lrf0_Rloc(nRstart:nRstop)

      !-- Local arrays
      real(cp) :: tmp_lrf0(2,n_r_max)
      complex(cp) :: lrf_r(lmStart:lmStop,n_r_max_3D)
      complex(cp) :: tmp_lrf(lm_max,nRstart3D:nRstop3D)
      !-- Local variables
      integer :: n_z, n_m_QG, n_m, m3D
      integer :: n_r_r, n_th_NHS, n_th_SHS
      integer :: n_r, lm, l, m, lmP
      real(cp) :: czavg


      !call allgather_from_rloc_3D(lrfLM(:,:,1),lrf_r,lm_max)
      !call get_dr( lrf_r, lrf_r, n_r_max_3D, rscheme_3D, nocopy=.true. )
   if( rank==0 ) print*, 'INSIDE LORENTZ_force - before get dr jxB_t',  lrf_Rloc(1,1,1)
      call transp_r2lm(r2lm_fields, lrfLM(:,:,1),lrf_r)
      call get_dr( lrf_r, lrf_r, lmStart, lmStop, &
           &       n_r_max_3D, rscheme_3D, nocopy=.true. )
      call transp_lm2r(lm2r_fields, lrf_r, lrfLM(:,:,1))

   if( rank==0 ) print*, 'INSIDE LORENTZ_force - before computing VxjxB_z',  lrf_Rloc(1,1,1)
      !-- Remaining term for the lorentz force
      !!$OMP PARALLEL DO default(shared) &
      !!$OMP& private(n_theta, n_phi)
      !-- Assemble curl jxB . e_z on the spherical grid
      !-- = 1/r * (dr(r jxB_t) - dth(jxB_r))
      do n_r=nRstart3D,nRstop3D
         do lm=2,lm_max
            l   =lm2l(lm)
            m   =lm2m(lm)
            lmP =lm2lmP(lm)
            tmp_lrf(lmP,n_r)=or1_3D(n_r)*(    &
            &                lrfLM(lmP,n_r,1)-&
            &                lrfLM(lmP,n_r,2))
         end do
      end do
      !!$OMP END PARALLEL DO

   if( rank==0 ) print*, 'INSIDE LORENTZ_force - before back to grid',  lrf_Rloc(1,1,1)
      !-- Back to the grid before the z-averaging
      do n_r=nRstart3D,nRstop3D
         call scal_to_spat(tmp_lrf(:,n_r), lrf_Rloc(:,:,n_r))
      end do

   if( rank==0 ) print*, 'INSIDE LORENTZ_force - before computing jxB_p z-avg',  lrf_Rloc(1,1,1)
      !-- z-average of the non-axisymmetric part (m==0)
      !-- = <jxB_p>
      do n_r=1,n_r_max
         do n_z=1,4*n_z_max
            n_th_NHS= this%interp_zt_mat(n_z,n_r)
            n_th_SHS= n_theta_max+1-n_th_NHS
            n_r_r = this%interp_zr_mat(n_z,n_r)
            czavg = this%interp_wt_mat(n_z,n_r)
            if ( n_r_r >= nRstart3D .and. n_r_r <= nRstop3D ) then
               tmp_lrf0(1,n_r) = tmp_lrf0(1,n_r) + czavg* &
               &                    (lrf0(n_th_NHS,n_r_r) &
               &                   + lrf0(n_th_SHS,n_r_r))
            end if
         end do
      end do

   if( rank==0 ) print*, 'INSIDE LORENTZ_force - before AllReduce',  lrf_Rloc(1,1,1)
      do n_r=1,n_r_max
         call MPI_Allreduce(MPI_IN_PLACE, tmp_lrf0(:,n_r), 2, MPI_DEF_REAL, &
              &             MPI_SUM, MPI_COMM_WORLD, ierr)
      end do

   if( rank==0 ) print*, 'INSIDE LORENTZ_force - before back to grid m=0',  lrf_Rloc(1,1,1)
      !-- Transforms back to spectral space
      do n_r=nRstart,nRstop
         !call fft(tmp_lrf0(n_r), lrf0_Rloc(n_r), l_3D=.true.)
         lrf0_Rloc(n_r)=cmplx(tmp_lrf0(1,n_r),0.0,kind=cp)
      end do

   end subroutine compute_lorentz_force
!--------------------------------------------------------------------------------
   subroutine fill_mat(this)

      class(zfunc_type) :: this

      !-- Local arrays
      real(cp) :: zz(0:n_r_max_3D)!(nRstop3D-nRstart3D)+1)

      !-- Local variables
      integer :: n_r, n_t, n_z, n_r_r, n_t_t
      integer :: n_start, n_p
      real(cp) :: r_r, z_r, c_t, h
      real(cp) :: norm, alpha_r, alpha_t
      real(cp) :: s_r, dz, th

      !-- Get z grid
      !-- for interpolation: n_r_max points on z-axis
      do n_r=1,n_r_max!-1 loop on all s
      !do n_r=nRstart,nRstop
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
            if ( n_z /= 0 ) then
               !-- Multiply by 2 because of the symmetries
               norm = two/(two*n_z_max+1)
            else
               norm = one/(two*n_z_max+1)
            endif
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
      !do n_r=2,n_r_max_3D
      !do n_r=max(1,nRstart3D-1),min(nRstop3D+1,n_r_max_3D)
      do n_r=max(2,nRstart3D),nRstop3D
         do n_t=1,n_theta_max/2
            s_r = r_3D(n_r)*sint(n_t)
            n_z = 0
            n_start = 2
            do n_r_r=n_start,n_r!-1 !-- Flip the loop since our radii are decreasing.
               if( n_r_r > n_r_max_3D) print*, 'w!! segFault in n_r_r loop!; n_r_r', n_r_r
               n_z = n_z+1
               th  = asin(s_r*or1_3D(n_r_r))
               if( th < 0.0_cp ) print*, 'w!! th<0; n_z, n_r, n_t =', n_z, n_r, n_t
               zz(n_z)=sqrt(r_3D(n_r_r)**2-s_r**2)
               if( th /= th ) print*, 'w!! nan in th; n_z, n_r, n_t =', th, n_z, n_r, n_t
               if( zz(n_z) /= zz(n_z) ) print*, 'w!! nan in zz; n_z, n_r, n_t =', zz(n_z), n_z, n_r, n_t
               if( th < theta(1) ) then
                  n_t_t = 1
                  this%interp_zp_thw(1,n_z,n_t,n_r)=n_r_r
                  this%interp_zp_thw(2,n_z,n_t,n_r)=n_t_t
                  this%interp_wt_thw(1,n_z,n_t,n_r)=one
                  this%interp_wt_thw(2,n_z,n_t,n_r)=0.0_cp
               else
                  n_t_t = 2
                  do while( .not.(th>=theta(n_t_t-1) .and. & 
                  &         th<=theta(n_t_t)) .and. n_t_t < n_theta_max/2 ) 
                     n_t_t=n_t_t+1
                  end do
                  if( n_t_t > n_theta_max/2 ) print*, 'w!! segFault in n_t_t loop!; n_t_t', n_t_t
                  if ( n_r_r==n_r ) then
                     this%interp_zp_thw(1,n_z,n_t,n_r)=n_r
                     this%interp_zp_thw(2,n_z,n_t,n_r)=n_t!heta_max/2
                     this%interp_wt_thw(1,n_z,n_t,n_r)=one
                     this%interp_wt_thw(2,n_z,n_t,n_r)=0.0_cp
                  else
                     this%interp_zp_thw(1,n_z,n_t,n_r)=n_r_r
                     this%interp_zp_thw(2,n_z,n_t,n_r)=n_t_t
                     this%interp_wt_thw(1,n_z,n_t,n_r)=(th-theta(n_t_t-1))/ &
                     &                          (theta(n_t_t)-theta(n_t_t-1))
                     this%interp_wt_thw(2,n_z,n_t,n_r)=(theta(n_t_t)-th)/   &
                     &                          (theta(n_t_t)-theta(n_t_t-1))
                  end if
               end if
            end do
            this%nzp_thw(n_t,n_r)=n_z
            if( n_z > n_r_max_3D) print*, 'w!! segFault in n_z loop!; n_z', n_z
            zz(0)=sqrt(r_3D(n_start-1)**2-s_r**2)
            do n_z=1,this%nzp_thw(n_t,n_r)
               dz=zz(n_z-1)-zz(n_z)
               if( dz /= dz ) print*, 'w!! nan in dz; n_z, n_r, n_t =', dz, n_z, n_r, n_t
               if( dz < 0.0_cp ) print*, 'w!! dz<0; n_z, n_r, n_t =', n_z, n_r, n_t
               this%interp_wt_thw(1,n_z,n_t,n_r)=dz* &
               &              this%interp_wt_thw(1,n_z,n_t,n_r)
               this%interp_wt_thw(2,n_z,n_t,n_r)=dz* &
               &              this%interp_wt_thw(2,n_z,n_t,n_r)
            end do
            !do n_p=0,n_procs-1
            !   if ( rank > n_p ) then
            !      n_r_r=radial_balance_3D(n_p)%nStop
            !      th  = asin(s_r*or1_3D(n_r_r))
            !      if( th < theta(1) ) then
            !         n_t_t = 1
            !         this%interp_zpb_thw(n_t,n_r,n_p)  =n_t_t
            !         this%interp_wtb_thw(1,n_t,n_r,n_p)=one
            !         this%interp_wtb_thw(2,n_t,n_r,n_p)=0.0_cp
            !      else
            !         n_t_t = 2
            !         do while (.not.(th>=theta(n_t_t-1) .and. th<=theta(n_t_t)))
            !            n_t_t=n_t_t+1
            !         end do
            !         this%interp_zpb_thw(n_t,n_r,n_p)  =n_t_t
            !         this%interp_wtb_thw(1,n_t,n_r,n_p)=(th-theta(n_t_t-1))/ &
            !         &                          (theta(n_t_t)-theta(n_t_t-1))
            !         this%interp_wtb_thw(2,n_t,n_r,n_p)=(theta(n_t_t)-th)/   &
            !         &                          (theta(n_t_t)-theta(n_t_t-1))
            !      end if
            !   end if
            !end do
         end do
      end do

   end subroutine fill_mat
!--------------------------------------------------------------------------------
end module z_functions
