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
   use communications, only: allgather_from_rloc, allgather_from_rloc_3D,      &
   &                         transp_lm2r,lm2r_fields, transp_r2lm,r2lm_fields, &
   &                         exchange_Nbound_from_Rloc_3D
   use constants, only: zero, half, third, one, two, four, ci, pi
   use blocking, only: nRstart, nRstop, nRstart3D, nRstop3D, lmStart, lmStop
   use blocking_lm, only: lm2l, lm2m, lm2lmP
   use truncation, only: n_r_max, n_z_max, n_m_max, minc, idx2m, m2idx
   use truncation_3D, only: n_r_max_3D, n_m_max_3D, n_theta_max,   &
       &                    minc_3D, idx2m3D, n_phi_max_3D, lm_max
   use namelists, only: r_icb, r_cmb, l_ek_pump, ktopv, CorFac, ek, ra, &
       &                BuoFac, l_heat_3D, l_thw_3D, l_cyl, l_mag_pump, &
       &                mag_pump_fac, l_QG_basis
   use horizontal, only: theta, cost, sint
   use radial_functions, only: r, r_3D, beta, oheight, ekpump, or1_3D, &
       &                       rgrav_3D, rscheme_3D

   implicit none

   private

   type, public :: zfunc_type

      !-- Bilinear zavg integration matrices
      integer, allocatable :: nzp_zavg(:)
      integer, allocatable :: interp_zr_mat(:,:)
      integer, allocatable :: interp_zt_mat(:,:)
      real(cp), allocatable :: interp_wt_mat(:,:)

      !-- 4th order zavg integration matrices
      integer, allocatable :: nzp_zcyl(:)
      integer, allocatable :: interp_zr_cyl(:,:)
      integer, allocatable :: interp_zt_cyl(:,:)
      real(cp), allocatable :: interp_wr_cyl(:,:,:)
      real(cp), allocatable :: interp_wt_cyl(:,:,:)

      !-- 2-neighbourgs zthw integration matrices
      integer, allocatable :: nzp_thw(:,:)
      integer, allocatable :: interp_zp_thw(:,:,:,:)
      !integer, allocatable :: interp_zpb_thw(:,:,:)
      real(cp), allocatable :: interp_wt_thw(:,:,:,:)
      !real(cp), allocatable :: interp_wtb_thw(:,:,:,:)

      !-- 2D->3D velocity field extrapolation matrices
      real(cp), allocatable :: us_phys_Rloc(:,:)
      real(cp), allocatable :: up_phys_Rloc(:,:)
      real(cp), allocatable :: ek_phys_Rloc(:,:)

      !-- Magnetic pumping matrices
      real(cp), allocatable :: up_magpump_Rloc(:,:)
      real(cp), allocatable :: uz_magpump_Rloc(:,:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: compute_zavg
      procedure :: cyl_avg
      procedure :: compute_zder
      procedure :: prep_extension_QGto3D
      procedure :: ext_QGto3D_vel
      procedure :: compute_thermal_wind
      procedure :: fill_zinterp_grid
   end type zfunc_type

contains

   subroutine initialize(this)

      class(zfunc_type) :: this

      !-- Local variable
      integer :: n_size

      !-- Bilinear zavg integration
      n_size = 4*(n_z_max+1)
      allocate( this%nzp_zavg(n_r_max) )!nRstart:nRstop) )!
      allocate( this%interp_zr_mat(n_size,n_r_max) )!nRstart:nRstop) )!
      allocate( this%interp_zt_mat(n_size,n_r_max) )!nRstart:nRstop) )!
      allocate( this%interp_wt_mat(n_size,n_r_max) )!nRstart:nRstop) )!

      this%nzp_zavg(:)=2 ! Minimum of 2 pts for Bilinear integration
      this%interp_zr_mat(:,:)=1
      this%interp_zt_mat(:,:)=1
      this%interp_wt_mat(:,:)=0.0_cp

      bytes_allocated = bytes_allocated+(2*n_size+1)*n_r_max*SIZEOF_INTEGER
      bytes_allocated = bytes_allocated+(n_size*n_r_max)*SIZEOF_DEF_REAL

      !-- 4th order zavg integration matrices
      if ( l_cyl ) then
         n_size = n_z_max+1!-- = 2*n_s_max+1
         allocate( this%nzp_zcyl(n_r_max) )!nRstart:nRstop) )!
         allocate( this%interp_zr_cyl(-n_r_max:n_r_max,n_r_max) )!nRstart:nRstop) )!
         allocate( this%interp_zt_cyl(-n_r_max:n_r_max,n_r_max) )!nRstart:nRstop) )!
         allocate( this%interp_wr_cyl(0:9,-n_r_max:n_r_max,n_r_max) )!nRstart:nRstop) )!
         allocate( this%interp_wt_cyl(0:9,-n_r_max:n_r_max,n_r_max) )!nRstart:nRstop) )!

         this%nzp_zcyl(:)=4 ! Minimum of 4 pts for Simpson integration
         this%interp_zr_cyl(:,:)=1
         this%interp_zt_cyl(:,:)=1
         this%interp_wr_cyl(:,:,:)=0.0_cp
         this%interp_wt_cyl(:,:,:)=0.0_cp

         bytes_allocated = bytes_allocated+(2*n_size+1)*n_r_max*SIZEOF_INTEGER
         bytes_allocated = bytes_allocated+2*(10*n_size*n_r_max)*SIZEOF_DEF_REAL
      else
         allocate( this%nzp_zcyl(1) )
         allocate( this%interp_zr_cyl(1,1) )
         allocate( this%interp_zt_cyl(1,1) )
         allocate( this%interp_wr_cyl(1,1,1) )
         allocate( this%interp_wt_cyl(1,1,1) )
      end if

      !-- 2-neighbourgs zthw integration
      if ( l_heat_3D .and. l_thw_3D ) then
         n_size=n_r_max_3D!!n_z_max/2!nRstop3D-nRstart3D+1
         allocate( this%nzp_thw(n_theta_max/2+1,nRstart3D:nRstop3D) )!n_r_max_3D) )!
         allocate( this%interp_zp_thw(2,n_size,n_theta_max/2+1,nRstart3D:nRstop3D) )!n_r_max_3D) )!
         !allocate( this%interp_zpb_thw(n_theta_max/2,nRstart3D:nRstop3D,0:n_procs-1) )
         allocate( this%interp_wt_thw(2,n_size,n_theta_max/2+1,nRstart3D:nRstop3D) )!n_r_max_3D) )!
         !allocate( this%interp_wtb_thw(2,n_theta_max/2,nRstart3D:nRstop3D,0:n_procs-1))

         this%nzp_thw(:,:)=1
         this%interp_zp_thw(:,:,:,:)=1
         !this%interp_zpb_thw(:,:,:)=1
         this%interp_wt_thw(:,:,:,:)=0.0_cp
         !this%interp_wtb_thw(:,:,:,:)=0.0_cp

         bytes_allocated = bytes_allocated+(2*n_size+1)*n_theta_max/2* &
         &                 (nRstop-nRstart+1)*SIZEOF_INTEGER
         !&                 n_r_max_3D)*SIZEOF_INTEGER
         bytes_allocated = bytes_allocated+2*(n_size*n_theta_max/2* &
         &                 (nRstop-nRstart+1))*SIZEOF_DEF_REAL
         !&                 n_r_max_3D)*SIZEOF_DEF_REAL
      else
         allocate( this%nzp_thw(1,1) )
         allocate( this%interp_zp_thw(1,1,1,1) )
         allocate( this%interp_wt_thw(1,1,1,1) )
      end if

      !-- 2D->3D velocity field extrapolation + Magnetic pumping
      allocate( this%us_phys_Rloc(n_phi_max_3D,nRstart:nRstop) )!n_r_max) )!
      allocate( this%up_phys_Rloc(n_phi_max_3D,nRstart:nRstop) )!n_r_max) )!
      allocate( this%ek_phys_Rloc(n_phi_max_3D,nRstart:nRstop) )!n_r_max) )!

      this%us_phys_Rloc(:,:)=0.0_cp
      this%up_phys_Rloc(:,:)=0.0_cp
      this%ek_phys_Rloc(:,:)=0.0_cp

      bytes_allocated = bytes_allocated+3*(n_phi_max_3D*      &
      &                 (nRstop-nRstart+1))*SIZEOF_DEF_REAL
      !&                 n_r_max)*SIZEOF_DEF_REAL

      if ( l_mag_pump ) then
         allocate( this%up_magpump_Rloc(n_phi_max_3D,nRstart:nRstop) )!n_r_max) )!
         allocate( this%uz_magpump_Rloc(n_phi_max_3D,nRstart:nRstop) )!n_r_max) )!

         this%up_magpump_Rloc(:,:)=0.0_cp
         this%uz_magpump_Rloc(:,:)=0.0_cp

         bytes_allocated = bytes_allocated+2*(n_phi_max_3D*      &
         &                 (nRstop-nRstart+1))*SIZEOF_DEF_REAL
         !&                 n_r_max)*SIZEOF_DEF_REAL
      else
         allocate( this%up_magpump_Rloc(1,1) )
         allocate( this%uz_magpump_Rloc(1,1) )
      end if

   end subroutine initialize
!--------------------------------------------------------------------------------
   subroutine finalize(this)

      class(zfunc_type) :: this

      deallocate( this%nzp_zavg, this%nzp_zcyl )
      deallocate( this%interp_zr_mat, this%interp_zt_mat, this%interp_wt_mat )
      deallocate( this%interp_zr_cyl, this%interp_zt_cyl, this%interp_wr_cyl, this%interp_wt_cyl )
      deallocate( this%nzp_thw, this%interp_zp_thw, this%interp_wt_thw )
      !deallocate( this%interp_zpb_thw, this%interp_wtb_thw )
      deallocate( this%us_phys_Rloc, this%up_phys_Rloc, this%ek_phys_Rloc )
      deallocate( this%up_magpump_Rloc, this%uz_magpump_Rloc )

   end subroutine finalize
!--------------------------------------------------------------------------------
   subroutine prep_extension_QGto3D(this, us_Rloc, up_Rloc, om_Rloc)

      class(zfunc_type) :: this

      !-- Input variables
      complex(cp), intent(in) :: us_Rloc(n_m_max,nRstart:nRstop)
      complex(cp), intent(in) :: up_Rloc(n_m_max,nRstart:nRstop)
      complex(cp), intent(in) :: om_Rloc(n_m_max,nRstart:nRstop)

      !-- Local variables
      complex(cp) :: usm3D_Rloc(n_m_max_3D,nRstart:nRstop)
      complex(cp) :: upm3D_Rloc(n_m_max_3D,nRstart:nRstop)
      complex(cp) :: ekpump_m3D(n_m_max_3D,nRstart:nRstop)

      complex(cp) :: mag_pump_upm3D(n_m_max_3D,nRstart:nRstop)
      complex(cp) :: mag_pump_uzm3D(n_m_max_3D,nRstart:nRstop)

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
                     &    5.0_cp*r_cmb*oheight(n_r) )*    us_Rloc(n_m,n_r))
                     if ( l_QG_basis ) then
                     !-- Ekman Pumping should not be modified (only dependant on BCs)
                        ekpump_m3D(n_m_3D,n_r)=ekpump_m3D(n_m_3D,n_r) + &
                        &  ekpump(n_r)*beta(n_r)*third*ci*real(m3D,cp)* &
                        &                               us_Rloc(n_m,n_r)
                     end if
                  else
                     ekpump_m3D(n_m_3D,n_r)=zero
                     !ekpump_m3D(n_m_3D,n_r)=-CorFac*                &
                     !&                                 ekpump(n_r)* &
                     !&                                 up_Rloc(n_m,n_r)
                  end if
               end if

               if ( l_mag_pump ) then
                  !-- Mag_pump u_phi porportional to r^2 u_s
                  mag_pump_upm3D(n_m_3D,n_r) = mag_pump_fac*us_Rloc(n_m,n_r)*r(n_r)**2.
                  !-- Mag_pump u_z porportional to omega_z
                  mag_pump_uzm3D(n_m_3D,n_r) = mag_pump_fac*om_Rloc(n_m,n_r)
               end if
            else
               usm3D_Rloc(n_m_3D,n_r) = zero
               upm3D_Rloc(n_m_3D,n_r) = zero
               ekpump_m3D(n_m_3D,n_r) = zero

               mag_pump_upm3D(n_m_3D,n_r) = zero
               mag_pump_uzm3D(n_m_3D,n_r) = zero
            end if
         end do
      end do

         !:: if needed:::
         !   do n_r=2,n_r_max
         !      if ( l_non_rot ) then
         !         h2 = one
         !      else
         !         h2 = r_cmb*r_cmb-r(n_r)*r(n_r)
         !      end if
         !      do n_m=nMstart,nMstop
         !         m = idx2m(n_m)
         !         if ( m > 0 ) then
         !            psi_Mloc(n_m, n_r) = -ci*r(n_r)/real(m,cp)/h2 * us_Mloc(n_m, n_r)
         !         else
         !            psi_Mloc(n_m, n_r) = 0.0_cp
         !         end if
         !      end do

      do n_r=nRstart,nRstop
         call ifft(usm3D_Rloc(:,n_r), this%us_phys_Rloc(:,n_r), l_3D=.true.)
         call ifft(upm3D_Rloc(:,n_r), this%up_phys_Rloc(:,n_r), l_3D=.true.)

         if( l_ek_pump ) &
         &   call ifft(ekpump_m3D(:,n_r), this%ek_phys_Rloc(:,n_r), l_3D=.true.)
         if( l_mag_pump ) then
            call ifft(mag_pump_upm3D(:,n_r), this%up_magpump_Rloc(:,n_r), l_3D=.true.)
            call ifft(mag_pump_uzm3D(:,n_r), this%uz_magpump_Rloc(:,n_r), l_3D=.true.)
         end if
      end do

      !-- Boundary point: fix Ek-pumping to zero
      if ( rank == 0 ) this%ek_phys_Rloc(:,1)=0.0_cp

   end subroutine prep_extension_QGto3D
!--------------------------------------------------------------------------------
   subroutine ext_QGto3D_vel(this, ur_Rloc, ut_Rloc, up_Rloc)

      class(zfunc_type) :: this

      !-- Output variables
      real(cp), intent(out) :: ur_Rloc(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp), intent(out) :: ut_Rloc(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp), intent(out) :: up_Rloc(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)

      !-- Local arrays
      real(cp) :: usr(n_phi_max_3D,n_r_max)
      real(cp) :: upp(n_phi_max_3D,n_r_max)
      real(cp) :: ekp(n_phi_max_3D,n_r_max)

      real(cp) :: upm(n_phi_max_3D,n_r_max)
      real(cp) :: uzm(n_phi_max_3D,n_r_max)
      !-- Local variables
      integer :: n_r, n_r_r, n_phi, n_th_NHS, n_th_SHS
      real(cp) :: s_r, z_r, z_eta
      real(cp) :: vs, vz, vrr, vth, vph
      real(cp) :: alpha_r1, alpha_r2
      real(cp) :: h_s, x_s, sfunc, hfunc, dhfunc, vpm, vzm

      call allgather_from_rloc(this%us_phys_Rloc,usr,n_phi_max_3D)
      call allgather_from_rloc(this%up_phys_Rloc,upp,n_phi_max_3D)

      if( l_ek_pump ) call allgather_from_rloc(this%ek_phys_Rloc,ekp,n_phi_max_3D)
      if( l_mag_pump ) then
         call allgather_from_rloc(this%up_magpump_Rloc,upm,n_phi_max_3D)
         call allgather_from_rloc(this%uz_magpump_Rloc,uzm,n_phi_max_3D)
      end if

      !do n_r=1,n_r_max
      !   do n_phi=1,n_phi_max_3D
      !      usr(n_phi,n_r) =(one-r(n_r))**2*r(n_r)**2
      !      upp(n_phi,n_r) =(one-r(n_r))**2*r(n_r)**2
      !      ekp(n_phi,n_r) = 1.0_cp
      !   end do
      !end do

      !-- Compute 3D velocity fields by a linear interpolation
      do n_r_r=nRstart3D,nRstop3D
         do n_th_NHS=1,n_theta_max/2
            n_th_SHS=n_theta_max+1-n_th_NHS
            s_r = r_3D(n_r_r)*sint(n_th_NHS)
            z_r = r_3D(n_r_r)*cost(n_th_NHS)
            if ( l_mag_pump ) then
               !-- Following the form of (Schaeffer, Silva, Pais, 2015)
               h_s = sqrt(r_cmb**2-s_r**2)
               x_s = z_r/h_s
               sfunc = 4.*(s_r - r_icb)*(r_cmb - s_r)
               hfunc = -(7./2.)*x_s*(1.-x_s)**2.*(1.+x_s)**2.
               dhfunc= -(7./2.)*(1./h_s)*( (1.-x_s)**2.*(1.+x_s)**2. + &
               &      2.*x_s*((1.-x_s)**2.*(1.+x_s) - (1.-x_s)*(1.+x_s)**2. ) )
            end if
            if ( s_r >= r_icb ) then !-- Outside TC
               n_r = 1
               do while ( r(n_r) >= s_r )! .and. n_r < n_r_max )
                  n_r = n_r+1
               end do
               alpha_r2 = (s_r-r(n_r))/(r(n_r-1)-r(n_r))
               alpha_r1 = one - alpha_r2
               z_eta = -s_r/(r_cmb*r_cmb-s_r*s_r)*z_r ! \beta * z
               do n_phi=1,n_phi_max_3D
                  vs = alpha_r1*usr(n_phi,n_r) + alpha_r2*usr(n_phi,n_r-1)
                  !-- vz = beta*z*vs
                  vz = z_eta*vs
                  if ( l_ek_pump ) then
                     !-- vz = beta*z*vs+ekpump
                     vz = vz + z_r*(alpha_r1*ekp(n_phi,n_r) + & 
                     &              alpha_r2*ekp(n_phi,n_r-1))
                  end if
                  vrr= vz*cost(n_th_NHS) + vs*sint(n_th_NHS)
                  vth= vs*cost(n_th_NHS) - vz*sint(n_th_NHS)
                  vph= alpha_r1*upp(n_phi,n_r) + alpha_r2*upp(n_phi,n_r-1)
                  if ( l_mag_pump ) then
                     vpm= alpha_r1*upm(n_phi,n_r) + alpha_r2*upm(n_phi,n_r-1)
                     vzm= alpha_r1*uzm(n_phi,n_r) + alpha_r2*uzm(n_phi,n_r-1)
                     vpm= sfunc*dhfunc*vpm!*dhfunc
                     vzm= sfunc* hfunc*vzm!*hfunc
                     !upmRloc(n_phi,n_th_NHS,n_r_r)= sfunc*dhfunc*vpm!*dhfunc
                     !upmRloc(n_phi,n_th_SHS,n_r_r)=-sfunc*dhfunc*vpm!*dhfunc
                     !uzmRloc(n_phi,n_th_NHS,n_r_r)= sfunc* hfunc*vzm!*hfunc
                     !uzmRloc(n_phi,n_th_SHS,n_r_r)=-sfunc* hfunc*vzm!*hfunc
                  else
                     vpm=0.0_cp
                     vzm=0.0_cp
                  end if
                  !if( rank == 0 .and. (n_r_r==1 .and. n_phi==1 .and. n_th_NHS==n_theta_max/2-1) ) &
                  !& print*, "mag_pump_ phi, z =", vpm, vzm
                  !if( rank == 0 .and. (n_r_r==1 .and. n_phi==1 .and. n_th_NHS==n_theta_max/2-1) ) &
                  !& print*, "UNmodified vel_ r, th, phi =", vrr, vth, vph
                  ur_Rloc(n_phi,n_th_NHS,n_r_r)= vrr + cost(n_th_NHS)*vzm
                  ur_Rloc(n_phi,n_th_SHS,n_r_r)= vrr - cost(n_th_NHS)*vzm
                  ut_Rloc(n_phi,n_th_NHS,n_r_r)= vth - sint(n_th_NHS)*vzm
                  ut_Rloc(n_phi,n_th_SHS,n_r_r)=-vth + sint(n_th_NHS)*vzm
                  up_Rloc(n_phi,n_th_NHS,n_r_r)= vph + vpm
                  up_Rloc(n_phi,n_th_SHS,n_r_r)= vph - vpm
                  !if( rank == 0 .and. (n_r_r==1 .and. n_phi==1 .and. n_th_NHS==n_theta_max/2-1) ) &
                  !& print*, "  MODified vel_ r, th, phi =", ur_Rloc(n_phi,n_th_NHS,n_r_r), &
                  !&         ut_Rloc(n_phi,n_th_NHS,n_r_r), up_Rloc(n_phi,n_th_NHS,n_r_r)
               end do

            else !-- Inside the tangent cylinder

               do n_phi=1,n_phi_max_3D
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

      !!print*, 'r_icb, what is wrong!!!!', r_icb, r_cmb
      !print*, 'eqjhn, what is wrong!!!!', rank, nRstart3D, nRstop3D

      !do n_r=nRstart3D,nRstop3D
      !   do n_th_NHS=1,n_theta_max
      !      n_th_SHS=n_theta_max+1-n_th_NHS
      !      s_r = r_3D(n_r)*sint(n_th_NHS)
      !!      ur_Rloc(:,n_th_NHS,n_r) = (one-r_3D(n_r))**2*r_3D(n_r)**2
      !!      ut_Rloc(:,n_th_NHS,n_r) = (one-r_3D(n_r))**2*r_3D(n_r)**2
      !!      up_Rloc(:,n_th_NHS,n_r) = (one-r_3D(n_r))**2*r_3D(n_r)**2
      !      if ( s_r < r_icb ) then !-- Inside TC
      !         if ( s_r > r_icb ) print*, 'n_th_NHS, s_r, r_icb', n_th_NHS, s_r, r_icb
      !         do n_phi=1,n_phi_max_3D
      !            ur_Rloc(n_phi,n_th_NHS,n_r)=0.0_cp
      !            ur_Rloc(n_phi,n_th_SHS,n_r)=0.0_cp
      !            ut_Rloc(n_phi,n_th_NHS,n_r)=0.0_cp
      !            ut_Rloc(n_phi,n_th_SHS,n_r)=0.0_cp
      !            up_Rloc(n_phi,n_th_NHS,n_r)=0.0_cp
      !            up_Rloc(n_phi,n_th_SHS,n_r)=0.0_cp
      !         end do!
      !      end if ! Inside/outside TC
      !   end do
      !end do

      !-- CMB values
      if ( nRstart3D == 1 ) then
         ur_Rloc(:,:,1) = 0.0_cp
         if ( ktopv == 2 ) then !-- Rigid boundaries
            ut_Rloc(:,:,1) = 0.0_cp
            up_Rloc(:,:,1) = 0.0_cp
         end if
      end if

#ifdef TOTO
      block
         integer :: n_r, file_handle

      if( rank == 0 ) then
         open(newunit=file_handle, file='extrapol_up_rloc.dat', status='new', form='formatted')
         write(file_handle, '(129es20.12)') 78374.0, theta(:n_theta_max)
         do n_r=1,n_r_max_3D
            write(file_handle, '(129es20.12)') r_3D(n_r), upmRloc(1,:n_theta_max,n_r)
         end do
         close(file_handle)

         open(newunit=file_handle, file='extrapol_uz_rloc.dat', status='new', form='formatted')
         write(file_handle, '(129es20.12)') 78374.0, theta(:n_theta_max)
         do n_r=1,n_r_max_3D
            write(file_handle, '(129es20.12)') r_3D(n_r), uzmRloc(1,:n_theta_max,n_r)
         end do
         close(file_handle)
      endif
      end block

      print*, 'ALL GOOD extrapolate!**'!

      stop
#endif

   end subroutine ext_QGto3D_vel
!--------------------------------------------------------------------------------
   subroutine compute_zavg(this,work_Rloc,zavg_Rloc)

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
      real(cp) :: work_Rext(n_theta_max,nRstart3D-2:nRstop3D+2)
      !real(cp) :: workall(n_theta_max,n_r_max_3D)
      real(cp) :: tmp(n_phi_max_3D,n_r_max)!nRstart:nRstop)!

      tmp(:,:)      =0.0_cp
      zavg_Rloc(:,:)=zero

      !-- z-average in spatial space
      if( l_cyl ) then !-- z-integration with a 4th-order scheme:: WARNING:: 3x slower!!
         do n_phi=1,n_phi_max_3D
            !call allgather_from_rloc_3D(work_Rloc(n_phi,:,:),workall,n_theta_max)
            call exchange_Nbound_from_Rloc_3D(work_Rloc(n_phi,:,:), work_Rext, 2, n_theta_max)
            call this%cyl_avg(work_Rext(:,:),tmp(n_phi,:))
         end do
      else !-- z-integration with a bilinear scheme
         do n_r=1,n_r_max
         !do n_r=nRstart,nRstop
            do n_z=1,4*n_z_max+1!this%nzp_zavg(n_r)!
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
      endif

         !-- TG: dirty fix: needs to be improved
         do n_r=1,n_r_max
            call MPI_Allreduce(MPI_IN_PLACE, tmp(:,n_r), n_phi_max_3D, MPI_DEF_REAL, &
                 &             MPI_SUM, MPI_COMM_WORLD, ierr)
         !call MPI_reduce(MPI_IN_PLACE, tmp(:,n_r), n_phi_max_3D, MPI_DEF_REAL, &
         !     &             MPI_SUM, MPI_COMM_WORLD, ierr)
         end do

#ifdef TOTO
      block
         integer :: n_r, file_handle

      if( rank == 0 ) then
         open(newunit=file_handle, file='zavg_rloc.dat', status='new', form='formatted')
         do n_r=1,n_r_max
            write(file_handle, '(2es20.12)') r(n_r), tmp(1,n_r)
         end do
         close(file_handle)
      endif
      end block

      if ( rank == 0) print*, 'ALL GOOD compute_avg!**'!

      stop
#endif

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

   end subroutine compute_zavg
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
      !real(cp) :: tmp(n_theta_max,n_r_max_3D,n_z_max)!nRstart3D:nRstop3D)!
      !real(cp) :: tmp(n_theta_max/2)
      !real(cp) :: Zwb(n_theta_max/2,0:n_procs-1)
      !-- Local variables
      integer :: n_th_NHS, n_z, n_z_r, n_z_t, n_th_SHS
      integer :: n_r!, n_t, n_r_r, n_t_t!, n_theta!, n_p
      real(cp) :: thwFac!, thwr, thwt, coefint
      !real(cp) :: r_i(n_r_max),r_i3D(n_r_max_3D)

      !open(1,file='./../../MagIC/test-build-out/urextra_com.dat')
      !do n_r=1,n_r_max
      !   read(1,*) r_i(n_r_max_3D+1-n_r)
      !enddo
      !close(1)
      !r_i3D(:) = r_i(:)

      thwFac=BuoFac/CorFac !Ra/Pr/(2/Ek) = RaEkPr/2

      call allgather_from_rloc_3D(dTdth_Rloc,dTdth,n_theta_max)

      !-- Remaining term for the temperature gradient
      dTzdt(:,:)=0.0_cp
      do n_r=1,n_r_max_3D!max(2,nRstart3D),nRstop3D
      !do n_r=max(2,nRstart3D),nRstop3D
         !do n_theta=1,n_theta_max
         do n_th_NHS=1,n_theta_max/2
            !dTzdt(n_th_NHS,n_r)=thwFac*r_3D(n_r)* dTdth(n_th_NHS,n_r)
            !-- TG I don't understand the r factor in the above equation
            !-- Th wind should be
            !-- duphi/dz = Ra/Pr  * g / r * dT/dtheta
            dTzdt(n_th_NHS,n_r)=thwFac*rgrav_3D(n_r)*or1_3D(n_r)* &
            !&                   cos(theta(n_th_NHS))
            &                   dTdth(n_th_NHS,n_r)
            !&                   dTdth_Rloc(n_th_NHS,n_r)
            !dTzdt(n_th_NHS,n_r)=(-r_3D(n_r)*cost(n_th_NHS))**2.!one!r_i3D(n_r)*theta(n_theta)
            !dTzdt(n_th_NHS,n_r)=exp(r_3D(n_r)*cost(n_th_NHS))!*sin(pi*r_3D(n_r)*cost(n_th_NHS))
         end do
      end do

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
         do n_th_NHS=2,n_theta_max/2!
            do n_z=1,this%nzp_thw(n_th_NHS,n_r)
               n_z_r = this%interp_zp_thw(1,n_z,n_th_NHS,n_r)
               n_z_t = this%interp_zp_thw(2,n_z,n_th_NHS,n_r)
               if ( n_z_t > 1 ) then
                  thw_Rloc(n_th_NHS,n_r)= thw_Rloc(n_th_NHS,n_r) -                &
                  &                      (this%interp_wt_thw(1,n_z,n_th_NHS,n_r)* &
                  & dTzdt(n_z_t,n_z_r) +  this%interp_wt_thw(2,n_z,n_th_NHS,n_r)* &
                  & dTzdt(n_z_t-1,n_z_r) )
               else !-- When n_z_t == 1; this%interp_wt_thw(2,n_z,n_theta,n_r)=0.0
                  thw_Rloc(n_th_NHS,n_r)=thw_Rloc(n_th_NHS,n_r) -               &
                  &                      this%interp_wt_thw(1,n_z,n_th_NHS,n_r)* &
                  & dTzdt(n_z_t,n_z_r)
               endif
            end do
         end do
      end do

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
   subroutine cyl_avg(this,work,cylavg)
      !-- Calculates on n_cyl_rad coaxial cylinders (outside
      !-- the tangentcylinder) the mean value of work(theta,r)
      !-- 
      !-- work: (input) array depending on phi,theta and r
      !-- cylavg: (output) mean value on cylinder surface

      class(zfunc_type) :: this

      !-- Input variables
      real(cp), intent(in) :: work(n_theta_max,nRStart3D-2:nRStop3D+2)
      !real(cp), intent(in) :: work(n_theta_max,n_r_max_3D)

      !-- Output variable
      real(cp), intent(out) :: cylavg(n_r_max)

      !-- Local variables
      integer :: n_z, nz, n_s, n_th, itr
      integer :: n_r0, n_r1, n_r2, n_r3, n_th0, n_th1, n_th2, n_th3
      real(cp) :: rr0, rr1, rr2, rr3, r10, r20, r30, r21, r31, r32
      real(cp) :: tt0, tt1, tt2, tt3, t10, t20, t30, t21, t31, t32
      real(cp) :: a01, a12, a23, a012, a123, tot
      real(cp) :: acyl(-n_r_max:n_r_max,n_r_max), aith(0:3)

      cylavg(:)=0.0_cp
      !-- Loop over axial cylinders starts here
      sLoop: do n_s=1,n_r_max

         nz = this%nzp_zcyl(n_s)
         !-- Loop over z starts
         do n_z=-nz,nz
            acyl(n_z,n_s)=0.0_cp
            !
            !  **** Interpolate values from (theta,r)-grid onto equidistant
            !  **** (z,rax)-grid using a fourth-order Lagrangian scheme
            !
            n_r2=this%interp_zr_cyl(n_z,n_s)
            n_r3=n_r2-1
            n_r1=n_r2+1
            n_r0=n_r2+2

            n_th1=this%interp_zt_cyl(n_z,n_s)
            n_th2=n_th1+1
            n_th3=n_th1+2
            n_th0=n_th1-1

            if ( n_r2 >= nRstart3D .and. n_r2 <= nRstop3D ) then

            !--  Calculate differences in r for 4th-order interpolation
            rr0=this%interp_wr_cyl(0,n_z,n_s)! rc-r_3D(n_r0)
            rr1=this%interp_wr_cyl(1,n_z,n_s)! rc-r_3D(n_r1)
            rr2=this%interp_wr_cyl(2,n_z,n_s)! rc-r_3D(n_r2)
            rr3=this%interp_wr_cyl(3,n_z,n_s)! rc-r_3D(n_r3)
            r10=this%interp_wr_cyl(4,n_z,n_s)! one/(r_3D(n_r1)-r_3D(n_r0))
            r20=this%interp_wr_cyl(5,n_z,n_s)! one/(r_3D(n_r2)-r_3D(n_r0))
            r30=this%interp_wr_cyl(6,n_z,n_s)! one/(r_3D(n_r3)-r_3D(n_r0))
            r21=this%interp_wr_cyl(7,n_z,n_s)! one/(r_3D(n_r2)-r_3D(n_r1))
            r31=this%interp_wr_cyl(8,n_z,n_s)! one/(r_3D(n_r3)-r_3D(n_r1))
            r32=this%interp_wr_cyl(9,n_z,n_s)! one/(r_3D(n_r3)-r_3D(n_r2))

            !--  Calculate differences in theta for 4th-order interpolation
            tt0=this%interp_wt_cyl(0,n_z,n_s)! thet-theta(n_th0)
            tt1=this%interp_wt_cyl(1,n_z,n_s)! thet-theta(n_th1)
            tt2=this%interp_wt_cyl(2,n_z,n_s)! thet-theta(n_th2)
            tt3=this%interp_wt_cyl(3,n_z,n_s)! thet-theta(n_th3)
            t10=this%interp_wt_cyl(4,n_z,n_s)! one/(theta(n_th1)-theta(n_th0))
            t20=this%interp_wt_cyl(5,n_z,n_s)! one/(theta(n_th2)-theta(n_th0))
            t30=this%interp_wt_cyl(6,n_z,n_s)! one/(theta(n_th3)-theta(n_th0))
            t21=this%interp_wt_cyl(7,n_z,n_s)! one/(theta(n_th2)-theta(n_th1))
            t31=this%interp_wt_cyl(8,n_z,n_s)! one/(theta(n_th3)-theta(n_th1))
            t32=this%interp_wt_cyl(9,n_z,n_s)! one/(theta(n_th3)-theta(n_th2))

            !-- Loop over 4 neighboring grid angles
            do itr=0,3
               n_th=n_th0+itr
               !-- Interpolation in r-direction
               a01=(rr0*work(n_th,n_r1) -    &
               &    rr1*work(n_th,n_r0))*r10
               a12=(rr1*work(n_th,n_r2) -    &
               &    rr2*work(n_th,n_r1))*r21
               a23=(rr2*work(n_th,n_r3) -    &
               &    rr3*work(n_th,n_r2))*r32

               a012=(rr0*a12-rr2*a01)*r20
               a123=(rr1*a23-rr3*a12)*r31

               aith(itr)=(rr0*a123-rr3*a012)*r30
            end do

            !-- Interpolation in theta-direction
            a01=(tt0*aith(1)-tt1*aith(0))*t10
            a12=(tt1*aith(2)-tt2*aith(1))*t21
            a23=(tt2*aith(3)-tt3*aith(2))*t32
        
            a012=(tt0*a12-tt2*a01)*t20
            a123=(tt1*a23-tt3*a12)*t31
            acyl(n_z,n_s)=acyl(n_z,n_s)+(tt0*a123-tt3*a012)*t30
            end if
         end do !-- end z-loop
         !
         !  *** interpolation completed
         !
         !  *** simpson integration
         !
         tot=acyl(-nz,n_s)+acyl(nz,n_s)
         do n_z=-nz+1,nz-1,2
            tot=tot+four*acyl(n_z,n_s)
         enddo
         do n_z=-nz+2,nz-2,2
            tot=tot+two*acyl(n_z,n_s)
         enddo
         cylavg(n_s)=tot/(6.0_cp*nz) !-- f_tot/2h *D/6 with D = h/N_pts
      end do sLoop

      !!--  special case s=rmax
      !cylavg(0) = half*( work(n_theta_max/2,1)  + &
      !&                  work(n_theta_max/2+1,1) )

   end subroutine cyl_avg
!--------------------------------------------------------------------------------
   subroutine compute_zder(this,work_Rloc,zder_Rloc)

      class(zfunc_type) :: this

      !-- Input variables
      real(cp), intent(in) :: work_Rloc(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)

      !-- Output variables
      real(cp), intent(out) :: zder_Rloc(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)

      !-- Local variables
      integer :: n_r, n_phi, n_th
      real(cp) :: work_Rextended(n_theta_max,nRstart3D-1:nRstop3D+1)

      zder_Rloc(:,:,:)=0.0_cp

      !-- z-derivation on the 3D grid
      !-- has to use finite difference::  dz = cos(theta)*dr - sin(theta)/r_3D*dtheta
      !--     but the precision is not a problem as zavg scheme order is =1.8

      !-- r-derivation first: dz = dr
      do n_phi=1,n_phi_max_3D
         call exchange_Nbound_from_Rloc_3D(work_Rloc(n_phi,:,:),           &
              &                            work_Rextended, 1, n_theta_max)
         do n_r=nRstart3D,nRstop3D
            if ( n_r == 1 ) then
               zder_Rloc(n_phi,:,n_r)=(work_Rextended(:,n_r+1) -           &
               &         work_Rextended(:,n_r))/(r_3D(n_r+1) - r_3D(n_r))
            else
               zder_Rloc(n_phi,:,n_r)=(work_Rextended(:,n_r) -             &
               &         work_Rextended(:,n_r-1))/(r_3D(n_r) - r_3D(n_r-1))
            endif
         end do
      end do

      !-- theta-derivation then: dz = cost*dz - sint*or1*dt
      do n_r=nRstart3D,nRstop3D
         do n_th=1,n_theta_max
            if ( n_th == 1 ) then
               zder_Rloc(:,n_th,n_r)=cost(n_th)*zder_Rloc(:,n_th,n_r) - &
               &        sint(n_th)*or1_3D(n_r)*(work_Rloc(:,n_th,n_r) - &
               &  work_Rloc(:,n_th+1,n_r))/(theta(n_th) - theta(n_th+1))
            else if ( n_th == n_theta_max ) then
               zder_Rloc(:,n_th,n_r)=cost(n_th)*zder_Rloc(:,n_th,n_r) -   &
               &        sint(n_th)*or1_3D(n_r)*(work_Rloc(:,n_th-1,n_r) - &
               &      work_Rloc(:,n_th,n_r))/(theta(n_th-1) - theta(n_th))
            else
               zder_Rloc(:,n_th,n_r)=cost(n_th)*zder_Rloc(:,n_th,n_r) -   &
               &        sint(n_th)*or1_3D(n_r)*(work_Rloc(:,n_th-1,n_r) - &
               &  work_Rloc(:,n_th+1,n_r))/(theta(n_th-1) - theta(n_th+1))
            endif
         end do
      end do

#ifdef TOTO
      block
         integer :: n_r, file_handle

      if( rank == 0 ) then
         open(newunit=file_handle, file='zder_rloc.dat', status='new', form='formatted')
         write(file_handle, '(257es20.12)') 78374.0, theta(:n_theta_max/2)
         do n_r=1,n_r_max_3D
            write(file_handle, '(257es20.12)') r_3D(n_r), zder_Rloc(1,:n_theta_max/2,n_r)
         end do
         close(file_handle)
      endif
      end block

      print*, 'Radius OUTPUT.dat', r_3D(n_r_max_3D/2), 'Latitude OUTPUT.dat', theta(n_theta_max/4)
      print*, 'ALL GOOD z_derivative!**'

      stop
#endif

   end subroutine compute_zder
!--------------------------------------------------------------------------------
   subroutine fill_zinterp_grid(this)

      !-- Input variables
      class(zfunc_type) :: this

      !-- Local arrays
      real(cp) :: zz(0:n_r_max_3D)!(nRstop3D-nRstart3D)+1)

      !-- Local variables
      integer :: n_r, n_t, n_z, n_r_r, n_t_t
      integer :: n_start!, n_p
      integer :: nz, n_s, n_th
      integer :: n_r0, n_r1, n_r2, n_r3, n_th0, n_th1, n_th2, n_th3
      real(cp) :: r_r, z_r, c_t, h
      real(cp) :: norm, alpha_r, alpha_t
      real(cp) :: s_r, dz, th
      real(cp) :: zmin, zmax, z, eps, rc, thet

#ifdef TOTO
      !-- Get z grid
      !-- for interpolation: n_r_max points on z-axis
      norm = one/(two*n_z_max+1)
      do n_r=1,n_r_max!-1! loop on all s
      !do n_r=nRstart,nRstop
         h = sqrt(r_cmb**2-r(n_r)**2)!half*height(n_r)! WARNING:: beta_shift in height!
         n_r_r = 2
         n_t_t = 2
         this%nzp_zavg(n_r)=0
         do n_z=n_z_max+1,1,-1
            z_r = h*(n_z-1)/n_z_max
            r_r = sqrt(r(n_r)**2 + z_r**2)
            c_t = z_r/r_r
            do while ( r_3D(n_r_r) > r_r  )
               n_r_r = n_r_r+1
            end do
            do while ( cost(n_t_t) > c_t  )
               n_t_t = n_t_t+1
            end do
            !-- Compute coeffs for bilinear interpolation
            alpha_r = (r_r-r_3D(n_r_r))/(r_3D(n_r_r-1)-r_3D(n_r_r))
            alpha_t = (c_t-cost(n_t_t))/(cost(n_t_t-1)-cost(n_t_t))
            !-- Coordinates neighbourg 1: n_r_r, n_t_t
            call add_new_point(this, n_r, n_r_r, n_t_t, &
                 &                   norm*(one-alpha_r)*(one-alpha_t))
            !-- Coordinates neighbourg 2: n_r_r, n_t_t-1
            call add_new_point(this, n_r, n_r_r, n_t_t-1, &
                 &                   norm*(one-alpha_r)*alpha_t)
            !-- Coordinates neighbourg 3: n_r_r-1, n_t_t
            call add_new_point(this, n_r, n_r_r-1, n_t_t, &
                 &                   norm*alpha_r*(one-alpha_t))
            !-- Coordinates neighbourg 4: n_r_r-1, n_t_t-1
            call add_new_point(this, n_r, n_r_r-1, n_t_t-1, &
                 &                   norm*alpha_r*alpha_t)
         end do
      end do
#endif
!#ifdef TOTO
      !-- Get z grid
      !-- for interpolation: n_r_max points on z-axis
      norm = one/(two*n_z_max+1)
      do n_r=1,n_r_max!-1 loop on all s
         h = sqrt(r_cmb**2-r(n_r)**2)!half*height(n_r)! WARNING:: beta_shift in height!
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
            !   norm = one/(two*n_z_max+1)
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
!#endif

      if ( l_cyl ) then
         !-- Get the z-avg interpolator 4th order method
         eps=10.0_cp*epsilon(one)
         !-- Loop over axial cylinders starts here
         sLoop: do n_s=1,n_r_max
            if ( r(n_s) < r_icb ) exit sLoop

            zmax = sqrt(r_cmb*r_cmb-r(n_s)*r(n_s)) ! zmax
            zmin = 0.0_cp
            nz = 2*int(n_r_max*(zmax-zmin)/(two*r_cmb)) ! Number of z points (one HS)
            nz = max(nz, 4)  ! Minimum to 4 for Simpson integration
            this%nzp_zcyl(n_s)=nz
            dz = (zmax-zmin)/real(nz,cp)
            !-- Loop over z starts
            do n_z=-nz,nz
               z=zmin+dz*n_z
               rc=sqrt(r(n_s)*r(n_s)+z*z)   ! radius from center
               if (rc >= r_cmb) rc=r_cmb-eps
               thet=half*pi-atan(z/r(n_s))  ! polar angle of point (rax,z)
               !
               !  **** Interpolate values from (theta,r)-grid onto equidistant
               !  **** (z,rax)-grid using a fourth-order Lagrangian scheme
               !
               !--  Find indices of radial grid levels that bracket rc
               rbracket: do n_r=n_r_max_3D-1,1,-1
                  if ( r_3D(n_r) >= rc ) then
                     n_r2 = n_r
                     exit rbracket
                  end if
               end do rbracket
               if(n_r2 == n_r_max_3D-1) n_r2=n_r_max_3D-2
               if(n_r2 == 1 ) n_r2=2
               this%interp_zr_cyl(n_z,n_s)=n_r2
               n_r3=n_r2-1
               n_r1=n_r2+1
               n_r0=n_r2+2

               !-- Find indices of angular grid levels that bracket thet
               tbracket: do n_th=n_theta_max,1,-1
                  if( theta(n_th) <= thet) then
                     n_th1=n_th
                     exit tbracket
                  end if
               end do tbracket
               if ( n_th1 == n_theta_max ) n_th1=n_theta_max-2
               if ( n_th1 == n_theta_max-1 ) n_th1=n_theta_max-2
               if ( n_th1 == 1 ) n_th1=2
               this%interp_zt_cyl(n_z,n_s)=n_th1
               n_th2=n_th1+1
               n_th3=n_th1+2
               n_th0=n_th1-1

               !--  Calculate differences in r for 4th-order interpolation
               this%interp_wr_cyl(0,n_z,n_s)=rc-r_3D(n_r0) !rr0
               this%interp_wr_cyl(1,n_z,n_s)=rc-r_3D(n_r1) !rr1
               this%interp_wr_cyl(2,n_z,n_s)=rc-r_3D(n_r2) !rr2
               this%interp_wr_cyl(3,n_z,n_s)=rc-r_3D(n_r3) !rr3
               this%interp_wr_cyl(4,n_z,n_s)= one/(r_3D(n_r1)-r_3D(n_r0)) !r10
               this%interp_wr_cyl(5,n_z,n_s)= one/(r_3D(n_r2)-r_3D(n_r0)) !r20
               this%interp_wr_cyl(6,n_z,n_s)= one/(r_3D(n_r3)-r_3D(n_r0)) !r30
               this%interp_wr_cyl(7,n_z,n_s)= one/(r_3D(n_r2)-r_3D(n_r1)) !r21
               this%interp_wr_cyl(8,n_z,n_s)= one/(r_3D(n_r3)-r_3D(n_r1)) !r31
               this%interp_wr_cyl(9,n_z,n_s)= one/(r_3D(n_r3)-r_3D(n_r2)) !r32

               !--  Calculate differences in theta for 4th-order interpolation
               this%interp_wt_cyl(0,n_z,n_s)=thet-theta(n_th0) !tt0
               this%interp_wt_cyl(1,n_z,n_s)=thet-theta(n_th1) !tt1
               this%interp_wt_cyl(2,n_z,n_s)=thet-theta(n_th2) !tt2
               this%interp_wt_cyl(3,n_z,n_s)=thet-theta(n_th3) !tt3
               this%interp_wt_cyl(4,n_z,n_s)=one/(theta(n_th1)-theta(n_th0)) !t10
               this%interp_wt_cyl(5,n_z,n_s)=one/(theta(n_th2)-theta(n_th0)) !t20
               this%interp_wt_cyl(6,n_z,n_s)=one/(theta(n_th3)-theta(n_th0)) !t30
               this%interp_wt_cyl(7,n_z,n_s)=one/(theta(n_th2)-theta(n_th1)) !t21
               this%interp_wt_cyl(8,n_z,n_s)=one/(theta(n_th3)-theta(n_th1)) !t31
               this%interp_wt_cyl(9,n_z,n_s)=one/(theta(n_th3)-theta(n_th2)) !t32
            end do !-- end z-loop
         end do sLoop
      end if

      if ( l_heat_3D .and. l_thw_3D ) then
         !-- Get theta weights for thermal wind calculation
         zz(:) = 0.0_cp
         !do n_r=2,n_r_max_3D
         !do n_r=max(1,nRstart3D-1),min(nRstop3D+1,n_r_max_3D)
         do n_r=max(2,nRstart3D),nRstop3D
            do n_t=1,n_theta_max/2+1
               s_r = r_3D(n_r)*sint(n_t)
               n_z = 0
               n_start = 2
               do n_r_r=n_start,n_r!-1 !-- Flip the loop since our radii are decreasing.
                  if( n_r_r > n_r_max_3D) print*, 'w!! segFault in n_r_r loop!; n_r_r', n_r_r
                  n_z = n_z+1
                  th  = asin(s_r*or1_3D(n_r_r))!/r_i3D(n_r_r))!
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
                     &         th<=theta(n_t_t)) .and. n_t_t < n_theta_max ) 
                        n_t_t=n_t_t+1
                     end do
                     if( n_t_t > n_theta_max ) print*, 'w!! segFault in n_t_t loop!; n_t_t', n_t_t
                     if ( n_r_r==n_r ) then
                        this%interp_zp_thw(1,n_z,n_t,n_r)=n_r
                        this%interp_zp_thw(2,n_z,n_t,n_r)=n_t
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
      end if

   end subroutine fill_zinterp_grid
!--------------------------------------------------------------------------------
   subroutine add_new_point(this, n_s, n_r, n_t, weight)

      !-- Input variables
      class(zfunc_type) :: this
      integer :: n_s, n_r, n_t
      real(cp) :: weight

      !-- Local variables
      integer :: n_p

      !-- Add a new point for the z-integral-interpolation
      do n_p=1,this%nzp_zavg(n_s)! number of points already stored
         !-- Scan to see if point is already used ...
         if ( n_r == this%interp_zr_mat(n_p,n_s) .and.  &
         &    n_t == this%interp_zt_mat(n_p,n_s) ) then
            this%interp_wt_mat(n_p,n_s) = this%  &
            &    interp_wt_mat(n_p,n_s) + weight
            return
         end if
      end do
      !-- Points is not stored yet: needs a new one!
      n_p=this%nzp_zavg(n_s)+1
      this%interp_zr_mat(n_p,n_s)=n_r
      this%interp_zt_mat(n_p,n_s)=n_t
      this%interp_wt_mat(n_p,n_s)=weight
      this%nzp_zavg(n_s)=n_p

   end subroutine add_new_point
!--------------------------------------------------------------------------------
end module z_functions
