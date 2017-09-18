module outputs

   use parallel_mod
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use namelists, only: tag, ra, pr, l_non_rot, l_vphi_balance, ek, &
       &                radratio, raxi, sc, tadvz_fac, kbotv, ktopv
   use communications, only: reduce_radial_on_rank
   use truncation, only: n_r_max, m2idx, n_m_max, idx2m
   use radial_functions, only: r, rscheme, rgrav, dtcond, height, tcond, &
       &                       beta, ekpump, or1
   use blocking, only: nMstart, nMstop, nRstart, nRstop, l_rank_has_m0, &
       &               nm_per_rank, m_balance
   use integration, only: rInt_R
   use useful, only: round_off, cc2real, cc22real, getMSD2, abortRun
   use constants, only: pi, two, four, surf, vol_otc, one
   use checkpoints, only: write_checkpoint
   use output_frames, only: write_snapshot_mloc
   use time_scheme, only: type_tscheme

   implicit none

   private

   integer :: frame_counter, n_calls, spec_counter
   real(cp) :: timeLast, timeAvg
   integer, public :: n_log_file
   integer :: n_rey_file, n_heat_file, n_kin_file, n_power_file
   integer :: n_vphi_bal_file
   character(len=144), public :: log_file
   character(len=144) :: kin_file, heat_file, power_file, rey_file

   real(cp), allocatable :: uphiR_mean(:), uphiR_SD(:)
   real(cp), allocatable :: tempR_mean(:), tempR_SD(:)

   type, public :: vp_bal_type
      real(cp), allocatable :: rey_stress(:)
      real(cp), allocatable :: dvpdt(:)
      real(cp), allocatable :: visc(:)
      real(cp), allocatable :: pump(:)
      integer :: n_calls
   end type vp_bal_type

   type(vp_bal_type), public :: vp_bal

   public :: initialize_outputs, finalize_outputs, get_time_series, &
   &         write_outputs

contains

   subroutine initialize_outputs

      logical :: log_exists

      if ( rank == 0 ) then
         log_file = 'log.'//tag
         inquire(file=log_file, exist=log_exists)
         if ( log_exists ) then
            call abortRun('! log file already exists, please change tag')
         end if
         open(newunit=n_log_file, file=log_file, status='new')
         kin_file = 'e_kin.'//tag
         open(newunit=n_kin_file, file=kin_file, status='new')
         power_file = 'power.'//tag
         open(newunit=n_power_file, file=power_file, status='new')
         rey_file = 'reynolds.'//tag
         open(newunit=n_rey_file,file=rey_file, status='new')
      end if 

      timeAvg = 0.0_cp

      if ( l_rank_has_m0 ) then
         heat_file = 'heat.'//tag
         open(newunit=n_heat_file, file=heat_file, status='new')

         allocate( uphiR_mean(n_r_max), uphiR_SD(n_r_max) )
         allocate( tempR_mean(n_r_max), tempR_SD(n_r_max) )
         bytes_allocated=bytes_allocated+4*n_r_max*SIZEOF_DEF_REAL

         uphiR_mean(:) = 0.0_cp
         uphiR_SD(:)   = 0.0_cp
         tempR_mean(:) = 0.0_cp
         tempR_SD(:)   = 0.0_cp
         n_calls       = 0
         timeLast      = 0.0_cp

         if ( l_vphi_balance ) then
            open(newunit=n_vphi_bal_file, file='vphi_bal.'//tag, &
            &    form='unformatted', status='new')

            allocate( vp_bal%rey_stress(n_r_max) )
            allocate( vp_bal%dvpdt(n_r_max) )
            allocate( vp_bal%visc(n_r_max) )
            allocate( vp_bal%pump(n_r_max) )
            vp_bal%n_calls = 0
            bytes_allocated = bytes_allocated+4*n_r_max*SIZEOF_DEF_REAL
         end if
      end if

      frame_counter = 1 ! For file suffix
      spec_counter = 1

   end subroutine initialize_outputs
!------------------------------------------------------------------------------
   subroutine finalize_outputs

      if ( l_rank_has_m0 ) then
         if ( l_vphi_balance ) then
            close(n_vphi_bal_file)
            deallocate( vp_bal%pump, vp_bal%visc )
            deallocate( vp_bal%dvpdt, vp_bal%rey_stress )
         end if
         deallocate( uphiR_mean, uphiR_SD, tempR_mean, tempR_SD )
         close(n_heat_file)
      end if

      if ( rank == 0 ) then
         close(n_rey_file)
         close(n_power_file)
         close(n_kin_file)
         close(n_log_file)
      end if

   end subroutine finalize_outputs
!------------------------------------------------------------------------------
   subroutine write_outputs(time, tscheme, n_time_step, l_log, l_rst,        &
              &             l_spec, l_frame, l_vphi_bal_write, l_stop_time,  &
              &             us_Mloc, up_Mloc, om_Mloc, temp_Mloc, dtemp_Mloc,&
              &             dtemp_exp_Mloc, dpsi_exp_Mloc)

      !-- Input variables
      real(cp),           intent(in) :: time
      type(type_tscheme), intent(in) :: tscheme
      integer,            intent(in) :: n_time_step
      logical,            intent(in) :: l_log
      logical,            intent(in) :: l_rst
      logical,            intent(in) :: l_spec
      logical,            intent(in) :: l_frame
      logical,            intent(in) :: l_vphi_bal_write
      logical,            intent(in) :: l_stop_time
      complex(cp),        intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),        intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),        intent(in) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),        intent(in) :: temp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),        intent(in) :: dtemp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),        intent(in) :: dtemp_exp_Mloc(tscheme%norder_exp,nMstart:nMstop,n_r_max)
      complex(cp),        intent(in) :: dpsi_exp_Mloc(tscheme%norder_exp,nMstart:nMstop,n_r_max)

      !-- Local variable
      character(len=144) :: frame_name
      character(len=144) :: spec_name

      timeAvg = timeAvg + tscheme%dt(1)

      if ( l_rst ) then
         call write_checkpoint(time, tscheme, n_time_step, n_log_file, &
              &                l_stop_time, temp_Mloc, us_Mloc, up_Mloc, &
              &                dtemp_exp_Mloc, dpsi_exp_Mloc)
      end if

      if ( l_spec ) then
         write(spec_name, '(A,I0,A,A)') 'spec_',spec_counter,'.',tag
         call write_spectra(spec_name, us_Mloc, up_Mloc, om_Mloc, temp_Mloc)
         spec_counter = spec_counter+1
      end if

      if ( l_vphi_bal_write ) then
         vp_bal%n_calls = vp_bal%n_calls+1
         call write_vphi_balance(time, up_Mloc)
      end if

      if ( l_frame ) then
         write(frame_name, '(A,I0,A,A)') 'frame_temp_',frame_counter,'.',tag
         call write_snapshot_mloc(frame_name, time, temp_Mloc)
         write(frame_name, '(A,I0,A,A)') 'frame_us_',frame_counter,'.',tag
         call write_snapshot_mloc(frame_name, time, us_Mloc)
         write(frame_name, '(A,I0,A,A)') 'frame_up_',frame_counter,'.',tag
         call write_snapshot_mloc(frame_name, time, up_Mloc)
         write(frame_name, '(A,I0,A,A)') 'frame_om_',frame_counter,'.',tag
         call write_snapshot_mloc(frame_name, time, om_Mloc)
         frame_counter = frame_counter+1
      end if

      if ( l_log ) then
         call get_time_series(time, us_Mloc, up_Mloc, om_Mloc, temp_Mloc, &
              &               dtemp_Mloc)

         call get_radial_averages(timeAvg, l_stop_time, up_Mloc, temp_Mloc)
      end if

   end subroutine write_outputs
!------------------------------------------------------------------------------
   subroutine get_radial_averages(timeAvg, l_stop_time, up_Mloc, temp_Mloc)

      !-- Input variables
      real(cp),    intent(in) :: timeAvg
      logical,     intent(in) :: l_stop_time
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(in) :: temp_Mloc(nMstart:nMstop, n_r_max)

      !-- Local variables
      real(cp) :: dtAvg
      integer :: n_r, idx, file_handle


      if ( l_rank_has_m0 ) then

         n_calls = n_calls+1
         dtAvg = timeAvg-timeLast

         idx = m2idx(0)
         do n_r=1,n_r_max
            call getMSD2(uphiR_mean(n_r), uphiR_SD(n_r), real(up_Mloc(idx,n_r)),&
                 &       n_calls, dtAvg, timeAvg)
            call getMSD2(tempR_mean(n_r), tempR_SD(n_r), real(temp_Mloc(idx,n_r)),&
                 &       n_calls, dtAvg, timeAvg)
         end do
         timeLast = timeAvg

         if ( l_stop_time ) then
            open(newunit=file_handle, file='radial_profiles.'//tag)
            do n_r=1,n_r_max
               uphiR_SD(n_r)=sqrt(uphiR_SD(n_r)/timeAvg)
               tempR_SD(n_r)=sqrt(tempR_SD(n_r)/timeAvg)
               write(file_handle, '(es20.12, 4es16.8)') r(n_r),            &
               &     round_off(uphiR_mean(n_r)), round_off(uphiR_SD(n_r)), &
               &     round_off(tempR_mean(n_r)+tcond(n_r)),                &
               &     round_off(tempR_SD(n_r))
            end do
            close(file_handle)
         end if

      end if

   end subroutine get_radial_averages
!------------------------------------------------------------------------------
   subroutine get_time_series(time, us_Mloc, up_Mloc, om_Mloc, temp_Mloc, &
              &               dtemp_Mloc)

      !-- Input variables
      real(cp),    intent(in) :: time
      complex(cp), intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: temp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: dtemp_Mloc(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_m, m, n_m0
      real(cp) :: tTop, tBot, visc, pow, pum
      real(cp) :: us2, up2, up2_axi
      real(cp) :: us2_r(n_r_max), up2_r(n_r_max)
      real(cp) :: up2_axi_r(n_r_max)
      real(cp) :: visc_diss(n_r_max), buo_power(n_r_max), pump(n_r_max)
      real(cp) :: NuTop, NuBot, beta_t, vol_fac
      real(cp) :: rey, rey_fluct, rey_zon

      !-- This is to determine the volume used to compute the Reynolds numbers
      if ( l_non_rot ) then
         vol_fac = surf
      else ! Rotating
         vol_fac = vol_otc
      end if

      do n_r=1,n_r_max
         us2_r(n_r)    =0.0_cp
         up2_r(n_r)    =0.0_cp
         up2_axi_r(n_r)=0.0_cp
         visc_diss(n_r)=0.0_cp
         buo_power(n_r)=0.0_cp
         pump(n_r)     =0.0_cp
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            us2_r(n_r)    =us2_r(n_r)+cc2real(us_Mloc(n_m,n_r),m)
            up2_r(n_r)    =up2_r(n_r)+cc2real(up_Mloc(n_m,n_r),m)
            visc_diss(n_r)=visc_diss(n_r)+cc2real(om_Mloc(n_m,n_r),m)
            buo_power(n_r)=buo_power(n_r)+cc22real(us_Mloc(n_m,n_r), &
            &              temp_Mloc(n_m,n_r), m)

            if ( m == 0 ) then
               up2_axi_r(n_r)=up2_axi_r(n_r)+cc2real(up_Mloc(n_m,n_r),m)
               pump(n_r)     =pump(n_r)+cc2real(up_Mloc(n_m,n_r),m)
            end if
         end do
         us2_r(n_r)    =us2_r(n_r)*r(n_r)*height(n_r)
         up2_r(n_r)    =up2_r(n_r)*r(n_r)*height(n_r)
         up2_axi_r(n_r)=up2_axi_r(n_r)*r(n_r)*height(n_r)
         visc_diss(n_r)=visc_diss(n_r)*r(n_r)*height(n_r)
         buo_power(n_r)=Ra/pr*buo_power(n_r)*rgrav(n_r)*r(n_r)*height(n_r)
         pump(n_r)     =pump(n_r)*ekpump(n_r)*height(n_r)*r(n_r)
      end do

      !-- MPI reductions to get the s-profiles on rank==0
      call reduce_radial_on_rank(us2_r, 0)
      call reduce_radial_on_rank(up2_r, 0)
      call reduce_radial_on_rank(up2_axi_r, 0)
      call reduce_radial_on_rank(visc_diss, 0)
      call reduce_radial_on_rank(buo_power, 0)
      call reduce_radial_on_rank(pump, 0)

      if ( rank == 0 ) then
         !-- Kinetic energy
         us2 = rInt_R(us2_r, r, rscheme)
         up2 = rInt_R(up2_r, r, rscheme)
         up2_axi = rInt_R(up2_axi_r, r, rscheme)
         us2 = round_off(pi*us2)
         up2 = round_off(pi*up2)
         up2_axi = round_off(pi*up2_axi)

         !-- Reynolds numbers
         rey       = sqrt(two*(us2+up2)/vol_fac)
         rey_zon   = sqrt(two*(up2_axi)/vol_fac)
         rey_fluct = sqrt(two*(us2+up2-up2_axi)/vol_fac)

         !-- Total viscous dissipation
         visc = rInt_R(visc_diss, r, rscheme)
         visc = round_off(two*pi*visc)

         !-- In case of stress-free we need to include the surface contributions
         if ( kbotv == 1 ) then
            visc = visc-four*or1(n_r_max)*up2_r(n_r_max)
         end if
         if ( ktopv == 1 ) then
            visc = visc-four*or1(1)*up2_r(1)
         end if

         !-- Total buoyancy power
         pow  = rInt_R(buo_power, r, rscheme)
         pow  = round_off(two*pi*pow)

         pum  = rInt_R(pump, r, rscheme)
         pum  = round_off(two*pi*pum)

         write(n_kin_file, '(1P, es20.12, 3es16.8)') time, us2, up2, &
         &                                           up2_axi
         write(n_power_file, '(1P, es20.12, 3es16.8)') time, pow, visc, pum

         write(n_rey_file, '(1P, es20.12, 3es16.8)') time, rey, rey_zon, &
         &                                           rey_fluct
      end if

      if ( l_rank_has_m0 ) then
         !-- Top and bottom temperatures
         n_m0 = m2idx(0)
         tTop = real(temp_Mloc(n_m0,1))+tcond(1)
         tBot = real(temp_Mloc(n_m0,n_r_max))+tcond(n_r_max)
         tTop = round_off(tTop)
         tBot = round_off(tBot)

         NuTop = one+real(dtemp_Mloc(n_m0,1))/dtcond(1)
         NuBot = one+real(dtemp_Mloc(n_m0,n_r_max))/dtcond(n_r_max)
         !&       (dtcond(n_r_max)-tadvz_fac*beta(n_r_max)*tcond(n_r_max))
         NuTop = round_off(NuTop)
         NuBot = round_off(NuBot)

         beta_t = dtcond(int(n_r_max/2))+real(dtemp_Mloc(n_m0,int(n_r_max/2)))
         beta_t = round_off(beta_t)

         write(n_heat_file, '(1P, ES20.12, 5ES16.8)') time, NuTop, NuBot, &
         &                                             tTop, tBot, beta_t
      end if

   end subroutine get_time_series
!------------------------------------------------------------------------------
   subroutine write_vphi_balance(time, up_Mloc)

      !-- Input variables
      real(cp),    intent(in) :: time
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)

      !-- Local variable
      integer :: idx_m0

      if ( l_rank_has_m0 ) then
         idx_m0 = m2idx(0)
         if ( vp_bal%n_calls == 1 ) then
            write(n_vphi_bal_file) ra, ek, pr, radratio, raxi, sc
            write(n_vphi_bal_file) r
         end if
         write(n_vphi_bal_file) time
         write(n_vphi_bal_file) real(up_Mloc(idx_m0,:))
         write(n_vphi_bal_file) vp_bal%dvpdt
         write(n_vphi_bal_file) vp_bal%rey_stress
         write(n_vphi_bal_file) vp_bal%pump
         write(n_vphi_bal_file) vp_bal%visc
      end if

   end subroutine write_vphi_balance
!------------------------------------------------------------------------------
   subroutine write_spectra(spec_name, us_Mloc, up_Mloc, om_Mloc, temp_Mloc)

      !-- Input variables
      character(len=*), intent(in) :: spec_name
      complex(cp),      intent(in) :: us_Mloc(nMstart:nMstop, n_r_max)
      complex(cp),      intent(in) :: up_Mloc(nMstart:nMstop, n_r_max)
      complex(cp),      intent(in) :: om_Mloc(nMstart:nMstop, n_r_max)
      complex(cp),      intent(in) :: temp_Mloc(nMstart:nMstop, n_r_max)

      !-- Local variables
      real(cp) :: ekin(n_r_max), enst(n_r_max), buo_power(n_r_max)
      real(cp) :: ekin_m(nMstart:nMstop), enstrophy_m(nMstart:nMstop)
      real(cp) :: ekin_m_global(n_m_max), enstrophy_m_global(n_m_max)
      real(cp) :: power_m(nMstart:nMstop), power_m_global(n_m_max)
      integer :: displs(0:n_procs-1), recvcounts(0:n_procs-1)
      integer :: n_m, m, n_r, n_p, file_handle


      !-- This is not cache-friendly but hopefully it's happening only
      !-- once in a while (otherwise we need (n_r, n_m) arrays
      do n_m=nMstart,nMstop
         m = idx2m(n_m)
         do n_r=1,n_r_max
            ekin(n_r)     =cc2real(us_Mloc(n_m,n_r),m) +  &
            &              cc2real(up_Mloc(n_m,n_r),m)
            enst(n_r)     =cc2real(om_Mloc(n_m,n_r),m)
            ekin(n_r)     =ekin(n_r)*r(n_r)*height(n_r)
            enst(n_r)     =enst(n_r)*r(n_r)*height(n_r)
            buo_power(n_r)=cc22real(us_Mloc(n_m,n_r),temp_Mloc(n_m,n_r), m)
            buo_power(n_r)=Ra/pr*buo_power(n_r)*rgrav(n_r)*r(n_r)*height(n_r)
         end do
         ekin_m(n_m)     =pi*rInt_R(ekin, r, rscheme)
         enstrophy_m(n_m)=pi*rInt_R(enst, r, rscheme)
         power_m(n_m)    =two*pi*rInt_R(buo_power, r, rscheme)
      end do

      do n_p=0,n_procs-1
         recvcounts(n_p)=m_balance(n_p)%n_per_rank
      end do
      displs(0)=0
      do n_p=1,n_procs-1
         displs(n_p)=displs(n_p-1)+recvcounts(n_p-1)
      end do
      call MPI_GatherV(ekin_m, nm_per_rank, MPI_DEF_REAL,  &
           &           ekin_m_global, recvcounts, displs,  &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(enstrophy_m, nm_per_rank, MPI_DEF_REAL,  &
           &           enstrophy_m_global, recvcounts, displs,  &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(power_m, nm_per_rank, MPI_DEF_REAL,  &
           &           power_m_global, recvcounts, displs,  &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)

      if ( rank == 0 ) then
         open(newunit=file_handle, file=spec_name, position='append')
         do n_m=1,n_m_max
            m = idx2m(n_m)
            write(file_handle, '(I4, 3es16.8)') m,              &
            &              round_off(ekin_m_global(n_m)),       &
            &              round_off(enstrophy_m_global(n_m)),  &
            &              round_off(power_m_global(n_m))
         end do
         close(file_handle)
      end if

   end subroutine write_spectra
!------------------------------------------------------------------------------
end module outputs
