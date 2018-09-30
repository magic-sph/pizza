module outputs
   !
   ! This module handles the non-binary I/O of pizza:
   !        - Time series (e_kin.TAG, length_scales.TAG, heat.TAG, power.TAG)
   !        - Time-averaged radial profiles (radial_profiles.TAG)
   !

   use parallel_mod
   use precision_mod
   use spectra, only: spectra_type
   use communications, only: my_allreduce_maxloc
   use mem_alloc, only: bytes_allocated
   use namelists, only: tag, BuoFac, ra, pr, l_non_rot, l_vphi_balance, ek, &
       &                radratio, raxi, sc, tadvz_fac, kbotv, ktopv,        &
       &                ViscFac, CorFac, time_scale, r_cmb, r_icb, TdiffFac,&
       &                l_vort_balance
   use communications, only: reduce_radial_on_rank
   use truncation, only: n_r_max, m2idx, n_m_max, idx2m, minc
   use radial_functions, only: r, rscheme, rgrav, dtcond, height, tcond, &
       &                       beta, ekpump, or1, oheight, or2
   use blocking, only: nMstart, nMstop, nRstart, nRstop, l_rank_has_m0, &
       &               nm_per_rank, m_balance
   use integration, only: rInt_R, simps
   use useful, only: round_off, cc2real, cc22real, getMSD2, abortRun
   use constants, only: pi, two, four, surf, vol_otc, one
   use checkpoints, only: write_checkpoint_mloc
   use output_frames, only: write_snapshot_mloc
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use char_manip, only: capitalize
   use vp_balance, only: vp_bal_type
   use vort_balance, only: vort_bal_type

   implicit none

   private

   integer :: frame_counter, n_calls
   real(cp) :: timeLast_rad, timeAvg_rad, timeAvg_vortbal
   real(cp) :: timeAvg_spec
   integer, public :: n_log_file
   integer :: n_rey_file_2D, n_heat_file, n_kin_file_2D, n_power_file_2D
   integer :: n_lscale_file, n_sig_file
   integer :: n_rey_file_3D, n_power_file_3D, n_kin_file_3D
   character(len=144), public :: log_file

   real(cp), allocatable :: uphiR_mean(:), uphiR_SD(:)
   real(cp), allocatable :: us2R_mean(:), us2R_SD(:), up2R_mean(:), up2R_SD(:)
   real(cp), allocatable :: enstrophyR_mean(:), enstrophyR_SD(:)
   real(cp), allocatable :: tempR_mean(:), tempR_SD(:)
   real(cp), allocatable :: fluxR_mean(:), fluxR_SD(:)

   type(vp_bal_type), public :: vp_bal
   type(spectra_type), public :: spec

   type(vort_bal_type), public :: vort_bal
   public :: initialize_outputs, finalize_outputs, get_time_series, &
   &         write_outputs, read_signal_file

contains

   subroutine initialize_outputs

      logical :: log_exists
      character(len=144) :: file_name

      if ( rank == 0 ) then
         log_file = 'log.'//tag
         inquire(file=log_file, exist=log_exists)
         if ( log_exists ) then
            call abortRun('! log file already exists, please change tag')
         end if
         open(newunit=n_log_file, file=log_file, status='new')
         file_name = 'e_kin.'//tag
         open(newunit=n_kin_file_2D, file=file_name, status='new')
         file_name = 'power.'//tag
         open(newunit=n_power_file_2D, file=file_name, status='new')
         if ( index(time_scale, 'ROT') /= 0 ) then
            file_name = 'rossby.'//tag
         else
            file_name = 'reynolds.'//tag
         end if
         open(newunit=n_rey_file_2D,file=file_name, status='new')
         file_name = 'length_scales.'//tag
         open(newunit=n_lscale_file, file=file_name, status='new')
         file_name = 'heat.'//tag
         open(newunit=n_heat_file, file=file_name, status='new')

         if ( .not. l_non_rot ) then
            file_name = 'e_kin_3D.' // tag
            open(newunit=n_kin_file_3D, file=file_name, status='new')
            file_name = 'power_3D.' // tag
            open(newunit=n_power_file_3D, file=file_name, status='new')
            if ( index(time_scale, 'ROT') /= 0 ) then
               file_name = 'rossby_3D.'//tag
            else
               file_name = 'reynolds_3D.'//tag
            end if
            open(newunit=n_rey_file_3D,file=file_name, status='new')
         end if

         file_name = 'signal.'//tag
         open(newunit=n_sig_file, file=file_name, status='unknown')
         write(n_sig_file,'(A3)') 'NOT'
         close(n_sig_file)
      end if 

      timeAvg_rad     = 0.0_cp
      timeAvg_spec    = 0.0_cp
      timeAvg_vortbal = 0.0_cp

      if ( l_rank_has_m0 ) then

         allocate( uphiR_mean(n_r_max), uphiR_SD(n_r_max) )
         allocate( tempR_mean(n_r_max), tempR_SD(n_r_max) )
         allocate( fluxR_mean(n_r_max), fluxR_SD(n_r_max) )
         allocate( us2R_mean(n_r_max), us2R_SD(n_r_max) )
         allocate( up2R_mean(n_r_max), up2R_SD(n_r_max) )
         allocate( enstrophyR_mean(n_r_max), enstrophyR_SD(n_r_max) )
         bytes_allocated=bytes_allocated+10*n_r_max*SIZEOF_DEF_REAL

         uphiR_mean(:)      = 0.0_cp
         uphiR_SD(:)        = 0.0_cp
         tempR_mean(:)      = 0.0_cp
         tempR_SD(:)        = 0.0_cp
         fluxR_mean(:)      = 0.0_cp
         fluxR_SD(:)        = 0.0_cp
         us2R_mean(:)       = 0.0_cp
         us2R_SD(:)         = 0.0_cp
         up2R_mean(:)       = 0.0_cp
         up2R_SD(:)         = 0.0_cp
         enstrophyR_mean(:) = 0.0_cp
         enstrophyR_SD(:)   = 0.0_cp
         n_calls            = 0
         timeLast_rad       = 0.0_cp
      end if

      frame_counter = 1 ! For file suffix

      call spec%initialize()

      if ( l_vphi_balance ) call vp_bal%initialize()
      if ( l_vort_balance ) call vort_bal%initialize()

   end subroutine initialize_outputs
!------------------------------------------------------------------------------
   subroutine finalize_outputs

      if ( rank == 0 ) then
         deallocate( uphiR_mean, uphiR_SD, tempR_mean, tempR_SD )
         deallocate( fluxR_mean, fluxR_SD, us2R_mean, us2R_SD )
         deallocate( up2R_mean, up2R_SD, enstrophyR_mean, enstrophyR_SD )
      end if

      call spec%finalize()

      if ( l_vort_balance ) call vort_bal%finalize()

      if ( l_vphi_balance ) call vp_bal%finalize()

      if ( rank == 0 ) then
         if ( .not. l_non_rot ) then
            close(n_rey_file_3D)
            close(n_power_file_3D)
            close(n_kin_file_3D)
         end if
         close(n_heat_file)
         close(n_lscale_file)
         close(n_rey_file_2D)
         close(n_power_file_2D)
         close(n_kin_file_2D)
         close(n_log_file)
      end if

   end subroutine finalize_outputs
!------------------------------------------------------------------------------
   subroutine read_signal_file(signals)

      !-- Outputs signals
      integer, intent(inout) :: signals(4)

      !-- Local variables:
      character(len=255) :: message
      character(len=76) :: SIG

      signals(:) = 0

      if ( rank == 0 ) then
         !----- Signalling via file signal:
         message='signal.'//tag
         open(newunit=n_sig_file, file=trim(message), status='old')
         read(n_sig_file,*) SIG
         close(n_sig_file)
         if ( len(trim(SIG)) > 0 ) then ! Non blank string ?
            call capitalize(SIG)

            if ( index(SIG,'END')/=0 ) signals(1)=1  !n_stop_signal=1

            if ( index(SIG,'FRA')/=0 ) then
               signals(2)=1
               open(newunit=n_sig_file, file=trim(message), status='unknown')
               write(n_sig_file,'(A3)') 'NOT'
               close(n_sig_file)
            end if
            if ( index(SIG,'RST')/=0 ) then
               signals(3)=1
               open(newunit=n_sig_file, file=trim(message), status='unknown')
               write(n_sig_file,'(A3)') 'NOT'
               close(n_sig_file)
            end if
            if ( index(SIG,'SPE')/=0 ) then
               signals(4)=1
               open(newunit=n_sig_file, file=trim(message), status='unknown')
               write(n_sig_file,'(A3)') 'NOT'
               close(n_sig_file)
            end if
         end if
      end if

      call MPI_Bcast(signals,4,MPI_Integer,0,MPI_COMM_WORLD,ierr)

   end subroutine read_signal_file
!------------------------------------------------------------------------------
   subroutine write_outputs(time, tscheme, n_time_step, l_log, l_rst, l_frame, &
              &             l_vphi_bal_write, l_stop_time,  us_Mloc, up_Mloc,  &
              &             om_Mloc, temp_Mloc, dtemp_Mloc, dpsidt, dTdt)

      !-- Input variables
      real(cp),            intent(in) :: time
      class(type_tscheme), intent(in) :: tscheme
      integer,             intent(in) :: n_time_step
      logical,             intent(in) :: l_log
      logical,             intent(in) :: l_rst
      logical,             intent(in) :: l_frame
      logical,             intent(in) :: l_vphi_bal_write
      logical,             intent(in) :: l_stop_time
      complex(cp),         intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(in) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(in) :: temp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(in) :: dtemp_Mloc(nMstart:nMstop,n_r_max)
      type(type_tarray),   intent(in) :: dpsidt
      type(type_tarray),   intent(in) :: dTdt

      !-- Local variable
      character(len=144) :: frame_name
      real(cp) :: us2_r(n_r_max), up2_r(n_r_max), enstrophy_r(n_r_max)
      real(cp) :: us2_m_Mloc(nMstart:nMstop), up2_m_Mloc(nMstart:nMstop)
      real(cp) :: enst_m_Mloc(nMstart:nMstop)
      real(cp) :: flux_r(n_r_max)

      timeAvg_rad  = timeAvg_rad  + tscheme%dt(1)
      timeAvg_spec = timeAvg_spec + tscheme%dt(1)

      !-- Write checkpoints
      if ( l_rst ) then
         call write_checkpoint_mloc(time, tscheme, n_time_step, n_log_file,   &
              &                     l_stop_time, temp_Mloc, us_Mloc, up_Mloc, &
              &                     dTdt, dpsidt)
      end if

      !-- Calculate spectra
      if ( l_log .or. spec%l_calc ) then
         call spec%calculate_spectra(timeAvg_spec, l_stop_time, us_Mloc, &
              &                      up_Mloc, om_Mloc, us2_m_Mloc,       &
              &                      up2_m_Mloc, enst_m_Mloc)
      end if

      !-- Write spectra
      if ( spec%l_calc ) then
         call spec%write_spectra(us2_m_Mloc, up2_m_Mloc, enst_m_Mloc)
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
         call get_time_series(time, us_Mloc, up_Mloc, om_Mloc, temp_Mloc,      &
              &               dtemp_Mloc, us2_m_Mloc, up2_m_Mloc, enst_m_Mloc, &
              &               us2_r, up2_r, enstrophy_r, flux_r)

         call get_radial_averages(timeAvg_rad, l_stop_time, up_Mloc, temp_Mloc, &
              &                   us2_r, up2_r, enstrophy_r, flux_r)
      end if

      if ( l_vphi_bal_write ) then
         call vp_bal%write_outputs(time, up_Mloc)
      end if

      if ( vort_bal%l_calc .or. (l_stop_time .and. l_vort_balance) ) then
         timeAvg_vortbal = timeAvg_vortbal + tscheme%dt(1)
         call vort_bal%calc_avg(timeAvg_vortbal,l_stop_time)
      end if

   end subroutine write_outputs
!------------------------------------------------------------------------------
   subroutine get_radial_averages(timeAvg_rad, l_stop_time, up_Mloc, temp_Mloc, &
              &                   us2_r, up2_r, enstrophy_r, flux_r)

      !-- Input variables
      real(cp),    intent(in) :: timeAvg_rad
      logical,     intent(in) :: l_stop_time
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(in) :: temp_Mloc(nMstart:nMstop, n_r_max)
      real(cp),    intent(in) :: us2_r(n_r_max)
      real(cp),    intent(in) :: up2_r(n_r_max)
      real(cp),    intent(in) :: enstrophy_r(n_r_max)
      real(cp),    intent(in) :: flux_r(n_r_max)

      !-- Local variables
      real(cp) :: dtAvg
      integer :: n_r, idx, file_handle


      if ( l_rank_has_m0 ) then

         n_calls = n_calls+1
         dtAvg = timeAvg_rad-timeLast_rad

         idx = m2idx(0)
         do n_r=1,n_r_max
            call getMSD2(uphiR_mean(n_r), uphiR_SD(n_r), real(up_Mloc(idx,n_r)),&
                 &       n_calls, dtAvg, timeAvg_rad)
            call getMSD2(tempR_mean(n_r), tempR_SD(n_r), real(temp_Mloc(idx,n_r)),&
                 &       n_calls, dtAvg, timeAvg_rad)
            call getMSD2(fluxR_mean(n_r), fluxR_SD(n_r), flux_r(n_r), n_calls, &
                 &       dtAvg, timeAvg_rad)
            call getMSD2(us2R_mean(n_r), us2R_SD(n_r), two*pi*us2_r(n_r), &
                 &       n_calls, dtAvg, timeAvg_rad)
            call getMSD2(up2R_mean(n_r), up2R_SD(n_r), two*pi*up2_r(n_r), &
                 &       n_calls, dtAvg, timeAvg_rad)
            call getMSD2(enstrophyR_mean(n_r), enstrophyR_SD(n_r),    &
                 &       enstrophy_r(n_r), n_calls, dtAvg, timeAvg_rad)
         end do
         timeLast_rad = timeAvg_rad

         if ( l_stop_time ) then
            open(newunit=file_handle, file='radial_profiles.'//tag)
            do n_r=1,n_r_max
               uphiR_SD(n_r)     =sqrt(uphiR_SD(n_r)/timeAvg_rad)
               tempR_SD(n_r)     =sqrt(tempR_SD(n_r)/timeAvg_rad)
               fluxR_SD(n_r)     =sqrt(fluxR_SD(n_r)/timeAvg_rad)
               us2R_SD(n_r)      =sqrt(us2R_SD(n_r)/timeAvg_rad)
               up2R_SD(n_r)      =sqrt(up2R_SD(n_r)/timeAvg_rad)
               enstrophyR_SD(n_r)=sqrt(enstrophyR_SD(n_r)/timeAvg_rad)
               write(file_handle, '(es20.12, 12es16.8)') r(n_r),           &
               &     round_off(us2R_mean(n_r)), round_off(us2R_SD(n_r)),   &
               &     round_off(up2R_mean(n_r)), round_off(up2R_SD(n_r)),   &
               &     round_off(enstrophyR_mean(n_r)),                      &
               &     round_off(enstrophyR_SD(n_r)),                        &
               &     round_off(uphiR_mean(n_r)), round_off(uphiR_SD(n_r)), &
               &     round_off(tempR_mean(n_r)+tcond(n_r)),                &
               &     round_off(tempR_SD(n_r)), round_off(fluxR_mean(n_r)), &
               &     round_off(fluxR_SD(n_r))
            end do
            close(file_handle)
         end if

      end if

   end subroutine get_radial_averages
!------------------------------------------------------------------------------
   subroutine get_time_series(time, us_Mloc, up_Mloc, om_Mloc, temp_Mloc, &
              &               dtemp_Mloc, us2_m, up2_m, enstrophy_m,      &
              &               us2_r, up2_r, enstrophy, flux_r)

      !-- Input variables
      real(cp),    intent(in) :: time
      complex(cp), intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: temp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: dtemp_Mloc(nMstart:nMstop,n_r_max)
      real(cp),    intent(in) :: us2_m(nMstart:nMstop)
      real(cp),    intent(in) :: up2_m(nMstart:nMstop)
      real(cp),    intent(in) :: enstrophy_m(nMstart:nMstop)

      !-- Output variables
      real(cp), intent(out) :: us2_r(n_r_max)
      real(cp), intent(out) :: up2_r(n_r_max)
      real(cp), intent(out) :: enstrophy(n_r_max)
      real(cp), intent(out) :: flux_r(n_r_max)

      !-- Local variables
      integer :: n_r, n_m, m, n_m0
      real(cp) :: dlus_peak, dlekin_peak, dlvort_peak
      real(cp) :: dl_diss, dlekin_int, pow_3D
      real(cp) :: tTop, tBot, visc_2D, pow_2D, pum, visc_3D
      real(cp) :: us2_2D, up2_2D, up2_axi_2D, us2_3D, up2_3D, up2_axi_3D, uz2_3D
      real(cp) :: up2_axi_r(n_r_max), nu_vol_r(n_r_max), nu_cond_r(n_r_max)
      real(cp) :: buo_power(n_r_max), pump(n_r_max), tmp(n_r_max)
      real(cp) :: theta(n_r_max)
      real(cp) :: NuTop, NuBot, beta_t, Nu_vol, Nu_int, fac
      real(cp) :: rey_2D, rey_fluct_2D, rey_zon_2D
      real(cp) :: rey_3D, rey_fluct_3D, rey_zon_3D

      do n_r=1,n_r_max
         us2_r(n_r)    =0.0_cp
         up2_r(n_r)    =0.0_cp
         up2_axi_r(n_r)=0.0_cp
         enstrophy(n_r)=0.0_cp
         buo_power(n_r)=0.0_cp
         pump(n_r)     =0.0_cp
         nu_vol_r(n_r) =0.0_cp
         flux_r(n_r)   =0.0_cp
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            us2_r(n_r)    =us2_r(n_r)+cc2real(us_Mloc(n_m,n_r),m)
            up2_r(n_r)    =up2_r(n_r)+cc2real(up_Mloc(n_m,n_r),m)
            enstrophy(n_r)=enstrophy(n_r)+cc2real(om_Mloc(n_m,n_r),m)
            buo_power(n_r)=buo_power(n_r)+cc22real(us_Mloc(n_m,n_r), &
            &              temp_Mloc(n_m,n_r), m)

            nu_vol_r(n_r) =nu_vol_r(n_r)+cc2real(dtemp_Mloc(n_m,n_r),m)&
            &              +real(m,cp)*real(m,cp)*or2(n_r)*            &
            &                             cc2real(temp_Mloc(n_m,n_r),m)

            flux_r(n_r)   =flux_r(n_r)+cc22real(us_Mloc(n_m,n_r), &
            &              temp_Mloc(n_m,n_r), m)

            if ( m == 0 ) then
               up2_axi_r(n_r)=up2_axi_r(n_r)+cc2real(up_Mloc(n_m,n_r),m)
               pump(n_r)     =pump(n_r)+cc2real(up_Mloc(n_m,n_r),m)
               flux_r(n_r)   =flux_r(n_r)-TdiffFac*real(dtemp_Mloc(n_m,n_r))-&
               &              TdiffFac*dtcond(n_r)
            end if
         end do
         nu_vol_r(n_r) =nu_vol_r(n_r)*r(n_r)*height(n_r)
      end do

      !-- MPI reductions to get the s-profiles on rank==0
      call reduce_radial_on_rank(nu_vol_r, 0)
      call reduce_radial_on_rank(flux_r, 0)
      call reduce_radial_on_rank(us2_r, 0)
      call reduce_radial_on_rank(up2_r, 0)
      call reduce_radial_on_rank(up2_axi_r, 0)
      call reduce_radial_on_rank(enstrophy, 0)
      call reduce_radial_on_rank(buo_power, 0)
      call reduce_radial_on_rank(pump, 0)

      !-------
      !-- Lengthscales
      !-------
      call get_lengthscales(us2_m, up2_m, enstrophy_m, dlus_peak, dlekin_peak,  &
           &                dlvort_peak, dlekin_int)

      if ( rank == 0 ) then

         !-----
         !-- Kinetic energy (2D and 3D)
         !-----
         tmp(:) = us2_r(:)*r(:)
         us2_2D = rInt_R(tmp, r, rscheme)
         us2_2D = round_off(pi*us2_2D)

         tmp(:) = up2_r(:)*r(:)
         up2_2D = rInt_R(tmp, r, rscheme)
         up2_2D = round_off(pi*up2_2D)

         tmp(:) = up2_axi_r(:)*r(:)
         up2_axi_2D = rInt_R(tmp, r, rscheme)
         up2_axi_2D = round_off(pi*up2_axi_2D)

         if ( .not. l_non_rot ) then
            tmp(:) = us2_r(:)*r(:)*height(:)
            us2_3D = rInt_R(tmp, r, rscheme)
            us2_3D = round_off(pi*us2_3D)

            tmp(:) = up2_r(:)*r(:)*height(:)
            up2_3D = rInt_R(tmp, r, rscheme)
            up2_3D = round_off(pi*up2_3D)

            tmp(:) = 4.0_cp/3.0_cp*us2_r(:)*r(:)*r(:)*r(:)*oheight(:)
            uz2_3D = rInt_R(tmp, r, rscheme)
            uz2_3D = round_off(pi*uz2_3D)

            tmp(:) = up2_axi_r(:)*r(:)*height(:)
            up2_axi_3D = rInt_R(tmp, r, rscheme)
            up2_axi_3D = round_off(pi*up2_axi_3D)
         end if

         !-----
         !-- Reynolds numbers
         !-----
         rey_2D       = sqrt(two*(us2_2D+up2_2D)/surf)
         rey_zon_2D   = sqrt(two*up2_axi_2D/surf)
         rey_fluct_2D = sqrt(two*(us2_2D+up2_2D-up2_axi_2D)/surf)
         
         if ( .not. l_non_rot ) then
            rey_3D       = sqrt(two*(us2_3D+up2_3D+uz2_3D)/vol_otc)
            rey_zon_3D   = sqrt(two*up2_axi_3D/vol_otc)
            rey_fluct_3D = sqrt(two*(us2_3D+up2_3D+uz2_3D-up2_axi_3D)/vol_otc)
         end if

         !-----
         !-- Power budget
         !-----
         tmp(:) = enstrophy(:)*r(:)
         visc_2D = rInt_R(tmp, r, rscheme)
         visc_2D = round_off(two*pi*visc_2D)
         !-- In case of stress-free we need to include the surface contributions
         if ( kbotv == 1 ) then
            visc_2D = visc_2D-four*or1(n_r_max)*up2_r(n_r_max)
         end if
         if ( ktopv == 1 ) then
            visc_2D = visc_2D-four*or1(1)*up2_r(1)
         end if

         tmp(:)=BuoFac*buo_power(:)*rgrav(:)*r(:)
         pow_2D = rInt_R(tmp, r, rscheme)
         pow_2D = round_off(two*pi*pow_2D)

         if ( .not. l_non_rot ) then
            tmp(:) = enstrophy(:)*r(:)*height(:)
            visc_3D = rInt_R(tmp, r, rscheme)
            visc_3D = round_off(two*pi*visc_3D)
            !-- In case of stress-free we need to include the surface contributions
            if ( kbotv == 1 ) then
               visc_3D = visc_3D-four*or1(n_r_max)*up2_r(n_r_max)
            end if
            if ( ktopv == 1 ) then
               visc_3D = visc_3D-four*or1(1)*up2_r(1)
            end if

            tmp(:)=BuoFac*buo_power(:)*rgrav(:)*r(:)*height(:)
            pow_3D = rInt_R(tmp, r, rscheme)
            pow_3D = round_off(two*pi*pow_3D)

            !-- \int\int (  g*T*(us*s/r+uz*z/r)*s dz ds )
            !-- =\int\int (  g*T*(us*s/r+beta*us*z**2/r)*s dz ds )
            !-- with r=sqrt(s**2+z**2)
            ! tmp(:)=-BuoFac*buo_power(:)*rgrav(:)*r(:)*beta(:)*( &
            ! &      log(or1(:)*(0.5_cp*height(:)+r_cmb))*(       &
            ! &      two*r_cmb*r_cmb-r(:)*r(:))-                  &
            ! &      0.5_cp*height(:)*r_cmb )
            ! pow_3D_b = rInt_R(tmp, r, rscheme)
            ! pow_3D_b = round_off(two*pi*pow_3D_b)

            tmp(:)=CorFac*pump(:)*ekpump(:)*height(:)*r(:)
            pum  = rInt_R(tmp, r, rscheme)
            pum  = round_off(two*pi*pum)
         end if


         write(n_kin_file_2D, '(1P, es20.12, 3es16.8)') time, us2_2D, up2_2D, &
         &                                              up2_axi_2D
         write(n_power_file_2D, '(1P, es20.12, 2es16.8)') time, pow_2D, visc_2D

         write(n_rey_file_2D, '(1P, es20.12, 3es16.8)') time, rey_2D, rey_zon_2D, &
         &                                              rey_fluct_2D

         if ( .not. l_non_rot ) then
            write(n_kin_file_3D, '(1P, es20.12, 4es16.8)') time, us2_3D, up2_3D, &
            &                                              uz2_3D, up2_axi_3D
            write(n_power_file_3D, '(1P, es20.12, 3es16.8)') time, pow_3D, &
            &                                                visc_3D, pum

            write(n_rey_file_3D, '(1P, es20.12, 3es16.8)') time, rey_3D, &
            &                                              rey_zon_3D,   &
            &                                              rey_fluct_3D
         end if

         if ( l_non_rot ) then
            !-- Dissipation lengthscale = \sqrt(2*Ekin/\omega^2)
            if ( abs(visc_2D) > 10.0_cp*epsilon(one) ) then
               dl_diss = sqrt(two*(us2_2D+up2_2D)*viscfac/visc_2D)
            else
               dl_diss = 0.0_cp
            end if
         else
            !-- Dissipation lengthscale = \sqrt(2*Ekin/\omega^2)
            if ( abs(visc_3D) > 10.0_cp*epsilon(one) ) then
               dl_diss = sqrt(two*(us2_3D+up2_3D)*viscfac/visc_3D)
            else
               dl_diss = 0.0_cp
            end if

         end if
         write(n_lscale_file, '(1P, es20.12, 5es16.8)') time, dlus_peak,          &
         &                                              dlekin_peak, dlvort_peak, &
         &                                              dlekin_int, dl_diss

         ! At this stage multiply us2_r, up2_r and enstrophy by r(:) and height(:)
         us2_r(:)    =us2_r(:)*r(:)*height(:)
         up2_r(:)    =up2_r(:)*r(:)*height(:)
         enstrophy(:)=enstrophy(:)*r(:)*height(:)

      end if

      if ( l_rank_has_m0 ) then
         !------
         !-- Heat transfer
         !------

         !-- Top and bottom temperatures
         n_m0 = m2idx(0)
         tTop = real(temp_Mloc(n_m0,1))+tcond(1)
         tBot = real(temp_Mloc(n_m0,n_r_max))+tcond(n_r_max)
         tTop = round_off(tTop)
         tBot = round_off(tBot)

         !-- Classical top and bottom Nusselt number
         NuTop = one+real(dtemp_Mloc(n_m0,1))/dtcond(1)
         NuBot = one+real(dtemp_Mloc(n_m0,n_r_max))/dtcond(n_r_max)
         !&       (dtcond(n_r_max)-tadvz_fac*beta(n_r_max)*tcond(n_r_max))
         NuTop = round_off(NuTop)
         NuBot = round_off(NuBot)

         !-- Volume-based Nusselt number
         Nu_vol = rInt_R(nu_vol_r, r, rscheme)
         nu_cond_r(:)=dtcond(:)*dtcond(:)*height(:)*r(:)
         Nu_vol = one+Nu_vol/rInt_R(nu_cond_r, r, rscheme)
         Nu_vol = round_off(Nu_vol)

         !-- Spherical shell surface-based Nusselt number
         if ( l_non_rot ) then
            tmp(:) = -TdiffFac*dtcond(:)*r(:)
            fac = rInt_R(tmp, r, rscheme)
            flux_r(:) = flux_r(:)*r(:)/fac
            Nu_int = rInt_R(flux_r, r, rscheme)
         else
            !-- The first sin(theta) comes from the projection of Flux_s to Flux_r
            !-- The second one from the surface integral on a spherical shell
            tmp(:) = -TdiffFac*dtcond(:) * (r(:)/r_cmb)**2 ! s/ro is sin(theta)
            theta(:) = asin(r(:)/r_cmb)
            !-- Integration over colatitudes (only simpson can work here)
            fac = simps(tmp, theta)

            !-- Again multiplication by sin(theta)**2 for the spherical surface
            tmp(:) = flux_r(:) * (r(:)/r_cmb)**2
            !-- Integration over colatitudes (only simpson can work here)
            Nu_int = simps(tmp, theta)
            Nu_int = Nu_int/fac

            !-- Finally transform the flux_r into a Nusselt(s) profile
            flux_r(:)=flux_r(:) * (r(:)/r_cmb)/fac
         end if
         Nu_int = round_off(Nu_int)

         !-- Mid-shell temperature gradient
         beta_t = dtcond(int(n_r_max/2))+real(dtemp_Mloc(n_m0,int(n_r_max/2)))
         beta_t = round_off(beta_t)

         write(n_heat_file, '(1P, ES20.12, 7ES16.8)') time, NuTop, NuBot,   &
         &                                            Nu_vol, Nu_int, tTop, &
         &                                            tBot, beta_t
      end if

   end subroutine get_time_series
!------------------------------------------------------------------------------
   subroutine get_lengthscales(us2_m, up2_m, enst_m, dlus_peak, dlekin_peak, &
              &                dlvort_peak, dlekin_int)

      !-- Input variables
      real(cp), intent(in) :: us2_m(nMstart:nMstop)
      real(cp), intent(in) :: up2_m(nMstart:nMstop)
      real(cp), intent(in) :: enst_m(nMstart:nMstop)

      !-- Output variables
      real(cp), intent(out) :: dlus_peak
      real(cp), intent(out) :: dlekin_peak
      real(cp), intent(out) :: dlvort_peak
      real(cp), intent(out) :: dlekin_int

      !-- Local variables:
      real(cp) :: E, Em
      integer :: lus_peak, lvort_peak, lekin_peak, m, n_m

      !-- Estimate the peak of the spectra (exclude m=0 for practical reason)
      lus_peak=my_allreduce_maxloc(us2_m(max(2,nMstart):nMstop))*minc
      if ( lus_peak > 0 ) then
         dlus_peak = pi/lus_peak
      else
         dlus_peak = 0.0_cp
      end if

      lekin_peak=my_allreduce_maxloc(us2_m(max(2,nMstart):nMstop)+ &
                 &                   up2_m(max(2,nMstart):nMstop))*minc
      if ( lekin_peak > 0 ) then
         dlekin_peak = pi/lekin_peak
      else
         dlekin_peak = 0.0_cp
      end if

      lvort_peak=my_allreduce_maxloc(enst_m(max(2,nMstart):nMstop))*minc
      if ( lvort_peak > 0 ) then
         dlvort_peak = pi/lvort_peak
      else
         dlvort_peak = 0.0_cp
      end if

      !-- Integral lengthscale
      E = 0.0_cp
      Em = 0.0_cp
      do n_m=nMstart,nMstop
         m=idx2m(n_m)
         E  = E  + (us2_m(n_m)+up2_m(n_m))
         Em = Em + real(m,cp)*(us2_m(n_m)+up2_m(n_m))
      end do

      !-- MPI reduction
      call MPI_AllReduce(MPI_IN_PLACE, E, 1, MPI_DEF_REAL, MPI_SUM, &
           &             MPI_COMM_WORLD, ierr)
      call MPI_AllReduce(MPI_IN_PLACE, Em, 1, MPI_DEF_REAL, MPI_SUM, &
           &             MPI_COMM_WORLD, ierr)

      if ( abs(E) > 10.0_cp*epsilon(one) ) then
         dlekin_int = pi*E/Em
      else
         dlekin_int = 0.0_cp
      end if

   end subroutine get_lengthscales
!------------------------------------------------------------------------------
end module outputs
