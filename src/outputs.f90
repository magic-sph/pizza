module outputs
   !
   ! This module handles the non-binary I/O of pizza:
   !        - Time series (e_kin.TAG, length_scales.TAG, heat.TAG, power.TAG)
   !        - Time-averaged radial profiles (radial_profiles.TAG)
   !

   use parallel_mod
   use precision_mod
   use spectra, only: spectra_type
   use communications, only: my_allreduce_maxloc, reduce_radial_on_rank
   use namelists, only: tag, BuoFac, ra, pr, l_non_rot, l_vphi_balance, ek, &
       &                radratio, raxi, sc, tadvz_fac, kbotv, ktopv,        &
       &                ViscFac, CorFac, time_scale, r_cmb, r_icb, TdiffFac,&
       &                l_vort_balance, l_corr, l_heat, l_chem, ChemFac
   use communications, only: reduce_radial_on_rank
   use truncation, only: n_r_max, m2idx, n_m_max, idx2m, minc
   use radial_functions, only: r, rscheme, rgrav, dtcond, height, tcond, &
       &                       beta, ekpump, or1, oheight, or2, xicond,  &
       &                       dxicond
   use blocking, only: nMstart, nMstop, nRstart, nRstop, l_rank_has_m0, &
       &               nm_per_rank, m_balance
   use integration, only: rInt_R, simps
   use useful, only: round_off, cc2real, cc22real, getMSD2, abortRun
   use constants, only: pi, two, four, surf, vol_otc, one, tiny_number
   use checkpoints, only: write_checkpoint_mloc
   use output_frames, only: write_snapshot_mloc
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use char_manip, only: capitalize
   use vp_balance, only: vp_bal_type
   use vort_balance, only: vort_bal_type
   use mean_sd, only: mean_sd_type

   implicit none

   private

   integer :: frame_counter, n_calls
   real(cp) :: timeLast_rad, timeAvg_rad, timeAvg_vortbal
   real(cp) :: timeAvg_spec
   integer, public :: n_log_file
   integer :: n_rey_file_2D, n_heat_file, n_kin_file_2D, n_power_file_2D
   integer :: n_lscale_file, n_sig_file, n_corr_file, n_chem_file
   integer :: n_rey_file_3D, n_power_file_3D, n_kin_file_3D
   character(len=144), public :: log_file

   type(mean_sd_type) :: uphiR, us2R, up2R, enstrophyR, tempR, fluxR, xiR

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

         if ( l_heat ) then
            file_name = 'heat.'//tag
            open(newunit=n_heat_file, file=file_name, status='new')
         end if

         if ( l_chem ) then
            file_name = 'composition.'//tag
            open(newunit=n_chem_file, file=file_name, status='new')
         end if

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

         if ( l_corr ) then
            file_name = 'corr.' // tag
            open(newunit=n_corr_file, file=file_name, status='new')
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
         call uphiR%initialize(1,n_r_max)
         call tempR%initialize(1,n_r_max)
         call xiR%initialize(1,n_r_max)
         call fluxR%initialize(1,n_r_max)
         call us2R%initialize(1,n_r_max)
         call up2R%initialize(1,n_r_max)
         call enstrophyR%initialize(1,n_r_max)
         n_calls     =0
         timeLast_rad=0.0_cp
      end if

      frame_counter = 1 ! For file suffix

      call spec%initialize()

      if ( l_vphi_balance ) call vp_bal%initialize()
      if ( l_vort_balance ) call vort_bal%initialize()

   end subroutine initialize_outputs
!------------------------------------------------------------------------------
   subroutine finalize_outputs

      if ( rank == 0 ) then
         call enstrophyR%finalize()
         call up2R%finalize()
         call us2R%finalize()
         call fluxR%finalize()
         call tempR%finalize()
         call xiR%finalize()
         call uphiR%finalize()
      end if

      call spec%finalize()

      if ( l_vort_balance ) call vort_bal%finalize()

      if ( l_vphi_balance ) call vp_bal%finalize()

      if ( rank == 0 ) then
         if ( l_corr ) close(n_corr_file)
         if ( .not. l_non_rot ) then
            close(n_rey_file_3D)
            close(n_power_file_3D)
            close(n_kin_file_3D)
         end if
         if ( l_chem ) close(n_chem_file)
         if ( l_heat ) close(n_heat_file)
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
              &             om_Mloc, temp_Mloc, dtemp_Mloc, xi_Mloc, dxi_Mloc, &
              &             dpsidt, dTdt, dxidt)

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
      complex(cp),         intent(in) :: xi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(in) :: dtemp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(in) :: dxi_Mloc(nMstart:nMstop,n_r_max)
      type(type_tarray),   intent(in) :: dpsidt
      type(type_tarray),   intent(in) :: dTdt
      type(type_tarray),   intent(in) :: dxidt

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
              &                     l_stop_time, temp_Mloc, xi_Mloc, us_Mloc, &
              &                     up_Mloc, dTdt, dxidt, dpsidt)
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
         if ( l_heat ) then
            write(frame_name, '(A,I0,A,A)') 'frame_temp_',frame_counter,'.',tag
            call write_snapshot_mloc(frame_name, time, temp_Mloc)
         end if
         if ( l_chem ) then
            write(frame_name, '(A,I0,A,A)') 'frame_xi_',frame_counter,'.',tag
            call write_snapshot_mloc(frame_name, time, xi_Mloc)
         end if
         write(frame_name, '(A,I0,A,A)') 'frame_us_',frame_counter,'.',tag
         call write_snapshot_mloc(frame_name, time, us_Mloc)
         write(frame_name, '(A,I0,A,A)') 'frame_up_',frame_counter,'.',tag
         call write_snapshot_mloc(frame_name, time, up_Mloc)
         write(frame_name, '(A,I0,A,A)') 'frame_om_',frame_counter,'.',tag
         call write_snapshot_mloc(frame_name, time, om_Mloc)
         frame_counter = frame_counter+1
      end if

      if ( l_log ) then
         call get_time_series(time, us_Mloc, up_Mloc, om_Mloc, temp_Mloc,  &
              &               dtemp_Mloc, xi_Mloc, dxi_Mloc, us2_m_Mloc,   &
              &               up2_m_Mloc, enst_m_Mloc, us2_r, up2_r,       &
              &               enstrophy_r, flux_r)

         call get_radial_averages(timeAvg_rad, l_stop_time, up_Mloc, temp_Mloc, &
              &                   xi_Mloc, us2_r, up2_r, enstrophy_r, flux_r)
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
              &                   xi_Mloc, us2_r, up2_r, enstrophy_r, flux_r)

      !-- Input variables
      real(cp),    intent(in) :: timeAvg_rad
      logical,     intent(in) :: l_stop_time
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(in) :: temp_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(in) :: xi_Mloc(nMstart:nMstop, n_r_max)
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
            call getMSD2(uphiR%mean(n_r), uphiR%SD(n_r), real(up_Mloc(idx,n_r)),&
                 &       n_calls, dtAvg, timeAvg_rad)
            if ( l_heat ) then
               call getMSD2(tempR%mean(n_r), tempR%SD(n_r), &
                    &       real(temp_Mloc(idx,n_r)), n_calls, dtAvg, timeAvg_rad)
            end if
            if ( l_chem ) then
               call getMSD2(xiR%mean(n_r), xiR%SD(n_r), &
                    &       real(xi_Mloc(idx,n_r)), n_calls, dtAvg, timeAvg_rad)
            end if
            call getMSD2(fluxR%mean(n_r), fluxR%SD(n_r), flux_r(n_r), n_calls, &
                 &       dtAvg, timeAvg_rad)
            call getMSD2(us2R%mean(n_r), us2R%SD(n_r), pi*us2_r(n_r), &
                 &       n_calls, dtAvg, timeAvg_rad)
            call getMSD2(up2R%mean(n_r), up2R%SD(n_r), pi*up2_r(n_r), &
                 &       n_calls, dtAvg, timeAvg_rad)
            call getMSD2(enstrophyR%mean(n_r), enstrophyR%SD(n_r),    &
                 &       enstrophy_r(n_r), n_calls, dtAvg, timeAvg_rad)
         end do
         timeLast_rad = timeAvg_rad

         if ( l_stop_time ) then
            open(newunit=file_handle, file='radial_profiles.'//tag)
            do n_r=1,n_r_max
               uphiR%SD(n_r)     =sqrt(uphiR%SD(n_r)/timeAvg_rad)
               fluxR%SD(n_r)     =sqrt(fluxR%SD(n_r)/timeAvg_rad)
               us2R%SD(n_r)      =sqrt(us2R%SD(n_r)/timeAvg_rad)
               up2R%SD(n_r)      =sqrt(up2R%SD(n_r)/timeAvg_rad)
               enstrophyR%SD(n_r)=sqrt(enstrophyR%SD(n_r)/timeAvg_rad)
               if ( l_heat ) tempR%SD(n_r)=sqrt(tempR%SD(n_r)/timeAvg_rad)
               if ( l_chem ) xiR%SD(n_r)=sqrt(xiR%SD(n_r)/timeAvg_rad)
               write(file_handle, '(es20.12, 14es16.8)') r(n_r),           &
               &     round_off(us2R%mean(n_r)), round_off(us2R%SD(n_r)),   &
               &     round_off(up2R%mean(n_r)), round_off(up2R%SD(n_r)),   &
               &     round_off(enstrophyR%mean(n_r)),                      &
               &     round_off(enstrophyR%SD(n_r)),                        &
               &     round_off(uphiR%mean(n_r)), round_off(uphiR%SD(n_r)), &
               &     round_off(tempR%mean(n_r)+tcond(n_r)),                &
               &     round_off(tempR%SD(n_r)),                             &
               &     round_off(xiR%mean(n_r)+xicond(n_r)),                 &
               &     round_off(xiR%SD(n_r)), round_off(fluxR%mean(n_r)),   &
               &     round_off(fluxR%SD(n_r))
            end do
            close(file_handle)
         end if

      end if

   end subroutine get_radial_averages
!------------------------------------------------------------------------------
   subroutine get_time_series(time, us_Mloc, up_Mloc, om_Mloc, temp_Mloc, &
              &               dtemp_Mloc, xi_Mloc, dxi_Mloc, us2_m, up2_m,&
              &               enstrophy_m, us2_r, up2_r, enstrophy, flux_r)

      !-- Input variables
      real(cp),    intent(in) :: time
      complex(cp), intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: temp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: dtemp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: xi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: dxi_Mloc(nMstart:nMstop,n_r_max)
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
      real(cp) :: chem_2D, chem_3D, ShTop, ShBot, beta_xi, Sh_vol, xiTop, xiBot
      real(cp) :: up2_axi_r(n_r_max), nu_vol_r(n_r_max), nu_cond_r(n_r_max)
      real(cp) :: buo_power(n_r_max), pump(n_r_max), tmp(n_r_max)
      real(cp) :: chem_power(n_r_max), sh_vol_r(n_r_max)
      real(cp) :: theta(n_r_max), sh_cond_r(n_r_max)
      real(cp) :: NuTop, NuBot, beta_t, Nu_vol, Nu_int, fac
      real(cp) :: rey_2D, rey_fluct_2D, rey_zon_2D
      real(cp) :: rey_3D, rey_fluct_3D, rey_zon_3D

      do n_r=1,n_r_max
         us2_r(n_r)     =0.0_cp
         up2_r(n_r)     =0.0_cp
         up2_axi_r(n_r) =0.0_cp
         enstrophy(n_r) =0.0_cp
         buo_power(n_r) =0.0_cp
         chem_power(n_r)=0.0_cp
         pump(n_r)      =0.0_cp
         nu_vol_r(n_r)  =0.0_cp
         sh_vol_r(n_r)  =0.0_cp
         flux_r(n_r)    =0.0_cp
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            us2_r(n_r)    =us2_r(n_r)+cc2real(us_Mloc(n_m,n_r),m)
            up2_r(n_r)    =up2_r(n_r)+cc2real(up_Mloc(n_m,n_r),m)
            enstrophy(n_r)=enstrophy(n_r)+cc2real(om_Mloc(n_m,n_r),m)

            if ( l_heat ) then
               buo_power(n_r)=buo_power(n_r)+cc22real(us_Mloc(n_m,n_r), &
               &                                      temp_Mloc(n_m,n_r), m)

               nu_vol_r(n_r) =nu_vol_r(n_r)+cc2real(dtemp_Mloc(n_m,n_r),m)&
               &              +real(m,cp)*real(m,cp)*or2(n_r)*            &
               &                             cc2real(temp_Mloc(n_m,n_r),m)

               flux_r(n_r)   =flux_r(n_r)+cc22real(us_Mloc(n_m,n_r), &
               &                                   temp_Mloc(n_m,n_r), m)
            end if

            if ( l_chem ) then
               chem_power(n_r)=chem_power(n_r)+cc22real(us_Mloc(n_m,n_r), &
               &                                        xi_Mloc(n_m,n_r), m)

               sh_vol_r(n_r) =sh_vol_r(n_r)+cc2real(dxi_Mloc(n_m,n_r),m)  &
               &              +real(m,cp)*real(m,cp)*or2(n_r)*            &
               &                             cc2real(xi_Mloc(n_m,n_r),m)
            end if

            if ( m == 0 ) then
               up2_axi_r(n_r)=up2_axi_r(n_r)+cc2real(up_Mloc(n_m,n_r),m)
               pump(n_r)     =pump(n_r)+cc2real(up_Mloc(n_m,n_r),m)
               if ( l_heat ) then
                  flux_r(n_r)   =flux_r(n_r)-TdiffFac*real(dtemp_Mloc(n_m,n_r))-&
                  &              TdiffFac*dtcond(n_r)
               end if
            end if
         end do
         if ( l_heat ) nu_vol_r(n_r)=nu_vol_r(n_r)*r(n_r)*height(n_r)
         if ( l_chem ) sh_vol_r(n_r)=sh_vol_r(n_r)*r(n_r)*height(n_r)
      end do

      !-- MPI reductions to get the s-profiles on rank==0
      call reduce_radial_on_rank(us2_r, 0)
      call reduce_radial_on_rank(up2_r, 0)
      call reduce_radial_on_rank(up2_axi_r, 0)
      call reduce_radial_on_rank(enstrophy, 0)
      call reduce_radial_on_rank(pump, 0)
      if ( l_heat ) then
         call reduce_radial_on_rank(buo_power, 0)
         call reduce_radial_on_rank(flux_r, 0)
         call reduce_radial_on_rank(nu_vol_r, 0)
      end if
      if ( l_chem ) then
         call reduce_radial_on_rank(sh_vol_r, 0)
         call reduce_radial_on_rank(chem_power, 0)
      end if

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
         !-- \int_v div ( \vec{\omega}\times\vec{u} ) = - \int_S \omega u_\phi dS
         !--                                          = - \int_S 2 u_\phi^2/s dS
         if ( kbotv == 1 ) then
            visc_2D = visc_2D+four*pi*up2_r(n_r_max)
         end if
         if ( ktopv == 1 ) then
            visc_2D = visc_2D-four*pi*up2_r(1)
         end if

         if ( l_heat ) then
            tmp(:)=BuoFac*buo_power(:)*rgrav(:)*r(:)
            pow_2D=rInt_R(tmp, r, rscheme)
            pow_2D=round_off(two*pi*pow_2D)
         else
            pow_2D=0.0_cp
         end if

         if ( l_chem ) then
            tmp(:)=ChemFac*chem_power(:)*rgrav(:)*r(:)
            chem_2D=rInt_R(tmp, r, rscheme)
            chem_2D=round_off(two*pi*chem_2D)
         else
            chem_2D=0.0_cp
         end if

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

            if ( l_heat ) then
               tmp(:)=BuoFac*buo_power(:)*rgrav(:)*r(:)*height(:)
               pow_3D=rInt_R(tmp, r, rscheme)
               pow_3D=round_off(two*pi*pow_3D)
            else
               pow_3D=0.0_cp
            end if

            if ( l_chem ) then
               tmp(:)=ChemFac*chem_power(:)*rgrav(:)*r(:)*height(:)
               chem_3D=rInt_R(tmp, r, rscheme)
               chem_3D=round_off(two*pi*chem_3D)
            else
               chem_3D=0.0_cp
            end if

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
         write(n_power_file_2D, '(1P, es20.12, 3es16.8)') time, pow_2D, chem_2D, &
         &                                                visc_2D

         write(n_rey_file_2D, '(1P, es20.12, 3es16.8)') time, rey_2D, rey_zon_2D,&
         &                                              rey_fluct_2D

         if ( .not. l_non_rot ) then
            write(n_kin_file_3D, '(1P, es20.12, 4es16.8)') time, us2_3D, up2_3D, &
            &                                              uz2_3D, up2_axi_3D
            write(n_power_file_3D, '(1P, es20.12, 4es16.8)') time, pow_3D, &
            &                                                chem_3D,      &
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
         write(n_lscale_file, '(1P, es20.12, 5es16.8)') time, dlus_peak,         &
         &                                              dlekin_peak, dlvort_peak,&
         &                                              dlekin_int, dl_diss

         ! At this stage multiply us2_r, up2_r and enstrophy by r(:) and height(:)
         us2_r(:)    =us2_r(:)*r(:)*height(:)
         up2_r(:)    =up2_r(:)*r(:)*height(:)
         enstrophy(:)=enstrophy(:)*r(:)*height(:)

      end if

      !-------
      !-- us uphi correlation
      !-------
      if ( l_corr ) call get_corr(time, us_Mloc, up_Mloc, us2_r, up2_r)


      if ( l_rank_has_m0 .and. l_heat ) then
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
         if (abs(dtcond(1)) > tiny_number) then
            NuTop = one+real(dtemp_Mloc(n_m0,1))/dtcond(1)
         else
            NuTop = one
         end if
         if (abs(dtcond(n_r_max)) > tiny_number) then
            NuBot = one+real(dtemp_Mloc(n_m0,n_r_max))/dtcond(n_r_max)
            !&       (dtcond(n_r_max)-tadvz_fac*beta(n_r_max)*tcond(n_r_max))
         else
            NuBot = one
         end if
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

      if ( l_rank_has_m0 .and. l_chem ) then
         !------
         !-- Chemical composition
         !------

         !-- Top and bottom temperatures
         n_m0 = m2idx(0)
         xiTop = real(xi_Mloc(n_m0,1))+xicond(1)
         xiBot = real(xi_Mloc(n_m0,n_r_max))+xicond(n_r_max)
         xiTop = round_off(xiTop)
         xiBot = round_off(xiBot)

         !-- Classical top and bottom Sherwood number
         ShTop = one+real(dxi_Mloc(n_m0,1))/dxicond(1)
         ShBot = one+real(dxi_Mloc(n_m0,n_r_max))/dxicond(n_r_max)
         !&       (dtcond(n_r_max)-tadvz_fac*beta(n_r_max)*tcond(n_r_max))
         ShTop = round_off(ShTop)
         ShBot = round_off(ShBot)

         !-- Volume-based Nusselt number
         Sh_vol = rInt_R(sh_vol_r, r, rscheme)
         sh_cond_r(:)=dxicond(:)*dxicond(:)*height(:)*r(:)
         Sh_vol = one+Sh_vol/rInt_R(sh_cond_r, r, rscheme)
         Sh_vol = round_off(Sh_vol)

         !-- Mid-shell temperature gradient
         beta_xi = dxicond(int(n_r_max/2))+real(dxi_Mloc(n_m0,int(n_r_max/2)))
         beta_xi = round_off(beta_xi)

         write(n_chem_file, '(1P, ES20.12, 6ES16.8)') time, ShTop, ShBot,   &
         &                                            Sh_vol, xiTop,        &
         &                                            xiBot, beta_xi
      end if

   end subroutine get_time_series
!------------------------------------------------------------------------------
   subroutine get_corr(time, us_Mloc, up_Mloc, us2_r, up2_r)
      !
      ! This routine calculates the following u_s u_\phi correlations:
      !    - stress = [ <u_s u_\phi> ]
      !    - corr = [<u_s u_\phi>] / [\sqrt{ <u_s^2> <u_\phi^2> } ]
      !
      ! where < > corresponds to a phi-average and [] corresponds to
      ! \int h*s ds
      !

      !-- Input variables:
      real(cp),    intent(in) :: time
      complex(cp), intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      real(cp),    intent(in) :: us2_r(n_r_max)
      real(cp),    intent(in) :: up2_r(n_r_max)

      !-- Local variables:
      integer :: n_m, m, n_r
      real(cp) :: tmp(n_r_max)
      real(cp) :: corr, stress, den


      do n_r=1,n_r_max
         tmp(n_r)=0.0_cp
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            tmp(n_r)=tmp(n_r)+r(n_r)*height(n_r)*&
            &        cc22real(us_Mloc(n_m,n_r),up_Mloc(n_m,n_r),m)
         end do
      end do

      !-- MPI reduction
      call reduce_radial_on_rank(tmp, 0)

      if ( rank == 0 ) then
         !-- Radial integration
         stress = rInt_R(tmp, r, rscheme)
         stress = two*pi*stress

         ! \int \sqrt{u_s^2*u_\phi^2} *h*s ds
         tmp(:)=sqrt(us2_r(:)*up2_r(:)*oheight(:)*oheight(:)*or2(:))*&
         &      r(:)*height(:)
         den = rInt_R(tmp, r, rscheme)
         den = two*pi*den

         if ( den > 0.0_cp ) then
            corr = stress/den
         else
            corr = one
         end if

         write(n_corr_file, '(1P, es20.12, 2es16.8)') time, stress, corr
      end if

   end subroutine get_corr
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
      lus_peak=my_allreduce_maxloc(us2_m(max(2,nMstart):nMstop),nMstart)
      lus_peak=idx2m(lus_peak)
      if ( lus_peak > 0 ) then
         dlus_peak = pi/lus_peak
      else
         dlus_peak = 0.0_cp
      end if

      lekin_peak=my_allreduce_maxloc(us2_m(max(2,nMstart):nMstop)+ &
                 &                   up2_m(max(2,nMstart):nMstop),nMstart)
      lekin_peak=idx2m(lekin_peak)
      if ( lekin_peak > 0 ) then
         dlekin_peak = pi/lekin_peak
      else
         dlekin_peak = 0.0_cp
      end if

      lvort_peak=my_allreduce_maxloc(enst_m(max(2,nMstart):nMstop),nMstart)
      lvort_peak=idx2m(lvort_peak)
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
