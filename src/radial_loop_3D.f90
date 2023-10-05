module rloop_3D

   use precision_mod
   use parallel_mod
   use constants, only: pi, ci, two, zero
   use nonlinear_lm_mod, only: nonlinear_lm_t
   use grid_space_arrays_mod, only: grid_space_arrays_t
   use blocking, only: nRstart3D, nRstop3D, nRstart, nRstop
   use namelists, only: BuoFac, DyMagFac, tag, l_heat_3D, l_thw_3D, l_mag_3D, l_mag_LF, &
   &                    l_cyl, l_QG_basis, r_cmb, r_icb, l_lin_solve, l_leibniz,        &
   &                    beta_shift, l_mag_B0, l_b_phiavg
   use truncation_3D, only: lm_max, lmP_max, n_phi_max_3D, n_theta_max, l_max, n_r_max_3D
   use truncation, only: idx2m, n_m_max, n_r_max
   use communications, only: scatter_from_rank0_to_rloc
   use courant_mod, only: courant_3D
#ifdef WITH_SHTNS
   use shtns, only: scal_to_SH, scal_to_spat, torpol_to_spat, spat_to_qst, &
       &            torpol_to_curl_spat, scal_axi_to_grad_spat
#endif
   use horizontal, only: cost, sint
   use radial_functions, only: r, r_3D, or1, beta
   use z_functions, only: zfunc_type
   use timers_mod, only: timers_type
   use time_schemes, only: type_tscheme
   use output_frames, only: open_snapshot_3D, close_snapshot_3D, &
       &                    write_bulk_snapshot_3D
   use fields, only: B0r_3D_Rloc, B0t_3D_Rloc, B0p_3D_Rloc!, curlB0r_3D_Rloc, &
   !    &             curlB0t_3D_Rloc, curlB0p_3D_Rloc

   implicit none

   private

   type(grid_space_arrays_t) :: gsa
   type(nonlinear_lm_t) :: nl_lm
   integer :: frame_counter=1
   integer :: n_bavg_file
   integer :: n_r_bavg, n_bavg_rank

   public :: radial_loop_3D, initialize_radial_loop_3D, finalize_radial_loop_3D

contains

   subroutine initialize_radial_loop_3D(lm_max, lmP_max)

      integer, intent(in) :: lm_max
      integer, intent(in) :: lmP_max

      character(len=144) :: file_name

      call gsa%initialize()
      call nl_lm%initialize(lm_max,lmP_max)

      n_r_bavg = int(n_r_max_3D*0.08) ! 92\% of the 3D-radius for butterfly diagram
      if( (n_r_bavg>=nrStart3D) .and. (n_r_bavg<=nrStop3D) ) then ! select processor
         n_bavg_rank = rank
      else
         n_bavg_rank = -1
      end if

      if ( rank == n_bavg_rank ) then
         if ( l_mag_3D .and. l_b_phiavg ) then
            file_name = 'b_phiavg.'//tag
            open(newunit=n_bavg_file, file=file_name, status='new')
         end if
       end if

   end subroutine initialize_radial_loop_3D
!------------------------------------------------------------------------------
   subroutine finalize_radial_loop_3D

      call nl_lm%finalize()
      call gsa%finalize()
      if ( rank == n_bavg_rank ) then
         if ( l_mag_3D .and. l_b_phiavg ) close(n_bavg_file)
      end if

   end subroutine finalize_radial_loop_3D
!------------------------------------------------------------------------------
   subroutine radial_loop_3D( time, ur, ut, up, temp, dtempdt, dVrTLM,  &
              &               buo_Rloc, b_3D, db_3D, ddb_3D, aj_3D,     &
              &               dj_3D, dbdt_3D, djdt_3D, lf_Rloc,         &
              &               djxB_Rloc, dVxBhLM, dpsidt_Rloc,          &
              &               dtr_3D_Rloc, dth_3D_Rloc, l_frame,        &
              &               zinterp, timers, tscheme )

      !-- Input variables
      complex(cp), intent(inout) :: temp(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(inout) :: b_3D(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(inout) :: db_3D(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(inout) :: ddb_3D(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(inout) :: aj_3D(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(inout) :: dj_3D(lm_max, nRstart3D:nRstop3D)
      real(cp),    intent(inout) :: ur(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp),    intent(inout) :: ut(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp),    intent(inout) :: up(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      logical,     intent(in) :: l_frame
      real(cp),    intent(in) :: time
      type(zfunc_type),    intent(in) :: zinterp
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp),       intent(out) :: dtempdt(lm_max, nRstart3D:nRstop3D)
      complex(cp),       intent(out) :: dVrTLM(lm_max, nRstart3D:nRstop3D)
      complex(cp),       intent(out) :: dbdt_3D(lm_max, nRstart3D:nRstop3D)
      complex(cp),       intent(out) :: djdt_3D(lm_max, nRstart3D:nRstop3D)
      complex(cp),       intent(out) :: dVxBhLM(lm_max, nRstart3D:nRstop3D)
      complex(cp),       intent(inout) :: buo_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),       intent(inout) :: lf_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),       intent(inout) :: djxB_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),       intent(inout) :: dpsidt_Rloc(n_m_max, nRstart:nRstop)
      real(cp),          intent(out) :: dtr_3D_Rloc(nRstart3D:nRstop3D)
      real(cp),          intent(out) :: dth_3D_Rloc(nRstart3D:nRstop3D)
      type(timers_type), intent(inout) :: timers

      !-- Local variables
      real(cp) :: buo_tmp(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp) :: jxBs(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp) :: jxBp(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp) :: jxBz(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp) :: djxBpdz(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp) :: djxBpds(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp) :: dTdth(n_theta_max,nRstart3D:nRstop3D)
      real(cp) :: runStart, runStop!, phi, theta!, rsint
      complex(cp) :: tmp_2D(n_m_max,n_r_max)!nRstart:nRstop)!
      complex(cp) :: buo_tmp_Rloc(n_m_max,nRstart:nRstop)
      complex(cp) :: lfs_tmp_Rloc(n_m_max,nRstart:nRstop)
      complex(cp) :: lfp_tmp_Rloc(n_m_max,nRstart:nRstop)
      complex(cp) :: lfz_tmp_Rloc(n_m_max,nRstart:nRstop)
      complex(cp) :: lfpds_tmp_Rloc(n_m_max,nRstart:nRstop)
      complex(cp) :: lfpdz_tmp_Rloc(n_m_max,nRstart:nRstop)
      complex(cp) :: lfp_bc_tmp_Rloc(n_m_max,nRstart:nRstop)
      complex(cp) :: lfp_h_tmp_Rloc(n_m_max,nRstart:nRstop)!n_r_max)
      integer :: n_r, n_m, m, fh_temp, info_temp, n_theta!, n_phi
      integer :: fh_ur, info_ur, fh_ut, info_ut, fh_up, info_up
      integer :: fh_br, info_br, fh_bt, info_bt, fh_bp, info_bp
      integer :: fh_jxbs, info_jxbs, fh_jxbp, info_jxbp, fh_jxbz, info_jxbz
      logical :: l_jxb_save
      character(len=144) :: frame_name

      l_jxb_save=.False.
      buo_tmp(:,:,:) = 0.0_cp
      !-- get thermal wind
      if ( l_heat_3D .and. l_thw_3D .and. (.not. l_lin_solve) ) then
         do n_r=nRstart3D,nRstop3D
            !-- Take the gradient of the axisymmetric part
            call transform_axi_to_dth_grid_space(temp(1:l_max+1,n_r), dTdth(:,n_r))
         end do

         !-- Compute thermal wind --> modify up_3D_Rloc
         call zinterp%compute_thermal_wind(dTdth, up)
      end if

      !-- Open the 3-D snapshots in physical space
      if ( l_frame .and. tscheme%istage==1 ) then
         write(frame_name, '(A,I0,A,A)') 'frame_ur_3D_',frame_counter,'.',tag 
         call open_snapshot_3D(frame_name, time, fh_ur, info_ur)
         write(frame_name, '(A,I0,A,A)') 'frame_ut_3D_',frame_counter,'.',tag
         call open_snapshot_3D(frame_name, time, fh_ut, info_ut)
         write(frame_name, '(A,I0,A,A)') 'frame_up_3D_',frame_counter,'.',tag
         call open_snapshot_3D(frame_name, time, fh_up, info_up)
         if ( l_heat_3D ) then
            write(frame_name, '(A,I0,A,A)') 'frame_temp_3D_',frame_counter,'.',tag
            call open_snapshot_3D(frame_name, time, fh_temp, info_temp)
         end if
         if ( l_mag_3D ) then
            write(frame_name, '(A,I0,A,A)') 'frame_br_3D_',frame_counter,'.',tag
            call open_snapshot_3D(frame_name, time, fh_br, info_br)
            write(frame_name, '(A,I0,A,A)') 'frame_bt_3D_',frame_counter,'.',tag
            call open_snapshot_3D(frame_name, time, fh_bt, info_bt)
            write(frame_name, '(A,I0,A,A)') 'frame_bp_3D_',frame_counter,'.',tag
            call open_snapshot_3D(frame_name, time, fh_bp, info_bp)
            if ( l_mag_LF .and. l_jxb_save ) then
               write(frame_name, '(A,I0,A,A)') 'frame_jxbs_3D_',frame_counter,'.',tag
               call open_snapshot_3D(frame_name, time, fh_jxbs, info_jxbs)
               write(frame_name, '(A,I0,A,A)') 'frame_jxbp_3D_',frame_counter,'.',tag
               call open_snapshot_3D(frame_name, time, fh_jxbp, info_jxbp)
               write(frame_name, '(A,I0,A,A)') 'frame_jxbz_3D_',frame_counter,'.',tag
               call open_snapshot_3D(frame_name, time, fh_jxbz, info_jxbz)
            end if
         end if
      end if

      !ur(:,:,:) = 0.0_cp
      !ut(:,:,:) = 0.0_cp
      !up(:,:,:) = 0.0_cp
      runStart = MPI_Wtime()
      do n_r=nRstart3D,nRstop3D
      !!if ( frame_counter == 1 ) then
      !      do n_theta=1,n_theta_max
      !         do n_phi=1,n_phi_max_3D
      !            phi = (n_phi-1)*two*pi/(n_phi_max_3D)
      !            !ur(n_phi,n_theta,n_r) = sin(pi*(r_3D(n_r)-r_icb))*sint(n_theta)*cos(4*phi+pi*0.25_cp)
      !            !ut(n_phi,n_theta,n_r) = sin(pi*(r_3D(n_r)-r_icb))*cost(n_theta)
      !            !up(n_phi,n_theta,n_r) = sin(pi*(r_3D(n_r)-r_icb))*sint(n_theta)*cos(4*phi)
      !            if ( r_3D(n_r)*sint(n_theta) > r_icb ) then
      !               !ur(n_phi,n_theta,n_r) = sin(pi*(r_3D(n_r)-r_icb))*sin(pi*(r_3D(n_r)*sint(n_theta)-r_icb))!*sint(n_theta)*cos(4*phi+pi*0.25_cp)
      !               !ut(n_phi,n_theta,n_r) = sin(pi*(r_3D(n_r)-r_icb))*sin(pi*(r_3D(n_r)*sint(n_theta)-r_icb))!*cost(n_theta)
      !               !up(n_phi,n_theta,n_r) = sin(pi*(r_3D(n_r)-r_icb))*sin(pi*(r_3D(n_r)*sint(n_theta)-r_icb))!*sint(n_theta)*cos(4*phi)
      !            ur(n_phi,n_theta,n_r) = 10*sint(n_theta)*sin(pi*(r_3D(n_r)-r_icb))*sin(5*pi*(r_3D(n_r)*sint(n_theta)-r_icb))* &
      !            &                       cos(8*phi+pi*0.25_cp)*r_3D(n_r)*sint(n_theta)*(r_cmb-r_3D(n_r)*sint(n_theta))/r_cmb**2
      !            ut(n_phi,n_theta,n_r) = 10*cost(n_theta)*sin(3*pi*(r_3D(n_r)*sint(n_theta)-r_icb))*cos(8*phi)* &
      !            &                       r_3D(n_r)*sint(n_theta)*(r_cmb-r_3D(n_r)*sint(n_theta))/r_cmb**2
      !            up(n_phi,n_theta,n_r) =-10*sin(2*pi*(r_3D(n_r)*sint(n_theta)-r_icb))*4*r_3D(n_r)*sint(n_theta) * &
      !            &                       ( r_cmb - r_3D(n_r)*sint(n_theta) )/r_cmb**2
      !            end if
      !        end do
      !      end do
      !!end if

         !-- Transform temperature and magnetic-field from (l,m) to (theta,phi)
         call transform_to_grid_space(temp(:,n_r), b_3D(:,n_r),   &
              &                      db_3D(:,n_r), ddb_3D(:,n_r), &
              &                      aj_3D(:,n_r),  dj_3D(:,n_r), &
              &                      n_r, gsa)

            !!-- Analytical B and VxB to benchmark Lorentz force if needed
            !do n_theta=1,n_theta_max
            !   do n_phi=1,n_phi_max_3D
            !      phi = (n_phi-1)*two*pi/(n_phi_max_3D)
            !!      gsa%Brc(n_phi,n_theta) = r_3D(n_r)*cos(phi)**2.*cost(n_theta)
            !!      gsa%Btc(n_phi,n_theta) = r_3D(n_r)*sin(phi)**2.
            !!      gsa%Bpc(n_phi,n_theta) =-r_3D(n_r)*sin(phi)**2.*cost(n_theta)
            !!      gsa%curlBrc(n_phi,n_theta) = ( (r_3D(n_r)*sin(phi)**2.*sint(n_theta)**2. &
            !!      &                             - r_3D(n_r)*sin(phi)**2.*cost(n_theta)**2.) &
            !!      &                            - 2.0_cp*r_3D(n_r)*sin(phi)*cos(phi) )/(r_3D(n_r)*sint(n_theta))
            !!      gsa%curlBtc(n_phi,n_theta) = ( (-r_3D(n_r)*2.0_cp*cos(phi)*sin(phi)*cost(n_theta))/sint(n_theta) &
            !!      &                            + 2.0_cp*r_3D(n_r)*sin(phi)**2.*cost(n_theta) )/r_3D(n_r)
            !!      gsa%curlBpc(n_phi,n_theta) = ( 2.0_cp*r_3D(n_r)*sin(phi)**2. &
            !!      &                            + r_3D(n_r)*cos(phi)**2.*sint(n_theta) )/r_3D(n_r)
            !      theta = acos(cost(n_theta))
            !      gsa%Brc(n_phi,n_theta) = 5.0_cp/8.0_cp*(8.0_cp*r_cmb-6.0_cp*r_3D(n_r)-2.0_cp*r_icb**4/r_3D(n_r)**3) &
            !      &                       *cost(n_theta)
            !      gsa%Btc(n_phi,n_theta) = 5.0_cp/8.0_cp*(9.0_cp*r_3D(n_r)-8.0_cp*r_cmb-r_icb**4/r_3D(n_r)**3) &
            !      &                       *sint(n_theta)
            !      gsa%Bpc(n_phi,n_theta) = 5.0_cp*sin(pi*(r_3D(n_r)-r_icb))*sin(2.0_cp*theta)
            !      gsa%curlBrc(n_phi,n_theta) = 5.0_cp/r_3D(n_r)*sin(pi*(r_3D(n_r)-r_icb))*(2.0_cp*cos(2.0_cp*theta) &
            !      &                           +sin(2.0_cp*theta)*cost(n_theta)/sint(n_theta))
            !      gsa%curlBtc(n_phi,n_theta) =-5.0_cp*(1.0_cp/r_3D(n_r)*sin(pi*(r_3D(n_r)-r_icb)) &
            !      &                           + pi*cos(pi*(r_3D(n_r)-r_icb)))*sin(2.0_cp*theta)
            !      gsa%curlBpc(n_phi,n_theta) = 15.0_cp/2.0_cp*sint(n_theta)
            !   end do
            !end do

         !print*, '!---------------------------', 'NEW LINE', '---------------------------!'
         !print*, 'alphaFac', gsa%Alphac(2,2), gsa%Alphac(5,7)
         !print*, 'br; btheta; bphi', gsa%Brc(2,2), '; ', gsa%Btc(2,2), '; ', gsa%Bpc(2,2)
         !print*, 'b0r; b0theta; b0phi', B0r_3D_Rloc(2,2,n_r), '; ', B0t_3D_Rloc(2,2,n_r), '; ', B0p_3D_Rloc(2,2,n_r)
         !-- Courant condition
         if ( l_mag_3D .and. tscheme%istage == 1 ) then
            if ( .not. l_mag_B0 ) then
            call courant_3D(n_r, dtr_3D_Rloc(n_r), dth_3D_Rloc(n_r),  &
                 &          ur(:,:,n_r),  ut(:,:,n_r),  up(:,:,n_r),  &
                 &          gsa%Brc(:,:), gsa%Btc(:,:), gsa%Bpc(:,:), &
                 &          gsa%Alphac(:,:), tscheme%courfac,         &
                 &          tscheme%alffac)
            else
               call courant_3D(n_r, dtr_3D_Rloc(n_r), dth_3D_Rloc(n_r),  &
                    &          ur(:,:,n_r),  ut(:,:,n_r),  up(:,:,n_r),  &
                    &          gsa%Brc(:,:)+B0r_3D_Rloc(:,:,n_r),        &
                    &          gsa%Btc(:,:)+B0t_3D_Rloc(:,:,n_r),        &
                    &          gsa%Bpc(:,:)+B0p_3D_Rloc(:,:,n_r),        &
                    &          gsa%Alphac(:,:), tscheme%courfac,         &
                    &          tscheme%alffac)
            end if
         end if


         !-- Construct non-linear terms in physical space
         call gsa%get_nl(ur(:,:,n_r), ut(:,:,n_r), up(:,:,n_r), n_r, &
              &          buo_tmp(:,:,n_r),jxBs(:,:,n_r),jxBp(:,:,n_r),jxBz(:,:,n_r))

         !-- Compute and Write phi-Average of the 3D B field just below CMB
         if ( tscheme%istage==1 ) then
            if ( l_b_phiavg .and. n_r==n_r_bavg ) then
               !call output_local_B_cmb(time, rank, gsa%Brc(:,:), gsa%Bpc(:,:))
               if (rank == n_bavg_rank ) call output_local_B_cmb(time, gsa%Brc(:,:), gsa%Bpc(:,:))
            end if
         end if

         !-- Write the snapshot on the grid (easier to handle)
         if ( l_frame .and. tscheme%istage==1 ) then
#ifdef aDEBUG
         block
            integer :: filehandle
            if ( n_r == 12 ) then
               !print*, r_3D(n_r)
               !         jxBs(n_phi,n_theta)=asign*delta_fac*sqrt(vpfluct(n_theta))*Bsavg(n_theta)
               !         jxBp(n_phi,n_theta)=asign*delta_fac*sqrt(vsfluct(n_theta))*Bpavg(n_theta)
               open(newunit=filehandle, file='test_VxB', form='unformatted', access='stream')
                  write(filehandle) jxBs(:,:,n_r), jxBp(:,:,n_r)
               close(filehandle)
            end if
         end block
#endif
            call write_bulk_snapshot_3D(fh_ur, ur(:,:,n_r))
            call write_bulk_snapshot_3D(fh_ut, ut(:,:,n_r))
            call write_bulk_snapshot_3D(fh_up, up(:,:,n_r))
            if ( l_heat_3D ) then
               call write_bulk_snapshot_3D(fh_temp, gsa%Tc(:,:))!buo_tmp(:,:,n_r))!
            end if
            if ( l_mag_3D ) then
               !if ( .not. l_mag_B0 ) then
                  call write_bulk_snapshot_3D(fh_br, gsa%Brc(:,:))!jxBs(:,:,n_r))!Vx!
                  call write_bulk_snapshot_3D(fh_bt, gsa%Btc(:,:))!Vx!
                  call write_bulk_snapshot_3D(fh_bp, gsa%Bpc(:,:))!jxBp(:,:,n_r))!Vx!
               !else
               !   call write_bulk_snapshot_3D(fh_br, B0r_3D_Rloc(:,:,n_r))!jxBs(:,:,n_r))!Vx!gsa%Brc(:,:)+
               !   call write_bulk_snapshot_3D(fh_bt, B0t_3D_Rloc(:,:,n_r))!Vx!gsa%Btc(:,:)+
               !   call write_bulk_snapshot_3D(fh_bp, B0p_3D_Rloc(:,:,n_r))!jxBp(:,:,n_r))!Vx!gsa%Bpc(:,:)+
               !end if
               if ( l_mag_LF .and. l_jxb_save ) then
                  call write_bulk_snapshot_3D(fh_jxbs, gsa%VxBr(:,:))!jxBs(:,:,n_r))
                  call write_bulk_snapshot_3D(fh_jxbp, gsa%VxBp(:,:))!jxBp(:,:,n_r))
                  call write_bulk_snapshot_3D(fh_jxbz, gsa%VxBt(:,:))!jxBz(:,:,n_r))
               end if
            end if
         end if

         !-- Transform back the non-linear terms to (l,m) space
         call transform_to_lm_space(gsa, nl_lm, dVrTLM(:,n_r))

         !-- Get theta and phi derivatives using recurrence relations
         call nl_lm%get_td(dtempdt(:,n_r),dVxBhLM(:,n_r),dbdt_3D(:,n_r), &
              &            djdt_3D(:,n_r), n_r )

            !!-- Analytical jxB to benchmark Lorentz force if needed
            !jxBp(:,:,n_r) = 0.0_cp
            !jxBs(:,:,n_r) = 0.0_cp
            !jxBz(:,:,n_r) = 0.0_cp
            !do n_theta=1,n_theta_max
            !   !if ( r_3D(n_r)*sint(n_theta) >= r_icb ) then
            !   do n_phi=1,n_phi_max_3D
            !      phi = (n_phi-1)*two*pi/(n_phi_max_3D)
            !      !jxBs(n_phi,n_theta,n_r) = sin(phi)**2 * (r_3D(n_r)*sint(n_theta))*(r_cmb-r_3D(n_r)*sint(n_theta))* &
            !      !&                        (r_3D(n_r)*sint(n_theta)-r_icb) * (r_3D(n_r)*cost(n_theta))**2
            !      jxBs(n_phi,n_theta,n_r) = sin(phi)**2*(r_3D(n_r)*sint(n_theta))*(r_3D(n_r)*cost(n_theta))**2
            !      !jxBs(n_phi,n_theta,n_r) = (r_3D(n_r)*sint(n_theta))*exp(-pi*r_3D(n_r)*cost(n_theta))
            !      !if ( l_QG_basis ) jxBs(n_phi,n_theta,n_r) = jxBs(n_phi,n_theta,n_r) -                               &
            !      if ( l_QG_basis ) jxBz(n_phi,n_theta,n_r) =                              &
            !      !&                 sin(phi)*(r_3D(n_r)*sint(n_theta))!sinphi*s*z/beta * beta/z !directly
            !      &                 sin(phi)**2 * (r_3D(n_r)*sint(n_theta))*(r_cmb-r_3D(n_r)*sint(n_theta))*          &
            !      !&  (r_3D(n_r)*sint(n_theta)-r_icb) * (r_3D(n_r)*cost(n_theta))* &!*beta(n_r)*(r_3D(n_r)*cost(n_theta))!
            !      !&  (-r_3D(n_r)*sint(n_theta)/(r_cmb**2 - (r_3D(n_r)*cost(n_theta))**2))*r_3D(n_r)*cost(n_theta)
            !      !!&  (-r_3D(n_r)*sint(n_theta)/((r_cmb+beta_shift)**2 - (r_3D(n_r)*cost(n_theta))**2))*r_3D(n_r)*cost(n_theta)
            !      &  (r_3D(n_r)*sint(n_theta)-r_icb) * (r_3D(n_r)*cost(n_theta))* (r_3D(n_r)*cost(n_theta))!z/beta * beta !directly
            !      jxBp(n_phi,n_theta,n_r) =-(0.5_cp+cos(phi)+sin(4.*phi)) * (r_cmb**2-(r_3D(n_r)*sint(n_theta))**2)**2 * &
            !      &                         cos(pi/2*r_3D(n_r)*cost(n_theta)/sqrt(r_cmb**2-(r_3D(n_r)*sint(n_theta))**2))* &!r_icb)* &!sqrt(r_cmb**2-r_icb**2))* &!
            !      &                         pi/(r_3D(n_r)*sint(n_theta)*sqrt(r_cmb**2-(r_3D(n_r)*sint(n_theta))**2))!r_icb)!sqrt(r_cmb**2-r_icb**2))!
            !      jxBp(n_phi,n_theta,n_r) = (0.5_cp+cos(phi)+sin(4.*phi))*(r_3D(n_r)*sint(n_theta))*exp(-pi*r_3D(n_r)*cost(n_theta))
            !      !jxBp(n_phi,n_theta,n_r) =  (r_3D(n_r)*cost(n_theta))*sin(pi*r_3D(n_r)*cost(n_theta))
            !   end do
            !   !end if
            !end do

         if ( n_r==1 .and. (l_mag_LF .and. (.not. l_lin_solve) .and. l_leibniz) ) then  !-- Surface terms at CMB? --> No need if potential field
            call zinterp%interp_1d(n_r,jxBp(:,:,n_r), tmp_2D)
         end if
      end do

      runStop = MPI_Wtime()
      if (runStop>runStart) then
         timers%n_r_loops_3D=timers%n_r_loops_3D+1
         timers%r_loop_3D   =timers%r_loop_3D+(runStop-runStart)
      end if

      !--Close the 3-D snapshots in physical space
      if ( l_frame .and. tscheme%istage==1 ) then
         call close_snapshot_3D(fh_ur, info_ur)
         call close_snapshot_3D(fh_ut, info_ut)
         call close_snapshot_3D(fh_up, info_up)
         if ( l_heat_3D ) then
            call close_snapshot_3D(fh_temp, info_temp)
         end if
         if ( l_mag_3D ) then
            call close_snapshot_3D(fh_br, info_br)
            call close_snapshot_3D(fh_bt, info_bt)
            call close_snapshot_3D(fh_bp, info_bp)
            if ( l_mag_LF .and. l_jxb_save ) then
               call close_snapshot_3D(fh_jxbs, info_jxbs)
               call close_snapshot_3D(fh_jxbp, info_jxbp)
               call close_snapshot_3D(fh_jxbz, info_jxbz)
            end if
         end if
         frame_counter = frame_counter+1
      end if

      !-- Compute z-derivative of the \phi-part of lorentz-force
      if (  l_mag_LF .and. l_QG_basis .and. (.not. l_lin_solve) ) then
         runStart = MPI_Wtime()
         call zinterp%compute_zder(jxBp, djxBpdz)
         runStop = MPI_Wtime()
         if (runStop>runStart) then
            timers%interp = timers%interp+(runStop-runStart)
            timers%n_zder=timers%n_zder+1
            timers%zder = timers%zder+(runStop-runStart)
         end if

         do n_r=nRstart3D,nRstop3D
            do n_theta=1,n_theta_max
               djxBpdz(:,n_theta,n_r) = r_3D(n_r)*cost(n_theta)*djxBpdz(:,n_theta,n_r)
            end do
         end do
         !if ( rank == 0 ) print*, "After z_der:djxBpdz",  djxBpdz(1,1:10,2)
      end if

      !-- Compute z-averaging of buoyancy and lorentz force
      runStart = MPI_Wtime()
      !!call zinterp%compute_zder(jxBs, djxBpds)
      !!call zinterp%compute_sder(jxBs, djxBpds)
      !!call zinterp%compute_zavg(jxBp, lfp_tmp_Rloc,2,.false.)
      if ( l_heat_3D ) call zinterp%compute_zavg(buo_tmp, buo_tmp_Rloc,2,.false.)
      if ( l_mag_LF .and. (.not. l_lin_solve) ) then  !-- Non-Linear terms?
         lf_Rloc(:,:)=zero !-- Initialisation of z-avg Lorentz force
         call zinterp%compute_zavg(jxBs, lfs_tmp_Rloc,2,.true.)
         call zinterp%compute_zavg(jxBp, lfp_tmp_Rloc,2,.false.)!0)!
         if ( l_leibniz ) then
            djxB_Rloc(:,:)=zero
            call scatter_from_rank0_to_rloc(tmp_2D, lfp_h_tmp_Rloc, n_m_max)!, rank)
            !print*, "After z_interp1d:: N_rank, lfp_h_tmp_Rloc", rank,"::", real(lfp_h_tmp_Rloc(1,:))
            if ( l_cyl ) then
               lfp_bc_tmp_Rloc(:,:) = lfp_tmp_Rloc(:,:)
            else !-- set the boundary points to 0 when using low precision z-avg scheme
               call zinterp%compute_zavg(jxBp, lfp_bc_tmp_Rloc,2,.true.)!0)!
            end if
         else !-- No leibniz rule
            do n_r=nRstart3D,nRstop3D
               do n_theta=1,n_theta_max
                  jxBp(:,n_theta,n_r) = r_3D(n_r)*sint(n_theta)*jxBp(:,n_theta,n_r)
               end do
            end do
            call zinterp%compute_sder(jxBp, djxBpds)
            call zinterp%compute_zavg(djxBpds, lfpds_tmp_Rloc,2,.true.)!lfp_tmp_Rloc,1)!
         end if
         if ( l_QG_basis ) then
            call zinterp%compute_zavg(jxBz, lfz_tmp_Rloc,2,.true.)
            call zinterp%compute_zavg(djxBpdz, lfpdz_tmp_Rloc,2,.true.)
         end if
      end if
      runStop = MPI_Wtime()
      if (runStop>runStart) then
         timers%interp = timers%interp+(runStop-runStart)
         timers%n_zavg=timers%n_zavg+1
         timers%zavg = timers%zavg+(runStop-runStart)
      end if

#ifdef aDEBUG
      block
         use truncation, only: n_r_max
         use truncation_3D, only: n_r_max_3D
         use radial_functions, only: r, r_3D

         integer :: file_handle, n_r

         open(newunit=file_handle, file='buo', status='new', form='formatted')
         do n_r=1,n_r_max
            write(file_handle, *) r(n_r),  real(buo_tmp_Rloc(:,n_r))
         end do

         !open(newunit=file_handle, file='buo_3D', status='new', form='unformatted',&
         !&    access='stream')
         !write(file_handle) buo_tmp
         !close(file_handle)
      end block
#endif

      !-- Finish assembling buoyancy and sum it with dpsidt
      !DyMagFac = 1.0_cp
      !lf_Rloc(:,:)=zero !--> should be encapsulated in l_mag_LF -> Done!
      !djxB_Rloc(:,:)=zero !--> should be encapsulated in l_leibniz -> Done!
      !dpsidt_Rloc(:,:)=zero
      do n_m=1,n_m_max
         m = idx2m(n_m)
         do n_r=nRstart,nRstop
            if ( l_heat_3D ) then!.and. ( .not. l_mag_B0 ) ) then
            !-- Finish assembly the buoyancy:: no 1/s term because a s terms arises from Vx[Tg/r(s.e_s + z.e_z)].e_z
               dpsidt_Rloc(n_m,n_r)= dpsidt_Rloc(n_m,n_r)-BuoFac*ci*real(m,cp)* &
               &                     buo_tmp_Rloc(n_m,n_r)!!*or1(n_r)
               !-- store z-avg 3D-buo for Force balance outputs! W: spurious at the moment
               buo_Rloc(n_m,n_r)=zero!-BuoFac*ci*real(m,cp)*buo_tmp_Rloc(n_m,n_r)
            end if
            if ( l_mag_LF .and. (.not. l_lin_solve) ) then
            !-- Finish assembling both r.h.s. parts of lorentz force and sum it with dpsidt
            !--   --> <Vx(jxB).e_z> = (+<d_s (s (jxB)_p)> - d_p <(jxB)_s>)/s
               if ( m == 0 ) then
                  !-- LF on m=0 --> <(jxB).e_p> = <(jxB)_p>
                  if ( n_r==1 .or. n_r==n_r_max ) lfp_tmp_Rloc(n_m,n_r) = zero
                  dpsidt_Rloc(n_m,n_r)=dpsidt_Rloc(n_m,n_r)+DyMagFac*lfp_tmp_Rloc(n_m,n_r)
                  !-- store z-avg 3D-lorentz for Force balance outputs on m=0
                  lf_Rloc(n_m,n_r)=DyMagFac*lfp_tmp_Rloc(n_m,n_r)
               else
                  !-- LF on m<>0 --> <Vx(jxB).e_z> = (+<d_s (s (jxB)_p)> - d_p <(jxB)_s>)/s
                  if ( l_leibniz ) then
                     djxB_Rloc(n_m,n_r)=DyMagFac*r(n_r)*lfp_tmp_Rloc(n_m,n_r)
                     !lf_Rloc(n_m,n_r)=DyMagFac*(beta(n_r)*r(n_r)*(lfp_bc_tmp_Rloc(n_m,n_r)-lfp_h_tmp_Rloc(n_m,n_r)))!))!
                     if ( n_r == 1 ) then !--> "unshifting beta" for a better Leibniz's rule accuracy
                        lf_Rloc(n_m,n_r)=DyMagFac*(-r(n_r)/(r_cmb**2-r(n_r)**2+epsilon(1.0_cp)*(r(n_r)-r_icb)/r_icb) &
                        &                   *r(n_r)*(lfp_bc_tmp_Rloc(n_m,n_r)-lfp_h_tmp_Rloc(n_m,n_r)))!))!
                     else
                        lf_Rloc(n_m,n_r)=DyMagFac*(-r(n_r)/(r_cmb**2-r(n_r)**2)*r(n_r)* &
                        &                  (lfp_bc_tmp_Rloc(n_m,n_r)-lfp_h_tmp_Rloc(n_m,n_r)))!))!
                     end if
                     lf_Rloc(n_m,n_r)=lf_Rloc(n_m,n_r) - DyMagFac*ci*real(m,cp)*lfs_tmp_Rloc(n_m,n_r)
                  else
                     lf_Rloc(n_m,n_r)=DyMagFac*(lfpds_tmp_Rloc(n_m,n_r)           &
                     &                         -ci*real(m,cp)*lfs_tmp_Rloc(n_m,n_r))*or1(n_r)
                     !!-- if LF benchmark on m<>0
                     !dpsidt_Rloc(n_m,n_r)=DyMagFac*(lfpds_tmp_Rloc(n_m,n_r)           &
                     !&                             -ci*real(m,cp)*lfs_tmp_Rloc(n_m,n_r))*or1(n_r)
                     !-- store z-avg 3D-lorentz for Force balance outputs on m<>0 (no need if NOT Leibniz because nothing directly added to dpsidt)
                  !   !djxB_Rloc(n_m,n_r)=DyMagFac*(lfpds_tmp_Rloc(n_m,n_r)-ci*real(m,cp)*lfs_tmp_Rloc(n_m,n_r))*or1(n_r)
                  end if
                  !!-- LF WRONG implementation before --> (--d_s (s <(jxB)_p>) - d_p <(jxB)_s>)/s
                  !!-- d_p <(jxB)_s>)/s part --> (d_s (-dVsom)/s is added later)
                  !lf_Rloc(n_m,n_r)= lf_Rloc(n_m,n_r)-DyMagFac*ci*real(m,cp)* &
                  !&                    lfs_tmp_Rloc(n_m,n_r)*or1(n_r)
                  !-- --d_s (s <(jxB)_p>) part without the LEIBNIZ extra term from switching d_s . with <.>
                  !djxB_Rloc(n_m,n_r)= -DyMagFac*r(n_r)*&
                  !&                   lfp_tmp_Rloc(n_m,n_r)
                  !-- store z-avg 3D-lorentz for Leibniz integration boundary in update_psi
                  !--   --> \beta s ( <(jxB)_p> - 0.5 [(jxB)_p(h) + (jxB)_p(-h)] ) / s
                  !djxB_Rloc(n_m,n_r)= DyMagFac*(beta(n_r)*(lfp_tmp_Rloc(n_m,n_r)))!-lfp_h_tmp_Rloc(n_m,n_r)))
                  if ( l_QG_basis ) then
                  !-- Finish assembling add. lorentz force from the QG basis projection
                  !-- --> \beta <zVx(jxB).e_s> = \beta (+(d_p <z (jxB)_z>)/s - <z d_z (jxB)_p>)
                  !-- But --> \beta <z (jxB)_z> has already been added to lfs_tmp_Rloc !-- PROBLEMS with this
                  !--     INSTEAD: --> only <z (jxB)_z> stored in lfz_tmp_Rloc
                  !--        + <z d_z (jxB)_p> stored in lfpdz_tmp_Rloc and everything need * \beta
                     if ( l_leibniz ) then
                        lf_Rloc(n_m,n_r)=lf_Rloc(n_m,n_r)+DyMagFac*beta(n_r)*r(n_r)*( -lfpdz_tmp_Rloc(n_m,n_r) &
      !lf_Rloc(n_m,n_r)=lf_Rloc(n_m,n_r)+DyMagFac*beta(n_r)*r(n_r)*( lfp_tmp_Rloc(n_m,n_r)-lfp_h_tmp_Rloc(n_m,n_r) &!!
                        &                                                + ci*real(m,cp)*lfz_tmp_Rloc(n_m,n_r) )
                        !if ( n_r == 1 ) then !--> "unshifting beta" for a better Leibniz's rule accuracy
                        !    lf_Rloc(n_m,n_r)=lf_Rloc(n_m,n_r)+DyMagFac*                                                  &
                        !    &                  r(n_r)/(r_cmb**2-r(n_r)**2+epsilon(1.0_cp)*(r(n_r)-r_icb)/r_icb)*r(n_r)*( &
                        !    &                lfp_tmp_Rloc(n_m,n_r)-lfp_h_tmp_Rloc(n_m,n_r) + ci*real(m,cp)*lfz_tmp_Rloc(n_m,n_r) )
                        !    !&               -lfpdz_tmp_Rloc(n_m,n_r) + ci*real(m,cp)*lfz_tmp_Rloc(n_m,n_r) )
                        !else
                        !    lf_Rloc(n_m,n_r)=lf_Rloc(n_m,n_r)+DyMagFac*r(n_r)/(r_cmb**2-r(n_r)**2)*r(n_r)*( &
                        !    &                lfp_tmp_Rloc(n_m,n_r)-lfp_h_tmp_Rloc(n_m,n_r) + ci*real(m,cp)*lfz_tmp_Rloc(n_m,n_r) )
                        !    !&               -lfpdz_tmp_Rloc(n_m,n_r) + ci*real(m,cp)*lfz_tmp_Rloc(n_m,n_r) )
                        !end if
                     else
                        lf_Rloc(n_m,n_r)=lf_Rloc(n_m,n_r)+DyMagFac*beta(n_r)*( &
                        &                          -lfpdz_tmp_Rloc(n_m,n_r) +  &
                        &               ci*real(m,cp)*lfz_tmp_Rloc(n_m,n_r)   )
                        !!-- if LF benchmark on m<>0
                        !dpsidt_Rloc(n_m,n_r)=dpsidt_Rloc(n_m,n_r) - DyMagFac* &
                        !&                   lfpdz_tmp_Rloc(n_m,n_r)*beta(n_r) &
                        !& + DyMagFac*ci*real(m,cp)*beta(n_r)*lfz_tmp_Rloc(n_m,n_r)
                     end if
                  end if
                  !if ( n_r==1 .or. n_r==n_r_max ) lf_Rloc(n_m,n_r) = zero !--> will be done after adding the d_s jxB_p part in all update_psi_* files
               end if
            end if
         end do
      end do
      !print*, l_mag_LF, l_lin_solve
      !lf_Rloc(:,:)= (4.0_cp,1.0_cp)

#ifdef TOTO
      block
         use truncation, only: n_phi_max!, n_r_max!, n_m_max
         use fourier, only: ifft
         !use communications, only: allgather_from_Rloc

         integer :: n_s, file_handle!, n_r, n_m
         !real(cp) :: dpsidt_hat(1,n_r_max)
         real(cp) :: dpsidt_grid(n_phi_max,nRstart:nRstop)

         !call allgather_from_Rloc(real(dpsidt_Rloc(1,:)), dpsidt_hat, n_m_max)

         !!-- if LF benchmark --> full LF (without leibniz) still stored in dpsidt_Rloc
         do n_s=nRstart,nRstop
         !-- Bring data on the grid
            call ifft(dpsidt_Rloc(:,n_s), dpsidt_grid(:,n_s))
         end do

      if ( rank == 0 ) then

         open(newunit=file_handle, file='zavgVxjxBez_Rloc.dat', status='new', form='formatted')
         do n_s=1,n_r_max
            write(file_handle, *) r(n_s), dpsidt_grid(:,n_s)
         end do
         close(file_handle)

         open(newunit=file_handle, file='Leibniz_extra.dat', status='new', form='formatted')
         do n_s=1,n_r_max
            !do n_m=1,n_m_max
            !   print*, lf_Rloc(n_m,n_s)
            !end do
            write(file_handle, *) r(n_s), real(lf_Rloc(:,n_s)), aimag(lf_Rloc(:,n_s))
         end do
         close(file_handle)

         !!-- if LF benchmark --> full LF (without leibniz) still stored in dpsidt_Rloc
         open(newunit=file_handle, file='dpsidt_Rloc.dat', status='new', form='formatted')
         !print*, n_m_max, size(dpsidt_Rloc), size(real(dpsidt_Rloc(:,4))), size(aimag(dpsidt_Rloc(:,4)))
         !print*, dpsidt_Rloc(:2,4), dpsidt_Rloc(65,4)
         !print*, real(dpsidt_Rloc(:2,4)), aimag(dpsidt_Rloc(:2,4))
         do n_s=1,n_r_max
            write(file_handle, *) r(n_s), real(dpsidt_Rloc(:,n_s)), aimag(dpsidt_Rloc(:,n_s))
         end do
         close(file_handle)

         open(newunit=file_handle, file='dslfpds_tmp_rloc.dat', status='new', form='formatted')
         do n_s=1,n_r_max
            write(file_handle, *) r(n_s), real(lfpds_tmp_rloc(:,n_s)), aimag(lfpds_tmp_rloc(:,n_s))
         end do
         close(file_handle)

         open(newunit=file_handle, file='dzlfpdz_tmp_rloc.dat', status='new', form='formatted')
         do n_s=1,n_r_max
            write(file_handle, *) r(n_s), real(lfpdz_tmp_rloc(:,n_s)), aimag(lfpdz_tmp_rloc(:,n_s))
         end do
         close(file_handle)

         open(newunit=file_handle, file='lfz_tmp_rloc.dat', status='new', form='formatted')
         do n_s=1,n_r_max
            write(file_handle, *) r(n_s), real(lfz_tmp_rloc(:,n_s)), aimag(lfz_tmp_rloc(:,n_s))
         end do
         close(file_handle)

         open(newunit=file_handle, file='lfp_tmp_rloc.dat', status='new', form='formatted')
         do n_s=1,n_r_max
            write(file_handle, *) r(n_s), real(lfp_tmp_rloc(:,n_s)), aimag(lfp_tmp_rloc(:,n_s))
         end do
         close(file_handle)

         open(newunit=file_handle, file='lfs_tmp_rloc.dat', status='new', form='formatted')
         do n_s=1,n_r_max
            write(file_handle, *) r(n_s), real(lfs_tmp_rloc(:,n_s)), aimag(lfs_tmp_rloc(:,n_s))
         end do
         close(file_handle)
      end if
      end block

      print*, 'ALL GOOD RAD_3D!**'
      !stop
#endif

   end subroutine radial_loop_3D
!-------------------------------------------------------------------------------
   subroutine transform_to_grid_space(temp, B, dB, ddB, aj, dj, n_r, gsa)

      complex(cp), intent(inout) :: temp(lm_max)
      complex(cp), intent(inout) :: B(lm_max), dB(lm_max), ddB(lm_max)
      complex(cp), intent(inout) :: aj(lm_max), dj(lm_max)
      integer,     intent(in) :: n_r
      type(grid_space_arrays_t) :: gsa

#ifdef WITH_SHTNS
      if ( l_heat_3D ) then
         call scal_to_spat(temp, gsa%Tc)
      end if

      if ( l_mag_3D .and. (.not. l_lin_solve) ) then
         call torpol_to_spat(B, dB, aj, n_r, gsa%Brc, gsa%Btc, gsa%Bpc)
         if (  l_mag_LF ) then
            call torpol_to_curl_spat(B, ddB, aj, dj, n_r, gsa%curlBrc, &
                 &                   gsa%curlBtc, gsa%curlBpc)
         end if
      end if
#endif

   end subroutine transform_to_grid_space
!-------------------------------------------------------------------------------
   subroutine transform_to_lm_space(gsa, nl_lm, dVrTLMR)

      complex(cp), intent(out) :: dVrTLMR(:)
      type(grid_space_arrays_t) :: gsa
      type(nonlinear_lm_t) :: nl_lm

#ifdef WITH_SHTNS
      if ( l_heat_3D ) then
         call spat_to_qst(gsa%VTr, gsa%VTt, gsa%VTp, dVrTLMR, nl_lm%VTtLM, &
              &           nl_lm%VTpLM)
      end if

      if ( l_mag_3D .and. (.not. l_lin_solve) ) then
         call scal_to_SH(gsa%VxBr, nl_lm%VxBrLM)
         call scal_to_SH(gsa%VxBt, nl_lm%VxBtLM)
         call scal_to_SH(gsa%VxBp, nl_lm%VxBpLM)
         !-- In get_td --> ifdef SHTNS could be used instead and replace the 3 above lines by:
         !-- spat_to_qst(gsa%VxBr,gsa%VxBt, gsa%VxBp, nl_lm%VxBrLM, nl_lm%VxBtLM, nl_lm%VxBpLM)
      end if
#endif

   end subroutine transform_to_lm_space
!-------------------------------------------------------------------------------
   subroutine transform_axi_to_dth_grid_space(work_l, grad_thc)

      !-- Input variable
      complex(cp), intent(inout) :: work_l(:)

      !-- Output variable
      real(cp), intent(out) :: grad_thc(n_theta_max)

#ifdef WITH_SHTNS
      call scal_axi_to_grad_spat(work_l, grad_thc)
#endif

   end subroutine transform_axi_to_dth_grid_space
!------------------------------------------------------------------------------
   subroutine output_local_B_cmb(time, Br_CMB, Bp_CMB)

      !-- Input variable
      real(cp), intent(in) :: time
      real(cp), intent(in) :: Br_CMB(n_phi_max_3D,n_theta_max)
      real(cp), intent(in) :: Bp_CMB(n_phi_max_3D,n_theta_max)

      !-- Local variable
      integer :: n_phi
      real(cp) :: Br_CMB_phiAvg(n_theta_max), Bp_CMB_phiAvg(n_theta_max)
      !character(len=10) :: length_theta

      if ( rank == n_bavg_rank ) then
         print*, 'Bp_CMB(4,4) =', Bp_CMB(4, 4)
         Br_CMB_phiAvg(:) = 0.0_cp
         Bp_CMB_phiAvg(:) = 0.0_cp
         do n_phi=1, n_phi_max_3D
            Br_CMB_phiAvg(:) = Br_CMB_phiAvg(:) + Br_CMB(n_phi,:)/n_phi_max_3D
            Bp_CMB_phiAvg(:) = Bp_CMB_phiAvg(:) + Bp_CMB(n_phi,:)/n_phi_max_3D
         end do

         !write(length_theta, *) 2*n_theta_max
         !write(n_bavg_file, '(1P, ES20.12, '//length_theta//'ES16.8)') time, Br_CMB_phiAvg, Bp_CMB_phiAvg
         write(n_bavg_file, *) time, Br_CMB_phiAvg, Bp_CMB_phiAvg
      end if

   end subroutine output_local_B_cmb
!------------------------------------------------------------------------------
end module rloop_3D
