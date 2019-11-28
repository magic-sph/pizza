module rloop_3D

   use precision_mod
   use parallel_mod
   use constants, only: pi, ci, two
   use nonlinear_lm_mod, only: nonlinear_lm_t
   use grid_space_arrays_mod, only: grid_space_arrays_t
   use blocking, only: nRstart3D, nRstop3D, nRstart, nRstop
   use namelists, only: BuoFac, DyMagFac, tag, l_mag_3D, l_mag_LF, r_cmb, r_icb
   use truncation_3D, only: lm_max, lmP_max, n_phi_max_3D, n_theta_max, l_max, n_r_max_3D
   use truncation, only: idx2m, n_m_max
#ifdef WITH_SHTNS
   use shtns, only: spat_to_SH, scal_to_spat, torpol_to_spat, &
       &            torpol_to_curl_spat, scal_axi_to_grad_spat
#endif
   use horizontal, only: cost, sint
   use radial_functions, only: r, r_3D, or1
   use z_functions, only: zfunc_type
   use timers_mod, only: timers_type
   use time_schemes, only: type_tscheme
   use output_frames, only: open_snapshot_3D, close_snapshot_3D, &
       &                    write_bulk_snapshot_3D

   implicit none

   private

   type(grid_space_arrays_t) :: gsa
   type(nonlinear_lm_t) :: nl_lm
   integer :: frame_counter=1

   public :: radial_loop_3D, initialize_radial_loop_3D, finalize_radial_loop_3D

contains

   subroutine initialize_radial_loop_3D(lmP_max)

      integer, intent(in) :: lmP_max

      call gsa%initialize()
      call nl_lm%initialize(lmP_max)

   end subroutine initialize_radial_loop_3D
!------------------------------------------------------------------------------
   subroutine finalize_radial_loop_3D

      call nl_lm%finalize()
      call gsa%finalize()

   end subroutine finalize_radial_loop_3D
!------------------------------------------------------------------------------
   subroutine radial_loop_3D( time, ur, ut, up, temp, dtempdt, dVrTLM, &
              &               b_3D, db_3D, ddb_3D, aj_3D, dj_3D,       &
              &               dbdt_3D, djdt_3D, dVxBhLM,               &
              &               dpsidt_Rloc, dVsOm_Rloc,                 &
              &               l_frame, zinterp, timers, tscheme )

      !-- Input variables
      complex(cp), intent(in) :: temp(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(in) :: b_3D(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(in) :: db_3D(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(in) :: ddb_3D(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(in) :: aj_3D(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(in) :: dj_3D(lm_max, nRstart3D:nRstop3D)
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
      complex(cp),       intent(inout) :: dpsidt_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),       intent(inout) :: dVsOm_Rloc(n_m_max, nRstart:nRstop)
      type(timers_type), intent(inout) :: timers

      !-- Local variables
      real(cp) :: buo_tmp(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp) :: dTdth(n_theta_max,nRstart3D:nRstop3D)
      real(cp) :: runStart, runStop, phi!, rsint
      complex(cp) :: buo_tmp_Rloc(n_m_max,nRstart:nRstop)
      complex(cp) :: lfs_tmp_Rloc(n_m_max,nRstart:nRstop)
      complex(cp) :: lfp_tmp_Rloc(n_m_max,nRstart:nRstop)
      integer :: n_r, n_m, m, fh_temp, info_temp, n_theta, n_phi
      integer :: fh_ur, info_ur, fh_ut, info_ut, fh_up, info_up
      integer :: fh_br, info_br, fh_bt, info_bt, fh_bp, info_bp
      character(len=144) :: frame_name

#ifdef aDEBUG
      block
      use blocking_lm, only: st_map

      print*, st_map%lm2(0,0)
      print*, st_map%lm2(l_max,0)

      end block
#endif

      !-- get thermal wind
      do n_r=nRstart3D,nRstop3D
         !-- Take the gradient of the axisymmetric part
         call transform_axi_to_dth_grid_space(temp(1:l_max+1,n_r), dTdth(:,n_r))
      end do

      !-- Compute thermal wind --> modify up_3D_Rloc
      call zinterp%compute_thermal_wind(dTdth, up)

      !-- Open the 3-D snapshots in physical space
      if ( l_frame .and. tscheme%istage==1 ) then
         write(frame_name, '(A,I0,A,A)') 'frame_ur_3D_',frame_counter,'.',tag
         call open_snapshot_3D(frame_name, time, fh_ur, info_ur)
         write(frame_name, '(A,I0,A,A)') 'frame_ut_3D_',frame_counter,'.',tag
         call open_snapshot_3D(frame_name, time, fh_ut, info_ut)
         write(frame_name, '(A,I0,A,A)') 'frame_up_3D_',frame_counter,'.',tag
         call open_snapshot_3D(frame_name, time, fh_up, info_up)
         write(frame_name, '(A,I0,A,A)') 'frame_temp_3D_',frame_counter,'.',tag
         call open_snapshot_3D(frame_name, time, fh_temp, info_temp)
         write(frame_name, '(A,I0,A,A)') 'frame_br_3D_',frame_counter,'.',tag
         call open_snapshot_3D(frame_name, time, fh_br, info_br)
         write(frame_name, '(A,I0,A,A)') 'frame_bt_3D_',frame_counter,'.',tag
         call open_snapshot_3D(frame_name, time, fh_bt, info_bt)
         write(frame_name, '(A,I0,A,A)') 'frame_bp_3D_',frame_counter,'.',tag
         call open_snapshot_3D(frame_name, time, fh_bp, info_bp)
      end if

      ur(:,:,:) = 0.0_cp
      ut(:,:,:) = 0.0_cp
      up(:,:,:) = 0.0_cp
      runStart = MPI_Wtime()
      do n_r=nRstart3D,nRstop3D
         if ( n_r .ne. 1 .or. n_r .ne. n_r_max_3D ) then
            do n_theta=1,n_theta_max
               !rsint = r_3D(n_r)*sint(n_theta)
               do n_phi=1,n_phi_max_3D
                  phi = (n_phi-1)*two*pi/(n_phi_max_3D)
                  ur(n_phi,n_theta,n_r) = 0.0_cp
                  ut(n_phi,n_theta,n_r) = sint(n_theta)
                  up(n_phi,n_theta,n_r) = cos(4*phi)
               end do
            end do
         end if

         !-- Transform temperature and magnetic-field from (l,m) to (theta,phi)
         call transform_to_grid_space(temp(:,n_r), b_3D(:,n_r),      &
              &                      db_3D(:,n_r), ddb_3D(:,n_r),    &
              &                      aj_3D(:,n_r), dj_3D(:,n_r), n_r, gsa)

         !-- Write the snapshot on the grid (easier to handle)
         if ( l_frame .and. tscheme%istage==1 ) then
            call write_bulk_snapshot_3D(fh_ur, ur(:,:,n_r))
            call write_bulk_snapshot_3D(fh_ut, ut(:,:,n_r))
            call write_bulk_snapshot_3D(fh_up, up(:,:,n_r))
            call write_bulk_snapshot_3D(fh_temp, gsa%Tc(:,:))
            call write_bulk_snapshot_3D(fh_br, gsa%Brc(:,:))
            call write_bulk_snapshot_3D(fh_bt, gsa%Btc(:,:))
            call write_bulk_snapshot_3D(fh_bp, gsa%Bpc(:,:))
         end if

         !-- Construct non-linear terms in physical space
         call gsa%get_nl(ur(:,:,n_r), ut(:,:,n_r), up(:,:,n_r), n_r, &
              &          buo_tmp(:,:,n_r))

         !-- Transform back the non-linear terms to (l,m) space
         call transform_to_lm_space(gsa, nl_lm)

         !-- Get theta and phi derivatives using recurrence relations
         call nl_lm%get_td(dVrTLM(:,n_r), dtempdt(:,n_r),                     &
              &            dVxBhLM(:,n_r),dbdt_3D(:,n_r),djdt_3D(:,n_r), n_r )
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
         call close_snapshot_3D(fh_temp, info_temp)
         call close_snapshot_3D(fh_br, info_br)
         call close_snapshot_3D(fh_bt, info_bt)
         call close_snapshot_3D(fh_bp, info_bp)
         frame_counter = frame_counter+1
      end if

      !-- Compute z-averaging of buoyancy and lorentz-force
      runStart = MPI_Wtime()
      call zinterp%compute_avg(buo_tmp, buo_tmp_Rloc)
      if (  l_mag_LF ) call zinterp%compute_avg(gsa%jxBs, lfs_tmp_Rloc)
      if (  l_mag_LF ) call zinterp%compute_avg(gsa%jxBp, lfp_tmp_Rloc)
      runStop = MPI_Wtime()
      if (runStop>runStart) then
         timers%interp = timers%interp+(runStop-runStart)
      end if

#ifdef DEBUG
      block

         use truncation, only: n_r_max
         use truncation_3D, only: n_r_max_3D
         use radial_functions, only: r, r_3D

         integer :: n_r, file_handle

         open(newunit=file_handle, file='thw', status='new')
         do n_r=1,n_r_max_3D
            write(file_handle, '(5ES20.12)') r_3D(n_r), work(5,1,n_r),work(5,4,n_r),work(5,8,n_r),work(5,16,n_r)
         end do
         close(file_handle)

         open(newunit=file_handle, file='up', status='new')
         do n_r=1,n_r_max_3D
            write(file_handle, '(5ES20.12)') r_3D(n_r), up(5,1,n_r),up(5,4,n_r),up(5,8,n_r),up(5,16,n_r)
         end do
         close(file_handle)

         open(newunit=file_handle, file='buo', status='new')
         do n_r=1,n_r_max
            write(file_handle, '(4ES20.12)') r(n_r),  real(buo_tmp_Rloc(5,n_r)), real(buo_tmp_Rloc(6,n_r)), &
           &               real(buo_tmp_Rloc(1,n_r))
         end do
         close(file_handle)

         !open(newunit=file_handle, file='buo_3D', status='new', form='unformatted',&
         !&    access='stream')
         !write(file_handle) buo_tmp
         !close(file_handle)

      end block
#endif

      !-- Finish assembling buoyancy and sum it with dpsidt
      do n_m=1,n_m_max
         m = idx2m(n_m)
         do n_r=nRstart,nRstop
            dpsidt_Rloc(n_m,n_r)= dpsidt_Rloc(n_m,n_r)-BuoFac*ci*m* &
            &                    buo_tmp_Rloc(n_m,n_r)
            if ( l_mag_LF ) then
            !-- Finish assembling both parts of lorentz force and sum it with dpsidt and dVsOm
               if ( m == 0 ) then
                  dpsidt_Rloc(n_m,n_r)=dpsidt_Rloc(n_m,n_r)+DyMagFac*lfp_tmp_Rloc(n_m,n_r)
               else
                  dVsOm_Rloc(n_m,n_r)= dVsOm_Rloc(n_m,n_r)-DyMagFac*r(n_r)*&
                  &                  lfp_tmp_Rloc(n_m,n_r)
                  dpsidt_Rloc(n_m,n_r)= dpsidt_Rloc(n_m,n_r)-DyMagFac*ci*m* &
                  &                    lfs_tmp_Rloc(n_m,n_r)*or1(n_r)
               end if
            end if
         end do
      end do

   end subroutine radial_loop_3D
!-------------------------------------------------------------------------------
   subroutine transform_to_grid_space(temp, B, dB, ddB, aj, dj, n_r, gsa)

      complex(cp), intent(in) :: temp(lm_max)
      complex(cp), intent(in) :: B(lm_max), dB(lm_max), ddB(lm_max)
      complex(cp), intent(in) :: aj(lm_max), dj(lm_max)
      integer, intent(in) :: n_r
      type(grid_space_arrays_t) :: gsa

#ifdef WITH_SHTNS
      call scal_to_spat(temp, gsa%Tc)

      if ( l_mag_3D ) then
         call torpol_to_spat(B, dB, aj, n_r, gsa%Brc, gsa%Btc, gsa%Bpc)
         if (  l_mag_LF ) then
            call torpol_to_curl_spat(B, ddB, aj, dj, n_r, gsa%curlBrc, &
                 &                   gsa%curlBtc, gsa%curlBpc)
         end if
      end if
#endif

   end subroutine transform_to_grid_space
!-------------------------------------------------------------------------------
   subroutine transform_to_lm_space(gsa, nl_lm)

      type(grid_space_arrays_t) :: gsa
      type(nonlinear_lm_t) :: nl_lm

#ifdef WITH_SHTNS
      call shtns_load_cfg(1)

      call spat_to_SH(gsa%VTr, nl_lm%VTrLM)
      call spat_to_SH(gsa%VTt, nl_lm%VTtLM)
      call spat_to_SH(gsa%VTp, nl_lm%VTpLM)

      if ( l_mag_3D ) then
         call spat_to_SH(gsa%VxBr, nl_lm%VxBrLM)
         call spat_to_SH(gsa%VxBt, nl_lm%VxBtLM)
         call spat_to_SH(gsa%VxBp, nl_lm%VxBpLM)
         !-- In get_td --> ifdef SHTNS could be used instead and replace the 3 above lines by:
         !-- spat_to_qst(gsa%VxBr,gsa%VxBt, gsa%VxBp, nl_lm%VxBrLM, nl_lm%VxBtLM, nl_lm%VxBpLM)
         !-- can also be done for the temperature!
      end if

      call shtns_load_cfg(0)
#endif

   end subroutine transform_to_lm_space
!-------------------------------------------------------------------------------
   subroutine transform_axi_to_dth_grid_space(work_l, grad_thc)

      !-- Input variable
      complex(cp), intent(in) :: work_l(:)

      !-- Output variable
      real(cp), intent(out) :: grad_thc(n_theta_max)

#ifdef WITH_SHTNS
      call scal_axi_to_grad_spat(work_l, grad_thc)
#endif

   end subroutine transform_axi_to_dth_grid_space
!------------------------------------------------------------------------------
end module rloop_3D
