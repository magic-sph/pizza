module rloop_3D

   use precision_mod
   use parallel_mod
   use constants, only: ci
   use nonlinear_lm_mod, only: nonlinear_lm_t
   use grid_space_arrays_mod, only: grid_space_arrays_t
   use blocking, only: nRstart3D, nRstop3D, nRstart, nRstop
   use namelists, only: BuoFac
   use truncation_3D, only: lm_max, lmP_max, n_phi_max_3D, n_theta_max
   use truncation, only: idx2m, n_m_max
#ifdef WITH_SHTNS
   use shtns, only: spat_to_SH, scal_to_spat
#endif
   use outputs_3D, only: write_snaps
   use z_functions, only: zfunc_type
   use timers_mod, only: timers_type
   use time_schemes, only: type_tscheme

   implicit none

   private

   type(grid_space_arrays_t) :: gsa
   type(nonlinear_lm_t) :: nl_lm

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
   subroutine radial_loop_3D( time, ur, ut, up, temp, dtempdt, dVrTLM, dpsidt_Rloc, &
              &               l_frame, zinterp, timers, tscheme)

      !-- Input variables
      complex(cp), intent(in) :: temp(lm_max, nRstart3D:nRstop3D)
      real(cp),    intent(in) :: ur(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp),    intent(in) :: ut(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp),    intent(in) :: up(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      logical,     intent(in) :: l_frame
      real(cp),    intent(in) :: time
      type(zfunc_type),    intent(in) :: zinterp
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp), intent(out) :: dtempdt(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(out) :: dVrTLM(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(inout) :: dpsidt_Rloc(n_m_max, nRstart:nRstop)
      type(timers_type), intent(inout) :: timers


      !-- Local variables
      real(cp) :: buo_tmp(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)
      real(cp) :: runStart, runStop
      complex(cp) :: buo_tmp_Rloc(n_m_max,nRstart:nRstop)
      integer :: n_r, n_m, m

      runStart = MPI_Wtime()
      do n_r=nRstart3D,nRstop3D

         !-- Transform temperature from (l,m) to (theta,phi)
         call transform_to_grid_space_shtns(temp(:,n_r), gsa)

         !-- Construct non-linear terms in physical space
         call gsa%get_nl(ur(:,:,n_r), ut(:,:,n_r), up(:,:,n_r), n_r, &
              &          buo_tmp(:,:,n_r))

         !-- Transform back the non-linear terms to (l,m) space
         call transform_to_lm_space_shtns(gsa, nl_lm)

         !-- Get theta and phi derivatives using recurrence relations
         call nl_lm%get_td(dVrTLM(:,n_r), dtempdt(:,n_r))

      end do
      runStop = MPI_Wtime()
      if (runStop>runStart) then
         timers%n_r_loops_3D=timers%n_r_loops_3D+1
         timers%r_loop_3D   =timers%r_loop_3D+(runStop-runStart)
      end if

      !-- Write the 3-D snapshots in physical space
      if ( tscheme%istage ==1 .and. l_frame ) then
         call write_snaps(time, ur, ut, up)
      end if

      !-- Compute z-averaging of buoyancy
      runStart = MPI_Wtime()
      call zinterp%compute_avg(buo_tmp, buo_tmp_Rloc)
      runStop = MPI_Wtime()
      if (runStop>runStart) then
         timers%interp = timers%interp+(runStop-runStart)
      end if

      !-- Finish assembling buoyancy and sum it with dpsidt
      do n_m=1,n_m_max
         m = idx2m(n_m)
         do n_r=nRstart,nRstop
            dpsidt_Rloc(n_m,n_r)=dpsidt_Rloc(n_m,n_r)-BuoFac*ci*m* &
            &                    buo_tmp_Rloc(n_m,n_r)
         end do
      end do

   end subroutine radial_loop_3D
!-------------------------------------------------------------------------------
   subroutine transform_to_grid_space_shtns(temp, gsa)

      complex(cp), intent(in) :: temp(lm_max)
      type(grid_space_arrays_t) :: gsa

#ifdef WITH_SHTNS
      call scal_to_spat(temp, gsa%Tc)
#endif

   end subroutine transform_to_grid_space_shtns
!-------------------------------------------------------------------------------
   subroutine transform_to_lm_space_shtns(gsa, nl_lm)

      type(grid_space_arrays_t) :: gsa
      type(nonlinear_lm_t) :: nl_lm

#ifdef WITH_SHTNS
      call shtns_load_cfg(1)

      call spat_to_SH(gsa%VTr, nl_lm%VTrLM)
      call spat_to_SH(gsa%VTt, nl_lm%VTtLM)
      call spat_to_SH(gsa%VTp, nl_lm%VTpLM)

      call shtns_load_cfg(0)
#endif

   end subroutine transform_to_lm_space_shtns
!------------------------------------------------------------------------------
end module rloop_3D
