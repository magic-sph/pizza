module step_time
   !
   ! This module controls the time advance of the code
   !

   use iso_fortran_env, only: output_unit
   use communications, only: m2r_fields, r2m_fields, r2m_single, m2r_single
   use fields, only: us_Mloc, us_Rloc, up_Mloc, up_Rloc, temp_Mloc,         &
       &             temp_Rloc, om_Rloc, om_Mloc, psi_Mloc, dtemp_Mloc,     &
       &             dom_Mloc, temp_hat_Mloc, psi_hat_Mloc, xi_Mloc,        &
       &             xi_Rloc, dxi_Mloc, xi_hat_Mloc, fields_Mloc_container, &
       &             fields_Rloc_container, dtemp_Rloc, dxi_Rloc
   use fieldsLast, only: dpsidt_Rloc, dtempdt_Rloc, dVsT_Rloc, dVsT_Mloc, &
       &                 dVsOm_Rloc, dVsOm_Mloc, dpsidt, dTdt, dxidt,     &
       &                 dVsXi_Mloc, dVsXi_Rloc, dxidt_Rloc,              &
       &                 dt_fields_Rloc_container, dt_fields_Mloc_container
   use truncation, only: n_m_max, idx2m
   use courant_mod, only: dt_courant
   use radial_der, only: exch_ghosts, bulk_to_ghost
   use blocking, only: nRstart, nRstop
   use constants, only: half, one, ci
   use radial_functions, only: r
   use update_temp_coll, only: get_temp_rhs_imp_coll
   use update_xi_coll, only: get_xi_rhs_imp_coll
   use update_temp_fd_mod, only: get_temp_rhs_imp_ghost, temp_ghost, fill_ghosts_temp
   use update_temp_integ, only: get_temp_rhs_imp_int
   use update_xi_integ, only: get_xi_rhs_imp_int
   use update_xi_fd_mod, only: get_xi_rhs_imp_ghost, xi_ghost, fill_ghosts_xi
   use update_psi_fd_mod, only: get_psi_rhs_imp_ghost, psi_ghost, up0_ghost, &
       &                        fill_ghosts_psi
   use update_psi_integ_smat, only: get_psi_rhs_imp_int_smat
   use update_psi_integ_dmat, only: get_psi_rhs_imp_int_dmat
   use update_psi_coll_dmat, only: get_psi_rhs_imp_coll_dmat
   use update_psi_coll_smat, only: get_psi_rhs_imp_coll_smat
   use mloop_mod, only: mloop, finish_explicit_assembly, assemble_stage
   use mloop_fd_mod, only: finish_explicit_assembly_Rdist, mloop_Rdist, &
       &                   assemble_stage_Rdist
   use rLoop, only: radial_loop
   use namelists, only: n_time_steps, alpha, dtMax, dtMin, l_bridge_step, &
       &                tEND, run_time_requested, n_log_step, n_frames,   &
       &                n_frame_step, n_checkpoints, n_checkpoint_step,   &
       &                n_spec_step, n_specs, l_vphi_balance, l_AB1,      &
       &                l_cheb_coll, l_direct_solve, l_vort_balance,      &
       &                l_heat, l_chem, l_packed_transp, l_finite_diff
   use outputs, only: n_log_file, write_outputs, vp_bal, vort_bal, &
       &              read_signal_file, spec
   use useful, only: logWrite, abortRun, formatTime, l_correct_step
   use time_schemes, only: type_tscheme
   use parallel_mod
   use precision_mod
   use timers_mod, only: timers_type

   implicit none

   private

   real(cp) :: tsig ! to measure the time between two signal files

   public :: time_loop

contains

   subroutine time_loop(time, tscheme, run_time_init)

      !-- Input variable
      real(cp),            intent(in) :: run_time_init

      !-- Output variables
      real(cp),            intent(inout) :: time
      class(type_tscheme), intent(inout) :: tscheme

      !-- Local variables
      integer :: n_time_step, n_time_steps_go, n_time_steps_run
      integer :: nPercent
      real(cp) :: tenth_n_time_steps
      character(len=255) :: message

      !-- Courant
      real(cp) :: dtr,dth,dt_new
      real(cp) :: dtr_Rloc(nRstart:nRstop), dth_Rloc(nRstart:nRstop)

      !-- Timings:
      type(timers_type) :: timers
      integer :: n_stage
      integer :: n_stop_signal, n_spec_signal, n_rst_signal, n_frame_signal
      integer :: signals(4)
      real(cp) :: run_time_passed
      real(cp) :: runStart, runStop, runStartT, runStopT

      logical :: l_new_dt, l_rst, l_frame, l_log, l_log_next
      logical :: l_vphi_bal_write, l_stop_time
      logical :: lMat, lMatNext

      tenth_n_time_steps=real(n_time_steps,kind=cp)/10.0_cp
      nPercent = 9

      l_new_dt        =.true.
      l_rst           =.false.
      l_frame         =.false.
      l_log           =.false.
      l_log_next      =.true.
      l_stop_time     =.false.
      l_vphi_bal_write=.false.
      lMatNext        =.true.

      tsig           =0.0_cp
      run_time_passed=0.0_cp

      n_stop_signal  = 0
      n_spec_signal  = 0
      n_rst_signal   = 0
      n_frame_signal = 0

      !-- Set the counters to zero
      call timers%initialize()

      !-- Dummy initial timings
      dtr_Rloc(:) = 1e10_cp
      dth_Rloc(:) = 1e10_cp

      !!!!! Time loop starts !!!!!!
      if ( n_time_steps == 1 ) then
         n_time_steps_run=1 ! Output only, for example G-file/movie etc.
      else if ( n_time_steps == 2 ) then
         n_time_steps_run=2 ! 
      else
         n_time_steps_run=n_time_steps+1  ! Last time step for output only !
      end if

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      n_time_steps_go = 0
      outer: do n_time_step=1,n_time_steps_run

         runStartT = MPI_Wtime()
         !-------------------
         !-- Determine whether we will need outputs at this time step
         !-------------------
         l_log = l_correct_step(n_time_step-1,n_time_steps,n_log_step,0)
         if ( n_time_step+1 <= n_time_steps+1 ) then
            l_log_next = l_correct_step(n_time_step,n_time_steps,n_log_step,0)
         end if
         l_rst = l_correct_step(n_time_step-1,n_time_steps,n_checkpoint_step, &
                 &              n_checkpoints) .or. n_rst_signal == 1
         l_frame = l_correct_step(n_time_step-1,n_time_steps,n_frame_step,n_frames) &
         &         .or. n_frame_signal == 1
         spec%l_calc = l_correct_step(n_time_step-1,n_time_steps,n_spec_step,n_specs) &
         &        .or. n_spec_signal == 1
         l_vphi_bal_write = (l_log .and. l_vphi_balance) .and.  (n_time_step>1)
         vp_bal%l_calc = l_log_next .and. l_vphi_balance
         vort_bal%l_calc = l_log_next .and. l_vort_balance

         !-----------------
         !-- Check SIGNALS
         !-----------------
         call check_signals(run_time_passed, signals)
         n_stop_signal  = signals(1)
         n_frame_signal = signals(2)
         n_rst_signal   = signals(3)
         n_spec_signal  = signals(4)

         !-------------------
         !-- Check whether the run is not getting out of time
         !-------------------
         call MPI_Allreduce(MPI_IN_PLACE,timers%tot,1,MPI_DEF_REAL, &
              &             MPI_MAX,MPI_COMM_WORLD,ierr)
         if ( timers%tot+run_time_init > run_time_requested ) then
            write(message,'("! Run time limit exeeded !")')
            call logWrite(message, n_log_file)
            l_stop_time=.true.
         end if

         !-- Some reasons to stop the run
         if ( n_stop_signal > 0 ) l_stop_time=.true.
         if ( n_time_step == n_time_steps_run ) l_stop_time=.true.
         if ( time > tEND .and. tEND /= 0.0_cp ) l_stop_time=.true.

         if ( n_time_step == 1 ) l_log=.true.

         if ( l_stop_time ) then             
            l_rst=.true.           
            l_log=.true.
         end if

         !-------------------
         !-- Outputs
         !-------------------
         !-- Get time series
         runStart = MPI_Wtime()
         if ( l_finite_diff .and. ( l_log .or. l_rst .or. l_frame ) ) then
            call transp_Rloc_to_Mloc_IO()
         end if
         call write_outputs(time, tscheme, n_time_step, l_log, l_rst, l_frame, &
              &             l_vphi_bal_write, l_stop_time, us_Mloc,  up_Mloc,  &
              &             om_Mloc, temp_Mloc, dtemp_Mloc, xi_Mloc, dxi_Mloc, &
              &             dpsidt, dTdt, dxidt)
         runStop = MPI_Wtime()
         if (runStop>runStart) then
            timers%n_io_calls=timers%n_io_calls+1
            timers%io        =timers%io+(runStop-runStart)
         end if

         !-- If this is running out of time, exit the loop and terminate the calculations
         if ( l_stop_time ) exit outer

         !-- Now loop over the stages
         tscheme%istage = 1
         do n_stage=1,tscheme%nstages

            if ( tscheme%l_exp_calc(n_stage) .or. vp_bal%l_calc .or. vort_bal%l_calc ) then
               !-------------------
               !-- MPI transpositions from m-distributed to r-distributed
               !-------------------
               call transp_Mloc_to_Rloc(timers%n_mpi_comms, timers%mpi_comms)

               !-------------------
               !-- Radial loop
               !-------------------
               runStart = MPI_Wtime()
               if ( l_finite_diff ) then
                  call radial_loop( us_Rloc, up_Rloc, om_Rloc, temp_Rloc, xi_Rloc, &
                       &            dTdt%expl(:,:,tscheme%istage), dVsT_Rloc,      &
                       &            dxidt%expl(:,:,tscheme%istage), dVsXi_Rloc,    &
                       &            dpsidt%expl(:,:,tscheme%istage), dVsOm_Rloc,   &
                       &            dtr_Rloc, dth_Rloc, timers, tscheme )
               else
                  call radial_loop( us_Rloc, up_Rloc, om_Rloc, temp_Rloc, xi_Rloc, &
                       &            dtempdt_Rloc, dVsT_Rloc, dxidt_Rloc,           &
                       &            dVsXi_Rloc, dpsidt_Rloc, dVsOm_Rloc, dtr_Rloc, &
                       &            dth_Rloc, timers, tscheme )
               end if
               runStop = MPI_Wtime()
               if (runStop>runStart) then
                  timers%n_r_loops=timers%n_r_loops+1
                  timers%r_loop   =timers%r_loop+(runStop-runStart)
               end if

               !------------------
               !-- MPI transpositions from r-distributed to m-distributed
               !------------------
               call transp_Rloc_to_Mloc(tscheme%istage, timers%mpi_comms)

               !--------------------
               !-- Finish assembling the explicit terms
               !--------------------
               runStart = MPI_Wtime()
               if ( l_finite_diff ) then
                  call finish_explicit_assembly_Rdist(us_Rloc,temp_Rloc,xi_Rloc,   &
                       &                              dVsT_Rloc,dVsXi_Rloc,        &
                       &                              dVsOm_Rloc,dTdt,dxidt,dpsidt,&
                       &                              tscheme,vp_bal,vort_bal)
               else
                  call finish_explicit_assembly(temp_Mloc, xi_Mloc, psi_Mloc,      &
                       &                        us_Mloc, up_Mloc, om_Mloc,         &
                       &                        dVsT_Mloc(:,:,tscheme%istage),     &
                       &                        dVsXi_Mloc(:,:,tscheme%istage),    &
                       &                        dVsOm_Mloc(:,:,tscheme%istage),    &
                       &                        dTdt, dxidt, dpsidt, tscheme,      &
                       &                        vp_bal, vort_bal)
               end if
               runStop = MPI_Wtime()
               if (runStop>runStart) then
                  timers%m_loop=timers%m_loop+(runStop-runStart)
               end if
            end if

            if ( tscheme%istage == 1 ) then
               !-------------------
               !------ Checking Courant criteria, l_new_dt and dt_new are output
               !-------------------
               call dt_courant(dtr,dth,l_new_dt,tscheme%dt(1),dt_new,dtMax, &
                    &          dtr_Rloc,dth_Rloc,time)

               call tscheme%set_dt_array(dt_new,dtMin,time,n_log_file,n_time_step,&
                    &                    l_new_dt)
               !-- Store the old weight factor of matrices
               !-- if it changes because of dt factors moving
               !-- matrix also needs to be rebuilt
               call tscheme%set_weights(lMatNext)

               !----- Advancing time:
               time=time+tscheme%dt(1) ! Update time
            end if

            lMat=lMatNext
            if ( (l_new_dt .or. lMat) .and. (tscheme%istage==1) ) then
               !----- Calculate matrices for new time step if dt /= dtLast
               lMat=.true.
               if ( rank == 0 ) then
                  write(*,'(1p,'' ! Building matrices at time step:'',   &
                       &              i8,ES16.6)') n_time_step,time
               end if
            end if
            lMatNext = .false.

            !-- If the scheme is a multi-step scheme that is not Crank-Nicolson 
            !-- we have to use a different starting scheme
            call start_from_another_scheme(l_bridge_step,n_time_step, tscheme)

            !--------------------
            !-- M-loop (update routines)
            !--------------------
            if ( (.not. tscheme%l_assembly) .or. (tscheme%istage/=tscheme%nstages) ) then
               runStart = MPI_Wtime()
               if ( l_finite_diff ) then
                  call mloop_Rdist(temp_Rloc, dtemp_Rloc, xi_Rloc, dxi_Rloc, &
                       &           om_Rloc, us_Rloc, up_Rloc, dTdt, dxidt,   &
                       &           dpsidt, vp_bal, vort_bal, tscheme, lMat,  &
                       &           timers)
               else
                  call mloop(temp_hat_Mloc, xi_hat_Mloc, temp_Mloc, dtemp_Mloc,       &
                       &     xi_Mloc, dxi_Mloc, psi_hat_Mloc, psi_Mloc, om_Mloc,      & 
                       &     dom_Mloc, us_Mloc, up_Mloc, dTdt, dxidt, dpsidt, vp_bal, &
                       &     vort_bal, tscheme, lMat, l_log_next, timers)
               end if
               runStop = MPI_Wtime()
               if ( .not. lMat ) then
                  if (runStop>runStart) then
                     timers%n_m_loops=timers%n_m_loops+1
                     timers%m_loop   =timers%m_loop+(runStop-runStart)
                  end if
               else
                  if (runStop>runStart) then
                     timers%n_m_loops_mat=timers%n_m_loops_mat+1
                     timers%m_loop_mat   =timers%m_loop_mat+(runStop-runStart)
                  end if
               end if

               ! Increment current stage
               tscheme%istage = tscheme%istage+1
            end if

         end do ! Finish loop over stages

         !----------------------------
         !-- Assembly stage of IMEX-RK (if needed)
         !----------------------------
         if ( tscheme%l_assembly ) then
            if ( l_finite_diff ) then
               call assemble_stage_Rdist(temp_Rloc, dtemp_Rloc, xi_Rloc, dxi_Rloc, &
                    &              us_Rloc, up_Rloc, om_Rloc, dTdt, dxidt, dpsidt, &
                    &              tscheme, vp_bal, vort_bal)
            else
               call assemble_stage(temp_hat_Mloc, xi_hat_Mloc, temp_Mloc, dtemp_Mloc, &
                    &              xi_Mloc, dxi_Mloc, psi_hat_Mloc, psi_Mloc,         &
                    &              us_Mloc, up_Mloc, om_Mloc, dTdt, dxidt, dpsidt,    &
                    &              tscheme, vp_bal, vort_bal, l_log_next, timers)
            end if
         end if

         !---------------------
         !-- Timings
         !---------------------
         runStopT = MPI_Wtime()
         if (runStopT>runStartT) then
            timers%tot=timers%tot+(runStopT-runStartT)
         end if
         n_time_steps_go = n_time_steps_go+1

         !---------------------
         !-- Info about run advance
         !---------------------
         run_time_passed=timers%tot
         run_time_passed=run_time_passed/n_time_steps_go
         if ( real(n_time_step,cp)+tenth_n_time_steps*real(nPercent,cp) >=  &
            & real(n_time_steps,cp)  .or. n_time_steps < 31 ) then
            write(message,'(" ! Time step finished:",i8)') n_time_step
            call logWrite(message, n_log_file)

            if ( real(n_time_step,cp)+tenth_n_time_steps*real(nPercent,cp) >= &
               & real(n_time_steps,cp) .and. n_time_steps >= 10 ) then
               write(message,'(" ! This is           :",i3,"%")') (10-nPercent)*10
               call logWrite(message, n_log_file)
               nPercent=nPercent-1
            end if
            if ( rank == 0 ) then
               call formatTime(output_unit,' ! Mean wall time for time step:',  &
               &               run_time_passed)
               call formatTime(n_log_file,' ! Mean wall time for time step:', &
               &               run_time_passed)
            end if

         end if

      end do outer ! end of time stepping !

      !-- Average timers over ranks and number of calls
      call timers%finalize(n_time_steps_go)

      !-- Write timers info in log file  and on display
      call timers%write_log(n_log_file)
      call timers%write_log(output_unit)

   end subroutine time_loop
!-------------------------------------------------------------------------------
   subroutine start_from_another_scheme(l_bridge_step, n_time_step, tscheme)

      !-- Input variables
      logical,             intent(in) :: l_bridge_step
      integer,             intent(in) :: n_time_step
      class(type_tscheme), intent(inout) :: tscheme

      !-- Local variables
      integer :: n_r, n_m, m
      complex(cp) :: psi_Rloc(n_m_max,nRstart:nRstop)

      !-- If the scheme is a multi-step scheme that is not Crank-Nicolson 
      !-- we have to use a different starting scheme
      if ( l_bridge_step .and. tscheme%time_scheme /= 'CNAB2' .and.  &
           n_time_step <= tscheme%nold-1 .and. tscheme%family=='MULTISTEP' ) then

         if ( l_finite_diff ) then
            if ( l_heat ) then
               call bulk_to_ghost(temp_Rloc, temp_ghost, 1, nRstart, nRstop, n_m_max,&
                    &             1, n_m_max)
               call exch_ghosts(temp_ghost, n_m_max, nRstart, nRstop, 1)
               call fill_ghosts_temp(temp_ghost)
               call get_temp_rhs_imp_ghost(temp_ghost, dtemp_Rloc, dTdt, 1, .true.)
            end if

            if ( l_chem ) then
               call bulk_to_ghost(xi_Rloc, xi_ghost, 1, nRstart, nRstop, n_m_max,&
                    &             1, n_m_max)
               call exch_ghosts(xi_ghost, n_m_max, nRstart, nRstop, 1)
               call fill_ghosts_xi(xi_ghost)
               call get_xi_rhs_imp_ghost(xi_ghost, dxi_Rloc, dxidt, 1, .true.)
            end if

            call bulk_to_ghost(up_Rloc(1,:), up0_ghost, 1, nRstart, nRstop, 1, 1, 1)
            call exch_ghosts(up0_ghost, 1, nRstart, nRstop, 1)
            !-- psi is unknown: rebuild it from us
            do n_r=nRstart,nRstop
               do n_m=1,n_m_max
                  m = idx2m(n_m)
                  if ( m == 0 ) cycle
                  psi_Rloc(n_m,n_r)=us_Rloc(n_m,n_r)/(ci*real(m,cp)) * r(n_r)
               end do
            end do
            call bulk_to_ghost(psi_Rloc, psi_ghost, 2, nRstart, nRstop, n_m_max, 1, &
                 &             n_m_max)
            call exch_ghosts(psi_ghost, 2, nRstart, nRstop, n_m_max)
            call fill_ghosts_psi(psi_ghost, up0_ghost)
            call get_psi_rhs_imp_ghost(psi_ghost, up0_ghost, us_Rloc, up_Rloc, om_Rloc, &
                 &                     dpsidt, 1, vp_bal, vort_bal, .true.)
         else
            if ( l_cheb_coll ) then

               if ( l_heat ) call get_temp_rhs_imp_coll(temp_Mloc,dtemp_Mloc, &
                                  &                     dTdt, 1, .true.)
               if ( l_chem ) call get_xi_rhs_imp_coll(xi_Mloc, dxi_Mloc, dxidt, 1, .true.)
               if ( l_direct_solve ) then
                  call get_psi_rhs_imp_coll_smat(us_Mloc, up_Mloc, om_Mloc,    &
                       &                         dom_Mloc, temp_Mloc, xi_Mloc, &
                       &                         dpsidt, 1, vp_bal, vort_bal,  &
                       &                         .true.)
               else
                  call get_psi_rhs_imp_coll_dmat(up_Mloc, om_Mloc, dom_Mloc,  &
                       &                         temp_Mloc, xi_Mloc, dpsidt,  &
                       &                         1, vp_bal, vort_bal, .true.)
               end if
            else
               if ( l_heat ) call get_temp_rhs_imp_int(temp_hat_Mloc, dTdt, 1, .true.)
               if ( l_chem ) call get_xi_rhs_imp_int(xi_hat_Mloc, dxidt, 1, .true.)
               if ( l_direct_solve ) then
                  call get_psi_rhs_imp_int_smat(psi_hat_Mloc,up_Mloc,temp_Mloc,psi_Mloc, &
                       &                        dpsidt, 1, vp_bal, .true.)
               else
                  call get_psi_rhs_imp_int_dmat(om_Mloc,up_Mloc,temp_Mloc,xi_Mloc,dpsidt,&
                       &                        1, vp_bal, .true.)
               end if
            end if
         end if

         call tscheme%bridge_with_cnab2()

      end if

      if ( l_AB1 .and. n_time_step == 1 ) then
         call tscheme%start_with_ab1()
         l_AB1 = .false.
      end if

   end subroutine start_from_another_scheme
!--------------------------------------------------------------------------------
   subroutine check_signals(run_time_passed, signals)

      !-- Input variable:
      real(cp), intent(in) :: run_time_passed

      !-- Output variables
      integer, intent(out) :: signals(4)

      !-- Local variables:
      logical :: l_check_signal

      tsig = tsig+run_time_passed
      if ( rank == 0 ) then
         if ( tsig > 2.0_cp ) then ! Only check signals every second
            l_check_signal = .true.
            tsig = 0.0_cp
         else
            l_check_signal = .false.
         end if
      end if

      call MPI_Bcast(l_check_signal,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

      if ( l_check_signal ) then
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         call read_signal_file(signals)
      else
         signals(:) = 0
      end if

   end subroutine check_signals
!--------------------------------------------------------------------------------
   subroutine transp_Mloc_to_Rloc(n_mpi_comms, mpi_comms)
      !
      !- MPI transposition from M-distributed to r-distributed
      !

      !-- Output variables
      integer,  intent(inout) :: n_mpi_comms
      real(cp), intent(inout) :: mpi_comms

      !-- Local variables
      real(cp) :: runStart, runStop

      runStart = MPI_Wtime()
      if ( .not. l_finite_diff ) then
         if ( l_packed_transp ) then
            call m2r_fields%transp_m2r(fields_Mloc_container, fields_Rloc_container)
         else
            call m2r_single%transp_m2r(us_Mloc, us_Rloc)
            call m2r_single%transp_m2r(up_Mloc, up_Rloc)
            call m2r_single%transp_m2r(om_Mloc, om_Rloc)
            if ( l_heat ) call m2r_single%transp_m2r(temp_Mloc, temp_Rloc)
            if ( l_chem ) call m2r_single%transp_m2r(xi_Mloc, xi_Rloc)
         end if
      end if
      runStop = MPI_Wtime()
      if (runStop>runStart) then
         n_mpi_comms=n_mpi_comms+1
         mpi_comms  =mpi_comms+runStop-runStart
      end if

   end subroutine transp_Mloc_to_Rloc
!--------------------------------------------------------------------------------
   subroutine transp_Rloc_to_Mloc(istage, mpi_comms)
      !
      !- MPI transposition from r-distributed to M-distributed
      !
      integer, intent(in) :: istage

      !-- Output variable
      real(cp), intent(inout) :: mpi_comms

      !-- Local variables
      real(cp) :: runStart, runStop

      runStart = MPI_Wtime()
      if ( .not. l_finite_diff ) then
         if ( l_packed_transp ) then
            call r2m_fields%transp_r2m(dt_fields_Rloc_container, &
                 &                     dt_fields_Mloc_container(:,:,:,istage))
         else
            if ( l_heat ) then
               call r2m_single%transp_r2m(dtempdt_Rloc,dTdt%expl(:,:,istage))
               call r2m_single%transp_r2m(dVsT_Rloc, dVsT_Mloc(:,:,istage))
            end if
            if ( l_chem ) then
               call r2m_single%transp_r2m(dxidt_Rloc,dxidt%expl(:,:,istage))
               call r2m_single%transp_r2m(dVsXi_Rloc, dVsXi_Mloc(:,:,istage))
            end if
            call r2m_single%transp_r2m(dpsidt_Rloc,dpsidt%expl(:,:,istage))
            call r2m_single%transp_r2m(dVsOm_Rloc, dVsOm_Mloc(:,:,istage))
         end if
      end if
      runStop = MPI_Wtime()
      if (runStop>runStart) then
         mpi_comms=mpi_comms+runStop-runStart
      end if

   end subroutine transp_Rloc_to_Mloc
!--------------------------------------------------------------------------------
   subroutine transp_Rloc_to_Mloc_IO

      call r2m_single%transp_r2m(us_Rloc, us_Mloc)
      call r2m_single%transp_r2m(up_Rloc, up_Mloc)
      call r2m_single%transp_r2m(om_Rloc, om_Mloc)
      if ( l_heat ) then
         call r2m_single%transp_r2m(temp_Rloc, temp_Mloc)
         call r2m_single%transp_r2m(dtemp_Rloc, dtemp_Mloc)
      end if

      if ( l_chem ) then
         call r2m_single%transp_r2m(xi_Rloc, xi_Mloc)
         call r2m_single%transp_r2m(dxi_Rloc, dxi_Mloc)
      end if

   end subroutine transp_Rloc_to_Mloc_IO
!--------------------------------------------------------------------------------
end module step_time
