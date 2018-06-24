module step_time

   use communications, only: transp_m2r, m2r_fields, transp_r2m, r2m_fields, &
       &                     gather_from_mloc_to_rank0, my_reduce_mean,      &
       &                     scatter_from_rank0_to_mloc
   use fields, only: us_Mloc, us_Rloc, up_Mloc, up_Rloc, temp_Mloc,     &
       &             temp_Rloc, om_Rloc, om_Mloc, psi_Mloc, dtemp_Mloc, &
       &             dom_Mloc
   use fieldsLast, only: dpsidt_Rloc, dtempdt_Rloc, dVsT_Rloc, dVsT_Mloc, &
       &                 dVsOm_Rloc, dVsOm_Mloc, buo_Mloc, dpsidt, dTdt
   use courant_mod, only: dt_courant
   use blocking, only: nRstart, nRstop
   use constants, only: half, one
   use update_temp_coll, only: update_temp_co, get_temp_rhs_imp_coll
   use update_temp_integ, only: update_temp_int, get_temp_rhs_imp_int
   use update_psi_integ_smat, only: get_psi_rhs_imp_int_smat
   use update_psi_integ_dmat, only: get_psi_rhs_imp_int_dmat
   use update_psi_coll_dmat, only: get_psi_rhs_imp_coll_dmat
   use update_psi_coll_smat, only: get_psi_rhs_imp_coll_smat
   use mloop_mod, only: mloop
   use rLoop, only: radial_loop
   use namelists, only: n_time_steps, alpha, dtMax, dtMin, l_bridge_step, &
       &                tEND, run_time_requested, n_log_step, n_frames,   &
       &                n_frame_step, n_checkpoints, n_checkpoint_step,   &
       &                n_spec_step, n_specs, l_vphi_balance, l_AB1,      &
       &                l_cheb_coll, l_direct_solve
   use outputs, only: n_log_file, write_outputs, vp_bal, terminate_vp_bal
   use useful, only: logWrite, abortRun, formatTime, l_correct_step
   use time_schemes, only: type_tscheme
   use parallel_mod
   use precision_mod

   implicit none

   private

   public :: time_loop

contains

   subroutine time_loop(time, tscheme)

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
      integer :: n_r_loops, n_mpi_comms, n_m_loops, n_m_loops_mat
      integer :: n_io_calls, n_fft_calls, n_solve_calls, n_dct_calls
      integer :: n_lu_calls, n_stage
      real(cp) :: run_time_r_loop, run_time_mpi_comms
      real(cp) :: run_time_m_loop, run_time_m_loop_mat, run_time_fft
      real(cp) :: run_time_tot, run_time_io, run_time_passed, run_time_solve
      real(cp) :: run_time_dct, run_time_lu
      real(cp) :: runStart, runStop, runStartT, runStopT

      logical :: l_new_dt
      logical :: l_rst
      logical :: l_frame
      logical :: l_spec
      logical :: l_log, l_log_next
      logical :: l_vphi_bal_calc, l_vphi_bal_write
      logical :: l_stop_time
      logical :: lMat, lMatNext

      tenth_n_time_steps=real(n_time_steps,kind=cp)/10.0_cp
      nPercent = 9

      l_new_dt        =.true.
      l_rst           =.false.
      l_frame         =.false.
      l_spec          =.false.
      l_log           =.false.
      l_log_next      =.true.
      l_stop_time     =.false.
      l_vphi_bal_calc =.false.
      l_vphi_bal_write=.false.
      lMatNext        =.true.

      !-- Dummy initial timings
      dtr_Rloc(:) = 1e10_cp
      dth_Rloc(:) = 1e10_cp

      n_r_loops     = 0
      n_m_loops     = 0
      n_m_loops_mat = 0
      n_mpi_comms   = 0
      n_io_calls    = 0
      n_fft_calls   = 0
      n_solve_calls = 0
      n_dct_calls   = 0
      n_lu_calls    = 0
      run_time_r_loop     = 0.0_cp
      run_time_m_loop     = 0.0_cp
      run_time_io         = 0.0_cp
      run_time_mpi_comms  = 0.0_cp
      run_time_m_loop_mat = 0.0_cp
      run_time_fft        = 0.0_cp
      run_time_lu         = 0.0_cp
      run_time_solve      = 0.0_cp
      run_time_dct        = 0.0_cp
      run_time_tot        = 0.0_cp

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

         !-------------------
         !-- Determine whether we will need outputs at this time step
         !-------------------
         l_log = l_correct_step(n_time_step-1,n_time_steps,n_log_step,0)
         if ( n_time_step+1 <= n_time_steps+1 ) then
            l_log_next = l_correct_step(n_time_step,n_time_steps,n_log_step,0)
         end if
         l_rst = l_correct_step(n_time_step-1,n_time_steps,n_checkpoint_step, &
                 &              n_checkpoints)
         l_frame = l_correct_step(n_time_step-1,n_time_steps,n_frame_step,n_frames)
         l_spec = l_correct_step(n_time_step-1,n_time_steps,n_spec_step,n_specs)
         l_vphi_bal_write = l_log .and. l_vphi_balance
         l_vphi_bal_calc = l_log_next .and. l_vphi_balance

         !-------------------
         !-- Check whether the run is not getting out of time
         !-------------------
         call MPI_Allreduce(MPI_IN_PLACE,run_time_tot,1,MPI_INTEGER8, &
              &             MPI_MAX,MPI_COMM_WORLD,ierr)
         if ( run_time_tot > run_time_requested ) then
            write(message,'("! Run time limit exeeded !")')
            call logWrite(message, n_log_file)
            l_stop_time=.true.
         end if
         !-- Some reasons to stop the run
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
         call write_outputs(time, tscheme, n_time_step, l_log, l_rst, l_spec,   &
              &             l_frame, l_vphi_bal_calc, l_vphi_bal_write,         &
              &             l_stop_time, us_Mloc,  up_Mloc, om_Mloc, temp_Mloc, &
              &             dtemp_Mloc, dpsidt, dTdt)
         runStop = MPI_Wtime()
         if (runStop>runStart) then
            n_io_calls  =n_io_calls+1
            run_time_io=run_time_io+(runStop-runStart)
         end if

         !-- If this is running out of time, exit the loop and terminate the calculations
         if ( l_stop_time ) exit outer

         !-- Now loop over the stages
         tscheme%istage = 1
         runStartT = MPI_Wtime()
         do n_stage=1,tscheme%nstages
            !-------------------
            !-- MPI transpositions from m-distributed to r-distributed
            !-------------------
            runStart = MPI_Wtime()
            call transp_m2r(m2r_fields, us_Mloc, us_Rloc)
            call transp_m2r(m2r_fields, up_Mloc, up_Rloc)
            call transp_m2r(m2r_fields, temp_Mloc, temp_Rloc)
            call transp_m2r(m2r_fields, om_Mloc, om_Rloc)
            runStop = MPI_Wtime()
            if (runStop>runStart) then
               n_mpi_comms  =n_mpi_comms+1
               run_time_mpi_comms=run_time_mpi_comms+(runStop-runStart)
            end if

            !-------------------
            !-- Radial loop
            !-------------------
            runStart = MPI_Wtime()
            call radial_loop( us_Rloc, up_Rloc, om_Rloc, temp_Rloc, dtempdt_Rloc, &
                 &            dVsT_Rloc, dpsidt_Rloc, dVsOm_Rloc, dtr_Rloc,       &
                 &            dth_Rloc, run_time_fft, n_fft_calls, tscheme )
            runStop = MPI_Wtime()
            if (runStop>runStart) then
               n_r_loops  =n_r_loops+1
               run_time_r_loop=run_time_r_loop+(runStop-runStart)
            end if

            !------------------
            !-- MPI transpositions from r-distributed to m-distributed
            !------------------
            runStart = MPI_Wtime()
            call transp_r2m(r2m_fields, dtempdt_Rloc, dTdt%expl(:,:,tscheme%istage))
            call transp_r2m(r2m_fields, dpsidt_Rloc, dpsidt%expl(:,:,tscheme%istage))
            call transp_r2m(r2m_fields, dVsT_Rloc, dVsT_Mloc)
            call transp_r2m(r2m_fields, dVsOm_Rloc, dVsOm_Mloc)
            runStop = MPI_Wtime()
            if (runStop>runStart) then
               run_time_mpi_comms=run_time_mpi_comms+(runStop-runStart)
            end if

            if ( tscheme%istage == 1 ) then

               !-------------------
               !------ Checking Courant criteria, l_new_dt and dt_new are output
               !-------------------
               call dt_courant(dtr,dth,l_new_dt,tscheme%dt(1),dt_new,dtMax, &
                    &          dtr_Rloc,dth_Rloc,time)

               call tscheme%set_dt_array(dt_new,dtMin,time,n_log_file,n_time_step, l_new_dt)
               !-- Store the old weight factor of matrices
               !-- it it changes because of dt factors moving
               !-- matrix also needs to be rebuilt
               call tscheme%set_weights(lMatNext)

               !----- Advancing time:
               time=time+tscheme%dt(1) ! Update time

               lMat=lMatNext
               if ( l_new_dt .or. lMat ) then
                  !----- Calculate matricies for new time step if dt /= dtLast
                  lMat=.true.
                  if ( rank == 0 ) then
                     write(*,'(1p,'' ! Building matricies at time step:'',   &
                          &              i8,ES16.6)') n_time_step,time
                  end if
               end if
               lMatNext = .false.

            end if

            !-- If the scheme is a multi-step scheme that is not Crank-Nicolson 
            !-- we have to use a different starting scheme
            call start_from_another_scheme(l_vphi_bal_calc, l_bridge_step, &
                 &                         n_time_step, tscheme)

            !--------------------
            !-- M-loop (update routines)
            !--------------------
            runStart = MPI_Wtime()
            call mloop(temp_Mloc, dtemp_Mloc, psi_Mloc, om_Mloc,  dom_Mloc, us_Mloc, up_Mloc,&
                 &     dVsT_Mloc, dVsOm_Mloc, buo_Mloc, dTdt, dpsidt, vp_bal, tscheme,       &
                 &     lMat, l_log_next, l_vphi_bal_calc, run_time_solve, n_solve_calls,     &
                 &     run_time_lu, n_lu_calls, run_time_dct,  n_dct_calls)
            runStop = MPI_Wtime()
            if ( .not. lMat ) then
               if (runStop>runStart) then
                  n_m_loops  =n_m_loops+1
                  run_time_m_loop=run_time_m_loop+(runStop-runStart)
               end if
            else
               if (runStop>runStart) then
                  n_m_loops_mat  =n_m_loops_mat+1
                  run_time_m_loop_mat=run_time_m_loop_mat+(runStop-runStart)
               end if
            end if

            ! Increment current stage
            tscheme%istage = tscheme%istage+1

         end do ! Finish loop over stages
         !---------------------
         !-- Timings
         !---------------------
         runStopT = MPI_Wtime()
         if (runStopT>runStartT) then
            run_time_tot=run_time_tot+(runStopT-runStartT)
         end if
         n_time_steps_go = n_time_steps_go+1

         if ( l_vphi_bal_calc ) call terminate_vp_bal(up_Mloc, vp_bal, tscheme)

         !---------------------
         !-- Info about run advance
         !---------------------
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
            run_time_passed=run_time_tot
            run_time_passed = run_time_passed/n_time_steps_go
            if ( rank == 0 ) then
               call formatTime(6,' ! Mean wall time for time step:',  &
               &               run_time_passed)
               call formatTime(n_log_file,' ! Mean wall time for time step:', &
               &               run_time_passed)
            end if

         end if

      end do outer ! end of time stepping !

      !--------------
      !-- Calculate wall time for different part of the code
      !-- and average over the different ranks
      !--------------
      run_time_io        = run_time_io/n_io_calls
      call my_reduce_mean(run_time_io, 0)
      run_time_r_loop    = run_time_r_loop/n_r_loops
      call my_reduce_mean(run_time_r_loop, 0)
      if ( n_m_loops /= 0 ) then
         run_time_m_loop    = run_time_m_loop/n_m_loops
         call my_reduce_mean(run_time_m_loop, 0)
      end if
      if ( n_m_loops_mat /= 0 ) then
         run_time_m_loop_mat= run_time_m_loop_mat/n_m_loops_mat
         call my_reduce_mean(run_time_m_loop_mat, 0)
      end if
      run_time_mpi_comms = run_time_mpi_comms/n_mpi_comms
      call my_reduce_mean(run_time_mpi_comms, 0)
      if ( n_time_steps_go /= 0 ) then
         run_time_tot       = run_time_tot/n_time_steps_go
         call my_reduce_mean(run_time_tot, 0)
      end if
      run_time_fft = run_time_fft/n_fft_calls
      call my_reduce_mean(run_time_fft, 0)
      run_time_lu = run_time_lu/n_lu_calls
      call my_reduce_mean(run_time_lu, 0)
      run_time_solve = run_time_solve/n_solve_calls
      call my_reduce_mean(run_time_solve, 0)
      run_time_dct = run_time_dct/n_dct_calls
      call my_reduce_mean(run_time_dct, 0)

      if ( rank == 0 ) then
         call formatTime(6, &
         &    '! Mean wall time for radial loop            :',run_time_r_loop)
         call formatTime(6, &
         &    '! Mean wall time for pure m loop            :',run_time_m_loop)
         call formatTime(6, &
         &    '! Mean wall time for m loop with matrix calc:',run_time_m_loop_mat)
         call formatTime(6, &
         &    '! Mean wall time for MPI communications     :',run_time_mpi_comms)
         call formatTime(6, &
         &    '! Mean wall time for output writting        :',run_time_io)
         call formatTime(6, &
         &    '! Mean wall time for one single FFT (rloop) :',run_time_fft)
         call formatTime(6, &
         &    '! Mean wall time for one single 2D-DCT      :',run_time_dct)
         call formatTime(6, &
         &    '! Mean wall time for one LU factor. (psi)   :',run_time_lu)
         call formatTime(6, &
         &    '! Mean wall time for one linear solve (psi) :',run_time_solve)
         call formatTime(6, &
         &    '! Mean wall time for one time step          :',run_time_tot)

         call formatTime(n_log_file, &
         &    '! Mean wall time for radial loop            :',run_time_r_loop)
         call formatTime(n_log_file,  &
         &    '! Mean wall time for pure m loop            :',run_time_m_loop)
         call formatTime(n_log_file, &
         &    '! Mean wall time for MPI communications     :',run_time_mpi_comms)
         call formatTime(n_log_file, &
         &    '! Mean wall time for m loop with matrix calc:',run_time_m_loop_mat)
         call formatTime(n_log_file, &
         &    '! Mean wall time for output writting        :',run_time_io)
         call formatTime(n_log_file, &
         &    '! Mean wall time for one single FFT (rloop) :',run_time_fft)
         call formatTime(n_log_file, &
         &    '! Mean wall time for one single 2D-DCT      :',run_time_dct)
         call formatTime(n_log_file, &
         &    '! Mean wall time for one LU factor. (psi)   :',run_time_lu)
         call formatTime(n_log_file, &
         &    '! Mean wall time for one linear solve (psi) :',run_time_solve)
         call formatTime(n_log_file,  &
         &    '! Mean wall time for one time step          :',run_time_tot)
      end if

   end subroutine time_loop
!--------------------------------------------------------------------------------
   subroutine start_from_another_scheme(l_vphi_bal_calc, l_bridge_step, &
              &                         n_time_step, tscheme)

      logical,             intent(in) :: l_vphi_bal_calc
      logical,             intent(in) :: l_bridge_step
      integer,             intent(in) :: n_time_step
      class(type_tscheme), intent(inout) :: tscheme

      !-- If the scheme is a multi-step scheme that is not Crank-Nicolson 
      !-- we have to use a different starting scheme
      if ( l_bridge_step .and. tscheme%time_scheme /= 'CNAB2' .and.  &
           n_time_step <= tscheme%norder_imp-2 .and.                 &
           tscheme%family=='MULTISTEP' ) then

         if ( l_cheb_coll ) then
            call get_temp_rhs_imp_coll(temp_Mloc,dtemp_Mloc, dTdt%old(:,:,1), &
                 &                     dTdt%impl(:,:,1),.true.)
            if ( l_direct_solve ) then
               call get_psi_rhs_imp_coll_smat(us_Mloc, up_Mloc, om_Mloc,   &
                    &                         dom_Mloc, dpsidt%old(:,:,1), &
                    &                         dpsidt%impl(:,:,1), vp_bal,  &
                    &                         l_vphi_bal_calc, .true.)
            else
               call get_psi_rhs_imp_coll_dmat(up_Mloc, om_Mloc, dom_Mloc,  &
                    &                         dpsidt%old(:,:,1),           &
                    &                         dpsidt%impl(:,:,1), vp_bal,  &
                    &                         l_vphi_bal_calc, .true.)
            end if
         else
            call get_temp_rhs_imp_int(temp_Mloc, dTdt%old(:,:,1), &
                 &                    dTdt%impl(:,:,1), .true.)
            if ( l_direct_solve ) then
               call get_psi_rhs_imp_int_smat(psi_Mloc,up_Mloc,dpsidt%old(:,:,1),&
                    &         dpsidt%impl(:,:,1), vp_bal, l_vphi_bal_calc,.true.)
            else
               call get_psi_rhs_imp_int_dmat(om_Mloc,up_Mloc,dpsidt%old(:,:,1), &
                    &         dpsidt%impl(:,:,1), vp_bal, l_vphi_bal_calc,.true.)
            end if
         end if

         call tscheme%bridge_with_cnab2()

      end if

      if ( l_AB1 .and. n_time_step == 1 ) then
         call tscheme%start_with_ab1()
         l_AB1 = .false.
      end if

   end subroutine start_from_another_scheme
!---------------------------------------------------------------------------------
end module step_time
