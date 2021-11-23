module mloop_fd_mod

   use iso_fortran_env, only: output_unit
   use precision_mod
   use parallel_mod
   use constants, only: one
   use useful, only: abortRun, logWrite
   use outputs, only: n_log_file
   use namelists, only: l_heat, l_chem
   use truncation, only: n_m_max
   use blocking, only: nRstart, nRstop
   use time_array, only: type_tarray
   use vp_balance, only: vp_bal_type
   use vort_balance, only: vort_bal_type
   use time_schemes, only: type_tscheme
   use timers_mod, only: timers_type, timer_type
   use update_temp_fd_mod
   use update_xi_fd_mod
   use update_psi_fd_mod

   implicit none

   private

   integer :: n_tri, n_penta ! Number of tri and pentadiagonal solvers
   integer :: block_sze, n_requests, nblocks
   integer, allocatable :: array_of_requests(:)

   public :: finish_explicit_assembly_Rdist, mloop_Rdist, initialize_mloop_fd, &
   &         finalize_mloop_fd, test_mloop, assemble_stage_Rdist

contains

   subroutine initialize_mloop_fd(tscheme)

      class(type_tscheme), intent(in) :: tscheme

      n_tri=0
      n_penta=1   
      if ( l_heat ) then
         call initialize_temp_fd()
         lTmat_FD(:)=.false.
         n_tri = n_tri+1
      end if
      if ( l_chem ) then
         call initialize_xi_fd()
         lXimat_FD(:)=.false.
         n_tri = n_tri+1
      end if
      call initialize_psi_fd(tscheme)
      lPsimat_FD(:)=.false.

      block_sze=50
      n_requests=10
      nblocks = n_m_max
      nblocks = set_block_number(nblocks)
      allocate( array_of_requests(n_requests))

   end subroutine initialize_mloop_fd
!------------------------------------------------------------------------------------
   subroutine finalize_mloop_fd(tscheme)

      class(type_tscheme), intent(in) :: tscheme

      deallocate( array_of_requests )
      if ( l_heat ) call finalize_temp_fd()
      if ( l_chem ) call finalize_xi_fd()
      call finalize_psi_fd(tscheme)

   end subroutine finalize_mloop_fd
!------------------------------------------------------------------------------------
   subroutine finish_explicit_assembly_Rdist(us_Rloc,temp_Rloc,xi_Rloc,dVsT_Rloc,    &
              &                              dVsXi_Rloc,dVsOm_Rloc,dTdt,dxidt,dpsidt,&
              &                              tscheme, vp_bal, vort_bal)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      complex(cp),         intent(in) :: us_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(in) :: temp_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(in) :: xi_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: dVsT_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: dVsXi_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: dVsOm_Rloc(n_m_max,nRstart:nRstop)

      !-- Output variables
      type(type_tarray),   intent(inout) :: dTdt
      type(type_tarray),   intent(inout) :: dxidt
      type(type_tarray),   intent(inout) :: dpsidt
      type(vp_bal_type),   intent(inout) :: vp_bal
      type(vort_bal_type), intent(inout) :: vort_bal

      call finish_exp_psi_Rdist(us_Rloc, dVsOm_Rloc, temp_Rloc, xi_Rloc,   &
           &                    dpsidt%expl(:,:,tscheme%istage), vp_bal,   &
           &                    vort_bal, tscheme)

      if ( l_heat ) then
         call finish_exp_temp_Rdist(us_Rloc, dVsT_Rloc, dTdt%expl(:,:,tscheme%istage))
      end if

      if ( l_chem ) then
         call finish_exp_xi_Rdist(us_Rloc, dVsXi_Rloc, dxidt%expl(:,:,tscheme%istage))
      end if

   end subroutine finish_explicit_assembly_Rdist
!------------------------------------------------------------------------------------
   subroutine test_mloop(tscheme)
      !
      ! This subroutine is used to solve dummy linear problem to estimate the best
      ! blocking size. This is done once at the initialisation stage of pizza.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme ! time scheme

      !-- Local variable
      type(type_tarray) :: dummy
      real(cp) :: dummy_time
      integer :: dummy_counter

      lPsimat_FD(:)=.false.
      if ( l_heat ) lTmat_FD(:) =.false.
      if ( l_chem ) lXimat_FD(:) =.false.

      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call dummy%initialize(1, n_m_max, nRstart, nRstop, tscheme%nold, tscheme%nexp,&
           &                tscheme%nimp, l_allocate_exp=.true.)

      if ( l_heat ) call prepare_temp_FD(tscheme, dummy)
      if ( l_chem ) call prepare_xi_FD(tscheme, dummy)
      call prepare_psi_fd(tscheme, dummy, dummy_time, dummy_counter)

      call find_faster_block() ! Find the fastest blocking

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call dummy%finalize()

   end subroutine test_mloop
!------------------------------------------------------------------------------------
   subroutine mloop_Rdist(temp_Rloc, dtemp_Rloc, xi_Rloc, dxi_Rloc, om_Rloc, &
              &           us_Rloc, up_Rloc, dTdt, dxidt, dpsidt, vp_bal,     &
              &           vort_bal, tscheme, lMat, timers)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lMat

      !-- Output variables
      complex(cp),         intent(out) :: temp_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(out) :: dtemp_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(out) :: xi_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(out) :: dxi_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(out) :: om_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: us_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: up_Rloc(n_m_max,nRstart:nRstop)
      type(vp_bal_type),   intent(inout) :: vp_bal
      type(vort_bal_type), intent(inout) :: vort_bal
      type(type_tarray),   intent(inout) :: dpsidt
      type(type_tarray),   intent(inout) :: dTdt
      type(type_tarray),   intent(inout) :: dxidt
      type(timers_type),   intent(inout) :: timers

      !-- Local variables
      real(cp) :: runStart, runStop

      if ( vp_bal%l_calc .and. tscheme%istage==1 ) then
         call vp_bal%initialize_dvpdt(up_Rloc, tscheme)
      end if

      if ( vort_bal%l_calc .and. tscheme%istage==1 ) then
         call vort_bal%initialize_domdt(om_Rloc, tscheme)
      end if

      if ( lMat ) then ! update matrices
         lPsimat_FD(:)=.false.
         if ( l_heat ) lTmat_FD(:)=.false.
         if ( l_chem ) lXimat_FD(:)=.false.
      end if

      !-- Mainly assemble the r.h.s. and rebuild the matrices if required
      if ( l_heat ) call prepare_temp_FD(tscheme, dTdt)
      if ( l_chem ) call prepare_xi_FD(tscheme, dxidt)
      call prepare_psi_FD(tscheme, dpsidt, timers%lu, timers%n_lu_calls)

      !-- solve the uphi0 equation on its own (one single m is involved)
      call upMat_FD%solver_single(up0_ghost, nRstart, nRstop)

      !-----------------------------------------------------------
      !--- This is where the matrices are solved
      !-- Here comes the real deal:
      runStart=MPI_Wtime()
      call parallel_solve(block_sze)
      runStop = MPI_Wtime()
      if ( runStop > runStart ) then
         timers%solve = timers%solve + (runStop-runStart)
         timers%n_solve_calls = timers%n_solve_calls+1
      end if
      !-----------------------------------------------------------

      !-- Now simply fill the ghost zones to ensure the boundary conditions
      if ( l_heat ) call fill_ghosts_temp(temp_ghost)
      if ( l_chem ) call fill_ghosts_xi(xi_ghost)
      call fill_ghosts_psi(psi_ghost, up0_ghost)

      !-- Finally build the radial derivatives and the arrays for next iteration
      if ( l_heat ) call update_temp_FD(temp_Rloc, dtemp_Rloc, dTdt, tscheme)
      if ( l_chem ) call update_xi_FD(xi_Rloc, dxi_Rloc, dxidt, tscheme)
      call update_psi_FD(us_Rloc, up_Rloc, om_Rloc, dpsidt, tscheme, vp_bal, &
           &             vort_bal)

      if ( vp_bal%l_calc .and. tscheme%istage==tscheme%nstages ) then
         call vp_bal%finalize_dvpdt(up_Rloc, tscheme)
      end if

      if ( vort_bal%l_calc .and. tscheme%istage==tscheme%nstages ) then
         call vort_bal%finalize_domdt(om_Rloc, tscheme)
      end if

   end subroutine mloop_Rdist
!------------------------------------------------------------------------------------
   subroutine assemble_stage_Rdist(temp_Rloc, dtemp_Rloc, xi_Rloc, dxi_Rloc, us_Rloc,&
              &                    up_Rloc, om_Rloc, dTdt, dxidt, dpsidt, tscheme,   &
              &                    vp_bal, vort_bal)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp),         intent(out) :: temp_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(out) :: dtemp_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(out) :: xi_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(out) :: dxi_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(out) :: om_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: us_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: up_Rloc(n_m_max,nRstart:nRstop)
      type(vp_bal_type),   intent(inout) :: vp_bal
      type(vort_bal_type), intent(inout) :: vort_bal
      type(type_tarray),   intent(inout) :: dpsidt
      type(type_tarray),   intent(inout) :: dTdt
      type(type_tarray),   intent(inout) :: dxidt

      if ( l_heat ) call assemble_temp_Rdist(temp_Rloc, dtemp_Rloc, dTdt, tscheme)
      if ( l_chem ) call assemble_xi_Rdist(xi_Rloc, dxi_Rloc, dxidt, tscheme)
      call assemble_psi_Rloc(block_sze, nblocks, us_Rloc, up_Rloc, om_Rloc, dpsidt, &
           &                 tscheme, vp_bal, vort_bal)

   end subroutine assemble_stage_Rdist
!------------------------------------------------------------------------------------
   subroutine parallel_solve(block_sze)
      !
      ! This subroutine handles the parallel solve of the time-advance matrices.
      ! This works with R-distributed arrays (finite differences).
      !
      integer, intent(in) :: block_sze ! Size ot the LM blocks

      !-- Local variables
      integer :: req
      integer :: start_m, stop_m, tag, n_m_block, n_ms_block

      array_of_requests(:)=MPI_REQUEST_NULL

      !$omp parallel default(shared) private(tag, req, start_m, stop_m)
      tag = 0
      req=1

      do n_ms_block=1,n_m_max,block_sze
         n_m_block = n_m_max-n_ms_block+1
         if ( n_m_block > block_sze ) n_m_block=block_sze
         start_m=n_ms_block; stop_m=n_ms_block+n_m_block-1
         call get_openmp_blocks(start_m,stop_m)
         !$omp barrier

         if ( l_heat ) then
            call tMat_FD%solver_up(temp_ghost, start_m, stop_m, nRstart, nRstop, tag, &
                 &                 array_of_requests, req, n_ms_block, n_m_block)
            tag = tag+1
         end if

         if ( l_chem ) then
            call xiMat_FD%solver_up(xi_ghost, start_m, stop_m, nRstart, nRstop, tag, &
                 &                  array_of_requests, req, n_ms_block, n_m_block)
            tag = tag+1
         end if

         call psiMat_FD%solver_up(psi_ghost, start_m, stop_m, nRstart, nRstop, tag, &
              &                   array_of_requests, req, n_ms_block, n_m_block)
         tag = tag+2
      end do

      do n_ms_block=1,n_m_max,block_sze
         n_m_block = n_m_max-n_ms_block+1
         if ( n_m_block > block_sze ) n_m_block=block_sze
         start_m=n_ms_block; stop_m=n_ms_block+n_m_block-1
         call get_openmp_blocks(start_m,stop_m)
         !$omp barrier

         if ( l_heat ) then
            call tMat_FD%solver_dn(temp_ghost, start_m, stop_m, nRstart, nRstop, tag, &
              &                    array_of_requests, req, n_ms_block, n_m_block)
            tag = tag+1
         end if

         if ( l_chem ) then
            call xiMat_FD%solver_dn(xi_ghost, start_m, stop_m, nRstart, nRstop, tag, &
                 &                  array_of_requests, req, n_ms_block, n_m_block)
            tag = tag+1
         end if

         call psiMat_FD%solver_dn(psi_ghost, start_m, stop_m, nRstart, nRstop, tag, &
              &                   array_of_requests, req, n_ms_block, n_m_block)
         tag = tag+2
      end do

      !$omp master
      do n_ms_block=1,n_m_max,block_sze
         n_m_block = n_m_max-n_ms_block+1
         if ( n_m_block > block_sze ) n_m_block=block_sze

         if ( l_heat ) then
            call tMat_FD%solver_finish(temp_ghost, n_ms_block, n_m_block, nRstart, nRstop, &
                 &                     tag, array_of_requests, req)
            tag = tag+1
         end if

         if ( l_chem ) then
            call xiMat_FD%solver_finish(xi_ghost, n_ms_block, n_m_block, nRstart, nRstop, &
                 &                      tag, array_of_requests, req)
            tag = tag+1
         end if

         call psiMat_FD%solver_finish(psi_ghost, n_ms_block, n_m_block, nRstart, nRstop, &
              &                       tag, array_of_requests, req)
         tag = tag+2
      end do

      call MPI_Waitall(req-1, array_of_requests(1:req-1), MPI_STATUSES_IGNORE, ierr)
      if ( ierr /= MPI_SUCCESS ) call abortRun('MPI_Waitall failed in LMLoop')
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      !$omp end master
      !$omp barrier

      !$omp end parallel

   end subroutine parallel_solve
!------------------------------------------------------------------------------------
   integer function set_block_number(nb) result(nbo)
      !
      ! This routine returns the number of lm-blocks for solving the LM loop.
      ! This is adapted from xshells.
      !
      integer, intent(inout) :: nb ! Number of blocks

      !-- Local variable
      integer :: nlm_4

      nlm_4=(n_m_max+3)/4 ! Don't split cache lines
      if ( nb < 1 ) nb=1
      if ( (nlm_4+nb-1)/nb < 16) nb=(nlm_4+15)/16 ! Size not less than 1KB (64 cplx)
      if ( nb > 2048 ) nb=2048 ! No more than 2048 blocks
      block_sze=((nlm_4+nb-1)/nb)*4 ! Block size (rounded)

      nblocks = (n_m_max+block_sze-1)/block_sze
      n_requests = ( n_tri*4+n_penta*8)*nblocks ! Max number of MPI request per block
      nbo = nb

      if ( allocated(array_of_requests) ) then
         deallocate(array_of_requests)
         allocate(array_of_requests(n_requests))
      end if

   end function set_block_number
!------------------------------------------------------------------------------------
   subroutine find_faster_block
      !
      ! This routine is used to find the best lm-block size for MPI communication.
      ! This is adapted from xshells.
      !
      real(cp), parameter :: alpha=0.61803_cp
      real(cp), parameter :: beta=one-alpha
      integer :: bmin, bmax, nb, nblk, b(0:2), nloops, bb
      real(cp) :: t(0:2), mul, tt
      character(len=72) :: str

      call logWrite('', n_log_file)
      call logWrite('! Time the LM Loop to get the optimal blocking:', n_log_file)

      t(:)=0.0_cp; b(:)=1
      nloops=1
      bmin = 1 ! Minimum number of blocks
      nb = n_m_max
      bmax =  set_block_number(nb) ! Maximum number of blocks
      nb = (n_procs+1)/2 ! Start with a default block number
      nblk = set_block_number(nb)

      if ( n_procs > 1 ) then
         nb=int(bmin**beta * bmax**alpha)
         nblk = set_block_number(nb)
      end if
      b(1) = set_block_number(nblk)

      !-- First estimate the number of loops to get accurate timing
      !-- To do so measure the walltime until this is stable to a 2% accuracy
      t(1) = time_M_loop(b(1), nloops) ! Starting value
      tt=t(1)+one
      do while( abs(t(1)-tt) > 0.02*abs(t(1)) ) ! 2% difference
         nloops = 2*nloops
         tt = t(1)
         t(1) = time_M_loop(b(1), nloops)
         if ( t(1)*(nloops/2) > 10 ) exit ! Longer than 10sec
      end do
      nloops=nloops/2
      nloops=max(nloops,2) ! At least 2 loops
      write(str,'(A,I5,A)') '! I am computing', nloops, ' calls'
      call logWrite(str, n_log_file)

      if ( n_procs > 1) then
         if ( bmax*bmin > b(1)*b(1) ) then
            mul = real(bmax,cp)
         else
            mul = real(bmin,cp)
         end if
         mul = (mul/b(1))**beta
         nb = int(b(1)*mul)
         b(0)=set_block_number(nb)
         t(0)=time_M_loop(b(0), nloops)
         if ( t(0) < t(1) ) then ! Exchange indices
            bb=b(0); b(0)=b(1); b(1)=bb
            tt=t(0); t(0)=t(1); t(1)=tt
         end if

         if ( b(1) > b(0) ) then
            mul = real(bmax,cp)
         else
            mul = real(bmin,cp)
         end if
         mul = (mul/b(1))**beta
         b(2)=b(1)

         do while ( (b(1) > bmin) .and. ( b(1) < bmax) )
            nb=int(b(2)*mul)
            b(2)=set_block_number(nb)
            t(2)=time_M_loop(b(2), nloops)
            if ( t(2) < t(1) ) then
               if ( t(1) > t(2)*1.02) then ! If we are within 2% of the minimum time
                  !-- Change the interval boundaries
                  t(0)=t(1); b(0)=b(1)
               end if
               t(1)=t(2); b(1)=b(2)
            end if
            if ( (b(2)==bmin) .or. (b(2)==bmax) .or. (t(2)>t(1)*1.02 ) ) exit
         end do
         write(str,'(A,1x,I0,A,1x,I0)') '! I am braketing the block number between', &
         &                              b(0), ' and', b(2)
         call logWrite(str, n_log_file)

         !-- At this stage, the minimum wall time has been stored in t(1) and is bracketed
         !-- between in t(0) and t(2)
         !if ( rank == 0 ) print*, 'b', 't', b, t

         !-- Determine the largest interval
         if ( abs(log(real(b(2)))-log(real(b(0)))) > abs(log(real(b(0)))-log(real(b(1)))) ) then
            nb=int(b(1)**alpha*b(2)**beta)
            bb=set_block_number(nb)
            tt=time_M_loop(bb, nloops)
         else
            bb=b(1); tt=t(1)
            nb=int(b(1)**alpha*b(0)**beta)
            t(1)=time_M_loop(b(1), nloops)
         end if

         !-- Refined block determination
         if ( rank == 0 ) write(output_unit,*) ''
         do while( (t(2)>min(t(1),tt)*1.02) .and. (t(0)>min(t(1),tt)*1.02) &
         &          .and. (abs(b(2)-b(0))>1) .and. (maxval(b) < 4) )
            if ( tt < t(1) ) then ! This is better than before
               t(0)=t(1); b(0)=b(1);
               t(1)=tt; b(1)=bb
               nb=int(b(1)**alpha*b(2)**beta)
               bb=set_block_number(nb)
               tt=time_M_loop(bb, nloops)
            else
               t(2)=tt; b(2)=bb
               tt=t(1); bb=b(1)
               nb=int(b(1)**alpha*b(0)**beta)
               b(1)=set_block_number(nb)
               t(1)=time_M_loop(b(1), nloops)
            end if
            if ( rank==0 ) write(output_unit, '(A,I0,1x,I0,1x,I0,1x,I0,A,4ES9.2,A)')  &
            &              'Searching for number of blocks b=(',b(0), b(1), bb, b(2), &
            &              '), t=(',t(0), t(1), tt, t(2),')'
         end do
         if ( rank == 0 ) write(output_unit,*) ''

         if ( tt < t(1) ) t(1)=tt; b(1)=bb
         call MPI_Barrier(MPI_COMM_WORLD, ierr)

         !- At this stage this should be done: b(1) should be the optimal block_number
         nb=b(1)
         nblk=set_block_number(nb)
      end if

      !print*, 'b,t=', b, t
      if ( nblk == 1 ) then
         write(str,'(A,1x,I0,A,ES9.2)') '! Found 1 block of size',  &
         &                                 block_sze, ', timing was:', t(1)
      else
         write(str,'(A,1x,I0,A,1x,I0,A,ES9.2)') '! Found', nblk, ' blocks of size',  &
         &                                      block_sze, ', timing was:', t(1)
      end if
      call logWrite(str, n_log_file)

   end subroutine find_faster_block
!------------------------------------------------------------------------------------
   real(cp) function time_M_loop(nblk, nloops) result(time)
      !
      ! This function is defined to measure the time of a call to the parallel solvers
      ! when a a number of block nblk is used. This takes the average over nloop calls.
      !

      integer, intent(inout) :: nblk ! Number of lm blocks
      integer, intent(in) :: nloops ! Number of calls

      !-- Local variables
      type(timer_type) :: tcount
      integer :: n_l

      time = 0.0_cp
      nblk = set_block_number(nblk)
      call tcount%initialize()

      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      do n_l=1,nloops
         call tcount%start_count()
         call parallel_solve(block_sze)
         call tcount%stop_count()
      end do
      call tcount%finalize()

      time=tcount%tTot

   end function time_M_loop
!------------------------------------------------------------------------------------
end module mloop_fd_mod
