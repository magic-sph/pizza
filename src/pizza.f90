program pizza

   use iso_fortran_env, only: output_unit
   use parallel_mod
   use precision_mod
   use mem_alloc
   use step_time, only: time_loop
   use courant_mod, only: initialize_courant, finalize_courant
   use radial_der, only: initialize_der_arrays, finalize_der_arrays
   use init_fields, only: get_start_fields
   use fields, only: initialize_fields, finalize_fields
   use fieldsLast, only: initialize_fieldsLast, finalize_fieldsLast
   use communications, only: initialize_communications, finalize_communications
   use blocking, only: set_mpi_domains, nMstart, nMstop, destroy_mpi_domains, &
       &               nRstart, nRstop
   use namelists, only: read_namelists, write_namelists, tag, time_scheme,    &  
       &                l_cheb_coll, l_rerror_fix, rerror_fac, l_direct_solve,&
       &                courfac, l_heat, l_chem, l_finite_diff
   use mloop_fd_mod, only: initialize_mloop_fd, finalize_mloop_fd, test_mloop
   use outputs, only: initialize_outputs, finalize_outputs, n_log_file
   use pre_calculations, only: preCalc
   use horizontal, only: initialize_mfunctions, finalize_mfunctions
   use radial_functions, only: initialize_radial_functions, &
       &                       finalize_radial_functions
   use truncation, only: initialize_truncation, n_r_max, n_phi_max, & 
       &                 finalize_truncation, n_m_max
   use fourier, only: initialize_fourier, finalize_fourier
   use rloop, only: initialize_radial_loop, finalize_radial_loop
   use update_temp_coll, only: initialize_temp_coll, finalize_temp_coll
   use update_xi_coll, only: initialize_xi_coll, finalize_xi_coll
   use update_temp_integ, only: initialize_temp_integ, finalize_temp_integ
   use update_xi_integ, only: initialize_xi_integ, finalize_xi_integ
   use update_psi_coll_smat, only: initialize_om_coll_smat, finalize_om_coll_smat
   use update_psi_coll_dmat, only: initialize_om_coll_dmat, finalize_om_coll_dmat
   use update_psi_integ_smat, only: initialize_psi_integ_smat, &
       &                            finalize_psi_integ_smat
   use update_psi_integ_dmat, only: initialize_psi_integ_dmat, &
       &                            finalize_psi_integ_dmat
   use select_time_scheme
   use time_schemes, only: type_tscheme
   use useful, only: formatTime
   use tests, only: solve_laplacian, test_radial_der, solve_biharmo, test_i4

   implicit none

   real(cp) :: time
   class(type_tscheme), pointer :: tscheme
   real(cp) :: run_stop, run_start, run_init
   integer(lip) :: local_bytes_used
   integer :: values(8)
   integer :: n, n_out
   character(len=72) :: date

   !-- Initialize MPI
   call initialize_mpi()

   run_start = MPI_Wtime()

   !-- Read input parameters
   call read_namelists()

   !-- Select the kind of time-integrator (multi-step or implicit R-K):
   call select_tscheme(time_scheme, tscheme)

   !-- Init memory counter
   call initialize_memory_counter(tag)

   !-- Initialize truncation
   call initialize_truncation()

   !-- Set the domain decomposition
   call set_mpi_domains()

   !-- Test radial derivatives
   !call test_i4()
   !call solve_laplacian(nMstart, nMstop)
   !call solve_biharmo(nMstart, nMstop)
   !call test_radial_der(nMstart,nMstop)
   !stop

   !-- Open output files
   local_bytes_used = bytes_allocated
   call initialize_outputs()
   local_bytes_used = bytes_allocated-local_bytes_used
   call memWrite('I/O', local_bytes_used)

   !-- Write header of log file
   if ( rank == 0 ) then
      do n=1,2
         if ( n==1 ) n_out=output_unit
         if ( n==2 ) n_out=n_log_file
         write(n_out,*)
         write(n_out,*) '!--- Program pizza  ---!'
         call date_and_time(values=values)
         write(date, '(i4,''/'',i0.2,''/'',i0.2,'' '', i0.2,'':'',i0.2,'':'',i0.2)') &
         &     values(1), values(2), values(3), values(5), values(6), values(7)
         write(n_out, *) '!  Start time:  ', date
      end do
   end if

   !-- Initialize time scheme
   call tscheme%initialize(time_scheme, courfac)

   !-- Initialize MPI communicators
   call initialize_communications(n_log_file)

   local_bytes_used = bytes_allocated
   call initialize_fields()
   call initialize_fieldsLast(nMstart, nMstop, n_m_max, nRstart, nRstop, n_r_max, &
        &                     tscheme%nold, tscheme%nexp, tscheme%nimp)
   local_bytes_used = bytes_allocated-local_bytes_used
   call memWrite('Fields', local_bytes_used)
   local_bytes_used = bytes_allocated
   call initialize_radial_functions()
   call initialize_der_arrays(n_r_max, nMstart, nMstop, l_rerror_fix, rerror_fac)
   local_bytes_used = bytes_allocated-local_bytes_used
   call memWrite('Radial functions', local_bytes_used)
   local_bytes_used = bytes_allocated
   call initialize_fourier(n_phi_max)
   call initialize_radial_loop(n_phi_max)
   local_bytes_used = bytes_allocated-local_bytes_used
   call memWrite('R loop', local_bytes_used)
   call initialize_mfunctions()


   if ( rank == 0 ) then
      call write_namelists(output_unit)
      call write_namelists(n_log_file)
      call tscheme%print_info(n_log_file)
   end if

   !-- Pre calculations has to be done before matrix initialisation
   call preCalc()

   local_bytes_used = bytes_allocated
   if ( l_finite_diff ) then
      call initialize_mloop_fd(tscheme)
      call test_mloop(tscheme) ! find faster block layout
   else
      if ( l_cheb_coll ) then
         if ( l_heat ) call initialize_temp_coll()
         if ( l_chem ) call initialize_xi_coll()
         if ( l_direct_solve ) then
            call initialize_om_coll_smat(tscheme)
         else
            call initialize_om_coll_dmat(tscheme)
         end if
      else
         if ( l_heat ) call initialize_temp_integ(tscheme)
         if ( l_chem ) call initialize_xi_integ(tscheme)
         if ( l_direct_solve ) then
            call initialize_psi_integ_smat(tscheme)
         else
            call initialize_psi_integ_dmat()
         end if
      end if
   end if
   local_bytes_used = bytes_allocated-local_bytes_used
   call memWrite('M loop', local_bytes_used)
   call finalize_memory_counter()

   !-- Start fields
   call get_start_fields(time, tscheme)

   !-- Open time step file
   call initialize_courant(time, tscheme%dt(1))

   !--- Write starting time to SDTOUT and logfile:
   if ( rank == 0 ) then
      do n=1,2
         if ( n == 1 ) n_out=output_unit
         if ( n == 2 ) n_out=n_log_file
         write(n_out,'(/,'' ! Starting time integration at:'')')
         write(n_out,'(''   start_time ='',1p,ES18.10)') time
         write(n_out,'(''   start dt   ='',1p,ES16.4)') tscheme%dt(1)
      end do
   end if
   run_init = MPI_Wtime()
   run_init = run_init - run_start
   call MPI_Allreduce(MPI_IN_PLACE,run_init,1,MPI_DEF_REAL,MPI_MAX, &
        &             MPI_COMM_WORLD, ierr)

   !-- Time integration
   call time_loop(time, tscheme, run_init)

   run_stop = MPI_Wtime()

   !--- Write stop time to SDTOUR and logfile:
   if ( rank == 0 ) then
      do n=1,2
         if ( n == 1 ) n_out=output_unit
         if ( n == 2 ) n_out=n_log_file
         write(n_out,'(/,'' ! Stopping time integration at:'')')
         write(n_out,'(''   stop time ='',1p,ES18.10)') time
         ! write(n_out,'(''   stop step ='',i10)') n_time_step
         ! write(n_out,'(''   steps gone='',i10)') (n_time_step-1)
         write(n_out,*)
         write(n_out,*)
         call formatTime(n_out,'! Total run time:', run_stop-run_start)
         write(n_out,*) '!-----------------------------------!'
         write(n_out,*) '!---- Ready for another pizza ? ----!'
         write(n_out,*) '!-----------------------------------!'
         write(n_out,*)
      end do
   end if

   !-- Close files
   if ( l_finite_diff ) then
      call finalize_mloop_fd(tscheme)
   else
      if ( l_cheb_coll ) then
         if ( l_heat ) call finalize_temp_coll()
         if ( l_chem ) call finalize_xi_coll()
         if ( l_direct_solve ) then
            call finalize_om_coll_smat(tscheme)
         else
            call finalize_om_coll_dmat(tscheme)
         end if
      else
         if ( l_heat ) call finalize_temp_integ(tscheme)
         if ( l_chem ) call finalize_xi_integ(tscheme)
         if ( l_direct_solve ) then
            call finalize_psi_integ_smat(tscheme)
         else
            call finalize_psi_integ_dmat()
         end if
      end if
   end if
   call finalize_mfunctions()
   call finalize_radial_loop()
   call finalize_fourier()
   call finalize_der_arrays()
   call finalize_radial_functions()
   call finalize_fieldsLast()
   call finalize_fields()
   call finalize_communications()
   call destroy_mpi_domains()
   call finalize_courant()
   call finalize_outputs()
   call tscheme%finalize()
   call finalize_truncation()

   !-- Initialize MPI
   call finalize_mpi()

end program pizza
