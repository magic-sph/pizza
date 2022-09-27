program pizza

   use parallel_mod
   use precision_mod
   use mem_alloc
   use iso_fortran_env, only: output_unit
   use step_time, only: time_loop
   use courant_mod, only: initialize_courant, finalize_courant
   use radial_der, only: initialize_der_arrays, finalize_der_arrays
   use init_fields, only: get_start_fields
   use fields, only: initialize_fields, finalize_fields
   use fieldsLast, only: initialize_fieldsLast, finalize_fieldsLast, &
       &                 initialize_fieldsLast_3D
   use communications, only: initialize_communications, finalize_communications
   use blocking, only: set_mpi_domains, nMstart, nMstop, destroy_mpi_domains, &
       &               nRstart, nRstop, nRstart3D, nRstop3D, lmStart, lmStop
   use blocking_lm, only: initialize_blocking, finalize_blocking
   use namelists, only: read_namelists, write_namelists, tag, time_scheme,    &  
       &                l_cheb_coll, l_rerror_fix, rerror_fac, l_direct_solve,&
       &                courfac, alffac, l_heat, l_chem, l_3D, l_heat_3D,     &
       &                l_mag_3D
   use outputs, only: initialize_outputs, finalize_outputs, n_log_file
   use outputs_3D, only: initialize_outputs_3D, finalize_outputs_3D
   use pre_calculations, only: preCalc
   use horizontal, only: initialize_mfunctions, finalize_mfunctions
   use radial_functions, only: initialize_radial_functions, r_3D, &
       &                       finalize_radial_functions
   use truncation, only: initialize_truncation, n_r_max, n_phi_max, & 
       &                 finalize_truncation, n_m_max
   use truncation_3D, only: initialize_truncation_3D, n_r_max_3D, lm_max, &
       &                    n_phi_max_3D, lmP_max
   !use z_functions, only: initialize_zfunctions, finalize_zfunctions
   use fourier, only: initialize_fourier, finalize_fourier
   use rloop, only: initialize_radial_loop, finalize_radial_loop
   use rloop_3D, only: initialize_radial_loop_3D, finalize_radial_loop_3D
   use update_temp_coll, only: initialize_temp_coll, finalize_temp_coll
   use update_xi_coll, only: initialize_xi_coll, finalize_xi_coll
   use update_temp_integ, only: initialize_temp_integ, finalize_temp_integ
   use update_temp_3D_mod, only: initialize_update_temp_3D, finalize_update_temp_3D
   use update_mag_3D_mod, only: initialize_update_mag_3D, finalize_update_mag_3D
   use update_xi_integ, only: initialize_xi_integ, finalize_xi_integ
   use update_psi_coll_smat, only: initialize_om_coll_smat, finalize_om_coll_smat
   use update_psi_coll_dmat, only: initialize_om_coll_dmat, finalize_om_coll_dmat
   use update_psi_integ_smat, only: initialize_psi_integ_smat, &
       &                            finalize_psi_integ_smat
   use update_psi_integ_dmat, only: initialize_psi_integ_dmat, &
       &                            finalize_psi_integ_dmat
   use select_time_scheme
   use time_schemes, only: type_tscheme
   use z_functions, only: zfunc_type
   use useful, only: formatTime
   use tests, only: solve_laplacian, test_radial_der, solve_biharmo, test_i4
#ifdef WITH_SHTNS
   use shtns, only: init_shtns
#endif

   implicit none

   real(cp) :: time
   class(type_tscheme), pointer :: tscheme
   type(zfunc_type) :: zinterp
   real(cp) :: run_stop, run_start, run_init
   integer(lip) :: local_bytes_used
   integer :: values(8)
   integer :: n, n_out
   character(len=72) :: date

   !-- Initialize MPI
   call initialize_mpi()

   run_start = MPI_Wtime()

   !--
   if ( rank == 0 ) then
      write(output_unit,*)
      write(output_unit,*) '!--- Program pizza  ---!'
      call date_and_time(values=values)
      write(date, '(i4,''/'',i0.2,''/'',i0.2,'' '', i0.2,'':'',i0.2,'':'',i0.2)') &
      &     values(1), values(2), values(3), values(5), values(6), values(7)
      write(output_unit, *) '!  Start time:  ', date

   end if

   !-- Read input parameters
   call read_namelists()

   !-- Select the kind of time-integrator (multi-step or implicit R-K):
   call select_tscheme(time_scheme, tscheme)

   !-- Init memory counter
   call initialize_memory_counter(tag)

   !-- Initialize truncation
   call initialize_truncation()

   !-- Initialize 3D truncation
   if ( l_3D ) call initialize_truncation_3D()

   !-- Set the domain decomposition
   call set_mpi_domains(l_3D)

   !-- Initialize 3D blocking
   if ( l_3D ) call initialize_blocking()

   if ( l_3D ) call zinterp%initialize()

   !-- Test radial derivatives
   !call test_i4()
   !call solve_laplacian(nMstart, nMstop)
   !call solve_biharmo(nMstart, nMstop)
   !call test_radial_der(nMstart,nMstop)
   !stop

   !-- Open output files
   local_bytes_used = bytes_allocated
   call initialize_outputs()
   if ( l_3D ) call initialize_outputs_3D()
   local_bytes_used = bytes_allocated-local_bytes_used
   call memWrite('I/O', local_bytes_used)


   !-- Initialize time scheme
   call tscheme%initialize(time_scheme, courfac, alffac)

   !-- Initialize MPI communicators
   call initialize_communications(l_3D)

   local_bytes_used = bytes_allocated
   call initialize_fields()
   call initialize_fieldsLast(nMstart, nMstop, n_m_max, nRstart, nRstop,       &
        &                     n_r_max, tscheme%norder_imp, tscheme%norder_exp, &
        &                     tscheme%norder_imp_lin)
   if ( l_3D ) then
      call initialize_fieldsLast_3D(lmStart, lmStop, lm_max, nRstart3D, nRstop3D, &
           &                        n_r_max_3D, tscheme%norder_imp,               &
           &                        tscheme%norder_exp, tscheme%norder_imp_lin)
   end if
   local_bytes_used = bytes_allocated-local_bytes_used
   call memWrite('Fields', local_bytes_used)
   local_bytes_used = bytes_allocated
   call initialize_radial_functions()
   local_bytes_used = bytes_allocated-local_bytes_used
   call memWrite('Radial functions', local_bytes_used)
   local_bytes_used = bytes_allocated
   call initialize_fourier(n_phi_max, n_phi_max_3D, l_3D)
   call initialize_radial_loop(n_phi_max)
   if ( l_3D ) call initialize_radial_loop_3D(lm_max,lmP_max)
   local_bytes_used = bytes_allocated-local_bytes_used
   call memWrite('R loop', local_bytes_used)
   call initialize_mfunctions()

#ifdef WITH_SHTNS
   if ( l_3D ) call init_shtns()
#endif

   if ( rank == 0 ) then
      call write_namelists(output_unit)
      call write_namelists(n_log_file)
      call tscheme%print_info(n_log_file)
   end if

   !-- Pre calculations has to be done before matrix initialisation
   call preCalc(zinterp)

   if ( l_3D ) then
      call initialize_der_arrays( l_rerror_fix, rerror_fac, r_3D)
   else
      call initialize_der_arrays( l_rerror_fix, rerror_fac)
   end if

   local_bytes_used = bytes_allocated
   if ( l_cheb_coll ) then
      if ( l_heat ) call initialize_temp_coll()
      if ( l_chem ) call initialize_xi_coll()
      if ( l_direct_solve ) then
         call initialize_om_coll_smat()
      else
         call initialize_om_coll_dmat()
      end if
   else
      if ( l_heat ) call initialize_temp_integ()
      if ( l_chem ) call initialize_xi_integ()
      if ( l_direct_solve ) then
         call initialize_psi_integ_smat()
      else
         call initialize_psi_integ_dmat()
      end if
   end if
   if ( l_heat_3D ) call initialize_update_temp_3D()
   if ( l_mag_3D ) call initialize_update_mag_3D()
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
   call time_loop(time, tscheme, run_init, zinterp)

   run_stop = MPI_Wtime()

   !--- Write stop time to SDTOUR and logfile:
   if ( rank == 0 ) then
      do n=1,2
         if ( n == 1 ) n_out=6
         if ( n == 2 )  n_out=n_log_file
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
   if ( l_mag_3D )  call finalize_update_mag_3D()
   if ( l_heat_3D )  call finalize_update_temp_3D()
   if ( l_cheb_coll ) then
      if ( l_heat ) call finalize_temp_coll()
      if ( l_chem ) call finalize_xi_coll()
      if ( l_direct_solve ) then
         call finalize_om_coll_smat()
      else
         call finalize_om_coll_dmat()
      end if
   else
      if ( l_heat ) call finalize_temp_integ()
      if ( l_chem ) call finalize_xi_integ()
      if ( l_direct_solve ) then
         call finalize_psi_integ_smat()
      else
         call finalize_psi_integ_dmat()
      end if
   end if
   call finalize_mfunctions()
   if ( l_3D ) call finalize_radial_loop_3D()
   call finalize_radial_loop()
   call finalize_fourier(l_3D)
   call finalize_radial_functions()
   call finalize_fieldsLast()
   call finalize_fields()
   if ( l_3D ) call finalize_blocking()
   call finalize_communications(l_3D)
   call destroy_mpi_domains(l_3D)
   call finalize_courant()
   if ( l_3D ) call finalize_outputs_3D()
   call finalize_outputs()
   if ( l_3D ) call finalize_der_arrays()
   call tscheme%finalize()
   if ( l_3D ) call zinterp%finalize()
   call finalize_truncation()

   !-- Initialize MPI
   call finalize_mpi()

end program pizza
