program pizza

   use parallel_mod
   use precision_mod
   use mem_alloc
   use step_time, only: time_loop
   use radial_der, only: initialize_der_arrays, finalize_der_arrays
   use init_fields, only: get_start_fields
   use fields, only: initialize_fields, finalize_fields
   use fieldsLast, only: initialize_fieldsLast, finalize_fieldsLast
   use communications, only: initialize_communications, finalize_communications
   use blocking, only: set_mpi_domains, nMstart, nMstop, destroy_mpi_domains
   use namelists, only: read_namelists, write_namelists, tag, imp_scheme, &
       &                exp_scheme
   use outputs, only: initialize_outputs, finalize_outputs, n_log_file
   use pre_calculations, only: preCalc
   use radial_functions, only: initialize_radial_functions, &
       &                       finalize_radial_functions
   use truncation, only: initialize_truncation, n_r_max, n_phi_max, & 
       &                 finalize_truncation
   use fourier, only: initialize_fourier, finalize_fourier
   use rloop, only: initialize_radial_loop, finalize_radial_loop
   use update_temperature, only: initialize_update_temp, finalize_update_temp
   use update_psi, only: initialize_update_om, finalize_update_om
   use time_scheme, only: type_tscheme
   use useful, only: formatTime

   implicit none

   real(cp) :: time
   type(type_tscheme) :: tscheme
   real(cp) :: runStop, runStart
   integer(lip) :: local_bytes_used
   integer :: values(8)
   integer :: n, n_out
   character(len=72) :: date

   !-- Initialize MPI
   call initialize_mpi()

   runStart = MPI_Wtime()

   !--
   if ( rank == 0 ) then
      write(*,*)
      write(*,*) '!--- Program pizza  ---!'
      call date_and_time(values=values)
      write(date, '(i4,''/'',i0.2,''/'',i0.2,'' '', i0.2,'':'',i0.2,'':'',i0.2)') &
      &     values(1), values(2), values(3), values(5), values(6), values(7)
      write(6, *) '!  Start time:  ', date

   end if

   !-- Read input parameters
   call read_namelists()

   !-- Init memory counter
   call initialize_memory_counter(tag)

   !-- Initialize truncation
   call initialize_truncation()

   !-- Set the domain decomposition
   call set_mpi_domains()

   !-- Open output files
   call initialize_outputs()

   !-- Initialize time scheme
   call tscheme%initialize(imp_scheme, exp_scheme)

   !-- Initialize MPI communicators
   call initialize_communications()

   local_bytes_used = bytes_allocated
   call initialize_fields()
   call initialize_fieldsLast(tscheme%norder_exp, tscheme%norder_imp)
   local_bytes_used = bytes_allocated-local_bytes_used
   call memWrite('Fields', local_bytes_used)
   call initialize_radial_functions()
   call initialize_der_arrays(n_r_max, nMstart, nMstop)
   local_bytes_used = bytes_allocated
   call initialize_fourier(n_phi_max)
   call initialize_radial_loop(n_phi_max)
   call memWrite('R loop', local_bytes_used)
   local_bytes_used = bytes_allocated-local_bytes_used
   local_bytes_used = bytes_allocated
   call initialize_update_om()
   call initialize_update_temp()
   local_bytes_used = bytes_allocated-local_bytes_used
   call memWrite('M loop', local_bytes_used)

   call finalize_memory_counter()

   if ( rank == 0 ) then
      call write_namelists(6)
      call write_namelists(n_log_file)
   end if

   !-- Pre calculations
   call preCalc()

   !-- Start fields
   call get_start_fields(time, tscheme)

   !--- Write starting time to SDTOUT and logfile:
   if ( rank == 0 ) then
      do n=1,2
         if ( n == 1 ) n_out=6
         if ( n == 2 ) n_out=n_log_file
         write(n_out,'(/,'' ! Starting time integration at:'')')
         write(n_out,'(''   start_time ='',1p,ES18.10)') time
         write(n_out,'(''   start dt   ='',1p,ES16.4)') tscheme%dt(1)
      end do
   end if

   !-- Time integration
   call time_loop(time, tscheme)

   runStop = MPI_Wtime()

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
         call formatTime(n_out,'! Total run time:', runStop-runStart)
         write(n_out,*) '!-----------------------------------!'
         write(n_out,*) '!---- Ready for another pizza ? ----!'
         write(n_out,*) '!-----------------------------------!'
         write(n_out,*)
      end do
   end if

   !-- Close files
   call finalize_update_temp()
   call finalize_update_om()
   call finalize_radial_loop()
   call finalize_fourier()
   call finalize_der_arrays()
   call finalize_radial_functions()
   call finalize_fieldsLast()
   call finalize_fields()
   call finalize_communications()
   call destroy_mpi_domains()
   call finalize_outputs()
   call finalize_truncation()

   !-- Initialize MPI
   call finalize_mpi()

end program pizza
