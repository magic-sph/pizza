module namelists

   use iso_c_binding
   use truncation, only: n_r_max, m_max, n_cheb_max, minc
   use parallel_mod, only: rank
   use useful, only: abortRun
   use char_manip, only: length_to_blank, capitalize
   use constants, only: one, half
   use precision_mod

   implicit none

   include 'fftw3.f03'

   private

   character(len=72), public :: tag
   real(cp), public :: ra            ! Rayleigh number
   real(cp), public :: ek            ! Ekman number
   real(cp), public :: pr            ! Prandtl number
   real(cp), public :: radratio      ! Radius ratio
   real(cp), public :: raxi          ! Compisitional Rayleigh number
   real(cp), public :: sc            ! Schmidt number
   real(cp), public :: tcond_fac     ! Rescaling of the conducting temperature
   logical,  public :: l_non_rot     ! Switch to do a non-rotatig annulus
   logical,  public :: l_temp_3D     ! 2D or 3D temperature background
   logical,  public :: l_ek_pump     ! With or without Ekman pumping
   logical,  public :: l_temp_advz   ! With or without vertical advection of temp
   integer,  public :: ktopt, kbott  ! Temperature boundary condition
   integer,  public :: ktopv, kbotv  ! Velocity boundary condition
   real(cp), public :: g0, g1, g2

   !-- For the nonlinear mapping)
   real(cp), public :: alph1  ! Input parameter for non-linear map to define degree of spacing (0.0:2.0)
   real(cp), public :: alph2  ! Input parameter for non-linear map to define central point of different spacing (-1.0:1.0)
   character(len=72), public :: map_function ! Mapping family: either tangent or arcsin
   character(len=72), public :: imp_scheme ! Implicit scheme
   character(len=72), public :: exp_scheme ! Implicit scheme
   logical, public :: l_newmap       ! Switch for non-linear mapping (see Bayliss and Turkel, 1990)
   real(cp), public :: dtMin           ! Minimum allowed time step
   real(cp), public :: dtMax           ! Maximum allowed time step
   real(cp), public :: alpha           ! Weight for implicit time step
   real(cp), public :: courfac   ! Courant factor
   real(cp), public :: tEND
   integer,  public :: n_time_steps
   integer :: runHours,runMinutes,runSeconds
   integer,  public :: run_time_requested
   integer :: n_fft_optim_lev ! FFTW planner flag
   integer(c_int), public :: fftw_plan_flag

   real(cp), public :: amp_t
   integer,  public :: init_t
   real(cp), public :: scale_t
   real(cp), public :: amp_u
   integer,  public :: init_u
   real(cp), public :: scale_u
   logical,  public :: l_start_file     ! taking fields from startfile ?
   character(len=72), public :: start_file  ! name of start_file           

   integer,  public :: n_log_step
   integer,  public :: n_frames
   integer,  public :: n_frame_step
   integer,  public :: n_checkpoints
   integer,  public :: n_checkpoint_step
   integer,  public :: n_specs
   integer,  public :: n_spec_step
   logical,  public :: l_vphi_balance    ! Calculate the vphi force balance
   logical,  public :: l_heat
   logical,  public :: l_chem
   logical,  public :: l_AB1
   logical,  public :: l_bridge_step
   real(cp), public :: tadvz_fac

   public :: read_namelists, write_namelists

contains

   subroutine read_namelists


      integer :: argument_count, input_handle, res
      character(len=100) :: input_filename
      logical :: nml_exist

      !-- Namelists:

      namelist/grid/n_r_max,n_cheb_max,m_max,minc
      namelist/control/tag,n_time_steps,alpha,l_newmap,map_function,&
      &                alph1,alph2,dtMax,courfac,tEnd,runHours,     &
      &                runMinutes,runSeconds,l_non_rot,             &
      &                n_fft_optim_lev,imp_scheme,exp_scheme
      namelist/phys_param/ra,ek,pr,raxi,sc,radratio,g0,g1,g2,  &
      &                   ktopt,kbott,ktopv,kbotv,l_ek_pump,   &
      &                   l_temp_3D,tcond_fac,l_temp_advz
      namelist/start_field/l_start_file,start_file,scale_t,init_t,amp_t, &
      &                    scale_u,init_u,amp_u
      namelist/output_control/n_log_step,n_checkpoints, n_checkpoint_step, &
      &                       n_frames, n_frame_step, n_specs, n_spec_step,&
      &                       l_vphi_balance

   !namelist/control/tag,n_times

      !-- Set default values of control parameters:
      call default_namelists()

      ! get the filename of the input file as first argument from the command line
      argument_count = command_argument_count()
      if (argument_count == 0) then
         call abortRun('The filename of the input file must be provided as first argument')
      else
         call get_command_argument(1,input_filename)

         inquire(file = input_filename, exist = nml_exist)

         if (.not. nml_exist) then
            call abortRun('! Input namelist file not found!')
         end if

         open(newunit=input_handle,file=trim(input_filename))

         if ( rank == 0 ) write(*,*) '!  Reading grid parameters!'
         read(input_handle,nml=grid,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(*,*) '! No grid namelist found!'
         end if
         close(input_handle)

         open(newunit=input_handle,file=trim(input_filename))
         !-- Reading control parameters from namelists in STDIN:
         if ( rank == 0 ) write(*,*) '!  Reading control parameters!'
         read(input_handle,nml=control,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(*,*) '! No control namelist found!'
         end if
         close(input_handle)

         open(newunit=input_handle,file=trim(input_filename))
         !-- Reading physical parameters from namelists in STDIN:
         if ( rank == 0 ) write(*,*) '!  Reading physical parameters!'
         read(input_handle,nml=phys_param,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(*,*) '! No phys_param namelist found!'
         end if
         close(input_handle)

         open(newunit=input_handle,file=trim(input_filename))
         !-- Reading start field info from namelists in STDIN:
         if ( rank == 0 ) write(*,*) '!  Reading start information!'
         read(input_handle,nml=start_field,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(*,*) '! No start_field namelist found!'
         end if
         close(input_handle)

         open(newunit=input_handle,file=trim(input_filename))
         !-- Reading start field info from namelists in STDIN:
         if ( rank == 0 ) write(*,*) '!  Reading output control information!'
         read(input_handle,nml=output_control,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(*,*) '! No output control namelist found!'
         end if
         close(input_handle)

      end if


      !--- Stuff for the radial non-linear mapping
      call capitalize(map_function)

      !-- This is the case of a tangent mapping (see Bayliss & Turkel, 1992)
      if ( index(map_function, 'TAN') /= 0 .or. index(map_function, 'BAY') /= 0 ) then
         !     alph1 can be any positive number, above 0
         !     alph2 has to be a value between -1 and 1 (interval in Chebyshev space)
         if ( (alph1 == 0.0_cp) .or. (alph2 < -one) .or. (alph2 > one) ) then
            call abortRun('! Chebyshev mapping parameter is not correct !')
         elseif ( l_newmap .and. (alph1 < 0.0_cp) ) then
            alph1=abs(alph1)
         end if
      !-- This is the case of the Kosloff & Tal-Ezer mapping (1993)
      else if ( index(map_function, 'ARCSIN') /= 0 .or. &
      &         index(map_function, 'KTL') /= 0 ) then
         if ( (alph1 < 0.0_cp) .or. (alph1 >= one) ) then
            call abortRun('! Chebyshev mapping parameter is not correct !')
         end if
      end if

      run_time_requested = runHours*3600+runMinutes*60+runSeconds

      if ( ek < 0.0_cp ) l_non_rot=.true.
      if ( l_non_rot ) then
         ek=-one ! used as a flag, not used for the calculation
      end if

      if ( ra == 0.0_cp ) then
         l_heat=.false.
      else
         l_heat=.true.
      end if

      if ( raxi == 0.0_cp ) then
         l_chem=.false.
      else
         l_heat=.true.
      end if

      !-- Determine whether we also consider vertical advection of temperature
      if ( l_temp_advz ) then
         tadvz_fac = one
      else
         tadvz_fac = 0.0_cp
      end if

      if ( n_fft_optim_lev == 0 ) then
         fftw_plan_flag = FFTW_ESTIMATE
      else if ( n_fft_optim_lev == 1 ) then
         fftw_plan_flag = FFTW_MEASURE
      else if ( n_fft_optim_lev == 2 ) then
         fftw_plan_flag = FFTW_PATIENT
      else if ( n_fft_optim_lev == 3 ) then
         fftw_plan_flag = FFTW_EXHAUSTIVE
      else
         fftw_plan_flag = FFTW_ESTIMATE
      end if

      !-- This is used when 1st order AB is used (i.e. Euler)
      l_AB1 = .false.
      l_bridge_step = .true.

   end subroutine read_namelists
!--------------------------------------------------------------------------------
   subroutine default_namelists

      !-- &grid namelist
      n_r_max          =32
      n_cheb_max       =32
      m_max            =32
      minc             =1

      !-- Control namelist
      n_time_steps     =100
      tag              ='test'
      !----- Non-linear mapping parameters (Bayliss, 1992):
      alpha            =half
      l_newmap         =.false.
      alph1            =0.8_cp
      alph2            =0.0_cp
      map_function     ='arcsin' ! By default Kosloff and Tal-Ezer mapping when l_newmap=.true.
      dtMax            =1.0e-5_cp
      courfac          =2.5_cp
      tEND             =0.0_cp    ! numerical time where run should end
      runHours         =0
      runMinutes       =10
      runSeconds       =0
      n_fft_optim_lev  =1 ! Optimisation level for FFT:
                          ! 0: FFTW_ESTIMATE, 1: FFTW_MEASURE, 2: FFTW_PATIENT
                          ! 3: FFTW_EXHAUSTIVE
      imp_scheme       ='CN'
      exp_scheme       ='AB2'

      !-- Physcal parameters
      l_non_rot        =.false.
      l_ek_pump        =.false.
      l_temp_advz      =.false.
      l_temp_3D        =.false.
      ra               =1.0e5_cp
      pr               =one
      ek               =1.0e-3_cp
      raxi             =0.0e5_cp
      sc               =0.0_cp
      radratio         =0.35_cp
      tcond_fac        =one
      !----- Gravity parameters: defaut value g propto r (i.e. g1=1)
      g0               =0.0_cp
      g1               =one
      g2               =0.0_cp
      !----- Boundary conditions        
      ktopt            =1
      kbott            =1
      ktopv            =2
      kbotv            =2

      !----- Namelist start_field:
      l_start_file     =.false.
      start_file       ="no_start_file"
      init_t           =0
      scale_t          =1.0_cp
      amp_t            =0.0_cp
      init_u           =0
      scale_u          =1.0_cp
      amp_u            =0.0_cp

      !----- Output namelist
      n_log_step       =50
      n_frames         =1
      n_frame_step     =0
      n_checkpoints    =1
      n_checkpoint_step=0
      n_specs          =0
      n_spec_step      =0
      l_vphi_balance   =.false.

   end subroutine default_namelists
!--------------------------------------------------------------------------------
   subroutine write_namelists(n_out)

      !-- Input file unit
      integer, intent(in) :: n_out

      !-- Local variables
      integer :: length

      write(n_out,*) " "
      write(n_out,*) "&grid"
      write(n_out,'(''  n_r_max         ='',i5,'','')') n_r_max
      write(n_out,'(''  n_cheb_max      ='',i5,'','')') n_cheb_max
      write(n_out,'(''  m_max           ='',i5,'','')') m_max
      write(n_out,'(''  minc            ='',i5,'','')') minc
      write(n_out,*) "/"

      write(n_out,*) "&control"
      length=length_to_blank(tag)
      write(n_out,*) " tag             = """,tag(1:length),""","
      write(n_out,'(''  n_time_steps    ='',i8,'','')') n_time_steps
      length=length_to_blank(imp_scheme)
      write(n_out,*) " imp_scheme      = """,imp_scheme(1:length),""","
      length=length_to_blank(exp_scheme)
      write(n_out,*) " exp_scheme      = """,exp_scheme(1:length),""","
      write(n_out,'(''  alpha           ='',ES14.6,'','')')   alpha
      write(n_out,'(''  l_newmap        ='',l3,'','')') l_newmap
      length=length_to_blank(map_function)
      write(n_out,*) " map_function    = """,map_function(1:length),""","
      write(n_out,'(''  alph1           ='',ES14.6,'','')') alph1
      write(n_out,'(''  alph2           ='',ES14.6,'','')') alph2
      write(n_out,'(''  dtMax           ='',ES14.6,'','')') dtMax
      write(n_out,'(''  courfac         ='',ES14.6,'','')') courfac
      write(n_out,'(''  runHours        ='',i4,'','')') runHours
      write(n_out,'(''  runMinutes      ='',i4,'','')') runMinutes
      write(n_out,'(''  runSeconds      ='',i4,'','')') runSeconds
      write(n_out,'(''  tEND            ='',ES14.6,'','')') tEND
      write(n_out,'(''  l_non_rot       ='',l3,'','')') l_non_rot
      write(n_out,'(''  n_fft_optim_lev ='',i4,'','')') n_fft_optim_lev
      write(n_out,*) "/"

      write(n_out,*) "&phys_param"
      write(n_out,'(''  ra              ='',ES14.6,'','')') ra
      write(n_out,'(''  pr              ='',ES14.6,'','')') pr
      write(n_out,'(''  ek              ='',ES14.6,'','')') ek
      write(n_out,'(''  raxi            ='',ES14.6,'','')') raxi
      write(n_out,'(''  sc              ='',ES14.6,'','')') sc
      write(n_out,'(''  radratio        ='',ES14.6,'','')') radratio
      write(n_out,'(''  tcond_fac       ='',ES14.6,'','')') tcond_fac
      write(n_out,'(''  g0              ='',ES14.6,'','')') g0
      write(n_out,'(''  g1              ='',ES14.6,'','')') g1
      write(n_out,'(''  g2              ='',ES14.6,'','')') g2
      write(n_out,'(''  l_ek_pump       ='',l3,'','')') l_ek_pump
      write(n_out,'(''  l_temp_advz     ='',l3,'','')') l_temp_advz
      write(n_out,'(''  l_temp_3D       ='',l3,'','')') l_temp_3D
      !--- Heat boundary condition:
      write(n_out,'(''  ktopt           ='',i3,'','')') ktopt
      write(n_out,'(''  kbott           ='',i3,'','')') kbott
      write(n_out,'(''  ktopv           ='',i3,'','')') ktopv
      write(n_out,'(''  kbotv           ='',i3,'','')') kbotv
      write(n_out,*) "/"

      write(n_out,*) "&start_field"
      write(n_out,'(''  l_start_file    ='',l3,'','')') l_start_file
      length=length_to_blank(start_file)
      write(n_out,*) " start_file      = """,start_file(1:length),""","
      write(n_out,'(''  scale_t         ='',ES14.6,'','')') scale_t
      write(n_out,'(''  init_t          ='',i7,'','')') init_t
      write(n_out,'(''  amp_t           ='',ES14.6,'','')') amp_t
      write(n_out,'(''  scale_u         ='',ES14.6,'','')') scale_u
      write(n_out,'(''  init_u          ='',i7,'','')') init_u
      write(n_out,'(''  amp_u           ='',ES14.6,'','')') amp_u
      write(n_out,*) "/"

      write(n_out,*) "&output_control"
      write(n_out,'(''  n_log_step      ='',i5,'','')') n_log_step
      write(n_out,'(''  n_checkpoints   ='',i5,'','')') n_checkpoints
      write(n_out,'('' n_checkpoint_step='',i5,'','')') n_checkpoint_step
      write(n_out,'(''  n_frames        ='',i5,'','')') n_frames
      write(n_out,'(''  n_frame_step    ='',i5,'','')') n_frame_step
      write(n_out,'(''  n_specs         ='',i5,'','')') n_specs
      write(n_out,'(''  n_spec_step     ='',i5,'','')') n_spec_step
      write(n_out,'(''  l_vphi_balance  ='',l3,'','')') l_vphi_balance
      write(n_out,*) "/"

   end subroutine write_namelists
!--------------------------------------------------------------------------------
end module namelists
