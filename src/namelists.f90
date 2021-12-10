module namelists
   !
   ! This module defines the different input namelists and their default
   ! values. This also controls the printing of the input values in the log.TAG
   ! file.
   !

   use iso_c_binding
   use iso_fortran_env, only: output_unit
   use truncation, only: n_r_max, m_max, n_cheb_max, minc
   use parallel_mod, only: rank
   use useful, only: abortRun
   use char_manip, only: length_to_blank, capitalize
   use constants, only: zero, one, half, two
   use precision_mod

   implicit none

   include 'fftw3.f03'

   private

   integer, parameter :: n_t_bounds=32 ! Number of parameters supported for t_top and t_bot

   character(len=72), public :: tag
   real(cp), public :: ra            ! Rayleigh number
   real(cp), public :: ek            ! Ekman number
   real(cp), public :: pr            ! Prandtl number
   real(cp), public :: radratio      ! Radius ratio
   real(cp), public :: raxi          ! Compisitional Rayleigh number
   real(cp), public :: sc            ! Schmidt number
   real(cp), public :: tcond_fac     ! Rescaling of the conducting temperature
   real(cp), public :: xicond_fac    ! Rescaling of the conducting composition
   real(cp), public :: beta_shift    ! Shift the upper bound of \beta
   real(cp), public :: beta_fac      ! Constant to modify the container (used for non-spherical containers)
   logical,  public :: l_finite_diff ! Shift to finite differences in radius
   logical,  public :: l_non_rot     ! Switch to do a non-rotatig annulus
   logical,  public :: l_temp_3D     ! 2D or 3D temperature background
   logical,  public :: l_xi_3D       ! 2D or 3D composition background
   logical,  public :: l_ek_pump     ! With or without Ekman pumping
   logical,  public :: l_temp_advz   ! With or without vertical advection of temp
   integer,  public :: ktopt, kbott  ! Temperature boundary condition
   integer,  public :: ktopxi, kbotxi! Boundary conditions for chemical composition
   integer,  public :: ktopv, kbotv  ! Velocity boundary condition
   real(cp), public :: t_bot(3*n_t_bounds), xi_bot(3*n_t_bounds)
   real(cp), public :: t_top(3*n_t_bounds), xi_top(3*n_t_bounds)
   real(cp), public :: g0, g1, g2
   real(cp), public :: h_temp, h_xi  ! Volumetric heating / volumetric chemical input

   !-- For the nonlinear mapping)
   real(cp), public :: alph1  ! Input parameter for non-linear map to define degree of spacing (0.0:2.0)
   real(cp), public :: alph2  ! Input parameter for non-linear map to define central point of different spacing (-1.0:1.0)
   character(len=72), public :: map_function ! Mapping family: either tangent or arcsin
   character(len=72), public :: time_scale  ! Time unit
   character(len=72), public :: time_scheme ! Time scheme
   character(len=72) :: bc_method    ! Galerkin or Tau-Lanczos method for BCs
   character(len=72), public :: container    ! shape of the container: sphere, expo, busse
   character(len=72) :: radial_scheme ! Chebyshev or finite differences
   character(len=72) :: cheb_method  ! Chebyshev method: collocation, integration
   character(len=72) :: matrix_solve ! Either direct or influence matrix
   character(len=72), public :: mpi_transp   ! 'AUTO', 'A2AV', 'A2AW'
   character(len=72) :: corio_term   ! Implicit or explicit treatment of Coriolis force
   character(len=72) :: buo_term    ! Implicit or explicit treatment of Buoyancy
   logical, public :: l_newmap      ! Switch for non-linear mapping (see Bayliss and Turkel, 1990)
   logical, public :: l_rerror_fix  ! Switch to fix round-off error in derivative calculation
   real(cp), public :: rerror_fac
   real(cp), public :: dtMin        ! Minimum allowed time step
   real(cp), public :: dtMax        ! Maximum allowed time step
   real(cp), public :: alpha        ! Weight for implicit time step
   real(cp), public :: courfac      ! Courant factor
   real(cp), public :: dt_fac       ! factor to control time step change
   real(cp), public :: tEND
   real(cp), public :: fd_stretch
   real(cp), public :: fd_ratio
   integer,  public :: n_time_steps
   integer :: runHours,runMinutes,runSeconds
   integer,  public :: run_time_requested
   integer :: n_fft_optim_lev ! FFTW planner flag
   integer(c_int), public :: fftw_plan_flag

   real(cp), public :: amp_t
   integer,  public :: init_t
   real(cp), public :: amp_xi
   integer,  public :: init_xi
   real(cp), public :: scale_t
   real(cp), public :: scale_xi
   real(cp), public :: amp_u
   integer,  public :: init_u
   real(cp), public :: scale_u
   logical,  public :: l_start_file     ! taking fields from startfile ?
   logical,  public :: l_reset_t ! Should we reset the time stored in the startfile?
   logical,  public :: l_packed_transp
   character(len=72), public :: start_file  ! name of start_file           

   integer,  public :: n_log_step
   integer,  public :: n_frames
   integer,  public :: n_frame_step
   integer,  public :: n_checkpoints
   integer,  public :: n_checkpoint_step
   integer,  public :: n_specs
   integer,  public :: n_spec_step
   logical,  public :: l_galerkin
   logical,  public :: l_vphi_balance    ! Calculate the vphi force balance
   logical,  public :: l_vort_balance    ! Calculate the vorticiy balance
   logical,  public :: l_2D_spectra      ! Calculate 2D spectra
   logical,  public :: l_2D_SD           ! Also store the standard deviation
   logical,  public :: l_corr            ! Calculte the correlation
   real(cp), public :: bl_cut            ! Cut-off boundary layers in the force balance
   logical,  public :: l_heat
   logical,  public :: l_chem
   logical,  public :: l_AB1
   logical,  public :: l_bridge_step
   logical,  public :: l_cheb_coll     ! Collocation method for Chebs
   logical,  public :: l_direct_solve  ! Direct solve or influence matrix method
   logical,  public :: l_coriolis_imp  ! Implicit treatment of Coriolis force
   logical,  public :: l_buo_imp       ! Implicit treatment of Buoyancy
   real(cp), public :: tadvz_fac
   real(cp), public :: r_cmb           ! Outer core radius
   real(cp), public :: r_icb           ! Inner core radius
   real(cp), public :: CorFac, BuoFac, TdiffFac, ViscFac, XiDiffFac, ChemFac
   real(cp), public :: hdif_temp       ! Hyperdiffusion amplitude on temperature
   real(cp), public :: hdif_comp       ! Hyperdiffusion amplitude on composition
   real(cp), public :: hdif_vel        ! Hyperdiffusion amplitude on velocity
   integer,  public :: hdif_exp        ! Exponent of the hyperdiffusion profile
   integer,  public :: hdif_m          ! Azimuthal wavenumber for hdif

   public :: read_namelists, write_namelists

contains

   subroutine read_namelists


      integer :: argument_count, input_handle, res
      character(len=100) :: input_filename, errmess
      logical :: nml_exist

      !-- Namelists:

      namelist/grid/n_r_max,n_cheb_max,m_max,minc,fd_ratio,fd_stretch
      namelist/control/tag,n_time_steps,alpha,l_newmap,map_function,&
      &                alph1,alph2,dtMax,courfac,tEND,runHours,     &
      &                runMinutes,runSeconds,l_non_rot,dt_fac,      &
      &                n_fft_optim_lev,time_scheme,cheb_method,     &
      &                l_rerror_fix, rerror_fac, time_scale,        &
      &                matrix_solve,corio_term,buo_term,bc_method,  &
      &                mpi_transp,l_packed_transp,radial_scheme
      namelist/hdif/hdif_temp,hdif_vel,hdif_exp,hdif_m,hdif_comp
      namelist/phys_param/ra,ek,pr,raxi,sc,radratio,g0,g1,g2,      &
      &                   ktopt,kbott,ktopv,kbotv,l_ek_pump,       &
      &                   l_temp_3D,tcond_fac,l_temp_advz,         &
      &                   beta_shift,ktopxi,kbotxi,t_bot,t_top,    &
      &                   xi_bot,xi_top, l_xi_3D, xicond_fac,      &
      &                   h_temp, h_xi, container, beta_fac
      namelist/start_field/l_start_file,start_file,scale_t,init_t,amp_t, &
      &                    scale_u,init_u,amp_u,l_reset_t,amp_xi,init_xi,&
      &                    scale_xi
      namelist/output_control/n_log_step,n_checkpoints, n_checkpoint_step, &
      &                       n_frames, n_frame_step, n_specs, n_spec_step,&
      &                       l_vphi_balance,l_vort_balance,bl_cut,        &
      &                       l_2D_spectra, l_2D_SD, l_corr

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
         read(input_handle,nml=grid,iostat=res,iomsg=errmess)
         if ( res > 0 .and. rank == 0 ) call abortRun(errmess)
         if ( res < 0 .and. rank == 0 ) then
            write(*,*) '!  No grid namelist found!'
         end if
         close(input_handle)

         open(newunit=input_handle,file=trim(input_filename))
         !-- Reading control parameters from namelists in STDIN:
         if ( rank == 0 ) write(*,*) '!  Reading control parameters!'
         read(input_handle,nml=control,iostat=res,iomsg=errmess)
         if ( res > 0 .and. rank == 0 ) call abortRun(errmess)
         if ( res < 0 .and. rank == 0 ) then
            write(*,*) '!  No control namelist found!'
         end if
         close(input_handle)

         open(newunit=input_handle,file=trim(input_filename))
         !-- Reading control parameters from namelists in STDIN:
         if ( rank == 0 ) write(*,*) '!  Reading hdif parameters!'
         read(input_handle,nml=hdif,iostat=res,iomsg=errmess)
         if ( res > 0 .and. rank == 0 ) call abortRun(errmess)
         if ( res < 0 .and. rank == 0 ) then
            write(*,*) '!  No hdif namelist found!'
         end if
         close(input_handle)

         open(newunit=input_handle,file=trim(input_filename))
         !-- Reading physical parameters from namelists in STDIN:
         if ( rank == 0 ) write(*,*) '!  Reading physical parameters!'
         read(input_handle,nml=phys_param,iostat=res,iomsg=errmess)
         if ( res > 0 .and. rank == 0 ) call abortRun(errmess)
         if ( res < 0 .and. rank == 0 ) then
            write(*,*) '!  No phys_param namelist found!'
         end if
         close(input_handle)

         open(newunit=input_handle,file=trim(input_filename))
         !-- Reading start field info from namelists in STDIN:
         if ( rank == 0 ) write(*,*) '!  Reading start information!'
         read(input_handle,nml=start_field,iostat=res,iomsg=errmess)
         if ( res > 0 .and. rank == 0 ) call abortRun(errmess)
         if ( res < 0 .and. rank == 0 ) then
            write(*,*) '! No start_field namelist found!'
         end if
         close(input_handle)

         open(newunit=input_handle,file=trim(input_filename))
         !-- Reading start field info from namelists in STDIN:
         if ( rank == 0 ) write(*,*) '!  Reading output control information!'
         read(input_handle,nml=output_control,iostat=res,iomsg=errmess)
         if ( res > 0 .and. rank == 0 ) call abortRun(errmess)
         if ( res < 0 .and. rank == 0 ) then
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
         l_chem=.true.
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

      r_cmb=one/(one-radratio)
      r_icb=r_cmb-one

      !-- Determine radial scheme
      call capitalize(radial_scheme)
      if ( index(radial_scheme, 'CHEB') /= 0 ) then
         l_finite_diff = .false.
      else
         l_finite_diff = .true.
      end if

      !-- Determine Cheb method
      call capitalize(cheb_method)
      if ( index(cheb_method, 'COLL') /= 0 ) then
         l_cheb_coll = .true.
      else
         l_cheb_coll = .false.
      end if

      !-- Determine BC method: Galerkin or Tau-Lanczos
      call capitalize(bc_method)
      if ( index(bc_method, 'GAL') /= 0 ) then
         l_galerkin = .true.
      else
         l_galerkin = .false.
      end if

      !-- Shape of the container (for QG computations)
      call capitalize(container)

      if ( index(container, 'SPHERE') == 0 .and. (.not. l_cheb_coll) .and. &
      &    (.not. l_finite_diff) ) then
         call abortRun('! This container is not supported by the sparse cheb form')
      end if

      !-- Determine solver method
      call capitalize(matrix_solve)
      if ( index(matrix_solve, 'DIRECT') /= 0 ) then
         l_direct_solve = .true.
      else
         l_direct_solve = .false.
      end if

      !-- Warning:: Galerkin method can not handle inhomogeneous BCs (yet)
      if ( l_galerkin .and. index(cheb_method, 'COLL') == 0 .and. &
      &    ( maxval(abs(t_top(:))) > 0.0_cp .or.                  &
      &      maxval(abs(xi_top(:))) > 0.0_cp) ) then
         call abortRun('! Inhomogeneous BCs not compatible with chosen chebyshev solver !')
      else if ( l_galerkin .and. index(cheb_method, 'COLL') == 0 .and. &
      &         ( maxval(abs(t_bot(:))) > 0.0_cp .or.                  &
      &           maxval(abs(xi_bot(:))) > 0.0_cp) ) then
         call abortRun('! Inhomogeneous BCs not compatible with chosen chebyshev solver !')
      end if

      !-- Implicit or Explicit treatment of Coriolis force
      call capitalize(corio_term)
      if ( index(corio_term, 'IMP') /= 0 ) then
         l_coriolis_imp = .true.
      else
         l_coriolis_imp = .false.
      end if

      !-- Implicit or Explicit treatment of Buoyancy
      call capitalize(buo_term)
      if ( index(buo_term, 'IMP') /= 0 ) then
         l_buo_imp = .true.
      else
         l_buo_imp = .false.
      end if

      if ( .not. l_direct_solve .and. l_coriolis_imp .and. (.not. l_non_rot) ) then
         l_coriolis_imp = .false.
         corio_term = 'EXP'
         if ( rank == 0 ) then
            write(output_unit, &
            &    '(" ! Implicit Coriolis term not compatible with influence matrix method")')
            write(output_unit, &
            &    '(" ! Coriolis term will be treated explicitly")')
         end if
      end if
      
      if ( l_finite_diff ) then
         l_buo_imp = .false.
         if ( rank == 0 ) then
            write(output_unit, &
            &    '(" ! Implicit buoyancy term not compatible with finite differences")')
            write(output_unit, &
            &    '(" ! Buoyancy term will be treated explicitly")')
         end if
      end if

      !-- Time unit
      call capitalize(time_scale) 
      if ( l_non_rot ) then
         CorFac = 0.0_cp
         if ( l_heat ) then
            BuoFac = ra/pr
            TdiffFac = one/pr
         else
            BuoFac = 0.0_cp
            TdiffFac = 0.0_cp
         end if
         if ( l_chem ) then
            ChemFac = raxi/sc
            XidiffFac = one/sc
         else
            ChemFac = 0.0_cp
            XidiffFac = 0.0_cp
         end if
         ViscFac = one
      else
         if ( index(time_scale, 'ROT') /= 0 ) then
            CorFac = two
            if ( l_heat ) then
               BuoFac = ra*ek*ek/pr
               TdiffFac = ek/pr
            else
               BuoFac = 0.0_cp
               TdiffFac = 0.0_cp
            end if
            if ( l_chem ) then
               ChemFac = raxi*ek*ek/sc
               XidiffFac = ek/sc
            else
               ChemFac = 0.0_cp
               XidiffFac = 0.0_cp
            end if
            ViscFac = ek
         else
            CorFac = two/ek
            if ( l_heat ) then
               BuoFac = ra/pr
               TdiffFac = one/pr
            else
               BuoFac = 0.0_cp
               TdiffFac = 0.0_cp
            end if
            if ( l_chem ) then
               ChemFac = raxi/sc
               XidiffFac = one/sc
            else
               ChemFac = 0.0_cp
               XidiffFac = 0.0_cp
            end if
            ViscFac = one
         end if
      end if

   end subroutine read_namelists
!--------------------------------------------------------------------------------
   subroutine default_namelists
      !
      !  Purpose of this subroutine is to set default parameters          
      !  for the namelists.                                               
      !

      !-- Local variable:
      integer :: n

      !-- &grid namelist
      n_r_max          =32
      n_cheb_max       =32
      m_max            =32
      minc             =1
      fd_stretch       =0.3_cp
      fd_ratio         =0.1_cp

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
      courfac          =1.0e3_cp
      dt_fac           =2.0_cp
      tEND             =0.0_cp    ! numerical time where run should end
      runHours         =0
      runMinutes       =10
      runSeconds       =0
      n_fft_optim_lev  =1 ! Optimisation level for FFT:
                          ! 0: FFTW_ESTIMATE, 1: FFTW_MEASURE, 2: FFTW_PATIENT
                          ! 3: FFTW_EXHAUSTIVE
      time_scheme      ='CNAB2'
      radial_scheme    ='cheb'
      cheb_method      ='colloc'
      bc_method        ='tau-lanczos'
      matrix_solve     ='DIRECT'
      l_rerror_fix     =.true.
      rerror_fac       =500.0_cp
      corio_term       ='IMPLICIT' ! Implicit treatment of Coriolis term
      buo_term         ='IMPLICIT' ! Implicit treatment of Buoyancy
      time_scale       ='VISC' ! viscous units
      mpi_transp       ='AUTO'
      l_packed_transp  =.false.

      !-- Hyperdiffusion
      hdif_vel         =0.0_cp
      hdif_temp        =0.0_cp
      hdif_comp        =0.0_cp
      hdif_m           =0
      hdif_exp         =0

      !-- Physcal parameters
      l_non_rot        =.false.
      l_ek_pump        =.false.
      l_temp_advz      =.false.
      l_temp_3D        =.false.
      l_xi_3D          =.false.
      ra               =1.0e5_cp
      pr               =one
      ek               =1.0e-3_cp
      raxi             =0.0e5_cp
      sc               =0.0_cp
      radratio         =0.35_cp
      tcond_fac        =one
      xicond_fac       =one
      beta_shift       =0.0_cp
      beta_fac         =0.4_cp
      !----- Gravity parameters: defaut value g propto r (i.e. g1=1)
      g0               =0.0_cp
      g1               =one
      g2               =0.0_cp
      !----- Volumetric sources
      h_temp           =0.0_cp
      h_xi             =0.0_cp
      !----- Boundary conditions        
      ktopt            =1
      kbott            =1
      ktopxi           =1
      kbotxi           =1
      ktopv            =2
      kbotv            =2
      !----- Parameters for temp. BCs
      do n=1,3*n_t_bounds
         t_bot(n) =0.0_cp
         t_top(n) =0.0_cp
         xi_bot(n)=0.0_cp
         xi_top(n)=0.0_cp
      end do
      container        ='sphere'

      !----- Namelist start_field:
      l_reset_t        =.false.
      l_start_file     =.false.
      start_file       ="no_start_file"
      init_t           =0
      scale_t          =1.0_cp
      amp_t            =0.0_cp
      init_xi          =0
      scale_xi         =1.0_cp
      amp_xi           =0.0_cp
      init_u           =0
      scale_u          =1.0_cp
      amp_u            =0.0_cp

      !----- Output namelist
      n_log_step       =50
      n_frames         =0
      n_frame_step     =0
      n_checkpoints    =1
      n_checkpoint_step=0
      n_specs          =0
      n_spec_step      =0
      l_vphi_balance   =.false.
      l_vort_balance   =.false.
      l_2D_spectra     =.false.
      l_2D_SD          =.false.
      l_corr           =.false.
      bl_cut           =1.0e-3_cp

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
      write(n_out,'(''  fd_stretch      ='',ES14.6,'','')') fd_stretch
      write(n_out,'(''  fd_ratio        ='',ES14.6,'','')') fd_ratio
      write(n_out,*) "/"

      write(n_out,*) "&control"
      length=length_to_blank(tag)
      write(n_out,*) " tag             = """,tag(1:length),""","
      write(n_out,'(''  n_time_steps    ='',i8,'','')') n_time_steps
      length=length_to_blank(time_scheme)
      write(n_out,*) " time_scheme     = """,time_scheme(1:length),""","
      length=length_to_blank(radial_scheme)
      write(n_out,*) " radial_scheme   = """,radial_scheme(1:length),""","
      length=length_to_blank(cheb_method)
      write(n_out,*) " cheb_method     = """,cheb_method(1:length),""","
      length=length_to_blank(bc_method)
      write(n_out,*) " bc_method       = """,bc_method(1:length),""","
      length=length_to_blank(matrix_solve)
      write(n_out,*) " matrix_solve    = """,matrix_solve(1:length),""","
      length=length_to_blank(mpi_transp)
      write(n_out,*) " mpi_transp      = """,mpi_transp(1:length),""","
      write(n_out,'(''  l_packed_transp ='',l3,'','')') l_packed_transp
      length=length_to_blank(corio_term)
      write(n_out,*) " corio_term      = """,corio_term(1:length),""","
      length=length_to_blank(buo_term)
      write(n_out,*) " buo_term        = """,buo_term(1:length),""","
      write(n_out,'(''  alpha           ='',ES14.6,'','')')   alpha
      write(n_out,'(''  l_newmap        ='',l3,'','')') l_newmap
      length=length_to_blank(map_function)
      write(n_out,*) " map_function    = """,map_function(1:length),""","
      write(n_out,'(''  alph1           ='',ES14.6,'','')') alph1
      write(n_out,'(''  alph2           ='',ES14.6,'','')') alph2
      write(n_out,'(''  dtMax           ='',ES14.6,'','')') dtMax
      write(n_out,'(''  dt_fac          ='',ES14.6,'','')') dt_fac
      write(n_out,'(''  runHours        ='',i4,'','')') runHours
      write(n_out,'(''  runMinutes      ='',i4,'','')') runMinutes
      write(n_out,'(''  runSeconds      ='',i4,'','')') runSeconds
      write(n_out,'(''  tEND            ='',ES14.6,'','')') tEND
      write(n_out,'(''  l_non_rot       ='',l3,'','')') l_non_rot
      write(n_out,'(''  l_rerror_fix    ='',l3,'','')') l_rerror_fix
      write(n_out,'(''  rerror_fac      ='',ES14.6,'','')') rerror_fac
      write(n_out,'(''  n_fft_optim_lev ='',i4,'','')') n_fft_optim_lev
      length=length_to_blank(time_scale)
      write(n_out,*) " time_scale      = """,time_scale(1:length),""","
      write(n_out,*) "/"

      write(n_out,*) "&hdif"
      length=length_to_blank(tag)
      write(n_out,'(''  hdif_vel        ='',ES14.6,'','')') hdif_vel
      write(n_out,'(''  hdif_temp       ='',ES14.6,'','')') hdif_temp
      write(n_out,'(''  hdif_comp       ='',ES14.6,'','')') hdif_comp
      write(n_out,'(''  hdif_exp        ='',i4,'','')') hdif_exp
      write(n_out,'(''  hdif_m          ='',i4,'','')') hdif_m
      write(n_out,*) "/"

      write(n_out,*) "&phys_param"
      write(n_out,'(''  ra              ='',ES14.6,'','')') ra
      write(n_out,'(''  pr              ='',ES14.6,'','')') pr
      write(n_out,'(''  ek              ='',ES14.6,'','')') ek
      write(n_out,'(''  raxi            ='',ES14.6,'','')') raxi
      write(n_out,'(''  sc              ='',ES14.6,'','')') sc
      write(n_out,'(''  radratio        ='',ES14.6,'','')') radratio
      write(n_out,'(''  tcond_fac       ='',ES14.6,'','')') tcond_fac
      write(n_out,'(''  xicond_fac      ='',ES14.6,'','')') xicond_fac
      length=length_to_blank(container)
      write(n_out,*) " container       = """,container(1:length),""","
      write(n_out,'(''  beta_shift      ='',ES14.6,'','')') beta_shift
      write(n_out,'(''  beta_fac        ='',ES14.6,'','')') beta_fac
      write(n_out,'(''  g0              ='',ES14.6,'','')') g0
      write(n_out,'(''  g1              ='',ES14.6,'','')') g1
      write(n_out,'(''  g2              ='',ES14.6,'','')') g2
      write(n_out,'(''  h_temp          ='',ES14.6,'','')') h_temp
      write(n_out,'(''  h_xi            ='',ES14.6,'','')') h_xi
      write(n_out,'(''  l_ek_pump       ='',l3,'','')') l_ek_pump
      write(n_out,'(''  l_temp_advz     ='',l3,'','')') l_temp_advz
      write(n_out,'(''  l_temp_3D       ='',l3,'','')') l_temp_3D
      write(n_out,'(''  l_xi_3D         ='',l3,'','')') l_xi_3D
      !--- Heat boundary condition:
      write(n_out,'(''  ktopt           ='',i3,'','')') ktopt
      write(n_out,'(''  kbott           ='',i3,'','')') kbott
      write(n_out,'(''  ktopxi          ='',i3,'','')') ktopxi
      write(n_out,'(''  kbotxi          ='',i3,'','')') kbotxi
      write(n_out,'(''  ktopv           ='',i3,'','')') ktopv
      write(n_out,'(''  kbotv           ='',i3,'','')') kbotv
      write(n_out,*) "/"

      write(n_out,*) "&start_field"
      write(n_out,'(''  l_start_file    ='',l3,'','')') l_start_file
      write(n_out,'(''  l_reset_t       ='',l3,'','')') l_reset_t
      length=length_to_blank(start_file)
      write(n_out,*) " start_file      = """,start_file(1:length),""","
      write(n_out,'(''  scale_t         ='',ES14.6,'','')') scale_t
      write(n_out,'(''  init_t          ='',i7,'','')') init_t
      write(n_out,'(''  amp_t           ='',ES14.6,'','')') amp_t
      write(n_out,'(''  scale_xi        ='',ES14.6,'','')') scale_xi
      write(n_out,'(''  init_xi         ='',i7,'','')') init_xi
      write(n_out,'(''  amp_xi          ='',ES14.6,'','')') amp_xi
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
      write(n_out,'(''  l_vort_balance  ='',l3,'','')') l_vort_balance
      write(n_out,'(''  bl_cut          ='',ES14.6,'','')') bl_cut
      write(n_out,'(''  l_2D_spectra    ='',l3,'','')') l_2D_spectra
      write(n_out,'(''  l_2D_SD         ='',l3,'','')') l_2D_SD
      write(n_out,'(''  l_corr          ='',l3,'','')') l_corr
      write(n_out,*) "/"

   end subroutine write_namelists
!--------------------------------------------------------------------------------
end module namelists
