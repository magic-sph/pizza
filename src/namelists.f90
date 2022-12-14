module namelists
   !
   ! This module defines the different input namelists and their default
   ! values. This also controls the printing of the input values in the log.TAG
   ! file.
   !

   use iso_c_binding
   use iso_fortran_env, only: output_unit
   use truncation, only: n_r_max, m_max, n_cheb_max, minc
   use truncation_3D, only: n_r_max_3D, n_cheb_max_3D, minc_3D, n_phi_tot_3D
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
   real(cp), public :: prmag         ! Magnetic Prandtl number
   real(cp), public :: radratio      ! Radius ratio
   real(cp), public :: raxi          ! Compositional Rayleigh number
   real(cp), public :: sc            ! Schmidt number
   real(cp), public :: alpha_fac     ! Alpha effect Amplitude control parameter
   real(cp), public :: mag_pump_fac  ! Magnetic pumping Amplitude control parameter
   real(cp), public :: delta_fac     ! Alpha^2 effect Amplitude control parameter (Davidson)
   real(cp), public :: epsrc0        ! Internal heat sources
   real(cp), public :: epsrc0xi      ! Internal chemical sources
   real(cp), public :: tcond_fac     ! Rescaling of the conducting temperature
   real(cp), public :: xicond_fac    ! Rescaling of the conducting composition
   real(cp), public :: beta_shift    ! Shift the upper bound of \beta
   logical,  public :: l_non_rot     ! Switch to do a non-rotatig annulus
   logical,  public :: l_tcond_3D    ! 2D or 3D temperature background
   logical,  public :: l_xi_3D       ! 2D or 3D composition background
   logical,  public :: l_ek_pump     ! With or without Ekman pumping
   logical,  public :: l_temp_advz   ! With or without vertical advection of temp
   integer,  public :: ktopt, kbott  ! Temperature boundary condition
   integer,  public :: ktopxi, kbotxi! Boundary conditions for chemical composition
   integer,  public :: ktopv, kbotv  ! Velocity boundary condition
   integer,  public :: ktopb, kbotb  ! Magnetic field boundary condition
   real(cp), public :: t_bot(3*n_t_bounds), xi_bot(3*n_t_bounds)
   real(cp), public :: t_top(3*n_t_bounds), xi_top(3*n_t_bounds)
   real(cp), public :: g0, g1, g2

   !-- For the nonlinear mapping)
   real(cp), public :: alph1  ! Input parameter for non-linear map to define degree of spacing (0.0:2.0)
   real(cp), public :: alph2  ! Input parameter for non-linear map to define central point of different spacing (-1.0:1.0)
   character(len=72), public :: map_function ! Mapping family: either tangent or arcsin
   character(len=72), public :: time_scale  ! Time unit
   character(len=72), public :: time_scheme ! Time scheme
   character(len=72) :: bc_method    ! Galerkin or Tau-Lanczos method for BCs
   character(len=72) :: cheb_method  ! Chebyshev method: collocation, integration
   character(len=72) :: matrix_solve ! Either direct or influence matrix
   character(len=72) :: corio_term   ! Implicit or explicit treatment of Coriolis force
   character(len=72) :: buo_term    ! Implicit or explicit treatment of Buoyancy
   logical, public :: l_newmap      ! Switch for non-linear mapping (see Bayliss and Turkel, 1990)
   logical, public :: l_rerror_fix  ! Switch to fix round-off error in derivative calculation
   real(cp), public :: rerror_fac
   real(cp), public :: dtMin        ! Minimum allowed time step
   real(cp), public :: dtMax        ! Maximum allowed time step
   real(cp), public :: alpha        ! Weight for implicit time step
   real(cp), public :: courfac      ! Courant factor
   real(cp), public :: alffac       ! Courant factor for Alfven waves
   real(cp), public :: dt_fac       ! factor to control time step change
   logical, public :: l_cour_alf_damp ! Modify Alven Courant condition based on Christensen et al., GJI, 1999 (.true. by default)
   real(cp), public :: tEND
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
   real(cp), public :: amp_B
   integer,  public :: init_B
   real(cp), public :: scale_B
   logical,  public :: l_start_file     ! taking fields from startfile ?
   logical,  public :: l_reset_t ! Should we reset the time stored in the startfile?
   character(len=72), public :: start_file  ! name of start_file
   character(len=72), public :: start_file_b0  ! name of start_file for Background field

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
   logical,  public :: l_b_phiavg        ! Calculate the phi-avg of B3D at r_CMB
   logical,  public :: l_corr            ! Calculate the correlation
   real(cp), public :: bl_cut            ! Cut-off boundary layers in the force balance
   logical,  public :: l_3D            ! Require 3-D functions
   logical,  public :: l_heat_3D       ! 3-D treatment of temperature
   logical,  public :: l_thw_3D        ! Computation of the 3D thermal wind contribution on u_phi-3D
   logical,  public :: l_mag_3D        ! Treatment of the magnetic field (3-D)
   logical,  public :: l_mag_B0        ! With or without a Background field to the problem (Magnetoconvection)
   logical,  public :: l_U0_3D         ! With or without a Velocity field to the problem (Magnetoconvection)
   logical,  public :: l_mag_LF        ! Treatment of the Lorentz-Force
   logical,  public :: l_mag_alpha     ! With or without magnetic Alpha effect in magnetic advection
   logical,  public :: l_mag_pump      ! With or without Magnetic pumping
   logical,  public :: l_mag_inertia   ! With or without magnetic Alpha^2 effect triggered by inertia waves
   logical,  public :: l_leibniz       ! With or without liebniz rule of integration for the LF (finite-differences if not)
   logical,  public :: l_cyl           ! Accuracy of the schemes for the z-integrations
   logical,  public :: l_heat
   logical,  public :: l_chem
   logical,  public :: l_AB1
   logical,  public :: l_bridge_step
   logical,  public :: l_cheb_coll     ! Collocation method for Chebs
   logical,  public :: l_direct_solve  ! Direct solve or influence matrix method
   logical,  public :: l_coriolis_imp  ! Implicit treatment of Coriolis force
   logical,  public :: l_buo_imp       ! Implicit treatment of Buoyancy
   logical,  public :: l_QG_basis      ! Switch on the extra terms from the projection on a QG basis
   logical,  public :: l_path_rescale  ! Switch on the rescaling of the checkpoint fields as in [Aubert et al., 2017]
   logical,  public :: l_lin_solve     ! Switch off the non-Linear terms for studying the Onset of Convection
   real(cp), public :: tadvz_fac
   real(cp), public :: r_cmb           ! Outer core radius
   real(cp), public :: r_icb           ! Inner core radius
   real(cp), public :: CorFac, BuoFac, TdiffFac, ViscFac, XiDiffFac, ChemFac, &
   &                   BdiffFac, DyMagFac
   real(cp), public :: hdif_temp       ! Hyperdiffusion amplitude on temperature
   real(cp), public :: hdif_comp       ! Hyperdiffusion amplitude on composition
   real(cp), public :: hdif_vel        ! Hyperdiffusion amplitude on velocity
   real(cp), public :: damp_zon        ! Damping factor on non-axisymmetric azimuthal velocity
   real(cp), public :: hdif_mag        ! Hyperdiffusion amplitude on magnetic field
   integer,  public :: hdif_exp        ! Exponent of the hyperdiffusion profile
   integer,  public :: hdif_l          ! Latitudinal wavenumber for 3D hyper diffusion
   integer,  public :: hdif_m          ! Azimuthal wavenumber for hdif

   public :: read_namelists, write_namelists

contains

   subroutine read_namelists


      integer :: argument_count, input_handle, res
      character(len=100) :: input_filename
      logical :: nml_exist

      !-- Namelists:

      namelist/grid/n_r_max,n_cheb_max,m_max,minc
      namelist/grid_3D/n_r_max_3D,n_cheb_max_3D,n_phi_tot_3D,minc_3D
      namelist/control/tag,n_time_steps,alpha,l_newmap,map_function,&
      &                alph1,alph2,dtMax,courfac,alffac,tEND,       &
      &                runHours,runMinutes,runSeconds,l_non_rot,    &
      &                dt_fac,l_cour_alf_damp,n_fft_optim_lev,      &
      &                time_scheme,cheb_method,l_rerror_fix,        &
      &                rerror_fac,time_scale,matrix_solve,          &
      &                corio_term,buo_term,bc_method,l_QG_basis,    &
      &                hdif_temp,hdif_vel,damp_zon,hdif_comp,       &
      &                hdif_mag,hdif_exp,hdif_l,hdif_m,             &
      &                l_path_rescale,l_lin_solve
      !namelist/hdif/hdif_temp,hdif_vel,damp_zon,hdif_comp,hdif_mag,&
      !&             hdif_exp,hdif_l,hdif_m
      namelist/phys_param/ra,ek,pr,prmag,raxi,sc,radratio,g0,g1,g2,&
      &                   ktopt,kbott,ktopv,kbotv,l_ek_pump,       &
      &                   l_tcond_3D,tcond_fac,l_temp_advz,        &
      &                   beta_shift,epsrc0,epsrc0xi,ktopxi,kbotxi,&
      &                   t_bot,t_top,xi_bot,xi_top, ktopb,kbotb,  &
      &                   xicond_fac,l_xi_3D,l_heat_3D,l_thw_3D,   &
      &                   l_mag_LF,l_mag_alpha,l_mag_pump,l_cyl,   &
      &                   l_mag_inertia,alpha_fac,mag_pump_fac,    &
      &                   delta_fac, l_mag_B0, l_U0_3D, l_leibniz
      namelist/start_field/l_start_file,start_file,start_file_b0,scale_t,&
      &                    init_t,amp_t,scale_u,init_u,amp_u,l_reset_t,  &
      &                    amp_xi,init_xi,scale_xi,scale_B,init_B,amp_B
      namelist/output_control/n_log_step,n_checkpoints, n_checkpoint_step, &
      &                       n_frames, n_frame_step, n_specs, n_spec_step,&
      &                       l_vphi_balance,l_vort_balance,bl_cut,        &
      &                       l_2D_spectra, l_2D_SD, l_b_phiavg, l_corr

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
         if ( rank == 0 ) write(output_unit,*) '!  Reading grid parameters!'
         read(input_handle,nml=grid,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(output_unit,*) '!  No grid namelist found!'
         end if
         close(input_handle)

         open(newunit=input_handle,file=trim(input_filename))
         if ( rank == 0 ) write(output_unit,*) '!  Reading grid_3D parameters!'
         read(input_handle,nml=grid_3D,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(output_unit,*) '!  No grid_3D namelist found!'
         end if
         close(input_handle)

         open(newunit=input_handle,file=trim(input_filename))
         !-- Reading control parameters from namelists in STDIN:
         if ( rank == 0 ) write(output_unit,*) '!  Reading control parameters!'
         read(input_handle,nml=control,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(output_unit,*) '!  No control namelist found!'
         end if
         close(input_handle)

         !open(newunit=input_handle,file=trim(input_filename))
         !!-- Reading hdif parameters from namelists in STDIN:
         !if ( rank == 0 ) write(output_unit,*) '!  Reading hdif parameters!'
         !read(input_handle,nml=param_overdif,iostat=res)
         !if ( res /= 0 .and. rank == 0 ) then
         !   write(output_unit,*) '!  No hdif namelist found!'
         !end if
         !close(input_handle)

         open(newunit=input_handle,file=trim(input_filename))
         !-- Reading physical parameters from namelists in STDIN:
         if ( rank == 0 ) write(output_unit,*) '!  Reading physical parameters!'
         read(input_handle,nml=phys_param,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(output_unit,*) '!  No phys_param namelist found!'
         end if
         close(input_handle)

         open(newunit=input_handle,file=trim(input_filename))
         !-- Reading start field info from namelists in STDIN:
         if ( rank == 0 ) write(output_unit,*) '!  Reading start information!'
         read(input_handle,nml=start_field,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(output_unit,*) '! No start_field namelist found!'
         end if
         close(input_handle)

         open(newunit=input_handle,file=trim(input_filename))
         !-- Reading start field info from namelists in STDIN:
         if ( rank == 0 ) write(output_unit,*) &
         &                      '!  Reading output control information!'
         read(input_handle,nml=output_control,iostat=res)
         if ( res /= 0 .and. rank == 0 ) then
            write(output_unit,*) '! No output control namelist found!'
         end if
         close(input_handle)

      end if

      !-- Default is 2-D
      l_3D = .false.

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

      if ( ra /= 0.0_cp .and. l_heat_3D ) then
         l_heat = .false.
         l_3D = .true.
      else
         if ( ra == 0.0_cp ) then
            l_heat=.false.
         else
            l_heat=.true.
         end if
      end if

      if ( prmag == 0.0_cp ) then
         l_mag_3D=.false.
         l_mag_LF=.false.
         l_mag_B0=.false.
         l_mag_alpha=.false.
         l_mag_pump=.false.
         l_mag_inertia=.false.
      else
         l_mag_3D=.true.
      end if

      if ( .not. l_mag_alpha ) alpha_fac=0.0_cp
      if ( .not. l_mag_pump )  mag_pump_fac=0.0_cp
      if ( .not. l_mag_inertia )  delta_fac=0.0_cp

      if ( l_mag_3D ) l_3D = .true.

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

      !-- Determine solver method
      call capitalize(matrix_solve)
      if ( index(matrix_solve, 'DIRECT') /= 0 ) then
         l_direct_solve = .true.
      else
         l_direct_solve = .false.
      end if

      !-- Control the add.terms from the projection of the vorticity eq. onto a QG basis
      !-- Default is no additional terms l_QG_basis = .false.
      if ( l_QG_basis .or. l_lin_solve ) then
         !write(output_unit, &
         !&    '(" ! New-QG-projection approach only compatible with Collocations and Direct solver methods")')
         l_cheb_coll=.true. ! Only implemented for coll. method yet
         l_direct_solve=.true. ! Only implemented for direct_solve method yet
         l_galerkin = .false. ! Incompatible with Galerkin BC method yet
      end if

      !-- Warning:: Galerkin method can not handle inhomogeneous BCs (yet)
      if ( l_galerkin .and. index(cheb_method, 'COLL') == 0 .and. &
      &    ( maxval(abs(t_top(:))) > 0.0_cp .or.                  &
      &      maxval(abs(xi_top(:))) > 0.0_cp) ) then
         write(output_unit,*) '! Inhomogeneous Top BCs not compatible with galerkin yet!'
         call abortRun('! Inhomogeneous BCs not compatible with chosen chebyshev solver !')
      else if ( l_galerkin .and. index(cheb_method, 'COLL') == 0 .and. &
      &         ( maxval(abs(t_bot(:))) > 0.0_cp .or.                  &
      &           maxval(abs(xi_bot(:))) > 0.0_cp) ) then
         write(output_unit,*) '! Inhomogeneous Bot BCs not compatible with galerkin yet!'
         call abortRun('! Inhomogeneous BCs not compatible with chosen chebyshev solver !')
      end if

      if ( .not. ( damp_zon .eq. 1.0_cp ) ) then
         if ( rank == 0 ) write(output_unit,*) '! WARNING!:: damp_zon is probably physically inconsistent!'
         if ( damp_zon < 0.0_cp ) damp_zon = abs(damp_zon)
         if ( damp_zon > 1.0_cp ) damp_zon = one/damp_zon
         if ( damp_zon == 0.0_cp ) then
            write(output_unit,*) '! damp_zon = 0.0 --> Zonal wind has been totally suppressed!'
            call abortRun('! damp_zon = 0.0 --> Zonal wind has been totally suppressed !')
         end if
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

      if ( .not. l_cheb_coll ) then
         call capitalize(time_scheme)
         if ( time_scheme == 'CNLF'.or. time_scheme=='MODCNAB' .or.      &
         &    time_scheme == 'ARS443' .or. time_scheme=='ARS222' .or.    &
         &    time_scheme == 'BPR353' .or. time_scheme=='PC2' .or.       &
         &    time_scheme == 'LZ453' .or. time_scheme=='TVB33' .or.      &
         &    time_scheme == 'LZ232' .or. time_scheme == 'CK232' ) then
            l_buo_imp = .false.
            buo_term = 'EXP'
            if ( rank == 0 ) then
               write(output_unit, &
               &    '(" ! Implicit Buoyancy term not compatible with chosen time scheme")')
               write(output_unit, &
               &    '(" ! Buoyancy term will be treated explicitly")')
            end if
         end if
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

      !-- Time unit
      if ( l_mag_3D ) then
         BdiffFac = one/prmag !--- To be modified!
         DyMagFac = one/ek/prmag     !--- To be modified!
      end if
      call capitalize(time_scale) 
      if ( l_non_rot ) then
         CorFac = 0.0_cp
         if ( l_heat .or. l_heat_3D ) then
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
            if ( l_heat .or. l_heat_3D ) then
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
            if ( l_heat .or. l_heat_3D ) then
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

      !-- &grid_3D namelist
      n_r_max_3D      =17
      n_cheb_max_3D   =15
      n_phi_tot_3D    =96
      minc_3D         =1

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
      l_cour_alf_damp  =.true. ! By default, use Christensen's (GJI, 1999) CFL
      alffac           =1.0e3_cp
      dt_fac           =2.0_cp
      tEND             =0.0_cp    ! numerical time where run should end
      runHours         =0
      runMinutes       =10
      runSeconds       =0
      n_fft_optim_lev  =1 ! Optimisation level for FFT:
                          ! 0: FFTW_ESTIMATE, 1: FFTW_MEASURE, 2: FFTW_PATIENT
                          ! 3: FFTW_EXHAUSTIVE
      time_scheme      ='CNAB2'
      cheb_method      ='colloc'
      bc_method        ='tau-lanczos'
      matrix_solve     ='DIRECT'
      l_rerror_fix     =.true.
      rerror_fac       =500.0_cp
      corio_term       ='IMPLICIT' ! Implicit treatment of Coriolis term
      buo_term         ='IMPLICIT' ! Implicit treatment of Buoyancy
      time_scale       ='VISC' ! viscous units
      l_QG_basis       =.false. ! No QG basis projection of the vorticity eq.
      l_path_rescale   =.false. ! No rescaling of the checkpoint fields
      l_lin_solve      =.false. ! No switching off of the non-linear terms

      !-- Hyperdiffusion
      hdif_vel         =0.0_cp
      damp_zon         =1.0_cp
      hdif_temp        =0.0_cp
      hdif_comp        =0.0_cp
      hdif_mag         =0.0_cp
      hdif_exp         =0
      hdif_l           =0
      hdif_m           =0

      !-- Physcal parameters
      l_non_rot        =.false.
      l_ek_pump        =.false.
      l_temp_advz      =.false.
      l_tcond_3D       =.false.
      l_xi_3D          =.false.
      ra               =1.0e5_cp
      pr               =one
      prmag            =0.0_cp
      ek               =1.0e-3_cp
      raxi             =0.0e5_cp
      sc               =0.0_cp
      radratio         =0.35_cp
      tcond_fac        =one
      xicond_fac       =one
      beta_shift       =0.0_cp
      !----- Gravity parameters: defaut value g propto r (i.e. g1=1)
      g0               =0.0_cp
      g1               =one
      g2               =0.0_cp
      !----- Heat/Chemical sources in the outer core
      epsrc0           =0.0_cp
      epsrc0xi         =0.0_cp
      !----- Boundary conditions
      ktopt            =1
      kbott            =1
      ktopxi           =1
      kbotxi           =1
      ktopv            =2
      kbotv            =2
      ktopb            =1
      kbotb            =1
      !----- Parameters for temp. BCs
      do n=1,3*n_t_bounds
         t_bot(n) =0.0_cp
         t_top(n) =0.0_cp
         xi_bot(n)=0.0_cp
         xi_top(n)=0.0_cp
      end do

      l_3D = .false.
      l_cyl= .false.
      l_thw_3D= .true.
      !-- 3D treatment of temperature
      l_heat_3D = .false.

      !-- treatment of magnetic field (3D)
      l_mag_3D = .false.
      l_mag_LF = .false.
      l_mag_B0 = .false.
      l_U0_3D  = .false.
      l_leibniz = .false.
      l_mag_alpha = .false.
      alpha_fac   = 0.0_cp
      l_mag_pump = .false.
      mag_pump_fac = 0.0_cp
      l_mag_inertia = .false.
      delta_fac   = 0.0_cp

      !----- Namelist start_field:
      l_reset_t        =.false.
      l_start_file     =.false.
      start_file       ="no_start_file"
      start_file_b0    ="no_b0_start_file"
      init_t           =0
      scale_t          =1.0_cp
      amp_t            =0.0_cp
      init_xi          =0
      scale_xi         =1.0_cp
      amp_xi           =0.0_cp
      init_u           =0
      scale_u          =1.0_cp
      amp_u            =0.0_cp
      init_B           =0
      scale_B          =1.0_cp
      amp_B            =0.0_cp

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
      l_b_phiavg       =.false.
      l_corr           =.false.
      bl_cut           =1.0e-3_cp

   end subroutine default_namelists
!--------------------------------------------------------------------------------
   subroutine write_namelists(n_out)

      !-- Input file unit
      integer, intent(in) :: n_out

      !-- Local variables
      integer :: length, lm

      write(n_out,*) " "
      write(n_out,*) "&grid"
      write(n_out,'(''  n_r_max         ='',i5,'','')') n_r_max
      write(n_out,'(''  n_cheb_max      ='',i5,'','')') n_cheb_max
      write(n_out,'(''  m_max           ='',i5,'','')') m_max
      write(n_out,'(''  minc            ='',i5,'','')') minc
      write(n_out,*) "/"

      write(n_out,*) " "
      write(n_out,*) "&grid_3D"
      write(n_out,'(''  n_r_max_3D      ='',i5,'','')') n_r_max_3D
      write(n_out,'(''  n_cheb_max_3D   ='',i5,'','')') n_cheb_max_3D
      write(n_out,'(''  n_phi_tot_3D    ='',i5,'','')') n_phi_tot_3D
      write(n_out,'(''  minc_3D         ='',i5,'','')') minc_3D
      write(n_out,*) "/"

      write(n_out,*) "&control"
      length=length_to_blank(tag)
      write(n_out,*) " tag             = """,tag(1:length),""","
      write(n_out,'(''  n_time_steps    ='',i8,'','')') n_time_steps
      length=length_to_blank(time_scheme)
      write(n_out,*) " time_scheme     = """,time_scheme(1:length),""","
      length=length_to_blank(cheb_method)
      write(n_out,*) " cheb_method     = """,cheb_method(1:length),""","
      length=length_to_blank(bc_method)
      write(n_out,*) " bc_method       = """,bc_method(1:length),""","
      length=length_to_blank(matrix_solve)
      write(n_out,*) " matrix_solve    = """,matrix_solve(1:length),""","
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
      write(n_out,'(''  l_cour_alf_damp ='',l3,'','')') l_cour_alf_damp
      write(n_out,'(''  runHours        ='',i4,'','')') runHours
      write(n_out,'(''  runMinutes      ='',i4,'','')') runMinutes
      write(n_out,'(''  runSeconds      ='',i4,'','')') runSeconds
      write(n_out,'(''  tEND            ='',ES14.6,'','')') tEND
      write(n_out,'(''  l_non_rot       ='',l3,'','')') l_non_rot
      write(n_out,'(''  l_rerror_fix    ='',l3,'','')') l_rerror_fix
      write(n_out,'(''  rerror_fac      ='',ES14.6,'','')') rerror_fac
      write(n_out,'(''  n_fft_optim_lev ='',i4,'','')') n_fft_optim_lev
      write(n_out,'(''  l_QG_basis      ='',l3,'','')') l_QG_basis
      write(n_out,'(''  l_path_rescale  ='',l3,'','')') l_path_rescale
      write(n_out,'(''  l_lin_solve     ='',l3,'','')') l_lin_solve
      length=length_to_blank(time_scale)
      write(n_out,*) " time_scale      = """,time_scale(1:length),""","
      write(n_out,*) "/"

      write(n_out,*) "&hdif"
      length=length_to_blank(tag)
      write(n_out,'(''  hdif_vel        ='',ES14.6,'','')') hdif_vel
      write(n_out,'(''  damp_zon        ='',ES14.6,'','')') damp_zon
      write(n_out,'(''  hdif_temp       ='',ES14.6,'','')') hdif_temp
      write(n_out,'(''  hdif_comp       ='',ES14.6,'','')') hdif_comp
      write(n_out,'(''  hdif_mag        ='',ES14.6,'','')') hdif_mag
      write(n_out,'(''  hdif_exp        ='',i4,'','')') hdif_exp
      write(n_out,'(''  hdif_l          ='',i4,'','')') hdif_l
      write(n_out,'(''  hdif_m          ='',i4,'','')') hdif_m
      write(n_out,*) "/"

      write(n_out,*) "&phys_param"
      write(n_out,'(''  ra              ='',ES14.6,'','')') ra
      write(n_out,'(''  pr              ='',ES14.6,'','')') pr
      write(n_out,'(''  ek              ='',ES14.6,'','')') ek
      write(n_out,'(''  prmag           ='',ES14.6,'','')') prmag
      write(n_out,'(''  raxi            ='',ES14.6,'','')') raxi
      write(n_out,'(''  sc              ='',ES14.6,'','')') sc
      write(n_out,'(''  radratio        ='',ES14.6,'','')') radratio
      write(n_out,'(''  tcond_fac       ='',ES14.6,'','')') tcond_fac
      write(n_out,'(''  xicond_fac      ='',ES14.6,'','')') xicond_fac
      write(n_out,'(''  beta_shift      ='',ES14.6,'','')') beta_shift
      write(n_out,'(''  g0              ='',ES14.6,'','')') g0
      write(n_out,'(''  g1              ='',ES14.6,'','')') g1
      write(n_out,'(''  g2              ='',ES14.6,'','')') g2
      write(n_out,'(''  epsrc0          ='',ES14.6,'','')') epsrc0
      write(n_out,'(''  epsrc0xi        ='',ES14.6,'','')') epsrc0xi
      write(n_out,'(''  l_ek_pump       ='',l3,'','')') l_ek_pump
      write(n_out,'(''  l_temp_advz     ='',l3,'','')') l_temp_advz
      write(n_out,'(''  l_tcond_3D      ='',l3,'','')') l_tcond_3D
      write(n_out,'(''  l_3D            ='',l3,'','')') l_3D
      write(n_out,'(''  l_cyl           ='',l3,'','')') l_cyl
      write(n_out,'(''  l_heat_3D       ='',l3,'','')') l_heat_3D
      write(n_out,'(''  l_xi_3D         ='',l3,'','')') l_xi_3D
      if ( l_heat_3D ) &
      & write(n_out,'(''  l_thw_3D        ='',l3,'','')') l_thw_3D
      write(n_out,'(''  l_mag_3D        ='',l3,'','')') l_mag_3D
      write(n_out,'(''  l_mag_LF        ='',l3,'','')') l_mag_LF
      write(n_out,'(''  l_mag_B0        ='',l3,'','')') l_mag_B0
      write(n_out,'(''  l_U0_3D         ='',l3,'','')') l_U0_3D
      write(n_out,'(''  l_leibniz       ='',l3,'','')') l_leibniz
      write(n_out,'(''  l_mag_alpha     ='',l3,'','')') l_mag_alpha
      if ( l_mag_alpha ) &
      & write(n_out,'(''  alpha_fac       ='',ES14.6,'','')') alpha_fac
      write(n_out,'(''  l_mag_pump      ='',l3,'','')') l_mag_pump
      if ( l_mag_pump ) &
      & write(n_out,'(''  mag_pump_fac    ='',ES14.6,'','')') mag_pump_fac
      write(n_out,'(''  l_mag_inertia   ='',l3,'','')') l_mag_inertia
      if ( l_mag_inertia ) &
      & write(n_out,'(''  delta_fac       ='',ES14.6,'','')') delta_fac
      !--- Heat boundary condition:
      write(n_out,'(''  ktopt           ='',i3,'','')') ktopt
      write(n_out,'(''  kbott           ='',i3,'','')') kbott
      write(n_out,'(''  ktopxi          ='',i3,'','')') ktopxi
      write(n_out,'(''  kbotxi          ='',i3,'','')') kbotxi
      write(n_out,'(''  ktopv           ='',i3,'','')') ktopv
      write(n_out,'(''  kbotv           ='',i3,'','')') kbotv
      write(n_out,'(''  ktopb           ='',i3,'','')') ktopb
      write(n_out,'(''  kbotb           ='',i3,'','')') kbotb
      if ( l_heat_3D .and. ra /= 0.0_cp ) then
         write(n_out,'(''  Top Boundary Condition l,m,T:'')')
         do lm=1,3*n_t_bounds,4
            if ( t_top(lm+2) /= 0.0_cp .or. t_top(lm+3) /= 0.0_cp ) &
            &  write(n_out,'(''     '',2i3,2ES14.6,'','')') int(t_top(lm:lm+1)), t_top(lm+2:lm+3)
         end do
         write(n_out,'(''  Bot Boundary Condition l,m,T:'')')
         do lm=1,3*n_t_bounds,4
            if ( t_bot(lm+2) /= 0.0_cp .or. t_bot(lm+3) /= 0.0_cp ) &
            &  write(n_out,'(''     '',2i3,2ES14.6,'','')') int(t_bot(lm:lm+1)), t_bot(lm+2:lm+3)
         end do
      else if ( ra /= 0.0_cp ) then
         write(n_out,'(''  Top Boundary Condition m,T:'')')
         do lm=1,3*n_t_bounds,3
            if ( t_top(lm+1) /= 0.0_cp .or. t_top(lm+2) /= 0.0_cp ) &
            &  write(n_out,'(''     '',i3,2ES14.6,'','')') int(t_top(lm)), t_top(lm+1:lm+2)
         end do
         write(n_out,'(''  Bot Boundary Condition m,T:'')')
         do lm=1,3*n_t_bounds,3
            if ( t_bot(lm+1) /= 0.0_cp .or. t_bot(lm+2) /= 0.0_cp ) &
            &  write(n_out,'(''     '',i3,2ES14.6,'','')') int(t_bot(lm)), t_bot(lm+1:lm+2)
         end do
      end if
      if ( l_xi_3D .and. raxi /= 0.0_cp ) then
         write(n_out,'(''  Top Boundary Condition l,m,Xi:'')')
         do lm=1,3*n_t_bounds,4
            if ( xi_top(lm+2) /= 0.0_cp .or. xi_top(lm+3) /= 0.0_cp ) &
            &  write(n_out,'(''     '',2i3,2ES14.6,'','')') int(xi_top(lm:lm+1)), xi_top(lm+2:lm+3)
         end do
         write(n_out,'(''  Bot Boundary Condition l,m,Xi:'')')
         do lm=1,3*n_t_bounds,4
            if ( xi_bot(lm+2) /= 0.0_cp .or. xi_bot(lm+3) /= 0.0_cp ) &
            &  write(n_out,'(''     '',2i3,2ES14.6,'','')') int(xi_bot(lm:lm+1)), xi_bot(lm+2:lm+3)
         end do
      else if ( raxi /= 0.0_cp ) then
         write(n_out,'(''  Top Boundary Condition m,Xi:'')')
         do lm=1,3*n_t_bounds,3
            if ( xi_top(lm+1) /= 0.0_cp .or. xi_top(lm+2) /= 0.0_cp ) &
            &  write(n_out,'(''     '',i3,2ES14.6,'','')') int(xi_top(lm)), xi_top(lm+1:lm+2)
         end do
         write(n_out,'(''  Bot Boundary Condition m,Xi:'')')
         do lm=1,3*n_t_bounds,3
            if ( xi_bot(lm+1) /= 0.0_cp .or. xi_bot(lm+2) /= 0.0_cp ) &
            &  write(n_out,'(''     '',i3,2ES14.6,'','')') int(xi_bot(lm)), xi_bot(lm+1:lm+2)
         end do
      end if
      write(n_out,*) "/"

      write(n_out,*) "&start_field"
      write(n_out,'(''  l_start_file    ='',l3,'','')') l_start_file
      write(n_out,'(''  l_reset_t       ='',l3,'','')') l_reset_t
      length=length_to_blank(start_file)
      write(n_out,*) " start_file      = """,start_file(1:length),""","
     length=length_to_blank(start_file_b0)
      write(n_out,*) " start_file_b0   = """,start_file_b0(1:length),""","
      write(n_out,'(''  scale_t         ='',ES14.6,'','')') scale_t
      write(n_out,'(''  init_t          ='',i7,'','')') init_t
      write(n_out,'(''  amp_t           ='',ES14.6,'','')') amp_t
      write(n_out,'(''  scale_xi        ='',ES14.6,'','')') scale_xi
      write(n_out,'(''  init_xi         ='',i7,'','')') init_xi
      write(n_out,'(''  amp_xi          ='',ES14.6,'','')') amp_xi
      write(n_out,'(''  scale_u         ='',ES14.6,'','')') scale_u
      write(n_out,'(''  init_u          ='',i7,'','')') init_u
      write(n_out,'(''  amp_u           ='',ES14.6,'','')') amp_u
      write(n_out,'(''  scale_B         ='',ES14.6,'','')') scale_B
      write(n_out,'(''  init_B          ='',i7,'','')') init_B
      write(n_out,'(''  amp_B           ='',ES14.6,'','')') amp_B
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
      write(n_out,'(''  l_b_phiavg      ='',l3,'','')') l_b_phiavg
      write(n_out,'(''  l_corr          ='',l3,'','')') l_corr
      write(n_out,*) "/"

   end subroutine write_namelists
!--------------------------------------------------------------------------------
end module namelists
