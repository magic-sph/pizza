module dirk_schemes

   use precision_mod
   use parallel_mod
   use constants, only: one, half, two
   use mem_alloc, only: bytes_allocated
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array

   implicit none

   private

   type, public, extends(type_tscheme) :: type_dirk
      real(cp), allocatable :: butcher_imp(:,:)
      real(cp), allocatable :: butcher_exp(:,:)
      real(cp), allocatable :: butcher_ass_imp(:) ! Implicit Assembly stage
      real(cp), allocatable :: butcher_ass_exp(:) ! Explicit Assembly stage
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_weights
      procedure :: set_dt_array
      procedure :: set_imex_rhs
      procedure :: set_imex_rhs_ghost
      procedure :: rotate_imex
      procedure :: bridge_with_cnab2
      procedure :: start_with_ab1
      procedure :: assemble_imex
   end type type_dirk

contains

   subroutine initialize(this, time_scheme, courfac_nml)

      class(type_dirk) :: this

      !-- Input/output variables
      real(cp),          intent(in) :: courfac_nml
      character(len=72), intent(inout) :: time_scheme

      !-- Local variables
      integer :: sizet
      real(cp) :: courfac_loc

      allocate ( this%dt(1) )
      this%dt(:)=0.0_cp
      allocate ( this%wimp_lin(1) )
      this%wimp_lin(1)=0.0_cp

      this%family='DIRK'

      this%nold=1
      if ( index(time_scheme, 'ARS222') /= 0 ) then
         this%time_scheme = 'ARS222'
         this%nimp = 2
         this%nexp = 2
         this%nstages = 2
         this%istage = 1
         courfac_loc = 1.35_cp
         this%l_assembly = .false. ! No assembly stage
      else if ( index(time_scheme, 'ARS232') /= 0 ) then
         this%time_scheme = 'ARS232'
         this%nimp = 3
         this%nexp = 3
         this%nstages = 3
         this%istage = 1
         courfac_loc = 1.35_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'ARS233') /= 0 ) then
         this%time_scheme = 'ARS233'
         this%nimp = 3
         this%nexp = 3
         this%nstages = 3
         this%istage = 1
         courfac_loc = 1.35_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'LZ232') /= 0 ) then
         this%time_scheme = 'LZ232'
         this%nimp = 2
         this%nexp = 2
         this%nstages = 2
         this%istage = 1
         courfac_loc = 1.25_cp
         this%l_assembly = .false. ! No assembly stage
      else if ( index(time_scheme, 'CK232') /= 0 ) then
         this%time_scheme = 'CK232'
         this%nimp = 2
         this%nexp = 2
         this%nstages = 2
         this%istage = 1
         courfac_loc = 1.25_cp
         this%l_assembly = .false. ! No assembly stage
      else if ( index(time_scheme, 'ARS443') /= 0 ) then
         this%time_scheme = 'ARS443'
         this%nimp = 4
         this%nexp = 4
         this%nstages = 4
         this%istage = 1
         courfac_loc = 0.9_cp
         this%l_assembly = .false. ! No assembly stage
      else if ( index(time_scheme, 'LZ453') /= 0 ) then
         this%time_scheme = 'LZ453'
         this%nimp = 4
         this%nexp = 4
         this%nstages = 4
         this%istage = 1
         courfac_loc = 1.15_cp
         this%l_assembly = .false. ! No assembly stage
      else if ( index(time_scheme, 'BPR353') /= 0 ) then
         this%time_scheme = 'BPR353'
         this%nimp = 4
         this%nexp = 4
         this%nstages = 4
         this%istage = 1
         courfac_loc = 1.0_cp
         this%l_assembly = .false. ! No assembly stage
      else if ( index(time_scheme, 'PC2') /= 0 ) then
         this%time_scheme = 'PC2'
         this%nimp = 3
         this%nexp = 3
         this%nstages = 3
         this%istage = 1
         courfac_loc = 0.8_cp
         this%l_assembly = .false. ! No assembly stage
      else if ( index(time_scheme, 'ARS343' ) /= 0 ) then
         this%time_scheme = 'ARS343'
         this%nimp = 4
         this%nexp = 4
         this%nstages = 4
         this%istage = 1
         courfac_loc = 0.7_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'MARS343') /= 0 ) then
         this%time_scheme = 'MARS343'
         this%nimp = 4
         this%nexp = 4
         this%nstages = 4
         this%istage = 1
         courfac_loc = 0.5_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'KC343') /= 0 ) then
         this%time_scheme = 'KC343'
         this%nimp = 4
         this%nexp = 4
         this%nstages = 4
         this%istage = 1
         courfac_loc = 0.5_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'KC564') /= 0 ) then
         this%time_scheme = 'KC564'
         this%nimp = 6
         this%nexp = 6
         this%nstages = 6
         this%istage = 1
         courfac_loc = 0.5_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'KC785') /= 0 ) then
         this%time_scheme = 'KC785'
         this%nimp = 8
         this%nexp = 8
         this%nstages = 8
         this%istage = 1
         courfac_loc = 0.5_cp
         this%l_assembly = .true.
      else if ( index(time_scheme, 'DBM453') /= 0 ) then
         this%time_scheme = 'DBM453'
         this%nimp = 5
         this%nexp = 5
         this%nstages = 5
         this%istage = 1
         courfac_loc = 0.5_cp
         this%l_assembly = .true.
      end if

      if ( abs(courfac_nml) >= 1.0e3_cp ) then
         this%courfac=courfac_loc
      else
         this%courfac=courfac_nml
      end if

      if ( .not. this%l_assembly ) then
         sizet = this%nstages+1
      else
         sizet = this%nstages
      end if

      allocate( this%butcher_imp(sizet,sizet), &
      &         this%butcher_exp(sizet,sizet) )
      this%butcher_imp(:,:)=0.0_cp
      this%butcher_exp(:,:)=0.0_cp
      bytes_allocated=bytes_allocated+2*sizet*sizet*SIZEOF_DEF_REAL

      allocate( this%butcher_ass_imp(sizet), &
      &         this%butcher_ass_exp(sizet) )
      this%butcher_ass_imp(:)=0.0_cp
      this%butcher_ass_exp(:)=0.0_cp
      bytes_allocated=bytes_allocated+2*sizet*SIZEOF_DEF_REAL

      allocate( this%l_exp_calc(this%nstages) )
      allocate( this%l_imp_calc_rhs(this%nstages) )
      this%l_exp_calc(:) = .true.
      this%l_imp_calc_rhs(:) = .true.
      bytes_allocated=bytes_allocated+2*this%nstages*SIZEOF_LOGICAL

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)

      class(type_dirk) :: this

      deallocate( this%dt, this%wimp_lin, this%butcher_exp, this%butcher_imp )
      deallocate( this%butcher_ass_imp, this%butcher_ass_exp )
      deallocate( this%l_exp_calc, this%l_imp_calc_rhs )

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine set_weights(this, lMatNext)

      class(type_dirk) :: this
      logical, intent(inout) :: lMatNext

      !-- Local variables
      real(cp) :: wimp_old, del, gam, b1, b2

      wimp_old = this%wimp_lin(1)

      select case ( this%time_scheme )
         case ('ARS222')
            gam = half * (two-sqrt(two))
            del = one-half/gam
            this%wimp_lin(1) = gam
            this%butcher_imp(:,:) = reshape([  0.0_cp,  0.0_cp, 0.0_cp,  &
                                    &          0.0_cp,  gam   , 0.0_cp,  &
                                    &          0.0_cp, one-gam, gam    ],&
                                    &          [3,3],order=[2,1])
            this%butcher_exp(:,:) = reshape([  0.0_cp,  0.0_cp, 0.0_cp,  &
                                    &             gam,  0.0_cp, 0.0_cp,  &
                                    &             del, one-del, 0.0_cp], &
                                    &          [3,3],order=[2,1])
            this%l_imp_calc_rhs(1)=.false.
         case ('ARS233')
            gam = (3.0_cp+sqrt(3.0_cp))/6.0_cp
            this%wimp_lin(1) = gam
            this%butcher_imp(:,:) = reshape([  0.0_cp,  0.0_cp,    0.0_cp,  &
                                    &          0.0_cp,  gam   ,    0.0_cp,  &
                                    &          0.0_cp, one-two*gam,   gam], &
                                    &          [3,3], order=[2,1])
            this%butcher_ass_imp(:) = [0.0_cp, half, half]
            this%butcher_exp(:,:) = reshape([  0.0_cp,        0.0_cp, 0.0_cp,  &
                                    &             gam,        0.0_cp, 0.0_cp,  &
                                    &         gam-one, two*(one-gam), 0.0_cp], &
                                    &          [3,3], order=[2,1])
            this%butcher_ass_exp(:) = [0.0_cp, half, half]
            this%l_imp_calc_rhs(1)=.false.
         case ('ARS232')
            gam = half * (two-sqrt(two))
            del = -two*sqrt(two)/3.0_cp
            this%wimp_lin(1) = gam
            this%butcher_imp(:,:) = reshape([  0.0_cp,  0.0_cp, 0.0_cp,  &
                                    &          0.0_cp,  gam   , 0.0_cp,  &
                                    &          0.0_cp, one-gam,    gam], &
                                    &          [3,3], order=[2,1])
            this%butcher_ass_imp(:) = [0.0_cp, one-gam, gam]
            this%butcher_exp(:,:) = reshape([  0.0_cp,  0.0_cp, 0.0_cp,  &
                                    &             gam,  0.0_cp, 0.0_cp,  &
                                    &             del, one-del, 0.0_cp], &
                                    &          [3,3], order=[2,1])
            this%butcher_ass_exp(:) = [0.0_cp, one-gam, gam]
            this%l_imp_calc_rhs(1)=.false.
         case ('LZ232')
            this%wimp_lin(1) = half
            this%butcher_imp(:,:) = reshape([  0.0_cp,  0.0_cp, 0.0_cp,  &
                                    &        -0.25_cp,    half, 0.0_cp,  &
                                    &            half,  0.0_cp, half ],  &
                                    &          [3,3],order=[2,1])
            this%butcher_exp(:,:) = reshape([  0.0_cp, 0.0_cp, 0.0_cp,  &
                                    &         0.25_cp, 0.0_cp, 0.0_cp,  &
                                    &            -one,    two, 0.0_cp], &
                                    &          [3,3],order=[2,1])
         case ('CK232')
            gam = one-half*sqrt(two)
            this%wimp_lin(1) = gam
            this%butcher_imp(:,:) = reshape([                                     &
            &                         0.0_cp,                     0.0_cp, 0.0_cp, &
            &  -1.0_cp/3.0_cp+half*sqrt(two),                        gam, 0.0_cp, &
            &      0.75_cp-0.25_cp*sqrt(two), -0.75_cp+0.75_cp*sqrt(two),    gam],&
            &                               [3,3],order=[2,1])
            this%butcher_exp(:,:) = reshape([        0.0_cp,  0.0_cp, 0.0_cp,  &
                                    &         2.0_cp/3.0_cp,  0.0_cp, 0.0_cp,  &
                                    &               0.25_cp, 0.75_cp, 0.0_cp], &
                                    &        [3,3],order=[2,1])
         case ('ARS443')
            this%wimp_lin(1) = half
            this%butcher_imp(:,:) = reshape(                               &
            &            [ 0.0_cp,       0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,   &
            &              0.0_cp,         half, 0.0_cp, 0.0_cp, 0.0_cp,   &
            &              0.0_cp,1.0_cp/6.0_cp,   half, 0.0_cp, 0.0_cp,   &
            &              0.0_cp,        -half,   half,   half, 0.0_cp,   &
            &              0.0_cp,       1.5_cp,-1.5_cp,   half,   half],  &
            &              [5,5],order=[2,1])
            this%butcher_exp(:,:) = reshape(                                     &
            &       [           0.0_cp,       0.0_cp,  0.0_cp,   0.0_cp, 0.0_cp, &
            &                    half,        0.0_cp,  0.0_cp,   0.0_cp, 0.0_cp, &
            &         11.0_cp/18.0_cp,1.0_cp/18.0_cp,  0.0_cp,   0.0_cp, 0.0_cp, &
            &           5.0_cp/6.0_cp,-5.0_cp/6.0_cp,    half,   0.0_cp, 0.0_cp, &
            &                 0.25_cp,       1.75_cp, 0.75_cp, -1.75_cp, 0.0_cp],&
            &         [5,5],order=[2,1])
            this%l_imp_calc_rhs(1)=.false.
         case ('BPR353')
            this%wimp_lin(1) = half
            this%butcher_imp(:,:) = reshape(                                   &
            &       [         0.0_cp,        0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,   &
            &                   half,          half, 0.0_cp, 0.0_cp, 0.0_cp,   &
            &         5.0_cp/18.0_cp,-1.0_cp/9.0_cp,   half, 0.0_cp, 0.0_cp,   &
            &                   half,        0.0_cp, 0.0_cp,   half, 0.0_cp,   &
            &                0.25_cp,        0.0_cp,0.75_cp,  -half,   half],  &
            &              [5,5],order=[2,1])
            this%butcher_exp(:,:) = reshape(                                   &
            &       [        0.0_cp,       0.0_cp,  0.0_cp,   0.0_cp, 0.0_cp,  &
            &                   one,       0.0_cp,  0.0_cp,   0.0_cp, 0.0_cp,  &
            &         4.0_cp/9.0_cp,2.0_cp/9.0_cp,  0.0_cp,   0.0_cp, 0.0_cp,  &
            &               0.25_cp,       0.0_cp, 0.75_cp,   0.0_cp, 0.0_cp,  &
            &               0.25_cp,       0.0_cp, 0.75_cp,   0.0_cp, 0.0_cp], &
            &         [5,5],order=[2,1])
            this%l_exp_calc(4)=.false. ! No need to calculte the explicit solve
         case ('LZ453')
            this%wimp_lin(1) = 1.2_cp
            this%butcher_imp(:,:) = reshape(                                                        &
            &[            0.0_cp,           0.0_cp,            0.0_cp,            0.0_cp, 0.0_cp,   &
            &   -44.0_cp/45.0_cp,           1.2_cp,            0.0_cp,            0.0_cp, 0.0_cp,   &
            &  -47.0_cp/300.0_cp,         -0.71_cp,            1.2_cp,            0.0_cp, 0.0_cp,   &
            &           3.375_cp,         -3.25_cp, -59.0_cp/120.0_cp,            1.2_cp, 0.0_cp,   &
            &    89.0_cp/50.0_cp,-486.0_cp/55.0_cp,            8.9_cp,-562.0_cp/275.0_cp, 1.2_cp],  &
            &[5,5],order=[2,1])
            this%butcher_exp(:,:) = reshape(                                   &
            & [            0.0_cp,            0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,  &
            &       2.0_cp/9.0_cp,            0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,  &
            &    71.0_cp/420.0_cp,  23.0_cp/140.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,  &
            &  -281.0_cp/336.0_cp, 187.0_cp/112.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,  &
            &             0.1_cp,             0.0_cp, 0.5_cp, 0.4_cp, 0.0_cp], &
            &  [5,5],order=[2,1])
         case ('PC2')
            this%wimp_lin(1) = half
            this%butcher_imp(:,:) = reshape([ 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, &
            &                                   half,   half, 0.0_cp, 0.0_cp, &
            &                                   half, 0.0_cp,   half, 0.0_cp, &
            &                                   half, 0.0_cp, 0.0_cp,   half],&
            &                               [4,4],order=[2,1])
            this%butcher_exp(:,:) = reshape([ 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, &
            &                                    one, 0.0_cp, 0.0_cp, 0.0_cp, &
            &                                   half,   half, 0.0_cp, 0.0_cp, &
            &                                   half, 0.0_cp,   half, 0.0_cp],&
            &                               [4,4],order=[2,1])
         case ('ARS343')
            gam = 0.435866521508459_cp
            b1 = -1.5_cp*gam*gam+4.0_cp*gam-0.25_cp
            b2 = 1.5_cp*gam*gam-5.0_cp*gam+1.25_cp
            this%wimp_lin(1) = gam
            this%butcher_imp(:,:) = reshape([                           &
            &                  0.0_cp,         0.0_cp, 0.0_cp,  0.0_cp, &
            &                  0.0_cp,            gam, 0.0_cp,  0.0_cp, &
            &                  0.0_cp, half*(one-gam),    gam,  0.0_cp, &
            &                  0.0_cp,             b1,     b2,    gam], &
            &                               [4,4], order=[2,1])
            this%butcher_ass_imp(:) = [0.0_cp, b1, b2, gam]
            this%butcher_ass_exp(:) = this%butcher_ass_imp(:)

            this%butcher_exp(:,:) = reshape([                                           &
            &                0.0_cp,               0.0_cp,               0.0_cp,0.0_cp, &
            &                   gam,               0.0_cp,               0.0_cp,0.0_cp, &
            & 0.3212788860286278_cp,0.3966543747256017_cp,               0.0_cp,0.0_cp, &
            &-0.1058582960718797_cp,0.5529291480359398_cp,0.5529291480359398_cp,0.0_cp],&
            &                               [4,4], order=[2,1])
            this%l_imp_calc_rhs(1)=.false.
         case ('KC343')
            gam = 1767732205903.0_cp/4055673282236.0_cp
            this%wimp_lin(1) = gam

            this%butcher_imp(:,:) = reshape([0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,        &
            &    gam, gam, 0.0_cp, 0.0_cp, 2746238789719.0_cp/10658868560708.0_cp,  &
            &    -640167445237.0_cp/6845629431997.0_cp,  gam, 0.0_cp,               &
            &    1471266399579.0_cp/7840856788654.0_cp, -4482444167858.0_cp/        &
            &    7529755066697.0_cp, 11266239266428.0_cp/11593286722821.0_cp, gam], &
            &                     [4,4],  order=[2,1])

            this%butcher_ass_imp(:) = [1471266399579.0_cp/7840856788654.0_cp,  &
            &                          -4482444167858.0_cp/7529755066697.0_cp, &
            &                          11266239266428.0_cp/11593286722821.0_cp,&
            &                          gam]

            this%butcher_exp(:,:) = reshape([0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,    &
            &    1767732205903.0_cp/2027836641118.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, &
            &    5535828885825.0_cp/10492691773637.0_cp, 788022342437.0_cp/     &
            &    10882634858940.0_cp,0.0_cp, 0.0_cp, 6485989280629.0_cp/        &
            &    16251701735622.0_cp,-4246266847089.0_cp/9704473918619.0_cp,    &
            &    10755448449292.0_cp/10357097424841.0_cp, 0.0_cp], [4,4],       &
            &    order=[2,1])
            this%butcher_ass_exp(:)=this%butcher_ass_imp(:)
         case ('MARS343')
            gam = 0.435866521508458_cp
            b1 =  1.20849664917601_cp
            b2 = -0.64436317068446_cp
            this%wimp_lin(1) = gam
            this%butcher_imp(:,:) = reshape([                           &
                            &  0.0_cp,         0.0_cp, 0.0_cp,  0.0_cp, &
                            &  0.0_cp,            gam, 0.0_cp,  0.0_cp, &
                            &  0.0_cp, half*(one-gam),    gam,  0.0_cp, &
                            &  0.0_cp,             b1,     b2,    gam], &
                            &          [4,4],order=[2,1])
            this%butcher_ass_imp(:) = [0.0_cp, b1, b2, gam]
            this%butcher_ass_exp(:) = this%butcher_ass_imp(:)
            this%butcher_exp(:,:) = reshape([                                    &
                   &          0.0_cp,          0.0_cp,         0.0_cp,  0.0_cp,  &
                   & 0.877173304301691_cp,     0.0_cp,         0.0_cp,  0.0_cp,  &
                   & 0.535396540307354_cp, 0.182536720446875_cp, 0.0_cp,  0.0_cp,&
                   & 0.63041255815287_cp, -0.83193390106308_cp, 1.20152134291021_cp,&
                   & 0.0_cp], [4,4], order=[2,1])
            this%l_imp_calc_rhs(1)=.false.
         case ('KC564')
            this%wimp_lin(1) = 0.25_cp

            this%butcher_imp(:,:) = reshape([0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,   &
            &    0.0_cp, 0.0_cp, 0.25_cp, 0.25_cp, 0.0_cp, 0.0_cp, 0.0_cp,     &
            &    0.0_cp, 8611.0_cp/62500.0_cp, -1743.0_cp/31250.0_cp, 0.25_cp, &
            &    0.0_cp, 0.0_cp, 0.0_cp, 5012029.0_cp/34652500.0_cp,           &
            &    -654441.0_cp/2922500.0_cp, 174375.0_cp/388108.0_cp, 0.25_cp,  &
            &    0.0_cp, 0.0_cp, 15267082809.0_cp/155376265600.0_cp,           &
            &    -71443401.0_cp/120774400.0_cp, 730878875.0_cp/902184768.0_cp, &
            &    2285395.0_cp/8070912.0_cp, 0.25_cp, 0.0_cp, 82889.0_cp/       &
            &    524892.0_cp, 0.0_cp, 15625.0_cp/83664.0_cp, 69875.0_cp/       &
            &    102672.0_cp, -2260.0_cp/8211.0_cp, 0.25_cp], [6,6], order=[2,1])

            this%butcher_ass_imp(:) = [82889.0_cp/524892.0_cp, 0.0_cp, &
            &                          15625.0_cp/83664.0_cp,          &
            &                          69875.0_cp/102672.0_cp,         &
            &                          -2260.0_cp/8211.0_cp, 0.25_cp]

            this%butcher_exp(:,:) = reshape([ 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,    &
            &    0.0_cp, 0.0_cp, 0.5_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, &
            &    13861.0_cp/62500.0_cp,6889.0_cp/62500.0_cp, 0.0_cp, 0.0_cp,     &
            &    0.0_cp, 0.0_cp, -116923316275.0_cp/2393684061468.0_cp,          &
            &    -2731218467317.0_cp/15368042101831.0_cp, 9408046702089.0_cp/    &
            &    11113171139209.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, -451086348788.0_cp/&
            &    2902428689909.0_cp, -2682348792572.0_cp/7519795681897.0_cp,     &
            &    12662868775082.0_cp/11960479115383.0_cp, 3355817975965.0_cp/    &
            &    11060851509271.0_cp, 0.0_cp, 0.0_cp, 647845179188.0_cp/         &
            &    3216320057751.0_cp, 73281519250.0_cp/8382639484533.0_cp,        &
            &    552539513391.0_cp/3454668386233.0_cp, 3354512671639.0_cp/       &
            &    8306763924573.0_cp, 4040.0_cp/17871.0_cp, 0.0_cp], [6,6],       &
            &    order=[2,1])

            this%butcher_ass_exp(:)=this%butcher_ass_imp(:)
         case ('DBM453')
            this%wimp_lin(1) =  0.32591194130117247_cp
            this%butcher_imp(:,:) = reshape(                                     &
            &       [         0.0_cp,        0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,     &
            &  -0.22284985318525410_cp, 0.32591194130117247_cp, 0.0_cp, 0.0_cp,  &
            &          0.0_cp, -0.46801347074080545_cp, 0.86349284225716961_cp,  &
            &   0.32591194130117247_cp, 0.0_cp, 0.0_cp, -0.46509906651927421_cp, &
            &   0.81063103116959553_cp, 0.61036726756832357_cp,                  &
            &   0.32591194130117247_cp, 0.0_cp, 0.87795339639076675_cp,          &
            &  -0.72692641526151547_cp, 0.75204137157372720_cp,                  &
            &  -0.22898029400415088_cp, 0.32591194130117247_cp], [5,5],order=[2,1])
            this%butcher_exp(:,:) = reshape(                                     &
            &  [ 0.0_cp,       0.0_cp,  0.0_cp,   0.0_cp, 0.0_cp,                &
            &    0.10306208811591838_cp,  0.0_cp,  0.0_cp,   0.0_cp, 0.0_cp,     &
            &   -0.94124866143519894_cp, 1.6626399742527356_cp, 0.0_cp, 0.0_cp,  &
            &    0.0_cp, -1.3670975201437765_cp,  1.3815852911016873_cp,         &
            &    1.2673234025619065_cp, 0.0_cp, 0.0_cp, -0.81287582068772448_cp, &
            &    0.81223739060505738_cp, 0.90644429603699305_cp,                 &
            &    0.094194134045674111_cp, 0.0_cp], [5,5],order=[2,1])

            this%butcher_ass_exp(:)=[0.87795339639076672_cp, -0.72692641526151549_cp, &
                                     0.7520413715737272_cp, -0.22898029400415090_cp,  &
                                     0.32591194130117247_cp]
            this%butcher_ass_imp(:)=this%butcher_ass_exp(:)
         case ('KC785')
            this%wimp_lin(1) = 2.0_cp/9.0_cp

            this%butcher_imp(:,:) = reshape([0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,          &
            &  0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 2.0_cp/9.0_cp, 2.0_cp/9.0_cp,          &
            &  0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 2366667076620.0_cp/    &
            &  8822750406821.0_cp, 2366667076620.0_cp/8822750406821.0_cp,             &
            &  2.0_cp/9.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,                 &
            &  -257962897183.0_cp/4451812247028.0_cp, -257962897183.0_cp/             &
            &  4451812247028.0_cp, 128530224461.0_cp/14379561246022.0_cp,             &
            &  2.0_cp/9.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, -486229321650.0_cp/     &
            &  11227943450093.0_cp, -486229321650.0_cp/11227943450093.0_cp,           &
            &  -225633144460.0_cp/6633558740617.0_cp, 1741320951451.0_cp/             &
            &  6824444397158.0_cp, 2.0_cp/9.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,             &
            &  621307788657.0_cp/4714163060173.0_cp, 621307788657.0_cp/               &
            &  4714163060173.0_cp, -125196015625.0_cp/3866852212004.0_cp,             &
            &  940440206406.0_cp/7593089888465.0_cp, 961109811699.0_cp/               &
            &  6734810228204.0_cp, 2.0_cp/9.0_cp, 0.0_cp, 0.0_cp, 2036305566805.0_cp/ &
            &  6583108094622.0_cp, 2036305566805.0_cp/6583108094622.0_cp,             &
            &  -3039402635899.0_cp/4450598839912.0_cp, -1829510709469.0_cp/           &
            &  31102090912115.0_cp, -286320471013.0_cp/6931253422520.0_cp,            &
            &  8651533662697.0_cp/9642993110008.0_cp, 2.0_cp/9.0_cp, 0.0_cp, 0.0_cp,  &
            &  0.0_cp, 3517720773327.0_cp/20256071687669.0_cp, 4569610470461.0_cp/    &
            &  17934693873752.0_cp,  2819471173109.0_cp/11655438449929.0_cp,          &
            &  3296210113763.0_cp/10722700128969.0_cp, -1142099968913.0_cp/           &
            &  5710983926999.0_cp, 2.0_cp/9.0_cp], [8,8], order=[2,1])

            this%butcher_ass_imp(:) = [0.0_cp, 0.0_cp, 3517720773327.0_cp/    &
            &    20256071687669.0_cp, 4569610470461.0_cp/17934693873752.0_cp, &
            &    2819471173109.0_cp/11655438449929.0_cp, 3296210113763.0_cp/  &
            &    10722700128969.0_cp, -1142099968913.0_cp/5710983926999.0_cp, &
            &    2.0_cp/9.0_cp]

            this%butcher_exp(:,:) = reshape([ 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,  &
            &  0.0_cp, 0.0_cp, 0.0_cp, 4.0_cp/9.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,  &
            &  0.0_cp, 0.0_cp, 0.0_cp, 1.0_cp/9.0_cp, 1183333538310.0_cp/              &
            &  1827251437969.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,     &
            &  895379019517.0_cp/9750411845327.0_cp, 477606656805.0_cp/                &
            &  13473228687314.0_cp, -112564739183.0_cp/9373365219272.0_cp, 0.0_cp,     &
            &  0.0_cp, 0.0_cp, 0.0_cp,0.0_cp, -4458043123994.0_cp/13015289567637.0_cp, &
            &  -2500665203865.0_cp/9342069639922.0_cp, 983347055801.0_cp/              &
            &  8893519644487.0_cp, 2185051477207.0_cp/2551468980502.0_cp, 0.0_cp,      &
            &  0.0_cp, 0.0_cp, 0.0_cp, -167316361917.0_cp/17121522574472.0_cp,         &
            &  1605541814917.0_cp/7619724128744.0_cp, 991021770328.0_cp/               &
            &  13052792161721.0_cp,  2342280609577.0_cp/11279663441611.0_cp,           &
            &  3012424348531.0_cp/12792462456678.0_cp, 0.0_cp, 0.0_cp, 0.0_cp,         &
            &  6680998715867.0_cp/14310383562358.0_cp, 5029118570809.0_cp/             &
            &  3897454228471.0_cp, 2415062538259.0_cp/6382199904604.0_cp,              &
            &  -3924368632305.0_cp/6964820224454.0_cp, -4331110370267.0_cp/            &
            &  15021686902756.0_cp, -3944303808049.0_cp/11994238218192.0_cp, 0.0_cp,   &
            &  0.0_cp, 2193717860234.0_cp/3570523412979.0_cp, 2193717860234.0_cp/      &
            &  3570523412979.0_cp, 5952760925747.0_cp/18750164281544.0_cp,             &
            &  -4412967128996.0_cp/6196664114337.0_cp, 4151782504231.0_cp/             &
            &  36106512998704.0_cp,  572599549169.0_cp/6265429158920.0_cp,             &
            &  -457874356192.0_cp/11306498036315.0_cp, 0.0_cp], [8,8], order=[2,1])

            this%butcher_ass_exp(:)=this%butcher_ass_imp(:)
      end select

      this%wimp_lin(1)      = this%dt(1)*this%wimp_lin(1)
      this%butcher_imp(:,:) = this%dt(1)*this%butcher_imp(:,:)
      this%butcher_exp(:,:) = this%dt(1)*this%butcher_exp(:,:)

      this%butcher_ass_imp(:) = this%dt(1)*this%butcher_ass_imp(:)
      this%butcher_ass_exp(:) = this%dt(1)*this%butcher_ass_exp(:)
         
   end subroutine set_weights
!------------------------------------------------------------------------------
   subroutine set_dt_array(this, dt_new, dt_min, time, n_log_file,  &
              &            n_time_step, l_new_dtNext)
      !
      ! This subroutine adjusts the time step
      !

      class(type_dirk) :: this

      !-- Input variables
      real(cp), intent(in) :: dt_new
      real(cp), intent(in) :: dt_min
      real(cp), intent(in) :: time
      integer,  intent(in) :: n_log_file
      integer,  intent(in) :: n_time_step
      logical,  intent(in) :: l_new_dtNext

      !-- Local variables
      real(cp) :: dt_old

      dt_old = this%dt(1)

      !-- Then overwrite the first element by the new timestep
      this%dt(1)=dt_new

      !----- Stop if time step has become too small:
      if ( dt_new < dt_min ) then
         if ( rank == 0 ) then
            write(*,'(1p,/,A,ES14.4,/,A)')             &
            &    " ! Time step too small, dt=",dt_new, &
            &    " ! I thus stop the run !"
            write(n_log_file,'(1p,/,A,ES14.4,/,A)')    &
            &    " ! Time step too small, dt=",dt_new, &
            &    " ! I thus stop the run !"
         end if
         call abortRun('Stop run in steptime!')
      end if

      if ( l_new_dtNext ) then
         !------ Writing info and getting new weights:
         if ( rank == 0 ) then
            write(*,'(1p,/,A,ES18.10,/,A,i9,/,A,ES15.8,/,A,ES15.8)')  &
            &    " ! Changing time step at time=",(time+this%dt(1)),  &
            &    "                 time step no=",n_time_step,        &
            &    "                      last dt=",dt_old,             &
            &    "                       new dt=",dt_new
            write(n_log_file,                                         &
            &    '(1p,/,A,ES18.10,/,A,i9,/,A,ES15.8,/,A,ES15.8)')     &
            &    " ! Changing time step at time=",(time+this%dt(1)),  &
            &    "                 time step no=",n_time_step,        &
            &    "                      last dt=",dt_old,             &
            &    "                       new dt=",dt_new
         end if
      end if

   end subroutine set_dt_array
!------------------------------------------------------------------------------
   subroutine set_imex_rhs(this, rhs, dfdt)
      !
      ! This subroutine assembles the right-hand-side of an IMEX scheme
      !

      class(type_dirk) :: this

      !-- Input variables:
      type(type_tarray), intent(in) :: dfdt

      !-- Output variable
      complex(cp), intent(out) :: rhs(dfdt%lm:dfdt%um,dfdt%nRstart:dfdt%nRstop)

      !-- Local variables
      integer :: n_stage, n_r, start_m, stop_m

      !$omp parallel default(shared) private(start_m, stop_m)
      start_m=dfdt%lm; stop_m=dfdt%um
      call get_openmp_blocks(start_m,stop_m)

      do n_r=dfdt%nRstart,dfdt%nRstop
         rhs(start_m:stop_m,n_r)=dfdt%old(start_m:stop_m,n_r,1)
      end do

      do n_stage=1,this%istage
         do n_r=dfdt%nRstart,dfdt%nRstop
            rhs(start_m:stop_m,n_r)=rhs(start_m:stop_m,n_r) +                &
            &                       this%butcher_exp(this%istage+1,n_stage)* &
            &                       dfdt%expl(start_m:stop_m,n_r,n_stage)
         end do
      end do

      do n_stage=1,this%istage
         do n_r=dfdt%nRstart,dfdt%nRstop
            rhs(start_m:stop_m,n_r)=rhs(start_m:stop_m,n_r) +                &
            &                       this%butcher_imp(this%istage+1,n_stage)* &
            &                       dfdt%impl(start_m:stop_m,n_r,n_stage)
         end do
      end do
      !$omp end parallel

   end subroutine set_imex_rhs
!------------------------------------------------------------------------------
   subroutine set_imex_rhs_ghost(this, rhs, dfdt, start_m, stop_m, ng)
      !
      ! This subroutine assembles the right-hand-side of an IMEX scheme in case
      ! an array with ghosts zones is provided
      !

      class(type_dirk) :: this

      !-- Input variables:
      type(type_tarray), intent(in) :: dfdt
      integer,           intent(in) :: start_m ! Starting m index
      integer,           intent(in) :: stop_m  ! Stopping m index
      integer,           intent(in) :: ng       ! Number of ghost zones

      !-- Output variable
      complex(cp), intent(out) :: rhs(dfdt%lm:dfdt%um,dfdt%nRstart-ng:dfdt%nRstop+ng)

      !-- Local variables
      integer :: n_stage, n_r

      do n_r=dfdt%nRstart,dfdt%nRstop
         rhs(start_m:stop_m,n_r)=dfdt%old(start_m:stop_m,n_r,1)
      end do

      do n_stage=1,this%istage
         do n_r=dfdt%nRstart,dfdt%nRstop
            rhs(start_m:stop_m,n_r)=rhs(start_m:stop_m,n_r) +                &
            &                       this%butcher_exp(this%istage+1,n_stage)* &
            &                       dfdt%expl(start_m:stop_m,n_r,n_stage)
         end do
      end do

      do n_stage=1,this%istage
         do n_r=dfdt%nRstart,dfdt%nRstop
            rhs(start_m:stop_m,n_r)=rhs(start_m:stop_m,n_r) +                &
            &                       this%butcher_imp(this%istage+1,n_stage)* &
            &                       dfdt%impl(start_m:stop_m,n_r,n_stage)
         end do
      end do

   end subroutine set_imex_rhs_ghost
!------------------------------------------------------------------------------
   subroutine rotate_imex(this, dfdt)
      !
      ! This subroutine is used to roll the time arrays from one time step
      !

      class(type_dirk) :: this

      !-- Output variables:
      type(type_tarray), intent(inout) :: dfdt

   end subroutine rotate_imex
!------------------------------------------------------------------------------
   subroutine bridge_with_cnab2(this)

      class(type_dirk) :: this

   end subroutine bridge_with_cnab2
!------------------------------------------------------------------------------
   subroutine start_with_ab1(this)

      class(type_dirk) :: this

   end subroutine start_with_ab1
!------------------------------------------------------------------------------
   subroutine assemble_imex(this, rhs, dfdt)

      class(type_dirk) :: this

      !-- Input variables
      type(type_tarray), intent(in) :: dfdt
      
      !-- Ouput variable
      complex(cp), intent(out) :: rhs(dfdt%lm:dfdt%um,dfdt%nRstart:dfdt%nRstop)

      !-- Local variables
      integer :: n_stage, n_r, start_m, stop_m

      !$omp parallel default(shared) private(start_m,stop_m,n_r)
      start_m=dfdt%lm; stop_m=dfdt%um
      call get_openmp_blocks(start_m,stop_m)

      do n_r=dfdt%nRstart,dfdt%nRstop
         rhs(start_m:stop_m,n_r)=dfdt%old(start_m:stop_m,n_r,1)
      end do

      do n_stage=1,this%nstages
         do n_r=dfdt%nRstart,dfdt%nRstop
            rhs(start_m:stop_m,n_r)=rhs(start_m:stop_m,n_r) +        &
            &                       this%butcher_ass_exp(n_stage)*       &
            &                       dfdt%expl(start_m:stop_m,n_r,n_stage)
         end do
      end do

      do n_stage=1,this%nstages
         do n_r=dfdt%nRstart,dfdt%nRstop
            rhs(start_m:stop_m,n_r)=rhs(start_m:stop_m,n_r) +       &
            &                       this%butcher_ass_imp(n_stage)*      &
            &                       dfdt%impl(start_m:stop_m,n_r,n_stage)
         end do
      end do
      !$omp end parallel


   end subroutine assemble_imex
!------------------------------------------------------------------------------
end module dirk_schemes
