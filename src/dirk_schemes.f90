module dirk_schemes

   use precision_mod
   use parallel_mod
   use namelists, only: alpha, l_cheb_coll
   use constants, only: one, half, two, ci, zero
   use mem_alloc, only: bytes_allocated
   use useful, only: abortRun
   use truncation, only: idx2m
   use radial_functions, only: or1
   use time_schemes, only: type_tscheme
   use time_array

   implicit none

   private

   type, public, extends(type_tscheme) :: type_dirk
      real(cp), allocatable :: butcher_imp(:,:)
      real(cp), allocatable :: butcher_exp(:,:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_weights
      procedure :: set_dt_array
      procedure :: set_imex_rhs
      procedure :: rotate_imex
      procedure :: assemble_implicit_buo
      procedure :: bridge_with_cnab2
      procedure :: start_with_ab1
   end type type_dirk

contains

   subroutine initialize(this, time_scheme, courfac_nml)

      class(type_dirk) :: this

      !-- Input/output variables
      real(cp),          intent(in) :: courfac_nml
      character(len=72), intent(inout) :: time_scheme

      !-- Local variables
      real(cp) :: courfac_loc

      allocate ( this%dt(1) )
      this%dt(:)=0.0_cp
      allocate ( this%wimp_lin(1) )
      this%wimp_lin(1)=0.0_cp

      this%family='DIRK'

      if ( index(time_scheme, 'ARS222') /= 0 ) then
         this%time_scheme = 'ARS222'
         this%norder_imp_lin = 3
         this%norder_imp = 3
         this%norder_exp = 2
         this%nstages = 2
         this%istage = 1
         courfac_loc = 1.35_cp
      else if ( index(time_scheme, 'LZ232') /= 0 ) then
         this%time_scheme = 'LZ232'
         this%norder_imp_lin = 3
         this%norder_imp = 3
         this%norder_exp = 2
         this%nstages = 2
         this%istage = 1
         courfac_loc = 1.25_cp
      else if ( index(time_scheme, 'ARS443') /= 0 ) then
         this%time_scheme = 'ARS443'
         this%norder_imp = 5
         this%norder_imp_lin = 5
         this%norder_exp = 4
         this%nstages = 4
         this%istage = 1
         courfac_loc = 0.9_cp
      else if ( index(time_scheme, 'LZ453') /= 0 ) then
         this%time_scheme = 'LZ453'
         this%norder_imp = 5
         this%norder_imp_lin = 5
         this%norder_exp = 4
         this%nstages = 4
         this%istage = 1
         courfac_loc = 1.15_cp
      else if ( index(time_scheme, 'BPR353') /= 0 ) then
         this%time_scheme = 'BPR353'
         this%norder_imp = 5
         this%norder_imp_lin = 5
         this%norder_exp = 4
         this%nstages = 4
         this%istage = 1
         courfac_loc = 1.0_cp
      else if ( index(time_scheme, 'PC2') /= 0 ) then
         this%time_scheme = 'PC2'
         this%norder_imp = 4
         this%norder_imp_lin = 4
         this%norder_exp = 3
         this%nstages = 3
         this%istage = 1
         courfac_loc = 0.8_cp
      end if

      if ( abs(courfac_nml) >= 1.0e3_cp ) then
         this%courfac=courfac_loc
      else
         this%courfac=courfac_nml
      end if

      allocate( this%butcher_imp(this%nstages+1,this%nstages+1), &
      &         this%butcher_exp(this%nstages+1,this%nstages+1) )
      this%butcher_imp(:,:)=0.0_cp
      this%butcher_exp(:,:)=0.0_cp
      bytes_allocated=bytes_allocated+2*(this%nstages+1)*(this%nstages+1)*&
      &               SIZEOF_DEF_REAL

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
      deallocate( this%l_exp_calc, this%l_imp_calc_rhs )

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine set_weights(this, lMatNext)

      class(type_dirk) :: this
      logical, intent(inout) :: lMatNext

      !-- Local variables
      real(cp) :: wimp_old, del, gam

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
      end select

      this%wimp_lin(1)      = this%dt(1)*this%wimp_lin(1)
      this%butcher_imp(:,:) = this%dt(1)*this%butcher_imp(:,:)
      this%butcher_exp(:,:) = this%dt(1)*this%butcher_exp(:,:)
         
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
   subroutine set_imex_rhs(this, rhs, dfdt, nMstart, nMstop, len_rhs)
      !
      ! This subroutine assembles the right-hand-side of an IMEX scheme
      !

      class(type_dirk) :: this

      !-- Input variables:
      integer,     intent(in) :: nMstart
      integer,     intent(in) :: nMstop
      integer,     intent(in) :: len_rhs
      type(type_tarray), intent(in) :: dfdt

      !-- Output variable
      complex(cp), intent(out) :: rhs(nMstart:nMstop,len_rhs)

      !-- Local variables
      integer :: n_stage, n_r, n_m

      do n_r=1,len_rhs
         do n_m=nMstart,nMstop
            rhs(n_m,n_r)=dfdt%old(n_m,n_r,1)
         end do
      end do

      do n_stage=1,this%istage
         do n_r=1,len_rhs
            do n_m=nMstart,nMstop
               rhs(n_m,n_r)=rhs(n_m,n_r)+this%butcher_exp(this%istage+1,n_stage)* &
               &            dfdt%expl(n_m,n_r,n_stage)
            end do
         end do
      end do

      do n_stage=1,this%istage
         do n_r=1,len_rhs
            do n_m=nMstart,nMstop
               rhs(n_m,n_r)=rhs(n_m,n_r)+this%butcher_imp(this%istage+1,n_stage)* &
               &                         dfdt%impl(n_m,n_r,n_stage)
            end do
         end do
      end do

   end subroutine set_imex_rhs
!------------------------------------------------------------------------------
   subroutine rotate_imex(this, dfdt, nMstart, nMstop, n_r_max)
      !
      ! This subroutine is used to roll the time arrays from one time step
      !

      class(type_dirk) :: this

      !-- Input variables:
      integer,     intent(in) :: nMstart
      integer,     intent(in) :: nMstop
      integer,     intent(in) :: n_r_max

      !-- Output variables:
      type(type_tarray), intent(inout) :: dfdt

   end subroutine rotate_imex
!------------------------------------------------------------------------------
   subroutine assemble_implicit_buo(this, buo, temp, dTdt, BuoFac, rgrav, &
              &                     nMstart, nMstop, n_r_max)
      !
      ! This subroutine is used to assemble Buoyancy
      !
      class(type_dirk) :: this

      !-- Input variables:
      integer,     intent(in) :: nMstart
      integer,     intent(in) :: nMstop
      integer,     intent(in) :: n_r_max
      real(cp),    intent(in) :: BuoFac
      real(cp),    intent(in) :: rgrav(n_r_max)
      complex(cp), intent(in) :: temp(nMstart:nMstop,n_r_max)
      type(type_tarray), intent(in) :: dTdt

      !-- Output variables:
      complex(cp), intent(out) :: buo(nMstart:nMstop,n_r_max)

      !-- Local variables:
      integer :: n_stage, n_m, n_r, m

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            buo(n_m,n_r)=zero
         end do
      end do

      do n_stage=1,this%istage
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               if ( m /= 0 ) then
                  buo(n_m,n_r)=buo(n_m,n_r)-                              &
                  &            this%butcher_imp(this%istage+1,n_stage)*   &
                  &            rgrav(n_r)*or1(n_r)*BuoFac*ci*real(m,cp)*  &
                  &            dTdt%old(n_m,n_r,n_stage)
               end if
            end do
         end do
      end do

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m /= 0 ) then
               buo(n_m,n_r)=buo(n_m,n_r)-this%wimp_lin(1)*rgrav(n_r)*  &
               &            or1(n_r)*BuoFac*ci*real(m,cp)*temp(n_m,n_r)
            end if
         end do
      end do

   end subroutine assemble_implicit_buo
!------------------------------------------------------------------------------
   subroutine bridge_with_cnab2(this)

      class(type_dirk) :: this

   end subroutine bridge_with_cnab2
!------------------------------------------------------------------------------
   subroutine start_with_ab1(this)

      class(type_dirk) :: this

   end subroutine start_with_ab1
!------------------------------------------------------------------------------
end module dirk_schemes
