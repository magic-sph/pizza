module time_scheme

   use precision_mod
   use parallel_mod
   use namelists, only: alpha
   use constants, only: one, half, two
   use mem_alloc, only: bytes_allocated
   use char_manip, only: capitalize
   use useful, only: abortRun

   implicit none

   private

   type, public :: type_tscheme
      integer :: norder_exp
      integer :: norder_imp
      character(len=4) :: imp_scheme
      character(len=4) :: exp_scheme
      real(cp), allocatable :: dt(:)
      real(cp), allocatable :: wimp(:)
      real(cp), allocatable :: wimp_lin(:)
      real(cp), allocatable :: wimp_buo(:)
      real(cp), allocatable :: wexp(:)

   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_weights
      procedure :: set_dt_array
   end type type_tscheme

contains

   subroutine initialize(this, imp_scheme, exp_scheme)

      class(type_tscheme) :: this

      !-- Input/output variables
      character(len=72), intent(inout) :: imp_scheme
      character(len=72), intent(inout) :: exp_scheme

      call capitalize(imp_scheme)

      if ( index(imp_scheme, 'CN') /= 0 ) then
         this%imp_scheme = 'CN'
         this%norder_imp = 2
      else if ( index(imp_scheme, 'BDF2') /= 0 ) then
         this%imp_scheme = 'BDF2'
         this%norder_imp = 3
      end if

      if ( index(exp_scheme, 'AB2') /= 0 ) then
         this%exp_scheme = 'AB2'
         this%norder_exp = 2
      else if ( index(exp_scheme, 'AB3') /= 0 ) then
         this%exp_scheme = 'AB3'
         this%norder_exp = 3
      end if


      allocate ( this%dt(this%norder_exp) )
      allocate ( this%wimp(this%norder_imp) )
      allocate ( this%wimp_lin(this%norder_imp) )
      allocate ( this%wimp_buo(this%norder_imp) )
      allocate ( this%wexp(this%norder_exp) )

      bytes_allocated = bytes_allocated+(2*this%norder_exp+3*this%norder_imp)* &
      &                 SIZEOF_DEF_REAL

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)

      class(type_tscheme) :: this

      deallocate( this%dt, this%wimp, this%wexp )

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine set_weights(this)

      class(type_tscheme) :: this

      real(cp) :: delta

      select case ( this%imp_scheme )
         case ('CN') 
            this%wimp(1)    =one
            this%wimp(2)    =one
            this%wimp_lin(1)=alpha*this%dt(1)
            this%wimp_lin(2)=(1-alpha)*this%dt(1)
            this%wimp_buo(1)=alpha*this%dt(1)
            this%wimp_buo(2)=(1-alpha)*this%dt(1)

            select case ( this%exp_scheme )
               case ('AB2')
               this%wexp(2)=-half*this%dt(1)/this%dt(2)
               this%wexp(1)=one-this%wexp(2)
            end select
         case ('BDF2')
            delta = this%dt(1)/this%dt(2)
            this%wimp_lin(1)=(one+delta)/(one+two*delta)*this%dt(1)
            this%wimp_lin(2)=0.0_cp
            this%wimp_lin(3)=0.0_cp
            this%wimp_buo(1)=(one+delta)/(one+two*delta)*this%dt(1)
            this%wimp_buo(2)=0.0_cp
            this%wimp_buo(3)=0.0_cp
            !this%wimp_buo(2)=-(one+delta)*(one+delta)/(one+two*delta)*this%dt(1)
            !this%wimp_buo(3)=delta*delta/(one+two*delta)*this%dt(1)
            this%wimp(1)=one
            this%wimp(2)=(one+delta)*(one+delta)/(one+two*delta)
            this%wimp(3)=-delta*delta/(one+two*delta)
            select case ( this%exp_scheme )
               case ('AB2')
                  this%wexp(1)=(one+delta)*(one+delta)/(one+two*delta)
                  this%wexp(2)=-delta*(one+delta)/(one+two*delta)
            end select
      end select

   end subroutine set_weights
!------------------------------------------------------------------------------
   subroutine set_dt_array(this, dt_new, dt_min, time, n_log_file,  &
              &            n_time_step, l_new_dtNext)

      class(type_tscheme) :: this

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

      !-- First roll the dt array
      this%dt   =cshift(this%dt,shift=this%norder_exp-1)
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
            write(*,'(1p,/,A,ES18.10,/,A,i9,/,A,ES15.8,/,A,ES15.8)')    &
            &    " ! Changing time step at time=",(time+this%dt(1)),    &
            &    "                 time step no=",n_time_step,          &
            &    "                      last dt=",dt_old,               &
            &    "                       new dt=",dt_new
            write(n_log_file,                                           &
            &    '(1p,/,A,ES18.10,/,A,i9,/,A,ES15.8,/,A,ES15.8)')       &
            &    " ! Changing time step at time=",(time+this%dt(1)),    &
            &    "                 time step no=",n_time_step,          &
            &    "                      last dt=",dt_old,               &
            &    "                       new dt=",dt_new
         end if
      end if

   end subroutine set_dt_array
!------------------------------------------------------------------------------
end module time_scheme
