module time_schemes

   use precision_mod
   use parallel_mod
   use namelists, only: alpha, l_cheb_coll
   use constants, only: one, half, two
   use mem_alloc, only: bytes_allocated
   use char_manip, only: capitalize
   use useful, only: abortRun
   use time_array

   implicit none

   private

   type, public :: type_tscheme
      integer :: norder_exp
      integer :: norder_imp
      integer :: norder_imp_lin
      character(len=8) :: time_scheme
      real(cp), allocatable :: dt(:)
      real(cp), allocatable :: wimp(:)
      real(cp), allocatable :: wimp_lin(:)
      real(cp), allocatable :: wexp(:)
      logical :: l_calc_lin_rhs ! Do we need to calculate the linear op in the past
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_weights
      procedure :: set_dt_array
      procedure :: set_imex_rhs
      procedure :: rotate_imex
   end type type_tscheme

contains

   subroutine initialize(this, time_scheme)

      class(type_tscheme) :: this

      !-- Input/output variables
      character(len=72), intent(inout) :: time_scheme

      call capitalize(time_scheme)

      if ( index(time_scheme, 'CNAB2') /= 0 ) then
         this%time_scheme = 'CNAB2'
         this%norder_imp_lin = 2
         this%norder_imp = 2
         this%norder_exp = 2
         this%l_calc_lin_rhs = .true.
      else if ( index(time_scheme, 'MODCNAB') /= 0 ) then
         this%time_scheme = 'MODCNAB'
         this%norder_imp = 3
         this%norder_imp_lin = 3
         this%norder_exp = 2
         this%l_calc_lin_rhs = .true.
      else if ( index(time_scheme, 'CNLF') /= 0 ) then
         this%time_scheme = 'CNLF'
         this%norder_imp = 3
         this%norder_imp_lin = 3
         this%norder_exp = 2
         this%l_calc_lin_rhs = .true.
      else if ( index(time_scheme, 'BDF2AB2') /= 0 ) then
         this%time_scheme = 'BDF2AB2'
         this%norder_imp = 3
         this%norder_imp_lin = 2 ! it should be one but we need to restart
         this%norder_exp = 2
         this%l_calc_lin_rhs = .false.
      else if ( index(time_scheme, 'BDF3AB3') /= 0 ) then
         this%time_scheme = 'BDF3AB3'
         this%norder_imp = 4
         this%norder_imp_lin = 2 ! it should be one but we need to restart
         this%norder_exp = 3
         this%l_calc_lin_rhs = .false.
      else if ( index(time_scheme, 'BDF4AB4') /= 0 ) then
         this%time_scheme = 'BDF4AB4'
         this%norder_imp = 5
         this%norder_imp_lin = 2 ! it should be one but we need to restart
         this%norder_exp = 4
         this%l_calc_lin_rhs = .false.
      end if

      if ( (this%time_scheme == 'CNLF'.or. this%time_scheme=='MODCNAB') .and. &
      &    (.not. l_cheb_coll) ) then
           call abortRun('! Integration method does not work with the chosen &
           &              time scheme. It would require to store buoyancy as well')
      end if

      allocate ( this%dt(this%norder_exp) )
      allocate ( this%wimp(this%norder_imp) )
      allocate ( this%wimp_lin(this%norder_imp_lin) )
      allocate ( this%wexp(this%norder_exp) )

      this%dt(:)       = 0.0_cp
      this%wimp(:)     = 0.0_cp
      this%wimp_lin(:) = 0.0_cp
      this%wexp(:)     = 0.0_cp

      bytes_allocated = bytes_allocated+(2*this%norder_exp+this%norder_imp+&
      &                 this%norder_imp_lin)*SIZEOF_DEF_REAL

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)

      class(type_tscheme) :: this

      deallocate( this%dt, this%wimp, this%wimp_lin, this%wexp )

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine set_weights(this)

      class(type_tscheme) :: this

      real(cp) :: delta, delta_n, delta_n_1, delta_n_2
      real(cp) :: a0, a1, a2, a3, a4, b0, b1, b2, b3, c1, c2, c3

      select case ( this%time_scheme )
         case ('CNAB2') 
            this%wimp(1)    =one
            this%wimp(2)    =one
            this%wimp_lin(1)=alpha*this%dt(1)
            this%wimp_lin(2)=(1-alpha)*this%dt(1)

            this%wexp(1)=(one+half*this%dt(1)/this%dt(2))*this%dt(1)
            this%wexp(2)=-half*this%dt(1)*this%dt(1)/this%dt(2)
         case ('CNLF')
            delta = this%dt(1)/this%dt(2)
            this%wimp(1)    =one
            this%wimp(2)    =(one-delta)*(one+delta)
            this%wimp(3)    =delta*delta
            this%wimp_lin(1)=half*(one+delta)/delta*this%dt(1)
            this%wimp_lin(2)=half*(one+delta)*(delta-one)/delta*this%dt(1)
            this%wimp_lin(3)=half*(one+delta)*this%dt(1)

            this%wexp(1)=(one+delta)*this%dt(1)
            this%wexp(2)=0.0_cp
         case ('MODCNAB') 
            delta = this%dt(1)/this%dt(2)
            this%wimp(1)    =one
            this%wimp(2)    =one
            this%wimp(3)    =0.0_cp
            this%wimp_lin(1)=(half+1.0_cp/delta/16.0_cp)*this%dt(1)
            this%wimp_lin(2)=(7.0_cp/16.0_cp-1.0_cp/delta/16.0_cp)*this%dt(1)
            this%wimp_lin(3)=1.0_cp/16.0_cp*this%dt(1)

            this%wexp(1)=(one+half*delta)*this%dt(1)
            this%wexp(2)=-half*delta*this%dt(1)
         case ('BDF2AB2')
            delta = this%dt(1)/this%dt(2)
            this%wimp_lin(1)=(one+delta)/(one+two*delta)*this%dt(1)
            this%wimp_lin(2)=0.0_cp
            this%wimp(1)=one
            this%wimp(2)=(one+delta)*(one+delta)/(one+two*delta)
            this%wimp(3)=-delta*delta/(one+two*delta)

            this%wexp(1)=(one+delta)*(one+delta)*this%dt(1)/(one+two*delta)
            this%wexp(2)=-delta*(one+delta)*this%dt(1)/(one+two*delta)
         case ('BDF3AB3')
            delta_n   = this%dt(2)/this%dt(1)
            delta_n_1 = this%dt(3)/this%dt(1)
            a0 = one+one/(one+delta_n)+one/(one+delta_n+delta_n_1)
            a1 = (one+delta_n)*(one+delta_n+delta_n_1)/ &
            &    (delta_n*(delta_n+delta_n_1))
            a2 = -(one+delta_n+delta_n_1)/(delta_n*delta_n_1* &
            &    (one+delta_n))
            a3 = (one+delta_n)/(delta_n_1*(delta_n+delta_n_1)* &
            &     (one+delta_n+delta_n_1))
            b0 = (one+delta_n)*(one+delta_n+delta_n_1)/(delta_n* &
            &    (delta_n+delta_n_1))
            b1 = -(one+delta_n+delta_n_1)/(delta_n*delta_n_1)
            b2 = (one+delta_n)/(delta_n_1*(delta_n+delta_n_1))

            this%wimp_lin(1)=one/a0 * this%dt(1)
            this%wimp_lin(2)=0.0_cp

            this%wimp(1)=one
            this%wimp(2)=a1/a0
            this%wimp(3)=a2/a0
            this%wimp(4)=a3/a0

            this%wexp(1)=b0/a0 * this%dt(1)
            this%wexp(2)=b1/a0 * this%dt(1)
            this%wexp(3)=b2/a0 * this%dt(1)
         case ('BDF4AB4')
            delta_n = this%dt(1)/this%dt(2)
            delta_n_1 = this%dt(2)/this%dt(3)
            delta_n_2 = this%dt(3)/this%dt(4)

            c1 = one+delta_n_2*(one+delta_n_1)
            c2 = one+delta_n_1*(one+delta_n)
            c3 = one+delta_n_2*c2

            a4 = (one+delta_n)/(one+delta_n_2) * c2/c1/c3 *(delta_n_2**4* &
            &    delta_n_1**3*delta_n**2)
            a3 = -delta_n_1**3*delta_n**2*(one+delta_n)/(one+delta_n_1) * c3/c2
            a2 = delta_n*(delta_n/(one+delta_n)+delta_n_1*delta_n * (c3+delta_n_2)/&
            &    (one+delta_n_2))
            a1 = -one-delta_n*(one+delta_n_1*(one+delta_n)*(one+delta_n_2*c2/c1)/ &
            &    (one+delta_n_1))
            a0 = one + delta_n/(one+delta_n)+delta_n_1*delta_n/c2+delta_n_2* &
            &    delta_n_1*delta_n/c3

            b3 = -delta_n_2**3*delta_n_1**2*delta_n*(one+delta_n)/(one+delta_n_2)* &
            &    c2/c1
            b2 = delta_n_1**2*delta_n*(one+delta_n)/(one+delta_n_1)*c3
            b1 = -c2*c3 * delta_n/(one+delta_n_2)
            b0 = delta_n_1*(one+delta_n)/(one+delta_n_1) * ((one+delta_n) * &
            &    (c3+delta_n_2)+(one+delta_n_2)/delta_n_1)/c1

            this%wimp_lin(1)=one/a0 * this%dt(1)
            this%wimp_lin(2)=0.0_cp

            this%wimp(1)=one
            this%wimp(2)=-a1/a0
            this%wimp(3)=-a2/a0
            this%wimp(4)=-a3/a0
            this%wimp(5)=-a4/a0

            this%wexp(1)=b0/a0 * this%dt(1)
            this%wexp(2)=b1/a0 * this%dt(1)
            this%wexp(3)=b2/a0 * this%dt(1)
            this%wexp(4)=b3/a0 * this%dt(1)
      end select

   end subroutine set_weights
!------------------------------------------------------------------------------
   subroutine set_dt_array(this, dt_new, dt_min, time, n_log_file,  &
              &            n_time_step, l_new_dtNext)
      !
      ! This subroutine adjusts the time step
      !

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

      class(type_tscheme) :: this

      !-- Input variables:
      integer,     intent(in) :: nMstart
      integer,     intent(in) :: nMstop
      integer,     intent(in) :: len_rhs
      type(type_tarray), intent(in) :: dfdt

      !-- Output variable
      complex(cp), intent(out) :: rhs(nMstart:nMstop,len_rhs)

      !-- Local variables
      integer :: n_o, n_r, n_m

      do n_o=1,this%norder_imp-1
         if ( n_o == 1 ) then
            do n_r=1,len_rhs
               do n_m=nMstart,nMstop
                  rhs(n_m,n_r)=this%wimp(n_o+1)*dfdt%old(n_m,n_r,n_o)
               end do
            end do
         else
            do n_r=1,len_rhs
               do n_m=nMstart,nMstop
                  rhs(n_m,n_r)=rhs(n_m,n_r)+this%wimp(n_o+1)*dfdt%old(n_m,n_r,n_o)
               end do
            end do
         end if
      end do

      do n_o=1,this%norder_imp_lin-1
         do n_r=1,len_rhs
            do n_m=nMstart,nMstop
               rhs(n_m,n_r)=rhs(n_m,n_r)+this%wimp_lin(n_o+1)*dfdt%impl(n_m,n_r,n_o)
            end do
         end do
      end do

      do n_o=1,this%norder_exp
         do n_r=1,len_rhs
            do n_m=nMstart,nMstop
               rhs(n_m,n_r)=rhs(n_m,n_r)+this%wexp(n_o)*dfdt%expl(n_m,n_r,n_o)
            end do
         end do
      end do

   end subroutine set_imex_rhs
!------------------------------------------------------------------------------
   subroutine rotate_imex(this, dfdt, nMstart, nMstop, n_r_max)
      !
      ! This subroutine is used to roll the time arrays from one time step
      !

      class(type_tscheme) :: this

      !-- Input variables:
      integer,     intent(in) :: nMstart
      integer,     intent(in) :: nMstop
      integer,     intent(in) :: n_r_max

      !-- Output variables:
      type(type_tarray), intent(inout) :: dfdt

      !-- Local variables:
      integer :: n_o, n_m, n_r

      do n_o=this%norder_exp,2,-1
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               dfdt%expl(n_m,n_r,n_o)=dfdt%expl(n_m,n_r,n_o-1)
            end do
         end do
      end do

      do n_o=this%norder_imp-1,2,-1
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               dfdt%old(n_m,n_r,n_o)=dfdt%old(n_m,n_r,n_o-1)
            end do
         end do
      end do

      do n_o=this%norder_imp_lin-1,2,-1
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               dfdt%impl(n_m,n_r,n_o)=dfdt%impl(n_m,n_r,n_o-1)
            end do
         end do
      end do

   end subroutine rotate_imex
!------------------------------------------------------------------------------
end module time_schemes
