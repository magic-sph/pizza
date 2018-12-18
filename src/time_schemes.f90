module time_schemes

   use time_array
   use precision_mod

   implicit none

   private 

   type, abstract, public :: type_tscheme

      character(len=10) :: family
      integer :: nstages
      integer :: istage
      integer :: norder_exp
      integer :: norder_imp
      integer :: norder_imp_lin
      character(len=8) :: time_scheme
      real(cp), allocatable :: dt(:)
      real(cp), allocatable :: wimp_lin(:)
      logical,  allocatable :: l_exp_calc(:)
      logical, allocatable :: l_imp_calc_rhs(:)
      real(cp) :: courfac ! Courant factor

   contains 

      procedure(initialize_if), deferred :: initialize
      procedure(finalize_if),  deferred :: finalize
      procedure(set_weights_if), deferred :: set_weights
      procedure(set_dt_array_if), deferred :: set_dt_array
      procedure(set_imex_rhs_if),  deferred :: set_imex_rhs
      procedure(rotate_imex_if), deferred :: rotate_imex
      procedure(assemble_implicit_buo_if), deferred :: assemble_implicit_buo
      procedure(bridge_with_cnab2_if), deferred :: bridge_with_cnab2
      procedure(start_with_ab1_if), deferred :: start_with_ab1
      procedure :: print_info

   end type type_tscheme

   interface

      subroutine initialize_if(this, time_scheme, courfac_nml)
         import
         class(type_tscheme) :: this
         real(cp),          intent(in)    :: courfac_nml
         character(len=72), intent(inout) :: time_scheme
      end subroutine initialize_if

      subroutine finalize_if(this)
         import
         class(type_tscheme) :: this
      end subroutine finalize_if

      subroutine set_weights_if(this, lMatNext)
         import
         class(type_tscheme) :: this
         logical, intent(inout) :: lMatNext
      end subroutine set_weights_if

      subroutine set_dt_array_if(this, dt_new, dt_min, time, n_log_file,  &
                 &                n_time_step, l_new_dtNext)
         import
         class(type_tscheme) :: this
         real(cp), intent(in) :: dt_new
         real(cp), intent(in) :: dt_min
         real(cp), intent(in) :: time
         integer,  intent(in) :: n_log_file
         integer,  intent(in) :: n_time_step
         logical,  intent(in) :: l_new_dtNext
      end subroutine set_dt_array_if

      subroutine set_imex_rhs_if(this, rhs, dfdt, nMstart, nMstop, len_rhs)
         import
         class(type_tscheme) :: this
         integer,     intent(in) :: nMstart
         integer,     intent(in) :: nMstop
         integer,     intent(in) :: len_rhs
         type(type_tarray), intent(in) :: dfdt
         complex(cp), intent(out) :: rhs(nMstart:nMstop,len_rhs)
      end subroutine set_imex_rhs_if

      subroutine rotate_imex_if(this, dfdt, nMstart, nMstop, n_r_max)
         import
         class(type_tscheme) :: this
         integer,     intent(in) :: nMstart
         integer,     intent(in) :: nMstop
         integer,     intent(in) :: n_r_max
         type(type_tarray), intent(inout) :: dfdt
      end subroutine rotate_imex_if

      subroutine assemble_implicit_buo_if(this, buo, temp, dTdt, BuoFac, rgrav, &
                 &                     nMstart, nMstop, n_r_max, l_init_buo)
         import
         class(type_tscheme) :: this
         integer,           intent(in) :: nMstart
         integer,           intent(in) :: nMstop
         integer,           intent(in) :: n_r_max
         real(cp),          intent(in) :: BuoFac
         real(cp),          intent(in) :: rgrav(n_r_max)
         complex(cp),       intent(in) :: temp(nMstart:nMstop,n_r_max)
         type(type_tarray), intent(in) :: dTdt
         complex(cp),       intent(out) :: buo(nMstart:nMstop,n_r_max)
         logical,           intent(in) :: l_init_buo
      end subroutine assemble_implicit_buo_if

      subroutine bridge_with_cnab2_if(this)
         import
         class(type_tscheme) :: this
      end subroutine bridge_with_cnab2_if

      subroutine start_with_ab1_if(this)
         import
         class(type_tscheme) :: this
      end subroutine start_with_ab1_if

   end interface

contains

      subroutine print_info(this, n_log_file)

         class(type_tscheme) :: this

         integer, intent(in) :: n_log_file

         !-- Local variables
         integer :: n, n_out

         do n=1,2
            if ( n == 1 ) n_out=6
            if ( n == 2 ) n_out=n_log_file
            write(n_out,*) ''
            write(n_out, '('' ! Time integrator  :'',1p,A10)') this%time_scheme
            write(n_out, '('' ! CFL value        :'',es12.4)') this%courfac
            write(n_out,*) ''
         end do

      end subroutine print_info
!------------------------------------------------------------------------------
end module time_schemes
