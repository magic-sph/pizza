module time_schemes

   use time_array
   use precision_mod

   implicit none

   private 

   type, abstract, public :: type_tscheme

      character(len=10) :: family
      integer :: nstages
      integer :: istage
      integer :: nexp
      integer :: nold
      integer :: nimp
      character(len=8) :: time_scheme
      logical :: l_assembly
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
      procedure(set_imex_rhs_ghost_if),  deferred :: set_imex_rhs_ghost
      procedure(rotate_imex_if), deferred :: rotate_imex
      procedure(bridge_with_cnab2_if), deferred :: bridge_with_cnab2
      procedure(start_with_ab1_if), deferred :: start_with_ab1
      procedure(assemble_imex_if), deferred :: assemble_imex
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

      subroutine set_imex_rhs_if(this, rhs, dfdt)
         import
         class(type_tscheme) :: this
         type(type_tarray), intent(in) :: dfdt
         complex(cp), intent(out) :: rhs(dfdt%lm:dfdt%um,dfdt%nRstart:dfdt%nRstop)
      end subroutine set_imex_rhs_if

      subroutine set_imex_rhs_ghost_if(this, rhs, dfdt, start_m, stop_m, ng)
         import
         class(type_tscheme) :: this
         type(type_tarray), intent(in) :: dfdt
         integer,           intent(in) :: start_m ! Starting index
         integer,           intent(in) :: stop_m  ! Stopping index
         integer,           intent(in) :: ng      ! Number of ghosts zones
         complex(cp), intent(out) :: rhs(dfdt%lm:dfdt%um,dfdt%nRstart-ng:dfdt%nRstop+ng)
      end subroutine set_imex_rhs_ghost_if

      subroutine rotate_imex_if(this, dfdt)
         import
         class(type_tscheme) :: this
         type(type_tarray), intent(inout) :: dfdt
      end subroutine rotate_imex_if

      subroutine assemble_imex_if(this, rhs, dfdt)
         import
         class(type_tscheme) :: this
         type(type_tarray), intent(in) :: dfdt
         complex(cp), intent(out) :: rhs(dfdt%lm:dfdt%um,dfdt%nRstart:dfdt%nRstop)
      end subroutine assemble_imex_if

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
