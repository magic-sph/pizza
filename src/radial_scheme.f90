module radial_scheme
   !
   ! This is an abstract type that defines the radial scheme
   !

   use precision_mod

   implicit none

   private

   type, abstract, public :: type_rscheme

      integer :: n_max
      integer :: order
      integer :: order_boundary
      integer :: nRmax
      real(cp) :: rnorm
      real(cp) :: boundary_fac
      real(cp), allocatable :: rMat(:,:)
      real(cp), allocatable :: drMat(:,:)
      real(cp), allocatable :: d2rMat(:,:)
      real(cp), allocatable :: dr(:,:)
      real(cp), allocatable :: ddr(:,:)
      real(cp), allocatable :: dddr(:,:)
      real(cp), allocatable :: dr_top(:,:)
      real(cp), allocatable :: dr_bot(:,:)
      real(cp), allocatable :: ddr_top(:,:)
      real(cp), allocatable :: ddr_bot(:,:)
      real(cp), allocatable :: dddr_top(:,:)
      real(cp), allocatable :: dddr_bot(:,:)
      real(cp), allocatable :: drx(:)   ! First derivative of non-linear mapping (see Bayliss and Turkel, 1990)
      real(cp), allocatable :: ddrx(:)  ! Second derivative of non-linear mapping

      character(len=72) :: version

   contains

      procedure(initialize_if),  deferred :: initialize
      procedure(empty_if),       deferred :: finalize
      procedure(get_der_mat_if), deferred :: get_der_mat
      procedure(get_grid_if),    deferred :: get_grid
      procedure :: costf1_complex_2d
      procedure :: costf1_real_1d
      generic :: costf1 => costf1_complex_2d, costf1_real_1d

   end type type_rscheme

   interface 

      subroutine empty_if(this, no_work_array)
         import
         class(type_rscheme) :: this
         logical, optional, intent(in) :: no_work_array
      end subroutine empty_if
      !------------------------------------------------------------------------
      subroutine get_grid_if(this,n_r_max,ricb,rcmb,ratio1,ratio2,r)
         import
         class(type_rscheme) :: this

         !-- Input quantities:
         integer,  intent(in) :: n_r_max    ! Number of grid points
         real(cp), intent(inout) :: ratio1  ! Nboudary/Nbulk
         real(cp), intent(in) :: ratio2     ! drMin/drMax
         real(cp), intent(in) :: ricb       ! inner boundary
         real(cp), intent(in) :: rcmb       ! outer boundary

         !-- Output quantities:
         real(cp), intent(out) :: r(n_r_max) ! radius

      end subroutine get_grid_if
      !------------------------------------------------------------------------
      subroutine initialize_if(this,n_r_max,order,order_boundary,l_cheb_coll, &
                 &             no_work_array)

         import
         class(type_rscheme) :: this
         integer, intent(in) :: n_r_max
         integer, intent(in) :: order
         integer, intent(in) :: order_boundary
         logical, intent(in) :: l_cheb_coll
         logical, optional, intent(in) :: no_work_array

      end subroutine initialize_if
      !------------------------------------------------------------------------
      subroutine get_der_mat_if(this,n_r_max,l_cheb_coll)

         import
         class(type_rscheme) :: this
         integer, intent(in) :: n_r_max
         logical, intent(in) :: l_cheb_coll

      end subroutine get_der_mat_if
      !------------------------------------------------------------------------
   end interface

contains

   subroutine costf1_complex_2d(this,f,nMstart,nMstop,n_r_max)

      class(type_rscheme) :: this

      !-- Input variables:
      integer,  intent(in) :: nMstart,nMstop
      integer,  intent(in) :: n_r_max

      !-- Output variables:
      complex(cp), intent(inout) :: f(nMstart:nMstop,n_r_max)

   end subroutine costf1_complex_2d
!------------------------------------------------------------------------------
   subroutine costf1_real_1d(this,f,n_r_max)

      class(type_rscheme) :: this

      !-- Input variable
      integer, intent(in) :: n_r_max

      !-- Output variable
      real(cp), intent(inout) :: f(n_r_max)

   end subroutine costf1_real_1d
!------------------------------------------------------------------------------
end module radial_scheme
