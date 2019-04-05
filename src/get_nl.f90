module general_arrays_mod
 
   implicit none
 
   private
 
   type, public, abstract :: general_arrays_t
 
   end type general_arrays_t
 
end module general_arrays_mod
!----------------------------------------------------------------------------
module grid_space_arrays_mod

   use general_arrays_mod
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation_3D, only: n_theta_max, n_phi_max_3D
   use radial_functions, only: or1_3D, r_3D
   use horizontal, only: osint1

   implicit none

   private

   type, public, extends(general_arrays_t) :: grid_space_arrays_t
      !----- Nonlinear terms in phi/theta space: 
      real(cp), allocatable :: VTr(:,:), VTt(:,:), VTp(:,:)

      !----- Fields calculated from these help arrays by legtf:
      real(cp), pointer :: Tc(:,:)
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: get_nl
   end type grid_space_arrays_t

contains

   subroutine initialize(this)

      class(grid_space_arrays_t) :: this

      allocate( this%VTr(n_phi_max_3D,n_theta_max) )
      allocate( this%VTt(n_phi_max_3D,n_theta_max) )
      allocate( this%VTp(n_phi_max_3D,n_theta_max) )

      !----- Fields calculated from these help arrays by legtf:
      allocate( this%Tc(n_phi_max_3D,n_theta_max) )
      bytes_allocated=bytes_allocated + 4*n_phi_max_3D*n_theta_max*SIZEOF_DEF_REAL

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(grid_space_arrays_t) :: this

      deallocate( this%VTr, this%VTt, this%VTp )
      deallocate( this%Tc )

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine get_nl(this, vr, vt, vp, n_r)
      !
      !  calculates non-linear products in grid-space for radial
      !  level n_r and returns them in arrays wnlr1-3, snlr1-3, bnlr1-3
      !
      !  temp_r: (input) velocity, and temperature on grid points
      !  n_r: (input) radial level
      !

      class(grid_space_arrays_t) :: this

      !-- Input of variables:
      real(cp), intent(in) :: vr(n_phi_max_3D,n_theta_max)
      real(cp), intent(in) :: vt(n_phi_max_3D,n_theta_max)
      real(cp), intent(in) :: vp(n_phi_max_3D,n_theta_max)
      integer,  intent(in) :: n_r

      !-- Local variables:
      integer :: n_theta
      integer :: n_phi
      real(cp) :: or1sn1, r2

      r2 = r_3D(n_r)*r_3D(n_r)

      !------ Get V T, the divergence of it is temperature advection:
      !$OMP PARALLEL DO default(shared) &
      !$OMP& private(n_theta, n_phi, or1sn1)
      do n_theta=1,n_theta_max
         or1sn1=or1_3D(n_r)*osint1(n_theta)
         do n_phi=1,n_phi_max_3D     ! calculate u*T components
            this%VTr(n_phi,n_theta) = &
            &    r2*vr(n_phi,n_theta)*this%Tc(n_phi,n_theta)
            this%VTt(n_phi,n_theta) = &
            &    or1sn1*vt(n_phi,n_theta)*this%Tc(n_phi,n_theta)
            this%VTp(n_phi,n_theta) = &
            &    or1sn1*vp(n_phi,n_theta)*this%Tc(n_phi,n_theta)
         end do
      end do  ! theta loop
      !$OMP END PARALLEL DO

   end subroutine get_nl
!----------------------------------------------------------------------------
end module grid_space_arrays_mod
