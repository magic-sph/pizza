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
   use constants, only: zero, two, pi
   !use namelists, only: DyMagFac
   use truncation_3D, only: n_r_max_3D, n_theta_max, n_phi_max_3D
   use radial_functions, only: or1_3D, or2_3D, r_3D, rgrav_3D
   use horizontal, only: cost, sint, osint1

   implicit none

   private

   type, public, extends(general_arrays_t) :: grid_space_arrays_t
      !----- Nonlinear terms in phi/theta space: 
      real(cp), allocatable :: VTr(:,:), VTt(:,:), VTp(:,:)
      real(cp), allocatable :: VxBr(:,:), VxBt(:,:), VxBp(:,:)
      real(cp), allocatable :: jxBr(:,:), jxBt(:,:), jxBp(:,:)

      !----- Fields calculated from these help arrays by legtf:
      real(cp), pointer :: Tc(:,:)
      real(cp), pointer :: Brc(:,:), Btc(:,:), Bpc(:,:)
      real(cp), pointer :: curlBrc(:,:), curlBtc(:,:), curlBpc(:,:)
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
      allocate( this%VxBr(n_phi_max_3D,n_theta_max) )
      allocate( this%VxBt(n_phi_max_3D,n_theta_max) )
      allocate( this%VxBp(n_phi_max_3D,n_theta_max) )
      allocate( this%jxBr(n_phi_max_3D,n_theta_max) )
      allocate( this%jxBt(n_phi_max_3D,n_theta_max) )
      allocate( this%jxBp(n_phi_max_3D,n_theta_max) )


      !----- Fields calculated from these help arrays by legtf:
      allocate( this%Tc(n_phi_max_3D,n_theta_max) )
      allocate( this%Brc(n_phi_max_3D,n_theta_max) )
      allocate( this%Btc(n_phi_max_3D,n_theta_max) )
      allocate( this%Bpc(n_phi_max_3D,n_theta_max) )
      allocate( this%curlBrc(n_phi_max_3D,n_theta_max) )
      allocate( this%curlBtc(n_phi_max_3D,n_theta_max) )
      allocate( this%curlBpc(n_phi_max_3D,n_theta_max) )
      bytes_allocated=bytes_allocated + 16*n_phi_max_3D*n_theta_max*SIZEOF_DEF_REAL

      this%VTr(:,:)=zero
      this%VTt(:,:)=zero
      this%VTp(:,:)=zero
      this%VxBr(:,:)=zero
      this%VxBt(:,:)=zero
      this%VxBp(:,:)=zero
      this%jxBr(:,:)=zero
      this%jxBt(:,:)=zero
      this%jxBp(:,:)=zero
      this%Tc(:,:)=zero
      this%Brc(:,:)=zero
      this%Btc(:,:)=zero
      this%Bpc(:,:)=zero
      this%curlBrc(:,:)=zero
      this%curlBtc(:,:)=zero
      this%curlBpc(:,:)=zero

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(grid_space_arrays_t) :: this

      deallocate( this%VTr, this%VTt, this%VTp )
      deallocate( this%VxBr, this%VxBt, this%VxBp )
      deallocate( this%jxBr, this%jxBt, this%jxBp )
      deallocate( this%Tc )
      deallocate( this%Brc, this%Btc, this%Bpc )
      deallocate( this%curlBrc, this%curlBtc, this%curlBpc )

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine get_nl(this, vr, vt, vp, n_r, buo, lfs)
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

      !-- Output of variable:
      real(cp), intent(out) :: buo(n_phi_max_3D,n_theta_max)
      real(cp), intent(out) :: lfs(n_phi_max_3D,n_theta_max)

      !-- Local variables:
      integer :: n_theta, n_phi
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

      !$OMP PARALLEL DO default(shared) &
      !$OMP& private(n_theta, n_phi)
      !-- Assemble g/r * T on the spherical grid
      do n_theta=1,n_theta_max
         do n_phi=1,n_phi_max_3D
            buo(n_phi,n_theta)=this%Tc(n_phi,n_theta)*rgrav_3D(n_r)* &
            &                  or1_3D(n_r)
         end do
      end do
      !$OMP END PARALLEL DO

      !------ Get (V x B) , the curl of this is the dynamo term:
      !$OMP PARALLEL DO default(shared) &
      !$OMP& private(n_theta, n_phi, or1sn1)
      do n_theta=1,n_theta_max
         or1sn1=or1_3D(n_r)*osint1(n_theta)!sint(n_theta)!
         do n_phi=1,n_phi_max_3D     ! calculate V x B components
            this%VxBr(n_phi,n_theta)=r2*(                   &
            &    vt(n_phi,n_theta)*this%Bpc(n_phi,n_theta)- &
            &    vp(n_phi,n_theta)*this%Btc(n_phi,n_theta) )

            this%VxBt(n_phi,n_theta)=or1sn1*(               &
            &    vp(n_phi,n_theta)*this%Brc(n_phi,n_theta)- &
            &    vr(n_phi,n_theta)*this%Bpc(n_phi,n_theta) )

            this%VxBp(n_phi,n_theta)=or1sn1*(         &
            &    vr(n_phi,n_theta)*this%Btc(n_phi,n_theta)- &
            &    vt(n_phi,n_theta)*this%Brc(n_phi,n_theta) )
         end do
      end do   ! theta loop
      !this%VxBr(:,:) = zero
      !this%VxBt(:,:) = zero
      !this%VxBp(:,:) = zero
      !$OMP END PARALLEL DO

      !------ Get the Lorentz force:
      !$OMP PARALLEL DO default(shared) &
      !$OMP& private(n_theta, n_phi, or1sn1, DyMagFac)
      do n_theta=1,n_theta_max
         !or1sn1=or1_3D(n_r)*osint1(n_theta)
         do n_phi=1,n_phi_max_3D
            !---- jxBr= r**2/(E*Pm) * ( curl(B)_t*B_p - curl(B)_p*B_t )
            this%jxBr(n_phi,n_theta)=1.0_cp*( &!DyMagFac*r2*(&
            &    this%curlBtc(n_phi,n_theta)*this%Bpc(n_phi,n_theta)- &
            &    this%curlBpc(n_phi,n_theta)*this%Btc(n_phi,n_theta) )
         end do
         !---- jxBt= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_p*B_r - curl(B)_r*B_p )
         do n_phi=1,n_phi_max_3D
            this%jxBt(n_phi,n_theta)=1.0_cp*( &!DyMagFac*or1sn1*(&
            &    this%curlBpc(n_phi,n_theta)*this%Brc(n_phi,n_theta)- &
            &    this%curlBrc(n_phi,n_theta)*this%Bpc(n_phi,n_theta) )
         end do
         !---- jxBp= 1/(E*Pm) * 1/(r*sin(theta)) * ( curl(B)_r*B_t - curl(B)_t*B_r )
         do n_phi=1,n_phi_max_3D
            this%jxBp(n_phi,n_theta)=1.0_cp*( &!DyMagFac*or1sn1*(&
            &    this%curlBrc(n_phi,n_theta)*this%Btc(n_phi,n_theta)- &
            &    this%curlBtc(n_phi,n_theta)*this%Brc(n_phi,n_theta) )
         end do
      end do   ! theta loop
      !$OMP END PARALLEL DO

      !-- Assemble parts of Lorentz-force on the spherical grid: 
      !-- curl( curl(B) x B )_z
      !-- (jxB) . e_s  and r . (jxB) . e_p
      do n_theta=1,n_theta_max
         do n_phi=1,n_phi_max_3D
            lfs(n_phi,n_theta)=sint(n_theta)*this%jxBr(n_phi,n_theta)+ &
            &      or1_3D(n_r)*cost(n_theta)*this%jxBt(n_phi,n_theta)
         end do
      end do
      !$OMP E

   end subroutine get_nl
!----------------------------------------------------------------------------
end module grid_space_arrays_mod
