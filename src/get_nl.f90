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
   use constants, only: pi
   use mem_alloc, only: bytes_allocated
   use namelists, only: l_heat_3D, l_mag_3D, l_mag_LF, r_icb, r_cmb
   use blocking, only: nRstart3D, nRstop3D
   use truncation_3D, only: n_theta_max, n_phi_max_3D
   use radial_functions, only: or1_3D, r_3D, rgrav_3D!, tcond_3D
   use horizontal, only: cost, sint, osint1

   implicit none

   private

   type, public, extends(general_arrays_t) :: grid_space_arrays_t
      !----- Nonlinear terms in phi/theta space: 
      real(cp), allocatable :: VTr(:,:), VTt(:,:), VTp(:,:)
      real(cp), allocatable :: VxBr(:,:), VxBt(:,:), VxBp(:,:)

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

      !----- Fields calculated from these help arrays by legtf:
      if ( l_heat_3D ) then
         allocate( this%VTr(n_phi_max_3D,n_theta_max) )
         allocate( this%VTt(n_phi_max_3D,n_theta_max) )
         allocate( this%VTp(n_phi_max_3D,n_theta_max) )
         allocate( this%Tc(n_phi_max_3D,n_theta_max) )
         bytes_allocated=bytes_allocated + 4*n_phi_max_3D*n_theta_max*SIZEOF_DEF_REAL
      end if
      if ( l_mag_3D ) then
         allocate( this%VxBr(n_phi_max_3D,n_theta_max) )
         allocate( this%VxBt(n_phi_max_3D,n_theta_max) )
         allocate( this%VxBp(n_phi_max_3D,n_theta_max) )
         allocate( this%Brc(n_phi_max_3D,n_theta_max) )
         allocate( this%Btc(n_phi_max_3D,n_theta_max) )
         allocate( this%Bpc(n_phi_max_3D,n_theta_max) )
         bytes_allocated=bytes_allocated + 6*n_phi_max_3D*n_theta_max*SIZEOF_DEF_REAL
         if ( l_mag_LF ) then
            allocate( this%curlBrc(n_phi_max_3D,n_theta_max) )
            allocate( this%curlBtc(n_phi_max_3D,n_theta_max) )
            allocate( this%curlBpc(n_phi_max_3D,n_theta_max) )
            bytes_allocated=bytes_allocated + 3*n_phi_max_3D*n_theta_max*SIZEOF_DEF_REAL
         end if
      end if

      if ( l_heat_3D ) then
         this%VTr(:,:)=0.0_cp
         this%VTt(:,:)=0.0_cp
         this%VTp(:,:)=0.0_cp
         this%Tc(:,:)=0.0_cp
      end if
      if ( l_mag_3D ) then
         this%VxBr(:,:)=0.0_cp
         this%VxBt(:,:)=0.0_cp
         this%VxBp(:,:)=0.0_cp
         this%Brc(:,:)=0.0_cp
         this%Btc(:,:)=0.0_cp
         this%Bpc(:,:)=0.0_cp
         if ( l_mag_LF ) then
            this%curlBrc(:,:)=0.0_cp
            this%curlBtc(:,:)=0.0_cp
            this%curlBpc(:,:)=0.0_cp
         end if
      end if

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(grid_space_arrays_t) :: this

      if ( l_heat_3D ) then
         deallocate( this%VTr, this%VTt, this%VTp )
          deallocate( this%Tc )
      end if
      if ( l_mag_3D ) then
         deallocate( this%VxBr, this%VxBt, this%VxBp )
         deallocate( this%Brc, this%Btc, this%Bpc )
         if ( l_mag_LF ) then
            deallocate( this%curlBrc, this%curlBtc, this%curlBpc )
         end if
      end if

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine get_nl(this, vr, vt, vp, n_r, buo, jxBs, jxBp, vpm, vzm)
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
      real(cp), intent(out) :: jxBs(n_phi_max_3D,n_theta_max)
      real(cp), intent(out) :: jxBp(n_phi_max_3D,n_theta_max)

      real(cp), intent(in) :: vpm(n_phi_max_3D,n_theta_max)
      real(cp), intent(in) :: vzm(n_phi_max_3D,n_theta_max)

      !-- Local variables:
      integer :: n_theta, n_phi
      real(cp) :: or1sn1, r2!, bamp2, rfunc

      r2 = r_3D(n_r)*r_3D(n_r)

      if ( l_heat_3D ) then
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
               buo(n_phi,n_theta)=this%Tc(n_phi,n_theta)*rgrav_3D(n_r)!* &
               !&                  or1_3D(n_r)
            end do
         end do
         !$OMP END PARALLEL DO
      end if ! l_heat_3D?

      if ( l_mag_3D ) then
         !------ Get (V x B) , the curl of this is the dynamo term:
         !$OMP PARALLEL DO default(shared) &
         !$OMP& private(n_theta, n_phi, or1sn1)
         do n_theta=1,n_theta_max
            or1sn1=or1_3D(n_r)*osint1(n_theta)
            do n_phi=1,n_phi_max_3D     ! calculate V x B components
               !bamp2 = this%Brc(n_phi,n_theta)*this%Brc(n_phi,n_theta) + &
               !&       this%Btc(n_phi,n_theta)*this%Btc(n_phi,n_theta) + &
               !&       this%Bpc(n_phi,n_theta)*this%Bpc(n_phi,n_theta)
               !rfunc = sin(pi*(r_3D(n_r)-r_icb))
               !rfunc = 4.*r_3D(n_r)*sint(n_theta) *                   &
               !&      ( r_cmb - r_3D(n_r)*sint(n_theta) ) / r_cmb**2.
               !rfunc = (r_3D(n_r)*sint(n_theta)* r_3D(n_r)*cost(n_theta)*   &
               !&( ( 2.*sqrt(r_cmb**2. - (r_3D(n_r)*sint(n_theta))**2.)**2. -&
               !&  (r_3D(n_r)*cost(n_theta))**2.) / (2.*sqrt(r_cmb**2. -     &
               !&  (r_3D(n_r)*sint(n_theta))**2.)**3.) ) ) / r_cmb
               this%VxBr(n_phi,n_theta)=r2*(                   &
               !&    vt(n_phi,n_theta)*this%Bpc(n_phi,n_theta)- &
               !&    vp(n_phi,n_theta)*this%Btc(n_phi,n_theta) &!)!&
               &   (vt(n_phi,n_theta) - vzm(n_phi,n_theta)*sint(n_theta))*this%Bpc(n_phi,n_theta)- &
               &   (vp(n_phi,n_theta) + vpm(n_phi,n_theta))*this%Btc(n_phi,n_theta) )!&&!
               !&    + 50.0*cost(n_theta)*rfunc*this%Brc(n_phi,n_theta)/(1. + bamp2))!sin(pi*(r_3D(n_r)-r_icb))
               !&    + vzm(n_phi,n_theta)*cost(n_theta)*this%Brc(n_phi,n_theta) )

               this%VxBt(n_phi,n_theta)=or1sn1*(               &
               !&    vp(n_phi,n_theta)*this%Brc(n_phi,n_theta)- &
               !&    vr(n_phi,n_theta)*this%Bpc(n_phi,n_theta) &!)!&
               &   (vp(n_phi,n_theta) + vpm(n_phi,n_theta))*this%Brc(n_phi,n_theta)- &
               &   (vr(n_phi,n_theta) + vzm(n_phi,n_theta)*cost(n_theta))*this%Bpc(n_phi,n_theta) )!&&!
               !&    + 50.0*cost(n_theta)*rfunc*this%Btc(n_phi,n_theta)/(1. + bamp2))!sin(pi*(r_3D(n_r)-r_icb))
               !&    - vzm(n_phi,n_theta)*sint(n_theta)*this%Btc(n_phi,n_theta) )

               this%VxBp(n_phi,n_theta)=or1sn1*(         &
               !&    vr(n_phi,n_theta)*this%Btc(n_phi,n_theta)- &
               !&    vt(n_phi,n_theta)*this%Brc(n_phi,n_theta) &!)!&
               &   (vr(n_phi,n_theta) + vzm(n_phi,n_theta)*cost(n_theta))*this%Btc(n_phi,n_theta)- &
               &   (vt(n_phi,n_theta) - vzm(n_phi,n_theta)*sint(n_theta))*this%Brc(n_phi,n_theta) )!&&!
               !&    + 50.0*cost(n_theta)*rfunc*this%Bpc(n_phi,n_theta)/(1. + bamp2))!sin(pi*(r_3D(n_r)-r_icb))
               !&    + vpm(n_phi,n_theta)*this%Bpc(n_phi,n_theta) )
            end do
         end do   ! theta loop
         !$OMP END PARALLEL DO

         if ( l_mag_LF ) then
            !------ Get the Lorentz force:
            !$OMP PARALLEL DO default(shared) &
            !$OMP& private(n_theta, n_phi, or1sn1)
            do n_theta=1,n_theta_max
               or1sn1=or1_3D(n_r)*osint1(n_theta)
               do n_phi=1,n_phi_max_3D
                  !---- jxBs= 1/(E*Pm) * ( curl(B)_z*B_p - curl(B)_p*B_z )
                  !--       = sint * jxBr + cost * jxBt
                  jxBs(n_phi,n_theta)=sint(n_theta)*r2*(                    &!(&!
                  &   this%curlBtc(n_phi,n_theta)*this%Bpc(n_phi,n_theta)-  &
                  &   this%curlBpc(n_phi,n_theta)*this%Btc(n_phi,n_theta) ) &
                  &                  +cost(n_theta)*or1sn1*(                &!(&!
                  &   this%curlBpc(n_phi,n_theta)*this%Brc(n_phi,n_theta)-  &
                  &   this%curlBrc(n_phi,n_theta)*this%Bpc(n_phi,n_theta) )

                  !---- jxBp= 1/(E*Pm) * ( curl(B)_r*B_t - curl(B)_t*B_r )
                  jxBp(n_phi,n_theta) =or1sn1*(                             &!(&!
                  &    this%curlBrc(n_phi,n_theta)*this%Btc(n_phi,n_theta)- &
                  &    this%curlBtc(n_phi,n_theta)*this%Brc(n_phi,n_theta) )
               end do
            end do   ! theta loop
            !$OMP END PARALLEL DO
         end if ! l_mag_LF?
      end if ! l_mag_3D?
      !this%VxBr(:,:) = 0.0_cp
      !this%VxBt(:,:) = 0.0_cp
      !this%VxBp(:,:) = 0.0_cp

   end subroutine get_nl
!----------------------------------------------------------------------------
end module grid_space_arrays_mod
