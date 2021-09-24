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
   use namelists, only: l_heat_3D, l_mag_3D, l_mag_LF, r_icb, r_cmb, &
   &                    l_mag_alpha, alpha_fac, l_mag_inertia,       &
   &                    delta_fac, l_QG_basis
   use blocking, only: nRstart3D, nRstop3D
   use truncation_3D, only: n_theta_max, n_phi_max_3D
   use radial_functions, only: or1_3D, r_3D, rgrav_3D, beta!, tcond_3D
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
   subroutine get_nl(this, vr, vt, vp, n_r, buo, jxBs, jxBp)
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

      !-- Local variables:
      integer :: n_theta, n_phi
      real(cp) :: r2, or1sn1
      real(cp) :: bamp2, rfunc, Afunc

      r2 = r_3D(n_r)*r_3D(n_r)

      if ( l_heat_3D ) then
         !------ Get V T, the divergence of it is temperature advection:
         !$OMP PARALLEL DO default(shared) &
         !$OMP& private(n_theta, n_phi)
         do n_theta=1,n_theta_max
            do n_phi=1,n_phi_max_3D     ! calculate u*T components
               this%VTr(n_phi,n_theta) = &
               &    r2*vr(n_phi,n_theta)*this%Tc(n_phi,n_theta)
               this%VTt(n_phi,n_theta) = &
               &    or1_3D(n_r)*vt(n_phi,n_theta)*this%Tc(n_phi,n_theta)
               this%VTp(n_phi,n_theta) = &
               &    or1_3D(n_r)*vp(n_phi,n_theta)*this%Tc(n_phi,n_theta)
            end do
         end do  ! theta loop
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO default(shared) &
         !$OMP& private(n_theta, n_phi)
         !-- Assemble g/r * T on the spherical grid (1/r comes from the decomposition of Tg.e_r = Tg(s/r.e_s + z/r.e_z))
         !-- and the s factor cancels afterward because of the z-compotent of the curl in cylindrical coord.
         do n_theta=1,n_theta_max
            do n_phi=1,n_phi_max_3D
               buo(n_phi,n_theta)=this%Tc(n_phi,n_theta)*rgrav_3D(n_r)* &
               &                  or1_3D(n_r)
               if ( l_QG_basis ) then
                  !-- Additional buoyancy term from the QG basis projection: -\beta/s z^2 Tg.r (.e_z part of the curl)
                  !-- Can be simplified with z = cost.r , s = sint.r so that z^2/(sr) = cost^2/sint
                  buo(n_phi,n_theta)=buo(n_phi,n_theta) - this%Tc(n_phi,n_theta)* &
                  &                                      rgrav_3D(n_r)*beta(n_r)* &
                  &                  cost(n_theta)*cost(n_theta)*osint1(n_theta)
               end if
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
               this%VxBr(n_phi,n_theta)=r2*(                   &
               &    vt(n_phi,n_theta)*this%Bpc(n_phi,n_theta)- &
               &    vp(n_phi,n_theta)*this%Btc(n_phi,n_theta) )!&&!
               !&   (vt(n_phi,n_theta) - vzm(n_phi,n_theta)*sint(n_theta))*this%Bpc(n_phi,n_theta)- &
               !&   (vp(n_phi,n_theta) + vpm(n_phi,n_theta))*this%Btc(n_phi,n_theta) )!&&!

               this%VxBt(n_phi,n_theta)=or1sn1*(               &
               &    vp(n_phi,n_theta)*this%Brc(n_phi,n_theta)- &
               &    vr(n_phi,n_theta)*this%Bpc(n_phi,n_theta) )!&&!
               !&   (vp(n_phi,n_theta) + vpm(n_phi,n_theta))*this%Brc(n_phi,n_theta)- &
               !&   (vr(n_phi,n_theta) + vzm(n_phi,n_theta)*cost(n_theta))*this%Bpc(n_phi,n_theta) )!&&!

               this%VxBp(n_phi,n_theta)=or1sn1*(         &
               &    vr(n_phi,n_theta)*this%Btc(n_phi,n_theta)- &
               &    vt(n_phi,n_theta)*this%Brc(n_phi,n_theta) )!&&!
               !&   (vr(n_phi,n_theta) + vzm(n_phi,n_theta)*cost(n_theta))*this%Btc(n_phi,n_theta)- &
               !&   (vt(n_phi,n_theta) - vzm(n_phi,n_theta)*sint(n_theta))*this%Brc(n_phi,n_theta) )!&&!
               if ( l_mag_alpha ) then
                  !-- Following the form of (Chan etal., 2001; eq.25)
                  bamp2 = this%Brc(n_phi,n_theta)*this%Brc(n_phi,n_theta) + &
                  &       this%Btc(n_phi,n_theta)*this%Btc(n_phi,n_theta) + &
                  &       this%Bpc(n_phi,n_theta)*this%Bpc(n_phi,n_theta)
                  rfunc = sin(pi*(r_3D(n_r)-r_icb))
                  !rfunc = 4.*r_3D(n_r)*sint(n_theta) *                   &
                  !&      ( r_cmb - r_3D(n_r)*sint(n_theta) ) / r_cmb**2.
                  Afunc = alpha_fac*cost(n_theta)*rfunc

                  this%VxBr(n_phi,n_theta)=this%VxBr(n_phi,n_theta) + r2*(             &
                  &                        Afunc*this%Brc(n_phi,n_theta)/(1. + bamp2) )

                  this%VxBt(n_phi,n_theta)=this%VxBt(n_phi,n_theta) + or1sn1*(         &
                  &                        Afunc*this%Btc(n_phi,n_theta)/(1. + bamp2) )

                  this%VxBp(n_phi,n_theta)=this%VxBp(n_phi,n_theta) + or1sn1*(         &
                  &                        Afunc*this%Bpc(n_phi,n_theta)/(1. + bamp2) )
               end if

               if ( l_mag_inertia ) then
                  !-- Following (Davidson and Ranjan, 2015; eq.53)
                  !--    <u x b>_s = \delta/Pm <u_p^2>.B_s  --> on VxBr and VxBt
                  !-- &  <u x b>_p = \delta/Pm <u_s^2>.B_p  --> on VxBp
                  !--              = terms coming from inertial waves
                  this%VxBr(n_phi,n_theta)=this%VxBr(n_phi,n_theta) + delta_fac*   &
                  &                        (sint(n_theta)*vp(n_phi,n_theta))**2.*( &
                  &                 or1sn1*this%Brc(n_phi,n_theta)*sint(n_theta) + &
                  &                     r2*this%Btc(n_phi,n_theta)*cost(n_theta) )

                  this%VxBt(n_phi,n_theta)=this%VxBt(n_phi,n_theta) + delta_fac*   &
                  &                        (cost(n_theta)*vp(n_phi,n_theta))**2.*( &
                  &                 or1sn1*this%Brc(n_phi,n_theta)*sint(n_theta) + &
                  &                     r2*this%Btc(n_phi,n_theta)*cost(n_theta) )

                  this%VxBp(n_phi,n_theta)=this%VxBp(n_phi,n_theta) + delta_fac*(  &
                  &                ( (or1sn1*vr(n_phi,n_theta)*sint(n_theta) +     &
                  &                    r2*vt(n_phi,n_theta)*cost(n_theta))**2. )*  &
                  &                                      this%Bpc(n_phi,n_theta) )
               end if
            end do
         end do   ! theta loop
         !$OMP END PARALLEL DO

         if ( l_mag_LF ) then
            !------ Get the Lorentz force:
            !------ We will only need to compute the z component of Vx(jxB)
            !$OMP PARALLEL DO default(shared) &
            !$OMP& private(n_theta, n_phi, or1sn1)
            do n_theta=1,n_theta_max
               or1sn1=or1_3D(n_r)*osint1(n_theta)
               do n_phi=1,n_phi_max_3D
                  !---- jxBs= 1/(E*Pm) * ( curl(B)_p*B_z - curl(B)_z*B_p )
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
                  if ( l_QG_basis ) then
                     !-- Additional buoyancy term from the QG basis projection: +\beta/s z Vx(jxB)_s
                     !-- Only need to compute the s component of Vx(jxB)
                     !-- Vx(jxB)_s = 1/s \partial_p jxB_z - \partial_z jxB_p (<- computed above)
                     !-- jxBz = 1/(E*Pm) * ( curl(B)_p*B_s - curl(B)_s*B_p )
                     !--      = cost * jxBr - sint * jxBt
                     !-- Can be added to jxBs because it will get 1/s \partial_\phi later
                     jxBs(n_phi,n_theta)=jxBs(n_phi,n_theta) -(cost(n_theta)*r2*( &
                     &   this%curlBtc(n_phi,n_theta)*this%Bpc(n_phi,n_theta)-     &
                     &   this%curlBpc(n_phi,n_theta)*this%Btc(n_phi,n_theta) )    &
                     &                                -sint(n_theta)*or1sn1*(     &
                     &   this%curlBpc(n_phi,n_theta)*this%Brc(n_phi,n_theta)-     &
                     &   this%curlBrc(n_phi,n_theta)*this%Bpc(n_phi,n_theta) ) )* &
                     &                          beta(n_r)*r_3D(n_r)*cost(n_theta)
                  end if
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
