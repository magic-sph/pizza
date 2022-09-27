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
   &                    delta_fac, l_QG_basis, l_lin_solve, l_mag_B0, beta_shift
   use blocking, only: nRstart3D, nRstop3D
   use truncation_3D, only: n_theta_max, n_phi_max_3D, n_r_max_3D
   use radial_functions, only: or1_3D, r_3D, rgrav_3D, beta!, tcond_3D
   use fields, only: B0r_3D_Rloc, B0t_3D_Rloc, B0p_3D_Rloc, curlB0r_3D_Rloc, &
       &             curlB0t_3D_Rloc, curlB0p_3D_Rloc
   use horizontal, only: cost, sint, osint1

   implicit none

   private

   type, public, extends(general_arrays_t) :: grid_space_arrays_t
      !----- Nonlinear terms in phi/theta space: 
      real(cp), allocatable :: VTr(:,:), VTt(:,:), VTp(:,:)
      real(cp), allocatable :: VxBr(:,:), VxBt(:,:), VxBp(:,:)
      real(cp), pointer :: Alphac(:,:)

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
         if ( l_mag_alpha ) then
            allocate( this%Alphac(n_phi_max_3D,n_theta_max) )
            bytes_allocated=bytes_allocated + n_phi_max_3D*n_theta_max*SIZEOF_DEF_REAL
         else
            allocate( this%Alphac(1,1) )
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
            this%Alphac(:,:)=0.0_cp !if ( l_mag_alpha )
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
            deallocate( this%Alphac )
         end if
      end if

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine get_nl(this, vr, vt, vp, n_r, buo, jxBs, jxBp, jxBz)
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
      real(cp), intent(out) :: jxBz(n_phi_max_3D,n_theta_max)

      !-- Local variables:
      integer :: n_theta, n_phi, n_count
      real(cp) :: r2, or1sn1
      real(cp) :: rfunc, asign!, bamp2!, Afunc
      real(cp) :: Bsavg(n_theta_max), Bpavg(n_theta_max)
      real(cp) :: vpavg(n_theta_max)
      real(cp) :: vsfluct2(n_theta_max), vpfluct2(n_theta_max)

      r2 = r_3D(n_r)*r_3D(n_r)

      if ( l_heat_3D ) then
         if ( .not. l_lin_solve ) then  !-- Non-Linear terms?
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
         end if

         !if ( .not. l_mag_B0 ) then
         !-- Assemble g/r * T on the spherical grid (1/r comes from the decomposition of Tg.e_r = Tg(s/r.e_s + z/r.e_z))
         !-- and the s factor cancels afterward because of the z-compotent of the curl in cylindrical coord.
         !$OMP PARALLEL DO default(shared) &
         !$OMP& private(n_theta, n_phi)
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
         !end if ! l_mag_B0?
      end if ! l_heat_3D?
      !this%VTr(:,:) = 0.0_cp
      !this%VTt(:,:) = 0.0_cp
      !this%VTp(:,:) = 0.0_cp
      !buo(:,:) = 0.0_cp

      if ( l_mag_inertia .and. (l_mag_3D .and. (.not. l_lin_solve)) ) then
         !------ Get V and B phi averaged:
         !$OMP PARALLEL DO default(shared) &
         !$OMP& private(n_theta)
         Bsavg(:) = 0.0_cp
         Bpavg(:) = 0.0_cp
         vpavg(:) = 0.0_cp
         do n_theta=1,n_theta_max
            do n_phi=1,n_phi_max_3D
               Bsavg(n_theta) = Bsavg(n_theta) + ( this%Brc(n_phi,n_theta)*sint(n_theta) +    &
               &                       this%Btc(n_phi,n_theta)*cost(n_theta) )
               Bpavg(n_theta) = Bpavg(n_theta) + this%Bpc(n_phi,n_theta)
               vpavg(n_theta) = vpavg(n_theta) + vp(n_phi,n_theta)
            end do
            Bsavg(n_theta) = Bsavg(n_theta)/n_phi_max_3D
            Bpavg(n_theta) = Bpavg(n_theta)/n_phi_max_3D
            vpavg(n_theta) = vpavg(n_theta)/n_phi_max_3D
         end do
         vsfluct2(:) = 0.0_cp
         vpfluct2(:) = 0.0_cp
         do n_theta=1,n_theta_max
            do n_phi=1,n_phi_max_3D
               vsfluct2(n_theta) = vsfluct2(n_theta) + ( vr(n_phi,n_theta)*sint(n_theta) +   &
               &                                         vt(n_phi,n_theta)*cost(n_theta) )**2
               vpfluct2(n_theta) = vpfluct2(n_theta) + (vp(n_phi,n_theta) - vpavg(n_theta))**2
            end do
            vsfluct2(n_theta) = vsfluct2(n_theta)/n_phi_max_3D
            vpfluct2(n_theta) = vpfluct2(n_theta)/n_phi_max_3D
         end do
         !$OMP END PARALLEL DO
      end if

      if ( l_mag_3D .and. (.not. l_lin_solve) ) then
         !------ Get (V x B) , the curl of this is the dynamo term:
         !$OMP PARALLEL DO default(shared) &
         !$OMP& private(n_theta, n_phi, or1sn1)
         do n_theta=1,n_theta_max
            or1sn1=or1_3D(n_r)*osint1(n_theta)
            do n_phi=1,n_phi_max_3D     ! calculate V x B components
               if ( .not. l_mag_B0 ) then
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
               else !-- Background field?
                  this%VxBr(n_phi,n_theta)=r2*(                   &
                  &    vt(n_phi,n_theta)*B0p_3D_Rloc(n_phi,n_theta,n_r)- &
                  &    vp(n_phi,n_theta)*B0t_3D_Rloc(n_phi,n_theta,n_r) )!

                  this%VxBt(n_phi,n_theta)=or1sn1*(               &
                  &    vp(n_phi,n_theta)*B0r_3D_Rloc(n_phi,n_theta,n_r)- &
                  &    vr(n_phi,n_theta)*B0p_3D_Rloc(n_phi,n_theta,n_r) )

                  this%VxBp(n_phi,n_theta)=or1sn1*(         &
                  &    vr(n_phi,n_theta)*B0t_3D_Rloc(n_phi,n_theta,n_r)- &
                  &    vt(n_phi,n_theta)*B0r_3D_Rloc(n_phi,n_theta,n_r) )
               end if

               if ( l_mag_alpha ) then !.and. r_3D(n_r)*sint(n_theta) >= r_icb ) then
                  !-- Following the form of (Chan etal., 2001; eq.25) --> only when nothing can prevent B from growing
                  !bamp2 = this%Brc(n_phi,n_theta)*this%Brc(n_phi,n_theta) + &
                  !&       this%Btc(n_phi,n_theta)*this%Btc(n_phi,n_theta) + &
                  !&       this%Bpc(n_phi,n_theta)*this%Bpc(n_phi,n_theta)
                  !rfunc = sin(pi*(r_3D(n_r)-r_icb))!*sin(pi*(r_3D(n_r)*sint(n_theta)-r_icb))
                  rfunc = 4.*r_3D(n_r)*sint(n_theta) *                   &
                  &      ( r_cmb - r_3D(n_r)*sint(n_theta) ) / r_cmb**2.
                  this%Alphac(n_phi,n_theta) = alpha_fac*rfunc*cost(n_theta)

                  this%VxBr(n_phi,n_theta)=this%VxBr(n_phi,n_theta) + r2*(              &
                  &    this%Alphac(n_phi,n_theta)*this%Brc(n_phi,n_theta))!/(1. + bamp2)

                  this%VxBt(n_phi,n_theta)=this%VxBt(n_phi,n_theta) + or1sn1*(          &
                  &    this%Alphac(n_phi,n_theta)*this%Btc(n_phi,n_theta))!/(1. + bamp2)

                  this%VxBp(n_phi,n_theta)=this%VxBp(n_phi,n_theta) + or1sn1*(          &
                  &    this%Alphac(n_phi,n_theta)*this%Bpc(n_phi,n_theta))!/(1. + bamp2)
               end if

               if ( l_mag_inertia ) then
                  !-- Following (Davidson and Ranjan, 2015; eq.53)
                  !--    <u x b>_s = \delta/Pm <u_p^2>.B_s  --> on VxBr and VxBt
                  !-- &  <u x b>_p = \delta/Pm <u_s^2>.B_p  --> on VxBp
                  !--              = terms coming from inertial waves
                  !--> Warning:: 1/r2 or 1/rsint is applied to the full product VxB;
                  !--            the different quantities have already been rescaled
                  asign = 1.0_cp!abs(cost(n_theta))/cost(n_theta) !-- +/-1 sign
                  if ( n_r == 1 .or. n_r == n_r_max_3D ) then
                     delta_fac = 1.0_cp
                  else
                     delta_fac = 1.0_cp/sqrt(vpfluct2(n_theta) + vsfluct2(n_theta)+epsilon(1.0_cp))
                  end if
                  this%VxBr(n_phi,n_theta)=asign*delta_fac*r2*(&!this%VxBr(n_phi,n_theta)+
                  &                        sint(n_theta)*vpfluct2(n_theta)*Bsavg(n_theta))!(  &
                  !&                        sint(n_theta)*sqrt((vp(n_phi,n_theta)-vpavg(n_theta))**2.)*Bsavg(n_theta))!(  &
                  !&                        this%Brc(n_phi,n_theta)*sint(n_theta) +  &
                  !&                        this%Btc(n_phi,n_theta)*cost(n_theta) ) )

                  this%VxBt(n_phi,n_theta)=asign*delta_fac*or1sn1*(&!this%VxBt(n_phi,n_theta)+
                  &                        cost(n_theta)*vpfluct2(n_theta)*Bsavg(n_theta))!(    &
                  !&                        cost(n_theta)*sqrt((vp(n_phi,n_theta)-vpavg(n_theta))**2.)*Bsavg(n_theta))!(    &
                  !&                        this%Brc(n_phi,n_theta)*sint(n_theta) +    &
                  !&                        this%Btc(n_phi,n_theta)*cost(n_theta) )  )

                  this%VxBp(n_phi,n_theta)=asign*delta_fac*or1sn1*(&!this%VxBp(n_phi,n_theta)+
                  &                       ( vsfluct2(n_theta) )*Bpavg(n_theta))!        &
                  !&                       ( sqrt((vr(n_phi,n_theta)*sint(n_theta) +        &
                  !&                       vt(n_phi,n_theta)*cost(n_theta))**2.) )*     &
                  !&                                      Bpavg(n_theta))!this%Bpc(n_phi,n_theta)    )
               end if
            end do
         end do   ! theta loop
         !$OMP END PARALLEL DO

         if ( l_mag_LF ) then
            !------ Get the Lorentz force:
            !------ We will only need to compute the z component of Vx(jxB)
            !$OMP PARALLEL DO default(shared) &
            !$OMP& private(n_theta, n_phi)
            n_count=0
            do n_theta=1,n_theta_max
               !--> Warning:: 1/rsint or 1/r2 is only there if back-transformed are done later
               !or1sn1=or1_3D(n_r)*osint1(n_theta)
               do n_phi=1,n_phi_max_3D
                  if ( .not. l_mag_B0 ) then
                     !---- jxBs= 1/(E*Pm) * ( curl(B)_p*B_z - curl(B)_z*B_p )
                     !--       = sint * jxBr + cost * jxBt
                     jxBs(n_phi,n_theta)=sint(n_theta)*(&!r2*(                    &!
                     &   this%curlBtc(n_phi,n_theta)*this%Bpc(n_phi,n_theta)-  &
                     &   this%curlBpc(n_phi,n_theta)*this%Btc(n_phi,n_theta) ) &
                     &                  +cost(n_theta)*(&!or1sn1*(                &!
                     &   this%curlBpc(n_phi,n_theta)*this%Brc(n_phi,n_theta)-  &
                     &   this%curlBrc(n_phi,n_theta)*this%Bpc(n_phi,n_theta) )

                     !---- jxBp= 1/(E*Pm) * ( curl(B)_r*B_t - curl(B)_t*B_r )
                     jxBp(n_phi,n_theta) =(&!or1sn1*(                             &!
                     &    this%curlBrc(n_phi,n_theta)*this%Btc(n_phi,n_theta)- &
                     &    this%curlBtc(n_phi,n_theta)*this%Brc(n_phi,n_theta) )
                  else !-- Background field?
                     !---- j0xBs= sint * j0xBr + cost * j0xBt
                     jxBs(n_phi,n_theta)=sint(n_theta)*(&!r2*(                     &!
                     &   curlB0t_3D_Rloc(n_phi,n_theta,n_r)*this%Bpc(n_phi,n_theta)-  &
                     &   curlB0p_3D_Rloc(n_phi,n_theta,n_r)*this%Btc(n_phi,n_theta) ) &
                     &                  +cost(n_theta)*(&!or1sn1*(                 &!
                     &   curlB0p_3D_Rloc(n_phi,n_theta,n_r)*this%Brc(n_phi,n_theta)-  &
                     &   curlB0r_3D_Rloc(n_phi,n_theta,n_r)*this%Bpc(n_phi,n_theta) )
                     !---- + jxB0s= sint * jxB0r + cost * jxB0t
                     jxBs(n_phi,n_theta)=jxBs(n_phi,n_theta)+sint(n_theta)*(&!r2*( &!
                     &   this%curlBtc(n_phi,n_theta)*B0p_3D_Rloc(n_phi,n_theta,n_r)-  &
                     &   this%curlBpc(n_phi,n_theta)*B0t_3D_Rloc(n_phi,n_theta,n_r) ) &
                     &                  +cost(n_theta)*(&!or1sn1*(                 &!
                     &   this%curlBpc(n_phi,n_theta)*B0r_3D_Rloc(n_phi,n_theta,n_r)-  &
                     &   this%curlBrc(n_phi,n_theta)*B0p_3D_Rloc(n_phi,n_theta,n_r) )
                     !if ( l_mag_B0LFsmall ) then
                        !---- + jxBs= sint * jxBr + cost * jxBt
                        jxBs(n_phi,n_theta)=jxBs(n_phi,n_theta)+sint(n_theta)*(&!r2*( &!
                        &   this%curlBtc(n_phi,n_theta)*this%Bpc(n_phi,n_theta)-  &
                        &   this%curlBpc(n_phi,n_theta)*this%Btc(n_phi,n_theta) ) &
                        &                  +cost(n_theta)*(&!or1sn1*(                &!
                        &   this%curlBpc(n_phi,n_theta)*this%Brc(n_phi,n_theta)-  &
                        &   this%curlBrc(n_phi,n_theta)*this%Bpc(n_phi,n_theta) )
                     !end if

                     !---- j0xBp= j0_r*B_t - j0_t*B_r
                     jxBp(n_phi,n_theta) =(&!or1sn1*(                              &!
                     &    curlB0r_3D_Rloc(n_phi,n_theta,n_r)*this%Btc(n_phi,n_theta)- &
                     &    curlB0t_3D_Rloc(n_phi,n_theta,n_r)*this%Brc(n_phi,n_theta) )
                     !---- + jxB0p= j_r*B0_t - j_t*B0_r
                     jxBp(n_phi,n_theta) =jxBp(n_phi,n_theta)+(&!or1sn1*(          &!
                     &    this%curlBrc(n_phi,n_theta)*B0t_3D_Rloc(n_phi,n_theta,n_r)- &
                     &    this%curlBtc(n_phi,n_theta)*B0r_3D_Rloc(n_phi,n_theta,n_r) )
                     !if ( l_mag_B0LFsmall ) then
                        !---- + jxBp= j_r*B_t - j_t*B_r
                        jxBp(n_phi,n_theta) =jxBp(n_phi,n_theta)+(&!or1sn1*(       &!
                        &    this%curlBrc(n_phi,n_theta)*this%Btc(n_phi,n_theta)- &
                        &    this%curlBtc(n_phi,n_theta)*this%Brc(n_phi,n_theta) )
                     !end if
                  end if

                  if ( l_QG_basis ) then
                     !-- Additional Lorentz force term from the QG basis projection: +\beta/s z Vx(jxB)_s
                     !-- Only need to compute the s component of Vx(jxB)
                     !-- Vx(jxB)_s = 1/s \partial_p jxB_z - \partial_z jxB_p (<- computed above)
                     if ( .not. l_mag_B0 ) then
                        !-- jxBz = 1/(E*Pm) * ( curl(B)_p*B_s - curl(B)_s*B_p )
                        !--      = cost * jxBr - sint * jxBt
                        !-- Can be added to jxBs because it will get -1/s \partial_\phi later --> but problematic
                        !jxBs(n_phi,n_theta)=jxBs(n_phi,n_theta) -(cost(n_theta)*(&!r2*( & !-- W:: Seems to work with B-pot but troubles with jxB_ana!!!
                        jxBz(n_phi,n_theta)=(cost(n_theta)*(&!
                        &   this%curlBtc(n_phi,n_theta)*this%Bpc(n_phi,n_theta)-     &
                        &   this%curlBpc(n_phi,n_theta)*this%Btc(n_phi,n_theta) )    &
                        &                                -sint(n_theta)*(&!or1sn1*(     &
                        &   this%curlBpc(n_phi,n_theta)*this%Brc(n_phi,n_theta)-     &
                        &   this%curlBrc(n_phi,n_theta)*this%Bpc(n_phi,n_theta) ) )* &
                        &                          r_3D(n_r)*cost(n_theta)!*beta(n_r)
                        !&   (-r_3D(n_r)*sint(n_theta)/(r_cmb**2 - (r_3D(n_r)*cost(n_theta))**2))*r_3D(n_r)*cost(n_theta)
                        !&   (-r_3D(n_r)*sint(n_theta)/((r_cmb+beta_shift)**2 - (r_3D(n_r)*cost(n_theta))**2))*r_3D(n_r)*cost(n_theta)
                     else
                        !-- j0xBz = cost * j0xBr - sint * j0xBt
                        jxBz(n_phi,n_theta)=(cost(n_theta)*(&!
                        &   curlB0t_3D_Rloc(n_phi,n_theta,n_r)*this%Bpc(n_phi,n_theta)-     &
                        &   curlB0p_3D_Rloc(n_phi,n_theta,n_r)*this%Btc(n_phi,n_theta) )    &
                        &                                -sint(n_theta)*(&!or1sn1*(     &
                        &   curlB0p_3D_Rloc(n_phi,n_theta,n_r)*this%Brc(n_phi,n_theta)-     &
                        &   curlB0r_3D_Rloc(n_phi,n_theta,n_r)*this%Bpc(n_phi,n_theta) ) )* &
                        &                          r_3D(n_r)*cost(n_theta)!*beta(n_r)
                        !-- + jxB0z = cost * jxB0r - sint * jxB0t
                        jxBz(n_phi,n_theta)=jxBz(n_phi,n_theta)+(cost(n_theta)*(&!
                        &   this%curlBtc(n_phi,n_theta)*B0p_3D_Rloc(n_phi,n_theta,n_r)-     &
                        &   this%curlBpc(n_phi,n_theta)*B0t_3D_Rloc(n_phi,n_theta,n_r) )    &
                        &                                -sint(n_theta)*(&!or1sn1*(     &
                        &   this%curlBpc(n_phi,n_theta)*B0r_3D_Rloc(n_phi,n_theta,n_r)-     &
                        &   this%curlBrc(n_phi,n_theta)*B0p_3D_Rloc(n_phi,n_theta,n_r) ) )* &
                        &                          r_3D(n_r)*cost(n_theta)!*beta(n_r)
                        !if ( l_mag_B0LFsmall ) then
                           !-- + jxBz = cost * jxBr - sint * jxBt
                           jxBz(n_phi,n_theta)=jxBz(n_phi,n_theta)+(cost(n_theta)*(&!
                           &   this%curlBtc(n_phi,n_theta)*this%Bpc(n_phi,n_theta)-     &
                           &   this%curlBpc(n_phi,n_theta)*this%Btc(n_phi,n_theta) )    &
                           &                                -sint(n_theta)*(&!or1sn1*(     &
                           &   this%curlBpc(n_phi,n_theta)*this%Brc(n_phi,n_theta)-     &
                           &   this%curlBrc(n_phi,n_theta)*this%Bpc(n_phi,n_theta) ) )* &
                           &                          r_3D(n_r)*cost(n_theta)!*beta(n_r)
                        !end if
                     end if
                  end if
               end do
            end do   ! theta loop
            !$OMP END PARALLEL DO
         end if ! l_mag_LF?
      end if ! l_mag_3D?
      !this%VxBr(:,:) = 0.0_cp
      !this%VxBt(:,:) = 0.0_cp
      !this%VxBp(:,:) = 0.0_cp

#ifdef aDEBUG
         block
            integer :: filehandle
            real(cp) :: VxBs(n_phi_max_3D,n_theta_max)
            real(cp) :: VxBp(n_phi_max_3D,n_theta_max)
            if ( n_r == 12 ) then
               print*, r_3D(n_r)
               do n_theta=1,n_theta_max
                  do n_phi=1,n_phi_max_3D     ! calculate V x B components
                        asign = abs(cost(n_theta))/cost(n_theta) !-- +/-1 sign
                        !VxBs(n_phi,n_theta)=asign*delta_fac*sqrt(vpfluct(n_theta))*Bsavg(n_theta)
                        jxBs(n_phi,n_theta)=asign*delta_fac*sqrt(vpfluct(n_theta))*Bsavg(n_theta)
                        !VxBp(n_phi,n_theta)=asign*delta_fac*sqrt(vsfluct(n_theta))*Bpavg(n_theta)
                        jxBp(n_phi,n_theta)=asign*delta_fac*sqrt(vsfluct(n_theta))*Bpavg(n_theta)
                  end do
               end do   ! theta loop

               !open(newunit=filehandle, file='test_VxB', form='unformatted', access='stream')
               !   write(filehandle) VxBs, VxBp
               !close(filehandle)
            end if

         end block
#endif

   end subroutine get_nl
!----------------------------------------------------------------------------
end module grid_space_arrays_mod
