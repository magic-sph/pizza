module nonlinear_lm_mod

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use blocking_lm, only: lm2l, lm2m, lm2lmP, lmP2lmPS, lmP2lmPA
   use truncation_3D, only: lm_max, lmP_max
   use radial_functions, only: r_3D, or2_3D
   use horizontal, only: dTheta1S, dTheta1A, dPhi, dTheta2S, dTheta2A
   use namelists, only: l_mag_3D, l_mag_LF
   use constants, only: zero, two

   implicit none
 
   type :: nonlinear_lm_t
      !----- Nonlinear terms in lm-space: 
      complex(cp), allocatable :: VTrLM(:),  VTtLM(:),  VTpLM(:)
      complex(cp), allocatable :: VxBrLM(:), VxBtLM(:), VxBpLM(:)
      complex(cp), allocatable :: jxBrLM(:), jxBtLM(:), jxBpLM(:)

 
   contains
 
      procedure :: initialize
      procedure :: finalize
      procedure :: get_td
 
   end type nonlinear_lm_t

contains

   subroutine initialize(this,lmP_max)

      class(nonlinear_lm_t) :: this
      integer, intent(in) :: lmP_max

      allocate( this%VTrLM(lmP_max), this%VTtLM(lmP_max), this%VTpLM(lmP_max) )
      allocate( this%VxBrLM(lmP_max), this%VxBtLM(lmP_max), this%VxBpLM(lmP_max) )
      allocate( this%jxBrLM(lmP_max), this%jxBtLM(lmP_max), this%jxBpLM(lmP_max) )
      bytes_allocated = bytes_allocated + 9*lmP_max*SIZEOF_DEF_COMPLEX

      this%VTrLM(:)=zero
      this%VTtLM(:)=zero
      this%VTpLM(:)=zero
      this%VxBrLM(:)=zero
      this%VxBtLM(:)=zero
      this%VxBpLM(:)=zero
      this%jxBrLM(:)=zero
      this%jxBtLM(:)=zero
      this%jxBpLM(:)=zero

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(nonlinear_lm_t) :: this

      deallocate( this%VTrLM, this%VTtLM, this%VTpLM )
      deallocate( this%VxBrLM, this%VxBtLM, this%VxBpLM )
      deallocate( this%jxBrLM, this%jxBtLM, this%jxBpLM )

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine get_td(this,dVTrLM,dTdt,lrfLM,dVxBhLM,dBdt,djdt,n_r)
      !
      !  Purpose of this to calculate time derivatives dTdt
      !  and auxiliary arrays dVTrLM from non-linear terms
      !  in spectral form,contained in flmw1-3,flms1-3, flmb1-3 (input)
      !
    
      !-- Input of variables:
      class(nonlinear_lm_t) :: this
      integer, intent(in) :: n_r

      !-- Output of variables:
      complex(cp), intent(out) :: dVTrLM(:)
      complex(cp), intent(out) :: dTdt(:)
      complex(cp), intent(out) :: lrfLM(:,:)
      complex(cp), intent(out) :: dVxBhLM(:)
      complex(cp), intent(out) :: dBdt(:), djdt(:)

      !-- Local variables:
      integer :: l,m,lm,lmP,lmPS,lmPA
      real(cp) :: dLh
      complex(cp) :: dTdt_loc
      !complex(cp) :: LF_pol(lm_max), LF_tor(lm_max)

      !-- Spherically-symmetric contribution:
      dVTrLM(1)=this%VTrLM(1)
      dTdt(1)  =zero
      lrfLM(:,:)=zero
      this%VxBrLM(:) = zero
      this%VxBtLM(:) = zero
      this%VxBpLM(:) = zero

      !PERFON('td_heat')
      !$OMP PARALLEL DEFAULT(shared) &
      !$OMP private(lm,l,m,lmP,lmPS,lmPA,dTdt_loc)
      !LIKWID_ON('td_heat')
      !$OMP DO
      do lm=2,lm_max
         l   =lm2l(lm)
         m   =lm2m(lm)
         lmP =lm2lmP(lm)
         lmPS=lmP2lmPS(lmP)
         lmPA=lmP2lmPA(lmP)
         !------ This is horizontal heat advection:
         if ( l > m ) then
            dTdt_loc=-dTheta1S(lm)*this%VTtLM(lmPS) &
            &        +dTheta1A(lm)*this%VTtLM(lmPA) &
            &            -dPhi(lm)*this%VTpLM(lmP)
         else if ( l == m ) then
            dTdt_loc= dTheta1A(lm)*this%VTtLM(lmPA) &
            &         -dPhi(lm)*this%VTpLM(lmP)
         end if
         dVTrLM(lm)=this%VTrLM(lmP)
         dTdt(lm)  =dTdt_loc
         end do
      !$OMP end do
      !LIKWID_OFF('td_heat')
      !$OMP END PARALLEL
      !PERFOFF

      if ( l_mag_3D  ) then
         !PERFON('td_magnl')
         !$OMP PARALLEL do default(none) &
         !$OMP private(lm,l,m,lmP,lmPS,lmPA) &
         !$OMP shared(lm_max,lm2l,lm2m,lm2lmP,lmP2lmPS,lmP2lmPA) &
         !$OMP shared(dBdt,djdt,dTheta1S,dTheta1A,dPhi) &
         !$OMP shared(LF_pol,LF_tor,dVxBhLM) &
         !$OMP shared(dLh,or2_3D,r_3D,n_r,this)
         do lm=1,lm_max
            if ( lm == 1 ) then
               lmP=1
               lmPA=lmP2lmPA(lmP)
               dVxBhLM(lm)=-r_3D(n_r)*r_3D(n_r)           &
               &           *dTheta1A(lm)*this%VxBtLM(lmPA)
               dBdt(lm)   =-dTheta1A(lm)*this%VxBpLM(lmPA)
               djdt(lm)   = zero
               cycle
               if ( l_mag_LF ) then! .and. nR>n_r_LCR ) then
                  !LF_pol(lm)= or2_3D(n_r) *this%jxBrLM(lm)
                  !LF_tor(lm)=-dTheta1A(lm)*this%jxBpLM(lmPA)
                  lrfLM(lm,2)=-dTheta2A(lm)*this%jxBrLM(lmPA)
               end if
            end if
            l   =lm2l(lm)
            m   =lm2m(lm)
            lmP =lm2lmP(lm)
            lmPS=lmP2lmPS(lmP)
            lmPA=lmP2lmPA(lmP)
            dLh =real(l*(l+1),kind=cp)

            if ( l_mag_LF ) then
               !------- lrf_1= d/d r part of curl( curl(B) x B )_z
               !LF_pol(lm)= or2_3D(n_r)*this%jxBrLM(lmP)
               lrfLM(lm,1)= r_3D(n_r)*this%jxBtLM(lmP)
            end if

            !------- This is the radial part of the dynamo terms \curl(VxB)
            !PERFON('td_mnl1')
            if ( l > m ) then
               dBdt(lm)= dTheta1S(lm)*this%VxBpLM(lmPS) &
               &        -dTheta1A(lm)*this%VxBpLM(lmPA) &
               &        -dPhi(lm)    *this%VxBtLM(lmP)
               if ( l_mag_LF ) then
                  !------ When RMS values are required, the Lorentz force is treated
                  !       separately:
                  !------- LFTor= 1/(E*Pm) * curl( curl(B) x B )_r
                  !------- lrf_2= d/d theta part of curl( curl(B) x B )_z
                  !LF_tor(lm)= dTheta1S(lm)*this%jxBpLM(lmPS) &
                  !&          -dTheta1A(lm)*this%jxBpLM(lmPA) &
                  !&          -dPhi(lm)    *this%jxBtLM(lmP)
                  lrfLM(lm,2)= dTheta2S(lm)*this%jxBrLM(lmPS) &
                  &           -dTheta2A(lm)*this%jxBrLM(lmPA)
               end if
            else if ( l == m ) then
               dBdt(lm)=-dTheta1A(lm)*this%VxBpLM(lmPA) &
               &        -dPhi(lm)    *this%VxBtLM(lmP)
               if ( l_mag_LF ) then
                  !LF_tor(lm)=-dTheta1A(lm)*this%jxBpLM(lmPA) &
                  !&          -dPhi(lm)    *this%jxBtLM(lmP)
                  lrfLM(lm,2)=-dTheta2A(lm)*this%jxBrLM(lmPA)
               end if
            end if
            !PERFOFF

            !------- Radial component of
            !           \curl\curl(UxB) = \grad\div(VxB) - \laplace(VxB)
            !------- This is the radial part of \laplace (VxB)
            djdt(lm)= dLh*or2_3D(n_r)*or2_3D(n_r) &
            &        *this%VxBrLM(lmP)

            !------- This is r^2 * horizontal divergence of (VxB)
            !        Radial derivative performed in get_dr_td
            !PERFON('td_mnl2')
            if ( l > m ) then
               dVxBhLM(lm)= r_3D(n_r)*r_3D(n_r)*(          &
               &             dTheta1S(lm)*this%VxBtLM(lmPS) &
               &            -dTheta1A(lm)*this%VxBtLM(lmPA) &
               &            +dPhi(lm)    *this%VxBpLM(lmP) )
            else if ( l == m ) then
               dVxBhLM(lm)= r_3D(n_r)*r_3D(n_r)*(          &
               &            -dTheta1A(lm)*this%VxBtLM(lmPA) &
               &            +dPhi(lm)    *this%VxBpLM(lmP) )
            end if
            !PERFOFF
         end do
         !$OMP END PARALLEL DO
         !PERFOFF
      else
         do lm=1,lm_max
            dBdt(lm)   =zero
            djdt(lm)   =zero
            dVxBhLM(lm)=zero
         end do
      end if

   end subroutine get_td
!-----------------------------------------------------------------------------
end module nonlinear_lm_mod
