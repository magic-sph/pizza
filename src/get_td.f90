module nonlinear_lm_mod

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use blocking_lm, only: lm2l, lm2m, lm2lmP, lmP2lmPS, lmP2lmPA, lo_map
   use truncation_3D, only: lm_max, lmP_max
   use radial_functions, only: r_3D, or2_3D
   use horizontal, only: dTheta1S, dTheta1A, dPhi, dLh
   use namelists, only: l_heat_3D, l_mag_3D
   use constants, only: zero

   implicit none
 
   type :: nonlinear_lm_t
      !----- Nonlinear terms in lm-space: 
      complex(cp), allocatable :: VTtLM(:),  VTpLM(:)
      complex(cp), allocatable :: VxBrLM(:), VxBtLM(:), VxBpLM(:)

 
   contains
 
      procedure :: initialize
      procedure :: finalize
      procedure :: get_td
 
   end type nonlinear_lm_t

contains

   subroutine initialize(this,lm_max,lmP_max)

      class(nonlinear_lm_t) :: this
      integer, intent(in) :: lm_max, lmP_max

      if ( l_heat_3D ) then
         allocate( this%VTtLM(lm_max), this%VTpLM(lm_max) )
         bytes_allocated = bytes_allocated + 2*lm_max*SIZEOF_DEF_COMPLEX
      endif
      if ( l_mag_3D ) then
         allocate( this%VxBrLM(lmP_max), this%VxBtLM(lmP_max), this%VxBpLM(lmP_max) )
         bytes_allocated = bytes_allocated + 3*lmP_max*SIZEOF_DEF_COMPLEX
      endif

      if ( l_heat_3D ) then
         this%VTtLM(:)=zero
         this%VTpLM(:)=zero
      endif
      if ( l_mag_3D ) then
         this%VxBrLM(:)=zero
         this%VxBtLM(:)=zero
         this%VxBpLM(:)=zero
      endif

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(nonlinear_lm_t) :: this

      if ( l_heat_3D ) deallocate( this%VTtLM, this%VTpLM )
      if ( l_mag_3D ) deallocate( this%VxBrLM, this%VxBtLM, this%VxBpLM )

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine get_td(this,dTdt,dVxBhLM,dBdt,djdt,n_r)
      !
      !  Purpose of this to calculate time derivatives dTdt
      !  and auxiliary arrays dVTrLM from non-linear terms
      !  in spectral form,contained in flmw1-3,flms1-3, flmb1-3 (input)
      !
    
      !-- Input of variables:
      class(nonlinear_lm_t) :: this
      integer, intent(in) :: n_r

      !-- Output of variables:
      complex(cp), intent(out) :: dTdt(:)
      complex(cp), intent(out) :: dVxBhLM(:)
      complex(cp), intent(out) :: dBdt(:), djdt(:)

      !-- Local variables:
      integer :: l,m,lm,lmP,lmPS,lmPA

      if ( l_heat_3D ) then
         !-- Spherically-symmetric contribution:
         dTdt(1)  =zero

         !$OMP PARALLEL DEFAULT(shared) &
         !$OMP private(lm)
         !$OMP DO
         do lm=2,lm_max
            dTdt(lm)=dLh(lm)*this%VTtLM(lm)
         end do
         !$OMP end do
         !$OMP END PARALLEL
      end if ! l_heat_3D?

      if ( l_mag_3D  ) then
         !PERFON('td_magnl')
         !$OMP PARALLEL do default(none) &
         !$OMP private(lm,l,m,lmP,lmPS,lmPA) &
         !$OMP shared(lm_max,lm2l,lm2m,lm2lmP,lmP2lmPS,lmP2lmPA) &
         !$OMP shared(dBdt,djdt,dTheta1S,dTheta1A,dPhi) &
         !$OMP shared(LF_pol,LF_tor,dVxBhLM) &
         !$OMP shared(or2_3D,r_3D,n_r,this)
         !if ( if qst_to_spat get_nl% ) then !#ifdef SHTNS
         ! do lm=1,lm_max
         !    l   =lm2l(lm)
         !    lmP =lm2lmP(lm)
         !    dBdt(lm)   = dLh(lm)*this%VxBpLM(lmP)
         !    dVxBhLM(lm)=-dLh(lm)*this%VxBtLM(lmP)*r_3D(n_r)!*r_3D(n_r)
         !    djdt(lm)   = dLh(lm)*this%VxBrLM(lmP)*or2_3D(n_r)!*or2_3D(n_r)
         ! end do
         !else !#endifdef SHTNS
         do lm=1,lm_max
            if ( lm == 1 ) then
               lmP=1
               lmPA=lmP2lmPA(lmP)
               dVxBhLM(lm)=-r_3D(n_r)*r_3D(n_r)           &
               &           *dTheta1A(lm)*this%VxBtLM(lmPA)
               dBdt(lm)   =-dTheta1A(lm)*this%VxBpLM(lmPA)
               djdt(lm)   = zero
               cycle
            end if
            l   =lm2l(lm)
            m   =lm2m(lm)
            lmP =lm2lmP(lm)
            lmPS=lmP2lmPS(lmP)
            lmPA=lmP2lmPA(lmP)

            !------- This is the radial part of the dynamo terms \curl(VxB)
            !PERFON('td_mnl1')
            if ( l > m ) then
               dBdt(lm)= dTheta1S(lm)*this%VxBpLM(lmPS) &
               &        -dTheta1A(lm)*this%VxBpLM(lmPA) &
               &        -dPhi(lm)    *this%VxBtLM(lmP)
            else if ( l == m ) then
               dBdt(lm)=-dTheta1A(lm)*this%VxBpLM(lmPA) &
               &        -dPhi(lm)    *this%VxBtLM(lmP)
            end if
            !PERFOFF

            !------- Radial component of
            !           \curl\curl(UxB) = \grad\div(VxB) - \laplace(VxB)
            !------- This is the radial part of \laplace (VxB)
            djdt(lm)= dLh(lm)*this%VxBrLM(lmP)*or2_3D(n_r)*or2_3D(n_r)!

            !------- This is r^2 * horizontal divergence of (VxB)
            !        Radial derivative performed in get_dr_td
            !PERFON('td_mnl2')
            if ( l > m ) then
               dVxBhLM(lm)= r_3D(n_r)*r_3D(n_r)*(           &
               &             dTheta1S(lm)*this%VxBtLM(lmPS) &
               &            -dTheta1A(lm)*this%VxBtLM(lmPA) &
               &            +dPhi(lm)    *this%VxBpLM(lmP) )
            else if ( l == m ) then
               dVxBhLM(lm)= r_3D(n_r)*r_3D(n_r)*(           &
               &            -dTheta1A(lm)*this%VxBtLM(lmPA) &
               &            +dPhi(lm)    *this%VxBpLM(lmP) )
            end if
            !PERFOFF
         end do
         !$OMP END PARALLEL DO
         !PERFOFF
      end if ! l_mag_3D?

   end subroutine get_td
!-----------------------------------------------------------------------------
end module nonlinear_lm_mod
