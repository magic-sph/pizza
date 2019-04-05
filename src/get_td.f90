module nonlinear_lm_mod

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use blocking_lm, only: lm2l, lm2m, lm2lmP, lmP2lmPS, lmP2lmPA
   use truncation_3D, only: lm_max, lmP_max
   use horizontal, only: dTheta1S, dTheta1A, dPhi
   use constants, only: zero, two

   implicit none
 
   type :: nonlinear_lm_t
      !----- Nonlinear terms in lm-space: 
      complex(cp), allocatable :: VTrLM(:),  VTtLM(:),  VTpLM(:)
 
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
      bytes_allocated = bytes_allocated + 3*lmP_max*SIZEOF_DEF_COMPLEX

      this%VTrLM(:)=zero
      this%VTtLM(:)=zero
      this%VTpLM(:)=zero

   end subroutine initialize
!----------------------------------------------------------------------------
   subroutine finalize(this)

      class(nonlinear_lm_t) :: this

      deallocate( this%VTrLM, this%VTtLM, this%VTpLM )

   end subroutine finalize
!----------------------------------------------------------------------------
   subroutine get_td(this,dVTrLM,dTdt)
      !
      !  Purpose of this to calculate time derivatives dTdt
      !  and auxiliary arrays dVTrLM from non-linear terms
      !  in spectral form,contained in flmw1-3,flms1-3, flmb1-3 (input)
      !
    
      !-- Input of variables:
      class(nonlinear_lm_t) :: this

      !-- Output of variables:
      complex(cp), intent(out) :: dVTrLM(:)
      complex(cp), intent(out) :: dTdt(:)
    
      !-- Local variables:
      integer :: l,m,lm,lmS,lmA,lmP,lmPS,lmPA
      complex(cp) :: dTdt_loc
    
      dVTrLM(1)=this%VTrLM(1)
      dTdt(1)  =zero

      !PERFON('td_heat')
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP private(lm,l,m,lmP,lmPS,lmPA,dTdt_loc) &
      !$OMP shared(lm2l,lm2m,lm2lmP,lmP2lmPS,lmP2lmPA) &
      !$OMP shared(lm_max,dTdt,dVTrLM,dTheta1S,dTheta1A,dPhi) &
      !$OMP shared(this)
      !LIKWID_ON('td_heat')
      !$OMP DO
      do lm=2,lm_max
         l   =lm2l(lm)
         m   =lm2m(lm)
         lmP =lm2lmP(lm)
         lmPS=lmP2lmPS(lmP)
         lmPA=lmP2lmPA(lmP)
         !------ This is horizontal heat advection:
         !PERFON('td_h1')
         if ( l > m ) then
            dTdt_loc= -dTheta1S(lm)*this%VTtLM(lmPS) &
            &         +dTheta1A(lm)*this%VTtLM(lmPA) &
            &         -dPhi(lm)*this%VTpLM(lmP)
         else if ( l == m ) then
            dTdt_loc=  dTheta1A(lm)*this%VTtLM(lmPA) &
            &          -dPhi(lm)*this%VTpLM(lmP)
         end if
            !PERFOFF
            !-----   simplified form for linear onset !
            !        not ds not saved in the current program form!
            !                 dTdt(lm)=
            !                    -dLh(lm)*w(lm,nR)*or2(nR)*dsR(1)
            dVTrLM(lm)=this%VTrLM(lmP)
            dTdt(lm) = dTdt_loc
         end do
      !$OMP end do
      !LIKWID_OFF('td_heat')
      !$OMP END PARALLEL
      !PERFOFF

   end subroutine get_td
!-----------------------------------------------------------------------------
end module nonlinear_lm_mod
