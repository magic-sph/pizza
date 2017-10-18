module update_psi_integ

   use precision_mod
   use parallel_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: one, zero, ci, half
   use outputs, only: vp_bal_type
   use namelists, only: kbotv, ktopv, alpha, r_cmb, r_icb, l_non_rot, CorFac, &
       &                l_ek_pump
   use radial_functions, only: rscheme, or1, or2, beta, dbeta, ekpump, oheight, r
   use blocking, only: nMstart, nMstop, l_rank_has_m0
   use truncation, only: n_r_max, idx2m, m2idx
   use radial_der, only: get_ddr, get_dr
   use fields, only: work_Mloc
   use time_schemes, only: type_tscheme
   use useful, only: abortRun, roll
   use matrix_types, only: type_bandmat_complex, type_bordmat_complex,        &
   &                       type_bandmat_real
   use chebsparselib, only: intcheb4rmult4lapl2, intcheb4rmult4lapl,          &
       &                    intcheb4rmult4, rmult2, intcheb1rmult1,           &
       &                    intcheb2rmult2, intcheb4rmult4hmult8laplrot2,     &
       &                    intcheb4rmult4hmult8laplrot, intcheb4hmult2,      &
       &                    intcheb4rmult4hmult6, intcheb2rmult2hmult4,       &
       &                    intcheb2rmult2hmult4laplaxi, intcheb2, intcheb4


   implicit none
   
   private

   logical,  allocatable :: lPsimat(:)
   complex(cp), allocatable :: rhs(:)

   type(type_bordmat_complex), allocatable :: LHS_mat(:)
   type(type_bandmat_complex), allocatable :: RHSIL_mat(:)
   type(type_bandmat_real), allocatable :: RHSI_mat(:)
   type(type_bandmat_real) :: RHSE_mat(2)

   public :: update_psi_int, initialize_psi_integ, finalize_psi_integ, &
   &         get_psi_rhs_imp_int

contains

   subroutine initialize_psi_integ
      !
      ! Memory allocation
      !


      !-- Local variables
      integer :: n_m, m

      !-- Allocate array of matrices (this is handy since it can store various bandwidths)
      allocate( RHSI_mat(nMstart:nMstop) )
      allocate( RHSIL_mat(nMstart:nMstop) )
      allocate( LHS_mat(nMstart:nMstop) )

      !-- Initialize matrices
      if ( l_non_rot ) then
         call RHSE_mat(1)%initialize(3, 3, n_r_max) ! This is m  = 0
         call RHSE_mat(2)%initialize(4, 4, n_r_max) ! This is m /= 0
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m == 0 ) then
               call RHSI_mat(n_m)%initialize(4, 4, n_r_max)
               call RHSIL_mat(n_m)%initialize(2, 2, n_r_max)
               call LHS_mat(n_m)%initialize(4, 4, 2, n_r_max)
            else
               call RHSI_mat(n_m)%initialize(6, 6, n_r_max)
               call RHSIL_mat(n_m)%initialize(4, 4, n_r_max)
               call LHS_mat(n_m)%initialize(6, 6, 4, n_r_max)
            end if
         end do
      else
         call RHSE_mat(1)%initialize(2, 2, n_r_max) ! This is m  = 0
         call RHSE_mat(2)%initialize(6, 6, n_r_max) ! This is m /= 0
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m == 0 ) then
               call RHSI_mat(n_m)%initialize(8, 8, n_r_max)
               call RHSIL_mat(n_m)%initialize(6, 6, n_r_max)
               call LHS_mat(n_m)%initialize(8, 8, 2, n_r_max)
            else
               call RHSI_mat(n_m)%initialize(14, 14, n_r_max)
               call RHSIL_mat(n_m)%initialize(14, 14, n_r_max)
               call LHS_mat(n_m)%initialize(14, 14, 4, n_r_max)
            end if
         end do
      end if

      !-- Fill matrices
      call get_rhs_exp_mat(RHSE_mat(1),1)
      call get_rhs_exp_mat(RHSE_mat(2),2)
      do n_m=nMstart,nMstop
         m = idx2m(n_m)
         call get_rhs_imp_mat(RHSI_mat(n_m), m)
         call get_rhs_imp_lin_mat(RHSIL_mat(n_m), m)
      end do

      allocate( lPsimat(nMstart:nMstop) )
      lPsimat(:)=.false.
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_LOGICAL

      allocate( rhs(n_r_max) )
      rhs(:)=zero
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_COMPLEX

   end subroutine initialize_psi_integ
!------------------------------------------------------------------------------
   subroutine finalize_psi_integ
      !
      ! Memory deallocation
      ! 

      !-- Local variable
      integer :: n_m

      call RHSE_mat(2)%finalize()
      call RHSE_mat(1)%finalize()
      do n_m=nMstart,nMstop
         call RHSIL_mat(n_m)%finalize()
         call RHSI_mat(n_m)%finalize()
         call LHS_mat(n_m)%finalize()
      end do
      deallocate( LHS_mat, RHSIL_mat, RHSI_mat )

      deallocate( rhs)
      deallocate( lPsimat )

   end subroutine finalize_psi_integ
!------------------------------------------------------------------------------
   subroutine update_psi_int(psi_Mloc, dpsi_Mloc, d2psi_Mloc, om_Mloc, us_Mloc, up_Mloc,            &
              &              dVsOm_Mloc, dpsi_exp_Mloc, dpsi_imp_Mloc,       &
              &              buo_imp_Mloc, vp_bal, tscheme, lMat, l_roll_imp,&
              &              l_vphi_bal_calc)

      !-- Input variables
      type(type_tscheme), intent(in) :: tscheme
      logical,            intent(in) :: lMat
      logical,            intent(in) :: l_vphi_bal_calc
      logical,            intent(in) :: l_roll_imp

      !-- Output variables
      complex(cp),       intent(out) :: psi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),       intent(out) :: dpsi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),       intent(out) :: d2psi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),       intent(out) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),       intent(out) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),       intent(out) :: up_Mloc(nMstart:nMstop,n_r_max)
      type(vp_bal_type), intent(inout) :: vp_bal
      complex(cp),       intent(inout) :: dpsi_exp_Mloc(nMstart:nMstop,n_r_max,tscheme%norder_exp)
      complex(cp),       intent(inout) :: dVsOm_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),       intent(inout) :: dpsi_imp_Mloc(nMstart:nMstop,n_r_max,tscheme%norder_imp-1)
      complex(cp),       intent(inout) :: buo_imp_Mloc(nMstart:nMstop,n_r_max)

      !-- Local variables
      real(cp) :: uphi0(n_r_max), om0(n_r_max)
      integer :: n_r, n_m, n_cheb, m, n_o

      if ( lMat ) lPsimat(:)=.false.

      !-- Finish calculation of advection
      call get_dr( dVsOm_Mloc, work_Mloc, nMstart, nMstop, n_r_max, &
           &       rscheme, nocopy=.true.)

      !-- Finish calculation of the explicit part for current time step
      do n_r=1,n_r_max
         do n_m=nMstart, nMstop
            m = idx2m(n_m)
            if ( m /= 0 ) then
               dpsi_exp_Mloc(n_m,n_r,1)=dpsi_exp_Mloc(n_m,n_r,1)-work_Mloc(n_m,n_r)
            end if
         end do
      end do

      !-- Add Ekman pumping as an explicit term if this is requested
      if ( l_ek_pump ) then

         do n_r=2,n_r_max-1
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               if ( m == 0 ) then
                  dpsi_exp_Mloc(n_m,n_r,1)=       dpsi_exp_Mloc(n_m,n_r,1) -     &
                  &                           ekpump(n_r)*up_Mloc(n_m,n_r)

               else 
                  dpsi_exp_Mloc(n_m,n_r,1)=       dpsi_exp_Mloc(n_m,n_r,1) +     &
                  &                        ekpump(n_r)*( -om_Mloc(n_m,n_r) +     &
                  &                     half*beta(n_r)*   up_Mloc(n_m,n_r) +     &
                  &       beta(n_r)*(-ci*real(m,cp)+5.0_cp*r_cmb*oheight(n_r))*  &
                  &                                       us_Mloc(n_m,n_r) )
               end if
            end do
         end do

      end if

      if ( l_rank_has_m0 .and. l_vphi_bal_calc ) then
         do n_r=2,n_r_max
            vp_bal%rey_stress(n_r)=real(dpsi_exp_Mloc(m2idx(0),n_r,1))*or2(n_r)/ &
            &                      (r_cmb**2-r(n_r)**2)**2
         end do
         vp_bal%rey_stress(1)=0.0_cp
      end if

      !-- Transform the explicit part to chebyshev space
      call rscheme%costf1(dpsi_exp_Mloc(:,:,1), nMstart, nMstop, n_r_max)

      !-- Matrix-vector multiplication by the operator \int\int\int\int h^2
      do n_m=nMstart,nMstop
         m = idx2m(n_m)

         do n_r=1,n_r_max
            rhs(n_r)=dpsi_exp_Mloc(n_m,n_r,1)
         end do

         if ( m == 0 ) then
            call RHSE_mat(1)%mat_vec_mul(rhs)
         else
            call RHSE_mat(2)%mat_vec_mul(rhs)
            rhs(3)=zero
            rhs(4)=zero
         end if
         rhs(1)=zero
         rhs(2)=zero

         do n_r=1,n_r_max
            dpsi_exp_Mloc(n_m,n_r,1)=rhs(n_r)
         end do
      end do

      !-- Transform buoyancy to Chebyshev space
      call rscheme%costf1(buo_imp_Mloc, nMstart, nMstop, n_r_max)

      !-- Matrix-vector multiplication by the operator \int\int\int\int r^4 h^8
      do n_m=nMstart,nMstop
         m = idx2m(n_m)

         if ( m > 0 ) then
            do n_r=1,n_r_max
               rhs(n_r)=buo_imp_Mloc(n_m,n_r)
            end do

            call RHSE_mat(2)%mat_vec_mul(rhs)

            rhs(1)=zero
            rhs(2)=zero
            rhs(3)=zero
            rhs(4)=zero
            do n_r=1,n_r_max
               buo_imp_Mloc(n_m,n_r)=rhs(n_r)
            end do
         end if
      end do


      !-- Calculation of the implicit part
      call get_psi_rhs_imp_int(psi_Mloc, up_Mloc, tscheme%wimp_lin(2), &
           &                   dpsi_imp_Mloc(:,:,1), vp_bal, l_vphi_bal_calc)


      do n_m=nMstart,nMstop

         m = idx2m(n_m)

         if ( .not. lPsimat(n_m) ) then
            call get_lhs_mat(tscheme, LHS_mat(n_m), m)
            lPsimat(n_m)=.true.
         end if

         !-- Store vphi force balance if needed
         if ( m == 0 .and. l_vphi_bal_calc ) then
            do n_r=1,n_r_max
               vp_bal%dvpdt(n_r)=real(up_Mloc(n_m,n_r))/tscheme%dt(1)
            end do
         end if

         !-- Assemble RHS
         do n_r=1,n_r_max
            do n_o=1,tscheme%norder_imp-1
               if ( n_o == 1 ) then
                  rhs(n_r)=tscheme%wimp(n_o+1)*dpsi_imp_Mloc(n_m,n_r,n_o)
               else
                  rhs(n_r)=rhs(n_r)+tscheme%wimp(n_o+1)* &
                  &        dpsi_imp_Mloc(n_m,n_r,n_o)
               end if
            end do
            do n_o=1,tscheme%norder_exp
               rhs(n_r)=rhs(n_r)+tscheme%wexp(n_o)*dpsi_exp_Mloc(n_m,n_r,n_o)
            end do
            !-- Add buoyancy
            if ( m /= 0 ) then
               rhs(n_r)=rhs(n_r)+buo_imp_Mloc(n_m,n_r)
            end if
         end do

         call LHS_mat(n_m)%solve(rhs, n_r_max)

         if ( m == 0 ) then
            do n_cheb=1,rscheme%n_max
               uphi0(n_cheb)=real(rhs(n_cheb))
            end do
         else
            do n_cheb=1,rscheme%n_max
               psi_Mloc(n_m,n_cheb)=rhs(n_cheb)
            end do
         end if

      end do

      !-- set cheb modes > rscheme%n_max to zero (dealiazing)
      if ( rscheme%n_max < n_r_max ) then ! fill with zeros !
         do n_cheb=rscheme%n_max+1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               if ( m == 0 ) then
                  uphi0(n_cheb)=0.0_cp
               else
                  psi_Mloc(n_m,n_cheb)=zero
               end if
            end do
         end do
      end if

      !-- Bring uphi0 to the physical space
      if ( l_rank_has_m0 ) then
         call rscheme%costf1(uphi0, n_r_max)
         call get_dr(uphi0, om0, n_r_max, rscheme)

         if ( l_vphi_bal_calc ) then
            do n_r=1,n_r_max
               vp_bal%dvpdt(n_r)=uphi0(n_r)/tscheme%dt(1)-vp_bal%dvpdt(n_r)
            end do
         end if
      end if

      !-- Bring psi back the physical space
      call rscheme%costf1(psi_Mloc, nMstart, nMstop, n_r_max)

      !-- Get the radial derivative of psi to calculate uphi, us and omega
      call get_ddr(psi_Mloc, dpsi_Mloc, d2psi_Mloc, nMstart, nMstop, n_r_max, &
           &       rscheme)

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)

            if ( m == 0 ) then
               us_Mloc(n_m,n_r)=0.0_cp
               up_Mloc(n_m,n_r)=uphi0(n_r)
               om_Mloc(n_m,n_r)=om0(n_r)+or1(n_r)*uphi0(n_r)

               !dpsi_Mloc(n_m,n_r)=-uphi0(n_r)
               !d2psi_Mloc(n_m,n_r)=0.0_cp
            else
               us_Mloc(n_m,n_r)=ci*real(m,cp)*or1(n_r)*psi_Mloc(n_m,n_r)
               up_Mloc(n_m,n_r)=-dpsi_Mloc(n_m,n_r)-beta(n_r)*psi_Mloc(n_m,n_r)
               om_Mloc(n_m,n_r)=-d2psi_Mloc(n_m,n_r)-(or1(n_r)+beta(n_r))*       &
               &               dpsi_Mloc(n_m,n_r)-(or1(n_r)*beta(n_r)+dbeta(n_r) &
               &               -real(m,cp)*real(m,cp)*or2(n_r))*psi_Mloc(n_m,n_r)
            end if
         end do
      end do


      !-- Roll the explicit arrays before filling again the first block
      call roll(dpsi_exp_Mloc, nMstart, nMstop, n_r_max, tscheme%norder_exp)
      if ( l_roll_imp ) then
         call roll(dpsi_imp_Mloc, nMstart, nMstop, n_r_max, tscheme%norder_imp-1)
      end if

   end subroutine update_psi_int
!------------------------------------------------------------------------------
   subroutine get_psi_rhs_imp_int(psi_Mloc, up_Mloc, wimp, dpsi_imp_Mloc_last, &
              &                   vp_bal, l_vphi_bal_calc)

      !-- Input variables
      complex(cp), intent(in) :: psi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      real(cp),    intent(in) :: wimp
      logical,     intent(in) :: l_vphi_bal_calc

      !-- Output variables
      complex(cp), intent(inout) :: dpsi_imp_Mloc_last(nMstart:nMstop,n_r_max)
      type(vp_bal_type), intent(inout) :: vp_bal

      !-- Local variables
      real(cp) :: duphi0(n_r_max), d2uphi0(n_r_max), uphi0(n_r_max)
      integer :: n_r, n_m, m, m0

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m == 0 ) then
               dpsi_imp_Mloc_last(n_m,n_r)=up_Mloc(n_m,n_r)
            else
               dpsi_imp_Mloc_last(n_m,n_r)=psi_Mloc(n_m,n_r)
            end if
         end do
      end do

      !-- Transform the implicit part to chebyshev space
      call rscheme%costf1(dpsi_imp_Mloc_last, nMstart, nMstop, n_r_max)

      !-- Matrix-vector multiplication by the operator -\int^4 r^4 \Delta .
      do n_m=nMstart,nMstop
         m = idx2m(n_m)

         do n_r=1,n_r_max
            rhs(n_r)= dpsi_imp_Mloc_last(n_m,n_r)
         end do

         call RHSI_mat(n_m)%mat_vec_mul(rhs)

         rhs(1)=zero ! vphi equation has only 2 BCs
         rhs(2)=zero
         if ( m > 0 ) then
            rhs(3)=zero
            rhs(4)=zero
         end if

         do n_r=1,n_r_max
            dpsi_imp_Mloc_last(n_m,n_r)=rhs(n_r)
         end do

      end do

      !-- Calculate and store vphi force balance if needed
      if ( l_rank_has_m0 .and. l_vphi_bal_calc ) then
         m0 = m2idx(0)
         do n_r=1,n_r_max
            uphi0(n_r)=real(up_Mloc(m0, n_r),kind=cp)
         end do
         call get_ddr(uphi0, duphi0, d2uphi0, n_r_max, rscheme)
         do n_r=1,n_r_max
            vp_bal%visc(n_r)=d2uphi0(n_r)+or1(n_r)*duphi0(n_r)-&
            &                or2(n_r)*uphi0(n_r)
            vp_bal%pump(n_r)=-ekpump(n_r)*uphi0(n_r)
         end do
      end if

      if ( wimp /= 0.0_cp ) then

         !-- Copy psi_Mloc into work_Mloc
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               if ( m == 0 ) then
                  work_Mloc(n_m,n_r)=up_Mloc(n_m,n_r)
               else
                  work_Mloc(n_m,n_r)=psi_Mloc(n_m,n_r)
               end if
            end do
         end do

         !-- Transform work_Mloc to Cheb space
         call rscheme%costf1(work_Mloc, nMstart, nMstop, n_r_max)

         !-- Matrix-vector multiplication by the LHS operator
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            do n_r=1,n_r_max
               rhs(n_r)=work_Mloc(n_m,n_r)
            end do
            call RHSIL_mat(n_m)%mat_vec_mul(rhs)
            rhs(1)=zero ! vphi equation has only 2 BCs
            rhs(2)=zero
            if ( m > 0 ) then
               rhs(3)=zero
               rhs(4)=zero
            end if
            do n_r=1,n_r_max
               work_Mloc(n_m,n_r)=rhs(n_r)
            end do
         end do

         !-- Finally assemble the right hand side
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               dpsi_imp_Mloc_last(n_m,n_r)=dpsi_imp_Mloc_last(n_m,n_r) &
               &                              +wimp*work_Mloc(n_m,n_r)
            end do
         end do

      end if ! if wimp /= 0

   end subroutine get_psi_rhs_imp_int
!------------------------------------------------------------------------------
   subroutine get_lhs_mat(tscheme, A_mat, m)

      !-- Input variables
      type(type_tscheme), intent(in) :: tscheme    ! time step
      integer,            intent(in) :: m          ! Azimuthal wavenumber

      !-- Output variables
      type(type_bordmat_complex), intent(inout) :: A_mat

      !-- Local variables
      real(cp) :: stencilA4(A_mat%nbands), CorSten(A_mat%nbands)
      integer :: n_r, i_r, n_band, n_b
      real(cp) :: a, b


      !-- We have to fill A3 with zeros again otherwise on the next iteration
      !-- with a different dt there might be some issues with spurious values
      do n_r=1,A_mat%nlines_band
         do n_b=1,A_mat%ntau
            A_mat%A3(n_r,n_b)=zero
         end do
      end do

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      !-- Fill A4 banded-block
      do n_r=1,A_mat%nlines_band
         i_r = n_r+A_mat%ntau

         !-- Define the equations
         if ( l_non_rot ) then
            if ( m == 0 ) then
               stencilA4 = intcheb2rmult2(a,b,i_r-1,A_mat%nbands)-               &
               &           tscheme%wimp_lin(1)*( rmult2(a,b,i_r-1,A_mat%nbands)- &
               &                  3.0_cp*intcheb1rmult1(a,b,i_r-1,A_mat%nbands) )
               CorSten   = 0.0_cp
            else
               stencilA4 = -intcheb4rmult4lapl(a,b,m,i_r-1,A_mat%nbands)+  &
               &           tscheme%wimp_lin(1)*                            &
               &           intcheb4rmult4lapl2(a,b,m,i_r-1,A_mat%nbands)  
               CorSten   = 0.0_cp
            end if
         else
            if ( m == 0 ) then
               stencilA4 = intcheb2rmult2hmult4(a,b,r_cmb,i_r-1,A_mat%nbands)-  &
               &           tscheme%wimp_lin(1)*(                                &
               &           intcheb2rmult2hmult4laplaxi(a,b,r_cmb,i_r-1,A_mat%nbands) )
               CorSten   = 0.0_cp
            else
               stencilA4 = -intcheb4rmult4hmult8laplrot(a,b,r_cmb,m,i_r-1, &
               &                                        A_mat%nbands)   +  &
               &           tscheme%wimp_lin(1)*                            &
               &           intcheb4rmult4hmult8laplrot2(a,b,r_cmb,m,i_r-1, &
               &                                        A_mat%nbands)  
               CorSten   = tscheme%wimp_lin(1)*CorFac*real(m,cp)*        &
               &           intcheb4rmult4hmult6(a,b,r_cmb,i_r-1,A_mat%nbands)
            end if
         end if

         !-- Roll the array for band storage
         do n_band=1,A_mat%nbands
            if ( n_r+A_mat%ku+1-n_band <= A_mat%nlines_band .and. n_r+A_mat%ku+1-n_band >= 1 ) then
               A_mat%A4(A_mat%kl+n_band,n_r+A_mat%ku+1-n_band) = rscheme%rnorm* &
               &             cmplx(stencilA4(n_band), CorSten(n_band), kind=cp)
            end if
         end do
      end do

      !-- Fill A3
      do n_r=1,A_mat%nlines_band
         i_r = n_r+A_mat%ntau

         if ( l_non_rot ) then
            if ( m == 0 ) then
               stencilA4 = intcheb2rmult2(a,b,i_r-1,A_mat%nbands)-               &
               &           tscheme%wimp_lin(1)*( rmult2(a,b,i_r-1,A_mat%nbands)- &
               &                  3.0_cp*intcheb1rmult1(a,b,i_r-1,A_mat%nbands) )
               CorSten   = 0.0_cp
            else
               stencilA4 = -intcheb4rmult4lapl(a,b,m,i_r-1,A_mat%nbands)+  &
               &           tscheme%wimp_lin(1)*                            &
               &           intcheb4rmult4lapl2(a,b,m,i_r-1,A_mat%nbands)  
               CorSten   = 0.0_cp
            end if
         else
            if ( m == 0 ) then
               stencilA4 = intcheb2rmult2hmult4(a,b,r_cmb,i_r-1,A_mat%nbands)-  &
               &           tscheme%wimp_lin(1)*(                                &
               &           intcheb2rmult2hmult4laplaxi(a,b,r_cmb,i_r-1,A_mat%nbands) )
               CorSten   = 0.0_cp
            else
               stencilA4 = -intcheb4rmult4hmult8laplrot(a,b,r_cmb,m,i_r-1, &
               &                                        A_mat%nbands)   +  &
               &           tscheme%wimp_lin(1)*                            &
               &           intcheb4rmult4hmult8laplrot2(a,b,r_cmb,m,i_r-1, &
               &                                        A_mat%nbands)  
               CorSten   = tscheme%wimp_lin(1)*CorFac*real(m,cp)*        &
               &           intcheb4rmult4hmult6(a,b,r_cmb,i_r-1,A_mat%nbands)
            end if
         end if

         !-- Only the lower bands can contribute to the matrix A3
         do n_band=1,A_mat%kl
            if ( n_r <= n_band .and. n_r+A_mat%ntau-n_band >= 1 ) then
               A_mat%A3(n_r,n_r+A_mat%ntau-n_band) = rscheme%rnorm*       &
               &      cmplx(stencilA4(A_mat%ku+1+n_band),                 &
               &              CorSten(A_mat%ku+1+n_band), kind=cp)
            end if
         end do
      end do

      !-- Add the tau Lines for boundary conditions into A1 and A2
      if ( m == 0 ) then

         do n_r=1,A_mat%nlines
            if ( n_r <= A_mat%ntau ) then
               if ( ktopv == 1 ) then ! Stress free
                  A_mat%A1(1,n_r)=rscheme%rnorm*( rscheme%drMat(1,n_r)-or1(1)*&
                  &                                rscheme%rMat(1,n_r) )
               else ! Rigid
                  A_mat%A1(1,n_r)=rscheme%rnorm*rscheme%rMat(1,n_r)
               end if
               if ( kbotv == 1 ) then ! Stress free
                  A_mat%A1(2,n_r)=rscheme%rnorm*( rscheme%drMat(2,n_r) &
                  &                  -or1(n_r_max)*rscheme%rMat(2,n_r) )
               else ! Rigid
                  A_mat%A1(2,n_r)=rscheme%rnorm*rscheme%rMat(2,n_r)
               end if
            else

               if ( ktopv == 1 ) then
                  A_mat%A2(1,n_r-A_mat%ntau)=rscheme%rnorm*(rscheme%drMat(1,n_r)-&
                  &                                   or1(1)*rscheme%rMat(1,n_r) )
               else
                  A_mat%A2(1,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%rMat(1,n_r)
               end if
               if ( kbotv == 1 ) then
                  A_mat%A2(2,n_r-A_mat%ntau)=rscheme%rnorm*(             &
                  &                               rscheme%drMat(2,n_r)-  &
                  &                   or1(n_r_max)*rscheme%rMat(2,n_r) )

               else
                  A_mat%A2(2,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%rMat(2,n_r)
               end if
            end if
         end do

      else

         do n_r=1,A_mat%nlines
            if ( n_r <= A_mat%ntau ) then
               A_mat%A1(1,n_r)=rscheme%rnorm*rscheme%rMat(1,n_r)
               A_mat%A1(2,n_r)=rscheme%rnorm*rscheme%rMat(2,n_r)
               if ( ktopv == 1 ) then !-- Stress-free
                  A_mat%A1(3,n_r)=rscheme%rnorm*(rscheme%d2rMat(1,n_r)-&
                  &                        or1(1)*rscheme%drMat(1,n_r) ) 
               else
                  A_mat%A1(3,n_r)=rscheme%rnorm*rscheme%drMat(1,n_r)
               end if
               if ( kbotv == 1 ) then
                  A_mat%A1(4,n_r)=rscheme%rnorm*(rscheme%d2rMat(2,n_r)-&
                  &                  or1(n_r_max)*rscheme%drMat(2,n_r) )
               else
                  A_mat%A1(4,n_r)=rscheme%rnorm*rscheme%drMat(2,n_r)
               end if
            else
               A_mat%A2(1,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%rMat(1,n_r)
               A_mat%A2(2,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%rMat(2,n_r)
               if ( ktopv == 1 ) then
                  A_mat%A2(3,n_r-A_mat%ntau)=rscheme%rnorm*(rscheme%d2rMat(1,n_r)-&
                  &                                   or1(1)*rscheme%drMat(1,n_r) )
               else
                  A_mat%A2(3,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%drMat(1,n_r)
               end if
               if ( kbotv == 1 ) then
                  A_mat%A2(4,n_r-A_mat%ntau)=rscheme%rnorm*(             &
                  &                              rscheme%d2rMat(2,n_r)-  &
                  &                  or1(n_r_max)*rscheme%drMat(2,n_r) )
               else
                  A_mat%A2(4,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%drMat(2,n_r)
               end if
            end if
         end do
      end if

      !-- Cheb factor for boundary conditions
      do n_r=1,A_mat%ntau
         A_mat%A1(n_r,1)                =rscheme%boundary_fac*A_mat%A1(n_r,1)
         A_mat%A2(n_r,A_mat%nlines_band)=rscheme%boundary_fac*A_mat%A2(n_r,A_mat%nlines_band)
      end do

      !-- LU factorisation
      call A_mat%prepare_LU()

   end subroutine get_lhs_mat
!------------------------------------------------------------------------------
   subroutine get_rhs_exp_mat(D_mat, m0)
      !
      ! This corresponds to the matrix that goes in front of the explicit terms
      !

      !-- Input variable
      integer, intent(in) :: m0

      !-- Output variable
      type(type_bandmat_real), intent(inout) :: D_mat

      !-- Local variables
      real(cp) :: stencilD(D_mat%nbands)
      real(cp) :: a, b
      integer :: n_band, n_r, i_r, n_bounds

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      if ( m0 == 1 ) then ! This is the axisymmetric part
         n_bounds = 2
      else                ! This is the non-axisymmetric part
         n_bounds = 4
      end if

      !-- Fill right-hand side matrix
      do n_r=1,D_mat%nlines
         i_r = n_r+n_bounds

         !-- Define right-hand side equations
         if ( l_non_rot ) then
            if ( m0 == 1) then
               stencilD = intcheb2(a,i_r-1,D_mat%nbands)
            else
               stencilD = intcheb4(a,i_r-1,D_mat%nbands)
            end if
         else
            if ( m0 == 1) then
               stencilD = intcheb2(a,i_r-1,D_mat%nbands)
            else
               stencilD = intcheb4hmult2(a,b,r_cmb,i_r-1,D_mat%nbands)
            end if
         end if

         !-- Roll array for band storage
         do n_band=1,D_mat%nbands
            if ( i_r+D_mat%ku+1-n_band <= D_mat%nlines .and. i_r+D_mat%ku+1-n_band >= 1 ) then
               D_mat%dat(n_band,i_r+D_mat%ku+1-n_band) = rscheme%rnorm*stencilD(n_band)
            end if
         end do
      end do

   end subroutine get_rhs_exp_mat
!------------------------------------------------------------------------------
   subroutine get_rhs_imp_mat(B_mat, m)
      !
      ! This matrix corresponds to the that goes in front of the d/dt term
      ! in the right-hand-side
      !

      !-- Input variable
      integer, intent(in) :: m

      !-- Output variable
      type(type_bandmat_real), intent(inout) :: B_mat

      !-- Local variables
      real(cp) :: stencilB(B_mat%nbands)
      real(cp) :: a, b
      integer :: n_band, n_r, i_r, n_bounds

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      if ( m == 0 ) then
         n_bounds = 2 ! For vphi equation, only 2 BCs
      else
         n_bounds = 4 ! For psi equation, 4 BCs
      end if

      !-- Fill right-hand side matrix
      do n_r=1,B_mat%nlines
         i_r = n_r+n_bounds

         !-- Define right-hand side equations
         if ( l_non_rot ) then
            if ( m == 0 ) then
               stencilB = intcheb2rmult2(a,b,i_r-1,B_mat%nbands)
            else
               stencilB = -intcheb4rmult4lapl(a,b,m,i_r-1,B_mat%nbands)
            end if
         else
            if ( m == 0 ) then
               stencilB = intcheb2rmult2hmult4(a,b,r_cmb,i_r-1,B_mat%nbands)
            else
               stencilB = -intcheb4rmult4hmult8laplrot(a,b,r_cmb,m,i_r-1,&
               &                                       B_mat%nbands)
            end if
         end if

         !-- Roll array for band storage
         do n_band=1,B_mat%nbands
            if ( i_r+B_mat%ku+1-n_band <= B_mat%nlines .and. i_r+B_mat%ku+1-n_band >= 1 ) then
               B_mat%dat(n_band,i_r+B_mat%ku+1-n_band) = rscheme%rnorm*stencilB(n_band)
            end if
         end do
      end do

   end subroutine get_rhs_imp_mat
!------------------------------------------------------------------------------
   subroutine get_rhs_imp_lin_mat(Cmat, m)
      !
      ! This matrix corresponds to the linear part that goes on the right-hand-side
      !

      !-- Input variable
      integer, intent(in) :: m

      !-- Output variable
      type(type_bandmat_complex), intent(inout) :: Cmat

      !-- Local variables
      real(cp) :: stencilC(Cmat%nbands), CorSten(Cmat%nbands)
      real(cp) :: a, b
      integer :: n_band, n_r, i_r, n_bounds

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      if ( m == 0 ) then
         n_bounds = 2 ! vphi equation has 2 BCs
      else
         n_bounds = 4 ! psi equation has 4 BCs
      end if

      !-- Fill right-hand side matrix
      do n_r=1,Cmat%nlines
         i_r = n_r+n_bounds

         !-- Define right-hand side equations
         if ( l_non_rot ) then
            if ( m == 0 ) then
               stencilC =               rmult2(a,b,i_r-1,Cmat%nbands)- &
               &         3.0_cp*intcheb1rmult1(a,b,i_r-1,Cmat%nbands)
               CorSten  = 0.0_cp
            else
               stencilC = -intcheb4rmult4lapl2(a,b,m,i_r-1,Cmat%nbands)
               CorSten  = 0.0_cp
            end if
         else
            if ( m == 0 ) then
               stencilC = intcheb2rmult2hmult4laplaxi(a,b,r_cmb,i_r-1,Cmat%nbands)
               CorSten  = 0.0_cp
            else
               stencilC = -intcheb4rmult4hmult8laplrot2(a,b,r_cmb,m,i_r-1,&
               &                                        Cmat%nbands)
               CorSten  = -CorFac*real(m,cp)*intcheb4rmult4hmult6(a,b,r_cmb, &
               &                                                  i_r-1,Cmat%nbands)
            end if
         end if

         !-- Roll array for band storage
         do n_band=1,Cmat%nbands
            if ( i_r+Cmat%ku+1-n_band <= Cmat%nlines .and. i_r+Cmat%ku+1-n_band >= 1 ) then
               Cmat%dat(n_band,i_r+Cmat%ku+1-n_band) = rscheme%rnorm*       &
               &            cmplx(stencilC(n_band), CorSten(n_band), kind=cp)
            end if
         end do
      end do

   end subroutine get_rhs_imp_lin_mat
!------------------------------------------------------------------------------
end module update_psi_integ
