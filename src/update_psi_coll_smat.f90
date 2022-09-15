module update_psi_coll_smat

   use precision_mod
   use parallel_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: one, zero, ci, half, third
   use horizontal, only: hdif_V
   use namelists, only: kbotv, ktopv, alpha, r_cmb, CorFac, ViscFac, &
       &                l_coriolis_imp, l_buo_imp, l_ek_pump,        &
       &                l_non_rot, l_mag_LF, l_leibniz, l_QG_basis, damp_zon
   use radial_functions, only: rscheme, or1, or2, beta, dbeta, ekpump, oheight
   use blocking, only: nMstart, nMstop, l_rank_has_m0
   use truncation, only: n_r_max, idx2m, m2idx
   use radial_der, only: get_ddr, get_dr
   use fields, only: work_Mloc
   use algebra, only: prepare_full_mat, solve_full_mat
   use time_schemes, only: type_tscheme
   use vort_balance, only: vort_bal_type
   use vp_balance, only: vp_bal_type
   use time_array
   use useful, only: abortRun
   use timers_mod, only: timers_type

   implicit none
   
   private

   logical,  allocatable :: lPsimat(:)
   complex(cp), allocatable :: psiMat(:,:,:)
   real(cp), allocatable :: uphiMat(:,:)
   integer,  allocatable :: psiPivot(:,:)
   real(cp), allocatable :: psiMat_fac(:,:,:)
   complex(cp), allocatable :: rhs(:)
   real(cp), allocatable :: rhs_m0(:)

   public :: update_om_coll_smat, initialize_om_coll_smat, finalize_om_coll_smat, &
   &         get_psi_rhs_imp_coll_smat, finish_exp_psi_coll_smat

contains

   subroutine initialize_om_coll_smat

      allocate( lPsimat(nMstart:nMstop) )
      lPsimat(:)=.false.
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_LOGICAL

      allocate( psiMat(2*n_r_max, 2*n_r_max, nMstart:nMstop) )
      allocate( psiPivot(2*n_r_max, nMstart:nMstop) )
      allocate( psiMat_fac(2*n_r_max, 2, nMstart:nMstop) )
      allocate( rhs(2*n_r_max), rhs_m0(n_r_max) )

      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*4*n_r_max*n_r_max* &
      &                 SIZEOF_DEF_COMPLEX+2*n_r_max*(nMstop-nMstart+1)*      &
      &                 SIZEOF_INTEGER+ n_r_max*(3+4*(nMstop-nMstart+1))*     &
      &                 SIZEOF_DEF_REAL

      allocate( uphiMat(n_r_max,n_r_max) )
      bytes_allocated = bytes_allocated+n_r_max*n_r_max*SIZEOF_DEF_REAL

   end subroutine initialize_om_coll_smat
!------------------------------------------------------------------------------
   subroutine finalize_om_coll_smat

      deallocate( rhs_m0, rhs, psiMat_fac )
      deallocate( lPsimat, psiMat, uphiMat, psiPivot )

   end subroutine finalize_om_coll_smat
!------------------------------------------------------------------------------
   subroutine update_om_coll_smat(psi_Mloc, om_Mloc, dom_Mloc, us_Mloc, up_Mloc, &
              &                   buo_Mloc, lf_Mloc, dpsidt, vp_bal, vort_bal, &
              &                   tscheme, lMat, timers)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lMat
      complex(cp),         intent(in) :: buo_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(in) :: lf_Mloc(nMstart:nMstop,n_r_max)

      !-- Output variables
      complex(cp),         intent(out) :: psi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: dom_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: up_Mloc(nMstart:nMstop,n_r_max)
      type(vp_bal_type),   intent(inout) :: vp_bal
      type(vort_bal_type), intent(inout) :: vort_bal
      type(type_tarray),   intent(inout) :: dpsidt
      type(timers_type),   intent(inout) :: timers

      !-- Local variables
      real(cp) :: uphi0(n_r_max), om0(n_r_max), runStart, runStop
      integer :: n_r, n_m, n_cheb, m

      if ( lMat ) lPsimat(:)=.false.


      !-- Calculation of the implicit part
      call get_psi_rhs_imp_coll_smat(us_Mloc, up_Mloc, om_Mloc, dom_Mloc,    &
           &                         dpsidt%old(:,:,tscheme%istage),         &
           &                         dpsidt%impl(:,:,tscheme%istage),        &
           &                         vp_bal, vort_bal,                       &
           &                         tscheme%l_imp_calc_rhs(tscheme%istage))

      !-- Calculate first part of time-derivative of \omega if needed
      if ( vort_bal%l_calc .and. tscheme%istage == 1  ) then
         call vort_bal%initialize_domdt(om_Mloc,tscheme)
      end if

      !-- Calculate first part of time-derivative of axisymmetric u_\phi if needed
      if ( vp_bal%l_calc .and. tscheme%istage == 1  ) then
         call vp_bal%initialize_dvpdt(up_Mloc,tscheme)
      end if

      !-- Now assemble the right hand side and store it in work_Mloc
      call tscheme%set_imex_rhs(work_Mloc, dpsidt, nMstart, nMstop, n_r_max)

      do n_m=nMstart,nMstop

         m = idx2m(n_m)

         if ( m == 0 ) then ! Axisymmetric component

            if ( .not. lPsimat(n_m) ) then
               call get_uphiMat(tscheme, uphiMat(:,:), psiPivot(1:n_r_max,n_m))
               lPsimat(n_m)=.true.
            end if

            rhs_m0(1)       = 0.0_cp
            rhs_m0(n_r_max) = 0.0_cp
            do n_r=2,n_r_max-1
               rhs_m0(n_r)=real(work_Mloc(n_m,n_r),kind=cp)
            end do

            if ( vp_bal%l_calc ) then
               do n_r=1,n_r_max
                  vp_bal%rey_stress(n_r)=real(dpsidt%expl(n_m,n_r,tscheme%istage))
                  if ( l_mag_LF ) then
                     !vp_bal%lorentz_force(n_r)=0.0_cp
                     vp_bal%lorentz_force(n_r)=real(lf_Mloc(n_m,n_r))
                     vp_bal%rey_stress(n_r)=vp_bal%rey_stress(n_r) - real(lf_Mloc(n_m,n_r))
                  end if
               end do
            end if

            call solve_full_mat(uphiMat(:,:), n_r_max, n_r_max,   &
                 &              psiPivot(1:n_r_max,n_m), rhs_m0(:))

            do n_cheb=1,rscheme%n_max
               uphi0(n_cheb)=rhs_m0(n_cheb)
            end do

         else ! Non-axisymmetric components
         
            if ( .not. lPsimat(n_m) ) then
               call get_psiMat(tscheme, m, psiMat(:,:,n_m), psiPivot(:,n_m), &
                    &          psiMat_fac(:,:,n_m), timers%lu, timers%n_lu_calls)
               lPsimat(n_m)=.true.
            end if

            if ( vort_bal%l_calc .and. tscheme%istage == tscheme%nstages ) then
               do n_r=1,n_r_max
                  vort_bal%buo(n_m,n_r)=buo_Mloc(n_m,n_r)/tscheme%dt(1)
               end do
            end if

            rhs(1)        =zero
            rhs(n_r_max)  =zero
            rhs(n_r_max+1)=zero
            rhs(2*n_r_max)=zero
            do n_r=2,n_r_max-1
               !-- Add buoyancy
               if ( l_buo_imp ) then
                  rhs(n_r)=work_Mloc(n_m,n_r)+buo_Mloc(n_m,n_r)
               else
                  rhs(n_r)=work_Mloc(n_m,n_r)
               end if
               !-- Second part is zero (no time-advance in the psi-block)
               rhs(n_r+n_r_max)=zero
            end do

            do n_r=1,2*n_r_max
               rhs(n_r) = rhs(n_r)*psiMat_fac(n_r,1,n_m)
            end do
            runStart = MPI_Wtime()
            call solve_full_mat(psiMat(:,:,n_m), 2*n_r_max, 2*n_r_max, &
                 &              psiPivot(:, n_m), rhs(:))
            runStop = MPI_Wtime()
            if ( runStop > runStart ) then
               timers%solve = timers%solve + (runStop-runStart)
               timers%n_solve_calls = timers%n_solve_calls+1
            end if
            do n_r=1,2*n_r_max
               rhs(n_r) = rhs(n_r)*psiMat_fac(n_r,2,n_m)
            end do

            do n_cheb=1,rscheme%n_max
               om_Mloc(n_m,n_cheb) =rhs(n_cheb)
               psi_Mloc(n_m,n_cheb)=rhs(n_cheb+n_r_max)
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
                  om_Mloc(n_m,n_cheb) =zero
                  psi_Mloc(n_m,n_cheb)=zero
               end if
            end do
         end do
      end if

      !-- Bring uphi0 to the physical space
      if ( l_rank_has_m0 ) then
         call get_dr(uphi0, om0, n_r_max, rscheme, l_dct=.false.)
         call rscheme%costf1(uphi0, n_r_max)
      end if

      !-- Get the radial derivative of psi to calculate uphi
      call get_dr(psi_Mloc, work_Mloc, nMstart, nMstop, n_r_max, &
           &      rscheme, l_dct_in=.false.)

      !-- Bring psi and omega to the physical space
      runStart = MPI_Wtime()
      call rscheme%costf1(psi_Mloc, nMstart, nMstop, n_r_max)
      call rscheme%costf1(om_Mloc, nMstart, nMstop, n_r_max)
      runStop = MPI_Wtime()
      if ( runStop > runStart ) then
         timers%dct = timers%dct + (runStop-runStart)
         timers%n_dct_calls = timers%n_dct_calls + 2
      end if

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)

            if ( m == 0 ) then
               us_Mloc(n_m,n_r)=0.0_cp
               up_Mloc(n_m,n_r)=uphi0(n_r)*damp_zon
               om_Mloc(n_m,n_r)=om0(n_r)+or1(n_r)*uphi0(n_r)*damp_zon
            else
               us_Mloc(n_m,n_r)=ci*real(m,cp)*or1(n_r)*psi_Mloc(n_m,n_r)
               up_Mloc(n_m,n_r)=-work_Mloc(n_m,n_r)-beta(n_r)*psi_Mloc(n_m,n_r)
            end if
         end do
      end do

      !-- Roll the time arrays before filling again the first block
      call tscheme%rotate_imex(dpsidt, nMstart, nMstop, n_r_max)

      !-- Finish calculation of du_\phi/dt if requested
      if ( vp_bal%l_calc .and. tscheme%istage == tscheme%nstages  ) then
         call vp_bal%finalize_dvpdt(up_Mloc, tscheme)
      end if

      !-- Finish calculation of d\omega/dt if requested
      if ( vort_bal%l_calc .and. tscheme%istage == tscheme%nstages  ) then
         call vort_bal%finalize_domdt(om_Mloc, tscheme)
      end if

   end subroutine update_om_coll_smat
!------------------------------------------------------------------------------
   subroutine finish_exp_psi_coll_smat(us_Mloc, dVsOm_Mloc, buo_Mloc, lf_Mloc, &
              &                        djxB_Mloc, dpsi_exp_last, vort_bal)

      !-- Input variables
      complex(cp), intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: dVsOm_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: buo_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: lf_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: djxB_Mloc(nMstart:nMstop,n_r_max)

      !-- Output variable
      complex(cp),         intent(inout) :: dpsi_exp_last(nMstart:nMstop,n_r_max)
      type(vort_bal_type), intent(inout) :: vort_bal

      !-- Local variables:
      integer :: n_r, n_m, m
      complex(cp) :: dslf_Mloc(nMstart:nMstop,n_r_max)

      !-- Finish calculation of advection
      call get_dr( dVsOm_Mloc, work_Mloc, nMstart, nMstop, n_r_max, &
           &       rscheme, nocopy=.true. )

      !-- Finish calculation of the z-avg Lorentz-force
      if ( l_mag_LF .and. l_leibniz ) then
         call get_dr( djxB_Mloc, dslf_Mloc, nMstart, nMstop, n_r_max, rscheme, nocopy=.true. )
      end if

#ifdef TOTO
      block
         integer :: n_s, file_handle!, n_m
      if( rank == 0  ) then

         open(newunit=file_handle, file='dslf_Mloc.dat', status='new', form='formatted')
         do n_s=1,n_r_max
            write(file_handle, *) 1./or1(n_s), real(dslf_Mloc(:,n_s)), aimag(dslf_Mloc(:,n_s))
         end do
         close(file_handle)

         open(newunit=file_handle, file='dpsi_exp_last.dat', status='new', form='formatted')
         do n_s=1,n_r_max
            if ( n_s==1 .or. n_s==n_r_max ) dslf_Mloc(:,n_s) = zero
            dpsi_exp_last(:,n_s) = or1(n_s)*(dslf_Mloc(:,n_s) + lf_Mloc(:,n_s)) !-- 1/s(ds<sjxBp> - dp<jxBs> + Leibniz)!-buo_Mloc(:,n_s)
            !if ( n_s==1 .or. n_s==n_r_max ) dpsi_exp_last(:,n_s) = zero
            write(file_handle, *) 1./or1(n_s), real(dpsi_exp_last(:,n_s)), aimag(dpsi_exp_last(:,n_s))
         end do
         close(file_handle)
         print*, 'ALL GOOD finish_exp_psi!**'!
         stop
      end if
      end block
#endif

      !-- Finish calculation of the explicit part for current time step
      do n_r=1,n_r_max
         do n_m=nMstart, nMstop
            m = idx2m(n_m)
            if ( m /= 0 ) then
               dpsi_exp_last(n_m,n_r)=     dpsi_exp_last(n_m,n_r)- &
               &                      or1(n_r)*work_Mloc(n_m,n_r)

               !-- Add z-avg Lorentz force on the r.h.s
               if ( l_mag_LF ) then
                  if ( l_leibniz ) then !-- if Leibniz rule of integration for LF
                     !if ( n_r==1 .or. n_r==n_r_max ) dslf_Mloc(n_m,n_r) = zero
                     lf_Mloc(n_m,n_r) = or1(n_r)*(dslf_Mloc(n_m,n_r) + lf_Mloc(n_m,n_r)) !-- 1/s(ds<sjxBp> + -dp<jxBs> + Leibniz)!
                  end if
                  if ( n_r==1 .or. n_r==n_r_max ) lf_Mloc(n_m,n_r) = zero
                  dpsi_exp_last(n_m,n_r)=dpsi_exp_last(n_m,n_r)+ &
                  &                            lf_Mloc(n_m,n_r)
               end if

               !-- If the force balance is requested get the advection here
               if ( vort_bal%l_calc ) then
                  vort_bal%adv(n_m,n_r)=dpsi_exp_last(n_m,n_r)
                  if ( l_mag_LF ) then
                     !vort_bal%lf(n_m,n_r)=zero
                     vort_bal%lf(n_m,n_r)=lf_Mloc(n_m,n_r)
                     vort_bal%adv(n_m,n_r)=vort_bal%adv(n_m,n_r) - lf_Mloc(n_m,n_r)
                  end if
               end if

               !-- If Coriolis is treated explicitly, add it here:
               if ( .not. l_coriolis_imp ) then
                  dpsi_exp_last(n_m,n_r)=       dpsi_exp_last(n_m,n_r) &
                  &                 +CorFac*beta(n_r)*us_Mloc(n_m,n_r)
               end if

               !if ( l_QG_basis ) then ! Missing NL part from the horizontal adv.
               !   dpsi_exp_last(n_m,n_r)=    dpsi_exp_last(n_m,n_r) &
               !   &            +beta(n_r)*or1(n_r)*dVsOm_Mloc(n_m,n_r)
               !   !-- WARNING!:: WHEN l_mag_LF THERE IS s-LORENTZ FORCE in dVsOm!!! WITH WRONG IMPLEM.
               !   !--        --> NOTHING MORE MAG. in dVsOm with RIGHT compute_sder()!!!
               !   !-- BUT WRONG aw!:: \beta \omega_z u_s COMES FROM the V.(u \omega_z) form
               !   !!-- Update June 2022!:: not even using dVsOm_Mloc to store the magnetic part anymore...
               !end if

               !-- If Buoyancy is treated explicitly, add it here:
               if ( .not. l_buo_imp ) then
                  dpsi_exp_last(n_m,n_r)=dpsi_exp_last(n_m,n_r)+buo_Mloc(n_m,n_r)
               end if
            end if
         end do
      end do

   end subroutine finish_exp_psi_coll_smat
!------------------------------------------------------------------------------
   subroutine get_psi_rhs_imp_coll_smat(us_Mloc, up_Mloc, om_Mloc, dom_Mloc,  &
              &                         psi_last, dpsi_imp_Mloc_last, vp_bal, &
              &                         vort_bal, l_calc_lin_rhs)

      !-- Input variables
      complex(cp), intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: om_Mloc(nMstart:nMstop,n_r_max)
      logical,     intent(in) :: l_calc_lin_rhs

      !-- Output variables
      complex(cp),         intent(out) :: dom_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(inout) :: psi_last(nMstart:nMstop,n_r_max)
      complex(cp),         intent(inout) :: dpsi_imp_Mloc_last(nMstart:nMstop,n_r_max)
      type(vp_bal_type),   intent(inout) :: vp_bal
      type(vort_bal_type), intent(inout) :: vort_bal

      !-- Local variables
      real(cp) :: duphi0(n_r_max), d2uphi0(n_r_max), uphi0(n_r_max)
      real(cp) :: dm2
      integer :: n_r, n_m, m, m0

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m == 0 ) then
               psi_last(n_m,n_r)=up_Mloc(n_m,n_r)
            else
               psi_last(n_m,n_r)=om_Mloc(n_m,n_r)
            end if
         end do
      end do

      if ( l_calc_lin_rhs .or. vp_bal%l_calc .or. vort_bal%l_calc ) then

         call get_ddr(om_Mloc, dom_Mloc, work_Mloc, nMstart, nMstop, &
              &       n_r_max, rscheme)

         m0 = m2idx(0)

         if ( l_rank_has_m0 ) then
            do n_r=1,n_r_max
               uphi0(n_r)=real(up_Mloc(m0, n_r),kind=cp)
            end do
            call get_ddr(uphi0, duphi0, d2uphi0, n_r_max, rscheme)
         end if

         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               if ( m == 0 ) then
                  dpsi_imp_Mloc_last(n_m,n_r)=ViscFac*hdif_V(n_m)*(        &
                  &                                      d2uphi0(n_r)+     &
                  &                        or1(n_r)*      duphi0(n_r)-     &
                  &                        or2(n_r)*       uphi0(n_r) ) -  &
                  &              CorFac*ekpump(n_r)*       uphi0(n_r)

                  if ( vp_bal%l_calc ) then
                     vp_bal%visc(n_r)=ViscFac*hdif_V(n_m)*(d2uphi0(n_r)+   &
                     &                or1(n_r)*duphi0(n_r)-or2(n_r)*uphi0(n_r))
                     if ( l_ek_pump ) then
                        vp_bal%pump(n_r)=-CorFac*ekpump(n_r)*uphi0(n_r)
                     end if
                  end if
               else
                  dm2 = real(m,cp)*real(m,cp)
                  dpsi_imp_Mloc_last(n_m,n_r)=ViscFac*hdif_V(n_m)*(         &
                  &                                    work_Mloc(n_m,n_r)   &
                  &                   +or1(n_r)*        dom_Mloc(n_m,n_r)   &
                  &                   -dm2*or2(n_r)*     om_Mloc(n_m,n_r) ) &
                  &             -CorFac*ekpump(n_r)*     om_Mloc(n_m,n_r)   &
                  & +half*CorFac*ekpump(n_r)*beta(n_r)*  up_Mloc(n_m,n_r)   &
                  & +CorFac*( ekpump(n_r)*beta(n_r)*(-ci*real(m,cp)+        &
                  &              5.0_cp*r_cmb*oheight(n_r)) )*              &
                  &                                      us_Mloc(n_m,n_r)

                  if ( l_coriolis_imp ) then
                     dpsi_imp_Mloc_last(n_m,n_r)=dpsi_imp_Mloc_last(n_m,n_r) &
                     &                   + CorFac*beta(n_r)*us_Mloc(n_m,n_r)
                  end if

                  if ( l_QG_basis ) then
                  !-- Ekman Pumping should not be modified (only dependant on BCs)
                     dpsi_imp_Mloc_last(n_m,n_r)=dpsi_imp_Mloc_last(n_m,n_r) + &
                     &                                     CorFac*ekpump(n_r)* &
                     &           beta(n_r)*third*ci*real(m,cp)*us_Mloc(n_m,n_r)
                  end if

                  !-- In case the force balance is requested:
                  if ( vort_bal%l_calc ) then
                     vort_bal%visc(n_m,n_r)=ViscFac*(     work_Mloc(n_m,n_r) &
                     &                     +or1(n_r)*      dom_Mloc(n_m,n_r) &
                     &                 -dm2*or2(n_r)*       om_Mloc(n_m,n_r) )
                     if ( .not. l_non_rot ) then
                        vort_bal%cor(n_m,n_r) =CorFac*beta(n_r)*us_Mloc(n_m,n_r)
                     end if
                     if ( l_ek_pump ) then
                        vort_bal%pump(n_m,n_r)=CorFac*ekpump(n_r)*(             &
                        &                                     -om_Mloc(n_m,n_r) &
                        &                    +half*beta(n_r)*  up_Mloc(n_m,n_r) &
                        & +beta(n_r)*(-ci*real(m,cp)+5.0_cp*r_cmb*oheight(n_r))*&
                        &                                      us_Mloc(n_m,n_r) )
                     end if
                  end if
               end if
            end do
         end do

      end if ! if wimp /= .or. vp_bal%l_calc .or. vort_bal%l_calc

   end subroutine get_psi_rhs_imp_coll_smat
!------------------------------------------------------------------------------
   subroutine get_psiMat(tscheme, m, psiMat, psiPivot, psiMat_fac, time_lu, &
              &          n_lu_calls)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      integer,             intent(in) :: m

      !-- Output variables
      complex(cp), intent(out) :: psiMat(2*n_r_max,2*n_r_max)
      integer,     intent(out) :: psiPivot(2*n_r_max)
      real(cp),    intent(out) :: psiMat_fac(2*n_r_max,2)
      real(cp),    intent(inout) :: time_lu
      integer,     intent(inout) :: n_lu_calls

      !-- Local variables
      integer :: nR_out, nR, nR_psi, nR_out_psi, info, n_m
      real(cp) :: dm2, runStart, runStop

      n_m = m2idx(m)
      dm2 = real(m,cp)*real(m,cp)

      !----- Boundary conditions:
      do nR_out=1,rscheme%n_max

         nR_out_psi = nR_out+n_r_max

         !-- Non-penetation condition
         psiMat(1,nR_out)          =0.0_cp
         psiMat(1,nR_out_psi)      =rscheme%rnorm*rscheme%rMat(1,nR_out)
         psiMat(n_r_max,nR_out)    =0.0_cp
         psiMat(n_r_max,nR_out_psi)=rscheme%rnorm*rscheme%rMat(n_r_max,nR_out)

         if ( ktopv == 1 ) then ! free-slip
            psiMat(n_r_max+1,nR_out)    =0.0_cp
            psiMat(n_r_max+1,nR_out_psi)=rscheme%rnorm*(                &
            &                                 rscheme%d2rMat(1,nR_out)- &
            &                           or1(1)*rscheme%drMat(1,nR_out) )
         else
            psiMat(n_r_max+1,nR_out)    =0.0_cp
            psiMat(n_r_max+1,nR_out_psi)=rscheme%rnorm*rscheme%drMat(1,nR_out)
         end if
         if ( kbotv == 1 ) then
            psiMat(2*n_r_max,nR_out)    =0.0_cp
            psiMat(2*n_r_max,nR_out_psi)=rscheme%rnorm*(                 &
            &                            rscheme%d2rMat(n_r_max,nR_out)- &
            &                or1(n_r_max)*rscheme%drMat(n_r_max,nR_out) )
         else
            psiMat(2*n_r_max,nR_out)    =0.0_cp
            psiMat(2*n_r_max,nR_out_psi)=rscheme%rnorm* &
            &                            rscheme%drMat(n_r_max,nR_out)
         end if
      end do


      if ( rscheme%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme%n_max+1,n_r_max
            nR_out_psi = nR_out+n_r_max
            psiMat(1,nR_out)            =0.0_cp
            psiMat(n_r_max,nR_out)      =0.0_cp
            psiMat(n_r_max+1,nR_out)    =0.0_cp
            psiMat(2*n_r_max,nR_out)    =0.0_cp
            psiMat(1,nR_out_psi)        =0.0_cp
            psiMat(n_r_max,nR_out_psi)  =0.0_cp
            psiMat(n_r_max+1,nR_out_psi)=0.0_cp
            psiMat(2*n_r_max,nR_out_psi)=0.0_cp
         end do
      end if

      !----- Other points:
      do nR_out=1,n_r_max
         nR_out_psi=nR_out+n_r_max
         do nR=2,n_r_max-1
            nR_psi=nR+n_r_max

            psiMat(nR,nR_out)= rscheme%rnorm * (                         &
            &                                  rscheme%rMat(nR,nR_out) - &
            &         tscheme%wimp_lin(1)*(ViscFac*hdif_V(n_m)*(         &
            &                              rscheme%d2rMat(nR,nR_out) +   &
            &            or1(nR)*           rscheme%drMat(nR,nR_out) -   &
            &            dm2*or2(nR)*        rscheme%rMat(nR,nR_out) ) - &
            &  CorFac*ekpump(nR)*            rscheme%rMat(nR,nR_out) ) )

            psiMat(nR,nR_out_psi)=-rscheme%rnorm*tscheme%wimp_lin(1)*(   &
            &-half*CorFac*ekpump(nR)*beta(nR)*rscheme%drMat(nR,nR_out)+  &
            &  CorFac*( -half*ekpump(nR)*beta(nR)*beta(nR)               &
            &   +ekpump(nR)*beta(nR)*or1(nR)*( dm2+                      &
            &              5.0_cp*r_cmb*oheight(nR)*ci*real(m,cp)) )*    &
            &                                  rscheme%rMat(nR,nR_out) ) 

            if ( l_QG_basis ) then
               !-- Ekman Pumping should not be modified (only dependant on BCs)
               !-- Added in psiMat(nR,nR_out) but the relation is between u_s and \psi
               psiMat(nR,nR_out_psi)= psiMat(nR,nR_out_psi) -             &
               &  tscheme%wimp_lin(1)*ViscFac*hdif_V(n_m)*rscheme%rnorm*( &
               &  CorFac*ekpump(nR)*third*dm2*or1(nR)*                    &
               &                                 rscheme%rMat(nR,nR_out) )
            end if

            if ( l_coriolis_imp ) then
               psiMat(nR,nR_out_psi) = psiMat(nR,nR_out_psi) -        &
               &                 rscheme%rnorm*tscheme%wimp_lin(1)*   &
               &               CorFac*beta(nR)*or1(nR)*ci*real(m,cp)* &
               &                rscheme%rMat(nR,nR_out)
            end if

            psiMat(nR_psi,nR_out)= rscheme%rnorm*rscheme%rMat(nR,nR_out)

            psiMat(nR_psi,nR_out_psi)= rscheme%rnorm * (              &
            &                             rscheme%d2rMat(nR,nR_out) + &
            &      (or1(nR)+beta(nR))*     rscheme%drMat(nR,nR_out) + &
            &  (or1(nR)*beta(nR)+dbeta(nR)-dm2*or2(nR))*              &
            &                               rscheme%rMat(nR,nR_out) )

            if ( l_QG_basis ) then ! Relation \omega_z / \psi above
               !-- Additional term of the projection is \beta/3 m^2/s \psi
               psiMat(nR_psi,nR_out_psi)= psiMat(nR_psi,nR_out_psi) +    &
               &                          rscheme%rnorm * (              &
               &   +beta(nR)*third*dm2*or1(nR)*                          &
               &                               rscheme%rMat(nR,nR_out) )
            end if

         end do
      end do

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         nR_psi = nR+n_r_max
         psiMat(nR,1)            =rscheme%boundary_fac*psiMat(nR,1)
         psiMat(nR,n_r_max)      =rscheme%boundary_fac*psiMat(nR,n_r_max)
         psiMat(nR,n_r_max+1)    =rscheme%boundary_fac*psiMat(nR,n_r_max+1)
         psiMat(nR,2*n_r_max)    =rscheme%boundary_fac*psiMat(nR,2*n_r_max)
         psiMat(nR_psi,1)        =rscheme%boundary_fac*psiMat(nR_psi,1)
         psiMat(nR_psi,n_r_max)  =rscheme%boundary_fac*psiMat(nR_psi,n_r_max)
         psiMat(nR_psi,n_r_max+1)=rscheme%boundary_fac*psiMat(nR_psi,n_r_max+1)
         psiMat(nR_psi,2*n_r_max)=rscheme%boundary_fac*psiMat(nR_psi,2*n_r_max)
      end do

      ! compute the linesum of each line
      do nR=1,2*n_r_max
         psiMat_fac(nR,1)=one/maxval(abs(psiMat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nR=1,2*n_r_max
         psiMat(nR,:) = psiMat(nR,:)*psiMat_fac(nR,1)
      end do

      ! also compute the rowsum of each column
      do nR=1,2*n_r_max
         psiMat_fac(nR,2)=one/maxval(abs(psiMat(:,nR)))
      end do
      ! now divide each row by the rowsum
      do nR=1,2*n_r_max
         psiMat(:,nR) = psiMat(:,nR)*psiMat_fac(nR,2)
      end do

      !----- LU decomposition:
      runStart = MPI_Wtime()
      call prepare_full_mat(psiMat,2*n_r_max,2*n_r_max,psiPivot,info)
      runStop = MPI_Wtime()
      if ( runStop > runStart ) then
         time_lu = time_lu+(runStop-runStart)
         n_lu_calls = n_lu_calls+1
      end if
      if ( info /= 0 ) then
         call abortRun('Singular matrix psiMat!')
      end if


   end subroutine get_psiMat
!------------------------------------------------------------------------------
   subroutine get_uphiMat(tscheme, uphiMat, uphiPivot)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      real(cp), intent(out) :: uphiMat(n_r_max,n_r_max)
      integer,  intent(out) :: uphiPivot(n_r_max)

      !-- Local variables
      integer :: nR_out, nR, info, n_m

      n_m = m2idx(0)

      !----- Boundary conditions:
      do nR_out=1,rscheme%n_max
         if ( ktopv == 1 ) then !-- Free-slip
            uphiMat(1,nR_out)=rscheme%rnorm*(rscheme%drMat(1,nR_out)-or1(1)* &
            &                                 rscheme%rMat(1,nR_out))
         else
            uphiMat(1,nR_out)=rscheme%rnorm*rscheme%rMat(1,nR_out)
         end if
         if ( kbotv == 1 ) then !-- Free-slip
            uphiMat(n_r_max,nR_out)=rscheme%rnorm*(                 &
            &                        rscheme%drMat(n_r_max,nR_out)  &
            &           -or1(n_r_max)*rscheme%rMat(n_r_max,nR_out))
         else
            uphiMat(n_r_max,nR_out)=rscheme%rnorm* &
            &                       rscheme%rMat(n_r_max,nR_out)
         end if
      end do


      if ( rscheme%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme%n_max+1,n_r_max
            uphiMat(1,nR_out)      =0.0_cp
            uphiMat(n_r_max,nR_out)=0.0_cp
         end do
      end if

      !----- Other points:
      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            uphiMat(nR,nR_out)= rscheme%rnorm * (                     &
            &                               rscheme%rMat(nR,nR_out) - &
            &tscheme%wimp_lin(1)*(ViscFac*hdif_V(n_m)*(               &
            &                           rscheme%d2rMat(nR,nR_out) +   &
            &            or1(nR)*        rscheme%drMat(nR,nR_out) -   &
            &            or2(nR)*         rscheme%rMat(nR,nR_out) ) - &
            &  CorFac*ekpump(nR)*         rscheme%rMat(nR,nR_out) ) )
         end do
      end do

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         uphiMat(nR,1)      =rscheme%boundary_fac*uphiMat(nR,1)
         uphiMat(nR,n_r_max)=rscheme%boundary_fac*uphiMat(nR,n_r_max)
      end do

      !----- LU decomposition:
      call prepare_full_mat(uphiMat,n_r_max,n_r_max,uphiPivot,info)
      if ( info /= 0 ) then
         call abortRun('Singular matrix uphiMat!')
      end if

   end subroutine get_uphiMat
!------------------------------------------------------------------------------
end module update_psi_coll_smat
