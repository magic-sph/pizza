module update_psi_coll_dmat

   use precision_mod
   use parallel_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: one, zero, ci, half
   use horizontal, only: hdif_V
   use namelists, only: kbotv, ktopv, alpha, r_cmb, CorFac, ViscFac, &
       &                l_non_rot, l_ek_pump, l_buo_imp
   use radial_functions, only: rscheme, or1, or2, beta, dbeta, ekpump, oheight
   use blocking, only: nMstart, nMstop, l_rank_has_m0
   use truncation, only: n_r_max, idx2m, m2idx
   use radial_der, only: get_ddr, get_dr
   use fields, only: work_Mloc
   use algebra, only: prepare_full_mat, solve_full_mat
   use time_schemes, only: type_tscheme
   use vort_balance, only: vort_bal_type
   use vp_balance, only: vp_bal_type
   use timers_mod, only: timers_type
   use time_array
   use useful, only: abortRun

   implicit none
   
   private

   logical,  allocatable :: lPsimat(:)
   real(cp), allocatable :: psiMat(:,:,:), omMat(:,:,:)
   real(cp), allocatable :: uphiMat(:,:)
   real(cp), allocatable :: om2_Mloc(:,:), om3_Mloc(:,:)
   real(cp), allocatable :: psi2_Mloc(:,:), psi3_Mloc(:,:)
   real(cp), allocatable :: influence_matrix_Mloc(:,:,:)
   integer,  allocatable :: psiPivot(:,:), omPivot(:,:)
   real(cp), allocatable :: psiMat_fac(:,:,:), omMat_fac(:,:,:)
   complex(cp), allocatable :: rhs(:)
   real(cp), allocatable :: rhs_m0(:)

   interface  get_bc_influence_matrix
      module procedure get_bc_influence_matrix_complex
      module procedure get_bc_influence_matrix_real
   end interface get_bc_influence_matrix

   public :: update_om_coll_dmat, initialize_om_coll_dmat, finalize_om_coll_dmat, &
   &         get_psi_rhs_imp_coll_dmat, finish_exp_psi_coll_dmat

contains

   subroutine initialize_om_coll_dmat

      !-- Local variables
      integer :: n_m, m

      allocate( lPsimat(nMstart:nMstop) )
      lPsimat(:)=.false.
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_LOGICAL

      allocate( psiMat(n_r_max, n_r_max, nMstart:nMstop) )
      allocate( psiPivot(n_r_max, nMstart:nMstop) )
      allocate( psiMat_fac(n_r_max, 2, nMstart:nMstop) )
      allocate( rhs(n_r_max), rhs_m0(n_r_max) )
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*n_r_max*n_r_max* &
      &                 SIZEOF_DEF_REAL+n_r_max*(nMstop-nMstart+1)*         &
      &                 SIZEOF_INTEGER+ n_r_max*(2+2*(nMstop-nMstart+1))*   &
      &                 SIZEOF_DEF_REAL

      allocate( omMat(n_r_max, n_r_max, nMstart:nMstop) )
      allocate( omPivot(n_r_max, nMstart:nMstop) )
      allocate( omMat_fac(n_r_max, 2, nMstart:nMstop) )
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*n_r_max*n_r_max* &
      &                 SIZEOF_DEF_REAL+n_r_max*(nMstop-nMstart+1)*         & 
      &                 SIZEOF_INTEGER+2*n_r_max*(nMstop-nMstart+1)*        &
      &                 SIZEOF_DEF_REAL

      allocate( uphiMat(n_r_max,n_r_max) )
      bytes_allocated = bytes_allocated+n_r_max*n_r_max*SIZEOF_DEF_REAL

      allocate( om2_Mloc(nMstart:nMstop,n_r_max),om3_Mloc(nMstart:nMstop,n_r_max) )
      allocate( psi2_Mloc(nMstart:nMstop,n_r_max),psi3_Mloc(nMstart:nMstop,n_r_max) )
      om2_Mloc(:,:) =0.0_cp
      om3_Mloc(:,:) =0.0_cp
      psi2_Mloc(:,:)=0.0_cp
      psi3_Mloc(:,:)=0.0_cp
      bytes_allocated = bytes_allocated+4*(nMstop-nMstart+1)*n_r_max*SIZEOF_DEF_REAL

      allocate( influence_matrix_Mloc(2,2,nMstart:nMstop) )
      influence_matrix_Mloc(:,:,:)=0.0_cp
      bytes_allocated = bytes_allocated+4*(nMstop-nMstart+1)*SIZEOF_DEF_REAL

      do n_m=nMstart,nMstop
         m = idx2m(n_m)
         if ( m /= 0 ) call get_psiMat(m, psiMat(:,:,n_m), psiPivot(:,n_m), &
                            &          psiMat_fac(:,:,n_m))
      end do

   end subroutine initialize_om_coll_dmat
!------------------------------------------------------------------------------
   subroutine finalize_om_coll_dmat

      deallocate( influence_matrix_Mloc )
      deallocate( psi2_Mloc, psi3_Mloc, om2_Mloc, om3_Mloc ) 
      deallocate( rhs_m0, rhs, psiMat_fac )
      deallocate( lPsimat, psiMat, uphiMat, psiPivot )
      deallocate( omMat_fac, omMat, omPivot )

   end subroutine finalize_om_coll_dmat
!------------------------------------------------------------------------------
   subroutine update_om_coll_dmat(psi_Mloc, om_Mloc, dom_Mloc, us_Mloc, up_Mloc, &
              &                   buo_Mloc, dpsidt, vp_bal, vort_bal, tscheme,   &
              &                   lMat, timers)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lMat
      complex(cp),         intent(in) :: buo_Mloc(nMstart:nMstop,n_r_max)

      !-- Output variables
      complex(cp),         intent(out) :: psi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(inout) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: dom_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(inout) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(inout) :: up_Mloc(nMstart:nMstop,n_r_max)
      type(vp_bal_type),   intent(inout) :: vp_bal
      type(vort_bal_type), intent(inout) :: vort_bal
      type(type_tarray),   intent(inout) :: dpsidt
      type(timers_type),   intent(inout) :: timers

      !-- Local variables
      real(cp) :: uphi0(n_r_max), om0(n_r_max), runStart, runStop
      integer :: n_r, n_m, n_cheb, m

      if ( lMat ) lPsimat(:)=.false.

      !-- Calculation of the implicit part
      call get_psi_rhs_imp_coll_dmat(up_Mloc, om_Mloc, dom_Mloc,           &
           &                         dpsidt%old(:,:,tscheme%istage),       &
           &                         dpsidt%impl(:,:,tscheme%istage),      &
           &                         vp_bal, vort_bal,                     &
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
               end do
            end if

            call solve_full_mat(uphiMat(:,:), n_r_max, n_r_max,   &
                 &              psiPivot(1:n_r_max,n_m), rhs_m0(:))

            do n_cheb=1,rscheme%n_max
               uphi0(n_cheb)=rhs_m0(n_cheb)
            end do

         else ! Non-axisymmetric components
         
            if ( .not. lPsimat(n_m) ) then
               call get_omMat(tscheme, m, omMat(:,:,n_m), omPivot(:,n_m), &
                    &          omMat_fac(:,:,n_m), timers%lu, timers%n_lu_calls)

               !-- Now solve the first intermediate problem
               !-- L_\omega \omega = 0; \omega(ro)=0; \omega(ri)=1
               call get_intermediate_om(omMat(:,:,n_m), m, omPivot(:,n_m), &
                    &                   omMat_fac(:,:,n_m), om2_Mloc, n_bc_set=1)
               !-- Now solve the second intermediate problem
               !-- L_\omega \omega = 0; \omega(ro)=1; \omega(ri)=0
               call get_intermediate_om(omMat(:,:,n_m), m, omPivot(:,n_m),  &
                    &                   omMat_fac(:,:,n_m), om3_Mloc, n_bc_set=2)

               lPsimat(n_m)=.true.
            end if

            if ( vort_bal%l_calc ) then
               do n_r=1,n_r_max
                  vort_bal%buo(n_m,n_r)=buo_Mloc(n_m,n_r)/tscheme%dt(1)
               end do
            end if

            rhs(1)        =zero
            rhs(n_r_max)  =zero
            do n_r=2,n_r_max-1
               !-- Add buoyancy
               if ( l_buo_imp ) then
                  rhs(n_r)=work_Mloc(n_m,n_r)+buo_Mloc(n_m,n_r)
               else
                  rhs(n_r)=work_Mloc(n_m,n_r)
               end if
            end do

            do n_r=1,n_r_max
               rhs(n_r) = rhs(n_r)*omMat_fac(n_r,1,n_m)
            end do
            runStart = MPI_Wtime()
            call solve_full_mat(omMat(:,:,n_m), n_r_max, n_r_max, &
                 &              omPivot(:, n_m), rhs(:))
            runStop = MPI_Wtime()
            if ( runStop > runStart ) then
               timers%solve = timers%solve + (runStop-runStart)
               timers%n_solve_calls = timers%n_solve_calls+1
            end if
            do n_r=1,n_r_max
               rhs(n_r) = rhs(n_r)*omMat_fac(n_r,2,n_m)
            end do

            do n_cheb=1,rscheme%n_max
               om_Mloc(n_m,n_cheb) =rhs(n_cheb)
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
                  if ( lMat ) then
                     om2_Mloc(n_m,n_cheb)=0.0_cp
                     om3_Mloc(n_m,n_cheb)=0.0_cp
                  end if
               end if
            end do
         end do
      end if

      !-- Bring uphi0 to the physical space
      if ( l_rank_has_m0 ) then
         call get_dr(uphi0, om0, n_r_max, rscheme, l_dct=.false.)
         call rscheme%costf1(uphi0, n_r_max)
      end if

      !-- Copy omega into work (it allows to keep omega in Cheb space)
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m /= 0 ) work_Mloc(n_m,n_r)=om_Mloc(n_m,n_r)
         end do
      end do

      !-- Bring  work_Mloc to the physical space
      !-- It will be used in the RHS of the psi equation
      runStart = MPI_Wtime()
      call rscheme%costf1(work_Mloc, nMstart, nMstop, n_r_max)
      runStop = MPI_Wtime()
      if ( runStop > runStart ) then
         timers%dct = timers%dct + (runStop-runStart)
         timers%n_dct_calls = timers%n_dct_calls + 1
      end if


      !-- Now that omega has been determined let's get psi
      do n_m=nMstart,nMstop
         m = idx2m(n_m)

         if ( m /= 0 ) then

            if ( lMat ) then
               !-- L_\psi \psi = -\omega
               do n_r=1,n_r_max
                  rhs_m0(n_r)=om2_Mloc(n_m,n_r)
               end do
               call get_intermediate_psi(psiMat(:,:,n_m), m, psiPivot(:,n_m),  &
                    &                    psiMat_fac(:,:,n_m), rhs_m0, psi2_Mloc)
               do n_r=1,n_r_max
                  rhs_m0(n_r)=om3_Mloc(n_m,n_r)
               end do
               call get_intermediate_psi(psiMat(:,:,n_m), m, psiPivot(:,n_m),  &
                    &                    psiMat_fac(:,:,n_m), rhs_m0, psi3_Mloc)

               lPsimat(n_m)=.true.
            end if

            rhs(1)      =zero
            rhs(n_r_max)=zero
            do n_r=2,n_r_max-1
               !-- L \psi = -\omega
               rhs(n_r)=-work_Mloc(n_m,n_r)
            end do

            do n_r=1,n_r_max
               rhs(n_r) = rhs(n_r)*psiMat_fac(n_r,1,n_m)
            end do
            runStart = MPI_Wtime()
            call solve_full_mat(psiMat(:,:,n_m), n_r_max, n_r_max, &
                 &              psiPivot(:, n_m), rhs(:))
            runStop = MPI_Wtime()
            if ( runStop > runStart ) then
               timers%solve = timers%solve + (runStop-runStart)
               timers%n_solve_calls = timers%n_solve_calls+1
            end if
            do n_r=1,n_r_max
               rhs(n_r) = rhs(n_r)*psiMat_fac(n_r,2,n_m)
            end do

            do n_cheb=1,rscheme%n_max
               psi_Mloc(n_m,n_cheb) =rhs(n_cheb)
            end do

         end if

      end do

      !-- set cheb modes > rscheme%n_max to zero (dealiazing)
      if ( rscheme%n_max < n_r_max ) then ! fill with zeros !
         do n_cheb=rscheme%n_max+1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               if ( m /= 0 ) then
                  psi_Mloc(n_m,n_cheb) =zero
                  if ( lMat ) then
                     psi2_Mloc(n_m,n_cheb)=0.0_cp
                     psi3_Mloc(n_m,n_cheb)=0.0_cp
                  end if
               end if
            end do
         end do
      end if


      !-- If the matrix has been rebuilt one also needs to assemble
      !-- again the influence matrix
      if ( lMat ) then
         call get_influence_matrix(psi2_Mloc, psi3_Mloc, influence_matrix_Mloc)
      end if

      !-- Finally assemble psi and omega from intermediate solutions:
      !-- \omega = \omega + c1*\omega_2 + c2*\omega_3
      !-- \psi = \psi + c1*\psi_2 + c2*\psi_3
      call solve_influence_matrix(psi_Mloc, psi2_Mloc, psi3_Mloc, om_Mloc, &
           &                      om2_Mloc, om3_Mloc, influence_matrix_Mloc)

      !-- Get the radial derivative of psi to calculate uphi
      call get_dr(psi_Mloc, work_Mloc, nMstart, nMstop, n_r_max, rscheme, &
           &      l_dct_in=.false.)

      !-- Finally bring psi and omega to the physical space
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
               up_Mloc(n_m,n_r)=uphi0(n_r)
               om_Mloc(n_m,n_r)=om0(n_r)+or1(n_r)*uphi0(n_r)
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

   end subroutine update_om_coll_dmat
!------------------------------------------------------------------------------
   subroutine finish_exp_psi_coll_dmat(us_Mloc, up_Mloc, om_Mloc, dVsOm_Mloc, &
              &                        buo_Mloc, dpsi_exp_last, vort_bal)

      !-- Input variables
      complex(cp), intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: dVsOm_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: buo_Mloc(nMstart:nMstop,n_r_max)

      !-- Output variable
      complex(cp),         intent(inout) :: dpsi_exp_last(nMstart:nMstop,n_r_max)
      type(vort_bal_type), intent(inout) :: vort_bal


      !-- Local variables:
      integer :: n_r, n_m, m
   
      !-- Finish calculation of advection
      call get_dr( dVsOm_Mloc, work_Mloc, nMstart, nMstop, n_r_max, &
           &       rscheme, nocopy=.true.)

      !-- Finish calculation of the explicit part for current time step
      do n_r=1,n_r_max
         do n_m=nMstart, nMstop
            m = idx2m(n_m)
            if ( m /= 0 ) then
               dpsi_exp_last(n_m,n_r)=dpsi_exp_last(n_m,n_r)-&
               &                        or1(n_r)*work_Mloc(n_m,n_r)

               !-- If the force balance is requested get the advection here
               if ( vort_bal%l_calc ) then
                  vort_bal%adv(n_m,n_r)=dpsi_exp_last(n_m,n_r)
                  if ( .not. l_non_rot ) then
                     vort_bal%cor(n_m,n_r) =CorFac*beta(n_r)*us_Mloc(n_m,n_r)
                  end if
               end if

               !-- If this is rotating, treat it explicitly in this case
               if ( .not. l_non_rot ) then
                  dpsi_exp_last(n_m,n_r)=            dpsi_exp_last(n_m,n_r) &
                  &                      +CorFac*beta(n_r)*us_Mloc(n_m,n_r)
               end if

               !-- If Buoyancy is treated explicitly, add it here:
               if ( .not. l_buo_imp ) then
                  dpsi_exp_last(n_m,n_r)=dpsi_exp_last(n_m,n_r)+buo_Mloc(n_m,n_r)
               end if
            end if
         end do
      end do

      if ( l_ek_pump ) then
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               if ( m /= 0 ) then
                  dpsi_exp_last(n_m,n_r)=       dpsi_exp_last(n_m,n_r)  - &
                  &                CorFac*ekpump(n_r)*om_Mloc(n_m,n_r)  + &
                  & half*CorFac*ekpump(n_r)*beta(n_r)*up_Mloc(n_m,n_r)  + &
                  &      CorFac*( ekpump(n_r)*beta(n_r)*(-ci*real(m,cp) + &
                  &    5.0_cp*r_cmb*oheight(n_r)) ) * us_Mloc(n_m,n_r)

                  !-- In case the force balance is requested:
                  if ( vort_bal%l_calc ) then
                     vort_bal%pump(n_m,n_r)=CorFac*ekpump(n_r)*(             &
                     &                                     -om_Mloc(n_m,n_r) &
                     &                    +half*beta(n_r)*  up_Mloc(n_m,n_r) &
                     & +beta(n_r)*(-ci*real(m,cp)+5.0_cp*r_cmb*oheight(n_r))*&
                     &                                      us_Mloc(n_m,n_r) )
                  end if
               end if
            end do
         end do
      end if

   end subroutine finish_exp_psi_coll_dmat
!------------------------------------------------------------------------------
   subroutine get_psi_rhs_imp_coll_dmat(up_Mloc, om_Mloc, dom_Mloc, psi_last, &
              &                         dpsi_imp_Mloc_last, vp_bal, vort_bal, &
              &                         l_calc_lin_rhs)

      !-- Input variables
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
                  dpsi_imp_Mloc_last(n_m,n_r)=ViscFac*hdif_V(n_m)*(      &
                  &                                      d2uphi0(n_r)+   &
                  &                        or1(n_r)*      duphi0(n_r)-   &
                  &                        or2(n_r)*       uphi0(n_r) )- &
                  &              CorFac*ekpump(n_r)*       uphi0(n_r)

                  if ( vp_bal%l_calc ) then
                     vp_bal%visc(n_r)=ViscFac*hdif_V(n_m)*(d2uphi0(n_r)+ &
                     &                or1(n_r)*duphi0(n_r)-or2(n_r)*uphi0(n_r))
                     if ( l_ek_pump ) then
                        vp_bal%pump(n_r)=-CorFac*ekpump(n_r)*uphi0(n_r)
                     end if
                  end if
               else
                  dm2 = real(m,cp)*real(m,cp)
                  dpsi_imp_Mloc_last(n_m,n_r)=ViscFac*hdif_V(n_m)*(       &
                  &                                    work_Mloc(n_m,n_r) &
                  &                   +or1(n_r)*        dom_Mloc(n_m,n_r) &
                  &               -dm2*or2(n_r)*         om_Mloc(n_m,n_r) )

                  !-- In case the force balance is requested:
                  if ( vort_bal%l_calc ) then
                     vort_bal%visc(n_m,n_r)=              work_Mloc(n_m,n_r) &
                     &                     +or1(n_r)*      dom_Mloc(n_m,n_r) &
                     &                 -dm2*or2(n_r)*       om_Mloc(n_m,n_r)
                  end if

               end if
            end do
         end do

      end if ! if wimp /= .or. vp_bal%l_calc .or. vort_bal%l_calc

   end subroutine get_psi_rhs_imp_coll_dmat
!------------------------------------------------------------------------------
   subroutine get_omMat(tscheme, m, omMat, omPivot, omMat_fac, time_lu, &
              &          n_lu_calls)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      integer,             intent(in) :: m

      !-- Output variables
      real(cp), intent(out) :: omMat(n_r_max,n_r_max)
      integer,  intent(out) :: omPivot(n_r_max)
      real(cp), intent(out) :: omMat_fac(n_r_max,2)
      real(cp), intent(inout) :: time_lu
      integer,  intent(inout) :: n_lu_calls

      !-- Local variables
      integer :: n_r, n_cheb, info, n_m
      real(cp) :: dm2, runStart, runStop

      dm2 = real(m,cp)*real(m,cp)
      n_m = m2idx(m)

      !----- Boundary conditions:
      do n_cheb=1,rscheme%n_max
         omMat(1,n_cheb)      =rscheme%rnorm*rscheme%rMat(1,n_cheb)
         omMat(n_r_max,n_cheb)=rscheme%rnorm*rscheme%rMat(n_r_max,n_cheb)
      end do

      if ( rscheme%n_max < n_r_max ) then ! fill with zeros !
         do n_cheb=rscheme%n_max+1,n_r_max
            omMat(1,n_cheb)      =0.0_cp
            omMat(n_r_max,n_cheb)=0.0_cp
         end do
      end if

      !----- Other points:
      do n_cheb=1,n_r_max
         do n_r=2,n_r_max-1

            omMat(n_r,n_cheb)= rscheme%rnorm * (                          &
            &                                  rscheme%rMat(n_r,n_cheb) - &
            &              tscheme%wimp_lin(1)*ViscFac*hdif_V(n_m)*(      &
            &                                rscheme%d2rMat(n_r,n_cheb) + &
            &            or1(n_r)*            rscheme%drMat(n_r,n_cheb) - &
            &            dm2*or2(n_r)*         rscheme%rMat(n_r,n_cheb) ) )

         end do
      end do

      !----- Factor for highest and lowest cheb:
      do n_r=1,n_r_max
         omMat(n_r,1)      =rscheme%boundary_fac*omMat(n_r,1)
         omMat(n_r,n_r_max)=rscheme%boundary_fac*omMat(n_r,n_r_max)
      end do

      ! compute the linesum of each line
      do n_r=1,n_r_max
         omMat_fac(n_r,1)=one/maxval(abs(omMat(n_r,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do n_r=1,n_r_max
         omMat(n_r,:) = omMat(n_r,:)*omMat_fac(n_r,1)
      end do

      ! also compute the rowsum of each column
      do n_r=1,n_r_max
         omMat_fac(n_r,2)=one/maxval(abs(omMat(:,n_r)))
      end do
      ! now divide each row by the rowsum
      do n_r=1,n_r_max
         omMat(:,n_r) = omMat(:,n_r)*omMat_fac(n_r,2)
      end do

      !----- LU decomposition:
      runStart = MPI_Wtime()
      call prepare_full_mat(omMat,n_r_max,n_r_max,omPivot,info)
      runStop = MPI_Wtime()
      if ( runStop > runStart ) then
         time_lu = time_lu+(runStop-runStart)
         n_lu_calls = n_lu_calls+1
      end if
      if ( info /= 0 ) then
         call abortRun('Singular matrix omMat!')
      end if

   end subroutine get_omMat
!------------------------------------------------------------------------------
   subroutine get_psiMat(m, psiMat, psiPivot, psiMat_fac)

      !-- Input variables
      integer, intent(in) :: m

      !-- Output variables
      real(cp), intent(out) :: psiMat(n_r_max,n_r_max)
      integer,  intent(out) :: psiPivot(n_r_max)
      real(cp), intent(out) :: psiMat_fac(n_r_max,2)

      !-- Local variables
      integer :: n_r, n_cheb, info
      real(cp) :: dm2

      dm2 = real(m,cp)*real(m,cp)

      !----- Boundary conditions:
      do n_cheb=1,rscheme%n_max
         psiMat(1,n_cheb)      =rscheme%rnorm*rscheme%rMat(1,n_cheb)
         psiMat(n_r_max,n_cheb)=rscheme%rnorm*rscheme%rMat(n_r_max,n_cheb)
      end do

      if ( rscheme%n_max < n_r_max ) then ! fill with zeros !
         do n_cheb=rscheme%n_max+1,n_r_max
            psiMat(1,n_cheb)      =0.0_cp
            psiMat(n_r_max,n_cheb)=0.0_cp
         end do
      end if

      !----- Other points:
      do n_cheb=1,n_r_max
         do n_r=2,n_r_max-1

            psiMat(n_r,n_cheb)= rscheme%rnorm * (                        &
            &                               rscheme%d2rMat(n_r,n_cheb) + &
            &      (or1(n_r)+beta(n_r))*     rscheme%drMat(n_r,n_cheb) + &
            &  (or1(n_r)*beta(n_r)+dbeta(n_r)-dm2*or2(n_r))*             &
            &                                 rscheme%rMat(n_r,n_cheb) )

         end do
      end do

      !----- Factor for highest and lowest cheb:
      do n_r=1,n_r_max
         psiMat(n_r,1)      =rscheme%boundary_fac*psiMat(n_r,1)
         psiMat(n_r,n_r_max)=rscheme%boundary_fac*psiMat(n_r,n_r_max)
      end do

      ! compute the linesum of each line
      do n_r=1,n_r_max
         psiMat_fac(n_r,1)=one/maxval(abs(psiMat(n_r,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do n_r=1,n_r_max
         psiMat(n_r,:) = psiMat(n_r,:)*psiMat_fac(n_r,1)
      end do

      ! also compute the rowsum of each column
      do n_r=1,n_r_max
         psiMat_fac(n_r,2)=one/maxval(abs(psiMat(:,n_r)))
      end do
      ! now divide each row by the rowsum
      do n_r=1,n_r_max
         psiMat(:,n_r) = psiMat(:,n_r)*psiMat_fac(n_r,2)
      end do

      !----- LU decomposition:
      call prepare_full_mat(psiMat,n_r_max,n_r_max,psiPivot,info)
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
            &        tscheme%wimp_lin(1)*(ViscFac*hdif_V(n_m)*(       &
            &                           rscheme%d2rMat(nR,nR_out) +   &
            &            or1(nR)*        rscheme%drMat(nR,nR_out) -   &
            &            or2(nR)*         rscheme%rMat(nR,nR_out) ) - &
            &          CorFac*ekpump(nR)* rscheme%rMat(nR,nR_out) ) )
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
   subroutine get_intermediate_om(om_mat, m, omPivot, omMat_fac, om_tmp, n_bc_set)

      !-- Input variables
      integer,  intent(in) :: m
      real(cp), intent(in) :: om_mat(n_r_max,n_r_max)
      integer,  intent(in) :: omPivot(n_r_max)
      real(cp), intent(in) :: omMat_fac(n_r_max,2)
      integer,  intent(in) :: n_bc_set

      !-- Output variables
      real(cp), intent(inout) :: om_tmp(nMstart:nMstop, n_r_max)

      !-- Local variables
      real(cp) :: rhs(n_r_max)
      integer :: n_r, n_m, n_cheb

      n_m = m2idx(m)

      do n_r=2,n_r_max-1
         rhs(n_r)=0.0_cp
      end do

      !-- Apply BCs on the 2 first points
      if ( n_bc_set == 1 ) then
         rhs(1)      =0.0_cp
         rhs(n_r_max)=omMat_fac(n_r_max,1)
      else if ( n_bc_set == 2 ) then
         rhs(1)      =omMat_fac(1,1)
         rhs(n_r_max)=0.0_cp
      end if

      call solve_full_mat(om_mat, n_r_max, n_r_max, omPivot, rhs)

      do n_r=1,n_r_max
         rhs(n_r)=rhs(n_r)*omMat_fac(n_r,2)
      end do

      do n_cheb=1,rscheme%n_max
         om_tmp(n_m,n_cheb) =rhs(n_cheb)
      end do

   end subroutine get_intermediate_om
!------------------------------------------------------------------------------
   subroutine get_intermediate_psi(psi_mat, m, psiPivot, psiMat_fac, om_tmp, &
              &                    psi_tmp)

      !-- Input variables
      integer,  intent(in) :: m
      real(cp), intent(in) :: psi_mat(n_r_max,n_r_max)
      integer,  intent(in) :: psiPivot(n_r_max)
      real(cp), intent(in) :: psiMat_fac(n_r_max,2)

      !-- Output variables
      real(cp), intent(inout) :: om_tmp(n_r_max)
      real(cp), intent(inout) :: psi_tmp(nMstart:nMstop, n_r_max)

      !-- Local variables
      real(cp) :: rhs(n_r_max)
      integer :: n_r, n_m, n_cheb

      n_m = m2idx(m)
      call rscheme%costf1(om_tmp,n_r_max)

      rhs(1)      =0.0_cp
      rhs(n_r_max)=0.0_cp
      do n_r=2,n_r_max-1
         rhs(n_r)=-om_tmp(n_r)
      end do

      do n_r=1,n_r_max
         rhs(n_r)=rhs(n_r)*psiMat_fac(n_r,1)
      end do

      call solve_full_mat(psi_mat, n_r_max, n_r_max, psiPivot, rhs)
      do n_r=1,n_r_max
         rhs(n_r)=rhs(n_r)*psiMat_fac(n_r,2)
      end do

      do n_cheb=1,rscheme%n_max
         psi_tmp(n_m,n_cheb)=rhs(n_cheb)
      end do

   end subroutine get_intermediate_psi
!------------------------------------------------------------------------------
   subroutine get_influence_matrix(psi2_Mloc, psi3_Mloc, influence_mat)

      !-- Input variables
      real(cp), intent(inout) :: psi2_Mloc(nMstart:nMstop,n_r_max)
      real(cp), intent(inout) :: psi3_Mloc(nMstart:nMstop,n_r_max)

      !-- Output variables
      real(cp), intent(out) :: influence_mat(2,2,nMstart:nMstop)

      !-- Local variables
      real(cp) :: bcTop2(nMstart:nMstop), bcBot2(nMstart:nMstop)
      real(cp) :: bcTop3(nMstart:nMstop), bcBot3(nMstart:nMstop)
      real(cp) :: det
      integer :: n_m, m

      call get_bc_influence_matrix(psi2_Mloc, bcBot2, bcTop2)
      call get_bc_influence_matrix(psi3_Mloc, bcBot3, bcTop3)

      do n_m=nMstart, nMstop

         m = idx2m(n_m)

         if ( m /= 0 ) then
            det = bcTop2(n_m)*bcBot3(n_m)-bcTop3(n_m)*bcBot2(n_m)

            influence_mat(1,1,n_m)= bcBot3(n_m)/det
            influence_mat(1,2,n_m)=-bcTop3(n_m)/det
            influence_mat(2,1,n_m)=-bcBot2(n_m)/det
            influence_mat(2,2,n_m)= bcTop2(n_m)/det

         end if

      end do

   end subroutine get_influence_matrix
!------------------------------------------------------------------------------
   subroutine solve_influence_matrix(psi1_Mloc, psi2_Mloc, psi3_Mloc, &
              &                      om1_Mloc, om2_Mloc, om3_Mloc,    &
              &                      influence_mat)

      !-- Input variables
      real(cp), intent(inout) :: psi2_Mloc(nMstart:nMstop,n_r_max)
      real(cp), intent(inout) :: psi3_Mloc(nMstart:nMstop,n_r_max)
      real(cp), intent(in) :: om2_Mloc(nMstart:nMstop,n_r_max)
      real(cp), intent(in) :: om3_Mloc(nMstart:nMstop,n_r_max)
      real(cp), intent(in) :: influence_mat(2,2,nMstart:nMstop)

      !-- Input/Output variables
      complex(cp), intent(inout) :: psi1_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: om1_Mloc(nMstart:nMstop,n_r_max)

      !-- Local variables
      complex(cp) :: bcTop1(nMstart:nMstop), bcBot1(nMstart:nMstop)
      complex(cp) :: c1, c2
      integer :: n_m, m, n_r

      call get_bc_influence_matrix(psi1_Mloc, bcBot1, bcTop1)

      do n_m=nMstart,nMstop

         m = idx2m(n_m)

         if ( m /= 0 ) then

            c1=-influence_mat(1,1,n_m)*bcTop1(n_m) &
            &  -influence_mat(1,2,n_m)*bcBot1(n_m)
            c2=-influence_mat(2,1,n_m)*bcTop1(n_m) &
            &  -influence_mat(2,2,n_m)*bcBot1(n_m)

            do n_r=1,n_r_max
               psi1_Mloc(n_m,n_r)=psi1_Mloc(n_m,n_r)+c1*psi2_Mloc(n_m,n_r)+&
               &                  c2*psi3_Mloc(n_m,n_r)
               om1_Mloc(n_m,n_r) =om1_Mloc(n_m,n_r)+c1*om2_Mloc(n_m,n_r)+&
               &                  c2*om3_Mloc(n_m,n_r)
            end do

         end if

      end do

   end subroutine solve_influence_matrix
!------------------------------------------------------------------------------
   subroutine get_bc_influence_matrix_real(psi_tmp, bcBot, bcTop)

      !-- Input variable
      real(cp), intent(inout) :: psi_tmp(nMstart:nMstop, n_r_max)

      !-- Output variables
      real(cp), intent(out) :: bcBot(nMstart:nMstop)
      real(cp), intent(out) :: bcTop(nMstart:nMstop)

      !-- Local variables
      real(cp) :: fac
      integer :: n_m, n_cheb, m

      !-- Initialize to zero
      do n_m=nMstart,nMstop
         bcTop(n_m)=0.0_cp
         bcBot(n_m)=0.0_cp
      end do

      !-- Direct computation of the derivatives (input in Cheb space, output in 
      !-- physical space)
      do n_cheb=1,n_r_max

         if ( (n_cheb==1) .or. (n_cheb==n_r_max) ) then
            fac = rscheme%boundary_fac
         else
            fac = one
         end if

         do n_m=nMstart,nMstop

            m = idx2m(n_m)

            if ( m /= 0 ) then

               if ( ktopv == 1 ) then ! Stress-free
                  bcTop(n_m)=bcTop(n_m)+fac*(rscheme%d2rMat(1,n_cheb)-or1(1)* &
                  &                     rscheme%drMat(1,n_cheb))*psi_tmp(n_m,n_cheb)
               else
                  bcTop(n_m)=bcTop(n_m)+fac*rscheme%drMat(1,n_cheb)* &
                  &          psi_tmp(n_m,n_cheb)
               end if

               if ( kbotv == 1 ) then
                  bcBot(n_m)=bcBot(n_m)+fac*(rscheme%d2rMat(n_r_max,n_cheb)-  &
                  &          or1(n_r_max)*rscheme%drMat(n_r_max,n_cheb))*     &
                  &          psi_tmp(n_m,n_cheb)
               else
                  bcBot(n_m)=bcBot(n_m)+fac*rscheme%drMat(n_r_max,n_cheb)*  &
                  &          psi_tmp(n_m,n_cheb)
               end if

            end if

         end do

      end do

      !-- Multiply by cheb factor
      do n_m=nMstart,nMstop
         bcTop(n_m)=rscheme%rnorm * bcTop(n_m)
         bcBot(n_m)=rscheme%rnorm * bcBot(n_m)
      end do

   end subroutine get_bc_influence_matrix_real
!------------------------------------------------------------------------------
   subroutine get_bc_influence_matrix_complex(psi_tmp, bcBot, bcTop)

      !-- Input variable
      complex(cp), intent(inout) :: psi_tmp(nMstart:nMstop, n_r_max)

      !-- Output variables
      complex(cp), intent(out) :: bcBot(nMstart:nMstop)
      complex(cp), intent(out) :: bcTop(nMstart:nMstop)

      !-- Local variables
      real(cp) :: fac
      integer :: n_m, n_cheb, m

      !-- Initialize to zero
      do n_m=nMstart,nMstop
         bcTop(n_m)=zero
         bcBot(n_m)=zero
      end do

      !-- Direct computation of the derivatives (input in Cheb space, output in 
      !-- physical space)
      do n_cheb=1,n_r_max

         if ( (n_cheb==1) .or. (n_cheb==n_r_max) ) then
            fac = rscheme%boundary_fac
         else
            fac = one
         end if

         do n_m=nMstart,nMstop

            m = idx2m(n_m)

            if ( m /= 0 ) then

               if ( ktopv == 1 ) then ! Stress-free
                  bcTop(n_m)=bcTop(n_m)+fac*(rscheme%d2rMat(1,n_cheb)-or1(1)* &
                  &                     rscheme%drMat(1,n_cheb))*psi_tmp(n_m,n_cheb)
               else
                  bcTop(n_m)=bcTop(n_m)+fac*rscheme%drMat(1,n_cheb)* &
                  &          psi_tmp(n_m,n_cheb)
               end if

               if ( kbotv == 1 ) then
                  bcBot(n_m)=bcBot(n_m)+fac*(rscheme%d2rMat(n_r_max,n_cheb)- &
                  &          or1(n_r_max)*rscheme%drMat(n_r_max,n_cheb))*    &
                  &          psi_tmp(n_m,n_cheb)
               else
                  bcBot(n_m)=bcBot(n_m)+fac*rscheme%drMat(n_r_max,n_cheb)*  &
                  &          psi_tmp(n_m,n_cheb)
               end if

            end if

         end do

      end do

      !-- Multiply by cheb factor
      do n_m=nMstart,nMstop
         bcTop(n_m)=rscheme%rnorm * bcTop(n_m)
         bcBot(n_m)=rscheme%rnorm * bcBot(n_m)
      end do

   end subroutine get_bc_influence_matrix_complex
!------------------------------------------------------------------------------
end module update_psi_coll_dmat
