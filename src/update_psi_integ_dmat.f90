module update_psi_integ_dmat

   use precision_mod
   use parallel_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: one, zero, ci, half
   use namelists, only: kbotv, ktopv, alpha, r_cmb, r_icb, l_non_rot, CorFac, &
       &                l_ek_pump, ViscFac, ek, l_buo_imp
   use horizontal, only: hdif_V
   use radial_functions, only: rscheme, or1, or2, beta, ekpump, oheight, r
   use blocking, only: nMstart, nMstop, l_rank_has_m0
   use truncation, only: n_r_max, idx2m, m2idx, n_cheb_max
   use radial_der, only: get_ddr, get_dr
   use fields, only: work_Mloc, dom_Mloc
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use useful, only: abortRun
   use vort_balance, only: vort_bal_type
   use vp_balance, only: vp_bal_type
   use band_matrix, only: type_bandmat_real
   use bordered_matrix, only: type_bordmat_real
   use timers_mod, only: timers_type
   use chebsparselib, only: rmult2, intcheb1rmult1, intcheb2rmult2,    &
       &                    intcheb2rmult2laplrot, intcheb2rmult2lapl


   implicit none
   
   private

   logical,  allocatable :: lPsimat(:) ! Do we need to rebuild the matrices?
   complex(cp), allocatable :: rhs(:)  ! Complex right-hand-side
   real(cp), allocatable :: omfac(:,:) ! Precondition matrix
   real(cp), allocatable :: psifac(:,:) ! Precondition matrix
   real(cp), allocatable :: om2_Mloc(:,:)
   real(cp), allocatable :: om3_Mloc(:,:)
   real(cp), allocatable :: psi2_Mloc(:,:)
   real(cp), allocatable :: psi3_Mloc(:,:)
   real(cp), allocatable :: influence_matrix_Mloc(:,:,:)

   type(type_bordmat_real), allocatable :: LHS_om_mat(:)
   type(type_bordmat_real), allocatable :: LHS_psi_mat(:)
   type(type_bandmat_real), allocatable :: RHSIL_mat(:)
   type(type_bandmat_real) :: RHSE_mat, RHSI_mat

   interface  get_bc_influence_matrix
      module procedure get_bc_influence_matrix_complex
      module procedure get_bc_influence_matrix_real
   end interface get_bc_influence_matrix

   public :: update_psi_int_dmat, initialize_psi_integ_dmat,    &
   &         finalize_psi_integ_dmat, get_psi_rhs_imp_int_dmat, &
   &         finish_exp_psi_int_dmat

contains

   subroutine initialize_psi_integ_dmat
      !
      ! Memory allocation
      !


      !-- Local variables
      integer :: n_m, m

      !-- Allocate arrays that are needed for the influence matrix method
      allocate( om2_Mloc(nMstart:nMstop,n_r_max), om3_Mloc(nMstart:nMstop,n_r_max) )
      allocate( psi2_Mloc(nMstart:nMstop,n_r_max), psi3_Mloc(nMstart:nMstop,n_r_max) )
      om2_Mloc(:,:) =0.0_cp
      om3_Mloc(:,:) =0.0_cp
      psi2_Mloc(:,:)=0.0_cp
      psi3_Mloc(:,:)=0.0_cp
      bytes_allocated = bytes_allocated+4*(nMstop-nMstart+1)*n_r_max* &
      &                 SIZEOF_DEF_REAL

      allocate( influence_matrix_Mloc(2,2,nMstart:nMstop) )
      influence_matrix_Mloc(:,:,:)=0.0_cp
      bytes_allocated = bytes_allocated+4*(nMstop-nMstart+1)*SIZEOF_DEF_REAL

      !-- Allocate array of matrices (this is handy since it can store various bandwidths)
      allocate( RHSIL_mat(nMstart:nMstop) )
      allocate( LHS_om_mat(nMstart:nMstop) )
      allocate( LHS_psi_mat(nMstart:nMstop) )

      call RHSE_mat%initialize(4, 4, n_r_max)

      call RHSI_mat%initialize(6, 6, n_r_max)
      do n_m=nMstart,nMstop
         m = idx2m(n_m)
         call RHSIL_mat(n_m)%initialize(4, 4, n_r_max)
         call LHS_om_mat(n_m)%initialize(6, 6, 2, n_r_max)
         if ( m /= 0 ) call LHS_psi_mat(n_m)%initialize(4, 4, 2, n_r_max)
      end do

      !-- Allocate prefactors
      allocate( omfac(n_r_max, nMstart:nMstop) )
      allocate( psifac(n_r_max, nMstart:nMstop) )
      omfac(:,:) =one
      psifac(:,:)=one
      bytes_allocated = bytes_allocated + 2*n_r_max*(nMstop-nMstart+1)* &
      &                 SIZEOF_DEF_REAL

      !-- Fill matrices
      call get_rhs_exp_mat(RHSE_mat)
      call get_rhs_imp_mat(RHSI_mat)
      do n_m=nMstart,nMstop
         m = idx2m(n_m)
         call get_rhs_imp_lin_mat(RHSIL_mat(n_m), m)
         if ( m /= 0 ) then
            call get_lhs_psi_mat(LHS_psi_mat(n_m), psifac(:, n_m), m)
         end if
      end do

      allocate( lPsimat(nMstart:nMstop) )
      lPsimat(:)=.false.
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_LOGICAL

      allocate( rhs(n_r_max) )
      rhs(:)=zero
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_COMPLEX


   end subroutine initialize_psi_integ_dmat
!------------------------------------------------------------------------------
   subroutine finalize_psi_integ_dmat
      !
      ! Memory deallocation
      ! 

      !-- Local variable
      integer :: n_m, m

      call RHSI_mat%finalize()
      call RHSE_mat%finalize()
      do n_m=nMstart,nMstop
         m = idx2m(n_m)
         call RHSIL_mat(n_m)%finalize()
         call LHS_om_mat(n_m)%finalize()
         if ( m /= 0 ) call LHS_psi_mat(n_m)%finalize()
      end do
      deallocate( LHS_psi_mat, LHS_om_mat, RHSIL_mat )
      deallocate( psifac, omfac, rhs, lPsimat )

      deallocate( influence_matrix_Mloc )
      deallocate( om2_Mloc, om3_Mloc, psi2_Mloc, psi3_Mloc )

   end subroutine finalize_psi_integ_dmat
!------------------------------------------------------------------------------
   subroutine update_psi_int_dmat(psi_Mloc, om_Mloc, us_Mloc, up_Mloc,         &
              &                   buo_Mloc, domdt, vp_bal, vort_bal, tscheme,  &
              &                   lMat, timers)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lMat

      !-- Output variables
      complex(cp),         intent(out) :: psi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: up_Mloc(nMstart:nMstop,n_r_max)
      type(vp_bal_type),   intent(inout) :: vp_bal
      type(vort_bal_type), intent(inout) :: vort_bal
      type(type_tarray),   intent(inout) :: domdt
      complex(cp),         intent(inout) :: buo_Mloc(nMstart:nMstop,n_r_max)
      type(timers_type),   intent(inout) :: timers

      !-- Local variables
      real(cp) :: uphi0(n_r_max), om0(n_r_max)
      real(cp) :: h2, runStart, runStop, dm2
      integer :: n_r, n_m, n_cheb, m

      if ( lMat ) lPsimat(:)=.false.

      !-- Calculate first part of time-derivative of \omega if needed
      if ( vort_bal%l_calc .and. tscheme%istage == 1  ) then
         call vort_bal%initialize_domdt(om_Mloc,tscheme)
      end if

      !-- Calculate first part of time-derivative of axisymmetric u_\phi if needed
      if ( vp_bal%l_calc .and. tscheme%istage == 1  ) then
         call vp_bal%initialize_dvpdt(up_Mloc,tscheme)
      end if

      !-- If Ekman pumping is requested, normalisation is different
      !-- Hence buoyancy has to be multiplied by h^2
      if ( l_buo_imp ) then

         !-- Store buoyancy for the force balance
         if ( vort_bal%l_calc .and. tscheme%istage == tscheme%nstages  ) then
            do n_r=1,n_r_max
               do n_m=nMstart,nMstop
                  m = idx2m(n_m)
                  if ( m /= 0 ) then
                     vort_bal%buo(n_m,n_r)=buo_Mloc(n_m,n_r)/tscheme%dt(1)
                  end if
               end do
            end do
         end if

         !-- Transform buoyancy to Chebyshev space
         runStart = MPI_Wtime()
         call rscheme%costf1(buo_Mloc, nMstart, nMstop, n_r_max)
         runStop = MPI_Wtime()
         if ( runStop > runStart ) then
            timers%dct = timers%dct + (runStop-runStart)
            timers%n_dct_calls = timers%n_dct_calls + 1
         end if

         !-- Matrix-vector multiplication by the operator \int\int r^2:
         do n_m=nMstart,nMstop
            m = idx2m(n_m)

            if ( m > 0 ) then
               do n_cheb=1,n_r_max
                  rhs(n_cheb)=buo_Mloc(n_m,n_cheb)
               end do

               call RHSE_mat%mat_vec_mul(rhs)

               rhs(1)=zero
               rhs(2)=zero
               do n_cheb=1,n_r_max
                  buo_Mloc(n_m,n_cheb)=rhs(n_cheb)
               end do

            end if
         end do
      end if

      !-- Calculation of the implicit part
      call get_psi_rhs_imp_int_dmat(om_Mloc, up_Mloc,                        &
           &                        domdt%old(:,:,tscheme%istage),           &
           &                        domdt%impl(:,:,tscheme%istage), vp_bal,  &
           &                        tscheme%l_imp_calc_rhs(tscheme%istage))


      !-- Now assemble the right hand side and store it in work_Mloc
      call tscheme%set_imex_rhs(work_Mloc, domdt, nMstart, nMstop, n_r_max)


      do n_m=nMstart,nMstop

         m = idx2m(n_m)

         if ( .not. lPsimat(n_m) ) then
            call get_lhs_om_mat(tscheme, LHS_om_mat(n_m), omfac(:, n_m), m, &
                 &              timers%lu, timers%n_lu_calls)

            if ( m /= 0 ) then
               !-- Now solve the first intermediate problem
               !-- L_\omega \omega = 0; \omega(ro)=0; \omega(ri)=1
               !-- L_\psi \psi = -\omega; \psi(ro)=0; \psi(ri)=1
               call get_intermediate_sol(LHS_om_mat(n_m), LHS_psi_mat(n_m), m, &
                    &                    omfac(:,n_m), psifac(:,n_m), om2_Mloc,&
                    &                    psi2_Mloc, n_bc_set=1)
               !-- Now solve the second intermediate problem
               !-- L_\omega \omega = 0; \omega(ro)=1; \omega(ri)=0
               !-- L_\psi \psi = -\omega; \psi(ro)=0; \psi(ri)=1
               call get_intermediate_sol(LHS_om_mat(n_m), LHS_psi_mat(n_m), m, &
                    &                    omfac(:,n_m), psifac(:,n_m), om3_Mloc,&
                    &                    psi3_Mloc, n_bc_set=2)
            end if
         end if

         !-- Assemble RHS
         do n_cheb=1,n_r_max
            rhs(n_cheb)=work_Mloc(n_m,n_cheb)
            !-- Add buoyancy
            if ( l_buo_imp ) then
               if ( m /=  0 ) rhs(n_cheb)=rhs(n_cheb)+buo_Mloc(n_m,n_cheb)
            end if
         end do

         !-- Multiply rhs by precond matrix
         do n_cheb=1,n_r_max
            rhs(n_cheb)=rhs(n_cheb)*omfac(n_cheb,n_m)
         end do

         runStart = MPI_Wtime()
         call LHS_om_mat(n_m)%solve(rhs,n_r_max)
         runStop = MPI_Wtime()
         if ( runStop > runStart ) then
            timers%solve = timers%solve + (runStop-runStart)
            timers%n_solve_calls = timers%n_solve_calls+1
         end if

         if ( m == 0 ) then
            do n_cheb=1,n_r_max
               uphi0(n_cheb)=real(rhs(n_cheb))
            end do
         else
            do n_cheb=1,n_r_max
               om_Mloc(n_m,n_cheb)=rhs(n_cheb)
            end do

            !-- Now solve for L_\psi \psi = -\omega
            !-- DGEMM to assemble right-hand side for psi equations
            call RHSE_mat%mat_vec_mul(rhs)
            !-- Boundary conditions for psi
            rhs(1)=0.0_cp
            rhs(2)=0.0_cp
            do n_cheb=1,n_r_max
               !-- Minus signs comes from \Del\psi = - \omega
               rhs(n_cheb)= -rhs(n_cheb)*psifac(n_cheb,n_m)
            end do

            call LHS_psi_mat(n_m)%solve(rhs,n_r_max)

            do n_cheb=1,n_r_max
               psi_Mloc(n_m,n_cheb)=rhs(n_cheb)
            end do
         end if  ! m /= 0

      end do

      !-- set cheb modes > n_cheb_max to zero (dealiasing)
      !if ( n_cheb_max < n_r_max ) then ! fill with zeros !
      !   do n_cheb=n_cheb_max+1,n_r_max
      !      do n_m=nMstart,nMstop
      !         m = idx2m(n_m)
      !         if ( m == 0 ) then
      !            uphi0(n_cheb)=0.0_cp
      !         else
      !            psi_Mloc(n_m,n_cheb)=zero
      !            om_Mloc(n_m,n_cheb) =zero
      !            if ( .not. lPsiMat(n_m) ) then
      !               om2_Mloc(n_m,n_cheb) =0.0_cp
      !               om3_Mloc(n_m,n_cheb) =0.0_cp
      !               psi2_Mloc(n_m,n_cheb)=0.0_cp
      !               psi3_Mloc(n_m,n_cheb)=0.0_cp
      !               lPsimat(n_m)=.true.
      !            end if
      !         end if ! m /= 0
      !      end do
      !   end do
      !end if

      !-- Bring uphi0 to the physical space
      if ( l_rank_has_m0 ) then
         call get_dr(uphi0, om0, n_r_max, rscheme, l_dct=.false.)
         call rscheme%costf1(uphi0, n_r_max)
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

      !-- Get the radial derivative of psi to calculate uphi from psi
      call get_dr(psi_Mloc, work_Mloc, nMstart, nMstop, n_r_max, rscheme, &
           &      l_dct_in=.false.)

      !-- Bring \psi and \omega back the physical space using DCTs
      runStart = MPI_Wtime()
      call rscheme%costf1(psi_Mloc, nMstart, nMstop, n_r_max)
      call rscheme%costf1(om_Mloc, nMstart, nMstop, n_r_max)
      runStop = MPI_Wtime()
      if ( runStop > runStart ) then
         timers%dct = timers%dct + (runStop-runStart)
         timers%n_dct_calls = timers%n_dct_calls + 2
      end if


      if ( l_non_rot ) then
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)

               if ( m == 0 ) then
                  us_Mloc(n_m,n_r)=0.0_cp
                  up_Mloc(n_m,n_r)=uphi0(n_r)
                  om_Mloc(n_m,n_r)=om0(n_r)+or1(n_r)*uphi0(n_r)
               else ! om_Mloc is already defined at this stage
                  us_Mloc(n_m,n_r)=ci*real(m,cp)*or1(n_r)*psi_Mloc(n_m,n_r)
                  up_Mloc(n_m,n_r)=-work_Mloc(n_m,n_r)
               end if
            end do
         end do
      else
         do n_r=1,n_r_max
            h2 = r_cmb*r_cmb-r(n_r)*r(n_r)
            do n_m=nMstart,nMstop
               m = idx2m(n_m)

               if ( m == 0 ) then
                  us_Mloc(n_m,n_r)=0.0_cp
                  up_Mloc(n_m,n_r)=uphi0(n_r)
                  om_Mloc(n_m,n_r)=om0(n_r)+or1(n_r)*uphi0(n_r)
               else ! om_Mloc is already known at this stage
                  us_Mloc(n_m,n_r)=ci*real(m,cp)*or1(n_r)*h2*   psi_Mloc(n_m,n_r)
                  up_Mloc(n_m,n_r)=-h2*                        work_Mloc(n_m,n_r) &
                  &                +3.0_cp*r(n_r)*              psi_Mloc(n_m,n_r)

               end if
            end do
         end do
      end if

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(domdt, nMstart, nMstop, n_r_max)

      !-- Finish calculation of du_\phi/dt if requested
      if ( vp_bal%l_calc .and. tscheme%istage == tscheme%nstages  ) then
         call vp_bal%finalize_dvpdt(up_Mloc, tscheme)
      end if

      !-- Finish calculation of d\omega/dt if requested
      !-- Compute the viscous term and Coriolis
      if ( vort_bal%l_calc .and. tscheme%istage == tscheme%nstages  ) then
         call get_ddr(om_Mloc, dom_Mloc, work_Mloc, nMstart, nMstop, n_r_max, &
              &       rscheme)
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               dm2 = real(m,cp)*real(m,cp)
               if ( m > 0 ) then
                  vort_bal%visc(n_m,n_r)=ViscFac*(  work_Mloc(n_m,n_r)+&
                  &                   or1(n_r)*      dom_Mloc(n_m,n_r)-&
                  &               dm2*or2(n_r)*       om_Mloc(n_m,n_r) )
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
            end do
         end do
         call vort_bal%finalize_domdt(om_Mloc, tscheme)
      end if

   end subroutine update_psi_int_dmat
!------------------------------------------------------------------------------
   subroutine finish_exp_psi_int_dmat(psi_Mloc, us_Mloc, up_Mloc, om_Mloc, &
              &                       dVsOm_Mloc, buo_Mloc, dom_exp_last,  &
              &                       vp_bal, vort_bal)

      !-- Input variables
      complex(cp), intent(in) :: psi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: dVsOm_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: buo_Mloc(nMstart:nMstop,n_r_max)

      !-- Output variables
      type(vp_bal_type),   intent(inout) :: vp_bal
      type(vort_bal_type), intent(inout) :: vort_bal
      complex(cp),         intent(inout) :: dom_exp_last(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_m, m, n_cheb
      real(cp) :: ekp_fac

      !-- Finish calculation of advection
      call rscheme%costf1( dVsOm_Mloc, nMstart, nMstop, n_r_max )
      do n_cheb=n_cheb_max+1,n_r_max
         do n_m=nMstart,nMstop
            dVsOm_Mloc(n_m,n_cheb)=zero
         end do
      end do
      call get_dr( dVsOm_Mloc, work_Mloc, nMstart, nMstop, n_r_max, &
           &       rscheme, nocopy=.true., l_dct_in=.false.)

      !-- Finish calculation of the explicit part for current time step
      do n_r=1,n_r_max
         do n_m=nMstart, nMstop
            m = idx2m(n_m)
            if ( m /= 0 ) then
               dom_exp_last(n_m,n_r)=      dom_exp_last(n_m,n_r)-   &
               &                     or1(n_r)*work_Mloc(n_m,n_r)

               !-- Store the advection term when vorticity balance is requested
               if ( vort_bal%l_calc ) then
                  vort_bal%adv(n_m,n_r)=dom_exp_last(n_m,n_r)
               end if

               !-- If Coriolis force is required it is added here:
               if ( .not. l_non_rot ) then
                  dom_exp_last(n_m,n_r)=                 dom_exp_last(n_m,n_r)-  &
                  &                     CorFac*ci*real(m,cp)*psi_Mloc(n_m,n_r)
               end if

               !-- If Buoyancy is treated exiplicitly:
               if ( .not. l_buo_imp ) then
                  !-- Store buoyancy for the vorticity balance:
                  if ( vort_bal%l_calc ) then
                     vort_bal%buo(n_m,n_r)=buo_Mloc(n_m,n_r)
                  end if
                  dom_exp_last(n_m,n_r)=dom_exp_last(n_m,n_r)+buo_Mloc(n_m,n_r)
               end if
            end if
         end do
      end do

      if ( l_rank_has_m0 .and. vp_bal%l_calc ) then
         do n_r=1,n_r_max
            vp_bal%rey_stress(n_r)=real(dom_exp_last(m2idx(0),n_r))
         end do
      end if


      !-- Add Ekman pumping as an explicit term if this is requested
      if ( l_ek_pump ) then
         do n_r=1,n_r_max
            ekp_fac = CorFac*ekpump(n_r)
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               if ( m == 0 ) then
                  dom_exp_last(n_m,n_r)= dom_exp_last(n_m,n_r) -  &
                  &                   ekp_fac*up_Mloc(n_m,n_r)
               else 
                  dom_exp_last(n_m,n_r)=   dom_exp_last(n_m,n_r) +           &
                  &                  ekp_fac*( -om_Mloc(n_m,n_r) +           &
                  &           half*beta(n_r)*   up_Mloc(n_m,n_r) +           &
                  &   beta(n_r)*(-ci*real(m,cp)+5.0_cp*r_cmb*oheight(n_r))*  &
                  &                             us_Mloc(n_m,n_r) )
               end if
            end do
         end do
      end if

      !-- Transform the explicit part to chebyshev space
      call rscheme%costf1(dom_exp_last, nMstart, nMstop, n_r_max)

      do n_cheb=n_cheb_max+1,n_r_max
         do n_m=nMstart,nMstop
            dom_exp_last(n_m,n_cheb)=zero
         end do
      end do


      !-- Matrix-vector multiplication by the operator \int\int r^2:
      do n_m=nMstart,nMstop
         m = idx2m(n_m)

         do n_cheb=1,n_r_max
            rhs(n_cheb)=dom_exp_last(n_m,n_cheb)
         end do

         call RHSE_mat%mat_vec_mul(rhs)
         rhs(1)=zero
         rhs(2)=zero

         do n_cheb=1,n_r_max
            dom_exp_last(n_m,n_cheb)=rhs(n_cheb)
         end do

      end do

   end subroutine finish_exp_psi_int_dmat
!------------------------------------------------------------------------------
   subroutine get_psi_rhs_imp_int_dmat(om_Mloc, up_Mloc, om_old,  &
              &                        dom_imp_Mloc_last, vp_bal, &
              &                        l_calc_lin_rhs)

      !-- Input variables
      complex(cp), intent(in) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      logical,     intent(in) :: l_calc_lin_rhs

      !-- Output variables
      complex(cp),       intent(inout) :: dom_imp_Mloc_last(nMstart:nMstop,n_r_max)
      complex(cp),       intent(inout) :: om_old(nMstart:nMstop,n_r_max)
      type(vp_bal_type), intent(inout) :: vp_bal

      !-- Local variables
      real(cp) :: duphi0(n_r_max), d2uphi0(n_r_max), uphi0(n_r_max)
      integer :: n_r, n_m, m, m0, n_cheb

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m == 0 ) then
               om_old(n_m,n_r)=up_Mloc(n_m,n_r)
            else
               om_old(n_m,n_r)=om_Mloc(n_m,n_r)
            end if
         end do
      end do

      !-- Transform the implicit part to chebyshev space
      call rscheme%costf1(om_old, nMstart, nMstop, n_r_max)

      !-- Matrix-vector multiplication by the operator \int^2 r^2.
      do n_m=nMstart,nMstop
         m = idx2m(n_m)

         do n_cheb=1,n_r_max
            rhs(n_cheb)= om_old(n_m,n_cheb)
         end do
         call RHSI_mat%mat_vec_mul(rhs)

         rhs(1)=zero 
         rhs(2)=zero

         do n_cheb=1,n_r_max
            om_old(n_m,n_cheb)=rhs(n_cheb)
         end do

      end do

      !-- Calculate and store vphi force balance if needed
      if ( l_rank_has_m0 .and. vp_bal%l_calc ) then
         m0 = m2idx(0)
         do n_r=1,n_r_max
            uphi0(n_r)=real(up_Mloc(m0, n_r),kind=cp)
         end do
         call get_ddr(uphi0, duphi0, d2uphi0, n_r_max, rscheme)
         do n_r=1,n_r_max
            vp_bal%visc(n_r)=ViscFac*hdif_V(m0)*(d2uphi0(n_r)+       &
            &                or1(n_r)*duphi0(n_r)-or2(n_r)*uphi0(n_r))
            vp_bal%pump(n_r)=-CorFac*ekpump(n_r)*uphi0(n_r)
         end do
      end if

      if ( l_calc_lin_rhs ) then

         !-- Copy om_Mloc into work_Mloc (or uphi(m=0) for axisymmetric equation)
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               if ( m == 0 ) then
                  work_Mloc(n_m,n_r)=up_Mloc(n_m,n_r)
               else
                  work_Mloc(n_m,n_r)=om_Mloc(n_m,n_r)
               end if
            end do
         end do

         !-- Transform work_Mloc to Cheb space
         call rscheme%costf1(work_Mloc, nMstart, nMstop, n_r_max)

         !-- Matrix-vector multiplication by the LHS operator
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            do n_cheb=1,n_r_max
               rhs(n_cheb)=work_Mloc(n_m,n_cheb)
            end do
            call RHSIL_mat(n_m)%mat_vec_mul(rhs)
            rhs(1)=zero ! vphi equation has only 2 BCs
            rhs(2)=zero
            do n_cheb=1,n_r_max
               work_Mloc(n_m,n_cheb)=rhs(n_cheb)
            end do
         end do

         !-- Finally assemble the right hand side
         do n_cheb=1,n_r_max
            do n_m=nMstart,nMstop
               dom_imp_Mloc_last(n_m,n_cheb)=work_Mloc(n_m,n_cheb)
            end do
         end do

      end if 

   end subroutine get_psi_rhs_imp_int_dmat
!------------------------------------------------------------------------------
   subroutine get_lhs_om_mat(tscheme, A_mat, omMat_fac, m, time_lu, n_lu_calls)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme    ! time step
      integer,             intent(in) :: m          ! Azimuthal wavenumber

      !-- Output variables
      type(type_bordmat_real), intent(inout) :: A_mat
      real(cp),                intent(inout) :: omMat_fac(A_mat%nlines)
      real(cp),                intent(inout) :: time_lu
      integer,                 intent(inout) :: n_lu_calls

      !-- Local variables
      real(cp) :: stencilA4(A_mat%nbands)
      integer :: n_r, i_r, n_band, n_b, n_m
      real(cp) :: a, b, runStart, runStop

      n_m = m2idx(m)

      !-- We have to fill A3 with zeros again otherwise on the next iteration
      !-- with a different dt there might be some issues with spurious values
      do n_r=1,A_mat%nlines_band
         do n_b=1,A_mat%ntau
            A_mat%A3(n_r,n_b)=0.0_cp
         end do
      end do

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      !-- Fill A4 banded-block
      do n_r=1,A_mat%nlines_band
         i_r = n_r+A_mat%ntau

         !-- Define the equations
         if ( m == 0 ) then
            stencilA4 = intcheb2rmult2(a,b,i_r-1,A_mat%nbands)-               &
            &           tscheme%wimp_lin(1)*ViscFac*hdif_V(n_m)* (            &
            &                                 rmult2(a,b,i_r-1,A_mat%nbands)- &
            &                  3.0_cp*intcheb1rmult1(a,b,i_r-1,A_mat%nbands) )
         else
            stencilA4 = intcheb2rmult2(a,b,i_r-1,A_mat%nbands)-         &
            &           tscheme%wimp_lin(1)*ViscFac*hdif_V(n_m)*        &
            &           intcheb2rmult2lapl(a,b,m,i_r-1,A_mat%nbands)
         end if

         !-- Roll the array for band storage
         do n_band=1,A_mat%nbands
            if ( n_r+A_mat%ku+1-n_band <= A_mat%nlines_band .and. n_r+A_mat%ku+1-n_band >= 1 ) then
               A_mat%A4(A_mat%kl+n_band,n_r+A_mat%ku+1-n_band) = rscheme%rnorm* &
               &             stencilA4(n_band)
            end if
         end do

         omMat_fac(n_r+A_mat%ntau)=one/rscheme%rnorm/maxval(abs(stencilA4))
      end do

      !-- Fill A3
      do n_r=1,A_mat%nlines_band
         i_r = n_r+A_mat%ntau

         if ( m == 0 ) then
            stencilA4 = intcheb2rmult2(a,b,i_r-1,A_mat%nbands)-               &
            &           tscheme%wimp_lin(1)*ViscFac*hdif_V(n_m)* (            &
            &                                 rmult2(a,b,i_r-1,A_mat%nbands)- &
            &                  3.0_cp*intcheb1rmult1(a,b,i_r-1,A_mat%nbands) )
         else
            stencilA4 = intcheb2rmult2(a,b,i_r-1,A_mat%nbands)-         &
            &           tscheme%wimp_lin(1)*ViscFac*hdif_V(n_m)*        & 
            &           intcheb2rmult2lapl(a,b,m,i_r-1,A_mat%nbands)  
         end if

         !-- Only the lower bands can contribute to the matrix A3
         do n_band=1,A_mat%kl
            if ( n_r <= n_band .and. n_r+A_mat%ntau-n_band >= 1 ) then
               A_mat%A3(n_r,n_r+A_mat%ntau-n_band) = rscheme%rnorm*       &
               &                                     stencilA4(A_mat%ku+1+n_band)
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
            else
               A_mat%A2(1,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%rMat(1,n_r)
               A_mat%A2(2,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%rMat(2,n_r)
            end if
         end do
      end if

      !-- Cheb factor for boundary conditions
      do n_r=1,A_mat%ntau
         A_mat%A1(n_r,1)                =rscheme%boundary_fac*A_mat%A1(n_r,1)
         if ( A_mat%nlines_band == n_r_max ) then
            A_mat%A2(n_r,A_mat%nlines_band)=rscheme%boundary_fac*A_mat%A2(n_r,A_mat%nlines_band)
         end if
      end do

      !-- Continue to assemble precond matrix
      do n_r=1,A_mat%ntau
         omMat_fac(n_r)=one/maxval(abs(A_mat%A2(n_r,:)))
      end do

      !-- Multiply the lines of the matrix by precond
      do n_r=1,A_mat%ntau
         A_mat%A1(n_r,:) = A_mat%A1(n_r,:)*omMat_fac(n_r)
         A_mat%A2(n_r,:) = A_mat%A2(n_r,:)*omMat_fac(n_r)
      end do

      do n_r=A_mat%ntau+1, A_mat%nlines
         A_mat%A3(n_r-A_mat%ntau,:) = A_mat%A3(n_r-A_mat%ntau,:)* omMat_fac(n_r)
      end do

      do n_r=1,A_mat%nlines_band
         i_r = n_r+A_mat%ntau

         do n_band=1,A_mat%nbands
            if ( n_r+A_mat%ku+1-n_band <= A_mat%nlines_band .and. &
            &                         n_r+A_mat%ku+1-n_band >= 1 ) then
               A_mat%A4(A_mat%kl+n_band,n_r+A_mat%ku+1-n_band) =  &
               & A_mat%A4(A_mat%kl+n_band,n_r+A_mat%ku+1-n_band)* &
               &                      omMat_fac(n_r+A_mat%ntau)
            end if
         end do

      end do

      ! if ( m == 10 ) call A_mat%write()

      !-- LU factorisation
      runStart = MPI_Wtime()
      call A_mat%prepare_LU()
      runStop = MPI_Wtime()
      if ( runStop > runStart ) then
         time_lu = time_lu + (runStop-runStart)
         n_lu_calls = n_lu_calls + 1
      end if

   end subroutine get_lhs_om_mat
!------------------------------------------------------------------------------
   subroutine get_lhs_psi_mat(A_mat, psiMat_fac, m)

      !-- Input variables
      integer, intent(in) :: m

      !-- Output variables
      type(type_bordmat_real), intent(inout) :: A_mat
      real(cp),                intent(inout) :: psiMat_fac(A_mat%nlines)

      !-- Local variables
      real(cp) :: stencilA4(A_mat%nbands)
      integer :: n_r, i_r, n_band, n_b
      real(cp) :: a, b, monen

      !-- We have to fill A3 with zeros again otherwise on the next iteration
      !-- with a different dt there might be some issues with spurious values
      do n_r=1,A_mat%nlines_band
         do n_b=1,A_mat%ntau
            A_mat%A3(n_r,n_b)=0.0_cp
         end do
      end do

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      !-- Fill A4 banded-block
      do n_r=1,A_mat%nlines_band
         i_r = n_r+A_mat%ntau

         !-- Define the equations
         if ( l_non_rot ) then
            stencilA4 = intcheb2rmult2lapl(a,b,m,i_r-1,A_mat%nbands)
         else
            stencilA4 = intcheb2rmult2laplrot(a,b,m,i_r-1,A_mat%nbands)
         end if

         !-- Roll the array for band storage
         do n_band=1,A_mat%nbands
            if ( n_r+A_mat%ku+1-n_band <= A_mat%nlines_band .and. n_r+A_mat%ku+1-n_band >= 1 ) then
               A_mat%A4(A_mat%kl+n_band,n_r+A_mat%ku+1-n_band) = rscheme%rnorm* &
               &             stencilA4(n_band)
            end if
         end do

         psiMat_fac(n_r+A_mat%ntau)=one/rscheme%rnorm/maxval(abs(stencilA4))
      end do

      !-- Fill A3
      do n_r=1,A_mat%nlines_band
         i_r = n_r+A_mat%ntau

         !-- Define the equations
         if ( l_non_rot ) then
            stencilA4 = intcheb2rmult2lapl(a,b,m,i_r-1,A_mat%nbands)  
         else
            stencilA4 = intcheb2rmult2laplrot(a,b,m,i_r-1,A_mat%nbands)  
         end if

         !-- Only the lower bands can contribute to the matrix A3
         do n_band=1,A_mat%kl
            if ( n_r <= n_band .and. n_r+A_mat%ntau-n_band >= 1 ) then
               A_mat%A3(n_r,n_r+A_mat%ntau-n_band) = rscheme%rnorm*       &
               &                                     stencilA4(A_mat%ku+1+n_band)
            end if
         end do
      end do

      !-- Add the tau Lines for boundary conditions into A1 and A2
      do n_r=1,A_mat%nlines
         if ( mod(n_r,2)==1 ) then
            monen = one
         else
            monen = -one
         end if

         if ( n_r <= A_mat%ntau ) then
            A_mat%A1(1,n_r)=rscheme%rnorm
            A_mat%A1(2,n_r)=rscheme%rnorm*monen
         else
            A_mat%A2(1,n_r-A_mat%ntau)=rscheme%rnorm
            A_mat%A2(2,n_r-A_mat%ntau)=rscheme%rnorm*monen
         end if
      end do

      !-- Cheb factor for boundary conditions
      do n_r=1,A_mat%ntau
         A_mat%A1(n_r,1)                =rscheme%boundary_fac*A_mat%A1(n_r,1)
         if ( A_mat%nlines_band == n_r_max ) then
            A_mat%A2(n_r,A_mat%nlines_band)=rscheme%boundary_fac*A_mat%A2(n_r,A_mat%nlines_band)
         end if
      end do

      !-- Continue to assemble precond matrix
      do n_r=1,A_mat%ntau
         psiMat_fac(n_r)=one/maxval(abs(A_mat%A2(n_r,:)))
      end do

      !-- Multiply the lines of the matrix by precond
      do n_r=1,A_mat%ntau
         A_mat%A1(n_r,:) = A_mat%A1(n_r,:)*psiMat_fac(n_r)
         A_mat%A2(n_r,:) = A_mat%A2(n_r,:)*psiMat_fac(n_r)
      end do

      do n_r=A_mat%ntau+1, A_mat%nlines
         A_mat%A3(n_r-A_mat%ntau,:) = A_mat%A3(n_r-A_mat%ntau,:)* psiMat_fac(n_r)
      end do

      do n_r=1,A_mat%nlines_band
         i_r = n_r+A_mat%ntau

         do n_band=1,A_mat%nbands
            if ( n_r+A_mat%ku+1-n_band <= A_mat%nlines_band .and. &
            &                         n_r+A_mat%ku+1-n_band >= 1 ) then
               A_mat%A4(A_mat%kl+n_band,n_r+A_mat%ku+1-n_band) =  &
               & A_mat%A4(A_mat%kl+n_band,n_r+A_mat%ku+1-n_band)* &
               &                      psiMat_fac(n_r+A_mat%ntau)
            end if
         end do

      end do

      !-- LU factorisation
      call A_mat%prepare_LU()

   end subroutine get_lhs_psi_mat
!------------------------------------------------------------------------------
   subroutine get_rhs_exp_mat(D_mat)
      !
      ! This corresponds to the matrix that goes in front of the explicit terms
      !

      !-- Output variable
      type(type_bandmat_real), intent(inout) :: D_mat

      !-- Local variables
      real(cp) :: stencilD(D_mat%nbands)
      real(cp) :: a, b
      integer :: n_band, n_r, i_r, n_bounds

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      n_bounds = 2 ! 2 boundary conditions

      !-- Fill right-hand side matrix
      do n_r=1,D_mat%nlines
         i_r = n_r+n_bounds

         !-- Define right-hand side equations
         stencilD = intcheb2rmult2(a,b,i_r-1,D_mat%nbands)

         !-- Roll array for band storage
         do n_band=1,D_mat%nbands
            if ( i_r+D_mat%ku+1-n_band <= D_mat%nlines .and. &
               & i_r+D_mat%ku+1-n_band >= 1 ) then
               D_mat%dat(n_band,i_r+D_mat%ku+1-n_band) = rscheme%rnorm*&
               &                                         stencilD(n_band)
            end if
         end do
      end do

   end subroutine get_rhs_exp_mat
!------------------------------------------------------------------------------
   subroutine get_rhs_imp_mat(B_mat)
      !
      ! This matrix corresponds to the that goes in front of the d/dt term
      ! in the right-hand-side
      !

      !-- Output variable
      type(type_bandmat_real), intent(inout) :: B_mat

      !-- Local variables
      real(cp) :: stencilB(B_mat%nbands)
      real(cp) :: a, b
      integer :: n_band, n_r, i_r, n_bounds

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      n_bounds = 2 

      !-- Fill right-hand side matrix
      do n_r=1,B_mat%nlines
         i_r = n_r+n_bounds

         !-- Define right-hand side equations
         stencilB = intcheb2rmult2(a,b,i_r-1,B_mat%nbands)

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
      type(type_bandmat_real), intent(inout) :: Cmat

      !-- Local variables
      real(cp) :: stencilC(Cmat%nbands)
      real(cp) :: a, b
      integer :: n_band, n_r, i_r, n_bounds, n_m

      n_m = m2idx(m)

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      n_bounds = 2 ! Only 2 BCs

      !-- Fill right-hand side matrix
      do n_r=1,Cmat%nlines
         i_r = n_r+n_bounds

         !-- Define right-hand side equations
         if ( m == 0 ) then
            stencilC =  ViscFac *hdif_V(n_m)* (                     &
            &                        rmult2(a,b,i_r-1,Cmat%nbands)- &
            &         3.0_cp*intcheb1rmult1(a,b,i_r-1,Cmat%nbands) )
         else
            stencilC = ViscFac*hdif_V(n_m)* &
            &          intcheb2rmult2lapl(a,b,m,i_r-1,Cmat%nbands)
         end if

         !-- Roll array for band storage
         do n_band=1,Cmat%nbands
            if ( i_r+Cmat%ku+1-n_band <= Cmat%nlines .and. i_r+Cmat%ku+1-n_band >= 1 ) then
               Cmat%dat(n_band,i_r+Cmat%ku+1-n_band) = rscheme%rnorm* &
               &                                       stencilC(n_band)
            end if
         end do
      end do

   end subroutine get_rhs_imp_lin_mat
!------------------------------------------------------------------------------
   subroutine get_intermediate_sol(om_mat, psi_mat, m, omMat_fac, psiMat_fac, &
                                   om_tmp, psi_tmp, n_bc_set)

      !-- Input variables
      integer,                 intent(in) :: m          ! Azimuthal wavenumber
      type(type_bordmat_real), intent(in) :: om_mat
      type(type_bordmat_real), intent(in) :: psi_mat
      real(cp),                intent(in) :: omMat_fac(om_mat%nlines)
      real(cp),                intent(in) :: psiMat_fac(psi_mat%nlines)
      integer,                 intent(in) :: n_bc_set

      !-- Output variables
      real(cp),  intent(inout) :: om_tmp(nMstart:nMstop, n_r_max)
      real(cp),  intent(inout) :: psi_tmp(nMstart:nMstop, n_r_max)

      !-- Local variables
      real(cp) :: rhs(om_mat%nlines)
      integer :: n_m, n_cheb

      n_m = m2idx(m)

      !-- First fill with zeros
      do n_cheb=1,om_mat%nlines
         rhs(n_cheb)=0.0_cp
      end do

      !-- Apply BCs on the 2 first points
      if ( n_bc_set == 1 ) then
         rhs(1)=0.0_cp
         rhs(2)=1.0_cp
      else if ( n_bc_set == 2 ) then
         rhs(1)=1.0_cp
         rhs(2)=0.0_cp
      end if

      !-- Apply precondition matrix
      do n_cheb=1,2
         rhs(n_cheb)=rhs(n_cheb)*omMat_fac(n_cheb)
      end do

      !-- Solve for omega an psi
      if ( m /= 0 ) then
         call om_mat%solve(rhs, om_mat%nlines)
         do n_cheb=1,om_mat%nlines
            om_tmp(n_m,n_cheb)=rhs(n_cheb)
         end do

         !-- DGEMM to assemble the rhs for psi equation
         call RHSE_mat%mat_vec_mul(rhs)
         !-- BC for psi
         rhs(1)=0.0_cp
         rhs(2)=0.0_cp

         do n_cheb=1,om_mat%nlines
            !-- Minus sign since \Del psi = -\omega
            rhs(n_cheb)=-rhs(n_cheb)*psiMat_fac(n_cheb)
         end do

         !--  Now solve for psi
         call psi_mat%solve(rhs, psi_mat%nlines)
         do n_cheb=1,psi_mat%nlines
            psi_tmp(n_m,n_cheb)=rhs(n_cheb)
         end do
      end if

   end subroutine get_intermediate_sol
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
      real(cp) :: fac, dn2
      real(cp) :: d4top(n_r_max)
      integer :: n_m, n_cheb, m

      do n_cheb=1,n_r_max
         dn2 = real(n_cheb-1,cp)*real(n_cheb-1,cp)
         d4top(n_cheb)=16.0_cp/105.0_cp*dn2*(dn2-one)*(dn2-4.0_cp)*(dn2-9.0_cp)
      end do

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
                  if ( l_non_rot ) then
                     bcTop(n_m)=bcTop(n_m)+fac*rscheme%drMat(1,n_cheb)* &
                     &          psi_tmp(n_m,n_cheb)
                  else
                     bcTop(n_m)=bcTop(n_m)+fac*d4top(n_cheb)*psi_tmp(n_m,n_cheb)
                  end if
               end if

               if ( kbotv == 1 ) then
                  bcBot(n_m)=bcBot(n_m)+fac*(rscheme%d2rMat(2,n_cheb)-or1(n_r_max)*&
                  &                     rscheme%drMat(2,n_cheb))*psi_tmp(n_m,n_cheb)
               else
                  bcBot(n_m)=bcBot(n_m)+fac*rscheme%drMat(2,n_cheb)*  &
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
      real(cp) :: fac, dn2
      integer :: n_m, n_cheb, m
      real(cp) :: d4top(n_r_max)

      do n_cheb=1,n_r_max
         dn2 = real(n_cheb-1,cp)*real(n_cheb-1,cp)
         d4top(n_cheb)=16.0_cp/105.0_cp*dn2*(dn2-one)*(dn2-4.0_cp)*(dn2-9.0_cp)
      end do


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
                  if ( l_non_rot ) then
                     bcTop(n_m)=bcTop(n_m)+fac*rscheme%drMat(1,n_cheb)* &
                     &          psi_tmp(n_m,n_cheb)
                  else
                     bcTop(n_m)=bcTop(n_m)+fac*d4top(n_cheb)*psi_tmp(n_m,n_cheb)
                  end if
               end if

               if ( kbotv == 1 ) then
                  bcBot(n_m)=bcBot(n_m)+fac*(rscheme%d2rMat(2,n_cheb)-or1(n_r_max)*&
                  &                     rscheme%drMat(2,n_cheb))*psi_tmp(n_m,n_cheb)
               else
                  bcBot(n_m)=bcBot(n_m)+fac*rscheme%drMat(2,n_cheb)*  &
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
end module update_psi_integ_dmat
