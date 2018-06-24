module update_psi_integ_smat

   use precision_mod
   use parallel_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: one, zero, ci, half
   use outputs, only: vp_bal_type
   use namelists, only: kbotv, ktopv, alpha, r_cmb, r_icb, l_non_rot, CorFac, &
       &                l_ek_pump, ViscFac, l_coriolis_imp, ek, l_buo_imp,    &
       &                l_galerkin
   use radial_functions, only: rscheme, or1, or2, beta, ekpump, oheight, r
   use blocking, only: nMstart, nMstop, l_rank_has_m0
   use truncation, only: n_r_max, idx2m, m2idx, n_cheb_max
   use radial_der, only: get_ddr, get_dr
   use fields, only: work_Mloc
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use useful, only: abortRun
   use galerkin
   use band_matrix, only: type_bandmat_complex, type_bandmat_real, band_band_product
   use bordered_matrix, only: type_bordmat_complex
   use chebsparselib, only: intcheb4rmult4lapl2, intcheb4rmult4lapl,    &
       &                    intcheb4rmult4, rmult2, intcheb1rmult1,     &
       &                    intcheb2rmult2, intcheb4rmult4laplrot2,     &
       &                    intcheb4rmult4laplrot, intcheb4rmult4hmult2,&
       &                    intcheb4rmult4hmult2laplrot,                &
       &                    intcheb4rmult4hmult2laplrot2,               &
       &                    intcheb2rmult2hmult2, intcheb2rmult2hmult2laplm1


   implicit none
   
   private

   logical,  allocatable :: lPsimat(:) ! Do we need to rebuild the matrices?
   complex(cp), allocatable :: rhs(:)  ! Complex right-hand-side
   real(cp), allocatable :: psifac(:,:) ! Precondition matrix

   type(type_bordmat_complex), allocatable :: LHS_mat_tau(:)
   type(type_bandmat_complex), allocatable :: LHS_mat_gal(:)
   type(type_bandmat_complex), allocatable :: RHSIL_mat(:)
   type(type_bandmat_real), allocatable :: RHSI_mat(:)
   type(type_bandmat_real) :: RHSE_mat(2), gal_sten(2)

   public :: update_psi_int_smat, initialize_psi_integ_smat,    &
   &         finalize_psi_integ_smat, get_psi_rhs_imp_int_smat, &
   &         finish_exp_psi_int_smat

contains

   subroutine initialize_psi_integ_smat
      !
      ! Memory allocation
      !


      !-- Local variables
      integer :: n_m, m

      if ( l_galerkin ) then
         !-- Define Galerkin basis and stencils
         if ( ktopv == 1 .and. kbotv == 1 ) then ! Stress-free stress-free
            call get_galerkin_stencil(gal_sten(1), n_r_max, 5)
            call get_galerkin_stencil(gal_sten(2), n_r_max, 9)
         else if ( ktopv /=1 .and. kbotv == 1 ) then
            call get_galerkin_stencil(gal_sten(1), n_r_max, 7)
            call get_galerkin_stencil(gal_sten(2), n_r_max, 11)
         else if ( ktopv ==1 .and. kbotv /= 1 ) then
            call get_galerkin_stencil(gal_sten(1), n_r_max, 6)
            call get_galerkin_stencil(gal_sten(2), n_r_max, 10)
         else if ( ktopv /=1 .and. kbotv /= 1 ) then ! Rigid/Rigid
            call get_galerkin_stencil(gal_sten(1), n_r_max, 1)
            if ( l_non_rot ) then
               call get_galerkin_stencil(gal_sten(2), n_r_max, 8)
            else
               call get_galerkin_stencil(gal_sten(2), n_r_max, 12)
            end if
         else
            call abortRun('This bc is not implemented yet')
         end if
      end if

      !-- Allocate array of matrices (this is handy since it can store various bandwidths)
      allocate( RHSI_mat(nMstart:nMstop) )
      allocate( RHSIL_mat(nMstart:nMstop) )
      if ( l_galerkin ) then
         allocate( LHS_mat_gal(nMstart:nMstop) )
      else
         allocate( LHS_mat_tau(nMstart:nMstop) )
      end if

      !-- Initialize matrices
      if ( l_non_rot ) then

         call RHSE_mat(1)%initialize(4, 4, n_r_max) ! This is m  = 0
         call RHSE_mat(2)%initialize(8, 8, n_r_max) ! This is m /= 0
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m == 0 ) then
               call RHSI_mat(n_m)%initialize(4, 4, n_r_max)
               call RHSIL_mat(n_m)%initialize(2, 2, n_r_max)
               if ( .not. l_galerkin ) then
                  call LHS_mat_tau(n_m)%initialize(4, 4, 2, n_r_max)
               end if
            else
               call RHSI_mat(n_m)%initialize(6, 6, n_r_max)
               call RHSIL_mat(n_m)%initialize(4, 4, n_r_max)
               if ( .not. l_galerkin ) then
                  call LHS_mat_tau(n_m)%initialize(6, 6, 4, n_r_max)
               end if
            end if
         end do

      else ! if this is rotating

         if ( l_ek_pump ) then

            call RHSE_mat(1)%initialize(4, 4, n_r_max) ! This is m  = 0
            call RHSE_mat(2)%initialize(8, 8, n_r_max) ! This is m /= 0
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               if ( m == 0 ) then
                  call RHSI_mat(n_m)%initialize(6, 6, n_r_max)
                  call RHSIL_mat(n_m)%initialize(4, 4, n_r_max)
                  if ( .not. l_galerkin ) then
                     call LHS_mat_tau(n_m)%initialize(6, 6, 2, n_r_max)
                  end if
               else
                  call RHSI_mat(n_m)%initialize(10, 10, n_r_max)
                  if ( .not. l_galerkin ) then
                     call LHS_mat_tau(n_m)%initialize(10, 10, 4, n_r_max)
                  end if
                  if ( l_coriolis_imp ) then
                     call RHSIL_mat(n_m)%initialize(10, 10, n_r_max)
                  else
                     call RHSIL_mat(n_m)%initialize(8, 8, n_r_max)
                  end if
               end if
            end do


         else ! if there's no Ekman pumping

            call RHSE_mat(1)%initialize(4, 4, n_r_max) ! This is m  = 0
            call RHSE_mat(2)%initialize(8, 8, n_r_max) ! This is m /= 0
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               if ( m == 0 ) then
                  call RHSI_mat(n_m)%initialize(4, 4, n_r_max)
                  call RHSIL_mat(n_m)%initialize(2, 2, n_r_max)
                  if ( .not. l_galerkin ) then
                     call LHS_mat_tau(n_m)%initialize(4, 4, 2, n_r_max)
                  end if
               else
                  call RHSI_mat(n_m)%initialize(8, 8, n_r_max)
                  if ( .not. l_galerkin ) then
                     call LHS_mat_tau(n_m)%initialize(8, 8, 4, n_r_max)
                  end if
                  if ( l_coriolis_imp ) then
                     call RHSIL_mat(n_m)%initialize(8, 8, n_r_max)
                  else
                     call RHSIL_mat(n_m)%initialize(6, 6, n_r_max)
                  end if
               end if
            end do

         end if
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

      allocate( psifac(n_r_max, nMstart:nMstop) )
      bytes_allocated = bytes_allocated + n_r_max*(nMstop-nMstart+1)* &
      &                 SIZEOF_DEF_REAL

   end subroutine initialize_psi_integ_smat
!------------------------------------------------------------------------------
   subroutine finalize_psi_integ_smat
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
         if ( l_galerkin ) then
            call LHS_mat_gal(n_m)%finalize()
         else
            call LHS_mat_tau(n_m)%finalize()
         end if
      end do
      if ( l_galerkin ) then
         call destroy_galerkin_stencil(gal_sten(1))
         call destroy_galerkin_stencil(gal_sten(2))
         deallocate( LHS_mat_gal )
      else
         deallocate( LHS_mat_tau )
      end if
      deallocate( RHSIL_mat, RHSI_mat )
      deallocate( psifac, rhs, lPsimat )

   end subroutine finalize_psi_integ_smat
!------------------------------------------------------------------------------
   subroutine update_psi_int_smat(psi_Mloc, om_Mloc, us_Mloc, up_Mloc,         &
              &                   buo_Mloc, dpsidt, vp_bal, tscheme, lMat,     &
              &                   l_vphi_bal_calc, time_solve, n_solve_calls,  &
              &                   time_lu, n_lu_calls, time_dct, n_dct_calls)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lMat
      logical,             intent(in) :: l_vphi_bal_calc

      !-- Output variables
      complex(cp),       intent(out) :: psi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),       intent(out) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),       intent(out) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),       intent(out) :: up_Mloc(nMstart:nMstop,n_r_max)
      type(vp_bal_type), intent(inout) :: vp_bal
      type(type_tarray), intent(inout) :: dpsidt
      complex(cp),       intent(inout) :: buo_Mloc(nMstart:nMstop,n_r_max)
      real(cp),          intent(inout) :: time_solve
      integer,           intent(inout) :: n_solve_calls
      real(cp),          intent(inout) :: time_lu
      integer,           intent(inout) :: n_lu_calls
      real(cp),          intent(inout) :: time_dct
      integer,           intent(inout) :: n_dct_calls

      !-- Local variables
      real(cp) :: uphi0(n_r_max), om0(n_r_max)
      real(cp) :: h2, runStart, runStop
      integer :: n_r, n_m, n_cheb, m

      if ( lMat ) lPsimat(:)=.false.

      !-- If Ekman pumping is needed then regularisation is different
      !-- Hence buoyancy has to be multiplied by h^2
      if ( l_buo_imp ) then
         if ( l_ek_pump ) then
            do n_r=1,n_r_max
               h2 = r_cmb*r_cmb-r(n_r)*r(n_r)
               do n_m=nMstart,nMstop
                  buo_Mloc(n_m,n_r)=h2*buo_Mloc(n_m,n_r)
               end do
            end do
         end if

         !-- Transform buoyancy to Chebyshev space
         runStart = MPI_Wtime()
         call rscheme%costf1(buo_Mloc, nMstart, nMstop, n_r_max)
         runStop = MPI_Wtime()
         if ( runStop > runStart ) then
            time_dct = time_dct + (runStop-runStart)
            n_dct_calls = n_dct_calls + 1
         end if

         !-- Matrix-vector multiplication by the operator \int\int\int\int r^4 .
         do n_m=nMstart,nMstop
            m = idx2m(n_m)

            if ( m > 0 ) then
               do n_cheb=1,n_r_max
                  rhs(n_cheb)=buo_Mloc(n_m,n_cheb)
               end do

               call RHSE_mat(2)%mat_vec_mul(rhs)

               rhs(1)=zero
               rhs(2)=zero
               rhs(3)=zero
               rhs(4)=zero
               do n_cheb=1,n_r_max
                  buo_Mloc(n_m,n_cheb)=rhs(n_cheb)
               end do

            end if
         end do
      end if

      !-- Calculation of the implicit part
      call get_psi_rhs_imp_int_smat(psi_Mloc, up_Mloc,                       &
           &                        dpsidt%old(:,:,tscheme%istage),          &
           &                        dpsidt%impl(:,:,tscheme%istage), vp_bal, &
           &                        l_vphi_bal_calc,                         &
           &                        tscheme%l_imp_calc_rhs(tscheme%istage))


      !-- Now assemble the right hand side and store it in work_Mloc
      call tscheme%set_imex_rhs(work_Mloc, dpsidt, nMstart, nMstop, n_r_max)

      do n_m=nMstart,nMstop

         m = idx2m(n_m)

         if ( .not. lPsimat(n_m) ) then
            if ( l_galerkin ) then
               call get_lhs_mat_gal(tscheme, LHS_mat_gal(n_m), psifac(:, n_m), m, &
                    &               time_lu, n_lu_calls)
            else
               call get_lhs_mat_tau(tscheme, LHS_mat_tau(n_m), psifac(:, n_m), m, &
                    &               time_lu, n_lu_calls)
            end if
            lPsimat(n_m)=.true.
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
            rhs(n_cheb)=rhs(n_cheb)*psifac(n_cheb,n_m)
         end do

         runStart = MPI_Wtime()
         if ( l_galerkin ) then
            if ( m == 0 ) then
               call LHS_mat_gal(n_m)%solve(rhs(3:n_r_max), n_r_max-2)
               !-- Put the two first zero lines at the end
               rhs = cshift(rhs, 2)
               !-- Transform from Galerkin space to Chebyshev space
               call galerkin2cheb(gal_sten(1), rhs)
            else
               call LHS_mat_gal(n_m)%solve(rhs(5:n_r_max), n_r_max-4)
               !-- Put the four first zero lines at the end
               rhs = cshift(rhs, 4)
               !-- Transform from Galerkin space to Chebyshev space
               call galerkin2cheb(gal_sten(2), rhs)
            end if
         else
            call LHS_mat_tau(n_m)%solve(rhs,n_r_max)
         end if
         runStop = MPI_Wtime()
         if ( runStop > runStart ) then
            time_solve = time_solve + (runStop-runStart)
            n_solve_calls = n_solve_calls+1
         end if

         if ( m == 0 ) then
            do n_cheb=1,n_r_max
               uphi0(n_cheb)=real(rhs(n_cheb))
            end do
         else
            do n_cheb=1,n_r_max
               psi_Mloc(n_m,n_cheb)=rhs(n_cheb)
            end do
         end if

      end do

      !-- set cheb modes > n_cheb_max to zero (dealiazing)
      !if ( n_cheb_max < n_r_max ) then ! fill with zeros !
      !   do n_cheb=n_cheb_max+1,n_r_max
      !      do n_m=nMstart,nMstop
      !         m = idx2m(n_m)
      !         if ( m == 0 ) then
      !            uphi0(n_cheb)=0.0_cp
      !         else
      !            psi_Mloc(n_m,n_cheb)=zero
      !         end if
      !      end do
      !   end do
      !end if

      !-- Bring uphi0 to the physical space
      if ( l_rank_has_m0 ) then
         call get_dr(uphi0, om0, n_r_max, rscheme, l_dct=.false.)
         call rscheme%costf1(uphi0, n_r_max)
      end if

      !-- Get the radial derivative of psi to calculate uphi, us and omega
      call get_ddr(psi_Mloc, work_Mloc, om_Mloc, nMstart, nMstop, n_r_max, &
           &       rscheme,l_dct=.false.)

      !-- Bring psi back the physical space
      runStart = MPI_Wtime()
      call rscheme%costf1(psi_Mloc, nMstart, nMstop, n_r_max)
      runStop = MPI_Wtime()
      if ( runStop > runStart ) then
         time_dct = time_dct + (runStop-runStart)
         n_dct_calls = n_dct_calls + 1
      end if

      if ( l_non_rot ) then
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)

               if ( m == 0 ) then
                  us_Mloc(n_m,n_r)=0.0_cp
                  up_Mloc(n_m,n_r)=uphi0(n_r)
                  om_Mloc(n_m,n_r)=om0(n_r)+or1(n_r)*uphi0(n_r)
               else
                  us_Mloc(n_m,n_r)=ci*real(m,cp)*or1(n_r)*psi_Mloc(n_m,n_r)
                  up_Mloc(n_m,n_r)=-work_Mloc(n_m,n_r)
                  om_Mloc(n_m,n_r)=-om_Mloc(n_m,n_r)-or1(n_r)*work_Mloc(n_m,n_r)+ &
                  &                real(m,cp)*real(m,cp)*or2(n_r)*psi_Mloc(n_m,n_r)
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
               else
                  us_Mloc(n_m,n_r)=ci*real(m,cp)*or1(n_r)*h2*   psi_Mloc(n_m,n_r)
                  up_Mloc(n_m,n_r)=-h2*                        work_Mloc(n_m,n_r) &
                  &                +3.0_cp*r(n_r)*              psi_Mloc(n_m,n_r)
                  om_Mloc(n_m,n_r)=-h2*                          om_Mloc(n_m,n_r) &
                  &     -(r_cmb*r_cmb*or1(n_r)-6.0_cp*r(n_r))* work_Mloc(n_m,n_r) &
                  &     +(6.0_cp+real(m,cp)*real(m,cp)*or2(n_r)*h2)*              &
                  &                                             psi_Mloc(n_m,n_r)
               end if
            end do
         end do
      end if

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dpsidt, nMstart, nMstop, n_r_max)

   end subroutine update_psi_int_smat
!------------------------------------------------------------------------------
   subroutine finish_exp_psi_int_smat(psi_Mloc, us_Mloc, up_Mloc, om_Mloc, &
              &                       dVsOm_Mloc, buo_Mloc, dpsi_exp_last, &
              &                       vp_bal, l_vphi_bal_calc)

      !-- Input variables
      complex(cp), intent(in) :: psi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: dVsOm_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: buo_Mloc(nMstart:nMstop,n_r_max)
      logical,     intent(in) :: l_vphi_bal_calc

      !-- Output variables
      type(vp_bal_type), intent(inout) :: vp_bal
      complex(cp),       intent(inout) :: dpsi_exp_last(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_m, m, n_cheb
      real(cp) :: h2, ekp_fac

      call rscheme%costf1( dVsOm_Mloc, nMstart, nMstop, n_r_max )
      do n_cheb=n_cheb_max+1,n_r_max
         do n_m=nMstart,nMstop
            dVsOm_Mloc(n_m,n_cheb)=zero
         end do
      end do
      !-- Finish calculation of advection
      call get_dr( dVsOm_Mloc, work_Mloc, nMstart, nMstop, n_r_max, &
           &       rscheme, nocopy=.true., l_dct=.false.)

      !-- Finish calculation of the explicit part for current time step
      do n_r=1,n_r_max
         do n_m=nMstart, nMstop
            m = idx2m(n_m)
            if ( m /= 0 ) then
               dpsi_exp_last(n_m,n_r)=dpsi_exp_last(n_m,n_r)-    &
               &                      or1(n_r)*work_Mloc(n_m,n_r)

               !-- If Coriolis force is treated explicitly it is added here:
               if ( .not. l_coriolis_imp ) then
                  dpsi_exp_last(n_m,n_r)=dpsi_exp_last(n_m,n_r)-   &
                  &                      CorFac*ci*real(m,cp)*psi_Mloc(n_m,n_r)
               end if

               !-- If Buoyancy is treated explicitly
               if ( .not. l_buo_imp ) then 
                  dpsi_exp_last(n_m,n_r)=dpsi_exp_last(n_m,n_r)+buo_Mloc(n_m,n_r)
               end if
            end if
         end do
      end do

      if ( l_rank_has_m0 .and. l_vphi_bal_calc ) then
         do n_r=1,n_r_max
            vp_bal%rey_stress(n_r)=real(dpsi_exp_last(m2idx(0),n_r))
         end do
      end if

      !-- Add Ekman pumping as an explicit term if this is requested
      if ( l_ek_pump ) then
         do n_r=1,n_r_max
            !-- The following statement is required to make sure
            !-- that h2 remains 0+ even when -O3 is used at the compilation
            if ( n_r == 1 ) then
               h2 = 0.0_cp
            else
               h2 = (r_cmb*r_cmb-r(n_r)*r(n_r))
            end if
            ekp_fac = CorFac*half*sqrt(ek*r_cmb)*h2**(0.25_cp)
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               if ( m == 0 ) then
                  dpsi_exp_last(n_m,n_r)=h2*dpsi_exp_last(n_m,n_r)- &
                  &                       ekp_fac*up_Mloc(n_m,n_r)
               else 
                  dpsi_exp_last(n_m,n_r)=  h2*dpsi_exp_last(n_m,n_r) +       &
                  &                      ekp_fac*( -om_Mloc(n_m,n_r) +       &
                  &               half*beta(n_r)*   up_Mloc(n_m,n_r) +       &
                  &   beta(n_r)*(-ci*real(m,cp)+5.0_cp*r_cmb*oheight(n_r))*  &
                  &                                 us_Mloc(n_m,n_r) )
               end if
            end do
         end do
      end if

      !-- Transform the explicit part to chebyshev space
      call rscheme%costf1(dpsi_exp_last, nMstart, nMstop, n_r_max)

      do n_cheb=n_cheb_max+1,n_r_max
         do n_m=nMstart,nMstop
            dpsi_exp_last(n_m,n_cheb)=zero
         end do
      end do

      !-- Matrix-vector multiplication by the operator \int\int\int\int r^4 .
      do n_m=nMstart,nMstop
         m = idx2m(n_m)

         do n_cheb=1,n_r_max
            rhs(n_cheb)=dpsi_exp_last(n_m,n_cheb)
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

         do n_cheb=1,n_r_max
            dpsi_exp_last(n_m,n_cheb)=rhs(n_cheb)
         end do

      end do

   end subroutine finish_exp_psi_int_smat
!------------------------------------------------------------------------------
   subroutine get_psi_rhs_imp_int_smat(psi_Mloc, up_Mloc, psi_old,   &
              &                        dpsi_imp_Mloc_last, vp_bal,   &
              &                        l_vphi_bal_calc, l_calc_lin_rhs)

      !-- Input variables
      complex(cp), intent(in) :: psi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      logical,     intent(in) :: l_vphi_bal_calc
      logical,     intent(in) :: l_calc_lin_rhs

      !-- Output variables
      complex(cp),       intent(inout) :: dpsi_imp_Mloc_last(nMstart:nMstop,n_r_max)
      complex(cp),       intent(inout) :: psi_old(nMstart:nMstop,n_r_max)
      type(vp_bal_type), intent(inout) :: vp_bal

      !-- Local variables
      real(cp) :: duphi0(n_r_max), d2uphi0(n_r_max), uphi0(n_r_max)
      integer :: n_r, n_m, m, m0, n_cheb

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

      !-- Transform the implicit part to chebyshev space
      call rscheme%costf1(work_Mloc, nMstart, nMstop, n_r_max)

      !-- Matrix-vector multiplication by the operator -\int^4 r^4 \Delta .
      do n_m=nMstart,nMstop
         m = idx2m(n_m)

         do n_cheb=1,n_r_max
            rhs(n_cheb)= work_Mloc(n_m,n_cheb)
         end do

         call RHSI_mat(n_m)%mat_vec_mul(rhs)

         rhs(1)=zero ! vphi equation has only 2 BCs
         rhs(2)=zero
         if ( m > 0 ) then
            rhs(3)=zero
            rhs(4)=zero
         end if

         do n_cheb=1,n_r_max
            psi_old(n_m,n_cheb)=rhs(n_cheb)
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
            vp_bal%visc(n_r)=ViscFac*(d2uphi0(n_r)+or1(n_r)*duphi0(n_r)-&
            &                or2(n_r)*uphi0(n_r))
            vp_bal%pump(n_r)=-CorFac*ekpump(n_r)*uphi0(n_r)
         end do
      end if

      if ( l_calc_lin_rhs ) then

         !-- Matrix-vector multiplication by the LHS operator
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            do n_cheb=1,n_r_max
               rhs(n_cheb)=work_Mloc(n_m,n_cheb)
            end do
            call RHSIL_mat(n_m)%mat_vec_mul(rhs)
            rhs(1)=zero ! vphi equation has only 2 BCs
            rhs(2)=zero
            if ( m > 0 ) then
               rhs(3)=zero
               rhs(4)=zero
            end if
            do n_cheb=1,n_r_max
               work_Mloc(n_m,n_cheb)=rhs(n_cheb)
            end do
         end do

         !-- Finally assemble the right hand side
         do n_cheb=1,n_r_max
            do n_m=nMstart,nMstop
               dpsi_imp_Mloc_last(n_m,n_cheb)=work_Mloc(n_m,n_cheb)
            end do
         end do

      end if 

   end subroutine get_psi_rhs_imp_int_smat
!------------------------------------------------------------------------------
   subroutine get_lhs_mat_gal(tscheme, Cmat, psiMat_fac, m, time_lu, n_lu_calls)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme    ! time step
      integer,             intent(in) :: m          ! Azimuthal wavenumber

      !-- Output variables
      type(type_bandmat_complex), intent(out) :: Cmat
      real(cp),                   intent(inout) :: psiMat_fac(:)
      real(cp),                   intent(inout) :: time_lu
      integer,                    intent(inout) :: n_lu_calls

      !-- Local variables
      type(type_bandmat_complex) :: Amat
      real(cp), allocatable :: stencilA(:), CorSten(:)
      integer :: n_r, i_r, n_band, klA, kuA, n_boundaries
      real(cp) :: a, b, runStart, runStop

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      if ( m == 0 ) then
         n_boundaries = 2
      else
         n_boundaries = 4
      end if

      if ( l_non_rot ) then
         if ( m == 0 ) then
            klA = 4
            kuA = 4
         else
            klA = 6
            kuA = 6
         end if
      else
         if ( l_ek_pump ) then
            if ( m == 0 ) then
               klA = 6
               kuA = 6
            else
               klA = 10
               kuA = 10
            end if
         else
            if ( m == 0 ) then
               klA = 4
               kuA = 4
            else
               klA = 8
               kuA = 8
            end if
         end if
      end if

      call Amat%initialize(klA, kuA, n_r_max)
      allocate( stencilA(Amat%nbands), CorSten(Amat%nbands) )


      !-- Fill A matrix
      do n_r=1,Amat%nlines
         i_r = n_r+n_boundaries

         !-- Define the equations
         if ( l_non_rot ) then

            if ( m == 0 ) then
               stencilA = intcheb2rmult2(a,b,i_r-1,Amat%nbands)-               &
               &  tscheme%wimp_lin(1)*ViscFac*( rmult2(a,b,i_r-1,Amat%nbands)- &
               &                 3.0_cp*intcheb1rmult1(a,b,i_r-1,Amat%nbands) )
               CorSten   = 0.0_cp
            else
               stencilA = -intcheb4rmult4lapl(a,b,m,i_r-1,Amat%nbands)+  &
               &          tscheme%wimp_lin(1)*ViscFac*                    &
               &          intcheb4rmult4lapl2(a,b,m,i_r-1,Amat%nbands)  
               CorSten   = 0.0_cp
            end if

         else ! this is rotating

            if ( l_ek_pump ) then

               if ( m == 0 ) then
                  stencilA = intcheb2rmult2hmult2(a,b,i_r-1,Amat%nbands)-      &
                  &          tscheme%wimp_lin(1)*ViscFac*                       &
                  &          intcheb2rmult2hmult2laplm1(a,b,i_r-1,Amat%nbands)
                  CorSten  = 0.0_cp
               else
                  stencilA = intcheb4rmult4hmult2laplrot(a,b,m,i_r-1,Amat%nbands)&
                  &          -tscheme%wimp_lin(1)*ViscFac*                        &
                  &          intcheb4rmult4hmult2laplrot2(a,b,m,i_r-1,Amat%nbands)  
                  if ( l_coriolis_imp ) then
                     CorSten   = tscheme%wimp_lin(1)*CorFac*real(m,cp)*        &
                     &           intcheb4rmult4hmult2(a,b,i_r-1,Amat%nbands)
                  else
                     CorSten   = 0.0_cp
                  end if
               end if

            else ! there's no Ekman pumping

               if ( m == 0 ) then
                  stencilA = intcheb2rmult2(a,b,i_r-1,Amat%nbands)-              &
                  &  tscheme%wimp_lin(1)*ViscFac*( rmult2(a,b,i_r-1,Amat%nbands)-&
                  &                 3.0_cp*intcheb1rmult1(a,b,i_r-1,Amat%nbands) )
                  CorSten   = 0.0_cp
               else
                  stencilA = intcheb4rmult4laplrot(a,b,m,i_r-1,Amat%nbands)   &
                  &          -tscheme%wimp_lin(1)*ViscFac*                     &
                  &          intcheb4rmult4laplrot2(a,b,m,i_r-1,Amat%nbands)  
                  if ( l_coriolis_imp ) then
                     CorSten   = tscheme%wimp_lin(1)*CorFac*real(m,cp)*     &
                     &           intcheb4rmult4(a,b,i_r-1,Amat%nbands)
                  else
                     CorSten   = 0.0_cp
                  end if
               end if

            end if

         end if

         !-- Roll the array for band storage
         do n_band=1,Amat%nbands
            if ( i_r+Amat%ku+1-n_band <= Amat%nlines .and. i_r+Amat%ku+1-n_band >= 1 ) then
               Amat%dat(n_band,i_r+Amat%ku+1-n_band) = rscheme%rnorm* &
               &             cmplx(stencilA(n_band), CorSten(n_band), kind=cp)
            end if
         end do

      end do

      if ( m == 0 ) then
         call band_band_product(Amat, gal_sten(1), Cmat, l_lhs=.true.)

         !-- Remove first blank rows (careful remapping of kl, ku and room for LU factorisation)
         call Cmat%remove_leading_blank_rows(2)

         !-- Truncate the last N lines
         call Cmat%remove_last_rows(2)
      else
         call band_band_product(Amat, gal_sten(2), Cmat, l_lhs=.true.)

         !-- Remove first blank rows (careful remapping of kl, ku and room for LU factorisation)
         call Cmat%remove_leading_blank_rows(4)

         !-- Truncate the last N lines
         call Cmat%remove_last_rows(4)
      end if

      deallocate( stencilA, CorSten )
      call Amat%finalize()

      !-- To be improved
      psiMat_fac(:)=one


      !do n_r=1,A_mat%nlines
      !   i_r = n_r+A_mat%ntau
!
!         do n_band=1,A_mat%nbands
!            if ( n_r+A_mat%ku+1-n_band <= A_mat%nlines .and. &
!            &                         n_r+A_mat%ku+1-n_band >= 1 ) then
!               A_mat%A4(A_mat%kl+n_band,n_r+A_mat%ku+1-n_band) =  &
!               & A_mat%A4(A_mat%kl+n_band,n_r+A_mat%ku+1-n_band)* &
!               &                      psiMat_fac(n_r+A_mat%ntau)
!            end if
!         end do

 !     end do

      ! if ( m == 10 ) call A_mat%write()

      !-- LU factorisation
      runStart = MPI_Wtime()
      call Cmat%prepare_LU()
      runStop = MPI_Wtime()
      if ( runStop > runStart ) then
         time_lu = time_lu + (runStop-runStart)
         n_lu_calls = n_lu_calls + 1
      end if

   end subroutine get_lhs_mat_gal
!------------------------------------------------------------------------------
   subroutine get_lhs_mat_tau(tscheme, A_mat, psiMat_fac, m, time_lu, n_lu_calls)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme    ! time step
      integer,             intent(in) :: m          ! Azimuthal wavenumber

      !-- Output variables
      type(type_bordmat_complex), intent(inout) :: A_mat
      real(cp),                   intent(inout) :: psiMat_fac(A_mat%nlines)
      real(cp),                   intent(inout) :: time_lu
      integer,                    intent(inout) :: n_lu_calls

      !-- Local variables
      real(cp) :: stencilA4(A_mat%nbands), CorSten(A_mat%nbands)
      integer :: n_r, i_r, n_band, n_b
      real(cp) :: a, b, dn2, runStart, runStop
      real(cp) :: d4top(A_mat%nlines)

      do n_r=1,A_mat%nlines
         dn2 = real(n_r-1,cp)*real(n_r-1,cp)
         d4top(n_r)=16.0_cp/105.0_cp*dn2*(dn2-one)*(dn2-4.0_cp)*(dn2-9.0_cp)
      end do


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
               &   tscheme%wimp_lin(1)*ViscFac*( rmult2(a,b,i_r-1,A_mat%nbands)- &
               &                  3.0_cp*intcheb1rmult1(a,b,i_r-1,A_mat%nbands) )
               CorSten   = 0.0_cp
            else
               stencilA4 = -intcheb4rmult4lapl(a,b,m,i_r-1,A_mat%nbands)+  &
               &           tscheme%wimp_lin(1)*ViscFac*                    &
               &           intcheb4rmult4lapl2(a,b,m,i_r-1,A_mat%nbands)  
               CorSten   = 0.0_cp
            end if

         else ! this is rotating

            if ( l_ek_pump ) then

               if ( m == 0 ) then
                  stencilA4 = intcheb2rmult2hmult2(a,b,i_r-1,A_mat%nbands)-      &
                  &           tscheme%wimp_lin(1)*ViscFac*                       &
                  &           intcheb2rmult2hmult2laplm1(a,b,i_r-1,A_mat%nbands)
                  CorSten   = 0.0_cp
               else
                  stencilA4 = intcheb4rmult4hmult2laplrot(a,b,m,i_r-1,A_mat%nbands)&
                  &           -tscheme%wimp_lin(1)*ViscFac*                        &
                  &           intcheb4rmult4hmult2laplrot2(a,b,m,i_r-1,A_mat%nbands)  
                  if ( l_coriolis_imp ) then
                     CorSten   = tscheme%wimp_lin(1)*CorFac*real(m,cp)*        &
                     &           intcheb4rmult4hmult2(a,b,i_r-1,A_mat%nbands)
                  else
                     CorSten   = 0.0_cp
                  end if
               end if

            else ! there's no Ekman pumping

               if ( m == 0 ) then
                  stencilA4 = intcheb2rmult2(a,b,i_r-1,A_mat%nbands)-              &
                  &   tscheme%wimp_lin(1)*ViscFac*( rmult2(a,b,i_r-1,A_mat%nbands)-&
                  &                  3.0_cp*intcheb1rmult1(a,b,i_r-1,A_mat%nbands) )
                  CorSten   = 0.0_cp
               else
                  stencilA4 = intcheb4rmult4laplrot(a,b,m,i_r-1,A_mat%nbands)   &
                  &           -tscheme%wimp_lin(1)*ViscFac*                     &
                  &           intcheb4rmult4laplrot2(a,b,m,i_r-1,A_mat%nbands)  
                  if ( l_coriolis_imp ) then
                     CorSten   = tscheme%wimp_lin(1)*CorFac*real(m,cp)*     &
                     &           intcheb4rmult4(a,b,i_r-1,A_mat%nbands)
                  else
                     CorSten   = 0.0_cp
                  end if
               end if

            end if

         end if

         !-- Roll the array for band storage
         do n_band=1,A_mat%nbands
            if ( n_r+A_mat%ku+1-n_band <= A_mat%nlines_band .and. n_r+A_mat%ku+1-n_band >= 1 ) then
               A_mat%A4(A_mat%kl+n_band,n_r+A_mat%ku+1-n_band) = rscheme%rnorm* &
               &             cmplx(stencilA4(n_band), CorSten(n_band), kind=cp)
            end if
         end do

         psiMat_fac(n_r+A_mat%ntau)=one/rscheme%rnorm/  &
         &                          maxval(abs(cmplx(stencilA4, CorSten, kind=cp)))
      end do

      !-- Fill A3
      do n_r=1,A_mat%nlines_band
         i_r = n_r+A_mat%ntau

         if ( l_non_rot ) then

            if ( m == 0 ) then
               stencilA4 = intcheb2rmult2(a,b,i_r-1,A_mat%nbands)-               &
               &   tscheme%wimp_lin(1)*ViscFac*( rmult2(a,b,i_r-1,A_mat%nbands)- &
               &                  3.0_cp*intcheb1rmult1(a,b,i_r-1,A_mat%nbands) )
               CorSten   = 0.0_cp
            else
               stencilA4 = -intcheb4rmult4lapl(a,b,m,i_r-1,A_mat%nbands)+  &
               &           tscheme%wimp_lin(1)*ViscFac*                    &
               &           intcheb4rmult4lapl2(a,b,m,i_r-1,A_mat%nbands)  
               CorSten   = 0.0_cp
            end if

         else ! this is rotating

            if ( l_ek_pump ) then 

               if ( m == 0 ) then
                  stencilA4 = intcheb2rmult2hmult2(a,b,i_r-1,A_mat%nbands)-    &
                  &           tscheme%wimp_lin(1)*ViscFac*                     &
                  &           intcheb2rmult2hmult2laplm1(a,b,i_r-1,A_mat%nbands)
                  CorSten   = 0.0_cp
               else
                  stencilA4 = intcheb4rmult4hmult2laplrot(a,b,m,i_r-1,A_mat%nbands)&
                  &           -tscheme%wimp_lin(1)*ViscFac*                        &
                  &           intcheb4rmult4hmult2laplrot2(a,b,m,i_r-1,A_mat%nbands)  
                  if ( l_coriolis_imp ) then
                     CorSten = tscheme%wimp_lin(1)*CorFac*real(m,cp)*        &
                     &         intcheb4rmult4hmult2(a,b,i_r-1,A_mat%nbands)
                  else
                     CorSten = 0.0_cp
                  end if
               end if

            else ! there's no Ekman pumping

               if ( m == 0 ) then
                  stencilA4 = intcheb2rmult2(a,b,i_r-1,A_mat%nbands)-              &
                  &   tscheme%wimp_lin(1)*ViscFac*( rmult2(a,b,i_r-1,A_mat%nbands)-&
                  &                  3.0_cp*intcheb1rmult1(a,b,i_r-1,A_mat%nbands) )
                  CorSten   = 0.0_cp
               else
                  stencilA4 = intcheb4rmult4laplrot(a,b,m,i_r-1,A_mat%nbands)   &
                  &           -tscheme%wimp_lin(1)*ViscFac*                     &
                  &           intcheb4rmult4laplrot2(a,b,m,i_r-1,A_mat%nbands)  
                  if ( l_coriolis_imp ) then
                     CorSten = tscheme%wimp_lin(1)*CorFac*real(m,cp)*        &
                     &         intcheb4rmult4(a,b,i_r-1,A_mat%nbands)
                  else
                     CorSten = 0.0_cp
                  end if
               end if

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
                  if ( l_non_rot ) then
                     A_mat%A1(3,n_r)=rscheme%rnorm*rscheme%drMat(1,n_r)
                  else
                     A_mat%A1(3,n_r)=rscheme%rnorm*d4top(n_r)
                  end if
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
                  if ( l_non_rot ) then
                     A_mat%A2(3,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%drMat(1,n_r)
                  else
                     A_mat%A2(3,n_r-A_mat%ntau)=rscheme%rnorm*d4top(n_r)
                  end if
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

      ! if ( m == 10 ) call A_mat%write()

      !-- LU factorisation
      runStart = MPI_Wtime()
      call A_mat%prepare_LU()
      runStop = MPI_Wtime()
      if ( runStop > runStart ) then
         time_lu = time_lu + (runStop-runStart)
         n_lu_calls = n_lu_calls + 1
      end if

   end subroutine get_lhs_mat_tau
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
         if ( m0 == 1) then
            stencilD = intcheb2rmult2(a,b,i_r-1,D_mat%nbands)
         else ! Non-axisymmetric terms
            stencilD = intcheb4rmult4(a,b,i_r-1,D_mat%nbands)
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

         else ! if this is rotating

            if ( l_ek_pump ) then

               if ( m == 0 ) then
                  stencilB = intcheb2rmult2hmult2(a,b,i_r-1,B_mat%nbands)
               else
                  stencilB = intcheb4rmult4hmult2laplrot(a,b,m,i_r-1,B_mat%nbands)
               end if

            else ! if there's no Ekman pumping

               if ( m == 0 ) then
                  stencilB = intcheb2rmult2(a,b,i_r-1,B_mat%nbands)
               else
                  stencilB = intcheb4rmult4laplrot(a,b,m,i_r-1,B_mat%nbands)
               end if

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
               stencilC =  ViscFac * (  rmult2(a,b,i_r-1,Cmat%nbands)- &
               &         3.0_cp*intcheb1rmult1(a,b,i_r-1,Cmat%nbands) )
               CorSten  = 0.0_cp
            else
               stencilC = -ViscFac*intcheb4rmult4lapl2(a,b,m,i_r-1,Cmat%nbands)
               CorSten  = 0.0_cp
            end if

         else ! if this is rotating

            if ( l_ek_pump ) then

               if ( m == 0 ) then
                  stencilC = ViscFac*intcheb2rmult2hmult2laplm1(a,b,i_r-1, &
                  &                                             Cmat%nbands)
                  CorSten  = 0.0_cp
               else
                  stencilC = ViscFac*intcheb4rmult4hmult2laplrot2(a,b,m,i_r-1, &
                  &                                               Cmat%nbands)
                  if ( l_coriolis_imp ) then
                     CorSten = -CorFac*real(m,cp)*intcheb4rmult4hmult2(a,b,i_r-1,&
                     &                                                 Cmat%nbands)
                  else
                     CorSten = 0.0_cp
                  end if
               end if

            else ! if there's no Ekman pumping

               if ( m == 0 ) then
                  stencilC =    ViscFac* ( rmult2(a,b,i_r-1,Cmat%nbands)- &
                  &         3.0_cp*intcheb1rmult1(a,b,i_r-1,Cmat%nbands))
                  CorSten  = 0.0_cp
               else
                  stencilC = ViscFac*intcheb4rmult4laplrot2(a,b,m,i_r-1,Cmat%nbands)
                  if ( l_coriolis_imp ) then
                     CorSten = -CorFac*real(m,cp)*intcheb4rmult4(a,b,i_r-1,&
                     &                                           Cmat%nbands)
                  else
                     CorSten = 0.0_cp
                  end if
               end if

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
end module update_psi_integ_smat
