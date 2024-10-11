module update_phi_coll

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: one, zero, ci
   use namelists, only: kbotphi, ktopphi, phaseDiffFac, l_full_disk, stef, pr, &
       &                phi_top, phi_bot
   use radial_functions, only: rscheme, or1, or2
   use blocking, only: nMstart, nMstop
   use truncation, only: n_r_max, idx2m, m2idx
   use radial_der, only: get_ddr
   use fields, only: work_Mloc
   use algebra, only: prepare_full_mat, solve_full_mat
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray

   implicit none
   
   private

   logical,  allocatable :: lPhimat(:)
   real(cp), allocatable :: phiMat(:, :, :)
   integer,  allocatable :: phiPivot(:, :)
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: phiMat_fac(:,:)
#endif
   complex(cp), allocatable :: rhs(:)

   public :: update_phi_co, initialize_phi_coll, finalize_phi_coll, &
   &         get_phi_rhs_imp_coll, assemble_phi_coll

contains

   subroutine initialize_phi_coll

      allocate( lPhimat(nMstart:nMstop) )
      lPhimat(:)=.false.
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_LOGICAL

      allocate( phiMat(n_r_max, n_r_max, nMstart:nMstop) )
      allocate( phiPivot(n_r_max, nMstart:nMstop) )
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*n_r_max*n_r_max* &
      &                 SIZEOF_DEF_REAL+n_r_max*(nMstop-nMstart+1)*SIZEOF_INTEGER
#ifdef WITH_PRECOND_S
      allocate( phiMat_fac(n_r_max, nMstart:nMstop) )
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*n_r_max*  &
      &                 SIZEOF_DEF_REAL
#endif
      allocate( rhs(n_r_max) )
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_COMPLEX

   end subroutine initialize_phi_coll
!------------------------------------------------------------------------------
   subroutine finalize_phi_coll

      deallocate( rhs )
#ifdef WITH_PRECOND_S
      deallocate( phiMat_fac )
#endif
      deallocate( lPhimat, phiMat, phiPivot )

   end subroutine finalize_phi_coll
!------------------------------------------------------------------------------
   subroutine update_phi_co(phi_Mloc, dphidt, tscheme, lMat)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lMat

      !-- Output variables
      complex(cp),       intent(out) :: phi_Mloc(nMstart:nMstop, n_r_max)
      type(type_tarray), intent(inout) :: dphidt

      !-- Local variables
      integer :: n_m, n_r_out, m

      if ( lMat ) lPhimat(:)=.false.

      !-- Now assemble the right hand side and store it in work_Mloc
      call tscheme%set_imex_rhs(work_Mloc, dphidt)

      do n_m=nMstart, nMstop

         m = idx2m(n_m)
         
         if ( .not. lPhimat(n_m) ) then
#ifdef WITH_PRECOND_S
            call get_phiMat(tscheme, m, phiMat(:,:,n_m), phiPivot(:,n_m), &
                 &          phiMat_fac(:,n_m))
#else
            call get_phiMat(tscheme, m, phiMat(:,:,n_m), phiPivot(:,n_m))
#endif
            lPhimat(n_m)=.true.
         end if

         !-- Inhomogeneous B.Cs (if not zero)
         if ( m == 0 ) then
            rhs(1)      =cmplx(phi_top,0.0_cp,cp)
            rhs(n_r_max)=cmplx(phi_bot,0.0_cp,cp)
         else
            rhs(1)      =zero
            rhs(n_r_max)=zero
         end if
         rhs(2:n_r_max-1)=work_Mloc(n_m,2:n_r_max-1)

#ifdef WITH_PRECOND_S
         rhs(:) = phiMat_fac(:,n_m)*rhs(:)
#endif

         call solve_full_mat(phiMat(:,:,n_m), n_r_max, n_r_max, phiPivot(:, n_m), &
              &              rhs(:))

         do n_r_out=1,rscheme%n_max
            phi_Mloc(n_m,n_r_out)=rhs(n_r_out)
         end do

      end do

      !-- set cheb modes > rscheme%n_max to zero (dealiazing)
      do n_r_out=rscheme%n_max+1,n_r_max
         phi_Mloc(:,n_r_out)=zero
      end do

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dphidt)

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_phi_rhs_imp_coll(phi_Mloc, dphidt, 1, tscheme%l_imp_calc_rhs(1), &
              &                    l_in_cheb_space=.true.)
      else
         call get_phi_rhs_imp_coll(phi_Mloc, dphidt, tscheme%istage+1,       &
              &                    tscheme%l_imp_calc_rhs(tscheme%istage+1), &
              &                    l_in_cheb_space=.true.)
      end if

   end subroutine update_phi_co
!------------------------------------------------------------------------------
   subroutine get_phi_rhs_imp_coll(phi_Mloc, dphidt, istage, l_calc_lin, &
              &                    l_in_cheb_space)

      !-- Input variables
      complex(cp),       intent(inout) :: phi_Mloc(nMstart:nMstop,n_r_max)
      integer,           intent(in) :: istage
      logical,           intent(in) :: l_calc_lin
      logical, optional, intent(in) :: l_in_cheb_space

      !-- Output variable
      type(type_tarray), intent(inout) :: dphidt

      !-- Local variables
      complex(cp) :: dphi_Mloc(nMstart:nMstop,n_r_max)
      logical :: l_in_cheb
      integer :: n_r, n_m, m
      real(cp) :: dm2

      if ( present(l_in_cheb_space) ) then
         l_in_cheb = l_in_cheb_space
      else
         l_in_cheb = .false.
      end if

      if ( l_calc_lin ) then
         call get_ddr(phi_Mloc, dphi_Mloc, work_Mloc, nMstart, nMstop, &
              &       n_r_max, rscheme, l_dct=.not. l_in_cheb)
         if ( l_in_cheb ) call rscheme%costf1(phi_Mloc, nMstart, nMstop, n_r_max)
      else
         if ( l_in_cheb ) call rscheme%costf1(phi_Mloc, nMstart, nMstop, n_r_max)
      end if

      if ( istage == 1 ) then
         do n_r=1,n_r_max
            dphidt%old(:,n_r,istage)=5.0_cp/6.0_cp*stef*pr*phi_Mloc(:,n_r)
         end do
      end if

      if ( l_calc_lin ) then
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               dm2 = real(m,cp)*real(m,cp)
               dphidt%impl(n_m,n_r,istage)=phaseDiffFac* (                &
               &                                      work_Mloc(n_m,n_r)  &
               &                         +or1(n_r)*   dphi_Mloc(n_m,n_r)  &
               &                     -dm2*or2(n_r)*    phi_Mloc(n_m,n_r) )
            end do
         end do
      end if

   end subroutine get_phi_rhs_imp_coll
!------------------------------------------------------------------------------
   subroutine assemble_phi_coll(phi_Mloc, dphidt, tscheme)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      type(type_tarray),   intent(inout) :: dphidt

      !-- Output variable
      complex(cp), intent(inout) :: phi_Mloc(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_m, m

      call tscheme%assemble_imex(work_Mloc, dphidt)

      do n_r=2,n_r_max-1
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m == 0 ) then
               phi_Mloc(n_m,n_r)=cmplx(real(work_Mloc(n_m,n_r)), 0.0_cp, cp)* &
               &                 6.0_cp/5.0_cp/stef/pr
            else
               phi_Mloc(n_m,n_r)=work_Mloc(n_m,n_r)*6.0_cp/5.0_cp/stef/pr
            end if
         end do
      end do

      if ( ktopphi==1 .and. kbotphi==1 ) then ! Dirichlet on both sides
         do n_m=nMstart,nMstop
            m=idx2m(n_m)
            if ( m == 0 ) then
               call rscheme%robin_bc(0.0_cp, one, cmplx(phi_top,0.0_cp,cp), &
                    &                one, 0.0_cp, cmplx(phi_bot,0.0_cp,cp), &
                    &                phi_Mloc(n_m,:))
            else
               call rscheme%robin_bc(0.0_cp, one, zero, one, 0.0_cp, zero, &
                    &                phi_Mloc(n_m,:))
            end if
         end do
      else if ( ktopphi==1 .and. kbotphi /= 1 ) then ! Dirichlet: top and Neumann: bot
         do n_m=nMstart,nMstop
            m=idx2m(n_m)
            if ( m == 0 ) then
               call rscheme%robin_bc(0.0_cp, one, cmplx(phi_top,0.0_cp,cp), &
                    &                one, 0.0_cp, zero, phi_Mloc(n_m,:))
            else
               call rscheme%robin_bc(0.0_cp, one, zero, one, 0.0_cp, &
                    &                zero, phi_Mloc(n_m,:))
            end if
         end do
      else if ( kbotphi==1 .and. ktopphi /= 1 ) then ! Dirichlet: bot and Neumann: top
         do n_m=nMstart,nMstop
            m=idx2m(n_m)
            if ( m == 0 ) then
               call rscheme%robin_bc(one, 0.0_cp, zero, 0.0_cp, one, &
                    &                cmplx(phi_bot,0.0_cp,cp), phi_Mloc(n_m,:))
            else
               call rscheme%robin_bc(one, 0.0_cp, zero, 0.0_cp, one, &
                    &                zero, phi_Mloc(n_m,:))
            end if
         end do
      else if ( kbotphi /= 1 .and. kbotphi /= 1 ) then ! Neumann on both sides
         do n_m=nMstart,nMstop
            call rscheme%robin_bc(one, 0.0_cp, zero, one, 0.0_cp, zero, &
                 &                phi_Mloc(n_m,:))
         end do
      end if

      call get_phi_rhs_imp_coll(phi_Mloc, dphidt, 1, tscheme%l_imp_calc_rhs(1),  &
           &                    l_in_cheb_space=.false.)

   end subroutine assemble_phi_coll
!------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_phiMat(tscheme, m, tMat, tPivot, tMat_fac)
#else
   subroutine get_phiMat(tscheme, m, tMat, tPivot)
#endif

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step
      integer,             intent(in) :: m

      !-- Output variables
      real(cp), intent(out) :: tMat(n_r_max,n_r_max)
      integer,  intent(out) :: tPivot(n_r_max)
#ifdef WITH_PRECOND_S
      real(cp),intent(out) :: tMat_fac(n_r_max)
#endif

      !-- Local variables
      integer :: nR_out, nR, info, n_m
      real(cp) :: dm2

      dm2 = real(m,cp)*real(m,cp)
      n_m = m2idx(m)

      !----- Boundary coditions:
      do nR_out=1,rscheme%n_max
         if ( ktopphi == 1 ) then
            tMat(1,nR_out)=rscheme%rnorm*rscheme%rMat(1,nR_out)
         else
            tMat(1,nR_out)=rscheme%rnorm*rscheme%drMat(1,nR_out)
         end if
         if ( l_full_disk ) then
            if ( m == 0 ) then
               tMat(n_r_max,nR_out)=rscheme%rnorm*rscheme%drMat(n_r_max,nR_out)
            else
               tMat(n_r_max,nR_out)=rscheme%rnorm*rscheme%rMat(n_r_max,nR_out)
            end if
         else
            if ( kbotphi == 1 ) then
               tMat(n_r_max,nR_out)=rscheme%rnorm*rscheme%rMat(n_r_max,nR_out)
            else
               tMat(n_r_max,nR_out)=rscheme%rnorm*rscheme%drMat(n_r_max,nR_out)
            end if
         end if
      end do

      if ( rscheme%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme%n_max+1,n_r_max
            tMat(1,nR_out)      =0.0_cp
            tMat(n_r_max,nR_out)=0.0_cp
         end do
      end if

      !----- Other points:
      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            tMat(nR,nR_out)= rscheme%rnorm * (                          &
            &         5.0_cp/6.0_cp*stef*pr*  rscheme%rMat(nR,nR_out) - &
            &                      tscheme%wimp_lin(1)*phaseDiffFac*(   &
            &                               rscheme%d2rMat(nR,nR_out) + &
            &          or1(nR)*              rscheme%drMat(nR,nR_out) - &
            &      dm2*or2(nR)*               rscheme%rMat(nR,nR_out) ) )
         end do
      end do

      !----- Factor for highest and lowest cheb:
      tMat(:,1)      =rscheme%boundary_fac*tMat(:,1)
      tMat(:,n_r_max)=rscheme%boundary_fac*tMat(:,n_r_max)

#ifdef WITH_PRECOND_S
      ! compute the linesum of each line
      do nR=1,n_r_max
         tMat_fac(nR)=one/maxval(abs(tMat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nR=1,n_r_max
         tMat(nR,:) = tMat(nR,:)*tMat_fac(nR)
      end do
#endif

      !----- LU decomposition:
      call prepare_full_mat(tMat,n_r_max,n_r_max,tPivot,info)
      if ( info /= 0 ) call abortRun('Singular matrix phiMat!')

   end subroutine get_phiMat
!------------------------------------------------------------------------------
end module update_phi_coll
