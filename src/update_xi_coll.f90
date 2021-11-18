module update_xi_coll

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: one, zero, four, ci
   use namelists, only: kbotxi, ktopxi, XidiffFac, ChemFac, l_heat
   use radial_functions, only: rscheme, or1, or2, dxicond, rgrav
   use horizontal, only: hdif_xi, botxi_Mloc, topxi_Mloc
   use blocking, only: nMstart, nMstop
   use truncation, only: n_r_max, idx2m, m2idx
   use radial_der, only: get_ddr, get_dr
   use fields, only: work_Mloc
   use algebra, only: prepare_full_mat, solve_full_mat
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray

   implicit none
   
   private

   logical,  allocatable :: lXimat(:)
   real(cp), allocatable :: xiMat(:, :, :)
   integer,  allocatable :: xiPivot(:, :)
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: xiMat_fac(:,:)
#endif
   complex(cp), allocatable :: rhs(:)

   public :: update_xi_co, initialize_xi_coll, finalize_xi_coll, &
   &         get_xi_rhs_imp_coll, finish_exp_xi_coll, assemble_xi_coll

contains

   subroutine initialize_xi_coll

      allocate( lXimat(nMstart:nMstop) )
      lXimat(:)=.false.
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_LOGICAL

      allocate( xiMat(n_r_max, n_r_max, nMstart:nMstop) )
      allocate( xiPivot(n_r_max, nMstart:nMstop) )
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*n_r_max*n_r_max* &
      &                 SIZEOF_DEF_REAL+n_r_max*(nMstop-nMstart+1)*SIZEOF_INTEGER
#ifdef WITH_PRECOND_S
      allocate( xiMat_fac(n_r_max, nMstart:nMstop) )
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*n_r_max*  &
      &                 SIZEOF_DEF_REAL
#endif
      allocate( rhs(n_r_max) )
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_COMPLEX

   end subroutine initialize_xi_coll
!------------------------------------------------------------------------------
   subroutine finalize_xi_coll

      deallocate( rhs )
#ifdef WITH_PRECOND_S
      deallocate( xiMat_fac )
#endif
      deallocate( lXimat, xiMat, xiPivot )

   end subroutine finalize_xi_coll
!------------------------------------------------------------------------------
   subroutine update_xi_co(xi_Mloc, dxi_Mloc, dxidt, tscheme, lMat, l_log_next)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lMat
      logical,             intent(in) :: l_log_next

      !-- Output variables
      complex(cp),       intent(out) :: xi_Mloc(nMstart:nMstop, n_r_max)
      complex(cp),       intent(out) :: dxi_Mloc(nMstart:nMstop, n_r_max)
      type(type_tarray), intent(inout) :: dxidt

      !-- Local variables
      integer :: n_r, n_m, n_r_out, m

      if ( lMat ) lXiMat(:)=.false.

      !-- Now assemble the right hand side and store it in work_Mloc
      call tscheme%set_imex_rhs(work_Mloc, dxidt)

      do n_m=nMstart, nMstop

         m = idx2m(n_m)
         
         if ( .not. lXimat(n_m) ) then
#ifdef WITH_PRECOND_S
            call get_xiMat(tscheme, m, xiMat(:,:,n_m), xiPivot(:,n_m), &
                 &         xiMat_fac(:,n_m))
#else
            call get_xiMat(tscheme, m, xiMat(:,:,n_m), xiPivot(:,n_m))
#endif
            lXimat(n_m)=.true.
         end if

         !-- Inhomogeneous B.Cs (if not zero)
         rhs(1)      =topxi_Mloc(n_m)
         rhs(n_r_max)=botxi_Mloc(n_m)
         do n_r=2,n_r_max-1
            rhs(n_r)=work_Mloc(n_m,n_r)
         end do

#ifdef WITH_PRECOND_S
         do n_r=1,n_r_max
            rhs(n_r) = xiMat_fac(n_r,n_m)*rhs(n_r)
         end do
#endif

         call solve_full_mat(xiMat(:,:,n_m), n_r_max, n_r_max, xiPivot(:, n_m), &
              &              rhs(:))

         do n_r_out=1,rscheme%n_max
            xi_Mloc(n_m, n_r_out)=rhs(n_r_out)
         end do

      end do

      !-- set cheb modes > rscheme%n_max to zero (dealiazing)
      do n_r_out=rscheme%n_max+1,n_r_max
         do n_m=nMstart,nMstop
            xi_Mloc(n_m,n_r_out)=zero
         end do
      end do

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dxidt)

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_xi_rhs_imp_coll(xi_Mloc, dxi_Mloc, dxidt, 1, &
              &                   tscheme%l_imp_calc_rhs(1), l_in_cheb_space=.true.)
      else
         call get_xi_rhs_imp_coll(xi_Mloc, dxi_Mloc, dxidt, tscheme%istage+1, &
              &                   tscheme%l_imp_calc_rhs(tscheme%istage+1),   &
              &                   l_in_cheb_space=.true.)
      end if

      !-- In case log is needed on the next iteration, recalculate dT/dr
      if ( l_log_next ) then
         call get_dr(xi_Mloc, dxi_Mloc, nMstart, nMstop, n_r_max, rscheme)
      end if

   end subroutine update_xi_co
!------------------------------------------------------------------------------
   subroutine finish_exp_xi_coll(us_Mloc, dVsXi_Mloc, dxi_exp_last)

      !-- Input variables
      complex(cp), intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: dVsXi_Mloc(nMstart:nMstop,n_r_max)

      !-- Output variables
      complex(cp), intent(inout) :: dxi_exp_last(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_m

      !-- Finish calculation of advection
      call get_dr( dVsXi_Mloc, work_Mloc, nMstart, nMstop, n_r_max, &
           &       rscheme, nocopy=.true. )

      !-- Finish calculation of the explicit part for current time step
      do n_r=1,n_r_max
         do n_m=nMstart, nMstop
            dxi_exp_last(n_m,n_r)=dxi_exp_last(n_m,n_r)         & 
            &                     -or1(n_r)*work_Mloc(n_m,n_r)  &
            &                     -us_Mloc(n_m,n_r)*dxicond(n_r)
         end do
      end do

   end subroutine finish_exp_xi_coll
!------------------------------------------------------------------------------
   subroutine get_xi_rhs_imp_coll(xi_Mloc, dxi_Mloc, dxidt, istage, l_calc_lin, &
              &                   l_in_cheb_space)

      !-- Input variables
      complex(cp),       intent(inout) :: xi_Mloc(nMstart:nMstop,n_r_max)
      integer,           intent(in) :: istage
      logical,           intent(in) :: l_calc_lin
      logical, optional, intent(in) :: l_in_cheb_space

      !-- Output variable
      complex(cp),       intent(out) :: dxi_Mloc(nMstart:nMstop,n_r_max)
      type(type_tarray), intent(inout) :: dxidt

      !-- Local variables
      logical :: l_in_cheb
      integer :: n_r, n_m, m
      real(cp) :: dm2

      if ( present(l_in_cheb_space) ) then
         l_in_cheb = l_in_cheb_space
      else
         l_in_cheb = .false.
      end if

      if ( l_calc_lin ) then
         call get_ddr(xi_Mloc, dxi_Mloc, work_Mloc, nMstart, nMstop, &
              &       n_r_max, rscheme, l_dct=.not. l_in_cheb)
         if ( l_in_cheb ) call rscheme%costf1(xi_Mloc, nMstart, nMstop, n_r_max)
      else
         if ( l_in_cheb ) call rscheme%costf1(xi_Mloc, nMstart, nMstop, n_r_max)
      end if

      if ( istage == 1 ) then
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               dxidt%old(n_m,n_r,istage)=xi_Mloc(n_m,n_r)
            end do
         end do
      end if

      if ( l_calc_lin ) then
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               dm2 = real(m,cp)*real(m,cp)
               dxidt%impl(n_m,n_r,istage)=XidiffFac*hdif_Xi(n_m)* (      &
               &                                     work_Mloc(n_m,n_r)  &
               &                         +or1(n_r)*   dxi_Mloc(n_m,n_r)  &
               &                     -dm2*or2(n_r)*    xi_Mloc(n_m,n_r) )
            end do
         end do
      end if

   end subroutine get_xi_rhs_imp_coll
!------------------------------------------------------------------------------
   subroutine assemble_xi_coll(xi_Mloc, dxi_Mloc, dxidt, tscheme, l_log_next)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: l_log_next
      type(type_tarray),   intent(inout) :: dxidt

      !-- Output variable
      complex(cp), intent(inout) :: xi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(out) :: dxi_Mloc(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_m, m

      call tscheme%assemble_imex(work_Mloc, dxidt)

      do n_r=2,n_r_max-1
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m == 0 ) then
               xi_Mloc(n_m,n_r) = cmplx(real(work_Mloc(n_m,n_r)), 0.0_cp, cp)
            else
               xi_Mloc(n_m,n_r) = work_Mloc(n_m,n_r)
            end if
         end do
      end do

      if ( ktopxi==1 .and. kbotxi==1 ) then ! Dirichlet on both sides
         do n_m=nMstart,nMstop
            call rscheme%robin_bc(0.0_cp, one, topxi_Mloc(n_m), 0.0_cp, one, &
                 &                botxi_Mloc(n_m), xi_Mloc(n_m,:))
         end do
      else if ( ktopxi==1 .and. kbotxi /= 1 ) then ! Dirichlet: top and Neumann: bot
         do n_m=nMstart,nMstop
            call rscheme%robin_bc(0.0_cp, one, topxi_Mloc(n_m), one, 0.0_cp, &
                 &                botxi_Mloc(n_m), xi_Mloc(n_m,:))
         end do
      else if ( kbotxi==1 .and. ktopxi /= 1 ) then ! Dirichlet: bot and Neumann: top
         do n_m=nMstart,nMstop
            call rscheme%robin_bc(one, 0.0_cp, topxi_Mloc(n_m), 0.0_cp, one, &
                 &                botxi_Mloc(n_m), xi_Mloc(n_m,:))
         end do
      else if ( kbotxi /= 1 .and. kbotxi /= 1 ) then ! Neumann on both sides
         do n_m=nMstart,nMstop
            call rscheme%robin_bc(one, 0.0_cp, topxi_Mloc(n_m), one, 0.0_cp, &
                 &                botxi_Mloc(n_m), xi_Mloc(n_m,:))
         end do
      end if

      call get_xi_rhs_imp_coll(xi_Mloc, dxi_Mloc, dxidt, 1,   &
           &                   tscheme%l_imp_calc_rhs(1),     &
           &                   l_in_cheb_space=.false.)

   end subroutine assemble_xi_coll
!------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_xiMat(tscheme, m, tMat, tPivot, tMat_fac)
#else
   subroutine get_xiMat(tscheme, m, tMat, tPivot)
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
         if ( ktopxi == 1 ) then
            tMat(1,nR_out)=rscheme%rnorm*rscheme%rMat(1,nR_out)
         else
            tMat(1,nR_out)=rscheme%rnorm*rscheme%drMat(1,nR_out)
         end if
         if ( kbotxi == 1 ) then
            tMat(n_r_max,nR_out)=rscheme%rnorm*rscheme%rMat(n_r_max,nR_out)
         else
            tMat(n_r_max,nR_out)=rscheme%rnorm*rscheme%drMat(n_r_max,nR_out)
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
            &                                 rscheme%rMat(nR,nR_out) - &
            &            tscheme%wimp_lin(1)*XidiffFac*hdif_Xi(n_m)*(   &
            &                               rscheme%d2rMat(nR,nR_out) + &
            &          or1(nR)*              rscheme%drMat(nR,nR_out) - &
            &      dm2*or2(nR)*               rscheme%rMat(nR,nR_out) ) )
         end do
      end do

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         tMat(nR,1)      =rscheme%boundary_fac*tMat(nR,1)
         tMat(nR,n_r_max)=rscheme%boundary_fac*tMat(nR,n_r_max)
      end do

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
      if ( info /= 0 ) then
         call abortRun('Singular matrix tMat!')
      end if

   end subroutine get_xiMat
!------------------------------------------------------------------------------
end module update_xi_coll
