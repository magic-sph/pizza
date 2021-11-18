module update_temp_coll

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: one, zero
   use namelists, only: kbott, ktopt, tadvz_fac, TdiffFac, BuoFac
   use radial_functions, only: rscheme, or1, or2, dtcond, tcond, beta, &
       &                       rgrav
   use horizontal, only: hdif_T, bott_Mloc, topt_Mloc
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

   logical,  allocatable :: lTmat(:)
   real(cp), allocatable :: tMat(:, :, :)
   integer,  allocatable :: tPivot(:, :)
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: tMat_fac(:,:)
#endif
   complex(cp), allocatable :: rhs(:)

   public :: update_temp_co, initialize_temp_coll, finalize_temp_coll, &
   &         get_temp_rhs_imp_coll, finish_exp_temp_coll, assemble_temp_coll

contains

   subroutine initialize_temp_coll

      allocate( lTmat(nMstart:nMstop) )
      lTmat(:)=.false.
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_LOGICAL

      allocate( tMat(n_r_max, n_r_max, nMstart:nMstop) )
      allocate( tPivot(n_r_max, nMstart:nMstop) )
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*n_r_max*n_r_max* &
      &                 SIZEOF_DEF_REAL+n_r_max*(nMstop-nMstart+1)*SIZEOF_INTEGER
#ifdef WITH_PRECOND_S
      allocate( tMat_fac(n_r_max, nMstart:nMstop) )
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*n_r_max*  &
      &                 SIZEOF_DEF_REAL
#endif
      allocate( rhs(n_r_max) )
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_COMPLEX

   end subroutine initialize_temp_coll
!------------------------------------------------------------------------------
   subroutine finalize_temp_coll

      deallocate( rhs )
#ifdef WITH_PRECOND_S
      deallocate( tMat_fac )
#endif
      deallocate( lTmat, tMat, tPivot )

   end subroutine finalize_temp_coll
!------------------------------------------------------------------------------
   subroutine update_temp_co(temp_Mloc, dtemp_Mloc, dTdt, tscheme, lMat, l_log_next)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lMat
      logical,             intent(in) :: l_log_next

      !-- Output variables
      complex(cp),       intent(out) :: temp_Mloc(nMstart:nMstop, n_r_max)
      complex(cp),       intent(out) :: dtemp_Mloc(nMstart:nMstop, n_r_max)
      type(type_tarray), intent(inout) :: dTdt

      !-- Local variables
      integer :: n_r, n_m, n_r_out, m

      if ( lMat ) lTMat(:)=.false.

      !-- Now assemble the right hand side and store it in work_Mloc
      call tscheme%set_imex_rhs(work_Mloc, dTdt)

      do n_m=nMstart, nMstop

         m = idx2m(n_m)
         
         if ( .not. lTmat(n_m) ) then
#ifdef WITH_PRECOND_S
            call get_tempMat(tscheme, m, tMat(:,:,n_m), tPivot(:,n_m), &
                 &           tMat_fac(:,n_m))
#else
            call get_tempMat(tscheme, m, tMat(:,:,n_m), tPivot(:,n_m))
#endif
            lTmat(n_m)=.true.
         end if

         !-- Inhomogeneous B.Cs (if not zero)
         rhs(1)      =topt_Mloc(n_m)
         rhs(n_r_max)=bott_Mloc(n_m)
         do n_r=2,n_r_max-1
            rhs(n_r)=work_Mloc(n_m,n_r)
         end do

#ifdef WITH_PRECOND_S
         do n_r=1,n_r_max
            rhs(n_r) = tMat_fac(n_r,n_m)*rhs(n_r)
         end do
#endif

         call solve_full_mat(tMat(:,:,n_m), n_r_max, n_r_max, tPivot(:, n_m), &
              &              rhs(:))

         do n_r_out=1,rscheme%n_max
            temp_Mloc(n_m, n_r_out)=rhs(n_r_out)
         end do

      end do

      !-- set cheb modes > rscheme%n_max to zero (dealiazing)
      do n_r_out=rscheme%n_max+1,n_r_max
         do n_m=nMstart,nMstop
            temp_Mloc(n_m,n_r_out)=zero
         end do
      end do

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dTdt)

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_temp_rhs_imp_coll(temp_Mloc, dtemp_Mloc, dTdt, 1,   &
              &                     tscheme%l_imp_calc_rhs(1),        &
              &                     l_in_cheb_space=.true.)
      else
         call get_temp_rhs_imp_coll(temp_Mloc, dtemp_Mloc, dTdt, tscheme%istage+1, &
              &                     tscheme%l_imp_calc_rhs(tscheme%istage+1),      &
              &                     l_in_cheb_space=.true.)
      end if


      !-- In case log is needed on the next iteration, recalculate dT/dr
      if ( l_log_next ) then
         call get_dr(temp_Mloc, dtemp_Mloc, nMstart, nMstop, n_r_max, rscheme)
      end if

   end subroutine update_temp_co
!------------------------------------------------------------------------------
   subroutine finish_exp_temp_coll(us_Mloc, dVsT_Mloc, dtemp_exp_last)

      !-- Input variables
      complex(cp), intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: dVsT_Mloc(nMstart:nMstop,n_r_max)

      !-- Output variables
      complex(cp), intent(inout) :: dtemp_exp_last(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_m

      !-- Finish calculation of advection
      call get_dr( dVsT_Mloc, work_Mloc, nMstart, nMstop, n_r_max, &
           &       rscheme, nocopy=.true. )

      !-- Finish calculation of the explicit part for current time step
      do n_r=1,n_r_max
         do n_m=nMstart, nMstop
            dtemp_exp_last(n_m,n_r)=dtemp_exp_last(n_m,n_r)         & 
            &                       -or1(n_r)*work_Mloc(n_m,n_r)    &
            &                       -us_Mloc(n_m,n_r)*(dtcond(n_r)- &
            &                       tadvz_fac*beta(n_r)*tcond(n_r))
         end do
      end do

   end subroutine finish_exp_temp_coll
!------------------------------------------------------------------------------
   subroutine get_temp_rhs_imp_coll(temp_Mloc, dtemp_Mloc, dTdt, istage, &
              &                     l_calc_lin, l_in_cheb_space)

      !-- Input variables
      complex(cp),       intent(inout) :: temp_Mloc(nMstart:nMstop,n_r_max)
      integer,           intent(in) :: istage
      logical,           intent(in) :: l_calc_lin
      logical, optional, intent(in) :: l_in_cheb_space

      !-- Output variable
      type(type_tarray), intent(inout) :: dTdt
      complex(cp),       intent(out) :: dtemp_Mloc(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_m, m
      logical :: l_in_cheb
      real(cp) :: dm2

      if ( present(l_in_cheb_space) ) then
         l_in_cheb = l_in_cheb_space
      else
         l_in_cheb = .false.
      end if

      if ( l_calc_lin ) then
         call get_ddr(temp_Mloc, dtemp_Mloc, work_Mloc, nMstart, nMstop, &
              &       n_r_max, rscheme, l_dct=.not. l_in_cheb)
         if ( l_in_cheb ) call rscheme%costf1(temp_Mloc, nMstart, nMstop, n_r_max)
      else
         if ( l_in_cheb ) call rscheme%costf1(temp_Mloc, nMstart, nMstop, n_r_max)
      end if

      if ( istage == 1 ) then
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               dTdt%old(n_m,n_r,istage)=temp_Mloc(n_m,n_r)
            end do
         end do
      end if

      if ( l_calc_lin ) then
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               dm2 = real(m,cp)*real(m,cp)
               dTdt%impl(n_m,n_r,istage)=TdiffFac*hdif_T(n_m)* (           &
               &                                       work_Mloc(n_m,n_r)  &
               &                         +or1(n_r)*   dtemp_Mloc(n_m,n_r)  &
               &                     -dm2*or2(n_r)*    temp_Mloc(n_m,n_r) )
            end do
         end do
      end if

   end subroutine get_temp_rhs_imp_coll
!------------------------------------------------------------------------------
   subroutine assemble_temp_coll(temp_Mloc, dtemp_Mloc, dTdt, tscheme, l_log_next)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: l_log_next
      type(type_tarray),   intent(inout) :: dTdt

      !-- Output variable
      complex(cp), intent(inout) :: temp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(out) :: dtemp_Mloc(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_m, m

      call tscheme%assemble_imex(work_Mloc, dTdt)

      do n_r=2,n_r_max-1
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m == 0 ) then
               temp_Mloc(n_m,n_r) = cmplx(real(work_Mloc(n_m,n_r)), 0.0_cp, cp)
            else
               temp_Mloc(n_m,n_r) = work_Mloc(n_m,n_r)
            end if
         end do
      end do

      if ( ktopt==1 .and. kbott==1 ) then ! Dirichlet on both sides
         do n_m=nMstart,nMstop
            call rscheme%robin_bc(0.0_cp, one, topt_Mloc(n_m), 0.0_cp, one, &
                 &                bott_Mloc(n_m), temp_Mloc(n_m,:))
         end do
      else if ( ktopt==1 .and. kbott /= 1 ) then ! Dirichlet: top and Neumann: bot
         do n_m=nMstart,nMstop
            call rscheme%robin_bc(0.0_cp, one, topt_Mloc(n_m), one, 0.0_cp, &
                 &                bott_Mloc(n_m), temp_Mloc(n_m,:))
         end do
      else if ( kbott==1 .and. ktopt /= 1 ) then ! Dirichlet: bot and Neumann: top
         do n_m=nMstart,nMstop
            call rscheme%robin_bc(one, 0.0_cp, topt_Mloc(n_m), 0.0_cp, one, &
                 &                bott_Mloc(n_m), temp_Mloc(n_m,:))
         end do
      else if ( kbott /= 1 .and. kbott /= 1 ) then ! Neumann on both sides
         do n_m=nMstart,nMstop
            call rscheme%robin_bc(one, 0.0_cp, topt_Mloc(n_m), one, 0.0_cp, &
                 &                bott_Mloc(n_m), temp_Mloc(n_m,:))
         end do
      end if

      !-- Finally call the construction of the implicit terms for the first stage
      !-- of next iteration
      call get_temp_rhs_imp_coll(temp_Mloc, dtemp_Mloc, dTdt, 1,   &
           &                     tscheme%l_imp_calc_rhs(1),        &
           &                     l_in_cheb_space=.false.)

      !-- In case log is needed on the next iteration, recalculate dT/dr
      if ( l_log_next ) then
         call get_dr(temp_Mloc, dtemp_Mloc, nMstart, nMstop, n_r_max, rscheme)
      end if

   end subroutine assemble_temp_coll
!------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_tempMat(tscheme, m, tMat, tPivot, tMat_fac)
#else
   subroutine get_tempMat(tscheme, m, tMat, tPivot)
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
         if ( ktopt == 1 ) then
            tMat(1,nR_out)=rscheme%rnorm*rscheme%rMat(1,nR_out)
         else
            tMat(1,nR_out)=rscheme%rnorm*rscheme%drMat(1,nR_out)
         end if
         if ( kbott == 1 ) then
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
            &              tscheme%wimp_lin(1)*TdiffFac*hdif_T(n_m)*(   &
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

   end subroutine get_tempMat
!------------------------------------------------------------------------------
end module update_temp_coll
