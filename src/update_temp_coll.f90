module update_temp_coll

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: one, zero, four, ci
   use namelists, only: kbott, ktopt, tadvz_fac, TdiffFac, BuoFac
   use radial_functions, only: rscheme, or1, or2, dtcond, tcond, beta, &
       &                       rgrav
   use blocking, only: nMstart, nMstop
   use truncation, only: n_r_max, idx2m
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
   &         get_temp_rhs_imp_coll

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
   subroutine update_temp_co(us_Mloc, temp_Mloc, dtemp_Mloc, dVsT_Mloc,      &
              &              buo_imp_Mloc, dTdt, tscheme, lMat, l_log_next)

      !-- Input variables
      type(type_tscheme), intent(in) :: tscheme
      logical,            intent(in) :: lMat
      logical,            intent(in) :: l_log_next
      complex(cp),        intent(in) :: us_Mloc(nMstart:nMstop, n_r_max)

      !-- Output variables
      complex(cp), intent(out) :: temp_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(out) :: dtemp_Mloc(nMstart:nMstop, n_r_max)
      type(type_tarray), intent(inout) :: dTdt
      complex(cp), intent(inout) :: buo_imp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: dVsT_Mloc(nMstart:nMstop, n_r_max)

      !-- Local variables
      integer :: n_r, n_m, n_r_out, m, n_o

      if ( lMat ) lTMat(:)=.false.

      !-- Finish calculation of advection
      call get_dr( dVsT_Mloc, work_Mloc, nMstart, nMstop, n_r_max, &
           &       rscheme, nocopy=.true. )

      !-- Finish calculation of the explicit part for current time step
      do n_r=1,n_r_max
         do n_m=nMstart, nMstop
            dTdt%expl(n_m,n_r,1)=dTdt%expl(n_m,n_r,1)             &
            &                     -or1(n_r)*work_Mloc(n_m,n_r)    &
            &                     -us_Mloc(n_m,n_r)*(dtcond(n_r)- &
            &                     tadvz_fac*beta(n_r)*tcond(n_r))
         end do
      end do

      !-- Calculation of the implicit part
      call get_temp_rhs_imp_coll(temp_Mloc, dtemp_Mloc, dTdt%old(:,:,1), &
           &                     dTdt%impl(:,:,1), tscheme%l_calc_lin_rhs)

      !-- Now assemble the right hand side and store it in work_Mloc
      call tscheme%set_imex_rhs(work_Mloc, dTdt, nMstart, nMstop, n_r_max)

      do n_m=nMstart, nMstop

         m = idx2m(n_m)
         
         if ( .not. lTmat(n_m) ) then
#ifdef WITH_PRECOND_S
            call get_tempMat(tscheme, m, tMat(:,:,n_m), tPivot(:,n_m), tMat_fac(:,n_m))
#else
            call get_tempMat(tscheme, m, tMat(:,:,n_m), tPivot(:,n_m))
#endif
            lTmat(n_m)=.true.
         end if

         rhs(1)      =zero
         rhs(n_r_max)=zero
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

      !-- Bring temperature back to physical space
      call rscheme%costf1(temp_Mloc, nMstart, nMstop, n_r_max)

      !-- Assemble buoyancy
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m /= 0 ) then
               buo_imp_Mloc(n_m,n_r)=-tscheme%wimp_lin(1)*rgrav(n_r)*   &
               &                     or1(n_r)*BuoFac*ci*real(m,cp)*     &
               &                     temp_Mloc(n_m,n_r)
               do n_o = 1,tscheme%norder_imp_lin-1
                  buo_imp_Mloc(n_m,n_r)=buo_imp_Mloc(n_m,n_r)-             &
                  &                     tscheme%wimp_lin(n_o+1)*rgrav(n_r)*&
                  &                     or1(n_r)*BuoFac*ci*real(m,cp)*     &
                  &                     dTdt%old(n_m,n_r,n_o)
               end do
            end if
         end do
      end do

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dTdt, nMstart, nMstop, n_r_max)

      !-- In case log is needed on the next iteration, recalculate dT/dr
      if ( l_log_next ) then
         call get_dr(temp_Mloc, dtemp_Mloc, nMstart, nMstop, n_r_max, rscheme)
      end if

   end subroutine update_temp_co
!------------------------------------------------------------------------------
   subroutine get_temp_rhs_imp_coll(temp_Mloc, dtemp_Mloc, temp_last, &
              &                     dtemp_imp_Mloc_last, l_calc_lin_rhs)

      !-- Input variables
      complex(cp), intent(in) :: temp_Mloc(nMstart:nMstop,n_r_max)
      logical,     intent(in) :: l_calc_lin_rhs

      !-- Output variable
      complex(cp), intent(out) :: dtemp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(out) :: temp_last(nMstart:nMstop,n_r_max)
      complex(cp), intent(out) :: dtemp_imp_Mloc_last(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_m, m
      real(cp) :: dm2

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            temp_last(n_m,n_r)=temp_Mloc(n_m,n_r)
         end do
      end do

      if ( l_calc_lin_rhs ) then

         call get_ddr(temp_Mloc, dtemp_Mloc, work_Mloc, nMstart, nMstop, &
              &       n_r_max, rscheme)
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               dm2 = real(m,cp)*real(m,cp)
               dtemp_imp_Mloc_last(n_m,n_r)=TdiffFac*( work_Mloc(n_m,n_r)  &
               &                         +or1(n_r)*   dtemp_Mloc(n_m,n_r)  &
               &                     -dm2*or2(n_r)*    temp_Mloc(n_m,n_r) )
            end do
         end do

      end if

   end subroutine get_temp_rhs_imp_coll
!------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_tempMat(tscheme, m, tMat, tPivot, tMat_fac)
#else
   subroutine get_tempMat(tscheme, m, tMat, tPivot)
#endif

      !-- Input variables
      type(type_tscheme), intent(in) :: tscheme        ! time step
      integer,            intent(in) :: m

      !-- Output variables
      real(cp), intent(out) :: tMat(n_r_max,n_r_max)
      integer,  intent(out) :: tPivot(n_r_max)
#ifdef WITH_PRECOND_S
      real(cp),intent(out) :: tMat_fac(n_r_max)
#endif

      !-- Local variables
      integer :: nR_out, nR, info
      real(cp) :: dm2

      dm2 = real(m,cp)*real(m,cp)

      !----- Boundary coditions:
      do nR_out=1,rscheme%n_max
         if ( ktopt == 1 ) then
            tMat(1,nR_out)=rscheme%rnorm*rscheme%rMat(1,nR_out)
         else
            tMat(1,nR_out)=rscheme%rnorm*rscheme%drMat(1,nR_out)
         end if
         if ( kbott == 1 ) then
            tMat(n_r_max,nR_out)=rscheme%rnorm* &
            &                     rscheme%rMat(n_r_max,nR_out)
         else
            tMat(n_r_max,nR_out)=rscheme%rnorm* &
            &                     rscheme%drMat(n_r_max,nR_out)
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
            &tscheme%wimp_lin(1)*TdiffFac*( rscheme%d2rMat(nR,nR_out) + &
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
