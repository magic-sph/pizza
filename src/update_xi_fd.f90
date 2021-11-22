module update_xi_fd_mod

   use precision_mod
   use parallel_mod, only: get_openmp_blocks
   use mem_alloc, only: bytes_allocated
   use constants, only: zero, one, two
   use truncation, only: n_m_max, n_r_max, idx2m
   use blocking, only: nRstart, nRstop
   use parallel_solvers, only: type_tri_par
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use namelists, only: XiDiffFac, kbotxi, ktopxi
   use radial_functions, only: or1, or2, rscheme, r, beta, dxicond
   use radial_der, only: get_ddr_ghost, get_dr_Rloc, exch_ghosts, bulk_to_ghost

   implicit none

   private

   type(type_tri_par), public :: xiMat_FD
   complex(cp), public, allocatable :: xi_ghost(:,:)
   logical, public, allocatable :: lXimat_FD(:)

   public :: initialize_xi_fd, finalize_xi_fd, finish_exp_xi_Rdist, &
   &         prepare_xi_fd, fill_ghosts_xi, update_xi_fd,           &
   &         get_xi_rhs_imp_ghost, assemble_xi_Rdist

contains

   subroutine initialize_xi_fd

      call xiMat_FD%initialize(nRstart,nRstop,1,n_m_max)
      allocate( xi_ghost(n_m_max,nRstart-1:nRstop+1) )
      bytes_allocated=bytes_allocated + n_m_max*(nRstop-nRstart+3)*SIZEOF_DEF_COMPLEX
      xi_ghost(:,:)=zero
      allocate( lXimat_FD(1:n_m_max) ) 

   end subroutine initialize_xi_fd
!---------------------------------------------------------------------------------
   subroutine finalize_xi_fd

      deallocate( lXimat_FD, xi_ghost)
      call xiMat_FD%finalize()

   end subroutine finalize_xi_fd
!---------------------------------------------------------------------------------
   subroutine prepare_xi_fd(tscheme, dxidt)
      !
      ! This subroutine is used to assemble the r.h.s. of the xierature equation
      ! when parallel F.D solvers are used. Boundary values are set here.
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      type(type_tarray), intent(inout) :: dxidt

      !-- Local variables:
      integer :: n_m_start, n_m_stop, nR, n_m

      if ( .not. lXimat_FD(1) ) then
         call get_xiMat_Rdist(tscheme, xiMat_FD)
         lXimat_FD(:)=.true.
      end if

      !$omp parallel default(shared) private(n_m_start,n_m_stop, nR, n_m)
      n_m_start=1; n_m_stop=n_m_max
      call get_openmp_blocks(n_m_start,n_m_stop)
      !$omp barrier

      !-- Now assemble the right hand side
      call tscheme%set_imex_rhs_ghost(xi_ghost, dxidt, n_m_start, n_m_stop, 1)

      !-- Set boundary conditions
      if ( nRstart == 1 ) then
         nR=1
         do n_m=n_m_start,n_m_stop
            if ( ktopxi == 1 ) then ! Fixed xierature
               xi_ghost(n_m,nR)=zero
            end if
            xi_ghost(n_m,nR-1)=zero ! Set ghost zone to zero
         end do
      end if

      if ( nRstop == n_r_max ) then
         nR=n_r_max
         do n_m=n_m_start,n_m_stop
            if ( kbotxi == 1 ) then ! Fixed xierature
                  xi_ghost(n_m,nR)=zero
            end if
            xi_ghost(n_m,nR+1)=zero ! Set ghost zone to zero
         end do
      end if
      !$omp end parallel

   end subroutine prepare_xi_fd
!---------------------------------------------------------------------------------
   subroutine fill_ghosts_xi(xig)
      !
      ! This subroutine is used to fill the ghosts zones that are located at
      ! nR=n_r_cmb-1 and nR=n_r_icb+1. This is used to properly set the Neuman
      ! boundary conditions. In case Dirichlet BCs are used, a simple first order
      ! extrapolation is employed. This is anyway only used for outputs (like Nusselt
      ! numbers).
      !
      complex(cp), intent(inout) :: xig(n_m_max,nRstart-1:nRstop+1)

      !-- Local variables
      integer :: n_m, n_m_start, n_m_stop

      !$omp parallel default(shared) private(n_m_start, n_m_stop, n_m)
      n_m_start=1; n_m_stop=n_m_max
      call get_openmp_blocks(n_m_start,n_m_stop)
      !$omp barrier

      !-- Handle upper boundary
      !dr = r(2)-r(1)
      if ( nRstart == 1 ) then
         do n_m=n_m_start,n_m_stop
            if ( ktopxi == 1 ) then
               xig(n_m,nRstart-1)=two*xig(n_m,nRstart)-xig(n_m,nRstart+1)
            else
               xig(n_m,nRstart-1)=xig(n_m,nRstart+1)!-two*dr*tops(l,m)
            end if
         end do
      end if

      !-- Handle Lower boundary
      !dr = r(n_r_max)-r(n_r_max-1)
      if ( nRstop == n_r_max ) then
         do n_m=n_m_start,n_m_stop
            if (kbotxi == 1) then ! Fixed xierature at bottom
               xig(n_m,nRstop+1)=two*xig(n_m,nRstop)-xig(n_m,nRstop-1)
            else
               xig(n_m,nRstop+1)=xig(n_m,nRstop-1)!+two*dr*bots(l,m)
            end if
         end do
      end if
      !$omp end parallel

   end subroutine fill_ghosts_xi
!---------------------------------------------------------------------------------
   subroutine update_xi_fd(xi, dxi, dxidt, tscheme)
      !
      ! This subroutine is called after the linear solves have been completed.
      ! This is then assembling the linear terms that will be used in the r.h.s.
      ! for the next iteration.
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      type(type_tarray), intent(inout) :: dxidt
      complex(cp),       intent(inout) :: xi(n_m_max,nRstart:nRstop) ! comp
      !-- Output: ds
      complex(cp),       intent(out) :: dxi(n_m_max,nRstart:nRstop) ! Radial derivative of xi

      !-- Local variables
      integer :: nR, n_m_start, n_m_stop, n_m

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dxidt)

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_xi_rhs_imp_ghost(xi_ghost, dxi, dxidt, 1, tscheme%l_imp_calc_rhs(1))
      else
         call get_xi_rhs_imp_ghost(xi_ghost, dxi, dxidt, tscheme%istage+1, &
              &                    tscheme%l_imp_calc_rhs(tscheme%istage+1))
      end if

      !$omp parallel default(shared) private(n_m_start,n_m_stop,nR,n_m)
      n_m_start=1; n_m_stop=n_m_max
      call get_openmp_blocks(n_m_start,n_m_stop)

      !-- Array copy from xi_ghost to xi
      do nR=nRstart,nRstop
         do n_m=n_m_start,n_m_stop
            xi(n_m,nR)=xi_ghost(n_m,nR)
         end do
      end do
      !$omp end parallel

   end subroutine update_xi_fd
!---------------------------------------------------------------------------------
   subroutine get_xi_rhs_imp_ghost(xig, dxi, dxidt, istage, l_calc_lin)
      !
      ! This subroutine computes the linear terms that enters the r.h.s.. This is
      ! used with R-distributed
      !

      !-- Input variables
      integer,     intent(in) :: istage
      logical,     intent(in) :: l_calc_lin
      complex(cp), intent(in) :: xig(n_m_max,nRstart-1:nRstop+1)

      !-- Output variables
      complex(cp),       intent(out) :: dxi(n_m_max,nRstart:nRstop)
      type(type_tarray), intent(inout) :: dxidt

      !-- Local variables
      complex(cp) :: work_Rloc(n_m_max,nRstart:nRstop)
      real(cp) :: dm2
      integer :: n_r, n_m, start_m, stop_m, m

      !$omp parallel default(shared)  private(start_m, stop_m, n_r, n_m, dm2, m)
      start_m=1; stop_m=n_m_max
      call get_openmp_blocks(start_m,stop_m)
      call get_ddr_ghost(xig, dxi, work_Rloc, n_m_max ,start_m, stop_m, &
           &             nRstart, nRstop, rscheme)

      if ( istage == 1 ) then
         do n_r=nRstart,nRstop
            do n_m=start_m,stop_m
               dxidt%old(n_m,n_r,istage)=xig(n_m,n_r)
            end do
         end do
      end if

      !-- Calculate implicit time step part:
      if ( l_calc_lin ) then
         do n_r=nRstart,nRstop
            do n_m=start_m,stop_m
               m = idx2m(n_m)
               dm2 = real(m*m,cp)
               dxidt%impl(n_m,n_r,istage)=XiDiffFac* (                     &
               &                                       work_Rloc(n_m,n_r)  &
               &                           +or1(n_r)*        dxi(n_m,n_r)  &
               &                     -dm2*or2(n_r)*          xig(n_m,n_r) )
            end do
         end do
      end if
      !$omp end parallel

   end subroutine get_xi_rhs_imp_ghost
!---------------------------------------------------------------------------------
   subroutine finish_exp_xi_Rdist(us_Rloc, dVsXiM, dxi_exp_last)
      !
      ! This subroutine completes the computation of the advection term by
      ! computing the radial derivative (R-distributed variant).
      !

      !-- Input variables
      complex(cp), intent(in) :: us_Rloc(n_m_max,nRstart:nRstop)
      complex(cp), intent(inout) :: dVsXiM(n_m_max,nRstart:nRstop)

      !-- Output variables
      complex(cp), intent(inout) :: dxi_exp_last(n_m_max,nRstart:nRstop)

      !-- Local variables
      complex(cp) :: work_Rloc(n_m_max,nRstart:nRstop)
      integer :: n_r, n_m, start_m, stop_m

      call get_dr_Rloc(dVsXiM, work_Rloc, n_m_max, nRstart, nRstop, n_r_max, &
           &           rscheme)

      !$omp parallel default(shared) private(n_r, n_m, start_m, stop_m)
      start_m=1; stop_m=n_m_max
      call get_openmp_blocks(start_m, stop_m)
      !$omp barrier

      !-- Finish calculation of the explicit part for current time step
      do n_r=nRstart,nRstop
         do n_m=start_m,stop_m
            dxi_exp_last(n_m,n_r)=dxi_exp_last(n_m,n_r)-or1(n_r)*work_Rloc(n_m,n_r)  &
            &                     -us_Rloc(n_m,n_r)*dxicond(n_r)
         end do
      end do
      !$omp end parallel

   end subroutine finish_exp_xi_Rdist
!---------------------------------------------------------------------------------
   subroutine assemble_xi_Rdist(xi, dxi, dxidt, tscheme)
      !
      ! This subroutine is used when an IMEX Runge-Kutta time scheme with an assembly
      ! stage is used. This is used when R is distributed.
      !

      !-- Input variable
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp),       intent(inout) :: xi(n_m_max,nRstart:nRstop)
      complex(cp),       intent(out) :: dxi(n_m_max,nRstart:nRstop)
      type(type_tarray), intent(inout) :: dxidt

      !-- Local variables
      integer :: m, n_m, n_r, start_m, stop_m
      complex(cp) :: work_Rloc(n_m_max,nRstart:nRstop)

      call tscheme%assemble_imex(work_Rloc, dxidt)

      !$omp parallel default(shared) private(start_m, stop_m, n_m, n_r, m)
      start_m=1; stop_m=n_m_max
      call get_openmp_blocks(start_m,stop_m)
      !$omp barrier

      do n_r=nRstart,nRstop
         do n_m=start_m,stop_m
            m = idx2m(n_m)
            if ( m == 0 ) then
               xi(n_m,n_r)=cmplx(real(work_Rloc(n_m,n_r)),0.0_cp,cp)
            else
               xi(n_m,n_r)=work_Rloc(n_m,n_r)
            end if
         end do
      end do

      if ( ktopxi == 1 .and. nRstart==1 ) then
         do n_m=start_m,stop_m
            xi(n_m,nRstart)=zero ! tops
         end do
      end if

      if ( kbotxi == 1 .and. nRstop==n_r_max ) then
         do n_m=start_m,stop_m
            xi(n_m,nRstop)=zero ! bots
         end do
      end if

      call bulk_to_ghost(xi, xi_ghost, 1, nRstart, nRstop, n_m_max, start_m, stop_m)
      !$omp end parallel

      call exch_ghosts(xi_ghost, n_m_max, nRstart, nRstop, 1)
      call fill_ghosts_xi(xi_ghost)

      !-- Finally call the construction of the implicit terms for the first stage
      !-- of next iteration
      call get_xi_rhs_imp_ghost(xi_ghost, dxi, dxidt, 1, tscheme%l_imp_calc_rhs(1))

   end subroutine assemble_xi_Rdist
!---------------------------------------------------------------------------------
   subroutine get_xiMat_Rdist(tscheme, xiMat)
      !
      !  This subroutine is used to construct the matrices when the parallel
      !  solver for F.D. is employed.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step

      !-- Output variables
      type(type_tri_par), intent(inout) :: xiMat

      !-- Local variables:
      real(cp) :: dm2
      integer :: nR, n_m, m

      !-- Bulk points: we fill all the points: this is then easier to handle
      !-- Neumann boundary conditions
      !$omp parallel default(shared) private(nR,m,n_m,dm2)
      !$omp do
      do nR=nRstart,nRstop
         do n_m=1,n_m_max
            m=idx2m(n_m)
            dm2 = real(m*m,cp)
            xiMat%diag(n_m,nR)=one - tscheme%wimp_lin(1)*XiDiffFac*   &
            &                (               rscheme%ddr(nR,1) +      &
            &                         or1(nR)*rscheme%dr(nR,1) -      &
            &                                        dm2*or2(nR)    )
            xiMat%low(n_m,nR)=   -tscheme%wimp_lin(1)*XiDiffFac*   &
            &                (               rscheme%ddr(nR,0) +   &
            &                         or1(nR)*rscheme%dr(nR,0) )
            xiMat%up(n_m,nR)=   -tscheme%wimp_lin(1)*XiDiffFac*    &
            &                (               rscheme%ddr(nR,2) +   &
            &                         or1(nR)*rscheme%dr(nR,2) )
         end do
      end do
      !$omp end do

      !----- Boundary conditions:
      if ( nRstart == 1 ) then
         !$omp do
         do n_m=1,n_m_max
            if ( ktopxi == 1 ) then
               xiMat%diag(n_m,1)=one
               xiMat%up(n_m,1)  =0.0_cp
               xiMat%low(n_m,1) =0.0_cp
            else
               xiMat%up(n_m,1)=xiMat%up(n_m,1)+xiMat%low(n_m,1)
               !fd_fac_top(n_m)=two*(r(2)-r(1))*xiMat%low(n_m,1)
            end if
         end do
         !$omp end do
      end if

      if ( nRstop == n_r_max ) then
         !$omp do
         do n_m=1,n_m_max
            if ( kbotxi == 1 ) then
               xiMat%diag(n_m,n_r_max)=one
               xiMat%up(n_m,n_r_max)  =0.0_cp
               xiMat%low(n_m,n_r_max) =0.0_cp
            else
               xiMat%low(n_m,n_r_max)=xiMat%up(n_m,n_r_max)+xiMat%low(n_m,n_r_max)
               !fd_fac_bot(n_m)=two*(r(n_r_max-1)-r(n_r_max))*xiMat%up(n_m,n_r_max)
            end if
         end do
         !$omp end do
      end if
      !$omp end parallel

      !-- LU decomposition:
      call xiMat%prepare_mat()

   end subroutine get_xiMat_Rdist
!---------------------------------------------------------------------------------
end module update_xi_fd_mod
