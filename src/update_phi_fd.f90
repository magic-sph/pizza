module update_phi_fd_mod

   use precision_mod
   use parallel_mod, only: get_openmp_blocks
   use mem_alloc, only: bytes_allocated
   use constants, only: zero, one, two
   use truncation, only: n_m_max, n_r_max, idx2m
   use blocking, only: nRstart, nRstop
   use parallel_solvers, only: type_tri_par
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use namelists, only: phaseDiffFac, kbotphi, ktopphi, l_full_disk, stef, pr, &
       &                phi_bot, phi_top
   use radial_functions, only: or1, or2, rscheme
   use radial_der, only: get_ddr_ghost, get_dr_Rloc, exch_ghosts, bulk_to_ghost

   implicit none

   private

   type(type_tri_par), public :: phiMat_FD
   complex(cp), public, allocatable :: phi_ghost(:,:)
   logical, public, allocatable :: lPhimat_FD(:)

   public :: initialize_phi_fd, finalize_phi_fd, prepare_phi_fd, fill_ghosts_phi, &
   &         update_phi_fd, get_phi_rhs_imp_ghost, assemble_phi_Rdist

contains

   subroutine initialize_phi_fd

      call phiMat_FD%initialize(nRstart,nRstop,1,n_m_max)
      allocate( phi_ghost(n_m_max,nRstart-1:nRstop+1) )
      bytes_allocated=bytes_allocated + n_m_max*(nRstop-nRstart+3)*SIZEOF_DEF_COMPLEX
      phi_ghost(:,:)=zero
      allocate( lPhimat_FD(1:n_m_max) )

   end subroutine initialize_phi_fd
!---------------------------------------------------------------------------------
   subroutine finalize_phi_fd

      deallocate( lPhimat_FD, phi_ghost)
      call phiMat_FD%finalize()

   end subroutine finalize_phi_fd
!---------------------------------------------------------------------------------
   subroutine prepare_phi_fd(tscheme, dphidt)
      !
      ! This subroutine is used to assemble the r.h.s. of the composition equation
      ! when parallel F.D solvers are used. Boundary values are set here.
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      type(type_tarray), intent(inout) :: dphidt

      !-- Local variables:
      integer :: n_m_start, n_m_stop, nR, n_m, m

      if ( .not. lPhimat_FD(1) ) then
         call get_phiMat_Rdist(tscheme, phiMat_FD)
         lPhimat_FD(:)=.true.
      end if

      !$omp parallel default(shared) private(n_m_start,n_m_stop, nR, n_m, m)
      n_m_start=1; n_m_stop=n_m_max
      call get_openmp_blocks(n_m_start,n_m_stop)
      !$omp barrier

      !-- Now assemble the right hand side
      call tscheme%set_imex_rhs_ghost(phi_ghost, dphidt, n_m_start, n_m_stop, 1)

      !-- Set boundary conditions
      if ( nRstart == 1 ) then
         nR=1
         do n_m=n_m_start,n_m_stop
            m = idx2m(n_m)
            if ( ktopphi == 1 ) then ! Fixed phase field
               if ( m == 0 ) then
                  phi_ghost(n_m,nR)=cmplx(phi_top,0.0_cp,cp)
               else
                  phi_ghost(n_m,nR)=zero
               end if
            end if
            phi_ghost(n_m,nR-1)=zero ! Set ghost zone to zero
         end do
      end if

      if ( nRstop == n_r_max ) then
         nR=n_r_max
         do n_m=n_m_start,n_m_stop
            m = idx2m(n_m)
            if ( l_full_disk ) then
               if ( m > 0 ) then
                  phi_ghost(n_m,nR)=zero
               end if
            else
               if ( kbotphi == 1 ) then ! Fixed phase field
                  if ( m == 0 ) then
                     phi_ghost(n_m,nR)=cmplx(phi_bot,0.0_cp,cp)
                  else
                     phi_ghost(n_m,nR)=zero
                  end if
               end if
            end if
            phi_ghost(n_m,nR+1)=zero ! Set ghost zone to zero
         end do
      end if
      !$omp end parallel

   end subroutine prepare_phi_fd
!---------------------------------------------------------------------------------
   subroutine fill_ghosts_phi(phig)
      !
      ! This subroutine is used to fill the ghosts zones that are located at
      ! nR=n_r_cmb-1 and nR=n_r_icb+1. This is used to properly set the Neuman
      ! boundary conditions. In case Dirichlet BCs are used, a simple first order
      ! extrapolation is employed. This is anyway only used for outputs (like Nusselt
      ! numbers).
      !
      complex(cp), intent(inout) :: phig(n_m_max,nRstart-1:nRstop+1)

      !-- Local variables
      integer :: n_m, n_m_start, n_m_stop, m

      !$omp parallel default(shared) private(n_m_start, n_m_stop, n_m, m)
      n_m_start=1; n_m_stop=n_m_max
      call get_openmp_blocks(n_m_start,n_m_stop)
      !$omp barrier

      !-- Handle upper boundary
      !dr = r(2)-r(1)
      if ( nRstart == 1 ) then
         do n_m=n_m_start,n_m_stop
            if ( ktopphi == 1 ) then
               phig(n_m,nRstart-1)=two*phig(n_m,nRstart)-phig(n_m,nRstart+1)
            else
               phig(n_m,nRstart-1)=phig(n_m,nRstart+1)!-two*dr*tops(l,m)
            end if
         end do
      end if

      !-- Handle Lower boundary
      !dr = r(n_r_max)-r(n_r_max-1)
      if ( nRstop == n_r_max ) then
         do n_m=n_m_start,n_m_stop
            if ( l_full_disk ) then
               m = idx2m(n_m)
               if ( m == 0 ) then
                  phig(n_m,nRstop+1)=phig(n_m,nRstop-1)!+two*dr*bots(l,m)
               else
                  phig(n_m,nRstop+1)=two*phig(n_m,nRstop)-phig(n_m,nRstop-1)
               end if
            else
               if (kbotphi == 1) then ! Fixed phase field at bottom
                  phig(n_m,nRstop+1)=two*phig(n_m,nRstop)-phig(n_m,nRstop-1)
               else
                  phig(n_m,nRstop+1)=phig(n_m,nRstop-1)!+two*dr*bots(l,m)
               end if
            end if
         end do
      end if
      !$omp end parallel

   end subroutine fill_ghosts_phi
!---------------------------------------------------------------------------------
   subroutine update_phi_fd(phi, dphidt, tscheme)
      !
      ! This subroutine is called after the linear solves have been completed.
      ! This is then assembling the linear terms that will be used in the r.h.s.
      ! for the next iteration.
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      type(type_tarray), intent(inout) :: dphidt
      complex(cp),       intent(inout) :: phi(n_m_max,nRstart:nRstop) ! comp

      !-- Local variables
      integer :: nR, n_m_start, n_m_stop, n_m

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dphidt)

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_phi_rhs_imp_ghost(phi_ghost, dphidt, 1, tscheme%l_imp_calc_rhs(1))
      else
         call get_phi_rhs_imp_ghost(phi_ghost, dphidt, tscheme%istage+1, &
              &                     tscheme%l_imp_calc_rhs(tscheme%istage+1))
      end if

      !$omp parallel default(shared) private(n_m_start,n_m_stop,nR,n_m)
      n_m_start=1; n_m_stop=n_m_max
      call get_openmp_blocks(n_m_start,n_m_stop)

      !-- Array copy from phi_ghost to phi
      do nR=nRstart,nRstop
         do n_m=n_m_start,n_m_stop
            phi(n_m,nR)=phi_ghost(n_m,nR)
         end do
      end do
      !$omp end parallel

   end subroutine update_phi_fd
!---------------------------------------------------------------------------------
   subroutine get_phi_rhs_imp_ghost(phig, dphidt, istage, l_calc_lin)
      !
      ! This subroutine computes the linear terms that enters the r.h.s.. This is
      ! used with R-distributed
      !

      !-- Input variables
      integer,     intent(in) :: istage
      logical,     intent(in) :: l_calc_lin
      complex(cp), intent(in) :: phig(n_m_max,nRstart-1:nRstop+1)

      !-- Output variables
      type(type_tarray), intent(inout) :: dphidt

      !-- Local variables
      complex(cp) :: dphi(n_m_max,nRstart:nRstop)
      complex(cp) :: work_Rloc(n_m_max,nRstart:nRstop)
      real(cp) :: dm2
      integer :: n_r, n_m, start_m, stop_m, m

      !$omp parallel default(shared)  private(start_m, stop_m, n_r, n_m, dm2, m)
      start_m=1; stop_m=n_m_max
      call get_openmp_blocks(start_m,stop_m)
      call get_ddr_ghost(phig, dphi, work_Rloc, n_m_max ,start_m, stop_m, &
           &             nRstart, nRstop, rscheme)

      if ( istage == 1 ) then
         do n_r=nRstart,nRstop
            do n_m=start_m,stop_m
               dphidt%old(n_m,n_r,istage)=5.0_cp/6.0_cp*stef*pr*phig(n_m,n_r)
            end do
         end do
      end if

      !-- Calculate implicit time step part:
      if ( l_calc_lin ) then
         do n_r=nRstart,nRstop
            do n_m=start_m,stop_m
               m = idx2m(n_m)
               dm2 = real(m*m,cp)
               dphidt%impl(n_m,n_r,istage)=PhaseDiffFac*(work_Rloc(n_m,n_r)  &
               &                           +or1(n_r)*         dphi(n_m,n_r)  &
               &                     -dm2*or2(n_r)*           phig(n_m,n_r) )
               if ( m == 0 .and. l_full_disk .and. n_r==n_r_max ) then
                  dphidt%impl(n_m,n_r,istage)=PhaseDiffFac*two*work_Rloc(n_m,n_r)
               end if
            end do
         end do
      end if
      !$omp end parallel

   end subroutine get_phi_rhs_imp_ghost
!---------------------------------------------------------------------------------
   subroutine assemble_phi_Rdist(phi, dphidt, tscheme)
      !
      ! This subroutine is used when an IMEX Runge-Kutta time scheme with an assembly
      ! stage is used. This is used when R is distributed.
      !

      !-- Input variable
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp),       intent(inout) :: phi(n_m_max,nRstart:nRstop)
      type(type_tarray), intent(inout) :: dphidt

      !-- Local variables
      integer :: m, n_m, n_r, start_m, stop_m
      complex(cp) :: work_Rloc(n_m_max,nRstart:nRstop)

      call tscheme%assemble_imex(work_Rloc, dphidt)

      !$omp parallel default(shared) private(start_m, stop_m, n_m, n_r, m)
      start_m=1; stop_m=n_m_max
      call get_openmp_blocks(start_m,stop_m)
      !$omp barrier

      do n_r=nRstart,nRstop
         do n_m=start_m,stop_m
            m = idx2m(n_m)
            if ( m == 0 ) then
               phi(n_m,n_r)=6.0_cp/5.0_cp*cmplx(real(work_Rloc(n_m,n_r)),0.0_cp,cp) &
               &            /stef/pr
            else
               phi(n_m,n_r)=6.0_cp/5.0_cp*work_Rloc(n_m,n_r)/stef/pr
            end if
         end do
      end do

      if ( ktopphi == 1 .and. nRstart==1 ) then
         do n_m=start_m,stop_m
            m = idx2m(n_m)
            if ( m == 0 ) then
               phi(n_m,nRstart)=cmplx(phi_top,0.0_cp,cp)
            else
               phi(n_m,nRstart)=zero
            end if
         end do
      end if

      if ( kbotphi == 1 .and. nRstop==n_r_max ) then
         do n_m=start_m,stop_m
            m = idx2m(n_m)
            if ( m == 0 ) then
               phi(n_m,nRstop)=cmplx(phi_bot,0.0_cp,cp)
            else
               phi(n_m,nRstop)=zero
            end if
         end do
      end if

      call bulk_to_ghost(phi, phi_ghost, 1, nRstart, nRstop, n_m_max, start_m, stop_m)
      !$omp end parallel

      call exch_ghosts(phi_ghost, n_m_max, nRstart, nRstop, 1)
      call fill_ghosts_phi(phi_ghost)

      !-- Finally call the construction of the implicit terms for the first stage
      !-- of next iteration
      call get_phi_rhs_imp_ghost(phi_ghost, dphidt, 1, tscheme%l_imp_calc_rhs(1))

   end subroutine assemble_phi_Rdist
!---------------------------------------------------------------------------------
   subroutine get_phiMat_Rdist(tscheme, phiMat)
      !
      !  This subroutine is used to construct the matrices when the parallel
      !  solver for F.D. is employed.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step

      !-- Output variables
      type(type_tri_par), intent(inout) :: phiMat

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
            phiMat%diag(n_m,nR)=5.0_cp/6.0_cp*stef*pr- tscheme%wimp_lin(1)* &
            &                   phaseDiffFac*(     rscheme%ddr(nR,1) +      &
            &                               or1(nR)*rscheme%dr(nR,1) -      &
            &                                            dm2*or2(nR)    )
            phiMat%low(n_m,nR)=   -tscheme%wimp_lin(1)*phaseDiffFac*   &
            &                    (               rscheme%ddr(nR,0) +   &
            &                             or1(nR)*rscheme%dr(nR,0) )
            phiMat%up(n_m,nR)=   -tscheme%wimp_lin(1)*phaseDiffFac*    &
            &                    (               rscheme%ddr(nR,2) +   &
            &                             or1(nR)*rscheme%dr(nR,2) )
            if ( m == 0 .and. l_full_disk .and. nR==n_r_max ) then
               phiMat%diag(n_m,nR)=5.0_cp/6.0_cp*stef*pr- tscheme%wimp_lin(1)* &
               &                   phaseDiffFac*two*rscheme%ddr(nR,1)
               phiMat%low(n_m,nR)=   -tscheme%wimp_lin(1)*phaseDiffFac*   &
               &                               two*rscheme%ddr(nR,0)
               phiMat%up(n_m,nR)=   -tscheme%wimp_lin(1)*phaseDiffFac*    &
               &                               two*rscheme%ddr(nR,2)
            end if
         end do
      end do
      !$omp end do

      !----- Boundary conditions:
      if ( nRstart == 1 ) then
         !$omp do
         do n_m=1,n_m_max
            if ( ktopphi == 1 ) then
               phiMat%diag(n_m,1)=one
               phiMat%up(n_m,1)  =0.0_cp
               phiMat%low(n_m,1) =0.0_cp
            else
               phiMat%up(n_m,1)=phiMat%up(n_m,1)+phiMat%low(n_m,1)
               !fd_fac_top(n_m)=two*(r(2)-r(1))*phiMat%low(n_m,1)
            end if
         end do
         !$omp end do
      end if

      if ( nRstop == n_r_max ) then
         !$omp do
         do n_m=1,n_m_max
            if ( l_full_disk ) then
               m = idx2m(n_m)
               if ( m == 0 ) then
                  phiMat%low(n_m,n_r_max)=phiMat%up(n_m,n_r_max)+phiMat%low(n_m,n_r_max)
                  !fd_fac_bot(n_m)=two*(r(n_r_max-1)-r(n_r_max))*phiMat%up(n_m,n_r_max)
               else
                  phiMat%diag(n_m,n_r_max)=one
                  phiMat%up(n_m,n_r_max)  =0.0_cp
                  phiMat%low(n_m,n_r_max) =0.0_cp
               end if
            else
               if ( kbotphi == 1 ) then
                  phiMat%diag(n_m,n_r_max)=one
                  phiMat%up(n_m,n_r_max)  =0.0_cp
                  phiMat%low(n_m,n_r_max) =0.0_cp
               else
                  phiMat%low(n_m,n_r_max)=phiMat%up(n_m,n_r_max)+phiMat%low(n_m,n_r_max)
                  !fd_fac_bot(n_m)=two*(r(n_r_max-1)-r(n_r_max))*phiMat%up(n_m,n_r_max)
               end if
            end if
         end do
         !$omp end do
      end if
      !$omp end parallel

      !-- LU decomposition:
      call phiMat%prepare_mat()

   end subroutine get_phiMat_Rdist
!---------------------------------------------------------------------------------
end module update_phi_fd_mod
