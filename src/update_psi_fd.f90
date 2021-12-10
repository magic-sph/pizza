module update_psi_fd_mod

   use precision_mod
   use parallel_mod
   use mem_alloc, only: bytes_allocated
   use useful, only: abortRun
   use constants, only: zero, ci, one, two, three, four, half
   use truncation, only: n_m_max, n_r_max, idx2m
   use blocking, only: nRstart, nRstop
   use parallel_solvers, only: type_tri_par, type_penta_par
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use namelists, only: kbotv, ktopv, BuoFac, ChemFac, l_chem, l_heat, CorFac, &
       &                l_coriolis_imp, ViscFac, r_cmb, l_ek_pump, l_non_rot
   use radial_functions, only: or1, rgrav, rscheme, beta, ekpump, or2, r, dbeta, &
       &                       d2beta, d3beta, or3, ekp_up, ekp_us, ekp_dusdp
   use radial_der, only: get_dr_Rloc, get_ddddr_ghost, get_ddr_ghost, exch_ghosts, &
       &                 bulk_to_ghost
   use vort_balance, only: vort_bal_type
   use vp_balance, only: vp_bal_type

   implicit none

   private

   type(type_penta_par), public :: psiMat_FD
   type(type_tri_par), public :: upMat_FD, ellMat_FD
   complex(cp), public, allocatable :: psi_ghost(:,:)
   complex(cp), public, allocatable :: up0_ghost(:)
   logical, public, allocatable :: lPsimat_FD(:), lEllmat_FD(:)

   public :: initialize_psi_fd, finalize_psi_fd, finish_exp_psi_Rdist, &
   &         prepare_psi_fd, fill_ghosts_psi, update_psi_FD,           &
   &         get_psi_rhs_imp_ghost, assemble_psi_Rloc

contains

   subroutine initialize_psi_fd(tscheme)

      !-- Input variable
      class(type_tscheme), intent(in) :: tscheme ! time scheme

      call upMat_FD%initialize(nRstart,nRstop,1,1)
      call psiMat_FD%initialize(nRstart,nRstop,1,n_m_max)
      allocate( psi_ghost(n_m_max,nRstart-2:nRstop+2) )
      allocate( up0_ghost(nRstart-1:nRstop+1) )
      bytes_allocated=bytes_allocated + n_m_max*(nRstop-nRstart+8)*SIZEOF_DEF_COMPLEX
      psi_ghost(:,:)=zero
      up0_ghost(:)=zero
      allocate( lPsimat_FD(1:n_m_max) )

      if ( tscheme%l_assembly ) then
         call ellMat_FD%initialize(nRstart,nRstop,1,n_m_max)
         allocate( lEllmat_FD(1:n_m_max) )
         lEllmat_FD(:)=.false.
      end if

   end subroutine initialize_psi_fd
!---------------------------------------------------------------------------------
   subroutine finalize_psi_fd(tscheme)

      !-- Input variable
      class(type_tscheme), intent(in) :: tscheme ! time scheme

      deallocate( lPsimat_FD, psi_ghost, up0_ghost)
      call psiMat_FD%finalize()
      call upMat_FD%finalize()
      if ( tscheme%l_assembly ) then
         deallocate( lEllmat_FD )
         call ellMat_FD%finalize()
      end if

   end subroutine finalize_psi_fd
!---------------------------------------------------------------------------------
   subroutine prepare_psi_fd(tscheme, dpsidt, lu_time, n_lu_calls)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      type(type_tarray), intent(inout) :: dpsidt
      real(cp),          intent(inout) :: lu_time
      integer,           intent(inout) :: n_lu_calls


      !-- Local variables
      integer :: n_m, m, n_m_start, n_m_stop, nR

      if ( .not. lPsimat_FD(1) ) then
         call get_psiMat_Rdist(tscheme, psiMat_FD, lu_time, n_lu_calls)
         call get_uphiMat_Rdist(tscheme, upMat_FD)
         lPsimat_FD(:)=.true.
      end if

      !$omp parallel default(shared) private(n_m_start,n_m_stop, nR, m, n_m)
      n_m_start=1; n_m_stop=n_m_max
      call get_openmp_blocks(n_m_start,n_m_stop)
      !$omp barrier

      call tscheme%set_imex_rhs_ghost(psi_ghost, dpsidt, n_m_start, n_m_stop, 2)

      do nR=nRstart,nRstop ! copy to fill up0 ! Does it work fine with OpenMP here?
         up0_ghost(nR)=psi_ghost(1,nR)
      end do

      !-- Set boundary conditions
      if ( nRstart == 1 ) then
         nR=1
         do n_m=n_m_start,n_m_stop
            m=idx2m(n_m)
            if ( m == 0 ) then
               if ( ktopv /= 1 ) up0_ghost(nR)=zero
               up0_ghost(nR-1)=zero 
            else
               psi_ghost(n_m,nR)  =zero ! Non-penetration condition
               psi_ghost(n_m,nR-1)=zero ! Ghost zones set to zero
               psi_ghost(n_m,nR-2)=zero
            end if
         end do
      end if

      if ( nRstop == n_r_max ) then
         nR=n_r_max
         do n_m=n_m_start,n_m_stop
            m=idx2m(n_m)
            if ( m == 0 ) then
               if ( kbotv /= 1 ) up0_ghost(nR)=zero
               up0_ghost(nR+1)=zero ! Set ghost zone to zero
            else
               psi_ghost(n_m,nR)  =zero ! Non-penetration condition
               psi_ghost(n_m,nR+1)=zero ! Ghost zones set to zero
               psi_ghost(n_m,nR+2)=zero
            end if
         end do
      end if
      !$omp end parallel

   end subroutine prepare_psi_fd
!---------------------------------------------------------------------------------
   subroutine fill_ghosts_psi(psig, up0g)

      !-- Input variables
      complex(cp), intent(inout) :: up0g(nRstart-1:nRstop+1)
      complex(cp), intent(inout) :: psig(n_m_max, nRstart-2:nRstop+2)

      !-- Local variables
      real(cp) :: dr
      integer :: n_m_start, n_m_stop, n_m, m

      !-- Boundary conditions for uphi0
      if ( nRstart == 1 ) then
         dr = r(2)-r(1)
         if ( ktopv == 1 ) then
            up0g(nRstart-1)=up0g(nRstart+1)-two*dr*or1(1)*up0g(nRstart)
         else
            up0g(nRstart-1)=two*up0g(nRstart)-up0g(nRstart+1)
         end if
      end if

      if ( nRstop == n_r_max ) then
         dr = r(n_r_max)-r(n_r_max-1)
         if ( kbotv == 1 ) then
            up0g(nRstop+1)=up0g(nRstop-1)+two*dr*or1(n_r_max)*up0g(nRstop)
         else
            up0g(nRstop+1)=two*up0g(nRstop)-up0g(nRstop-1)
         end if
      end if

      !$omp parallel default(shared) private(n_m_start, n_m_stop, n_m, m)
      n_m_start=1; n_m_stop=n_m_max
      call get_openmp_blocks(n_m_start,n_m_stop)
      !$omp barrier

      !-- Upper boundary
      dr = r(2)-r(1)
      if ( nRstart == 1 ) then ! Rank with n_r_mcb
         do n_m=n_m_start,n_m_stop
            m = idx2m(n_m)
            if ( m == 0 ) cycle
            if ( ktopv == 1 ) then  ! Stress-free
               psig(n_m,nRstart-1)=-(one-half*or1(1)*dr)/ &
               &                    (one+half*or1(1)*dr) * psig(n_m,nRstart+1)
            else ! Rigid boundary condition
               psig(n_m,nRstart-1)=psig(n_m,nRstart+1) ! dpsi/dr=0
            end if
            psig(n_m,nRstart-2)=zero
         end do
      end if

      !-- Lower boundary
      dr = r(n_r_max)-r(n_r_max-1)
      if ( nRstop == n_r_max ) then
         do n_m=n_m_start,n_m_stop
            m = idx2m(n_m)
            if ( m == 0 ) cycle
            if ( kbotv == 1 ) then ! Stress-free
               psig(n_m,nRstop+1)=-(one+half*or1(n_r_max)*dr)/ &
               &                (one-half*or1(n_r_max)*dr) * psig(n_m,nRstop-1)
            else
               psig(n_m,nRstop+1)=psig(n_m,nRstop-1) ! dw=0
            end if
            psig(n_m,nRstop+2)=zero
         end do
      end if
      !$omp end parallel

   end subroutine fill_ghosts_psi
!---------------------------------------------------------------------------------
   subroutine update_psi_FD(us_Rloc, up_Rloc, om_Rloc, dpsidt, tscheme, &
              &             vp_bal, vort_bal)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp),         intent(inout) :: us_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: up_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: om_Rloc(n_m_max,nRstart:nRstop)
      type(type_tarray),   intent(inout) :: dpsidt
      type(vp_bal_type),   intent(inout) :: vp_bal
      type(vort_bal_type), intent(inout) :: vort_bal

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dpsidt)

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_psi_rhs_imp_ghost(psi_ghost, up0_ghost, us_Rloc, up_Rloc, &
              &                     om_Rloc, dpsidt, 1, vp_bal, vort_bal,   &
              &                     tscheme%l_imp_calc_rhs(1))
      else
         call get_psi_rhs_imp_ghost(psi_ghost, up0_ghost, us_Rloc, up_Rloc, om_Rloc, &
              &                     dpsidt, tscheme%istage+1, vp_bal, vort_bal,      &
              &                     tscheme%l_imp_calc_rhs(tscheme%istage+1))
      end if

   end subroutine update_psi_FD
!---------------------------------------------------------------------------------
   subroutine get_psi_rhs_imp_ghost(psig, up0g, us_Rloc, up_Rloc, om_Rloc,     &
              &                     dpsidt, istage, vp_bal, vort_bal, l_calc_lin)

      !-- Input variables
      integer,     intent(in) :: istage
      logical,     intent(in) :: l_calc_lin
      complex(cp), intent(in) :: up0g(nRstart-1:nRstop+1)
      complex(cp), intent(in) :: psig(n_m_max,nRstart-2:nRstop+2)

      !-- Output variables
      type(type_tarray),   intent(inout) :: dpsidt
      complex(cp),         intent(inout) :: us_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: up_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: om_Rloc(n_m_max,nRstart:nRstop)
      type(vp_bal_type),   intent(inout) :: vp_bal
      type(vort_bal_type), intent(inout) :: vort_bal

      !-- Local variables
      complex(cp) :: dup0(nRstart:nRstop), ddup0(nRstart:nRstop)
      complex(cp) :: work_Rloc(n_m_max,nRstart:nRstop)
      complex(cp) :: dddpsi_Rloc(n_m_max,nRstart:nRstop)
      real(cp) :: dm2
      integer :: n_m_start, n_m_stop, n_m, nR, m

      call get_ddr_ghost(up0g, dup0, ddup0, 1, 1, 1, nRstart, nRstop, rscheme)

      !$omp parallel default(shared)  private(n_m_start, n_m_stop, nR, n_m, m, dm2)
      n_m_start=1; n_m_stop=n_m_max
      call get_openmp_blocks(n_m_start,n_m_stop)
      call get_ddddr_ghost(psig, up_Rloc, om_Rloc, dddpsi_Rloc, work_Rloc, n_m_max, &
           &               n_m_start, n_m_stop, nRstart, nRstop, rscheme)
      !$omp barrier
      !-- Warning: at this stage, us, up and om are used as temp arrays to store the
      !-- derivatives of psi 

      if ( istage == 1 ) then
         do nR=nRstart,nRstop
            do n_m=n_m_start,n_m_stop
               m = idx2m(n_m)
               dm2 = real(m*m,cp)
               if ( m == 0 ) then
                  dpsidt%old(n_m,nR,istage)=up0g(nR)
               else
                  dpsidt%old(n_m,nR,istage)=-(om_Rloc(n_m,nR)+(or1(nR)+beta(nR))* &
                  &                           up_Rloc(n_m,nR)+(or1(nR)*beta(nR)+  &
                  &                           dbeta(nR)-dm2*or2(nR))*psig(n_m,nR))
               end if
            end do
         end do
      end if

      if ( l_calc_lin  .or. vp_bal%l_calc .or. vort_bal%l_calc ) then
         do nR=nRstart,nRstop
            do n_m=n_m_start,n_m_stop
               m = idx2m(n_m)
               dm2 = real(m*m,cp)
               if ( m == 0 ) then
                  dpsidt%impl(n_m,nR,istage)=ViscFac*(ddup0(nR)+or1(nR)*dup0(nR)-&
                  &                                   or2(nR)*up0g(nR)) -        &
                  &                          CorFac*ekpump(nR)*up0g(nR)

                  if ( vp_bal%l_calc ) then
                     vp_bal%visc(nR)=ViscFac*(real(ddup0(nR))+or1(nR)*real(dup0(nR))-&
                     &                        or2(nR)*real(up0g(nR)))
                     if ( l_ek_pump ) then
                        vp_bal%pump(nR)=-CorFac*ekpump(nR)*real(up0g(nR))
                     end if
                  end if
               else
                  dpsidt%impl(n_m,nR,istage)=-ViscFac* (         work_Rloc(n_m, nR)+ &
                  &        ( beta(nR)+two*or1(nR) )*            dddpsi_Rloc(n_m,nR)+ &
                  &        ( three*dbeta(nR)+two*beta(nR)*or1(nR)-or2(nR)-two*dm2*   &
                  &          or2(nR) ) *                            om_Rloc(n_m,nR)+ &
                  &        ( three*d2beta(nR)+four*dbeta(nR)*or1(nR)-beta(nR)*       &
                  &          or2(nR)+or3(nR)+dm2*or2(nR)*(two*or1(nR)-beta(nR)) ) *  &
                  &                                                 up_Rloc(n_m,nR)+ &
                  &        ( d3beta(nR)+two*d2beta(nR)*or1(nR)-dbeta(nR)*or2(nR) +   &
                  &          beta(nR)*or3(nR)-dm2*or2(nR)*(four*or2(nR)+beta(nR)*    &
                  &          or1(nR)+dbeta(nR)-dm2*or2(nR) ) ) *       psig(n_m,nR) )

                  if ( l_coriolis_imp .and. (.not. l_non_rot) ) then
                     dpsidt%impl(n_m,nR,istage)=dpsidt%impl(n_m,nR,istage)+CorFac* &
                     &                          beta(nR)*ci*real(m,cp)*or1(nR)*    &
                     &                          psig(n_m,nR)
                  end if

                  if ( vort_bal%l_calc ) then
                     vort_bal%visc(n_m,nR)=dpsidt%impl(n_m,nR,istage)
                     if ( .not. l_non_rot ) then
                        vort_bal%cor(n_m,nR)=CorFac*beta(nR)*ci*real(m,cp)*or1(nR)* &
                        &                    psig(n_m,nR)
                     end if
                  end if

                  if ( l_ek_pump ) then
                     !-- Careful again that om is 2nd derivative of psi, up first
                     ! derivative
                     dpsidt%impl(n_m,nR,istage)=dpsidt%impl(n_m,nR,istage)+CorFac*   &
                     &                    ekpump(nR) * (          om_Rloc(n_m,nR) +  &
                     &            (or1(nR)+beta(nR)+ekp_up(nR))  *up_Rloc(n_m,nR) +  &
                     &  (dbeta(nR)+beta(nR)*or1(nR)-dm2*or2(nR)+beta(nR)*ekp_up(nR)- &
                     &          ci*real(m,cp)*or1(nR)*(ekp_us(nR)+ci*real(m,cp)*     &
                     &          ekp_dusdp(nR)))*                     psig(n_m,nR))
                     if ( vort_bal%l_calc ) then
                        vort_bal%pump(n_m,nR)=CorFac*ekpump(nR)* (om_Rloc(n_m,nR) +  &
                        &         (or1(nR)+beta(nR)+ekp_up(nR))  *up_Rloc(n_m,nR) +  &
                        &      (dbeta(nR)+beta(nR)*or1(nR)-dm2+beta(nR)*ekp_up(nR)-  &
                        &       ci*real(m,cp)*or1(nR)*(ekp_us(nR)+ci*real(m,cp)*     &
                        &       ekp_dusdp(nR)))*                     psig(n_m,nR))
                     end if
                  end if
               end if
            end do
         end do
      end if

      !-- Now properly set us, up and om
      do nR=nRstart,nRstop
         do n_m=n_m_start, n_m_stop
            m = idx2m(n_m)
            dm2 = real(m*m,cp)
            if (m == 0) then
               us_Rloc(n_m,nR)=zero
               up_Rloc(n_m,nR)=up0g(nR)
               om_Rloc(n_m,nR)=dup0(nR)+or1(nR)*up0g(nR)
            else
               us_Rloc(n_m,nR)=ci*real(m,cp)*or1(nR)*psig(n_m,nR)
               om_Rloc(n_m,nR)=-(om_Rloc(n_m,nR)+(or1(nR)+beta(nR))*up_Rloc(n_m,nR)+&
               &                 (beta(nR)*or1(nR)+dbeta(nR)-dm2*or2(nR))*psig(n_m,nR))
               up_Rloc(n_m,nR)=-up_Rloc(n_m,nR)-beta(nR)*psig(n_m,nR)
            end if
         end do
      end do

      !$omp end parallel

   end subroutine get_psi_rhs_imp_ghost
!---------------------------------------------------------------------------------
   subroutine finish_exp_psi_Rdist(us_Rloc, dVsOmM, temp_Rloc, xi_Rloc, &
              &                    dpsi_exp_last, vp_bal, vort_bal, tscheme)
      !
      ! This subroutine completes the computation of the advection term by
      ! computing the radial derivative (R-distributed variant).
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      complex(cp),         intent(in) :: us_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: dVsOmM(n_m_max,nRstart:nRstop)
      complex(cp),         intent(in) :: temp_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(in) :: xi_Rloc(n_m_max,nRstart:nRstop)

      !-- Output variables
      complex(cp),         intent(inout) :: dpsi_exp_last(n_m_max,nRstart:nRstop)
      type(vp_bal_type),   intent(inout) :: vp_bal
      type(vort_bal_type), intent(inout) :: vort_bal

      !-- Local variables
      complex(cp) :: work_Rloc(n_m_max,nRstart:nRstop)
      integer :: n_r, n_m, m, start_m, stop_m

      call get_dr_Rloc(dVsOmM, work_Rloc, n_m_max, nRstart, nRstop, n_r_max, &
           &           rscheme)

      !$omp parallel default(shared) private(n_r, n_m, m, start_m, stop_m)
      start_m=1; stop_m=n_m_max
      call get_openmp_blocks(start_m, stop_m)
      !$omp barrier

      !-- Finish calculation of the explicit part for current time step
      do n_r=nRstart,nRstop
         do n_m=start_m,stop_m
            m = idx2m(n_m)
            if ( m == 0 ) then
               if ( vp_bal%l_calc .and. tscheme%istage==1 ) then
                  vp_bal%rey_stress(n_r)=real(dpsi_exp_last(n_m,n_r))
               end if
            else
               dpsi_exp_last(n_m,n_r)=dpsi_exp_last(n_m,n_r)-or1(n_r)*work_Rloc(n_m,n_r)
               if ( vort_bal%l_calc .and. tscheme%istage==1 ) then
                  vort_bal%adv(n_m,n_r)=dpsi_exp_last(n_m,n_r)
               end if
               if ( .not. l_coriolis_imp ) then
                  dpsi_exp_last(n_m,n_r)=       dpsi_exp_last(n_m,n_r) &
                  &                 +CorFac*beta(n_r)*us_Rloc(n_m,n_r)
               end if
               !-- Buoyancy is handled explicitly
               if ( l_heat ) then
                  dpsi_exp_last(n_m,n_r)=dpsi_exp_last(n_m,n_r)-BuoFac*rgrav(n_r)*&
                  &                      or1(n_r)*ci*real(m,cp)*temp_Rloc(n_m,n_r)
                  if ( vort_bal%l_calc .and. tscheme%istage==1 ) then
                     vort_bal%buo(n_m,n_r)=-BuoFac*rgrav(n_r)*&
                     &                     or1(n_r)*ci*real(m,cp)*temp_Rloc(n_m,n_r)
                  end if
               end if
               if ( l_chem ) then
                  dpsi_exp_last(n_m,n_r)=dpsi_exp_last(n_m,n_r)-ChemFac*rgrav(n_r)*&
                  &                      or1(n_r)*ci*real(m,cp)*xi_Rloc(n_m,n_r)
                  if ( vort_bal%l_calc .and. tscheme%istage==1 ) then
                     vort_bal%buo(n_m,n_r)=vort_bal%buo(n_m,n_r)-ChemFac*rgrav(n_r)*&
                     &                     or1(n_r)*ci*real(m,cp)*xi_Rloc(n_m,n_r)
                  end if
               end if
            end if
         end do
      end do
      !$omp end parallel

   end subroutine finish_exp_psi_Rdist
!---------------------------------------------------------------------------------
   subroutine assemble_psi_Rloc(block_sze, nblocks, us_Rloc, up_Rloc, om_Rloc, &
              &                 dpsidt, tscheme, vp_bal, vort_bal)

      !-- Input variables
      integer,             intent(in) :: block_sze
      integer,             intent(in) :: nblocks
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      type(type_tarray),   intent(inout) :: dpsidt
      complex(cp),         intent(inout) :: us_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: up_Rloc(n_m_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: om_Rloc(n_m_max,nRstart:nRstop)
      type(vp_bal_type),   intent(inout) :: vp_bal
      type(vort_bal_type), intent(inout) :: vort_bal

      !-- Local variables
      integer :: n_m_start, n_m_stop, tag, req, nm_block, n_m, nms_block, n_r
      integer :: array_of_requests(4*nblocks)
      complex(cp) :: up0_Rloc(nRstart:nRstop)
      complex(cp) :: work_Rloc(n_m_max, nRstart:nRstop)
      complex(cp) :: work_ghost(n_m_max, nRstart-1:nRstop+1)

      !-- LU factorisation if needed
      if ( .not. lEllmat_FD(1) ) then
         call get_elliptic_mat_Rdist(ellMat_FD)
         lEllmat_FD(:)=.true.
      end if

      !-- First assemble IMEX to get an r.h.s. stored in work_Rloc
      call tscheme%assemble_imex(work_Rloc, dpsidt)

      !-- Handle the uphi0 part
      up0_Rloc(:)=work_Rloc(1,:)

      if ( ktopv /= 1 .and. nRstart==1 ) then
         up0_Rloc(nRstart)=zero 
      end if

      if ( kbotv /= 1 .and. nRstop==n_r_max ) then
         up0_Rloc(nRstop)=zero
      end if
      call bulk_to_ghost(up0_Rloc, up0_ghost, 1, nRstart, nRstop, 1, 1, 1)
      call exch_ghosts(up0_ghost, 1, nRstart, nRstop, 1)

      !--------------------------
      ! Now handle the -\L \psi = \omega equation
      !--------------------------
      array_of_requests(:)=MPI_REQUEST_NULL

      !-- Now solve to finally get psi
      !$omp parallel default(shared) private(tag, req, n_m_start, n_m_stop)
      n_m_start=1; n_m_stop=n_m_max
      call get_openmp_blocks(n_m_start,n_m_stop)

      !-- Non-penetration boundary condition
      if ( nRstart==1 ) then
         do n_m=n_m_start,n_m_stop
            work_Rloc(n_m,1)=zero
         end do
      end if
      if ( nRstop==n_r_max ) then
         do n_m=n_m_start,n_m_stop
            work_Rloc(n_m,n_r_max)=zero
         end do
      end if

      !-- Now copy into an array with proper ghost zones
      call bulk_to_ghost(work_Rloc, work_ghost, 1, nRstart, nRstop, n_m_max, n_m_start, &
           &             n_m_stop)

      tag = 0
      req=1

      do nms_block=1,n_m_max,block_sze
         nm_block = n_m_max-nms_block+1
         if ( nm_block > block_sze ) nm_block=block_sze
         n_m_start=nms_block; n_m_stop=nms_block+nm_block-1
         call get_openmp_blocks(n_m_start,n_m_stop)
         !$omp barrier

         call ellMat_FD%solver_up(work_ghost, n_m_start, n_m_stop, nRstart, nRstop, tag, &
              &                   array_of_requests, req, nms_block, nm_block)
         tag = tag+1
      end do

      do nms_block=1,n_m_max,block_sze
         nm_block = n_m_max-nms_block+1
         if ( nm_block > block_sze ) nm_block=block_sze
         n_m_start=nms_block; n_m_stop=nms_block+nm_block-1
         call get_openmp_blocks(n_m_start,n_m_stop)
         !$omp barrier

         call ellMat_FD%solver_dn(work_ghost, n_m_start, n_m_stop, nRstart, nRstop, tag, &
              &                   array_of_requests, req, nms_block, nm_block)
         tag = tag+1
      end do

      !$omp master
      do nms_block=1,n_m_max,block_sze
         nm_block = n_m_max-nms_block+1
         if ( nm_block > block_sze ) nm_block=block_sze

         call ellMat_FD%solver_finish(work_ghost, nms_block, nm_block, nRstart, nRstop, &
                 &                    tag, array_of_requests, req)
         tag = tag+1
      end do

      call MPI_Waitall(req-1, array_of_requests(1:req-1), MPI_STATUSES_IGNORE, ierr)
      if ( ierr /= MPI_SUCCESS ) call abortRun('MPI_Waitall failed in assemble_pol_Rloc')
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      !$omp end master
      !$omp barrier

      n_m_start=1; n_m_stop=n_m_max
      call get_openmp_blocks(n_m_start,n_m_stop)
      do n_r=nRstart-1,nRstop+1
         do n_m=n_m_start,n_m_stop
            psi_ghost(n_m,n_r)=work_ghost(n_m,n_r)
         end do
      end do
      !$omp end parallel

      ! nRstart-1 and nRstop+1 are already known, only the next one is not known
      !call exch_ghosts(w_ghost, n_m_max, nRstart-1, nRstop+1, 1)
      ! Apparently it yields some problems, not sure why yet
      call exch_ghosts(psi_ghost, n_m_max, nRstart, nRstop, 2)

      !-- Fill ghost zones: this will ensure the second part of the psi bc
      call fill_ghosts_psi(psi_ghost, up0_ghost)

      call get_psi_rhs_imp_ghost(psi_ghost, up0_ghost, us_Rloc, up_Rloc, om_Rloc, &
           &                     dpsidt, 1, vp_bal, vort_bal, tscheme%l_imp_calc_rhs(1))

   end subroutine assemble_psi_Rloc
!---------------------------------------------------------------------------------
   subroutine get_psimat_Rdist(tscheme, psiMat, time_lu, n_lu_calls)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      type(type_penta_par), intent(inout) :: psiMat
      real(cp),             intent(inout) :: time_lu
      integer,              intent(inout) :: n_lu_calls

      !-- Local variables:
      real(cp) :: dm2, dr, fac, runStart, runStop
      integer :: nR, n_m, m

      !-- Bulk points: we fill all the points: this is then easier to handle
      !-- Neumann boundary conditions
      !$omp parallel default(shared) private(nR,m,n_m,dm2)
      !$omp do
      do nR=nRstart,nRstop
         do n_m=1,n_m_max
            m=idx2m(n_m)
            dm2 = real(m*m,cp)
            if ( m == 0 ) cycle
            psiMat%diag(n_m,nR)=-(                     rscheme%ddr(nR,1)  + &
            &                     (or1(nR)+beta(nR))*   rscheme%dr(nR,1)  + &
            &                  (beta(nR)*or1(nR)+dbeta(nR)-dm2*or2(nR)) )   &
            &                      +tscheme%wimp_lin(1)*(ViscFac*         ( &
            &                                        rscheme%ddddr(nR,2)    &
            &          +   (two*or1(nR)+beta(nR))*    rscheme%dddr(nR,2)    &
            &    +( three*dbeta(nR)+two*beta(nR)*or1(nR)-or2(nR)-two*dm2*   &
            &       or2(nR) ) *                        rscheme%ddr(nR,1)    &
            &    +( three*d2beta(nR)+four*dbeta(nR)*or1(nR)-beta(nR)*       &
            &       or2(nR)+or3(nR)+dm2*or2(nR)*(two*or1(nR)-beta(nR))) *   &
            &                                           rscheme%dr(nR,1)    &
            &    +( d3beta(nR)+two*d2beta(nR)*or1(nR)-dbeta(nR)*or2(nR) +   &
            &       beta(nR)*or3(nR)-dm2*or2(nR)*(four*or2(nR)+beta(nR)*    &
            &       or1(nR)+dbeta(nR)-dm2*or2(nR) ) ) ) -                   &
            &                                    CorFac*ekpump(nR)* (       &
            &                                          rscheme%ddr(nR,1)+   &
            &          (beta(nR)+or1(nR)+ekp_up(nR))  * rscheme%dr(nR,1)+   &
            &        (dbeta(nR)+beta(nR)*or1(nR)-dm2*or2(nR)+beta(nR)*      &
            &         ekp_up(nR)-ci*real(m,cp)*or1(nR)*(ekp_us(nR)+ci*      &
            &         real(m,cp)*ekp_dusdp(nR)))) -                         &
            &                   CorFac*beta(nR)*ci*real(m,cp)*or1(nR) )
            psiMat%low1(n_m,nR)=-(                     rscheme%ddr(nR,0)  + &
            &                     (or1(nR)+beta(nR))*   rscheme%dr(nR,0) )  &
            &                      +tscheme%wimp_lin(1)*(ViscFac*         ( &
            &                                        rscheme%ddddr(nR,1)    &
            &          +   (two*or1(nR)+beta(nR))*    rscheme%dddr(nR,1)    &
            &    +( three*dbeta(nR)+two*beta(nR)*or1(nR)-or2(nR)-two*dm2*   &
            &       or2(nR) ) *                        rscheme%ddr(nR,0)    &
            &    +( three*d2beta(nR)+four*dbeta(nR)*or1(nR)-beta(nR)*       &
            &       or2(nR)+or3(nR)+dm2*or2(nR)*(two*or1(nR)-beta(nR))) *   &
            &                                           rscheme%dr(nR,0) )- &
            &                                    CorFac*ekpump(nR)* (       &
            &                                          rscheme%ddr(nR,0)+   &
            &     (beta(nR)+or1(nR)+ekp_up(nR))  *      rscheme%dr(nR,0) ) )
            psiMat%up1(n_m,nR) =-(                     rscheme%ddr(nR,2)  + &
            &                     (or1(nR)+beta(nR))*   rscheme%dr(nR,2) )  &
            &                      +tscheme%wimp_lin(1)*(ViscFac*         ( &
            &                                        rscheme%ddddr(nR,3)    &
            &          +   (two*or1(nR)+beta(nR))*    rscheme%dddr(nR,3)    &
            &    +( three*dbeta(nR)+two*beta(nR)*or1(nR)-or2(nR)-two*dm2*   &
            &       or2(nR) ) *                        rscheme%ddr(nR,2)    &
            &    +( three*d2beta(nR)+four*dbeta(nR)*or1(nR)-beta(nR)*       &
            &       or2(nR)+or3(nR)+dm2*or2(nR)*(two*or1(nR)-beta(nR))) *   &
            &                                           rscheme%dr(nR,2) )- &
            &                                    CorFac*ekpump(nR)* (       &
            &                                          rscheme%ddr(nR,2)+   &
            &     (beta(nR)+or1(nR)+ekp_up(nR))  *      rscheme%dr(nR,2) ) )
            psiMat%low2(n_m,nR)=tscheme%wimp_lin(1)*ViscFac*         (   &
            &                                        rscheme%ddddr(nR,0) &
            &          +   (two*or1(nR)+beta(nR))*    rscheme%dddr(nR,0))
            psiMat%up2(n_m,nR) =tscheme%wimp_lin(1)*ViscFac*         (   &
            &                                        rscheme%ddddr(nR,4) &
            &          +   (two*or1(nR)+beta(nR))*    rscheme%dddr(nR,4))
         end do
      end do
      !$omp end do
      !$omp end parallel

      !-- Boundary conditions
      if ( nRstart == 1) then
         do n_m=1,n_m_max
            m = idx2m(n_m)
            if ( m == 0 ) cycle
            !-- Non-penetration condition at both boundaries
            psiMat%diag(n_m,1)=one
            psiMat%low1(n_m,1)=zero
            psiMat%low2(n_m,1)=zero
            psiMat%up1(n_m,1) =zero
            psiMat%up2(n_m,1) =zero

            !-- Second part of B.C.
            if ( ktopv == 1 ) then ! free slip
               dr=r(2)-r(1)
               fac=(one-half*or1(1)*dr)/(one+half*or1(1)*dr)
               psiMat%diag(n_m,2)=psiMat%diag(n_m,2)-fac*psiMat%low2(n_m,2)
               psiMat%low1(n_m,2)=psiMat%low1(n_m,2)+two/(one+half*or1(1)*dr)*&
               &              psiMat%low2(n_m,2)
            else ! No slip
               psiMat%diag(n_m,2)=psiMat%diag(n_m,2)+psiMat%low2(n_m,2)
            end if
         end do
      end if

      if ( nRstop == n_r_max ) then
         do n_m=1,n_m_max
            m = idx2m(n_m)
            if ( m == 0 ) cycle
            psiMat%diag(n_m,n_r_max)=one
            psiMat%low1(n_m,n_r_max)=zero
            psiMat%low2(n_m,n_r_max)=zero
            psiMat%up1(n_m,n_r_max) =zero
            psiMat%up2(n_m,n_r_max) =zero

            if ( kbotv == 1 ) then ! free-slip
               dr=r(n_r_max)-r(n_r_max-1)
               fac=(one+half*or1(n_r_max)*dr)/(one-half*or1(n_r_max)*dr)
               psiMat%diag(n_m,n_r_max-1)=psiMat%diag(n_m,n_r_max-1) - &
               &                          fac*psiMat%up2(n_m,n_r_max-1)
               psiMat%up1(n_m,n_r_max-1)=psiMat%up1(n_m,n_r_max-1) +   &
               &                         two*psiMat%up2(n_m,n_r_max-1) &
               &                         /(one-half*or1(n_r_max)*dr)
            else ! no slip
               psiMat%diag(n_m,n_r_max-1)=psiMat%diag(n_m,n_r_max-1) + &
               &                          psiMat%up2(n_m,n_r_max-1)
            end if
         end do
      end if

      runStart = MPI_Wtime()
      !-- LU decomposition of the pentadiagonal matrix for \psi
      call psiMat%prepare_mat()
      runStop = MPI_Wtime()
      if ( runStop > runStart ) then
         time_lu = time_lu+(runStop-runStart)
         n_lu_calls = n_lu_calls+1
      end if

   end subroutine get_psimat_Rdist
!---------------------------------------------------------------------------------
   subroutine get_uphimat_Rdist(tscheme, upMat)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      type(type_tri_par), intent(inout) :: upMat

      !-- Local variables:
      integer :: nR, n_m

      n_m=1

      !-- Bulk points: we fill all the points: this is then easier to handle
      !-- Neumann boundary conditions
      do nR=nRstart,nRstop
         upMat%diag(n_m,nR)=one - tscheme%wimp_lin(1)*(ViscFac*  &
         &                (               rscheme%ddr(nR,1) +    &
         &                         or1(nR)*rscheme%dr(nR,1) -    &
         &                                        or2(nR))  -    &
         &                              CorFac*ekpump(nR) )
         upMat%low(n_m,nR)=   -tscheme%wimp_lin(1)*ViscFac*   &
         &                (               rscheme%ddr(nR,0) + &
         &                         or1(nR)*rscheme%dr(nR,0) )
         upMat%up(n_m,nR)=   -tscheme%wimp_lin(1)*ViscFac*    &
         &                (               rscheme%ddr(nR,2) + &
         &                         or1(nR)*rscheme%dr(nR,2) )
      end do

      !-- Boundary conditions
      if ( nRstart== 1 ) then
         if ( ktopv == 1 ) then ! free slip -> duphi/dr - uphi/r = 0
            upMat%up(n_m,1)  =upMat%up(n_m,1)+upMat%low(n_m,1)
            upMat%diag(n_m,1)=upMat%diag(n_m,1)-two*(r(2)-r(1))*or1(1)*upMat%low(n_m,1)
         else ! no slip
            upMat%diag(n_m,1)=1.0_cp
            upMat%low(n_m,1) =0.0_cp
            upMat%up(n_m,1)  =0.0_cp
         end if
      end if

      if ( nRstop == n_r_max ) then
         if ( kbotv == 1 ) then ! free slip
            upMat%low(n_m,n_r_max) =upMat%low(n_m,n_r_max)+upMat%up(n_m,n_r_max)
            upMat%diag(n_m,n_r_max)=upMat%diag(n_m,n_r_max)+two*(r(n_r_max)-r(n_r_max-1))* &
            &                       or1(n_r_max)*upMat%up(n_m,n_r_max)
         else ! no slip
            upMat%diag(n_m,n_r_max)=1.0_cp
            upMat%low(n_m,n_r_max) =0.0_cp
            upMat%up(n_m,n_r_max)  =0.0_cp
         end if
      end if

      !-- LU decomposition:
      call upMat%prepare_mat()

   end subroutine get_uphimat_Rdist
!---------------------------------------------------------------------------------
   subroutine get_elliptic_mat_Rdist(ellMat)
      !
      !  Purpose of this subroutine is to contruct the matrix needed
      !  for the derivation of w for the time advance of the poloidal equation
      !  if the double curl form is used. This is the R-dist version.
      !

      !-- Output variables:
      type(type_tri_par), intent(inout) :: ellMat

      !-- Local variables:
      real(cp) :: dm2
      integer :: nR, m, n_m

      !----- Bulk points:
      do nR=nRstart,nRstop
         do n_m=1,n_m_max
            m = idx2m(n_m)
            if ( m == 0 ) cycle
            dm2 = real(m*m,cp)
            ellMat%diag(n_m,nR)=-(                   rscheme%ddr(nR,1) + &
            &                      (beta(nR)+or1(nR))*rscheme%dr(nR,1) + &
            &                     (dbeta(nR)+beta(nR)*or1(nR)-dm2*or2(nR)) )
            ellMat%up(n_m,nR)  =-(                   rscheme%ddr(nR,2) + &
            &                      (beta(nR)+or1(nR))*rscheme%dr(nR,2) )
            ellMat%low(n_m,nR) =-(                   rscheme%ddr(nR,0) + &
            &                      (beta(nR)+or1(nR))*rscheme%dr(nR,0) )
         end do
      end do

      !-- Non penetrative boundary condition
      if ( nRstart == 1 ) then
         do n_m=1,n_m_max
            ellMat%diag(n_m,1)      =one
            ellMat%up(n_m,1)        =0.0_cp
            ellMat%low(n_m,1)       =0.0_cp
         end do
      end if

      if ( nRstop == n_r_max ) then
         do n_m=1,n_m_max
            ellMat%diag(n_m,n_r_max)=one
            ellMat%up(n_m,n_r_max)  =0.0_cp
            ellMat%low(n_m,n_r_max) =0.0_cp
         end do
      end if

      !-- Lu factorisation
      call ellMat%prepare_mat()

   end subroutine get_elliptic_mat_Rdist
!---------------------------------------------------------------------------------
end module update_psi_fd_mod
