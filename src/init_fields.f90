module init_fields
   !
   ! This module sets the initial fields either by reading in a checkpoint file
   ! or by setting a given perturbation in temperature and/or velocity specified
   ! in the input namelist
   !

   use constants, only: zero, one, two, three, ci, pi, half
   use blocking, only: nRstart, nRstop
   use communications, only: r2m_fields, r2m_single, m2r_single
   use radial_functions, only: r, rscheme, or1, or2, beta, dbeta
   use namelists, only: l_start_file, dtMax, init_t, amp_t, init_u, amp_u, &
       &                radratio, r_cmb, r_icb, l_cheb_coll, l_non_rot,    &
       &                l_reset_t, l_chem, l_heat, amp_xi, init_xi,        &
       &                l_direct_solve, l_finite_diff
   use outputs, only: n_log_file, vp_bal, vort_bal
   use parallel_mod, only: rank
   use blocking, only: nMstart, nMstop, nM_per_rank
   use truncation, only: m_max, n_r_max, minc, m2idx, idx2m, n_phi_max, n_m_max
   use useful, only: logWrite, abortRun, gausslike_compact_center, &
       &             gausslike_compact_middle, gausslike_compact_edge
   use radial_der, only: get_dr, exch_ghosts, bulk_to_ghost, get_dr_Rloc
   use fourier, only: fft
   use checkpoints, only: read_checkpoint
   use time_schemes, only: type_tscheme
   use update_temp_coll, only: get_temp_rhs_imp_coll
   use update_temp_integ, only: get_temp_rhs_imp_int
   use update_temp_fd_mod, only: get_temp_rhs_imp_ghost, temp_ghost, fill_ghosts_temp
   use update_xi_coll, only: get_xi_rhs_imp_coll
   use update_xi_integ, only: get_xi_rhs_imp_int
   use update_xi_fd_mod, only: get_xi_rhs_imp_ghost, xi_ghost, fill_ghosts_xi
   use update_psi_coll_smat, only: get_psi_rhs_imp_coll_smat
   use update_psi_coll_dmat, only: get_psi_rhs_imp_coll_dmat
   use update_psi_integ_smat, only: get_psi_rhs_imp_int_smat
   use update_psi_integ_dmat, only: get_psi_rhs_imp_int_dmat
   use update_psi_fd_mod, only: get_psi_rhs_imp_ghost, psi_ghost, fill_ghosts_psi, &
       &                        up0_ghost
   use fields
   use fieldsLast
   use precision_mod

   implicit none

   private

   public :: get_start_fields

contains

   subroutine get_start_fields(time, tscheme)

      !-- Output variables
      real(cp),            intent(out) :: time
      class(type_tscheme), intent(inout) :: tscheme

      !-- Local variables
      complex(cp) :: psi_Rloc(n_m_max, nRstart:nRstop)
      integer :: m, n_r, n_m
      logical :: lMat
      real(cp) :: h2
      character(len=76) :: message

      if ( l_start_file ) then
         call read_checkpoint(us_Mloc, up_Mloc, temp_Mloc, xi_Mloc, dpsidt, &
              &               dTdt, dxidt, time, tscheme)

         if ( l_reset_t ) time = 0.0_cp

         !-- If integration method is used, since u_s is stored, one needs to
         !-- reconstruct psi(m/=0)
         if ( .not. l_cheb_coll ) then
            do n_r=2,n_r_max
               if ( l_non_rot ) then
                  h2 = one
               else
                  h2 = r_cmb*r_cmb-r(n_r)*r(n_r)
               end if
               do n_m=nMstart,nMstop
                  m = idx2m(n_m)
                  if ( m > 0 ) then
                     psi_Mloc(n_m, n_r) = -ci*r(n_r)/real(m,cp)/h2 * us_Mloc(n_m, n_r)
                  else
                     psi_Mloc(n_m, n_r) = 0.0_cp
                  end if
               end do
            end do
            !-- Boundary point (since h2 is singular there)
            do n_m=nMstart, nMstop
               psi_Mloc(n_m,1)=zero
            end do
         end if

      else

         if ( l_heat ) temp_Mloc(:,:)=zero
         if ( l_chem ) xi_Mloc(:,:)  =zero
         us_Mloc(:,:)  =zero
         up_Mloc(:,:)  =zero
         psi_Mloc(:,:) =zero
         call dpsidt%set_initial_values()
         if ( l_heat ) call dTdt%set_initial_values()
         if ( l_chem ) call dxidt%set_initial_values()

         time=0.0_cp
         tscheme%dt(:)=dtMax

         if (rank == 0) write(message,'(''! Using dtMax time step:'',ES16.6)') dtMax
         call logWrite(message, n_log_file)

      end if

      !-- Initialize the weights of the time scheme
      call tscheme%set_weights(lMat)

      if ( init_t /= 0 .and. l_heat ) call initProp(temp_Mloc, init_t, amp_t)

      if ( init_xi /= 0 .and. l_chem ) call initProp(xi_Mloc, init_xi, amp_xi)

      if ( init_u /= 0 ) call initU(us_Mloc, up_Mloc)

      if ( l_finite_diff ) then ! Finite diff. are used in radius

         if ( l_heat ) call m2r_single%transp_m2r(temp_Mloc, temp_Rloc)
         if ( l_chem ) call m2r_single%transp_m2r(xi_Mloc, xi_Rloc)
         call m2r_single%transp_m2r(us_Mloc, us_Rloc)
         call m2r_single%transp_m2r(up_Mloc, up_Rloc)

         if ( l_heat ) then
            call bulk_to_ghost(temp_RLoc, temp_ghost, 1, nRstart, nRstop, n_m_max, &
                 &             1, n_m_max)
            call exch_ghosts(temp_ghost, n_m_max, nRstart, nRstop, 1)
            call fill_ghosts_temp(temp_ghost)
            call get_temp_rhs_imp_ghost(temp_ghost, dtemp_Rloc, dTdt, 1, .true.)
         end if

         if ( l_chem ) then
            call bulk_to_ghost(xi_RLoc, xi_ghost, 1, nRstart, nRstop, n_m_max, &
                 &             1, n_m_max)
            call exch_ghosts(xi_ghost, n_m_max, nRstart, nRstop, 1)
            call fill_ghosts_xi(xi_ghost)
            call get_xi_rhs_imp_ghost(xi_ghost, dxi_Rloc, dxidt, 1, .true.)
         end if

         !-- psi is unknown: rebuild it from us
         do n_r=nRstart,nRstop
            do n_m=2,n_m_max
               m = idx2m(n_m)
               psi_Rloc(n_m,n_r)=us_Rloc(n_m,n_r)/(ci*real(m,cp)) * r(n_r)
            end do
         end do
         call bulk_to_ghost(up_Rloc(1,:), up0_ghost, 1, nRstart, nRstop, 1, 1, 1)
         call exch_ghosts(up0_ghost, 1, nRstart, nRstop, 1)
         call bulk_to_ghost(psi_RLoc, psi_ghost, 2, nRstart, nRstop, n_m_max, &
              &             1, n_m_max)
         call exch_ghosts(psi_ghost, n_m_max, nRstart, nRstop, 2)
         call fill_ghosts_psi(psi_ghost, up0_ghost)
         call get_psi_rhs_imp_ghost(psi_ghost, up0_ghost, us_Rloc, up_Rloc, om_Rloc, &
              &                     dpsidt, 1, vp_bal, vort_bal, .true.)

         !-- Construct vorticity (use psi_Rloc for temp array)
         call get_dr_Rloc(up_Rloc, psi_Rloc, n_m_max, nRstart, nRstop, n_r_max, rscheme)
         do n_r=nRstart,nRstop
            do n_m=1,n_m_max
               m = idx2m(n_m)
               om_Rloc(n_m,n_r)=psi_Rloc(n_m,n_r)+or1(n_r)*up_Rloc(n_m,n_r)- &
               &                    or1(n_r)*ci*real(m,cp)*us_Rloc(n_m,n_r)
            end do
         end do

      else ! Chebyshev polynomials are used for the radial representation

         !-- Reconstruct missing fields, dtemp_Mloc, om_Mloc
         if ( l_heat ) call get_dr(temp_Mloc, dtemp_Mloc, nMstart, nMstop, &
                            &      n_r_max, rscheme)
         if ( l_chem ) call get_dr(xi_Mloc, dxi_Mloc, nMstart, nMstop, &
                            &      n_r_max, rscheme)
         call get_dr(up_Mloc, work_Mloc, nMstart, nMstop, n_r_max, rscheme)
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               om_Mloc(n_m,n_r)=work_Mloc(n_m,n_r)+or1(n_r)*up_Mloc(n_m,n_r)- &
               &                     or1(n_r)*ci*real(m,cp)*us_Mloc(n_m,n_r)
            end do
         end do

         !-- When not using Collocation also store temp_hat and psi_hat
         !-- This saves 2 DCTs per iteration
         if ( .not. l_cheb_coll ) then
            do n_r=1,n_r_max
               do n_m=nMstart,nMstop
                  if ( l_heat ) temp_hat_Mloc(n_m,n_r)=temp_Mloc(n_m,n_r)
                  if ( l_chem ) xi_hat_Mloc(n_m,n_r)=xi_Mloc(n_m,n_r)
                  psi_hat_Mloc(n_m,n_r) =psi_Mloc(n_m,n_r)
               end do
            end do
            if ( l_heat ) call rscheme%costf1(temp_hat_Mloc, nMstart, nMstop, &
                               &              n_r_max)
            if ( l_chem ) call rscheme%costf1(xi_hat_Mloc, nMstart, nMstop, &
                               &              n_r_max)
            call rscheme%costf1(psi_hat_Mloc, nMstart, nMstop, n_r_max)
         end if

         !-- Compute the first implicit state
         if ( l_cheb_coll ) then
            if ( l_heat ) call get_temp_rhs_imp_coll(temp_Mloc, dtemp_Mloc, dTdt, 1, .true.)
            if ( l_chem ) call get_xi_rhs_imp_coll(xi_Mloc, dxi_Mloc, dxidt, 1, .true.)

            if ( l_direct_solve ) then
               call get_psi_rhs_imp_coll_smat(us_Mloc, up_Mloc, om_Mloc, dom_Mloc,   &
                    &                         temp_Mloc, xi_Mloc, dpsidt, 1, vp_bal, &
                    &                         vort_bal, .true.)
            else
               call get_psi_rhs_imp_coll_dmat(up_Mloc, om_Mloc, dom_Mloc, temp_Mloc, &
                    &                         xi_Mloc, dpsidt, 1, vp_bal, vort_bal, .true.)
            end if
         else
            if ( l_heat ) call get_temp_rhs_imp_int(temp_hat_Mloc, dTdt, 1, .true.)
            if ( l_chem ) call get_xi_rhs_imp_int(xi_hat_Mloc, dxidt, 1, .true.)
            if ( l_direct_solve ) then
               call get_psi_rhs_imp_int_smat(psi_hat_Mloc,up_Mloc,temp_Mloc,psi_Mloc, &
                    &                        dpsidt, 1, vp_bal, .true.)
            else
               call get_psi_rhs_imp_int_dmat(om_Mloc,up_Mloc,temp_Mloc,xi_Mloc,dpsidt,&
                    &                        1, vp_bal, .true.)
            end if
         end if
      end if

   end subroutine get_start_fields
!----------------------------------------------------------------------------------
   subroutine initProp(prop_Mloc, init_prop, amp_prop)

      !-- Output variables
      complex(cp), intent(inout) :: prop_Mloc(nMstart:nMstop, n_r_max)
      real(cp),    intent(in) :: amp_prop
      integer,     intent(in) :: init_prop

      !-- Local variables
      integer :: m_pertu, n_r, idx, n_phi, n_m, ir
      real(cp) :: x, c_r, rdm, c1, c2, rc, L, sigma_r
      real(cp) :: t1(n_r_max), gasp(n_r_max)
      real(cp) :: phi, phi0
      real(cp) :: phi_func(n_phi_max)

      !-- Radial dependence of perturbation in t1:
      do n_r=1,n_r_max
         x=two*r(n_r)-r_cmb-r_icb
         t1(n_r)=sin(pi*(r(n_r)-r_icb))
      end do

      if ( init_prop > 0 ) then ! Initialize a peculiar m mode
         
         m_pertu = init_prop

         if ( mod(m_pertu,minc) /= 0 ) then
            write(*,*) '! Wave number of mode for temperature initialisation'
            write(*,*) '! not compatible with phi-symmetry:',m_pertu
            call abortRun('Stop run in init')
         end if
         if ( m_pertu > m_max ) then
            write(*,*) '! Degree of mode for temperature initialisation'
            write(*,*) '! > m_max  !',m_pertu
            call abortRun('Stop run in init')
         end if

         idx = m2idx(m_pertu)
         if ( idx >= nMstart .and. idx <= nMstop ) then
            do n_r=1,n_r_max
               c_r=t1(n_r)*amp_prop
               prop_Mloc(idx,n_r)=prop_Mloc(idx,n_r)+cmplx(c_r,0.0_cp,kind=cp)
            end do
         end if

      else if ( init_prop == -1 ) then ! bubble = Gaussian in r and phi

         sigma_r = 0.1_cp/sqrt(two)
         rc = half*(r_cmb+r_icb)
         !-- Find the closest point to the middle radius
         ir=minloc(abs(r-rc),1)
         !-- Overwrite rc to be on the grid
         rc = r(ir)
         c1 = half*(r(1)-rc)
         c2 = half*(rc-r(n_r_max))
         L = half * sigma_r

         !-- Piecewise-definition of a compact-support Gaussian-like profile
         !-- From Eq. (4.7) from Gaspari et al. (1999)
         !-- This ensures that the BCs will always be fulfilled
         do n_r=1,n_r_max
            if ( n_r == ir ) then ! Middle point
               gasp(n_r)=one
            else ! Middle and Edge parts
               if ( r(n_r) < rc-c2 ) then
                  gasp(n_r) = gausslike_compact_edge(rc-r(n_r),c2,L)
               else if ( r(n_r) >= rc-c2 .and. r(n_r) < rc ) then
                  gasp(n_r) = gausslike_compact_middle(rc-r(n_r),c2,L)
               else if ( r(n_r) > rc .and. r(n_r) <= rc+c1 ) then
                  gasp(n_r) = gausslike_compact_middle(r(n_r)-rc,c1,L)
               else if ( r(n_r) > rc+c1 ) then
                  gasp(n_r) = gausslike_compact_edge(r(n_r)-rc,c1,L)
               end if
            end if
         end do

         !-- Normalisation of the two branches by the middle value
         do n_r=1,n_r_max
            if ( n_r < ir ) then
               gasp(n_r) = gasp(n_r)/gausslike_compact_center(c1,L)
            else if ( n_r > ir ) then
               gasp(n_r) = gasp(n_r)/gausslike_compact_center(c2,L)
            end if
         end do

         !-- Now finally define the bubble
         phi0 = pi/minc
         do n_r=nRstart,nRstop
            c_r = amp_prop*gasp(n_r)
            do n_phi=1,n_phi_max
               phi = (n_phi-1)*two*pi/minc/(n_phi_max)
               phi_func(n_phi)=c_r*exp(-(phi-phi0)**2/0.2_cp**2)
            end do

            !-- us_Rloc is used as a work r-distributed array here
            call fft(phi_func, us_Rloc(:,n_r))
         end do

         !-- MPI transpose is needed here
         call r2m_single%transp_r2m(us_Rloc, work_Mloc)

         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               prop_Mloc(n_m,n_r) = prop_Mloc(n_m,n_r)+work_Mloc(n_m,n_r)
            end do
         end do

      else ! random noise

         do n_r=1,n_r_max
            t1(n_r)=sin(pi*(r(n_r)-r_icb))
         end do

         do n_r=1,n_r_max
            do n_m=nMstart, nMstop
               m_pertu = idx2m(n_m)
               if ( m_pertu > 0 ) then
                  call random_number(rdm)
                  prop_Mloc(n_m, n_r) = amp_prop*rdm*m_pertu**(-1.5_cp)*t1(n_r)
               end if
            end do
         end do

      end if

   end subroutine initProp
!----------------------------------------------------------------------------------
   subroutine initU(us_Mloc, up_Mloc)

      !-- Output variables
      complex(cp), intent(inout) :: us_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(inout) :: up_Mloc(nMstart:nMstop, n_r_max)

      !-- Local variables
      integer :: m_pertu, n_r, idx, m, n_m
      real(cp) :: c_r
      real(cp) :: u1(n_r_max)


      !-- Radial dependence of perturbation in t1:
      do n_r=1,n_r_max
         u1(n_r)=sin(pi*(r(n_r)-r_icb))
      end do

      if ( init_u > 0 ) then
         
         m_pertu = init_u

         if ( mod(m_pertu,minc) /= 0 ) then
            write(*,*) '! Wave number of mode for velocity initialisation'
            write(*,*) '! not compatible with phi-symmetry:',m_pertu
            call abortRun('Stop run in init')
         end if
         if ( m_pertu > m_max ) then
            write(*,*) '! Degree of mode for velocity initialisation'
            write(*,*) '! > m_max  !',m_pertu
            call abortRun('Stop run in init')
         end if

         idx = m2idx(m_pertu)
         if ( idx >= nMstart .and. idx <= nMstop ) then
            do n_r=1,n_r_max
               c_r=u1(n_r)*amp_u
               us_Mloc(idx,n_r)=us_Mloc(idx,n_r)+cmplx(c_r,0.0_cp,kind=cp)
            end do
         end if

         !-- Get the corresponding vorticity
         do n_r=1,n_r_max
            do n_m=nMstart, nMstop
               m = idx2m(n_m)
               om_Mloc(n_m,n_r)=-ci*m*or1(n_r)*us_Mloc(n_m,n_r)
            end do
         end do

      else ! initialize an axisymmetric vphi

         idx = m2idx(0)
         if ( idx >= nMstart .and. idx <= nMstop ) then
            do n_r=1,n_r_max
               c_r=amp_u*sin(pi*(r(n_r)-r_icb))
               up_Mloc(idx,n_r)=up_Mloc(idx,n_r)+cmplx(c_r,0.0_cp,kind=cp)
            end do
         end if

      end if

   end subroutine initU
!----------------------------------------------------------------------------------
end module init_fields
