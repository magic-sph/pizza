module init_fields
   !
   ! This module sets the initial fields either by reading in a checkpoint file
   ! or by setting a given perturbation in temperature and/or velocity specified
   ! in the input namelist
   !

   use constants, only: zero, one, two, three, ci, pi, half, sq4pi
   use blocking, only: nRstart, nRstop, nRstart3D, nRstop3D
   use communications, only: transp_r2m, r2m_fields, transp_r2lm, r2lm_fields
   use radial_functions, only: r, rscheme, or1, or2, beta, dbeta, rscheme_3D, &
       &                       r_3D, tcond_3D
   use horizontal, only: theta
   use namelists, only: l_start_file, dtMax, init_t, amp_t, init_u, amp_u, &
       &                radratio, r_cmb, r_icb, l_cheb_coll, l_non_rot,    &
       &                l_reset_t, l_chem, l_heat, amp_xi, init_xi,        &
       &                l_heat_3D
   use outputs, only: n_log_file
   use parallel_mod, only: rank
   use blocking, only: nMstart, nMstop, nM_per_rank, lmStart, lmStop
   use blocking_lm, only: lo_map
   use truncation, only: m_max, n_r_max, minc, m2idx, idx2m, n_phi_max
   use truncation_3D, only: n_r_max_3D, minc_3D, l_max, n_theta_max, &
       &                    n_phi_max_3D, lm_max
   use useful, only: logWrite, abortRun, gausslike_compact_center, &
       &             gausslike_compact_middle, gausslike_compact_edge
   use radial_der, only: get_dr
   use fourier, only: fft
   use checkpoints, only: read_checkpoint
   use time_schemes, only: type_tscheme
   use fields
   use fieldsLast
   use precision_mod
#ifdef WITH_SHTNS
   use shtns, only: spat_to_SH
#endif

   implicit none

   private

   public :: get_start_fields

contains

   subroutine get_start_fields(time, tscheme)

      !-- Output variables
      real(cp),            intent(out) :: time
      class(type_tscheme), intent(inout) :: tscheme

      !-- Local variables
      integer :: m, n_r, n_m
      logical :: lMat
      real(cp) :: h2
      character(len=76) :: message

      if ( l_start_file ) then
         call read_checkpoint(us_Mloc, up_Mloc, temp_Mloc, xi_Mloc, dpsidt, &
              &               dTdt, dxidt, temp_3D_LMloc, dTdt_3D, time, tscheme)

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

         if ( l_heat_3D ) temp_3D_LMloc(:,:)=zero
         if ( l_heat ) temp_Mloc(:,:)=zero
         if ( l_chem ) xi_Mloc(:,:)  =zero
         us_Mloc(:,:) =zero
         up_Mloc(:,:) =zero
         psi_Mloc(:,:)=zero
         call dpsidt%set_initial_values()
         if ( l_heat_3D ) call dTdt_3D%set_initial_values()
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

      if ( init_t /= 0 .and. l_heat_3D ) call initT_3D(temp_3D_LMloc, init_t, amp_t)

      if ( init_u /= 0 ) call initU(us_Mloc, up_Mloc)

      !-- Reconstruct missing fields, dtemp_Mloc, om_Mloc
      if ( l_heat ) call get_dr(temp_Mloc, dtemp_Mloc, nMstart, nMstop, &
                         &      n_r_max, rscheme)
      if ( l_chem ) call get_dr(xi_Mloc, dxi_Mloc, nMstart, nMstop, &
                         &      n_r_max, rscheme)
      if ( l_heat_3D ) call get_dr(temp_3D_LMloc, dtemp_3D_LMloc, lmStart, lmStop,&
                            &      n_r_max_3D, rscheme_3D)
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

   end subroutine get_start_fields
!----------------------------------------------------------------------------------
   subroutine initProp(prop_Mloc, init_prop, amp_prop)

      !-- Input variables
      real(cp),    intent(in) :: amp_prop
      integer,     intent(in) :: init_prop

      !-- Output variable
      complex(cp), intent(inout) :: prop_Mloc(nMstart:nMstop, n_r_max)

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

            !-- temp_Rloc is used as a work r-distributed array here
            call fft(phi_func, temp_Rloc(:,n_r))
         end do

         !-- MPI transpose is needed here
         call transp_r2m(r2m_fields, temp_Rloc, work_Mloc)

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
   subroutine initT_3D(temp_LMloc, init_t, amp_t)

      !-- Input variables
      real(cp),    intent(in) :: amp_t
      integer,     intent(in) :: init_t

      !-- Output variable
      complex(cp), intent(inout) :: temp_LMloc(lmStart:lmStop, n_r_max_3D)

      !-- Local variables
      real(cp) :: tpert(n_r_max_3D), gasp(n_r_max_3D)
      real(cp) :: x, ra1, ra2, c_r, c_i, rdm, rc, c2, c1, sigma_r, LL
      real(cp) :: theta0, phi0, phi
      real(cp) :: ang_func(n_phi_max,n_theta_max)
      integer :: n_r, lm00, lm, l, m, ir, n_th, n_phi
      logical :: rank_has_l0m0

      lm00 = lo_map%lm2(0,0)
      rank_has_l0m0=.false.

      if ( lm00 >= lmStart .and. lm00 <= lmStop ) then
         rank_has_l0m0=.true.
      end if

      if ( (.not. l_start_file) .and. ( rank_has_l0m0) ) then
         !-- Conducting state on the spherically-symmetric part
         do n_r=1,n_r_max_3D
            temp_LMloc(lm00,n_r)=tcond_3D(n_r)*sq4pi
         end do
      end if

      !-- Radial dependence of perturbation in s1:
      do n_r=1,n_r_max_3D
         x=two*r_3D(n_r)-r_cmb-r_icb
         tpert(n_r)=one-three*x**2+three*x**4-x**6
      end do

      if ( init_t < 100 .and. init_t > 0 ) then
         !-- Random noise

         do lm=max(lmStart,2),lmStop
            call random_number(rdm)
            m = lo_map%lm2m(lm)
            l = lo_map%lm2l(lm)
            ra1=(-one+two*rdm)*amp_t/(real(l,cp))**(init_t-1)
            ra2=(-one+two*rdm)*amp_t/(real(l,cp))**(init_t-1)
            do n_r=1,n_r_max_3D
               c_r=ra1*tpert(n_r)
               c_i=ra2*tpert(n_r)
               if ( m > 0 ) then  ! non axisymmetric modes
                  temp_LMloc(lm,n_r)=temp_LMloc(lm,n_r)+cmplx(c_r,c_i,kind=cp)
               else
                  temp_LMloc(lm,n_r)=temp_LMloc(lm,n_r)+cmplx(c_r,0.0_cp,kind=cp)
               end if
            end do
         end do

      else if ( init_t == -1 ) then ! Bubble

         sigma_r = 0.1_cp/sqrt(two)
         rc = half*(r_cmb+r_icb)
         !-- Find the closest point to the middle radius
         ir=minloc(abs(r_3D-rc),1)
         !-- Overwrite rc to be on the grid
         rc = r_3D(ir)
         c1 = half*(r_3D(1)-rc)
         c2 = half*(rc-r_3D(n_r_max_3D))
         LL = half * sigma_r

         !-- Piecewise-definition of a compact-support Gaussian-like profile
         !-- From Eq. (4.7) from Gaspari et al. (1999)
         !-- This ensures that the BCs will always be fulfilled
         do n_r=1,n_r_max_3D
            if ( n_r == ir ) then ! Middle point
               gasp(n_r)=one
            else ! Middle and Edge parts
               if ( r_3D(n_r) < rc-c2 ) then
                  gasp(n_r) = gausslike_compact_edge(rc-r_3D(n_r),c2,LL)
               else if ( r_3D(n_r) >= rc-c2 .and. r_3D(n_r) < rc ) then
                  gasp(n_r) = gausslike_compact_middle(rc-r_3D(n_r),c2,LL)
               else if ( r_3D(n_r) > rc .and. r_3D(n_r) <= rc+c1 ) then
                  gasp(n_r) = gausslike_compact_middle(r_3D(n_r)-rc,c1,LL)
               else if ( r_3D(n_r) > rc+c1 ) then
                  gasp(n_r) = gausslike_compact_edge(r_3D(n_r)-rc,c1,LL)
               end if
            end if
         end do

         !-- Normalisation of the two branches by the middle value
         do n_r=1,n_r_max_3D
            if ( n_r < ir ) then
               gasp(n_r) = gasp(n_r)/gausslike_compact_center(c1,LL)
            else if ( n_r > ir ) then
               gasp(n_r) = gasp(n_r)/gausslike_compact_center(c2,LL)
            end if
         end do

         !-- Now finally define the bubble
         phi0 = pi/minc_3D
         theta0 = half*pi
         do n_r=nRstart3D,nRstop3D
            c_r = amp_t*gasp(n_r)
            do n_th=1,n_theta_max
               do n_phi=1,n_phi_max_3D
                  phi = (n_phi-1)*two*pi/minc_3D/(n_phi_max_3D)
                  ang_func(n_phi,n_th)= c_r*exp(-(phi-phi0)**2/0.2_cp**2)* &
                  &                   exp(-(theta(n_th)-theta0)**2/0.2_cp**2)
               end do
            end do

#ifdef WITH_SHTNS
            call spat_to_SH(ang_func, temp_3D_Rloc(:,n_r))
#endif
         end do

         call transp_r2lm(r2lm_fields, temp_3D_Rloc, work_LMloc)

         do n_r=1,n_r_max_3D
            do lm=lmStart,lmStop
               temp_LMloc(lm,n_r)=temp_LMloc(lm,n_r)+work_LMloc(lm,n_r)
            end do
         end do

      else  if ( init_t >= 100 ) then
         !-- Initialize one specific mode

         l = init_t/100
         if ( l > 99 ) l = init_t/1000
         m = mod(init_t, 100)
         if ( l > 99 ) m = mod(init_t, 1000)
         if ( mod(m,minc_3D) /= 0 ) then
            write(*,*) '! Wave number of mode for 3D temperature initialisation'
            write(*,*) '! not compatible with phi-symmetry:',m
            call abortRun('Stop run in init')
         end if
         if ( l > l_max .or. l < m ) then
            write(*,*) '! Degree of mode for 3D temperature initialisation'
            write(*,*) '! > l_max or < m !',l
            call abortRun('Stop run in init')
         end if

         lm=lo_map%lm2(l,m)
         if( (lm>=lmStart) .and. (lm<=lmStop) ) then
            do n_r=1,n_r_max_3D
               c_r=tpert(n_r)*amp_t
               temp_LMloc(lm,n_r)=temp_LMloc(lm,n_r)+cmplx(c_r,0.0_cp,kind=cp)
            end do

            write(*,'(/'' ! Temperature (3D) initialized at mode:'', &
            &      '' l='',i4,'' m='',i4,'' Ampl='',f8.5)') l,m,amp_t
         end if

      end if

   end subroutine initT_3D
!----------------------------------------------------------------------------------
end module init_fields
