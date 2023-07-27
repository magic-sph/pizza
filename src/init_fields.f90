module init_fields
   !
   ! This module sets the initial fields either by reading in a checkpoint file
   ! or by setting a given perturbation in temperature and/or velocity specified
   ! in the input namelist
   !

   use iso_fortran_env, only: output_unit
   use constants, only: zero, one, two, three, four, ci, pi, half, third, sq4pi
   use blocking, only: nRstart, nRstop, nRstart3D, nRstop3D
   use communications, only: transp_r2m, r2m_fields, transp_r2lm, r2lm_fields, &
       &               transp_lm2r, lm2r_fields, scatter_from_rank0_to_lmloc
   use radial_functions, only: r, rscheme, or1, or2, beta, rscheme_3D, &
       &                       r_3D, or1_3D, tcond_3D
   use horizontal, only: theta
   use namelists, only: l_start_file, dtMax, init_t, amp_t, init_u, amp_u,     &
       &                radratio, r_cmb, r_icb, l_cheb_coll, l_QG_basis,       &
       &                l_non_rot, l_reset_t, l_chem, l_heat, amp_xi, init_xi, &
       &                l_heat_3D, l_mag_3D, amp_B, init_B, BdiffFac, l_mag_B0,&
       &                amp_B0, init_B0, lmax_trunc_B0, start_file_b0,         &
       &                l_mag_inertia
   use outputs, only: n_log_file
   use parallel_mod, only: rank
   use algebra, only: prepare_full_mat, solve_full_mat
   use blocking, only: nMstart, nMstop, nM_per_rank, lmStart, lmStop
   use blocking_lm, only: lo_map
   use truncation, only: m_max, n_r_max, minc, m2idx, idx2m, n_phi_max
   use truncation_3D, only: n_r_max_3D, minc_3D, l_max, n_theta_max, &
       &                    n_phi_max_3D, m_max_3D, lm_max
   use useful, only: logWrite, abortRun, gausslike_compact_center, &
       &             gausslike_compact_middle, gausslike_compact_edge
   use radial_der, only: get_dr, get_ddr
   use fourier, only: fft
   use checkpoints, only: read_checkpoint
   use time_schemes, only: type_tscheme
   use fields
   use fieldsLast
   use precision_mod
#ifdef WITH_SHTNS
   use shtns, only: scal_to_SH, torpol_to_spat, torpol_to_curl_spat
#endif

   implicit none

   private

   public :: get_start_fields

   !----- Peak values for magnetic field:
   real(cp), public :: bpeakbot, bpeaktop

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
         call read_checkpoint(us_Mloc, up_Mloc, temp_Mloc, xi_Mloc, dpsidt,   &
              &               dTdt, dxidt, temp_3D_LMloc, dTdt_3D, b_3D_LMloc,&
              &               dBdt_3D, aj_3D_LMloc, djdt_3D, time, tscheme)

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

         if ( l_mag_3D ) then
            b_3D_LMloc(:,:) =zero
            aj_3D_LMloc(:,:)=zero
            db_3D_LMloc(:,:)=zero
            dj_3D_LMloc(:,:)=zero
            ddb_3D_LMloc(:,:)=zero
         end if
         if ( l_heat_3D ) temp_3D_LMloc(:,:)=zero
         if ( l_heat ) temp_Mloc(:,:)=zero
         if ( l_chem ) xi_Mloc(:,:)  =zero
         us_Mloc(:,:) =zero
         up_Mloc(:,:) =zero
         psi_Mloc(:,:)=zero
         call dpsidt%set_initial_values()
         if ( l_mag_3D ) then
            call dBdt_3D%set_initial_values()
            call djdt_3D%set_initial_values()
         end if
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

      if ( l_mag_3D ) then
         if ( init_B /= 0 ) call initB_3D(b_3D_LMloc, aj_3D_LMloc, init_B, amp_B)
         if ( init_B0 /= 0 .and. l_mag_B0 ) call initB0_3D(b_3D_LMloc, aj_3D_LMloc, &
                                                 &         init_B0, amp_B0, lmax_trunc_B0)
      end if

      !-- Reconstruct missing fields, dtemp_Mloc, om_Mloc, dB, dj, etc.
      if ( l_heat ) call get_dr(temp_Mloc, dtemp_Mloc, nMstart, nMstop, &
                         &      n_r_max, rscheme)
      if ( l_chem ) call get_dr(xi_Mloc, dxi_Mloc, nMstart, nMstop, &
                         &      n_r_max, rscheme)
      if ( l_heat_3D ) call get_dr(temp_3D_LMloc, dtemp_3D_LMloc, lmStart, lmStop,&
                            &      n_r_max_3D, rscheme_3D)
      if ( l_mag_3D ) then
         call get_ddr(b_3D_LMloc, db_3D_LMloc, ddb_3D_LMloc, lmStart, lmStop,&
              &      n_r_max_3D, rscheme_3D)
         call get_dr(aj_3D_LMloc, dj_3D_LMloc, lmStart, lmStop,&
              &      n_r_max_3D, rscheme_3D)
      end if
      call get_dr(up_Mloc, work_Mloc, nMstart, nMstop, n_r_max, rscheme)
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            om_Mloc(n_m,n_r)=work_Mloc(n_m,n_r)+or1(n_r)*up_Mloc(n_m,n_r)- &
            &                     or1(n_r)*ci*real(m,cp)*us_Mloc(n_m,n_r)
            if ( l_QG_basis .and. m>0 ) then ! No phi derivative at m=0
               om_Mloc(n_m,n_r)=om_Mloc(n_m,n_r) + beta(n_r)*third*        &
               &                           ci*real(m,cp)*us_Mloc(n_m,n_r)
            end if
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
            write(output_unit,*) &
            &            '! Wave number of mode for temperature initialisation'
            write(output_unit,*) '! not compatible with phi-symmetry:',m_pertu
            call abortRun('Stop run in init')
         end if
         if ( m_pertu > m_max ) then
            write(output_unit,*) '! Degree of mode for temperature initialisation'
            write(output_unit,*) '! > m_max  !',m_pertu
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
            write(output_unit,*) '! Wave number of mode for velocity initialisation'
            write(output_unit,*) '! not compatible with phi-symmetry:',m_pertu
            call abortRun('Stop run in init')
         end if
         if ( m_pertu > m_max ) then
            write(output_unit,*) '! Degree of mode for velocity initialisation'
            write(output_unit,*) '! > m_max  !',m_pertu
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
            call scal_to_SH(ang_func, temp_3D_Rloc(:,n_r))
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
            write(output_unit,*)  &
            &      '! Wave number of mode for 3D temperature initialisation'
            write(output_unit,*) '! not compatible with phi-symmetry:',m
            call abortRun('Stop run in init')
         end if
         if ( l > l_max .or. l < m ) then
            write(output_unit,*) &
            &      '! Degree of mode for 3D temperature initialisation'
            write(output_unit,*) '! > l_max or < m !',l
            call abortRun('Stop run in init')
         end if

         lm=lo_map%lm2(l,m)
         if( (lm>=lmStart) .and. (lm<=lmStop) ) then
            do n_r=1,n_r_max_3D
               c_r=tpert(n_r)*amp_t
               temp_LMloc(lm,n_r)=temp_LMloc(lm,n_r)+cmplx(c_r,0.0_cp,kind=cp)
            end do

            write(output_unit,'(/'' ! Temperature (3D) initialized at mode:'', &
            &      '' l='',i4,'' m='',i4,'' Ampl='',f8.5)') l,m,amp_t
         end if

      end if

   end subroutine initT_3D
!----------------------------------------------------------------------------------
   subroutine initB_3D(b_LMloc, aj_LMloc, init_B, amp_B)
      !
      ! Purpose of this subroutine is to initialize the magnetic field
      ! according to the control parameters amp_B and init_B.
      !

      !-- Input variables
      real(cp),    intent(in) :: amp_B
      integer,     intent(in) :: init_B

      !-- Output variables:
      complex(cp), intent(inout) :: b_LMloc(lmStart:lmStop, n_r_max_3D)
      complex(cp), intent(inout) :: aj_LMloc(lmStart:lmStop, n_r_max_3D)

      !-- Local variables:
      integer :: n_r
      real(cp) :: b_pol, b_tor
      !-- poloidal and toroidal arrays for grid_space quantities computation
      complex(cp) :: db_LMloc(lmStart:lmStop, n_r_max_3D)
      complex(cp) :: b_Rloc(lm_max,nRstart3D:nRstop3D)
      complex(cp) :: db_Rloc(lm_max,nRstart3D:nRstop3D)
      complex(cp) :: aj_Rloc(lm_max,nRstart3D:nRstop3D)
      real(cp) :: Br_shell(n_phi_max_3D,n_theta_max)
      real(cp) :: Bt_shell(n_phi_max_3D,n_theta_max)
      real(cp) :: Bp_shell(n_phi_max_3D,n_theta_max)
      !-- arrays and helpers to read complex b(t=0) from any code
      integer :: llmm2lm(l_max+1,m_max_3D+1)
      real(cp) :: load_head(8)
      real(cp), allocatable :: work(:,:)
      complex(cp) :: load_b0pol(lm_max, n_r_max_3D)
      complex(cp) :: load_b0tor(lm_max, n_r_max_3D)
      integer :: l_max_old, l_max_f, lm_max_old, lm_max_f
      integer :: n_start_file, n_r_max_old, n_r_max_f
      logical :: startfile_does_exist

      real(cp) :: b1(n_r_max_3D)
      real(cp) :: bR, bI, rdm
      real(cp) :: x
      integer :: bExp

      integer :: lm, l1, m1, k
      integer :: lm10, lm20, lm30, lm11, lm44

      lm10 = lo_map%lm2(1,0)
      lm20 = lo_map%lm2(2,0)
      lm30 = lo_map%lm2(3,0)
      lm11 = lo_map%lm2(1,1)
      lm44 = lo_map%lm2(4,4)


      if ( init_B < 0 ) then  ! l,m mixture, random init
      !-- Random noise initialization of all (l,m) modes exept (l=0,m=0):
         bExp=abs(init_B)
         do n_r=1,n_r_max_3D
            !x=two*r_3D(n_r)-r_cmb-r_icb
            !b1(n_r)=one-three*x**2+three*x**4-x**6
            b1(n_r)=(r_3D(n_r)/r_cmb)**2 * ( one-three/5.0_cp*(r_3D(n_r)/r_cmb)**2 )
         end do

         do lm=max(lmStart,2),lmStop
            call random_number(rdm)
            l1=lo_map%lm2l(lm)
            m1=lo_map%lm2m(lm)
            if ( l1 > 0 ) then
               bR=(-one+two*rdm)*amp_B/(real(l1,cp))**(bExp-1)
               bI=(-one+two*rdm)*amp_B/(real(l1,cp))**(bExp-1)
            else
               bR=0.0_cp
               bI=0.0_cp
            end if
            if ( m1 == 0 ) bI=0.0_cp
            do n_r=1,n_r_max_3D
               b_LMloc(lm,n_r)=b_LMloc(lm,n_r) + cmplx(bR*b1(n_r),bI*b1(n_r),kind=cp)
            end do
         end do

      else if ( init_B == 1 ) then  ! l=1,m0 poloidal field
      ! which is potential field at r_cmb
         if( (lm10>=lmStart) .and. (lm10<=lmStop) ) then ! select processor
!         if ( lmStartB(rank+1) <= lm10 .and. lmStopB(rank+1) >= lm10 ) then ! select processor
            b_pol=amp_B*5.0_cp*half*sqrt(third*pi)*r_icb**2
            do n_r=1,n_r_max_3D
               b_LMloc(lm10,n_r)=b_LMloc(lm10,n_r)+b_pol*(r_3D(n_r)/r_icb)**2 * &
               &               ( one - three/5.0_cp*(r_3D(n_r)/r_cmb)**2 )
            end do
         end if

      else if ( init_B == 2 ) then  ! l=1,m=0 analytical toroidal field
      ! with a maximum of amp_B at mid-radius
      ! between r_icb and r_cmb for an insulating
      ! inner core and at r_cmb/2 for a conducting
      ! inner core
         if( (lm10>=lmStart) .and. (lm10<=lmStop) ) then ! select processor
!         if ( lmStartB(rank+1) <= lm10 .and. lmStopB(rank+1) >= lm10 ) then ! select processor
            b_tor=-two*amp_B*sqrt(third*pi)  ! minus sign makes phi comp. > 0
            do n_r=1,n_r_max_3D
               aj_LMloc(lm10,n_r)=aj_LMloc(lm10,n_r) + b_tor*r_3D(n_r)*sin(pi*(r_3D(n_r)-r_icb))
            end do
         end if

      else if ( init_B == 3 ) then
         ! l=2,m=0 toroidal field and l=1,m=0 poloidal field
         ! toroidal field has again its maximum of amp_B
         ! at mid-radius between r_icb and r_cmb for an
         ! insulating inner core and at r_cmb/2 for a
         ! conducting inner core
         ! The outer core poloidal field is defined by
         ! a homogeneous  current density, its maximum at
         ! the ICB is set to amp_B.
         ! The inner core poloidal field is chosen accordingly.
         if( (lm10>=lmStart) .and. (lm10<=lmStop) ) then ! select processor
!         if ( lmStartB(rank+1) <= lm10 .and. lmStopB(rank+1) >= lm10 ) then ! select processor
            b_tor=-four*third*amp_B*sqrt(pi/5.0_cp)
            b_pol= amp_B*sqrt(three*pi)/four
            do n_r=1,n_r_max_3D
               b_LMloc(lm10,n_r)=b_LMloc(lm10,n_r)+b_pol * (       &
               &     r_3D(n_r)**3 - four*third*r_cmb*r_3D(n_r)**2  &
               &     + third*r_icb**4/r_3D(n_r)                )
            end do

            !b_tor=-four*third*amp_B*sqrt(pi/20.0_cp)
            !b_pol= amp_B*sqrt(three*pi)/10.0_cp
            !do n_r=1,n_r_max_3D
            !   b_LMloc(lm44,n_r)=b_LMloc(lm44,n_r)+b_pol * (       &
            !   &     (r_3D(n_r)-r_icb)**3 - four*third*r_cmb*r_3D(n_r)**2  &
            !   &     + third*r_icb**4/r_3D(n_r)                )
            !end do
         end if

         if( (lm20>=lmStart) .and. (lm20<=lmStop) ) then ! select processor
!         if ( lmStartB(rank+1) <= lm20 .and. lmStopB(rank+1) >= lm20 ) then ! select processor
            b_tor=-four*third*amp_B*sqrt(pi/5.0_cp)
            b_pol=amp_B*sqrt(three*pi)/four
            do n_r=1,n_r_max_3D
               aj_LMloc(lm20,n_r)=aj_LMloc(lm20,n_r) + b_tor*r_3D(n_r)*sin(pi*(r_3D(n_r)-r_icb))
            end do
         end if

      else if ( init_B == 4 ) then  ! l=1,m0 poloidal field
      ! with max field amplitude amp_B at r_icb
         if( (lm10>=lmStart) .and. (lm10<=lmStop) ) then ! select processor
!         if ( lmStartB(rank+1) <= lm10 .and. lmStopB(rank+1) >= lm10 ) then ! select processor
            b_pol=-amp_B*r_icb**3*sqrt(third*pi)
            do n_r=1,n_r_max_3D
               b_LMloc(lm10,n_r)=b_LMloc(lm10,n_r)+b_pol*or1_3D(n_r)
            end do
         end if

     else if ( init_B == 5 ) then  ! l=1,m=0 poloidal field , constant in r !
      ! no potential at r_cmb but simple
         if( (lm10>=lmStart) .and. (lm10<=lmStop) ) then ! select processor
!         if ( lmStartB(rank+1) <= lm10 .and. lmStopB(rank+1) >= lm10 ) then ! select processor
            b_pol=amp_B
            do n_r=1,n_r_max_3D
               b_LMloc(lm10,n_r)=b_LMloc(lm10,n_r)+b_pol*r_3D(n_r)**2!*(r_cmb+r_icb)/2**2!
            end do
         end if

      else if ( init_B == 6 ) then  ! l=1,m0 poloidal field
      ! which is potential field at r_cmb
      ! SAME as init_B == 1: kept this one because for easier comparisons with MagIC 
         if( (lm10>=lmStart) .and. (lm10<=lmStop) ) then ! select processor
!         if ( lmStartB(rank+1) <= lm10 .and. lmStopB(rank+1) >= lm10 ) then ! select processor
            b_pol=amp_B*5.0_cp*half*sqrt(third*pi)*r_icb**2
            do n_r=1,n_r_max_3D
               b_LMloc(lm10,n_r)=b_LMloc(lm10,n_r)+b_pol*(r_3D(n_r)/r_icb)**2 * &
               &               ( one - three/5.0_cp*(r_3D(n_r)/r_cmb)**2 )
            end do
         end if

      else if ( init_B == -1 ) then  ! Read a file containing a starting field of any complexity from any code
      ! b0 format:: Spherical harmonics and chebychev polynomials
      !             (2*lm_max,2*n_r_max)=(complex stored in real array, bpol+btor)
      ! !WARNING!:: Fortran format =' formatted' ('unformatted' is for binary files; not implemented yet)
         if ( rank == 0 ) then
            inquire(file=start_file_b0, exist=startfile_does_exist)

            if ( startfile_does_exist ) then
               open(newunit=n_start_file, file=start_file_b0, status='old', &
               &    form='formatted', access='stream')
               print*, "! Init_B Starting from a B field stored in ", start_file_b0
            else
               call abortRun('! The restart file for B does not exist !')
            end if

            read(n_start_file, *) load_head ! unknown length
            l_max_old = int(load_head(1)) ! l_max from checkpoint
            lm_max_old = int(load_head(2)) ! lm_max from checkpoint
            n_r_max_old = int(load_head(4))! n_r_max from checkpoint
            !-- extracting bpol/btor from the file
            allocate(work(2*lm_max_old, 2*n_r_max_old))
            do lm=1,2*lm_max_old
               read(n_start_file, *) work(lm,:)
            end do
            close(n_start_file)

            !-- lookup table for converting lm to l,m from checkpoint file
            lm=0
            do m1=0,m_max_3D
               do l1=m1,l_max
                  lm=lm+1
                  llmm2lm(l1+1,m1+1)=lm
               end do
            end do

            !-- adapting to the current resolution
            l_max_f = min(l_max, l_max_old)
            lm_max_f = min(lm_max, lm_max_old)
            n_r_max_f = min(n_r_max_3D, n_r_max_old)
            !-- De-sorting the coefficient for the Back transforms
            k=1
            do l1=0,l_max_f!-- WARNING: must start at 0 because (l,m)=(0,0) is stored as well
               do m1=0,l1
                  lm = llmm2lm(l1+1,m1+1)
                  !lm = lo_map%lm2(l,m)
                  load_b0pol(lm,:n_r_max_f) = amp_B*cmplx(work(k  ,             1:            n_r_max_f), &
                  &                                       work(k+1,             1:            n_r_max_f), kind=cp)
                  load_b0tor(lm,:n_r_max_f) = amp_B*cmplx(work(k  , n_r_max_old+1:n_r_max_old+n_r_max_f), &
                  &                                       work(k+1, n_r_max_old+1:n_r_max_old+n_r_max_f), kind=cp)
                  k=k+2
               end do
            end do
            deallocate(work)

            !!WARNING:: costf1 function is distributed in lm!call rscheme_3D%costf1(load_b0pol,1, lm_max, n_r_max_3D)
         end if
         call scatter_from_rank0_to_lmloc(load_b0pol,  b_LMloc)
         call scatter_from_rank0_to_lmloc(load_b0tor, aj_LMloc)

         !-- bringing the loaded b0 to the physical radial grid
         call rscheme_3D%costf1(b_LMloc, lmStart, lmStop, n_r_max_3D)
         call rscheme_3D%costf1(aj_LMloc, lmStart, lmStop, n_r_max_3D)

      else if ( init_B == 1010 ) then  ! l=1,m=0 decay mode , radial bessel function of the first kind and order 0!!
         ! WARNING:: the order of the bessel function used should match l
         !           --> can build complex (l,m) using different bessel functions an p-roots: see Gillet 2011
         if( (lm10>=lmStart) .and. (lm10<=lmStop) ) then ! select processor
            do n_r=1,n_r_max_3D
               x = pi*r_3D(n_r)/r_cmb!pi*(r_3D(n_r)-r_icb+1e-4_cp)/(r_cmb-r_icb)!
               b1(n_r) = (sin(x)/x)!/x - cos(x)/x)!
               b_LMloc(lm10,n_r)=amp_B*b1(n_r)
            end do
         end if
      !-- Bessel function initialization of all (l,m) modes exept (l=0,m=0):
         do lm=max(lmStart,2),lmStop
            l1=lo_map%lm2l(lm)
            m1=lo_map%lm2m(lm)
            if ( l1 > 0 ) then
               bR=amp_B/(real(l1,cp))**8.
               bI=amp_B/(real(l1,cp))**8.
            else
               bR=0.0_cp
               bI=0.0_cp
            end if
            if ( m1 == 0 ) bI=0.0_cp
            do n_r=1,n_r_max_3D
               b_LMloc(lm,n_r)=b_LMloc(lm,n_r) + cmplx(bR*b1(n_r),bI*b1(n_r),kind=cp)
            end do
         end do
      end if

      if ( l_mag_inertia ) then
         call get_dr(b_LMloc, db_LMloc, lmStart, lmStop,&
              &       n_r_max_3D, rscheme_3D)

         call transp_lm2r(lm2r_fields, b_LMloc, b_Rloc)
         call transp_lm2r(lm2r_fields, db_LMloc, db_Rloc)
         call transp_lm2r(lm2r_fields, aj_LMloc, aj_Rloc)

#ifdef WITH_SHTNS
         do n_r=nRstart3D,nRstop3D
            call torpol_to_spat(b_Rloc(:,n_r), db_Rloc(:,n_r), aj_Rloc(:,n_r), n_r, &
                 &              Br_shell(:,:), Bt_shell(:,:), Bp_shell(:,:))
            B0_3D_Rloc(:,:,n_r) = sqrt(Br_shell(:,:)**2 + Bt_shell(:,:)**2 + Bp_shell(:,:)**2)
         end do
#endif
      end if

   end subroutine initB_3D
!----------------------------------------------------------------------------------
   subroutine initB0_3D(b_LMloc, aj_LMloc, init_B0, amp_B0, lmax_trunc_B0)
      !
      ! Purpose of this subroutine is to initialize a background magnetic field
      ! according to the control parameters init_B0.
      ! Can be used for magneto convection or benchmark purposes
      !

      !-- Input variables
      real(cp),    intent(in) :: amp_B0
      integer,     intent(in) :: init_B0
      integer,     intent(in) :: lmax_trunc_B0

      !-- Output variables:
      complex(cp), intent(inout) :: b_LMloc(lmStart:lmStop, n_r_max_3D)
      complex(cp), intent(inout) :: aj_LMloc(lmStart:lmStop, n_r_max_3D)

      !-- Local variables:
      integer :: n_r
      real(cp) :: b_pol, b_tor
      !-- poloidal and toroidal arrays for grid_space quantities computation
      complex(cp) :: b0_LMloc(lmStart:lmStop, n_r_max_3D)
      complex(cp) :: db0_LMloc(lmStart:lmStop, n_r_max_3D)
      complex(cp) :: ddb0_LMloc(lmStart:lmStop, n_r_max_3D)
      complex(cp) :: aj0_LMloc(lmStart:lmStop, n_r_max_3D)
      complex(cp) :: dj0_LMloc(lmStart:lmStop, n_r_max_3D)
      complex(cp) :: b0_Rloc(lm_max,nRstart3D:nRstop3D)
      complex(cp) :: db0_Rloc(lm_max,nRstart3D:nRstop3D)
      complex(cp) :: ddb0_Rloc(lm_max,nRstart3D:nRstop3D)
      complex(cp) :: aj0_Rloc(lm_max,nRstart3D:nRstop3D)
      complex(cp) :: dj0_Rloc(lm_max,nRstart3D:nRstop3D)
      !-- arrays to read complex b0
      integer :: llmm2lm(l_max+1,m_max_3D+1)
      real(cp) :: load_head(8)
      real(cp), allocatable :: work(:,:)!, work_left(:)
      complex(cp) :: load_b0pol(lm_max, n_r_max_3D)
      complex(cp) :: load_b0tor(lm_max, n_r_max_3D)
      integer ::  k, l, m, lm
      integer :: l_max_old, l_max_f, lm_max_old
      integer :: n_start_file, n_r_max_old, n_r_max_f

      integer :: lm10, lm20
      logical :: startfile_does_exist

      !-- Set the values of the perturbation field to zero
      if ( .not. l_start_file ) then
         b_LMloc(:,:) = zero
         aj_LMloc(:,:) = zero
      end if

      lm10 = lo_map%lm2(1,0)
      lm20 = lo_map%lm2(2,0)

      b0_LMloc(:,:) = zero
      aj0_LMloc(:,:) = zero
      if ( init_B0 == 1 ) then  ! Simple Background field
      ! l=1,m0 poloidal field
      ! which is potential field at r_cmb
         if( (lm10>=lmStart) .and. (lm10<=lmStop) ) then ! select processor
            b_pol=amp_B0*5.0_cp*half*sqrt(third*pi)*r_icb**2
            do n_r=1,n_r_max_3D
               b0_LMloc(lm10,n_r)=b0_LMloc(lm10,n_r)+b_pol*(r_3D(n_r)/r_icb)**2 * &
               &                     ( one - three/5.0_cp*(r_3D(n_r)/r_cmb)**2 )
            end do
         end if

      else if ( init_B0 == 2 ) then  ! l=1,m=0 analytical toroidal field
      ! with a maximum of amp_B0 at mid-radius
      ! between r_icb and r_cmb for an insulating
      ! inner core and at r_cmb/2 for a conducting
      ! inner core
         if( (lm10>=lmStart) .and. (lm10<=lmStop) ) then ! select processor
!         if ( lmStartB(rank+1) <= lm10 .and. lmStopB(rank+1) >= lm10 ) then ! select processor
            b_tor=-two*amp_B0*sqrt(third*pi)  ! minus sign makes phi comp. > 0
            do n_r=1,n_r_max_3D
               aj0_LMloc(lm10,n_r)=aj0_LMloc(lm10,n_r) + b_tor*r_3D(n_r)*sin(pi*(r_3D(n_r)-r_icb))
            end do
         end if

      else if ( init_B0 == 3 ) then
         ! l=2,m=0 toroidal field and l=1,m=0 poloidal field
         ! toroidal field has again its maximum of amp_B
         ! at mid-radius between r_icb and r_cmb for an
         ! insulating inner core and at r_cmb/2 for a
         ! conducting inner core
         ! The outer core poloidal field is defined by
         ! a homogeneous  current density, its maximum at
         ! the ICB is set to amp_B.
         ! The inner core poloidal field is chosen accordingly.
         if( (lm10>=lmStart) .and. (lm10<=lmStop) ) then ! select processor
            b_tor=-four*third*amp_B*sqrt(pi/5.0_cp)
            b_pol= amp_B*sqrt(three*pi)/four
            do n_r=1,n_r_max_3D
               b0_LMloc(lm10,n_r)=b0_LMloc(lm10,n_r)+b_pol * (       &
               &     r_3D(n_r)**3 - four*third*r_cmb*r_3D(n_r)**2  &
               &     + third*r_icb**4/r_3D(n_r)                )
            end do
         end if

         if( (lm20>=lmStart) .and. (lm20<=lmStop) ) then ! select processor
            b_tor=-four*third*amp_B*sqrt(pi/5.0_cp)
            b_pol=amp_B*sqrt(three*pi)/four
            do n_r=1,n_r_max_3D
               aj0_LMloc(lm20,n_r)=aj0_LMloc(lm20,n_r) + b_tor*r_3D(n_r)*sin(pi*(r_3D(n_r)-r_icb))
            end do
         end if

      else if ( init_B0 == -1 ) then  ! Read a file containing a background field of any complexity
      ! b0 format:: Spherical harmonics and chebychev polynomials
      !             (2*lm_max,2*n_r_max)=(complex stored in real array, bpol+btor)
      ! !WARNING!:: Fortran format =' formatted' ('unformatted' is for binary files; not implemented yet)
         if ( rank == 0 ) then
            inquire(file=start_file_b0, exist=startfile_does_exist)

            if ( startfile_does_exist ) then
               open(newunit=n_start_file, file=start_file_b0, status='old', &
               &    form='formatted', access='stream')
               print*, "! Reading background field B0 stored in ", start_file_b0
            else
               call abortRun('! The restart file for B_0 does not exist !')
            end if

            read(n_start_file, *) load_head ! unknown length
            l_max_old = int(load_head(1)) ! l_max from checkpoint
            lm_max_old = int(load_head(2)) ! lm_max from checkpoint
            n_r_max_old = int(load_head(4))! n_r_max from checkpoint
            !print*, load_head(:8)
            !print*, l_max_old, lm_max_old, n_r_max_old
            !-- no need for reading the rest of the line apparently...
            !allocate(work_left(2*n_r_max_old-8))
            !read(n_start_file, *) work_left

            !-- extracting bpol/btor from the file
            allocate(work(2*lm_max_old, 2*n_r_max_old))
            !read(n_start_file, *) work ! apparently not possible to get the full block at once
            do lm=1,2*lm_max_old
               read(n_start_file, *) work(lm,:)
            end do
            !print*, work(2,:8)
            close(n_start_file)

            !-- lookup table for converting lm to l,m from checkpoint file
            lm=0
            do m=0,m_max_3D
               do l=m,l_max
                  lm=lm+1
                  llmm2lm(l+1,m+1)=lm
               end do
            end do

            !-- adapting to the current resolution
            l_max_f = min(lmax_trunc_B0, l_max_old)
            n_r_max_f = min(n_r_max_3D, n_r_max_old)
            !-- De-sorting the coefficient for the Back transforms
            k=1
            do l=0,l_max_f !-- WARNING: must start at 0 because (l,m)=(0,0) is stored as well
               do m=0,l
                  lm = llmm2lm(l+1,m+1)
                  !lm = lo_map%lm2(l,m)
                  load_b0pol(lm,:n_r_max_f) = amp_B0*cmplx(work(k  ,             1:            n_r_max_f), &
                  &                                        work(k+1,             1:            n_r_max_f), kind=cp)
                  load_b0tor(lm,:n_r_max_f) = amp_B0*cmplx(work(k  , n_r_max_old+1:n_r_max_old+n_r_max_f), &
                  &                                        work(k+1, n_r_max_old+1:n_r_max_old+n_r_max_f), kind=cp)
                  k=k+2
               end do
            end do
            deallocate(work)

            !!WARNING:: costf1 function is distributed in lm!call rscheme_3D%costf1(load_b0pol,1, lm_max, n_r_max_3D)
         end if
         call scatter_from_rank0_to_lmloc(load_b0pol,  b0_LMloc)
         call scatter_from_rank0_to_lmloc(load_b0tor, aj0_LMloc)

         !-- bringing the loaded b0 to the physical radial grid
         call rscheme_3D%costf1(b0_LMloc, lmStart, lmStop, n_r_max_3D)
         call rscheme_3D%costf1(aj0_LMloc, lmStart, lmStop, n_r_max_3D)

      else  ! Simple Background field !--> Same as init_B0 = 1 for the moment!!!
      ! l=1,m0 poloidal field
      ! which is potential field at r_cmb
         if( (lm10>=lmStart) .and. (lm10<=lmStop) ) then ! select processor
            b_pol=amp_B0*5.0_cp*half*sqrt(third*pi)*r_icb**2
            do n_r=1,n_r_max_3D
               b0_LMloc(lm10,n_r)=b0_LMloc(lm10,n_r)+b_pol*(r_3D(n_r)/r_icb)**2 * &
               &                     ( one - three/5.0_cp*(r_3D(n_r)/r_cmb)**2 )
            end do
         end if

      end if

      call get_ddr(b0_LMloc, db0_LMloc, ddb0_LMloc, lmStart, lmStop,&
           &       n_r_max_3D, rscheme_3D)
      call get_dr(aj0_LMloc, dj0_LMloc, lmStart, lmStop,&
           &      n_r_max_3D, rscheme_3D)

      call transp_lm2r(lm2r_fields, b0_LMloc, b0_Rloc)
      call transp_lm2r(lm2r_fields, db0_LMloc, db0_Rloc)
      call transp_lm2r(lm2r_fields, ddb0_LMloc, ddb0_Rloc)
      call transp_lm2r(lm2r_fields, aj0_LMloc, aj0_Rloc)
      call transp_lm2r(lm2r_fields, dj0_LMloc, dj0_Rloc)

#ifdef WITH_SHTNS
      do n_r=nRstart3D,nRstop3D
         call torpol_to_spat(b0_Rloc(:,n_r), db0_Rloc(:,n_r), aj0_Rloc(:,n_r), n_r, &
              &              B0r_3D_Rloc(:,:,n_r), B0t_3D_Rloc(:,:,n_r), B0p_3D_Rloc(:,:,n_r))
         call torpol_to_curl_spat(b0_Rloc(:,n_r), ddb0_Rloc(:,n_r), aj0_Rloc(:,n_r), dj0_Rloc(:,n_r), &
              &                   n_r, curlB0r_3D_Rloc(:,:,n_r), curlB0t_3D_Rloc(:,:,n_r), curlB0p_3D_Rloc(:,:,n_r))
      end do
#endif

   end subroutine initB0_3D
!----------------------------------------------------------------------------------
end module init_fields
