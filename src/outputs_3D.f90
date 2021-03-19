module outputs_3D
   !
   ! This module computes the outputs that are related to fiels expressed in 3-D
   !

   use parallel_mod
   use precision_mod
   use spectra_3D, only: spectra_3D_type
   use constants, only: zero, two, half, pi, osq4pi
   use communications, only: reduce_radial_on_rank
   use truncation_3D, only: n_r_max_3D, n_theta_max, n_phi_max_3D
   use blocking, only: lmStart, lmStop, nRstart3D, nRstop3D
   use blocking_lm, only: lo_map
   use namelists, only: l_heat_3D, l_mag_3D, tag, DyMagFac
   use radial_functions, only: rscheme_3D, dtcond_3D, tcond_3D, r_3D, or2_3D
   use radial_der, only: get_dr
   use mean_sd, only: mean_sd_type
   use time_schemes, only: type_tscheme
   use integration, only: rInt_R
   use useful, only: round_off, getMSD2, cc2real

   implicit none

   private

   type(mean_sd_type) :: tempR
   type(mean_sd_type) :: bpolR
   type(mean_sd_type) :: btorR
   type(spectra_3D_type), public :: spec_3D
   integer :: n_heat_file, n_mag_file, n_calls
   real(cp) :: timeLast_rad, timeAvg_rad
   real(cp) :: timeAvg_spec

   public :: initialize_outputs_3D, finalize_outputs_3D, write_outputs_3D

contains

   subroutine initialize_outputs_3D

      character(len=144) :: file_name

      timeAvg_rad = 0.0_cp
      timeAvg_spec= 0.0_cp

      if ( rank == 0 ) then
         if ( l_heat_3D ) then
            file_name = 'heat_3D.'//tag
            open(newunit=n_heat_file, file=file_name, status='new')
         end if
         if ( l_mag_3D ) then
            file_name = 'e_mag_3D.'//tag
            open(newunit=n_mag_file, file=file_name, status='new')
         end if
         call tempR%initialize(1,n_r_max_3D)
         call bpolR%initialize(1,n_r_max_3D)
         call btorR%initialize(1,n_r_max_3D)
         n_calls     =0
         timeLast_rad=0.0_cp
      end if

      call spec_3D%initialize()

   end subroutine initialize_outputs_3D
!---------------------------------------------------------------------------------
   subroutine finalize_outputs_3D

      if ( rank == 0 ) then
         call btorR%finalize()
         call bpolR%finalize()
         call tempR%finalize()
         if ( l_heat_3D ) close(n_heat_file)
         if ( l_mag_3D ) close(n_mag_file)
      end if

      call spec_3D%finalize()

   end subroutine finalize_outputs_3D
!---------------------------------------------------------------------------------
   subroutine write_outputs_3D(time, tscheme, l_log, l_stop_time, temp_3D, &
   &                           b_3D, db_3D, aj_3D)

      !-- Input variables
      real(cp),            intent(in) :: time
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: l_log
      logical,             intent(in) :: l_stop_time
      complex(cp),         intent(in) :: temp_3D(lmStart:lmStop,n_r_max_3D)
      complex(cp),         intent(in) :: b_3D(lmStart:lmStop,n_r_max_3D)
      complex(cp),         intent(in) :: db_3D(lmStart:lmStop,n_r_max_3D)
      complex(cp),         intent(in) :: aj_3D(lmStart:lmStop,n_r_max_3D)

      !-- Local variables
      real(cp) :: bpol(n_r_max_3D), btor(n_r_max_3D)

      timeAvg_rad  = timeAvg_rad  + tscheme%dt(1)
      timeAvg_spec = timeAvg_spec + tscheme%dt(1)

      if ( l_log ) then
         call write_time_series(time, temp_3D, b_3D, db_3D, aj_3D, bpol, btor)
         call get_radial_averages(timeAvg_rad, l_stop_time, temp_3D, bpol, btor)
      end if

      if ( l_mag_3D ) then
      !-- Calculate spectra
      if ( l_log .or. spec_3D%l_calc ) then
         call spec_3D%compute_spectra_3D(timeAvg_spec, l_stop_time, b_3D, &
              &                         db_3D, aj_3D)
      end if
      end if

      !!-- Write spectra
      !if ( spec%l_calc ) then
      !   call spec%write_spectra(E_pol_l, E_tor_l)
      !end if

   end subroutine write_outputs_3D
!---------------------------------------------------------------------------------
   subroutine write_time_series(time, temp_3D, b_3D, db_3D, aj_3D, bpol, btor)

      !-- Input variables
      real(cp),    intent(in) :: time
      complex(cp), intent(in) :: temp_3D(lmStart:lmStop,n_r_max_3D)
      complex(cp), intent(in) :: b_3D(lmStart:lmStop,n_r_max_3D)
      complex(cp), intent(in) :: db_3D(lmStart:lmStop,n_r_max_3D)
      complex(cp), intent(in) :: aj_3D(lmStart:lmStop,n_r_max_3D)

      !-- Output variables
      real(cp),    intent(out) :: bpol(n_r_max_3D)
      real(cp),    intent(out) :: btor(n_r_max_3D)

      !-- Local variables
      real(cp) :: temp0(n_r_max_3D), dtemp0(n_r_max_3D)
      real(cp) :: NuTop, NuBot, tTop, tBot, NuDelta, beta_t
      !logical :: rank_has_lm10, rank_has_lm11
      !integer :: sr_tag, request1, request2
      integer :: n_r, lm, l, m
      integer :: lm10, lm11
      complex(cp) :: b10, b11
      real(cp) :: dLh, e_pol_tmp, e_tor_tmp, E_dip_axis
      real(cp) :: rad, theta_dip, phi_dip, Efac
      real(cp) :: E_pol      ! Volume averaged poloidal magnetic energy
      real(cp) :: E_tor      ! Volume averaged toroidal magnetic energy
      real(cp) :: E_pol_axis ! Volume averaged axisymmetric poloidal magnetic energy
      real(cp) :: E_tor_axis ! Volume averaged axisymmetric toroidal magnetic energy
      real(cp) :: E_cmb      ! Magnetic energy at the CMB
      real(cp) :: EDip       ! Relative magnetic energy of axial dipole
      real(cp) :: e_pol_r(n_r_max_3D), e_tor_r(n_r_max_3D)
      real(cp) :: e_pol_axis_r(n_r_max_3D), e_tor_axis_r(n_r_max_3D), e_dipole_axis_r(n_r_max_3D)

      if ( l_heat_3D .and. rank == 0  ) then
         temp0(:) = osq4pi*real(temp_3D(1, :))
         call get_dr(temp0, dtemp0, n_r_max_3D, rscheme_3D)

         !print*, 'toptemp, bottemp', temp0(1), temp0(n_r_max_3D)
         !print*, 'topgrad, botgrad', dtemp0(1), dtemp0(n_r_max_3D)
         !!print*, 'topdtes, botdtes',  osq4pi*real(dtemp_3D(1,1)),  osq4pi*real(dtemp_3D(1,n_r_max_3D))
         !print*, 'topcond, botcond', dtcond_3D(1), dtcond_3D(n_r_max_3D)
         !print*, ''

         !-- Nusselt at top and bottom
         NuTop = dtemp0(1)/(dtcond_3D(1)+epsilon(1.0_cp))
         NuBot = dtemp0(n_r_max_3D)/(dtcond_3D(n_r_max_3D)+epsilon(1.0_cp)) !-> Divide by 0 when T_cond = cste

         !-- Nusset based on the ratio of temperature contrast
         NuDelta = (tcond_3D(n_r_max_3D)-tcond_3D(1)) / &
         &         (temp0(n_r_max_3D)-temp0(1)+epsilon(1.0_cp)) !-> Divide by 0 when 0 Temperature

         !-- Temperature gradient at mid-shell
         beta_t = dtemp0(n_r_max_3D/2)

         !-- Top and bottom temperature
         tTop = temp0(1)
         tBot = temp0(n_r_max_3D)

         write(n_heat_file, '(2P, ES20.12, 7ES16.8)') time, NuBot, NuTop,   &
         &                                            NuDelta, tBot, tTop,  &
         &                                            beta_t

      end if

      Efac     = 0.0_cp
      E_pol     = 0.0_cp
      E_tor     = 0.0_cp
      E_pol_axis= 0.0_cp
      E_tor_axis= 0.0_cp
      E_cmb     = 0.0_cp
      EDip      = 0.0_cp
      if ( l_mag_3D ) then
         do n_r=1,n_r_max_3D
            e_pol_r(n_r)        = 0.0_cp
            e_tor_r(n_r)        = 0.0_cp
            e_pol_axis_r(n_r)   = 0.0_cp
            e_tor_axis_r(n_r)   = 0.0_cp
            e_dipole_axis_r(n_r)= 0.0_cp
            do lm=max(2,lmStart),lmStop
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               dLh = real(l*(l+1),cp)

               e_pol_tmp= dLh*(dLh*or2_3D(n_r)*cc2real( b_3D(lm,n_r),m) &
               &                              +cc2real(db_3D(lm,n_r),m) )
               e_tor_tmp= dLh*cc2real(aj_3D(lm,n_r),m)

               if ( m == 0 ) then  ! axisymmetric part
                  if ( l == 1 ) then
                     e_dipole_axis_r(n_r)=e_pol_tmp
                  end if
                  e_pol_axis_r(n_r)=e_pol_axis_r(n_r) + e_pol_tmp
                  e_tor_axis_r(n_r)=e_tor_axis_r(n_r) + e_tor_tmp
                  !if ( mod(l,2) == 1 ) then
                  !   e_pol_eaxis_r(n_r)=e_pol_eaxis_r(n_r)+e_pol_tmp
                  !else
                  !   e_tor_eaxis_r(n_r)=e_tor_eaxis_r(n_r)+e_tor_tmp
                  !end if
               else
                  !if ( l == 1 ) e_dipole_r(n_r)=e_dipole_r(n_r)+e_pol_tmp
                  e_pol_r(n_r)=e_pol_r(n_r) + e_pol_tmp
                  e_tor_r(n_r)=e_tor_r(n_r) + e_tor_tmp
               end if
               !if ( mod(l+m,2) == 1 ) then
               !   e_pol_es_r(n_r)=e_pol_es_r(n_r) + e_pol_tmp
               !else
               !   e_tor_es_r(n_r)=e_tor_es_r(n_r) + e_tor_tmp
               !end if

            end do    ! do loop over lms in block
            e_pol_r(n_r)=e_pol_r(n_r)+e_pol_axis_r(n_r)
            e_tor_r(n_r)=e_tor_r(n_r)+e_tor_axis_r(n_r)
         end do    ! radial grid points
         !-- A reduction is needed here
         call reduce_radial_on_rank(e_pol_r, n_r_max_3D, 0)
         call reduce_radial_on_rank(e_tor_r, n_r_max_3D, 0)
         call reduce_radial_on_rank(e_pol_axis_r, n_r_max_3D, 0)
         call reduce_radial_on_rank(e_tor_axis_r, n_r_max_3D, 0)
         call reduce_radial_on_rank(e_dipole_axis_r, n_r_max_3D, 0)
         bpol(:) = e_pol_r(:)
         btor(:) = e_tor_r(:)

         if ( rank == 0 ) then
            !-- Get Values at CMB: n_r=1 == n_r_cmb
            E_cmb      =e_pol_r(1)+e_tor_r(1)
            E_pol      =rInt_R(e_pol_r,r_3D,rscheme_3D)
            E_tor      =rInt_R(e_tor_r,r_3D,rscheme_3D)
            E_pol_axis =rInt_R(e_pol_axis_r,r_3D,rscheme_3D)
            E_tor_axis =rInt_R(e_tor_axis_r,r_3D,rscheme_3D)
            E_dip_axis =rInt_R(e_dipole_axis_r,r_3D,rscheme_3D)

            Efac=half*DyMagFac!*eScale
            !print*, Efac, E_cmb, E_pol, E_tor, E_pol_axis, E_tor_axis, E_dip_axis
            E_cmb      =Efac*E_cmb
            E_pol      =Efac*E_pol
            E_tor      =Efac*E_tor
            E_pol_axis =Efac*E_pol_axis
            E_tor_axis =Efac*E_tor_axis
            E_dip_axis =Efac*E_dip_axis

            lm10=lo_map%lm2(1,0)
            lm11=lo_map%lm2(1,1)

#ifdef TOTO
            ! some arbitrary send recv tag
            sr_tag=18657
            rank_has_lm10=.false.
            rank_has_lm11=.false.
            if ( (lm10 >= lmStart) .and. (lm10 <= lmStop) ) then
               b10=b_3D(lm10,1)
               if (rank /= 0) then
                  call MPI_Send(b10,1,MPI_DEF_COMPLEX,0,sr_tag,MPI_COMM_WORLD,ierr)
               end if
               rank_has_lm10=.true.
            end if
            if ( lm11 > 0 ) then
               if ( (lm11 >= lmStart) .and. (lm11 <= lmStop) ) then
                  b11=b_3D(lm11,1)
                  if (rank /= 0) then
                     call MPI_Send(b11,1,MPI_DEF_COMPLEX,0,sr_tag+1, &
                          &        MPI_COMM_WORLD,ierr)
                  end if
                  rank_has_lm11=.true.
               end if
            else
               b11=zero
               rank_has_lm11=.true.
            end if
            rad =180.0_cp/pi
            if (.not.rank_has_lm10) then
               call MPI_IRecv(b10,1,MPI_DEF_COMPLEX,MPI_ANY_SOURCE,&
                    &         sr_tag,MPI_COMM_WORLD,request1, ierr)
            end if
            if ( .not. rank_has_lm11 ) then
               call MPI_IRecv(b11,1,MPI_DEF_COMPLEX,MPI_ANY_SOURCE,&
                    &         sr_tag+1,MPI_COMM_WORLD,request2,ierr)
            end if
            if ( .not. rank_has_lm10 ) then
               call MPI_Wait(request1,status,ierr)
            end if
            if ( .not. rank_has_lm11 ) then
               call MPI_Wait(request2,status,ierr)
            end if
#endif

            theta_dip= rad*atan2(sqrt(two)*abs(b11),real(b10))
            if ( theta_dip < 0.0_cp ) theta_dip=180.0_cp+theta_dip
            if ( abs(b11) < 1.e-20_cp ) then
               phi_dip=0.0_cp
            else
               phi_dip=-rad*atan2(aimag(b11),real(b11))
            end if
            EDip      =E_dip_axis/(E_pol+E_tor+epsilon(1.0_cp)) !-> Divide by 0 when 0 field

            write(n_mag_file, '(1P, ES20.12, 4ES16.8)') time, E_pol, E_tor,   & 
            &                                           E_pol_axis, E_tor_axis
         end if

      end if ! l_mag_3D?

   end subroutine write_time_series
!---------------------------------------------------------------------------------
   subroutine get_radial_averages(timeAvg_rad, l_stop_time, temp_3D, bpol, btor)

      !-- Input variables
      real(cp),    intent(in) :: timeAvg_rad
      logical,     intent(in) :: l_stop_time
      complex(cp), intent(in) :: temp_3D(lmStart:lmStop,n_r_max_3D)
      real(cp),    intent(in) :: bpol(n_r_max_3D)
      real(cp),    intent(in) :: btor(n_r_max_3D)

      !-- Local variables
      real(cp) :: dtAvg, dat
      integer :: n_r, file_handle

      if ( rank == 0 ) then

         n_calls = n_calls+1
         dtAvg = timeAvg_rad-timeLast_rad

         do n_r=1,n_r_max_3D
            if ( l_heat_3D ) then
               dat = real(temp_3D(1,n_r))*osq4pi
               call getMSD2(tempR%mean(n_r), tempR%SD(n_r), dat, &
                    &       n_calls, dtAvg, timeAvg_rad)
            end if
            if ( l_mag_3D ) then
               dat = bpol(n_r)
               call getMSD2(bpolR%mean(n_r), bpolR%SD(n_r), dat, &
                    &       n_calls, dtAvg, timeAvg_rad)
               dat = btor(n_r)
               call getMSD2(btorR%mean(n_r), btorR%SD(n_r), dat, &
                    &       n_calls, dtAvg, timeAvg_rad)
            end if
         end do
         timeLast_rad = timeAvg_rad

         if ( l_stop_time ) then
            open(newunit=file_handle, file='radial_profiles_3D.'//tag)
            do n_r=1,n_r_max_3D
               if ( l_heat_3D ) then
                  tempR%SD(n_r)=sqrt(tempR%SD(n_r)/timeAvg_rad)
               end if
               if ( l_mag_3D ) then
                  bpolR%SD(n_r)=sqrt(bpolR%SD(n_r)/timeAvg_rad)
                  btorR%SD(n_r)=sqrt(btorR%SD(n_r)/timeAvg_rad)
               end if
               write(file_handle, '(es20.12, 6es16.8)') r_3D(n_r),        &
               &     round_off(tempR%mean(n_r)), round_off(tempR%SD(n_r)),&
               &     round_off(bpolR%mean(n_r)), round_off(bpolR%SD(n_r)),&
               &     round_off(btorR%mean(n_r)), round_off(btorR%SD(n_r))
            end do
            close(file_handle)
         end if

      end if

   end subroutine get_radial_averages
!---------------------------------------------------------------------------------
end module outputs_3D
