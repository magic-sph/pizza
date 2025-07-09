module rloop

   use precision_mod
   use parallel_mod
   use constants, only: ci, one, half, pi, two, zero
   use mem_alloc, only: bytes_allocated
   use namelists, only: ek, tadvz_fac, CorFac, l_heat, l_chem, n_rings,    &
       &                radratio, amp_forcing, radius_forcing, dx_forcing, &
       &                forcing_type, tag, l_full_disk, l_phase_field,     &
       &                epsPhase, penaltyFac, phaseDiffFac, tmelt
   use radial_functions, only: or1, r, beta, dtcond, ekpump, m_R, tcond
   use outputs_phase, only: calc_out_phase
   use blocking, only: nRstart, nRstop
   use truncation, only: n_m_max, n_phi_max, idx2m, m2idx, n_r_max
   use courant_mod, only: courant
   use fourier, only: fft, ifft
   use useful, only: cc22real
   use time_schemes, only: type_tscheme
   use timers_mod, only: timers_type
   use output_frames, only: write_snapshot_rloc

   implicit none

   private

   real(cp), allocatable :: us_grid(:), up_grid(:)
   real(cp), allocatable :: om_grid(:), temp_grid(:), xi_grid(:)
   real(cp), allocatable :: phi_grid(:), phiTerms(:)
   real(cp), allocatable :: usT_grid(:), upT_grid(:)
   real(cp), allocatable :: usXi_grid(:), upXi_grid(:)
   real(cp), allocatable :: usOm_grid(:), upOm_grid(:)
   complex(cp), allocatable :: forcing_Rloc(:,:)
   real(cp), public, allocatable :: amp_pump(:), x_pump(:), y_pump(:)

   public :: radial_loop, initialize_radial_loop, finalize_radial_loop

contains

   subroutine initialize_radial_loop(n_phi_max)

      !-- Input variables
      integer,  intent(in) :: n_phi_max

      !-- Local variable:
      real(cp), allocatable :: r_rings(:), phi_ring(:)
      real(cp) :: force(n_phi_max)
      real(cp) :: ricb, rcmb, dr, dphi, xgrid, ygrid, phi, dy_forcing, rad
      integer, allocatable :: n_phi_ring(:)
      integer :: n_nphi_arrays, nr, np, Ntot, npump, Npumps, i, j, Nx, Ny

      allocate( us_grid(n_phi_max), up_grid(n_phi_max), om_grid(n_phi_max) )
      allocate( usOm_grid(n_phi_max), upOm_grid(n_phi_max) )
      n_nphi_arrays = 5
      if ( l_heat ) then
         allocate( usT_grid(n_phi_max), upT_grid(n_phi_max) )
         allocate( temp_grid(n_phi_max) )
         n_nphi_arrays = n_nphi_arrays+3
      end if
      if ( l_chem ) then
         allocate( usXi_grid(n_phi_max), upXi_grid(n_phi_max) )
         allocate( xi_grid(n_phi_max) )
         n_nphi_arrays = n_nphi_arrays+3
      end if
      if ( l_phase_field ) then
         allocate( phi_grid(n_phi_max), phiTerms(n_phi_max) )
         n_nphi_arrays = n_nphi_arrays+2
      end if
      bytes_allocated = bytes_allocated+n_nphi_arrays*n_phi_max*SIZEOF_DEF_REAL

      us_grid(:)  =0.0_cp
      up_grid(:)  =0.0_cp
      om_grid(:)  =0.0_cp
      usOm_grid(:)=0.0_cp
      upOm_grid(:)=0.0_cp
      if ( l_heat ) then
         temp_grid(:)=0.0_cp
         usT_grid(:) =0.0_cp
         upT_grid(:) =0.0_cp
      end if
      if ( l_chem ) then
         xi_grid(:)  =0.0_cp
         usXi_grid(:)=0.0_cp
         upXi_grid(:)=0.0_cp
      end if
      if ( l_phase_field ) then
         phi_grid(:)=0.0_cp
         phiTerms(:)=0.0_cp
      end if

      if ( n_rings > 0 .or. dx_forcing > 0.0_cp ) then
         allocate(forcing_Rloc(n_m_max,nRstart:nRstop) )
         forcing_Rloc(:,:)=zero
         bytes_allocated = bytes_allocated+n_m_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
         ricb=radratio/(one-radratio)
         rcmb=one/(one-radratio)

         if ( index(forcing_type, 'POLAR') == 1 ) then

            allocate(r_rings(n_rings), n_phi_ring(n_rings))
            dr=(rcmb-ricb)/(n_rings+1)
            !-- Define rings
            do nr=1,n_rings
               r_rings(nr)=rcmb-dr-(nr-1)*dr
            end do

            !-- Determine the number of vortices
            Ntot = 0
            do nr=1,n_rings
               n_phi_ring(nr) = int(r_rings(nr)*two*pi/dr)
               !-- Make sure we get an even number per ring
               if ( mod(n_phi_ring(nr),2)/=0 ) n_phi_ring(nr)=n_phi_ring(nr)+(-1)**(nr-1)
               Ntot = Ntot + n_phi_ring(nr)
            end do

            !-- Get the coordinates
            allocate(x_pump(Ntot), y_pump(Ntot), amp_pump(Ntot))
            Npumps=0
            do nr=1,n_rings
               allocate(phi_ring(n_phi_ring(nr)))
               dphi = two*pi/n_phi_ring(nr)
               do np=1,n_phi_ring(nr)
                  phi_ring(np)=(np-1)*dphi
               end do
               x_pump(Npumps+1:Npumps+n_phi_ring(nr))=r_rings(nr)*cos(phi_ring(:))
               y_pump(Npumps+1:Npumps+n_phi_ring(nr))=r_rings(nr)*sin(phi_ring(:))
               deallocate(phi_ring)
               Npumps=Npumps+n_phi_ring(nr)
            end do

            !-- Set amplitude array
            do npump=1,Ntot
               amp_pump(npump)=(-1)**npump
            end do

            deallocate(r_rings, n_phi_ring)

         else if ( index(forcing_type, 'CARTESIAN') == 1 ) then

            Nx=int(two*rcmb/dx_forcing+1)
            Ny=Nx
            !-- Redefine dx (because of integer)
            dx_forcing=two*rcmb/(Nx-1)
            dy_forcing=dx_forcing

            !-- Determine the number of vortices
            Ntot = 0
            do i=1,Nx
               xgrid=-rcmb+(i-1)*dx_forcing
               do j=1,Ny
                  ygrid=-rcmb+(j-1)*dy_forcing
                  rad=sqrt(xgrid*xgrid+ygrid*ygrid)
                  if ( rad >= ricb+half*dx_forcing .and. rad <= rcmb-half*dx_forcing) then
                     Ntot=Ntot+1
                  end if
               end do
            end do

            !-- Get the coordinates
            allocate(x_pump(Ntot), y_pump(Ntot), amp_pump(Ntot))
            npump=1
            do i=1,Nx
               xgrid=-rcmb+(i-1)*dx_forcing
               do j=1,Ny
                  ygrid=-rcmb+(j-1)*dy_forcing
                  rad = sqrt(xgrid*xgrid+ygrid*ygrid)
                  if ( rad >= ricb+half*dx_forcing .and. rad <= rcmb-half*dx_forcing) then
                     x_pump(npump)  =xgrid
                     y_pump(npump)  =ygrid
                     amp_pump(npump)=(-1)**i*(-1)**j
                     npump=npump+1
                  end if
               end do
            end do

            !-- If number of vortices is odd remove one to make sur no bias is introduced
            if ( mod(Ntot,2)/=0 ) Ntot=Ntot-1
            !print*, 'Number of vortices:', Ntot
         end if


         do nr=nRstart,nRstop
            force(:)=0.0_cp
            do np=1,n_phi_max
               phi = (np-1)*two*pi/n_phi_max
               xgrid=r(nr)*cos(phi)
               ygrid=r(nr)*sin(phi)

               do npump=1,Ntot
                  force(np)=force(np)+amp_forcing*amp_pump(npump)            * &
                  &         exp(-(x_pump(npump)-xgrid)**2/radius_forcing**2) * &
                  &         exp(-(y_pump(npump)-ygrid)**2/radius_forcing**2)
               end do
               !if ( np == 1) print*, nr, force(np)
            end do

            call fft(force, forcing_Rloc(:,nr), n_m_max)
         end do

         !-- Write a snapshot which contain the layout of the forcing
         call write_snapshot_rloc('forcing.'//tag, 0.0_cp, forcing_Rloc)

      end if

   end subroutine initialize_radial_loop
!------------------------------------------------------------------------------
   subroutine finalize_radial_loop

      if ( n_rings > 0 ) deallocate( forcing_Rloc, x_pump, y_pump, amp_pump )
      if ( l_heat ) deallocate (usT_grid, upT_grid,  temp_grid )
      if ( l_chem ) deallocate (usXi_grid, upXi_grid,  xi_grid )
      if ( l_phase_field ) deallocate( phi_grid, phiTerms )
      deallocate( usOm_grid, upOm_grid, us_grid, up_grid, om_grid )

   end subroutine finalize_radial_loop
!------------------------------------------------------------------------------
   subroutine radial_loop(us_Rloc, up_Rloc, om_Rloc, temp_Rloc, xi_Rloc,   &
              &           phi_Rloc, dtempdt_Rloc, dVsT_Rloc, dxidt_Rloc,   &
              &           dVsXi_Rloc, dphidt_Rloc, dpsidt_Rloc, dVsOm_Rloc,&
              &           dtr_Rloc, dth_Rloc, timers, tscheme, l_log)

      !-- Input variables
      complex(cp),         intent(in) :: us_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),         intent(in) :: up_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),         intent(in) :: om_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),         intent(in) :: temp_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),         intent(in) :: phi_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),         intent(in) :: xi_Rloc(n_m_max, nRstart:nRstop)
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: l_log

      !-- Output variables
      complex(cp),       intent(out) :: dpsidt_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),       intent(out) :: dtempdt_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),       intent(out) :: dVsT_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),       intent(out) :: dxidt_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),       intent(out) :: dphidt_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),       intent(out) :: dVsXi_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),       intent(out) :: dVsOm_Rloc(n_m_max, nRstart:nRstop)
      real(cp),          intent(out) :: dtr_Rloc(nRstart:nRstop)
      real(cp),          intent(out) :: dth_Rloc(nRstart:nRstop)
      type(timers_type), intent(inout) :: timers

      !-- Local variables
      real(cp) :: usom, runStart, runStop, upPhi
      complex(cp) :: us_fluct
      integer :: n_r, n_m, m, idx_m0

      idx_m0 = m2idx(0)

      do n_r=nRstart,nRstop

         !-- Calculate Reynolds stress for axisymmetric equation
         usom = 0.0_cp
         do n_m=1,n_m_max
            m = idx2m(n_m)

            if ( m == 0 ) then ! Add first order contribution for this term
                               ! us(m=0) * vortz(m=0)
               us_fluct = half*ek*CorFac*ekpump(n_r)*up_Rloc(n_m,n_r)
            else
               us_fluct = us_Rloc(n_m,n_r)
            end if
            usom = usom+cc22real(us_fluct,om_Rloc(n_m,n_r),m)
         end do

         !-- Penalty term for the axisymmetric equation
         if ( l_phase_field ) then
            upPhi = 0.0_cp
            do n_m=1,n_m_max
               upPhi = upPhi+cc22real(up_Rloc(n_m,n_r),phi_Rloc(n_m,n_r),m)
            end do
         end if

         !-----------------
         !-- Bring data on the grid
         !-----------------
         runStart = MPI_Wtime()
         call ifft(us_Rloc(:,n_r), us_grid, m_R(n_r))
         call ifft(up_Rloc(:,n_r), up_grid, m_R(n_r))
         call ifft(om_Rloc(:,n_r), om_grid, m_R(n_r))
         runStop = MPI_Wtime()
         if ( l_heat ) call ifft(temp_Rloc(:,n_r), temp_grid, m_R(n_r))
         if ( l_chem ) call ifft(xi_Rloc(:,n_r), xi_grid, m_R(n_r))
         if ( l_phase_field ) call ifft(phi_Rloc(:,n_r), phi_grid, m_R(n_r))
         if ( runStop > runStart ) then
            timers%fft = timers%fft + (runStop-runStart)
            timers%n_fft_calls = timers%n_fft_calls + 3
         end if

         !-- Compute diagnostics related to phase field
         if ( l_phase_field .and. l_log ) then
            call calc_out_phase(us_grid, up_grid, phi_grid, temp_grid, n_r)
         end if

         !-- Courant condition
         if ( tscheme%istage == 1 .and. (.not. l_full_disk .or. n_r /= n_r_max) ) then
            call courant(n_r, dtr_Rloc(n_r), dth_Rloc(n_r), us_grid, up_grid, &
                 &       tscheme%courfac)
         end if

         !-- Get nonlinear products
         if ( l_heat ) then
            usT_grid(:) =r(n_r)*us_grid(:)*temp_grid(:)
            upT_grid(:) =       up_grid(:)*temp_grid(:)
         end if
         if ( l_chem ) then
            usXi_grid(:) =r(n_r)*us_grid(:)*xi_grid(:)
            upXi_grid(:) =       up_grid(:)*xi_grid(:)
         end if
         usOm_grid(:)=r(n_r)*us_grid(:)*om_grid(:)
         upOm_grid(:)=       up_grid(:)*om_grid(:)

         if ( l_phase_field ) then
            phiTerms(:)=-one/epsPhase**2* phi_grid(:)*(one-phi_grid(:))*( &
            &           phaseDiffFac*(one-two*phi_grid(:))+temp_grid(:)+  &
            &           tcond(n_r)-tmelt)
            !usOm_grid(:)=usOm_grid(:)+r(n_r)*us_grid(:)*phi_grid(:)/ &
            !&            epsPhase**2/penaltyFac**2
            !upOm_grid(:)=upOm_grid(:)+up_grid(:)*phi_grid(:)/ &
            !&            epsPhase**2/penaltyFac**2
            usOm_grid(:)=usOm_grid(:)+r(n_r)*up_grid(:)*phi_grid(:)/ &
            &            epsPhase**2/penaltyFac**2
            upOm_grid(:)=upOm_grid(:)-us_grid(:)*phi_grid(:)/ &
            &            epsPhase**2/penaltyFac**2
         end if

         !-- Bring data back on the spectral domain
         runStart = MPI_Wtime()
         call fft(upOm_grid, dpsidt_Rloc(:,n_r), m_R(n_r))
         call fft(usOm_grid, dVsOm_Rloc(:,n_r), m_R(n_r))
         runStop = MPI_Wtime()
         if ( l_heat ) then
            call fft(upT_grid, dtempdt_Rloc(:,n_r), m_R(n_r))
            call fft(usT_grid, dVsT_Rloc(:,n_r), m_R(n_r))
         end if
         if ( l_chem ) then
            call fft(upXi_grid, dxidt_Rloc(:,n_r), m_R(n_r))
            call fft(usXi_grid, dVsXi_Rloc(:,n_r), m_R(n_r))
         end if
         if ( l_phase_field ) call fft(phiTerms, dphidt_Rloc(:,n_r), m_R(n_r))
         if ( runStop > runStart ) then
            timers%fft = timers%fft + (runStop-runStart)
            timers%n_fft_calls = timers%n_fft_calls + 2
         end if

         if ( l_heat ) then
            do n_m=1,n_m_max
               m = idx2m(n_m)
               dtempdt_Rloc(n_m,n_r)=-or1(n_r)*ci*m*dtempdt_Rloc(n_m,n_r) &
               &                     -(one-tadvz_fac)*beta(n_r)*or1(n_r)* &
               &                     dVsT_Rloc(n_m,n_r)
            end do
         end if
         if ( l_chem ) then
            do n_m=1,n_m_max
               m = idx2m(n_m)
               dxidt_Rloc(n_m,n_r)=-or1(n_r)*ci*m*dxidt_Rloc(n_m,n_r)
            end do
         end if
         do n_m=1,n_m_max
            m = idx2m(n_m)
            if ( m == 0 ) then
               dpsidt_Rloc(n_m,n_r)=-usom
               if ( l_phase_field ) then ! Penalty term for the <uphi> equation
                  dpsidt_Rloc(n_m,n_r)=dpsidt_Rloc(n_m,n_r)-upPhi / &
                  &                    epsPhase**2/penaltyFac**2
               end if
            else
               dpsidt_Rloc(n_m,n_r)=-or1(n_r)*ci*m*dpsidt_Rloc(n_m,n_r)
            end if
         end do

         if ( amp_forcing > 0.0_cp ) then
            do n_m=1,n_m_max
               m = idx2m(n_m)
               if ( m > 0 ) dpsidt_Rloc(n_m,n_r)=dpsidt_Rloc(n_m,n_r) + &
               &                                 forcing_Rloc(n_m,n_r)
            end do
         end if
      end do

   end subroutine radial_loop
!------------------------------------------------------------------------------
end module rloop
