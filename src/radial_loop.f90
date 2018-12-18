module rloop

   use precision_mod
   use parallel_mod
   use constants, only: ci, one, half
   use mem_alloc, only: bytes_allocated
   use namelists, only: ek, tadvz_fac, CorFac, l_heat, l_chem
   use radial_functions, only: or1, r, beta, oheight, dtcond, ekpump
   use blocking, only: nRstart, nRstop
   use truncation, only: n_m_max, n_phi_max, idx2m, m2idx
   use courant_mod, only: courant
   use fourier, only: fft, ifft
   use useful, only: cc22real
   use time_schemes, only: type_tscheme
   use timers_mod, only: timers_type

   implicit none

   private

   real(cp), allocatable :: us_grid(:), up_grid(:)
   real(cp), allocatable :: om_grid(:), temp_grid(:), xi_grid(:)
   real(cp), allocatable :: usT_grid(:), upT_grid(:)
   real(cp), allocatable :: usXi_grid(:), upXi_grid(:)
   real(cp), allocatable :: usOm_grid(:), upOm_grid(:)

   public :: radial_loop, initialize_radial_loop, finalize_radial_loop

contains

   subroutine initialize_radial_loop(n_phi_max)

      !-- Input variable
      integer, intent(in) :: n_phi_max

      !-- Local variable:
      integer :: n_nphi_arrays

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

      end subroutine initialize_radial_loop
   !------------------------------------------------------------------------------
      subroutine finalize_radial_loop

         if ( l_heat ) deallocate (usT_grid, upT_grid,  temp_grid )
         if ( l_chem ) deallocate (usXi_grid, upXi_grid,  xi_grid )
         deallocate( usOm_grid, upOm_grid, us_grid, up_grid, om_grid )

      end subroutine finalize_radial_loop
   !------------------------------------------------------------------------------
      subroutine radial_loop(us_Rloc, up_Rloc, om_Rloc, temp_Rloc, xi_Rloc,   &
                 &           dtempdt_Rloc, dVsT_Rloc, dxidt_Rloc, dVsXi_Rloc, &
                 &           dpsidt_Rloc, dVsOm_Rloc, dtr_Rloc, dth_Rloc,     &
                 &           timers, tscheme)

         !-- Input variables
         complex(cp),         intent(in) :: us_Rloc(n_m_max, nRstart:nRstop)
         complex(cp),         intent(in) :: up_Rloc(n_m_max, nRstart:nRstop)
         complex(cp),         intent(in) :: om_Rloc(n_m_max, nRstart:nRstop)
         complex(cp),         intent(in) :: temp_Rloc(n_m_max, nRstart:nRstop)
         complex(cp),         intent(in) :: xi_Rloc(n_m_max, nRstart:nRstop)
         class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp),       intent(out) :: dpsidt_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),       intent(out) :: dtempdt_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),       intent(out) :: dVsT_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),       intent(out) :: dxidt_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),       intent(out) :: dVsXi_Rloc(n_m_max, nRstart:nRstop)
      complex(cp),       intent(out) :: dVsOm_Rloc(n_m_max, nRstart:nRstop)
      real(cp),          intent(out) :: dtr_Rloc(nRstart:nRstop)
      real(cp),          intent(out) :: dth_Rloc(nRstart:nRstop)
      type(timers_type), intent(inout) :: timers

      !-- Local variables
      real(cp) :: usom, runStart, runStop
      complex(cp) :: us_fluct
      integer :: n_r, n_phi, n_m, m, idx_m0

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

         !-----------------
         !-- Bring data on the grid
         !-----------------
         runStart = MPI_Wtime()
         call ifft(us_Rloc(:,n_r), us_grid)
         call ifft(up_Rloc(:,n_r), up_grid)
         call ifft(om_Rloc(:,n_r), om_grid)
         runStop = MPI_Wtime()
         if ( l_heat ) call ifft(temp_Rloc(:,n_r), temp_grid)
         if ( l_chem ) call ifft(xi_Rloc(:,n_r), xi_grid)
         if ( runStop > runStart ) then
            timers%fft = timers%fft + (runStop-runStart)
            timers%n_fft_calls = timers%n_fft_calls + 3
         end if

         !-- Courant condition
         if ( tscheme%istage == 1 ) then
            call courant(n_r, dtr_Rloc(n_r), dth_Rloc(n_r), us_grid, up_grid, &
                 &       tscheme%courfac)
         end if

         !-- Get nonlinear products
         if ( l_heat ) then
            do n_phi=1,n_phi_max
               usT_grid(n_phi) =r(n_r)*us_grid(n_phi)*temp_grid(n_phi)
               upT_grid(n_phi) =       up_grid(n_phi)*temp_grid(n_phi)
            end do
         end if
         if ( l_chem ) then
            do n_phi=1,n_phi_max
               usXi_grid(n_phi) =r(n_r)*us_grid(n_phi)*xi_grid(n_phi)
               upXi_grid(n_phi) =       up_grid(n_phi)*xi_grid(n_phi)
            end do
         end if
         do n_phi=1,n_phi_max
            usOm_grid(n_phi)=r(n_r)*us_grid(n_phi)*  om_grid(n_phi)
            upOm_grid(n_phi)=       up_grid(n_phi)*  om_grid(n_phi)
         end do

         !-- Bring data back on the spectral domain
         runStart = MPI_Wtime()
         call fft(upOm_grid, dpsidt_Rloc(:,n_r))
         call fft(usOm_grid, dVsOm_Rloc(:,n_r))
         runStop = MPI_Wtime()
         if ( l_heat ) then
            call fft(upT_grid, dtempdt_Rloc(:,n_r))
            call fft(usT_grid, dVsT_Rloc(:,n_r))
         end if
         if ( l_chem ) then
            call fft(upXi_grid, dxidt_Rloc(:,n_r))
            call fft(usXi_grid, dVsXi_Rloc(:,n_r))
         end if
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
            else
               dpsidt_Rloc(n_m,n_r)=-or1(n_r)*ci*m*dpsidt_Rloc(n_m,n_r)
            end if
         end do

      end do

   end subroutine radial_loop
!------------------------------------------------------------------------------
end module rloop
