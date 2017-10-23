module rloop_integ

   use precision_mod
   use constants, only: ci, one, half, two, three, zero
   use mem_alloc, only: bytes_allocated
   use namelists, only: ek, tadvz_fac, r_cmb, l_non_rot
   use radial_functions, only: r, ekpump, or1
   use blocking, only: nRstart, nRstop
   use truncation, only: n_m_max, n_phi_max, idx2m, m2idx
   use courant_mod, only: courant
   use fourier, only: fft, ifft
   use useful, only: cc22real

   implicit none

   private

   real(cp), allocatable :: us_grid(:), up_grid(:)
   real(cp), allocatable :: om_grid(:), temp_grid(:)
   real(cp), allocatable :: usT_grid(:), upT_grid(:)
   real(cp), allocatable :: usOm_grid(:), upOm_grid(:)

   complex(cp), allocatable :: omMod(:), upMod(:), usMod(:)
   real(cp), allocatable :: psi_grid(:), omMod_grid(:), upMod_grid(:)
   real(cp), allocatable :: usMod_grid(:)

   public :: radial_loop_integ, initialize_radial_loop_integ, &
   &         finalize_radial_loop_integ

contains

   subroutine initialize_radial_loop_integ(n_phi_max)

      integer, intent(in) :: n_phi_max

      allocate( us_grid(n_phi_max), up_grid(n_phi_max) )
      allocate( temp_grid(n_phi_max), om_grid(n_phi_max)  )
      allocate( usT_grid(n_phi_max), upT_grid(n_phi_max) )
      allocate( usOm_grid(n_phi_max), upOm_grid(n_phi_max) )
      bytes_allocated = bytes_allocated+8*n_phi_max*SIZEOF_DEF_REAL

      us_grid(:)     =0.0_cp
      up_grid(:)     =0.0_cp
      temp_grid(:)   =0.0_cp
      om_grid(:)     =0.0_cp
      usT_grid(:)    =0.0_cp
      upT_grid(:)    =0.0_cp
      usOm_grid(:)   =0.0_cp
      upOm_grid(:)   =0.0_cp

      allocate( omMod(n_m_max), upMod(n_m_max), usMod(n_m_max) )
      allocate( psi_grid(n_phi_max), omMod_grid(n_phi_max) )
      allocate( usMod_grid(n_phi_max), upMod_grid(n_phi_max) )
      omMod(:)     =zero
      upMod(:)     =zero
      usMod(:)     =zero
      omMod_grid(:)=0.0_cp
      psi_grid(:)  =0.0_cp
      upMod_grid(:)=0.0_cp
      usMod_grid(:)=0.0_cp

   end subroutine initialize_radial_loop_integ
!------------------------------------------------------------------------------
   subroutine finalize_radial_loop_integ

      deallocate( usT_grid, upT_grid, usOm_grid, upOm_grid )
      deallocate( us_grid, up_grid, temp_grid, om_grid )

      deallocate( usMod_grid, upMod_grid, omMod_grid, psi_grid )
      deallocate( usMod, upMod, omMod )

   end subroutine finalize_radial_loop_integ
!------------------------------------------------------------------------------
   subroutine radial_loop_integ(psi_Rloc, dpsi_Rloc, d2psi_Rloc, us_Rloc,    &
              &                 up_Rloc, om_Rloc, temp_Rloc, dtempdt_Rloc,   &
              &           dVsT_Rloc, dpsidt_Rloc, dVsOm_Rloc, dtr_Rloc,      &
              &           dth_Rloc)

      !-- Input variables
      complex(cp), intent(in) :: us_Rloc(n_m_max, nRstart:nRstop)
      complex(cp), intent(in) :: up_Rloc(n_m_max, nRstart:nRstop)
      complex(cp), intent(in) :: om_Rloc(n_m_max, nRstart:nRstop)
      complex(cp), intent(in) :: psi_Rloc(n_m_max, nRstart:nRstop)
      complex(cp), intent(in) :: dpsi_Rloc(n_m_max, nRstart:nRstop)
      complex(cp), intent(in) :: d2psi_Rloc(n_m_max, nRstart:nRstop)
      complex(cp), intent(in) :: temp_Rloc(n_m_max, nRstart:nRstop)

      !-- Output variables
      complex(cp), intent(out) :: dpsidt_Rloc(n_m_max, nRstart:nRstop)
      complex(cp), intent(out) :: dtempdt_Rloc(n_m_max, nRstart:nRstop)
      complex(cp), intent(out) :: dVsT_Rloc(n_m_max, nRstart:nRstop)
      complex(cp), intent(out) :: dVsOm_Rloc(n_m_max, nRstart:nRstop)
      real(cp),    intent(out) :: dtr_Rloc(nRstart:nRstop)
      real(cp),    intent(out) :: dth_Rloc(nRstart:nRstop)

      !-- Local variables
      real(cp) :: usom, ro2ms2
      complex(cp) :: us_fluct
      integer :: n_r, n_phi, n_m, m, idx_m0

      idx_m0 = m2idx(0)

      do n_r=nRstart,nRstop

         !-- (r_cmb^2-s^2)
         ro2ms2 = r_cmb**2-r(n_r)**2

         usom = 0.0_cp
         do n_m=1,n_m_max
            m = idx2m(n_m)

            if ( l_non_rot ) then
               upMod(n_m) =-dpsi_Rloc(n_m,n_r)
               omMod(n_m) =-r(n_r)**2*            d2psi_Rloc(n_m,n_r)    &
               &           -r(n_r)*                dpsi_Rloc(n_m,n_r)    &
               &           +real(m,cp)*real(m,cp)*  psi_Rloc(n_m,n_r)
            else
               if ( m == 0 ) then
                  upMod(n_m)=ro2ms2*up_Rloc(n_m,n_r)
                  omMod(n_m)=ro2ms2**2*r(n_r)**2*om_Rloc(n_m,n_r)
               else
                  upMod(n_m)=-ro2ms2*dpsi_Rloc(n_m,n_r)+r(n_r)*psi_Rloc(n_m,n_r)
                  omMod(n_m)=-ro2ms2**2*r(n_r)**2*    d2psi_Rloc(n_m,n_r)    &
                  &          -ro2ms2*r(n_r)*(r_cmb**2-two*r(n_r)**2)*        &
                  &                                    dpsi_Rloc(n_m,n_r)    &
                  &          +( two*r_cmb**2*r(n_r)**2+ro2ms2**2*real(m,cp)* &
                  &             real(m,cp) )*           psi_Rloc(n_m,n_r)
               end if
            end if
            usMod(n_m) = ci*real(m,cp)*psi_Rloc(n_m,n_r)

            if ( m == 0 ) then ! Add first order contribution for this term
                               ! us(m=0) * vortz(m=0)
               us_fluct = half*ek*ekpump(n_r)*up_Rloc(n_m,n_r)
            else
               us_fluct = us_Rloc(n_m,n_r)
            end if

            !-- Calculate Reynolds stress for axisymmetric equation
            usom = usom+cc22real(us_fluct,omMod(n_m),m)

         end do

         !-----------------
         !-- Bring data on the grid
         !-----------------
         call ifft(us_Rloc(:,n_r), us_grid)
         call ifft(up_Rloc(:,n_r), up_grid)
         call ifft(temp_Rloc(:,n_r), temp_grid)

         call ifft(omMod, omMod_grid)
         call ifft(usMod, usMod_grid)
         call ifft(upMod, upMod_grid)
         call ifft(psi_Rloc(:,n_r), psi_grid)

         !-- Courant condition
         call courant(n_r, dtr_Rloc(n_r), dth_Rloc(n_r), us_grid, up_grid)

         !-- Get nonlinear products
         do n_phi=1,n_phi_max
            usT_grid(n_phi) =usMod_grid(n_phi)*temp_grid(n_phi)
            upT_grid(n_phi) =upMod_grid(n_phi)*temp_grid(n_phi)

            !usOm_grid(n_phi)=ro2ms2*r(n_r)**3*omMod_grid(n_phi)*usMod_grid(n_phi)
            usOm_grid(n_phi)=       usMod_grid(n_phi)*omMod_grid(n_phi)
            upOm_grid(n_phi)=r(n_r)*upMod_grid(n_phi)*omMod_grid(n_phi)
         end do

         !-- Bring data back on the spectral domain
         call fft(upT_grid, dtempdt_Rloc(:,n_r))
         call fft(usT_grid, dVsT_Rloc(:,n_r))
         call fft(upOm_grid, dpsidt_Rloc(:,n_r))
         call fft(usOm_grid, dVsOm_Rloc(:,n_r))

         do n_m=1,n_m_max
            m = idx2m(n_m)
            if ( l_non_rot ) then
               dtempdt_Rloc(n_m,n_r)=-ci*m*dtempdt_Rloc(n_m,n_r) 
               if ( m == 0 ) then
                  dpsidt_Rloc(n_m,n_r)=-usom
               else
                  dpsidt_Rloc(n_m,n_r)=-ci*m*dpsidt_Rloc(n_m,n_r)+ &
                  &                     three*dVsOm_Rloc(n_m,n_r)
                  dVsOm_Rloc(n_m, n_r)=r(n_r)*dVsOm_Rloc(n_m,n_r)
               end if

            else
               dtempdt_Rloc(n_m,n_r)=-ci*m*dtempdt_Rloc(n_m,n_r)              &
               &                     +(one-tadvz_fac)*r(n_r)*dVsT_Rloc(n_m,n_r)
               if ( m == 0 ) then
                  dpsidt_Rloc(n_m,n_r)=-usom
               else
                  dpsidt_Rloc(n_m,n_r)=-ci*m*dpsidt_Rloc(n_m,n_r)+ &
                  &         three*(r_cmb**2-three*r(n_r)**2)*dVsOm_Rloc(n_m,n_r)
                  dVsOm_Rloc(n_m, n_r)=ro2ms2*r(n_r)*dVsOm_Rloc(n_m,n_r)
               end if
            end if
         end do

      end do

   end subroutine radial_loop_integ
!------------------------------------------------------------------------------
end module rloop_integ
