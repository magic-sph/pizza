module fields
   !
   !  This module contains the potential fields and their radial
   !  derivatives
   !
   use precision_mod
   use constants, only: zero
   use namelists, only: l_cheb_coll, l_heat, l_chem, l_heat_3D
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_m_max, n_r_max
   use truncation_3D, only: lm_max, n_r_max_3D, n_theta_max, n_phi_max_3D
   use blocking, only: nMstart, nMstop, nRstart, nRstop, nRstart3D, nRstop3D
   use blocking_lm, only: llm, ulm

 
   implicit none

   private
 
   !-- Velocity potentials:
   complex(cp), allocatable, public :: temp_Mloc(:,:), psi_Mloc(:,:)
   complex(cp), allocatable, public :: temp_Rloc(:,:), dtemp_Mloc(:,:)
   complex(cp), allocatable, public :: xi_Mloc(:,:), dxi_Mloc(:,:)
   complex(cp), allocatable, public :: om_Rloc(:,:), om_Mloc(:,:)
   complex(cp), allocatable, public :: dom_Mloc(:,:), xi_Rloc(:,:)
   complex(cp), allocatable, public :: work_Mloc(:,:)

   complex(cp), allocatable, public :: us_Mloc(:,:), up_Mloc(:,:)
   complex(cp), allocatable, public :: us_Rloc(:,:), up_Rloc(:,:)

   complex(cp), allocatable, public :: psi_hat_Mloc(:,:)
   complex(cp), allocatable, public :: temp_hat_Mloc(:,:)
   complex(cp), allocatable, public :: xi_hat_Mloc(:,:)

   !-- 3-D Temperature:
   complex(cp), allocatable, public :: temp_3D_LMloc(:,:)
   complex(cp), allocatable, public :: dtemp_3D_LMloc(:,:)
   complex(cp), allocatable, public :: work_LMloc(:,:)
   complex(cp), allocatable, public :: temp_3D_Rloc(:,:)
   complex(cp), allocatable, public :: work_3D_Rloc(:,:)

   !-- 3-D velocity in the physical space
   real(cp), allocatable, public :: ur_3D_Rloc(:,:,:)
   real(cp), allocatable, public :: ut_3D_Rloc(:,:,:)
   real(cp), allocatable, public :: up_3D_Rloc(:,:,:)

 
   public :: initialize_fields, finalize_fields

contains

   subroutine initialize_fields

      integer :: n_mloc_fields, n_rloc_fields

      n_mloc_fields = 0
      if ( l_heat ) then
         allocate( temp_Mloc(nMStart:nMstop,n_r_max) )
         allocate( dtemp_Mloc(nMStart:nMstop,n_r_max) )
         n_mloc_fields = n_mloc_fields + 2
      else
         allocate( temp_Mloc(0,0), dtemp_Mloc(0,0) )
      end if
      if ( l_chem ) then
         allocate( xi_Mloc(nMStart:nMstop,n_r_max) )
         allocate( dxi_Mloc(nMStart:nMstop,n_r_max) )
         n_mloc_fields = n_mloc_fields + 2
      else
         allocate(xi_Mloc(0,0), dxi_Mloc(0,0) )
      end if
      allocate( psi_Mloc(nMStart:nMstop,n_r_max) )
      allocate( om_Mloc(nMStart:nMstop,n_r_max) )
      allocate( dom_Mloc(nMStart:nMstop,n_r_max) )
      allocate( us_Mloc(nMStart:nMstop,n_r_max) )
      allocate( up_Mloc(nMStart:nMstop,n_r_max) )
      allocate( work_Mloc(nMStart:nMstop,n_r_max) )
      n_mloc_fields = n_mloc_fields + 6
      bytes_allocated = bytes_allocated + n_mloc_fields* &
      &                 (nMstop-nMstart+1)*n_r_max*SIZEOF_DEF_COMPLEX

      if ( l_heat ) then
         temp_Mloc(:,:) =zero
         dtemp_Mloc(:,:)=zero
      end if
      if ( l_chem ) then
         xi_Mloc(:,:) =zero
         dxi_Mloc(:,:)=zero
      end if
      psi_Mloc(:,:) =zero
      om_Mloc(:,:)  =zero
      dom_Mloc(:,:) =zero
      us_Mloc(:,:)  =zero
      up_Mloc(:,:)  =zero
      work_Mloc(:,:)=zero

      if ( .not. l_cheb_coll ) then
         allocate( psi_hat_Mloc(nMstart:nMstop,n_r_max) )
         n_mloc_fields = 1
         if ( l_heat ) then
            allocate( temp_hat_Mloc(nMstart:nMstop,n_r_max) )
            n_mloc_fields = n_mloc_fields + 1
         else
            allocate( temp_hat_Mloc(0,0) )
         end if
         if ( l_chem ) then
            allocate( xi_hat_Mloc(nMstart:nMstop,n_r_max) )
            n_mloc_fields = n_mloc_fields + 1
         else
            allocate( xi_hat_Mloc(0,0) )
         end if
         bytes_allocated = bytes_allocated+n_mloc_fields* &
         &                 (nMstop-nMstart+1)*n_r_max*SIZEOF_DEF_COMPLEX
         psi_hat_Mloc(:,:) =zero
         if ( l_heat ) temp_hat_Mloc(:,:)=zero
         if ( l_chem ) xi_hat_Mloc(:,:)=zero
      else
         ! We have to allocate the arrays to make the code runs when
         ! debug flags are turned on
         allocate (psi_hat_Mloc(0,0), temp_hat_Mloc(0,0), xi_hat_Mloc(0,0) )
      end if


      n_rloc_fields = 0
      if ( l_heat ) then
         allocate( temp_Rloc(n_m_max,nRstart:nRstop) )
         n_rloc_fields = n_rloc_fields+1
      else
         allocate( temp_Rloc(0,0) )
      end if
      if ( l_chem ) then
         allocate( xi_Rloc(n_m_max,nRstart:nRstop) )
         n_rloc_fields = n_rloc_fields+1
      else
         allocate( xi_Rloc(0,0) )
      end if
      allocate( us_Rloc(n_m_max,nRstart:nRstop) )
      allocate( up_Rloc(n_m_max,nRstart:nRstop) )
      allocate( om_Rloc(n_m_max,nRstart:nRstop) )
      n_rloc_fields = n_rloc_fields+3
      bytes_allocated = bytes_allocated + n_rloc_fields* &
      &                 n_m_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

      if ( l_heat ) temp_Rloc(:,:)=zero
      if ( l_chem ) xi_Rloc(:,:)=zero
      us_Rloc(:,:)=zero
      up_Rloc(:,:)=zero
      om_Rloc(:,:)=zero

      !-- 3-D temperature:
      if ( l_heat_3D ) then
         allocate( temp_3D_LMloc(llm:ulm,n_r_max_3D) )
         allocate( dtemp_3D_LMloc(llm:ulm,n_r_max_3D) )
         allocate( work_LMloc(llm:ulm,n_r_max_3D) )
         allocate( temp_3D_Rloc(lm_max,nRstart3D:nRstop3D) )
         allocate( work_3D_Rloc(lm_max,nRstart3D:nRstop3D) )
         bytes_allocated = bytes_allocated + &
         &                 2*(ulm-llm+1)*n_r_max_3D*SIZEOF_DEF_COMPLEX
         bytes_allocated = bytes_allocated + &
         &                 2*lm_max*(nRstop3D-nRstart3D+1)*SIZEOF_DEF_COMPLEX
         temp_3D_LMloc(:,:) =zero
         dtemp_3D_LMloc(:,:)=zero
         work_LMloc(:,:)    =zero
         temp_3D_Rloc(:,:)  =zero
         work_3D_Rloc(:,:)  =zero

         allocate( ur_3D_Rloc(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D) )
         allocate( ut_3D_Rloc(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D) )
         allocate( up_3D_Rloc(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D) )
         bytes_allocated = bytes_allocated + 3*n_phi_max_3D*n_theta_max*&
         &                 (nRstop3D-nRstart3D+1)*SIZEOF_DEF_REAL
         ur_3D_Rloc(:,:,:)=0.0_cp
         ut_3D_Rloc(:,:,:)=0.0_cp
         up_3D_Rloc(:,:,:)=0.0_cp
      else
         allocate( temp_3D_Rloc(1,1), temp_3D_LMloc(1,1) )
      end if

   end subroutine initialize_fields
!----------------------------------------------------------------------------
   subroutine finalize_fields

      if ( l_heat_3D ) then
         deallocate( ur_3D_Rloc, ut_3D_Rloc, up_3D_Rloc )
         deallocate( work_LMloc, temp_3D_LMloc, dtemp_3D_LMloc )
         deallocate( work_3D_Rloc, temp_3D_Rloc )
      end if
      if ( l_heat ) then
         deallocate( temp_Mloc, dtemp_Mloc, temp_Rloc )
         if ( .not. l_cheb_coll ) deallocate( temp_hat_Mloc )
      end if
      if ( l_chem ) then
         deallocate( xi_Mloc, dxi_Mloc, xi_Rloc )
         if ( .not. l_cheb_coll ) deallocate( xi_hat_Mloc )
      end if
      deallocate( us_Rloc, up_Rloc, us_Mloc, up_Mloc, work_Mloc )
      deallocate( om_Rloc, om_Mloc, dom_Mloc )
      if ( .not. l_cheb_coll ) deallocate( psi_hat_Mloc )

   end subroutine finalize_fields
!----------------------------------------------------------------------------
end module fields
