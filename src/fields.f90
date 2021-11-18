module fields
   !
   !  This module contains the potential fields and their radial
   !  derivatives
   !
   use precision_mod
   use constants, only: zero
   use namelists, only: l_cheb_coll, l_heat, l_chem,  l_finite_diff
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_m_max, n_r_max
   use blocking, only: nMstart, nMstop, nRstart, nRstop
 
   implicit none

   private
 
   !-- Velocity potentials:
   complex(cp), allocatable, public :: psi_Mloc(:,:), dtemp_Mloc(:,:)
   complex(cp), allocatable, public :: dxi_Mloc(:,:), dom_Mloc(:,:)
   complex(cp), allocatable, public :: work_Mloc(:,:)
   complex(cp), allocatable, public :: dtemp_Rloc(:,:), dxi_Rloc(:,:)
   complex(cp), pointer, public :: temp_Mloc(:,:), om_Mloc(:,:), xi_Mloc(:,:)
   complex(cp), pointer, public :: us_Mloc(:,:), up_Mloc(:,:)
   complex(cp), allocatable, target, public :: fields_Mloc_container(:,:,:)

   complex(cp), pointer, public :: us_Rloc(:,:), up_Rloc(:,:)
   complex(cp), pointer, public :: om_Rloc(:,:), temp_Rloc(:,:), xi_Rloc(:,:)
   complex(cp), allocatable, target, public :: fields_Rloc_container(:,:,:)

   complex(cp), allocatable, public :: psi_hat_Mloc(:,:)
   complex(cp), allocatable, public :: temp_hat_Mloc(:,:)
   complex(cp), allocatable, public :: xi_hat_Mloc(:,:)
 
   public :: initialize_fields, finalize_fields

contains

   subroutine initialize_fields

      integer :: n_mloc_fields, n_rloc_fields

      n_mloc_fields = 0
      if ( l_heat ) then
         allocate( dtemp_Mloc(nMStart:nMstop,n_r_max) )
         n_mloc_fields = n_mloc_fields + 1
      else
         allocate( temp_Mloc(0,0), dtemp_Mloc(0,0) )
      end if
      if ( l_chem ) then
         allocate( dxi_Mloc(nMStart:nMstop,n_r_max) )
         n_mloc_fields = n_mloc_fields + 1
      else
         allocate( xi_Mloc(0,0), dxi_Mloc(0,0) )
      end if

      allocate( fields_Mloc_container(nMstart:nMstop, n_r_max, 3+n_mloc_fields) )
      us_Mloc(nMstart:,1:) => fields_Mloc_container(nMstart:nMstop,1:n_r_max,1)
      up_Mloc(nMstart:,1:) => fields_Mloc_container(nMstart:nMstop,1:n_r_max,2)
      om_Mloc(nMstart:,1:) => fields_Mloc_container(nMstart:nMstop,1:n_r_max,3)

      if ( l_heat ) then
         temp_Mloc(nMstart:,1:) => fields_Mloc_container(nMstart:nMstop,1:n_r_max,4)
         if ( l_chem ) then ! Heat and Chem
            xi_Mloc(nMstart:,1:) => fields_Mloc_container(nMstart:nMstop,1:n_r_max,5)
         end if
      else
         if ( l_chem ) then
            xi_Mloc(nMstart:,1:) => fields_Mloc_container(nMstart:nMstop,1:n_r_max,4)
         end if
      end if
      allocate( psi_Mloc(nMStart:nMstop,n_r_max) )
      allocate( dom_Mloc(nMStart:nMstop,n_r_max) )
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
         psi_hat_Mloc(:,:)=zero
         if ( l_heat ) temp_hat_Mloc(:,:)=zero
         if ( l_chem ) xi_hat_Mloc(:,:)=zero
      else
         ! We have to allocate the arrays to make the code runs when
         ! debug flags are turned on
         allocate (psi_hat_Mloc(0,0), temp_hat_Mloc(0,0), xi_hat_Mloc(0,0) )
      end if

      n_rloc_fields = 0
      if ( l_heat ) then
         n_rloc_fields = n_rloc_fields+1
      else
         allocate( temp_Rloc(0,0) )
      end if
      if ( l_chem ) then
         n_rloc_fields = n_rloc_fields+1
      else
         allocate( xi_Rloc(0,0) )
      end if
      allocate( fields_Rloc_container(n_m_max, nRstart:nRstop, 3+n_rloc_fields) )
      us_Rloc(1:,nRstart:) => fields_Rloc_container(1:n_m_max,nRstart:nRstop,1)
      up_Rloc(1:,nRstart:) => fields_Rloc_container(1:n_m_max,nRstart:nRstop,2)
      om_Rloc(1:,nRstart:) => fields_Rloc_container(1:n_m_max,nRstart:nRstop,3)
      if ( l_heat ) then
         temp_Rloc(1:,nRstart:) => fields_Rloc_container(1:n_m_max,nRstart:nRstop,4)
         if ( l_chem ) then
            xi_Rloc(1:,nRstart:) => fields_Rloc_container(1:n_m_max,nRstart:nRstop,5)
         end if
      else
         if ( l_chem ) then
            xi_Rloc(1:,nRstart:) => fields_Rloc_container(1:n_m_max,nRstart:nRstop,4)
         end if
      end if

      n_rloc_fields = n_rloc_fields+3
      bytes_allocated = bytes_allocated + n_rloc_fields* &
      &                 n_m_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

      fields_Rloc_container(:,:,:)=zero

      if ( l_finite_diff ) then
         n_rloc_fields=0
         if ( l_heat ) then
            n_rloc_fields=n_rloc_fields+1
            allocate( dtemp_Rloc(n_m_max,nRstart:nRstop) )
         else
            allocate( dtemp_Rloc(1,1) )
         end if
         if ( l_chem ) then
            n_rloc_fields=n_rloc_fields+1
            allocate( dxi_Rloc(n_m_max,nRstart:nRstop) )
         else
            allocate( dxi_Rloc(1,1) )
         end if
         bytes_allocated = bytes_allocated + n_rloc_fields* &
         &                 n_m_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
      else
         allocate( dtemp_Rloc(1,1) )
         allocate( dxi_Rloc(1,1) )
      end if

   end subroutine initialize_fields
!----------------------------------------------------------------------------
   subroutine finalize_fields

      deallocate( work_Mloc, dom_Mloc, psi_Mloc )
      if ( .not. l_cheb_coll ) deallocate( psi_hat_Mloc )
      if ( l_chem ) then
         deallocate( dxi_Mloc )
         if ( .not. l_cheb_coll ) deallocate( xi_hat_Mloc )
      end if
      if ( l_heat ) then
         deallocate( dtemp_Mloc )
         if ( .not. l_cheb_coll ) deallocate( temp_hat_Mloc )
      end if

      deallocate(fields_Rloc_container, fields_Mloc_container)
      deallocate( dtemp_Rloc, dxi_Rloc )

   end subroutine finalize_fields
!----------------------------------------------------------------------------
end module fields
