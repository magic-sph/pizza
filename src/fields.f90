module fields
   !
   !  This module contains the potential fields and their radial
   !  derivatives
   !
   use precision_mod
   use constants, only: zero
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_m_max, n_r_max
   use blocking, only: nMstart, nMstop, nRstart, nRstop
 
   implicit none

   private
 
   !-- Velocity potentials:
   complex(cp), allocatable, public :: temp_Mloc(:,:), psi_Mloc(:,:)
   complex(cp), allocatable, public :: temp_Rloc(:,:), dtemp_Mloc(:,:)
   complex(cp), allocatable, public :: om_Rloc(:,:), om_Mloc(:,:)
   complex(cp), allocatable, public :: dom_Mloc(:,:)
   complex(cp), allocatable, public :: work_Mloc(:,:)

   complex(cp), allocatable, public :: us_Mloc(:,:), up_Mloc(:,:)
   complex(cp), allocatable, public :: us_Rloc(:,:), up_Rloc(:,:)
 
   public :: initialize_fields, finalize_fields

contains

   subroutine initialize_fields

      allocate( temp_Mloc(nMStart:nMstop,n_r_max) )
      allocate( dtemp_Mloc(nMStart:nMstop,n_r_max) )
      allocate( psi_Mloc(nMStart:nMstop,n_r_max) )
      allocate( om_Mloc(nMStart:nMstop,n_r_max) )
      allocate( dom_Mloc(nMStart:nMstop,n_r_max) )
      allocate( us_Mloc(nMStart:nMstop,n_r_max) )
      allocate( up_Mloc(nMStart:nMstop,n_r_max) )
      allocate( work_Mloc(nMStart:nMstop,n_r_max) )
      bytes_allocated = bytes_allocated + &
      &                 8*(nMstop-nMstart+1)*n_r_max*SIZEOF_DEF_COMPLEX

      temp_Mloc(:,:) =zero
      dtemp_Mloc(:,:)=zero
      psi_Mloc(:,:)  =zero
      om_Mloc(:,:)   =zero
      dom_Mloc(:,:)  =zero
      us_Mloc(:,:)   =zero
      up_Mloc(:,:)   =zero
      work_Mloc(:,:) =zero

      allocate( temp_Rloc(n_m_max,nRstart:nRstop) )
      allocate( us_Rloc(n_m_max,nRstart:nRstop) )
      allocate( up_Rloc(n_m_max,nRstart:nRstop) )
      allocate( om_Rloc(n_m_max,nRstart:nRstop) )
      bytes_allocated = bytes_allocated + &
      &                 4*n_m_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

      temp_Rloc(:,:) =zero
      us_Rloc(:,:)   =zero
      up_Rloc(:,:)   =zero
      om_Rloc(:,:)   =zero

   end subroutine initialize_fields
!----------------------------------------------------------------------------
   subroutine finalize_fields

      deallocate( us_Rloc, up_Rloc, us_Mloc, up_Mloc, work_Mloc )
      deallocate( om_Rloc, om_Mloc, dom_Mloc )
      deallocate( temp_Mloc, dtemp_Mloc, psi_Mloc, temp_Rloc )

   end subroutine finalize_fields
!----------------------------------------------------------------------------
end module fields
