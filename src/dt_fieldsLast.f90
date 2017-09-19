module fieldsLast
   !
   ! This module contains time-derivatives array of the previous time-step
   ! They are needed in the time-stepping scheme.
   !
   ! The variables labeled with a suffix 'Last' are provided
   ! by the restart file for the first time step or
   ! calculated here or by the update routines for the
   ! following time step.
   ! These fields remain in the M-distributed space 
 
   use precision_mod
   use constants, only: zero
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, n_m_max
   use blocking, only: nMStart, nMstop, nRstart, nRstop

   implicit none

   private

   complex(cp), public, allocatable :: dVsT_Mloc(:,:)
   complex(cp), public, allocatable :: dVsOm_Mloc(:,:)
   complex(cp), public, allocatable :: dtemp_imp_Mloc(:,:,:)
   complex(cp), public, allocatable :: dtemp_exp_Mloc(:,:,:)
   complex(cp), public, allocatable :: dpsi_imp_Mloc(:,:,:)
   complex(cp), public, allocatable :: buo_imp_Mloc(:,:)
   complex(cp), public, allocatable :: dpsi_exp_Mloc(:,:,:)
   complex(cp), public, allocatable :: dpsidt_Rloc(:,:)
   complex(cp), public, allocatable :: dtempdt_Rloc(:,:)
   complex(cp), public, allocatable :: dVsT_Rloc(:,:)
   complex(cp), public, allocatable :: dVsOm_Rloc(:,:)

   public :: initialize_fieldsLast, finalize_fieldsLast

contains

   subroutine initialize_fieldsLast(norder_exp,norder_imp)

      integer, intent(in) :: norder_exp
      integer, intent(in) :: norder_imp

      allocate( dtemp_imp_Mloc(nMStart:nMstop,n_r_max,norder_imp-1) )
      allocate( dtemp_exp_Mloc(nMStart:nMstop,n_r_max,norder_exp) )
      allocate( dpsi_imp_Mloc(nMStart:nMstop,n_r_max,norder_imp-1) )
      allocate( buo_imp_Mloc(nMStart:nMstop,n_r_max) )
      allocate( dpsi_exp_Mloc(nMStart:nMstop,n_r_max,norder_exp) )
      allocate( dVsT_Mloc(nMStart:nMstop,n_r_max) )
      allocate( dVsOm_Mloc(nMStart:nMstop,n_r_max) )
      bytes_allocated = bytes_allocated + (3+2*norder_exp+2*(norder_exp-1))*&
      &                 (nMstop-nMStart+1)*n_r_max*SIZEOF_DEF_COMPLEX

      dtemp_imp_Mloc(:,:,:)=zero
      dtemp_exp_Mloc(:,:,:)=zero
      dpsi_imp_Mloc(:,:,:) =zero
      buo_imp_Mloc(:,:)    =zero
      dpsi_exp_Mloc(:,:,:) =zero
      dVsT_Mloc(:,:)       =zero
      dVsOm_Mloc(:,:)      =zero

      allocate( dpsidt_Rloc(n_m_max,nRstart:nRstop) )
      allocate( dtempdt_Rloc(n_m_max,nRstart:nRstop) )
      allocate( dVsT_Rloc(n_m_max,nRstart:nRstop) )
      allocate( dVsOm_Rloc(n_m_max,nRstart:nRstop) )
      bytes_allocated = bytes_allocated + &
      &                 4*(nRstop-nRStart+1)*n_m_max*SIZEOF_DEF_COMPLEX

      dpsidt_Rloc(:,:) =zero
      dtempdt_Rloc(:,:)=zero
      dVsT_Rloc(:,:)   =zero
      dVsOm_Rloc(:,:)  =zero

   end subroutine initialize_fieldsLast
!-------------------------------------------------------------------------------
   subroutine finalize_fieldsLast

      deallocate( dVsOm_Rloc, dVsOm_Mloc )
      deallocate( dVsT_Rloc, dVsT_Mloc )
      deallocate( dtemp_imp_Mloc, dtemp_exp_Mloc )
      deallocate( dpsi_imp_Mloc, dpsi_exp_Mloc, buo_imp_Mloc )
      deallocate( dpsidt_Rloc, dtempdt_Rloc )

   end subroutine finalize_fieldsLast
!-------------------------------------------------------------------------------
end module fieldsLast
