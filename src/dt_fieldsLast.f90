module fieldsLast
   !
   ! This module contains time-derivatives array of the previous time-step
   ! They are needed in the time-stepping scheme.
   !
 
   use precision_mod
   use constants, only: zero
   use mem_alloc, only: bytes_allocated
   use time_array

   implicit none

   private

   complex(cp), public, allocatable :: dVsT_Mloc(:,:)
   complex(cp), public, allocatable :: dVsOm_Mloc(:,:)
   complex(cp), public, allocatable :: buo_imp_Mloc(:,:)
   complex(cp), public, allocatable :: dpsidt_Rloc(:,:)
   complex(cp), public, allocatable :: dtempdt_Rloc(:,:)
   complex(cp), public, allocatable :: dVsT_Rloc(:,:)
   complex(cp), public, allocatable :: dVsOm_Rloc(:,:)
   type(type_tarray), public :: dpsidt
   type(type_tarray), public :: dTdt

   public :: initialize_fieldsLast, finalize_fieldsLast

contains

   subroutine initialize_fieldsLast(nMstart,nMstop,n_m_max,nRstart,nRstop,n_r_max,&
              &                     norder_imp, norder_exp, norder_imp_lin)

      !-- Input variables
      integer, intent(in) :: nMstart
      integer, intent(in) :: nMstop
      integer, intent(in) :: n_m_max
      integer, intent(in) :: nRstart
      integer, intent(in) :: nRstop
      integer, intent(in) :: n_r_max
      integer, intent(in) :: norder_imp
      integer, intent(in) :: norder_exp
      integer, intent(in) :: norder_imp_lin

      call dpsidt%initialize(nMstart, nMstop, n_r_max, norder_imp, norder_exp, &
           &                 norder_imp_lin)
      call dTdt%initialize(nMstart, nMstop, n_r_max, norder_imp, norder_exp, &
           &               norder_imp_lin)

      allocate( buo_imp_Mloc(nMStart:nMstop,n_r_max) )
      allocate( dVsT_Mloc(nMStart:nMstop,n_r_max) )
      allocate( dVsOm_Mloc(nMStart:nMstop,n_r_max) )
      bytes_allocated = bytes_allocated + 3*(nMstop-nMStart+1)*n_r_max*&
      &                 SIZEOF_DEF_COMPLEX

      buo_imp_Mloc(:,:)=zero
      dVsT_Mloc(:,:)   =zero
      dVsOm_Mloc(:,:)  =zero

      allocate( dpsidt_Rloc(n_m_max,nRstart:nRstop) )
      allocate( dtempdt_Rloc(n_m_max,nRstart:nRstop) )
      allocate( dVsT_Rloc(n_m_max,nRstart:nRstop) )
      allocate( dVsOm_Rloc(n_m_max,nRstart:nRstop) )
      bytes_allocated = bytes_allocated + 4*(nRstop-nRStart+1)*n_m_max* &
      &                 SIZEOF_DEF_COMPLEX

      dpsidt_Rloc(:,:) =zero
      dtempdt_Rloc(:,:)=zero
      dVsT_Rloc(:,:)   =zero
      dVsOm_Rloc(:,:)  =zero

   end subroutine initialize_fieldsLast
!-------------------------------------------------------------------------------
   subroutine finalize_fieldsLast

      call dTdt%finalize()
      call dpsidt%finalize()
      deallocate( dVsOm_Rloc, dVsOm_Mloc, buo_imp_Mloc )
      deallocate( dVsT_Rloc, dVsT_Mloc )
      deallocate( dpsidt_Rloc, dtempdt_Rloc )

   end subroutine finalize_fieldsLast
!-------------------------------------------------------------------------------
end module fieldsLast
