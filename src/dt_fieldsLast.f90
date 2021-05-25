module fieldsLast
   !
   ! This module contains time-derivatives array of the previous time-step
   ! They are needed in the time-stepping scheme.
   !
 
   use precision_mod
   use constants, only: zero
   use mem_alloc, only: bytes_allocated
   use namelists, only: l_heat, l_chem
   use time_array

   implicit none

   private

   complex(cp), allocatable, target, public :: dt_fields_Rloc_container(:,:,:)
   complex(cp), allocatable, target, public :: dt_fields_Mloc_container(:,:,:,:)
   complex(cp), pointer, public :: dVsT_Mloc(:,:,:)
   complex(cp), pointer, public :: dVsXi_Mloc(:,:,:)
   complex(cp), pointer, public :: dVsOm_Mloc(:,:,:)
   complex(cp), pointer, public :: dpsidt_Rloc(:,:)
   complex(cp), pointer, public :: dtempdt_Rloc(:,:)
   complex(cp), pointer, public :: dxidt_Rloc(:,:)
   complex(cp), pointer, public :: dVsT_Rloc(:,:)
   complex(cp), pointer, public :: dVsXi_Rloc(:,:)
   complex(cp), pointer, public :: dVsOm_Rloc(:,:)
   type(type_tarray), public :: dpsidt, dTdt, dXidt

   public :: initialize_fieldsLast, finalize_fieldsLast

contains

   subroutine initialize_fieldsLast(nMstart,nMstop,n_m_max,nRstart,nRstop,n_r_max,&
              &                     norder_imp, nexp, norder_imp_lin)

      !-- Input variables
      integer, intent(in) :: nMstart
      integer, intent(in) :: nMstop
      integer, intent(in) :: n_m_max
      integer, intent(in) :: nRstart
      integer, intent(in) :: nRstop
      integer, intent(in) :: n_r_max
      integer, intent(in) :: norder_imp
      integer, intent(in) :: nexp
      integer, intent(in) :: norder_imp_lin

      !-- Local variables
      integer :: n_mloc_fields, n_rloc_fields


      call dpsidt%initialize(nMstart, nMstop, n_r_max, norder_imp, nexp, &
           &                 norder_imp_lin)
      if ( l_heat ) then
         call dTdt%initialize(nMstart, nMstop, n_r_max, norder_imp, nexp, &
              &               norder_imp_lin)
      end if
      if ( l_chem ) then
         call dxidt%initialize(nMstart, nMstop, n_r_max, norder_imp, nexp, &
              &                norder_imp_lin)
      end if

      n_mloc_fields = 0
      if ( l_heat ) then
         n_mloc_fields = n_mloc_fields + 1
      else
         allocate( dVsT_Mloc(0,0,0) )
      end if
      if ( l_chem ) then
         n_mloc_fields = n_mloc_fields + 1
      else
         allocate( dVsXi_Mloc(0,0,0) )
      end if
      allocate(dt_fields_Mloc_container(nMStart:nMstop,n_r_max,2+2*n_mloc_fields,1:nexp))
      dpsidt%expl(nMstart:,1:,1:) => dt_fields_Mloc_container(nMstart:nMstop,1:n_r_max,1,&
      &                                                       1:nexp)
      dVsOm_Mloc(nMstart:,1:,1:) => dt_fields_Mloc_container(nMstart:nMstop,1:n_r_max,2, &
      &                                                      1:nexp)
      if ( l_heat ) then
         dTdt%expl(nMstart:,1:,1:) => dt_fields_Mloc_container(nMstart:nMstop,1:n_r_max, &
         &                                                     3,1:nexp)
         dVsT_Mloc(nMstart:,1:,1:) => dt_fields_Mloc_container(nMstart:nMstop,1:n_r_max, &
         &                                                     4,1:nexp)
         if ( l_chem ) then
            dxidt%expl(nMstart:,1:,1:) => dt_fields_Mloc_container(nMstart:nMstop, &
            &                                                     1:n_r_max,5,1:nexp)
            dVsXi_Mloc(nMstart:,1:,1:) => dt_fields_Mloc_container(nMstart:nMstop, &
            &                                                      1:n_r_max,6,1:nexp)
         end if
      else
         if ( l_chem ) then
            dxidt%expl(nMstart:,1:,1:) => dt_fields_Mloc_container(nMstart:nMstop, &
            &                                                     1:n_r_max,3,1:nexp)
            dVsXi_Mloc(nMstart:,1:,1:) => dt_fields_Mloc_container(nMstart:nMstop, &
            &                                                      1:n_r_max,4,1:nexp)
         end if
      end if

      bytes_allocated = bytes_allocated + (2*n_mloc_fields+2)*nexp*(nMstop-nMStart+1)* &
      &                 n_r_max*SIZEOF_DEF_COMPLEX

      dt_fields_Mloc_container(:,:,:,:)=zero

      n_rloc_fields = 2
      if ( l_heat ) then
         n_rloc_fields = n_rloc_fields+2
      else
         allocate( dtempdt_Rloc(0,0), dVsT_Rloc(0,0) )
      end if
      if ( l_chem ) then
         n_rloc_fields = n_rloc_fields+2
      else
         allocate( dxidt_Rloc(0,0), dVsXi_Rloc(0,0) )
      end if
      allocate( dt_fields_Rloc_container(n_m_max,nRstart:nRstop,n_rloc_fields) )
      dpsidt_Rloc(1:,nRstart:) => dt_fields_Rloc_container(1:n_m_max,nRstart:nRstop,1)
      dVsOm_Rloc(1:,nRstart:) => dt_fields_Rloc_container(1:n_m_max,nRstart:nRstop,2)
      if ( l_heat ) then
         dtempdt_Rloc(1:,nRstart:) => dt_fields_Rloc_container(1:n_m_max,nRstart:nRstop,3)
         dVsT_Rloc(1:,nRstart:) => dt_fields_Rloc_container(1:n_m_max,nRstart:nRstop,4)
         if ( l_chem ) then
            dxidt_Rloc(1:,nRstart:) => dt_fields_Rloc_container(1:n_m_max,nRstart:nRstop,5)
            dVsXi_Rloc(1:,nRstart:) => dt_fields_Rloc_container(1:n_m_max,nRstart:nRstop,6)
         end if
      else
         if ( l_chem ) then
            dxidt_Rloc(1:,nRstart:) => dt_fields_Rloc_container(1:n_m_max,nRstart:nRstop,3)
            dVsXi_Rloc(1:,nRstart:) => dt_fields_Rloc_container(1:n_m_max,nRstart:nRstop,4)
         end if
      end if
      bytes_allocated = bytes_allocated + n_rloc_fields*(nRstop-nRStart+1)* &
      &                 n_m_max*SIZEOF_DEF_COMPLEX

      dt_fields_Rloc_container(:,:,:)=zero

   end subroutine initialize_fieldsLast
!-------------------------------------------------------------------------------
   subroutine finalize_fieldsLast

      deallocate( dt_fields_Rloc_container, dt_fields_Mloc_container )
      if ( l_heat ) call dTdt%finalize()
      if ( l_chem ) call dxidt%finalize()
      call dpsidt%finalize()

   end subroutine finalize_fieldsLast
!-------------------------------------------------------------------------------
end module fieldsLast
