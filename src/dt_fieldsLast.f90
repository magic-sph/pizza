module fieldsLast
   !
   ! This module contains time-derivatives array of the previous time-step
   ! They are needed in the time-stepping scheme.
   !
 
   use precision_mod
   use constants, only: zero
   use mem_alloc, only: bytes_allocated
   use namelists, only: l_heat, l_chem, l_heat_3D, l_mag_3D, l_mag_LF
   use time_array

   implicit none

   private

   complex(cp), public, allocatable :: dVsT_Mloc(:,:)
   complex(cp), public, allocatable :: dVsXi_Mloc(:,:)
   complex(cp), public, allocatable :: dVsOm_Mloc(:,:)
   complex(cp), public, allocatable :: djxB_Mloc(:,:)
   complex(cp), public, allocatable :: buo_Mloc(:,:)
   complex(cp), public, allocatable :: dpsidt_Rloc(:,:)
   complex(cp), public, allocatable :: dtempdt_Rloc(:,:)
   complex(cp), public, allocatable :: dxidt_Rloc(:,:)
   complex(cp), public, allocatable :: dVsT_Rloc(:,:)
   complex(cp), public, allocatable :: dVsXi_Rloc(:,:)
   complex(cp), public, allocatable :: dVsOm_Rloc(:,:)
   complex(cp), public, allocatable :: buo_Rloc(:,:)
   complex(cp), public, allocatable :: djxB_Rloc(:,:)

   complex(cp), public, allocatable :: dVrT_3D_LMloc(:,:)
   complex(cp), public, allocatable :: dtempdt_3D_Rloc(:,:)
   complex(cp), public, allocatable :: dVrT_3D_Rloc(:,:)
   complex(cp), public, allocatable :: dVxBh_3D_LMloc(:,:)
   complex(cp), public, allocatable :: djdt_3D_Rloc(:,:)
   complex(cp), public, allocatable :: dbdt_3D_Rloc(:,:)
   complex(cp), public, allocatable :: dVxBh_3D_Rloc(:,:)

   type(type_tarray), public :: dpsidt
   type(type_tarray), public :: dTdt, dXidt
   type(type_tarray), public :: dTdt_3D
   type(type_tarray), public :: dBdt_3D
   type(type_tarray), public :: djdt_3D

   public :: initialize_fieldsLast, finalize_fieldsLast, initialize_fieldsLast_3D

contains

   subroutine initialize_fieldsLast(nMstart,nMstop,n_m_max,nRstart,nRstop, &
              &                     n_r_max,norder_imp, norder_exp,        &
              &                     norder_imp_lin)

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

      !-- Local variables
      integer :: n_mloc_fields, n_rloc_fields


      call dpsidt%initialize(nMstart, nMstop, n_r_max, norder_imp, norder_exp, &
           &                 norder_imp_lin)
      if ( l_heat ) then
         call dTdt%initialize(nMstart, nMstop, n_r_max, norder_imp, norder_exp, &
              &               norder_imp_lin)
      end if
      if ( l_chem ) then
         call dxidt%initialize(nMstart, nMstop, n_r_max, norder_imp, norder_exp, &
              &                norder_imp_lin)
      end if

      allocate( buo_Mloc(nMStart:nMstop,n_r_max) )
      n_mloc_fields = 1
      if ( l_heat ) then
         allocate( dVsT_Mloc(nMStart:nMstop,n_r_max) )
         n_mloc_fields = n_mloc_fields + 1
      else
         allocate( dVsT_Mloc(0,0) )
      end if
      if ( l_chem ) then
         allocate( dVsXi_Mloc(nMStart:nMstop,n_r_max) )
         n_mloc_fields = n_mloc_fields + 1
      else
         allocate( dVsXi_Mloc(0,0) )
      end if
      if ( l_mag_LF ) then
         allocate( djxB_Mloc(nMStart:nMstop,n_r_max) )
         n_mloc_fields = n_mloc_fields + 1
      else
         allocate( djxB_Mloc(0,0) )
      end if
      allocate( dVsOm_Mloc(nMStart:nMstop,n_r_max) )
      n_mloc_fields = n_mloc_fields + 1
      bytes_allocated = bytes_allocated + n_mloc_fields*(nMstop-nMStart+1)* &
      &                 n_r_max*SIZEOF_DEF_COMPLEX

      buo_Mloc(:,:)  =zero
      dVsOm_Mloc(:,:)=zero
      if ( l_heat ) dVsT_Mloc(:,:)=zero
      if ( l_chem ) dVsXi_Mloc(:,:)=zero
      if ( l_mag_LF ) djxB_Mloc(:,:)=zero

      allocate( dpsidt_Rloc(n_m_max,nRstart:nRstop) )
      allocate( dVsOm_Rloc(n_m_max,nRstart:nRstop) )
      n_rloc_fields = 2
      if ( l_heat ) then
         allocate( dtempdt_Rloc(n_m_max,nRstart:nRstop) )
         allocate( dVsT_Rloc(n_m_max,nRstart:nRstop) )
         n_rloc_fields = n_rloc_fields+2
      else
         allocate( dtempdt_Rloc(0,0), dVsT_Rloc(0,0) )
      end if
      if ( l_chem ) then
         allocate( dxidt_Rloc(n_m_max,nRstart:nRstop) )
         allocate( dVsXi_Rloc(n_m_max,nRstart:nRstop) )
         n_rloc_fields = n_rloc_fields+2
      else
         allocate( dxidt_Rloc(0,0), dVsXi_Rloc(0,0) )
      end if
      if ( l_heat_3D ) then
         allocate( buo_Rloc(n_m_max,nRstart:nRstop) )
         n_rloc_fields = n_rloc_fields+1
      else
         allocate( buo_Rloc(0,0) )
      end if
      if ( l_mag_LF ) then
         allocate( djxB_Rloc(n_m_max,nRstart:nRstop) )
         n_rloc_fields = n_rloc_fields+1
      else
         allocate( djxB_Rloc(0,0) )
      end if
      bytes_allocated = bytes_allocated + n_rloc_fields*(nRstop-nRStart+1)* &
      &                 n_m_max*SIZEOF_DEF_COMPLEX

      dpsidt_Rloc(:,:)=zero
      dVsOm_Rloc(:,:) =zero
      if ( l_heat ) then
         dVsT_Rloc(:,:)   =zero
         dtempdt_Rloc(:,:)=zero
      end if
      if ( l_chem ) then
         dVsXi_Rloc(:,:)=zero
         dxidt_Rloc(:,:)=zero
      end if
      if ( l_heat_3D ) buo_Rloc(:,:)=zero
      if ( l_mag_LF ) djxB_Rloc(:,:)=zero

   end subroutine initialize_fieldsLast
!-------------------------------------------------------------------------------
   subroutine initialize_fieldsLast_3D(lmStart,lmStop,lm_max,nRstart,nRstop,   &
              &                        n_r_max,norder_imp, norder_exp,  &
              &                        norder_imp_lin)

      !-- Input variables
      integer, intent(in) :: lmStart
      integer, intent(in) :: lmStop
      integer, intent(in) :: lm_max
      integer, intent(in) :: nRstart
      integer, intent(in) :: nRstop
      integer, intent(in) :: n_r_max
      integer, intent(in) :: norder_imp
      integer, intent(in) :: norder_exp
      integer, intent(in) :: norder_imp_lin

      if ( l_heat_3D ) then
         call dTdt_3D%initialize(lmStart, lmStop, n_r_max, norder_imp, norder_exp, &
              &                 norder_imp_lin)
         allocate( dVrT_3D_LMloc(lmStart:lmStop,n_r_max) )
         bytes_allocated = bytes_allocated + 1*(lmStop-lmStart+1)* &
         &                 n_r_max*SIZEOF_DEF_COMPLEX
         dVrT_3D_LMloc(:,:)=zero
      else
         allocate( dVrT_3D_LMloc(0,0) )
      end if

      if ( l_mag_3D ) then
         call dBdt_3D%initialize(lmStart, lmStop, n_r_max, norder_imp, norder_exp, &
              &                 norder_imp_lin)
         call djdt_3D%initialize(lmStart, lmStop, n_r_max, norder_imp, norder_exp, &
              &                 norder_imp_lin)
         allocate( dVxBh_3D_LMloc(lmStart:lmStop,n_r_max) )
         bytes_allocated = bytes_allocated + 1*(lmStop-lmStart+1)* &
         &                 n_r_max*SIZEOF_DEF_COMPLEX
         dVxBh_3D_LMloc(:,:)=zero
      else
         allocate( dVxBh_3D_LMloc(0,0) )
      end if

      if ( l_heat_3D ) then
         allocate( dtempdt_3D_Rloc(lm_max,nRstart:nRstop) )
         allocate( dVrT_3D_Rloc(lm_max,nRstart:nRstop) )
         bytes_allocated = bytes_allocated + 2*(nRstop-nRStart+1)* &
         &                 lm_max*SIZEOF_DEF_COMPLEX
         dVrT_3D_Rloc(:,:)   =zero
         dtempdt_3D_Rloc(:,:)=zero
      else
         allocate( dtempdt_3D_Rloc(0,0), dVrT_3D_Rloc(0,0) )
      end if

      if ( l_mag_3D ) then
         allocate( dbdt_3D_Rloc(lm_max,nRstart:nRstop) )
         allocate( djdt_3D_Rloc(lm_max,nRstart:nRstop) )
         allocate( dVxBh_3D_Rloc(lm_max,nRstart:nRstop) )
         bytes_allocated = bytes_allocated + 3*(nRstop-nRStart+1)* &
         &                 lm_max*SIZEOF_DEF_COMPLEX
         dVxBh_3D_Rloc(:,:)=zero
         djdt_3D_Rloc(:,:) =zero
         dbdt_3D_Rloc(:,:) =zero
      else
         allocate( dbdt_3D_Rloc(0,0), djdt_3D_Rloc(0,0), dVxBh_3D_Rloc(0,0) )
      end if

   end subroutine initialize_fieldsLast_3D
!-------------------------------------------------------------------------------
   subroutine finalize_fieldsLast

      if ( l_mag_3D ) then
         call dBdt_3D%finalize()
         call djdt_3D%finalize()
         deallocate( dVxBh_3D_LMloc, dVxBh_3D_Rloc, djdt_3D_Rloc, dbdt_3D_Rloc )
         if ( l_mag_LF ) then
            deallocate( djxB_Mloc, djxB_Rloc )
         end if
      end if
      if ( l_heat_3D ) then
         call dTdt_3D%finalize()
         deallocate( buo_Rloc )
         deallocate( dVrT_3D_LMloc, dVrT_3D_Rloc, dtempdt_3D_Rloc )
      end if
      if ( l_heat ) then
         call dTdt%finalize()
         deallocate( dVsT_Mloc, dVsT_Rloc, dtempdt_Rloc )
      end if
      if ( l_chem ) then
         call dxidt%finalize()
         deallocate( dVsXi_Mloc, dVsXi_Rloc, dxidt_Rloc )
      end if
      call dpsidt%finalize()
      deallocate( dVsOm_Rloc, dVsOm_Mloc, buo_Mloc, dpsidt_Rloc )

   end subroutine finalize_fieldsLast
!-------------------------------------------------------------------------------
end module fieldsLast
