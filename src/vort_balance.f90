module vort_balance
   !
   ! This module is used to compute the balance of the terms that enter
   ! the vorticity equations.
   !

   use precision_mod
   use parallel_mod
   use constants, only: zero, vol_otc
   use time_schemes, only: type_tscheme
   use mem_alloc, only: bytes_allocated
   use communications, only: r2m_single
   use truncation, only: n_r_max, m_max, minc, n_m_max, idx2m
   use blocking, only: nMstart, nMstop, nm_per_rank, m_balance, nRstart, nRstop
   use namelists, only: ra, pr, ek, radratio, sc, raxi, tag, r_cmb, r_icb, &
       &                bl_cut, l_2D_SD, l_finite_diff
   use radial_functions, only: r, height
   use integration, only: simps
   use useful, only: cc2real, getMSD2
   use mean_sd, only: mean_sd_type, mean_sd_2D_type

   implicit none

   private

   type, public :: vort_bal_type
      integer :: nRmin, nRmax, mMin, mMax
      complex(cp), allocatable :: visc(:,:)
      complex(cp), allocatable :: cor(:,:)
      complex(cp), allocatable :: adv(:,:)
      complex(cp), allocatable :: pump(:,:)
      complex(cp), allocatable :: dwdt(:,:)
      complex(cp), allocatable :: buo(:,:)
      type(mean_sd_2D_type) :: visc2D
      type(mean_sd_2D_type) :: cor2D
      type(mean_sd_2D_type) :: pump2D
      type(mean_sd_2D_type) :: adv2D
      type(mean_sd_2D_type) :: dwdt2D
      type(mean_sd_2D_type) :: buo2D
      type(mean_sd_2D_type) :: iner2D
      type(mean_sd_2D_type) :: thwind2D
      type(mean_sd_2D_type) :: cia2D
      type(mean_sd_type) :: viscM
      type(mean_sd_type) :: pumpM
      type(mean_sd_type) :: buoM
      type(mean_sd_type) :: corM
      type(mean_sd_type) :: advM
      type(mean_sd_type) :: dwdtM
      type(mean_sd_type) :: inerM
      type(mean_sd_type) :: thwindM
      type(mean_sd_type) :: ciaM
      real(cp) :: timeLast
      real(cp) :: dt
      integer :: n_calls
      logical :: l_calc
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: calc_avg
      procedure :: mean_sd_loc
      procedure :: initialize_domdt
      procedure :: finalize_domdt
      procedure :: write_outputs
   end type vort_bal_type

contains

   subroutine initialize(this)
      !
      ! Memory allocation
      !

      class(vort_bal_type) :: this

      if ( l_finite_diff ) then
         this%mMin=1
         this%mMax=n_m_max
         this%nRmin=nRstart
         this%nRmax=nRstop
      else
         this%mMin=nMstart
         this%mMax=nMstop
         this%nRmin=1
         this%nRmax=n_r_max
      end if

      allocate( this%visc(this%mMin:this%mMax,this%nRmin:this%nRmax) )
      allocate( this%cor(this%mMin:this%mMax,this%nRmin:this%nRmax) )
      allocate( this%adv(this%mMin:this%mMax,this%nRmin:this%nRmax) )
      allocate( this%dwdt(this%mMin:this%mMax,this%nRmin:this%nRmax) )
      allocate( this%pump(this%mMin:this%mMax,this%nRmin:this%nRmax) )
      allocate( this%buo(this%mMin:this%mMax,this%nRmin:this%nRmax) )
      bytes_allocated = bytes_allocated + 6*(this%nRmax-this%nRmin+1)* &
      &                 (this%mMax-this%mMin+1)*SIZEOF_DEF_COMPLEX

      call this%visc2D%initialize(nMstart,nMstop,n_r_max,l_2D_SD)
      call this%cor2D%initialize(nMstart,nMstop,n_r_max,l_2D_SD)
      call this%dwdt2D%initialize(nMstart,nMstop,n_r_max,l_2D_SD)
      call this%pump2D%initialize(nMstart,nMstop,n_r_max,l_2D_SD)
      call this%adv2D%initialize(nMstart,nMstop,n_r_max,l_2D_SD)
      call this%buo2D%initialize(nMstart,nMstop,n_r_max,l_2D_SD)
      call this%iner2D%initialize(nMstart,nMstop,n_r_max,l_2D_SD)
      call this%thwind2D%initialize(nMstart,nMstop,n_r_max,l_2D_SD)
      call this%cia2D%initialize(nMstart,nMstop,n_r_max,l_2D_SD)

      this%visc(:,:)=zero
      this%cor(:,:) =zero
      this%adv(:,:) =zero
      this%dwdt(:,:)=zero
      this%buo(:,:) =zero
      this%pump(:,:)=zero

      this%n_calls = 0
      this%l_calc = .false.
      this%timeLast = 0.0_cp
      this%dt = 0.0_cp

      call this%viscM%initialize(nMstart,nMstop)
      call this%pumpM%initialize(nMstart,nMstop)
      call this%corM%initialize(nMstart,nMstop)
      call this%advM%initialize(nMstart,nMstop)
      call this%buoM%initialize(nMstart,nMstop)
      call this%dwdtM%initialize(nMstart,nMstop)
      call this%inerM%initialize(nMstart,nMstop)
      call this%thwindM%initialize(nMstart,nMstop)
      call this%ciaM%initialize(nMstart,nMstop)

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation
      !

      class(vort_bal_type) :: this

      call this%ciaM%finalize()
      call this%thwindM%finalize()
      call this%inerM%finalize()
      call this%dwdtM%finalize()
      call this%buoM%finalize()
      call this%advM%finalize()
      call this%corM%finalize()
      call this%pumpM%finalize()
      call this%viscM%finalize()
      call this%cor2D%finalize()
      call this%buo2D%finalize()
      call this%adv2D%finalize()
      call this%dwdt2D%finalize()
      call this%pump2D%finalize()
      call this%visc2D%finalize()
      call this%thwind2D%finalize()
      call this%iner2D%finalize()
      call this%cia2D%finalize()
      deallocate( this%visc, this%cor, this%adv, this%dwdt, this%pump, this%buo )

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine initialize_domdt(this, om, tscheme)

      class(vort_bal_type) :: this

      !-- Input variables
      complex(cp),         intent(in) :: om(this%mMin:this%mMax,this%nRmin:this%nRmax)
      class(type_tscheme), intent(in) :: tscheme

      this%dwdt(:,:)=om(:,:)/tscheme%dt(1)

   end subroutine initialize_domdt
!------------------------------------------------------------------------------
   subroutine finalize_domdt(this, om, tscheme)

      class(vort_bal_type) :: this

      !-- Input variables
      complex(cp),         intent(in) :: om(this%mMin:this%mMax,this%nRmin:this%nRmax)
      class(type_tscheme), intent(in) :: tscheme

      this%dwdt(:,:)=om(:,:)/tscheme%dt(1)-this%dwdt(:,:)

   end subroutine finalize_domdt
!------------------------------------------------------------------------------
   subroutine mean_sd_loc(this, time, input, output, outM)

      class(vort_bal_type) :: this

      !-- Input variables
      real(cp),    intent(in) :: time
      complex(cp), intent(in) :: input(this%mMin:this%mMax,this%nRmin:this%nRmax)

      !-- Output variables
      type(mean_sd_2D_type), intent(inout) :: output
      type(mean_sd_type),    intent(inout) :: outM

      !-- Local variables
      complex(cp) :: tmp(nMstart:nMstop, n_r_max)
      real(cp) :: dat(n_r_max), sd, val_m
      integer :: n_r, n_m, m

      if ( l_finite_diff ) then
         call r2m_single%transp_r2m(input, tmp)
      else
         tmp(:,:)=input(:,:)
      end if

      !-- This loop is not very cache-friendly but only done for some outputing
      do n_m=nMstart,nMstop
         m = idx2m(n_m)
         do n_r=1,n_r_max
            dat(n_r)=cc2real(tmp(n_m,n_r), m)
            dat(n_r)=dat(n_r)*r(n_r)*height(n_r)
            if ( output%l_SD ) then
               call getMSD2(output%mean(n_m,n_r), output%SD(n_m,n_r), dat(n_r), &
                    &       this%n_calls, this%dt, time)
            else
               call getMSD2(output%mean(n_m,n_r), sd, dat(n_r), this%n_calls, &
                    &       this%dt, time)
            end if

            if ( (r(n_r) >= r_cmb-bl_cut) .or. (r(n_r) <= r_icb+bl_cut) ) then
               dat(n_r)=0.0_cp
            end if
         end do
         val_m = simps(dat,r)
         !val_m = sqrt(val_m) ! I rather save the squared forces for a coherent 
         ! definition of forces: intcheb(mean_2D) = mean_1D; intcheb(SD_2D= SD_1D
         call getMSD2(outM%mean(n_m), outM%SD(n_m), val_m, this%n_calls,&
              &       this%dt, time)
      end do

   end subroutine mean_sd_loc
!------------------------------------------------------------------------------
   subroutine calc_avg(this, time, l_stop_time)

      !-- Input variables:
      real(cp), intent(in) :: time
      logical,  intent(in) :: l_stop_time
      class(vort_bal_type) :: this

      this%n_calls = this%n_calls+1
      this%dt = time-this%timeLast

      !------
      !-- Buoyancy term
      !------
      call this%mean_sd_loc(time,this%buo,this%buo2D,this%buoM)

      !------
      !-- Coriolis term
      !------
      call this%mean_sd_loc(time,this%cor,this%cor2D,this%corM)

      !------
      !-- Advection term
      !------
      call this%mean_sd_loc(time,this%adv,this%adv2D,this%advM)

      !------
      !-- d\omega/dt term
      !------
      call this%mean_sd_loc(time,this%dwdt,this%dwdt2D,this%dwdtM)

      !------
      !-- Viscous term
      !------
      call this%mean_sd_loc(time,this%visc,this%visc2D,this%viscM)

      !------
      !-- Ekman pumping term
      !------
      call this%mean_sd_loc(time,this%pump,this%pump2D,this%pumpM)

      !------
      !-- Thermal wind balance
      !------
      !-- Make use of pump as a temporary array here
      this%pump(:,:)=this%buo(:,:)+this%cor(:,:)
      call this%mean_sd_loc(time,this%pump,this%thwind2D,this%thwindM)

      !------
      !-- Inertial term: d\omega/dt + div( u \omega )
      !------
      this%pump(:,:)=-this%dwdt(:,:)+this%adv(:,:)
      call this%mean_sd_loc(time,this%pump,this%iner2D,this%inerM)

      !------
      !-- Coriolis-Inertia-Archimedian force balance
      !------
      this%pump(:,:)=-this%dwdt(:,:)+this%adv(:,:)+this%buo(:,:)+this%cor(:,:)
      call this%mean_sd_loc(time,this%pump,this%cia2D,this%ciaM)

      !-- Put zeros in this%pump again!
      this%pump(:,:)=zero

      this%timeLast = time

      if ( l_stop_time ) call this%write_outputs()

   end subroutine calc_avg
!------------------------------------------------------------------------------
   subroutine write_outputs(this)

      class(vort_bal_type) ::this

      !-- Local variables
      integer :: n_m,file_handle,m,n_p,n_r
      integer :: info, fh, version, header_size, filetype
      integer :: istat(MPI_STATUS_SIZE)
      integer :: arr_size(2), arr_loc_size(2), arr_start(2)
      integer(kind=MPI_OFFSET_KIND) :: disp
      character(len=100) :: forces_file
      integer :: displs(0:n_procs-1), recvcounts(0:n_procs-1)
      real(cp) :: buoM_global(n_m_max), buoSD_global(n_m_max)
      real(cp) :: corM_global(n_m_max), corSD_global(n_m_max)
      real(cp) :: advM_global(n_m_max), advSD_global(n_m_max)
      real(cp) :: dwdtM_global(n_m_max), dwdtSD_global(n_m_max)
      real(cp) :: viscM_global(n_m_max), viscSD_global(n_m_max)
      real(cp) :: pumpM_global(n_m_max), pumpSD_global(n_m_max)
      real(cp) :: ciaM_global(n_m_max), ciaSD_global(n_m_max)
      real(cp) :: thwindM_global(n_m_max), thwindSD_global(n_m_max)
      real(cp) :: inerM_global(n_m_max), inerSD_global(n_m_max)

      !----------
      !-- Force balance integrated over radii (ascii file written by rank0)
      !----------
      do n_p=0,n_procs-1
         recvcounts(n_p)=m_balance(n_p)%n_per_rank
      end do
      displs(0)=0
      do n_p=1,n_procs-1
         displs(n_p)=displs(n_p-1)+recvcounts(n_p-1)
      end do
      call MPI_GatherV(this%buoM%mean, nm_per_rank, MPI_DEF_REAL,   &
           &           buoM_global, recvcounts, displs,             &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%buoM%SD, nm_per_rank, MPI_DEF_REAL,     &
           &           buoSD_global, recvcounts, displs,            &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%corM%mean, nm_per_rank, MPI_DEF_REAL,   &
           &           corM_global, recvcounts, displs,             &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%corM%SD, nm_per_rank, MPI_DEF_REAL,     &
           &           corSD_global, recvcounts, displs,            &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%advM%mean, nm_per_rank, MPI_DEF_REAL,   &
           &           advM_global, recvcounts, displs,             &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%advM%SD, nm_per_rank, MPI_DEF_REAL,     &
           &           advSD_global, recvcounts, displs,            &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%dwdtM%mean, nm_per_rank, MPI_DEF_REAL,  &
           &           dwdtM_global, recvcounts, displs,            &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%dwdtM%SD, nm_per_rank, MPI_DEF_REAL,    &
           &           dwdtSD_global, recvcounts, displs,           &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%viscM%mean, nm_per_rank, MPI_DEF_REAL,  &
           &           viscM_global, recvcounts, displs,            &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%viscM%SD, nm_per_rank, MPI_DEF_REAL,    &
           &           viscSD_global, recvcounts, displs,           &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%pumpM%mean, nm_per_rank, MPI_DEF_REAL,  &
           &           pumpM_global, recvcounts, displs,            &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%pumpM%SD, nm_per_rank, MPI_DEF_REAL,    &
           &           pumpSD_global, recvcounts, displs,           &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%thwindM%mean, nm_per_rank, MPI_DEF_REAL,&
           &           thwindM_global, recvcounts, displs,          &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%thwindM%SD, nm_per_rank, MPI_DEF_REAL,  &
           &           thwindSD_global, recvcounts, displs,         &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%ciaM%mean, nm_per_rank, MPI_DEF_REAL,   &
           &           ciaM_global, recvcounts, displs,             &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%ciaM%SD, nm_per_rank, MPI_DEF_REAL,     &
           &           ciaSD_global, recvcounts, displs,            &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%inerM%mean, nm_per_rank, MPI_DEF_REAL,  &
           &           inerM_global, recvcounts, displs,            &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%inerM%SD, nm_per_rank, MPI_DEF_REAL,    &
           &           inerSD_global, recvcounts, displs,           &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)

      if ( rank == 0 ) then
         open(newunit=file_handle, file='vort_terms_avg.'//tag)
         do n_m=1,n_m_max
            m = idx2m(n_m)
            if ( m >= 1 ) then
               corSD_global(n_m)   =sqrt(corSD_global(n_m)/this%timeLast)
               buoSD_global(n_m)   =sqrt(buoSD_global(n_m)/this%timeLast)
               advSD_global(n_m)   =sqrt(advSD_global(n_m)/this%timeLast)
               dwdtSD_global(n_m)  =sqrt(dwdtSD_global(n_m)/this%timeLast)
               viscSD_global(n_m)  =sqrt(viscSD_global(n_m)/this%timeLast)
               pumpSD_global(n_m)  =sqrt(pumpSD_global(n_m)/this%timeLast)
               thwindSD_global(n_m)=sqrt(thwindSD_global(n_m)/this%timeLast)
               inerSD_global(n_m)  =sqrt(inerSD_global(n_m)/this%timeLast)
               ciaSD_global(n_m)   =sqrt(ciaSD_global(n_m)/this%timeLast)

               write(file_handle, '(I5,18es16.8)') m,                      &
               &                buoM_global(n_m), buoSD_global(n_m),       &
               &                corM_global(n_m), corSD_global(n_m),       &
               &                advM_global(n_m), advSD_global(n_m),       &
               &                dwdtM_global(n_m), dwdtSD_global(n_m),     &
               &                viscM_global(n_m), viscSD_global(n_m),     &
               &                pumpM_global(n_m), pumpSD_global(n_m),     &
               &                thwindM_global(n_m), thwindSD_global(n_m), &
               &                inerM_global(n_m), inerSD_global(n_m),     &
               &                ciaM_global(n_m), ciaSD_global(n_m)
            end if
         end do
         close(file_handle)
      end if

      !-----------
      !-- 2D force balance file (make use of MPI-IO)
      !-----------
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop

            if ( this%cor2D%l_SD ) then
               this%cor2D%SD(n_m,n_r)   =sqrt(this%cor2D%SD(n_m,n_r)/this%timeLast)
               this%buo2D%SD(n_m,n_r)   =sqrt(this%buo2D%SD(n_m,n_r)/this%timeLast)
               this%adv2D%SD(n_m,n_r)   =sqrt(this%adv2D%SD(n_m,n_r)/this%timeLast)
               this%dwdt2D%SD(n_m,n_r)  =sqrt(this%dwdt2D%SD(n_m,n_r)/this%timeLast)
               this%visc2D%SD(n_m,n_r)  =sqrt(this%visc2D%SD(n_m,n_r)/this%timeLast)
               this%pump2D%SD(n_m,n_r)  =sqrt(this%pump2D%SD(n_m,n_r)/this%timeLast)
               this%iner2D%SD(n_m,n_r)  =sqrt(this%iner2D%SD(n_m,n_r)/this%timeLast)
               this%thwind2D%SD(n_m,n_r)=sqrt(this%thwind2D%SD(n_m,n_r)/this%timeLast)
               this%cia2D%SD(n_m,n_r)   =sqrt(this%cia2D%SD(n_m,n_r)/this%timeLast)
            end if
         end do
      end do

      forces_file='vort_bal.'//tag
      if ( this%cor2D%l_SD ) then
         version = 2
      else
         version = 1
      end if

      header_size = SIZEOF_INTEGER+6*SIZEOF_DEF_REAL+4*SIZEOF_INTEGER   &
      &             +n_r_max*SIZEOF_DEF_REAL

      call MPI_Info_create(info, ierr)

      !-- Enable collective buffering
      call MPI_Info_set(info, "romio_cb_write", "automatic", ierr)
      call MPI_Info_set(info, "romio_cb_read", "automatic", ierr)

      !-- Disable data sieving (let the filesystem handles it)
      call MPI_Info_set(info, "romio_ds_write", "disable", ierr)
      call MPI_Info_set(info, "romio_ds_read", "disable", ierr)

      !-- Set the stripping unit to 4M
      call MPI_Info_set(info, "stripping_unit", "4194304", ierr)

      !-- Set the buffer size to 4M
      call MPI_Info_set(info,"cb_buffer_size","4194304", ierr)

      !-- Open file
      call MPI_File_Open(MPI_COMM_WORLD, forces_file, ior(MPI_MODE_WRONLY, &
           &             MPI_MODE_CREATE), info, fh, ierr)

      !-- Only rank=0 writes the header of the file
      if ( rank == 0 ) then
         call MPI_File_Write(fh, version, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, ra, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, pr, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, raxi, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, sc, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, ek, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, radratio, 1, MPI_DEF_REAL, istat, ierr)

         call MPI_File_Write(fh, n_r_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, n_m_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, m_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, minc, 1, MPI_INTEGER, istat, ierr)

         call MPI_File_Write(fh, r, n_r_max, MPI_DEF_REAL, istat, ierr)
      end if

      arr_size(1) = n_m_max
      arr_size(2) = n_r_max
      arr_loc_size(1) = nm_per_rank
      arr_loc_size(2) = n_r_max
      arr_start(1) = nMstart-1
      arr_start(2) = 0
      call MPI_Type_Create_Subarray(2,arr_size, arr_loc_size, arr_start, &
           &                        MPI_ORDER_FORTRAN, MPI_DEF_REAL,     &
           &                        filetype, ierr)
      call MPI_Type_Commit(filetype, ierr)

      !-- Set the view after the header
      disp = header_size
      call MPI_File_Set_View(fh, disp, MPI_DEF_REAL, filetype, "native", &
           &                 info, ierr)

      !-- Now finally write the fields
      call MPI_File_Write_all(fh, this%buo2D%mean, nm_per_rank*n_r_max,   &
           &                  MPI_DEF_REAL, istat, ierr)
      if ( this%buo2D%l_SD ) then
         call MPI_File_Write_all(fh, this%buo2D%SD, nm_per_rank*n_r_max,  &
              &                  MPI_DEF_REAL, istat, ierr)
      end if
      call MPI_File_Write_all(fh, this%cor2D%mean, nm_per_rank*n_r_max,   &
           &                  MPI_DEF_REAL, istat, ierr)
      if ( this%cor2D%l_SD ) then
         call MPI_File_Write_all(fh, this%cor2D%SD, nm_per_rank*n_r_max,  &
              &                  MPI_DEF_REAL, istat, ierr)
      end if
      call MPI_File_Write_all(fh, this%adv2D%mean, nm_per_rank*n_r_max,   &
           &                  MPI_DEF_REAL, istat, ierr)
      if ( this%adv2D%l_SD ) then
         call MPI_File_Write_all(fh, this%adv2D%SD, nm_per_rank*n_r_max,  &
              &                  MPI_DEF_REAL, istat, ierr)
      end if
      call MPI_File_Write_all(fh, this%dwdt2D%mean, nm_per_rank*n_r_max,  &
           &                  MPI_DEF_REAL, istat, ierr)
      if ( this%dwdt2D%l_SD ) then
         call MPI_File_Write_all(fh, this%dwdt2D%SD, nm_per_rank*n_r_max, &
              &                  MPI_DEF_REAL, istat, ierr)
      end if
      call MPI_File_Write_all(fh, this%visc2D%mean, nm_per_rank*n_r_max,  &
           &                  MPI_DEF_REAL, istat, ierr)
      if ( this%visc2D%l_SD ) then
         call MPI_File_Write_all(fh, this%visc2D%SD, nm_per_rank*n_r_max, &
              &                  MPI_DEF_REAL, istat, ierr)
      end if
      call MPI_File_Write_all(fh, this%pump2D%mean, nm_per_rank*n_r_max,  &
           &                  MPI_DEF_REAL, istat, ierr)
      if ( this%pump2D%l_SD ) then
         call MPI_File_Write_all(fh, this%pump2D%SD, nm_per_rank*n_r_max, &
              &                  MPI_DEF_REAL, istat, ierr)
      end if
      call MPI_File_Write_all(fh, this%thwind2D%mean, nm_per_rank*n_r_max,&
           &                  MPI_DEF_REAL, istat, ierr)
      if ( this%thwind2D%l_SD ) then
         call MPI_File_Write_all(fh, this%thwind2D%SD,nm_per_rank*n_r_max,&
              &                  MPI_DEF_REAL, istat, ierr)
      end if
      call MPI_File_Write_all(fh, this%iner2D%mean, nm_per_rank*n_r_max,  &
           &                  MPI_DEF_REAL, istat, ierr)
      if ( this%iner2D%l_SD ) then
         call MPI_File_Write_all(fh, this%iner2D%SD, nm_per_rank*n_r_max, &
              &                  MPI_DEF_REAL, istat, ierr)
      end if
      call MPI_File_Write_all(fh, this%cia2D%mean, nm_per_rank*n_r_max,   &
           &                  MPI_DEF_REAL, istat, ierr)
      if ( this%cia2D%l_SD ) then
         call MPI_File_Write_all(fh, this%cia2D%SD, nm_per_rank*n_r_max,  &
              &                  MPI_DEF_REAL, istat, ierr)
      end if

      call MPI_Info_free(info, ierr)
      call MPI_File_close(fh, ierr)

   end subroutine write_outputs
!------------------------------------------------------------------------------
end module vort_balance
