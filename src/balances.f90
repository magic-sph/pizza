module balances
   !
   ! This module is used to compute the balance of the terms that enter
   ! the vorticity equations.
   !

   use precision_mod
   use parallel_mod
   use constants, only: zero
   use time_schemes, only: type_tscheme
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, m_max, minc, n_m_max, idx2m
   use blocking, only: nMstart, nMstop, nm_per_rank, m_balance
   use namelists, only: ra, pr, ek, radratio, sc, raxi, tag, r_cmb, r_icb, &
       &                bl_cut
   use radial_functions, only: r, height
   use integration, only: simps
   use useful, only: cc2real, getMSD2

   implicit none

   private

   type, public :: vort_bal_type
      complex(cp), allocatable :: visc(:,:)
      complex(cp), allocatable :: cor(:,:)
      complex(cp), allocatable :: adv(:,:)
      complex(cp), allocatable :: pump(:,:)
      complex(cp), allocatable :: dwdt(:,:)
      complex(cp), allocatable :: buo(:,:)
      real(cp), allocatable :: visc_mean(:,:)
      real(cp), allocatable :: cor_mean(:,:)
      real(cp), allocatable :: pump_mean(:,:)
      real(cp), allocatable :: adv_mean(:,:)
      real(cp), allocatable :: dwdt_mean(:,:)
      real(cp), allocatable :: buo_mean(:,:)
      real(cp), allocatable :: iner_mean(:,:)
      real(cp), allocatable :: thwind_mean(:,:)
      real(cp), allocatable :: cia_mean(:,:)
      real(cp), allocatable :: viscM_mean(:)
      real(cp), allocatable :: viscM_SD(:)
      real(cp), allocatable :: pumpM_mean(:)
      real(cp), allocatable :: pumpM_SD(:)
      real(cp), allocatable :: buoM_mean(:)
      real(cp), allocatable :: buoM_SD(:)
      real(cp), allocatable :: corM_mean(:)
      real(cp), allocatable :: corM_SD(:)
      real(cp), allocatable :: advM_mean(:)
      real(cp), allocatable :: advM_SD(:)
      real(cp), allocatable :: dwdtM_mean(:)
      real(cp), allocatable :: dwdtM_SD(:)
      real(cp), allocatable :: thwindM_mean(:)
      real(cp), allocatable :: thwindM_SD(:)
      real(cp), allocatable :: inerM_SD(:)
      real(cp), allocatable :: inerM_mean(:)
      real(cp), allocatable :: ciaM_SD(:)
      real(cp), allocatable :: ciaM_mean(:)
      real(cp) :: timeLast
      real(cp) :: dt
      integer :: n_calls
      logical :: l_calc
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: calc_avg
      procedure :: mean_sd
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

      !-- Local variables
      integer :: n_r, n_m

      allocate( this%visc(nMstart:nMstop,n_r_max) )
      allocate( this%cor(nMstart:nMstop,n_r_max) )
      allocate( this%adv(nMstart:nMstop,n_r_max) )
      allocate( this%dwdt(nMstart:nMstop,n_r_max) )
      allocate( this%pump(nMstart:nMstop,n_r_max) )
      allocate( this%buo(nMstart:nMstop,n_r_max) )
      allocate( this%visc_mean(nMstart:nMstop,n_r_max) )
      allocate( this%cor_mean(nMstart:nMstop,n_r_max) )
      allocate( this%dwdt_mean(nMstart:nMstop,n_r_max) )
      allocate( this%pump_mean(nMstart:nMstop,n_r_max) )
      allocate( this%adv_mean(nMstart:nMstop,n_r_max) )
      allocate( this%buo_mean(nMstart:nMstop,n_r_max) )
      allocate( this%iner_mean(nMstart:nMstop,n_r_max) )
      allocate( this%cia_mean(nMstart:nMstop,n_r_max) )
      allocate( this%thwind_mean(nMstart:nMstop,n_r_max) )
      bytes_allocated = bytes_allocated + 6*n_r_max*(nMstop-nMstart+1)* &
      &                 SIZEOF_DEF_COMPLEX+9*n_r_max*(nMstop-nMstart+1)*&
      &                 SIZEOF_DEF_REAL

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            this%visc(n_m,n_r)=zero
            this%cor(n_m,n_r) =zero
            this%adv(n_m,n_r) =zero
            this%dwdt(n_m,n_r)=zero
            this%buo(n_m,n_r) =zero
            this%pump(n_m,n_r)=zero
            this%visc_mean(n_m,n_r)  =0.0_cp
            this%buo_mean(n_m,n_r)   =0.0_cp
            this%cor_mean(n_m,n_r)   =0.0_cp
            this%adv_mean(n_m,n_r)   =0.0_cp
            this%pump_mean(n_m,n_r)  =0.0_cp
            this%dwdt_mean(n_m,n_r)  =0.0_cp
            this%iner_mean(n_m,n_r)  =0.0_cp
            this%thwind_mean(n_m,n_r)=0.0_cp
            this%cia_mean(n_m,n_r)   =0.0_cp
         end do
      end do

      this%n_calls = 0
      this%l_calc = .false.
      this%timeLast = 0.0_cp
      this%dt = 0.0_cp

      allocate( this%viscM_mean(nMstart:nMstop), this%viscM_SD(nMstart:nMstop) )
      allocate( this%pumpM_mean(nMstart:nMstop), this%pumpM_SD(nMstart:nMstop) )
      allocate( this%corM_mean(nMstart:nMstop), this%corM_SD(nMstart:nMstop) )
      allocate( this%advM_mean(nMstart:nMstop), this%advM_SD(nMstart:nMstop) )
      allocate( this%buoM_mean(nMstart:nMstop), this%buoM_SD(nMstart:nMstop) )
      allocate( this%dwdtM_mean(nMstart:nMstop), this%dwdtM_SD(nMstart:nMstop) )
      allocate( this%thwindM_mean(nMstart:nMstop), this%thwindM_SD(nMstart:nMstop) )
      allocate( this%ciaM_mean(nMstart:nMstop), this%ciaM_SD(nMstart:nMstop) )
      allocate( this%inerM_mean(nMstart:nMstop), this%inerM_SD(nMstart:nMstop) )
      bytes_allocated=bytes_allocated+18*(nMstop-nMstart+1)*SIZEOF_DEF_REAL
      this%viscM_mean(:)  =0.0_cp
      this%viscM_SD(:)    =0.0_cp
      this%pumpM_mean(:)  =0.0_cp
      this%pumpM_SD(:)    =0.0_cp
      this%corM_mean(:)   =0.0_cp
      this%corM_SD(:)     =0.0_cp
      this%buoM_mean(:)   =0.0_cp
      this%buoM_SD(:)     =0.0_cp
      this%advM_mean(:)   =0.0_cp
      this%advM_SD(:)     =0.0_cp
      this%dwdtM_mean(:)  =0.0_cp
      this%dwdtM_SD(:)    =0.0_cp
      this%thwindM_mean(:)=0.0_cp
      this%thwindM_SD(:)  =0.0_cp
      this%inerM_mean(:)  =0.0_cp
      this%inerM_SD(:)    =0.0_cp
      this%ciaM_mean(:)   =0.0_cp
      this%ciaM_SD(:)     =0.0_cp

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation
      !

      class(vort_bal_type) :: this

      deallocate( this%viscM_mean, this%viscM_SD )
      deallocate( this%pumpM_mean, this%pumpM_SD )
      deallocate( this%buoM_mean, this%buoM_SD )
      deallocate( this%advM_mean, this%advM_SD )
      deallocate( this%dwdtM_mean, this%dwdtM_SD )
      deallocate( this%corM_mean, this%corM_SD )
      deallocate( this%thwindM_mean, this%thwindM_SD )
      deallocate( this%ciaM_mean, this%ciaM_SD )
      deallocate( this%inerM_mean, this%inerM_SD )
      deallocate( this%visc_mean, this%cor_mean, this%adv_mean, this%dwdt_mean )
      deallocate( this%pump_mean, this%buo_mean, this%iner_mean, this%cia_mean )
      deallocate( this%thwind_mean)
      deallocate( this%visc, this%cor, this%adv, this%dwdt, this%pump, this%buo )

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine initialize_domdt(this, om_Mloc, tscheme)

      class(vort_bal_type) :: this

      !-- Input variables
      complex(cp),         intent(in) :: om_Mloc(nMstart:nMstop,n_r_max)
      class(type_tscheme), intent(in) :: tscheme

      !-- Local variables:
      integer :: n_r, n_m

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            this%dwdt(n_m,n_r)=om_Mloc(n_m,n_r)/tscheme%dt(1)
         end do
      end do

   end subroutine initialize_domdt
!------------------------------------------------------------------------------
   subroutine finalize_domdt(this, om_Mloc, tscheme)

      class(vort_bal_type) :: this

      !-- Input variables
      complex(cp),         intent(in) :: om_Mloc(nMstart:nMstop,n_r_max)
      class(type_tscheme), intent(in) :: tscheme

      !-- Local variables:
      integer :: n_r, n_m

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            this%dwdt(n_m,n_r)=om_Mloc(n_m,n_r)/tscheme%dt(1)-this%dwdt(n_m,n_r)
         end do
      end do

   end subroutine finalize_domdt
!------------------------------------------------------------------------------
   subroutine mean_sd(this, input, output, outM_mean, outM_SD)

      class(vort_bal_type) :: this

      !-- Input variables
      complex(cp), intent(in) :: input(nMstart:nMstop,n_r_max)

      !-- Output variables
      real(cp), intent(out) :: output(nMstart:nMstop,n_r_max)
      real(cp), intent(out) :: outM_mean(nMstart:nMstop)
      real(cp), intent(out) :: outM_SD(nMstart:nMstop)

      !-- Local variables
      real(cp) :: dat(n_r_max), sd, val_m
      integer :: n_r, n_m, m

      !-- This loop is not very cache-friendly but only done for some outputing
      do n_m=nMstart,nMstop
         m = idx2m(n_m)
         do n_r=1,n_r_max
            dat(n_r) = cc2real(input(n_m,n_r), m)
            call getMSD2(output(n_m,n_r), sd, dat(n_r), this%n_calls, &
                 &       this%dt, this%timeLast)

            if ( (r(n_r) >= r_cmb-bl_cut) .or. (r(n_r) <= r_icb+bl_cut) ) then
               dat(n_r)=0.0_cp
            else
               dat(n_r)=dat(n_r)*r(n_r)*height(n_r)
            end if
         end do
         val_m = simps(dat,r)
         val_m = sqrt(val_m)
         call getMSD2(outM_mean(n_m), outM_SD(n_m), val_m, this%n_calls,&
              &       this%dt, this%timeLast)
      end do

   end subroutine mean_sd
!------------------------------------------------------------------------------
   subroutine calc_avg(this, time, l_stop_time)

      !-- Input variables:
      real(cp), intent(in) :: time
      logical,  intent(in) :: l_stop_time
      class(vort_bal_type) :: this

      !-- Local variables:
      integer :: n_r, n_m

      this%n_calls = this%n_calls+1
      this%dt = time-this%timeLast

      !------
      !-- Buoyancy term
      !------
      call this%mean_sd(this%buo,this%buo_mean,this%buoM_mean,this%buoM_SD)

      !------
      !-- Coriolis term
      !------
      call this%mean_sd(this%cor,this%cor_mean,this%corM_mean,this%corM_SD)

      !------
      !-- Advection term
      !------
      call this%mean_sd(this%adv,this%adv_mean,this%advM_mean,this%advM_SD)

      !------
      !-- d\omega/dt term
      !------
      call this%mean_sd(this%dwdt,this%dwdt_mean,this%dwdtM_mean,this%dwdtM_SD)

      !------
      !-- Viscous term
      !------
      call this%mean_sd(this%visc,this%visc_mean,this%viscM_mean,this%viscM_SD)

      !------
      !-- Ekman pumping term
      !------
      call this%mean_sd(this%pump,this%pump_mean,this%pumpM_mean,this%pumpM_SD)

      !------
      !-- Thermal wind balance
      !------
      !-- Make use of pump as a temporary array here
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            this%pump(n_m,n_r)=this%buo(n_m,n_r)+this%cor(n_m,n_r)
         end do
      end do
      call this%mean_sd(this%pump,this%thwind_mean,this%thwindM_mean,this%thwindM_SD)

      !------
      !-- Inertial term: d\omega/dt + div( u \omega )
      !------
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            this%pump(n_m,n_r)=-this%dwdt(n_m,n_r)+this%adv(n_m,n_r)
         end do
      end do
      call this%mean_sd(this%pump,this%iner_mean,this%inerM_mean,this%inerM_SD)

      !------
      !-- Coriolis-Inertia-Archimedian force balance
      !------
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            this%pump(n_m,n_r)=-this%dwdt(n_m,n_r)+this%adv(n_m,n_r)+ &
            &                   this%buo(n_m,n_r)+this%cor(n_m,n_r)
         end do
      end do
      call this%mean_sd(this%pump,this%cia_mean,this%ciaM_mean,this%ciaM_SD)

      !-- Put zeros in this%pump again!
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            this%pump(n_m,n_r)=zero
         end do
      end do

      this%timeLast = time

      if ( l_stop_time ) call this%write_outputs()

   end subroutine calc_avg
!------------------------------------------------------------------------------
   subroutine write_outputs(this)

      class(vort_bal_type) ::this

      !-- Local variables
      integer :: n_m,file_handle,m,n_p
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
      call MPI_GatherV(this%buoM_mean, nm_per_rank, MPI_DEF_REAL,   &
           &           buoM_global, recvcounts, displs,             &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%buoM_SD, nm_per_rank, MPI_DEF_REAL,     &
           &           buoSD_global, recvcounts, displs,            &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%corM_mean, nm_per_rank, MPI_DEF_REAL,   &
           &           corM_global, recvcounts, displs,             &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%corM_SD, nm_per_rank, MPI_DEF_REAL,     &
           &           corSD_global, recvcounts, displs,            &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%advM_mean, nm_per_rank, MPI_DEF_REAL,   &
           &           advM_global, recvcounts, displs,             &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%advM_SD, nm_per_rank, MPI_DEF_REAL,     &
           &           advSD_global, recvcounts, displs,            &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%dwdtM_mean, nm_per_rank, MPI_DEF_REAL,  &
           &           dwdtM_global, recvcounts, displs,            &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%dwdtM_SD, nm_per_rank, MPI_DEF_REAL,    &
           &           dwdtSD_global, recvcounts, displs,           &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%viscM_mean, nm_per_rank, MPI_DEF_REAL,  &
           &           viscM_global, recvcounts, displs,            &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%viscM_SD, nm_per_rank, MPI_DEF_REAL,    &
           &           viscSD_global, recvcounts, displs,           &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%pumpM_mean, nm_per_rank, MPI_DEF_REAL,  &
           &           pumpM_global, recvcounts, displs,            &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%pumpM_SD, nm_per_rank, MPI_DEF_REAL,    &
           &           pumpSD_global, recvcounts, displs,           &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%thwindM_mean, nm_per_rank, MPI_DEF_REAL,&
           &           thwindM_global, recvcounts, displs,          &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%thwindM_SD, nm_per_rank, MPI_DEF_REAL,  &
           &           thwindSD_global, recvcounts, displs,         &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%ciaM_mean, nm_per_rank, MPI_DEF_REAL,   &
           &           ciaM_global, recvcounts, displs,             &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%ciaM_SD, nm_per_rank, MPI_DEF_REAL,     &
           &           ciaSD_global, recvcounts, displs,            &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%inerM_mean, nm_per_rank, MPI_DEF_REAL,  &
           &           inerM_global, recvcounts, displs,            &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%inerM_SD, nm_per_rank, MPI_DEF_REAL,    &
           &           inerSD_global, recvcounts, displs,           &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)

      if ( rank == 0 ) then
         open(newunit=file_handle, file='vort_terms_avg.'//tag)
         do n_m=1,n_m_max
            m = idx2m(n_m)
            if ( m >= 1 ) then
               corSD_global(n_m) = sqrt(corSD_global(n_m)/this%timeLast)

               write(file_handle, '(I4,18es16.8)') m,                      &
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
      forces_file='vort_bal.'//tag
      version = 1

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
      call MPI_File_Write_all(fh, this%buo_mean, nm_per_rank*n_r_max, &
           &                  MPI_DEF_REAL, istat, ierr)
      call MPI_File_Write_all(fh, this%cor_mean, nm_per_rank*n_r_max, &
           &                  MPI_DEF_REAL, istat, ierr)
      call MPI_File_Write_all(fh, this%adv_mean, nm_per_rank*n_r_max, &
           &                  MPI_DEF_REAL, istat, ierr)
      call MPI_File_Write_all(fh, this%dwdt_mean, nm_per_rank*n_r_max, &
           &                  MPI_DEF_REAL, istat, ierr)
      call MPI_File_Write_all(fh, this%visc_mean, nm_per_rank*n_r_max, &
           &                  MPI_DEF_REAL, istat, ierr)
      call MPI_File_Write_all(fh, this%pump_mean, nm_per_rank*n_r_max, &
           &                  MPI_DEF_REAL, istat, ierr)
      call MPI_File_Write_all(fh, this%thwind_mean, nm_per_rank*n_r_max, &
           &                  MPI_DEF_REAL, istat, ierr)
      call MPI_File_Write_all(fh, this%iner_mean, nm_per_rank*n_r_max, &
           &                  MPI_DEF_REAL, istat, ierr)
      call MPI_File_Write_all(fh, this%cia_mean, nm_per_rank*n_r_max, &
           &                  MPI_DEF_REAL, istat, ierr)

      call MPI_Info_free(info, ierr)
      call MPI_File_close(fh, ierr)

   end subroutine write_outputs
!------------------------------------------------------------------------------
end module balances
