module spectra
   !
   ! This module handles the calculation and the writing of spectra.
   !   -spec_#.TAG files for snapshot spectra
   !   -spec_avg.TAG files for time-averaged spectra (and standard deviation)
   !   -2D_spec_avg.TAG files which correspond to 2-D spectra in a (r,m) plane

   use precision_mod
   use parallel_mod
   use constants, only: pi, one
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, idx2m, n_m_max, minc, m_max
   use blocking, only: nMstart, nMstop, nm_per_rank, m_balance
   use namelists, only: tag, ra, ek, pr, ek, raxi, sc, radratio, l_2D_spectra
   use useful, only: cc2real, getMSD2, round_off
   use radial_functions, only: r, rscheme, height
   use integration, only: rInt_R

   implicit none

   private

   type, public :: spectra_type
      real(cp), allocatable :: us2_mean(:,:)
      real(cp), allocatable :: up2_mean(:,:)
      real(cp), allocatable :: enst_mean(:,:)
      real(cp), allocatable :: us2M_mean(:)
      real(cp), allocatable :: us2M_SD(:)
      real(cp), allocatable :: up2M_mean(:)
      real(cp), allocatable :: up2M_SD(:)
      real(cp), allocatable :: enstM_mean(:)
      real(cp), allocatable :: enstM_SD(:)
      integer :: n_calls
      integer :: ispec_counter
      real(cp) :: dt
      real(cp) :: timeLast
      logical :: l_calc
   contains 
      procedure :: initialize
      procedure :: finalize
      procedure :: write_spectra
      procedure :: write_spectra_avg
      procedure :: calculate_spectra
   end type spectra_type

contains

   subroutine initialize(this)
      !
      ! Memory allocation and initial values
      !
      class(spectra_type) :: this

      !-- Local variables:
      integer :: n_r, n_m

      if ( l_2D_spectra ) then
         allocate( this%us2_mean(nMstart:nMstop,n_r_max) )
         allocate( this%up2_mean(nMstart:nMstop,n_r_max) )
         allocate( this%enst_mean(nMstart:nMstop,n_r_max) )
         bytes_allocated = bytes_allocated+3*(nMstop-nMstart+1)*n_r_max*&
         &                 SIZEOF_DEF_REAL

         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               this%us2_mean(n_m,n_r) =0.0_cp
               this%up2_mean(n_m,n_r) =0.0_cp
               this%enst_mean(n_m,n_r)=0.0_cp
            end do
         end do
      end if

      allocate( this%us2M_mean(nMstart:nMstop), this%up2M_mean(nMstart:nMstop) )
      allocate( this%enstM_mean(nMstart:nMstop), this%us2M_SD(nMstart:nMstop) )
      allocate( this%up2M_SD(nMstart:nMstop), this%enstM_SD(nMstart:nMstop) )
      bytes_allocated=bytes_allocated+6*(nMstop-nMstart+1)*SIZEOF_DEF_REAL

      this%us2M_mean(:)  = 0.0_cp
      this%up2M_mean(:)  = 0.0_cp
      this%enstM_mean(:) = 0.0_cp
      this%us2M_SD(:)    = 0.0_cp
      this%up2M_SD(:)    = 0.0_cp
      this%enstM_SD(:)   = 0.0_cp

      this%ispec_counter = 1
      this%n_calls = 0
      this%dt = 0.0_cp
      this%timeLast = 0.0_cp
      this%l_calc = .false.

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)
      !
      ! Memory deallocation
      !
      class(spectra_type) :: this

      deallocate(this%us2M_mean, this%us2M_SD, this%up2M_mean)
      deallocate(this%up2M_SD, this%enstM_mean, this%enstM_SD)
      if ( l_2D_spectra ) then
         deallocate(this%enst_mean, this%us2_mean, this%up2_mean)
      end if

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine calculate_spectra(this, time, l_stop_time, us_Mloc, up_Mloc, &
              &                 om_Mloc, us2_m, up2_m, enst_m)

      class(spectra_type) :: this

      !-- Input variables:
      real(cp),    intent(in) :: time
      logical,     intent(in) :: l_stop_time
      complex(cp), intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: om_Mloc(nMstart:nMstop,n_r_max)

      !-- Output variables
      real(cp), intent(out) :: us2_m(nMstart:nMstop)
      real(cp), intent(out) :: up2_m(nMstart:nMstop)
      real(cp), intent(out) :: enst_m(nMstart:nMstop)

      !-- Local variables:
      real(cp) :: us2(n_r_max), up2(n_r_max), enst(n_r_max)
      real(cp) :: sd
      integer :: n_r, n_m, m

      this%n_calls = this%n_calls+1
      this%dt = time-this%timeLast

      !-- This is not cache-friendly but hopefully it's happening only
      !-- once in a while (otherwise we need (n_r, n_m) arrays
      do n_m=nMstart,nMstop
         m = idx2m(n_m)
         do n_r=1,n_r_max
            us2(n_r) =cc2real(us_Mloc(n_m,n_r),m)
            up2(n_r) =cc2real(up_Mloc(n_m,n_r),m)
            enst(n_r)=cc2real(om_Mloc(n_m,n_r),m)
            us2(n_r) =pi*us2(n_r)*r(n_r)*height(n_r)
            up2(n_r) =pi*up2(n_r)*r(n_r)*height(n_r)
            enst(n_r)=pi*enst(n_r)*r(n_r)*height(n_r)
            if ( l_2D_spectra ) then
               call getMSD2(this%us2_mean(n_m,n_r), sd, us2(n_r), this%n_calls, &
                    &       this%dt, time)
               call getMSD2(this%up2_mean(n_m,n_r), sd, up2(n_r), this%n_calls, &
                    &       this%dt, time)
               call getMSD2(this%enst_mean(n_m,n_r), sd, enst(n_r), this%n_calls, &
                    &       this%dt, time)
            end if
         end do
         !--Radial integration
         us2_m(n_m) =rInt_R(us2, r, rscheme)
         up2_m(n_m) =rInt_R(up2, r, rscheme)
         enst_m(n_m)=rInt_R(enst, r, rscheme)

         !-- Mean and SD of m-spectra
         call getMSD2(this%us2M_mean(n_m), this%us2M_SD(n_m), us2_m(n_m), &
              &       this%n_calls, this%dt, time)
         call getMSD2(this%up2M_mean(n_m), this%up2M_SD(n_m), up2_m(n_m), &
              &       this%n_calls, this%dt, time)
         call getMSD2(this%enstM_mean(n_m), this%enstM_SD(n_m), enst_m(n_m), &
              &       this%n_calls, this%dt, time)
      end do
      this%timeLast = time

      if ( l_stop_time ) call this%write_spectra_avg()

   end subroutine calculate_spectra
!----------------------------------------------------------------------
   subroutine write_spectra_avg(this)
      !
      ! This subroutine writes the time-averages spectra: spec_avg.TAG
      !

      class(spectra_type) :: this

      !-- Local variables
      real(cp) :: us2_mean_global(n_m_max), us2_SD_global(n_m_max)
      real(cp) :: up2_mean_global(n_m_max), up2_SD_global(n_m_max)
      real(cp) :: enst_mean_global(n_m_max), enst_SD_global(n_m_max)
      integer :: file_handle, n_m, m, n_p
      integer :: displs(0:n_procs-1), recvcounts(0:n_procs-1)

      !----------------
      !- First write the time-average m spectra
      !----------------

      !-- Finally gather everything on rank=0
      do n_p=0,n_procs-1
         recvcounts(n_p)=m_balance(n_p)%n_per_rank
      end do
      displs(0)=0
      do n_p=1,n_procs-1
         displs(n_p)=displs(n_p-1)+recvcounts(n_p-1)
      end do
      call MPI_GatherV(this%us2M_mean, nm_per_rank, MPI_DEF_REAL,  &
           &           us2_mean_global, recvcounts, displs,        &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%us2M_SD, nm_per_rank, MPI_DEF_REAL,    &
           &           us2_SD_global, recvcounts, displs,          &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%up2M_mean, nm_per_rank, MPI_DEF_REAL,  &
           &           up2_mean_global, recvcounts, displs,        &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%up2M_SD, nm_per_rank, MPI_DEF_REAL,    &
           &           up2_SD_global, recvcounts, displs,          &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%enstM_mean, nm_per_rank, MPI_DEF_REAL, &
           &           enst_mean_global, recvcounts, displs,       &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%enstM_SD, nm_per_rank, MPI_DEF_REAL,   &
           &           enst_SD_global, recvcounts, displs,         &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)

      !-- Only rank==0 writes the spec_avg.TAG file
      if ( rank == 0 ) then
         open(newunit=file_handle, file='spec_avg.'//tag)
         do n_m=1,n_m_max
            m = idx2m(n_m)
            print*, us2_SD_global(n_m)
            us2_SD_global(n_m) =sqrt(us2_SD_global(n_m)/this%timeLast)
            up2_SD_global(n_m) =sqrt(up2_SD_global(n_m)/this%timeLast)
            enst_SD_global(n_m)=sqrt(enst_SD_global(n_m)/this%timeLast)
            write(file_handle, '(I4, 6es16.8)') m,                                &
            &     round_off(us2_mean_global(n_m)), round_off(us2_SD_global(n_m)), &
            &     round_off(up2_mean_global(n_m)), round_off(up2_SD_global(n_m)), &
            &     round_off(enst_mean_global(n_m)), round_off(enst_SD_global(n_m))
         end do
         close(file_handle)
      end if

      !-------------------
      !- 2nd: store the time-average 2-D spectra in a (r,m) plane)
      !-------------------
      if ( l_2D_spectra ) then
         call write_2D_spectra(this%us2_mean, this%up2_mean, this%enst_mean)
      end if

   end subroutine write_spectra_avg
!----------------------------------------------------------------------
   subroutine write_spectra(this, us2_m, up2_m, enst_m)
      !
      ! This subroutine writes one snapshot spectrum
      !

      class(spectra_type) :: this

      !-- Input variables
      real(cp), intent(in) :: us2_m(nMstart:nMstop)
      real(cp), intent(in) :: up2_m(nMstart:nMstop)
      real(cp), intent(in) :: enst_m(nMstart:nMstop)

      !-- Local variables
      real(cp) :: us2_m_global(n_m_max), up2_m_global(n_m_max)
      real(cp) :: enst_m_global(n_m_max)
      integer :: n_p, m, n_m
      integer :: displs(0:n_procs-1), recvcounts(0:n_procs-1)
      character(len=144) :: spec_name
      integer :: file_handle

      do n_p=0,n_procs-1
         recvcounts(n_p)=m_balance(n_p)%n_per_rank
      end do
      displs(0)=0
      do n_p=1,n_procs-1
         displs(n_p)=displs(n_p-1)+recvcounts(n_p-1)
      end do
      call MPI_GatherV(us2_m, nm_per_rank, MPI_DEF_REAL,        &
           &           us2_m_global, recvcounts, displs,        &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(up2_m, nm_per_rank, MPI_DEF_REAL,        &
           &           up2_m_global, recvcounts, displs,        &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(enst_m, nm_per_rank, MPI_DEF_REAL,       &
           &           enst_m_global, recvcounts, displs,       &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)

      if ( rank == 0 ) then
         write(spec_name, '(A,I0,A,A)') 'spec_',this%ispec_counter,'.',tag

         open(newunit=file_handle, file=spec_name, position='append')
         do n_m=1,n_m_max
            m = idx2m(n_m)
            write(file_handle, '(I4, 3es16.8)') m,       &
            &              round_off(us2_m_global(n_m)), &
            &              round_off(up2_m_global(n_m)), &
            &              round_off(enst_m_global(n_m))
         end do
         close(file_handle)

         this%ispec_counter = this%ispec_counter+1

      end if

   end subroutine write_spectra
!----------------------------------------------------------------------
   subroutine write_2D_spectra(us2_mean, up2_mean, enst_mean)

      !-- Input variables:
      real(cp), intent(in) :: us2_mean(nMstart:nMstop,n_r_max)
      real(cp), intent(in) :: up2_mean(nMstart:nMstop,n_r_max)
      real(cp), intent(in) :: enst_mean(nMstart:nMstop,n_r_max)

      !-- Local variables:
      integer :: info, fh, version, header_size
      integer :: istat(MPI_STATUS_SIZE), filetype
      integer :: arr_size(2), arr_loc_size(2), arr_start(2)
      integer(kind=MPI_OFFSET_KIND) :: disp
      character(len=100) :: spec_file

      spec_file='2D_spec_avg.'//tag
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
      call MPI_File_Open(MPI_COMM_WORLD, spec_file, ior(MPI_MODE_WRONLY, &
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
      call MPI_File_Write_all(fh, us2_mean, nm_per_rank*n_r_max,  & 
           &                  MPI_DEF_REAL, istat, ierr)
      call MPI_File_Write_all(fh, up2_mean, nm_per_rank*n_r_max,  &
           &                  MPI_DEF_REAL, istat, ierr)
      call MPI_File_Write_all(fh, enst_mean, nm_per_rank*n_r_max, &
           &                  MPI_DEF_REAL, istat, ierr)

      call MPI_Info_free(info, ierr)
      call MPI_File_close(fh, ierr)

   end subroutine write_2D_spectra
!----------------------------------------------------------------------
end module spectra
