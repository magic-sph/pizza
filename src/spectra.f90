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
   use namelists, only: tag, ra, ek, pr, ek, raxi, sc, radratio, l_2D_spectra, &
       &                l_2D_SD
   use useful, only: cc2real, getMSD2, round_off
   use radial_functions, only: r, rscheme, height
   use integration, only: rInt_R
   use mean_sd, only: mean_sd_type, mean_sd_2D_type

   implicit none

   private

   type, public :: spectra_type
      type(mean_sd_2D_type) :: us2
      type(mean_sd_2D_type) :: up2
      type(mean_sd_2D_type) :: enst
      type(mean_sd_type) :: us2M
      type(mean_sd_type) :: up2M
      type(mean_sd_type) :: enstM
      integer :: n_calls
      integer :: ispec_counter
      real(cp) :: dt
      real(cp) :: timeLast
      logical :: l_calc
   contains 
      procedure :: initialize
      procedure :: finalize
      procedure :: write_spectra
      procedure :: write_2D_spectra
      procedure :: write_spectra_avg
      procedure :: calculate_spectra
   end type spectra_type

contains

   subroutine initialize(this)
      !
      ! Memory allocation and initial values
      !
      class(spectra_type) :: this

      if ( l_2D_spectra ) then
         call this%us2%initialize(nMstart,nMstop,n_r_max,l_2D_SD)
         call this%up2%initialize(nMstart,nMstop,n_r_max,l_2D_SD)
         call this%enst%initialize(nMstart,nMstop,n_r_max,l_2D_SD)
      end if

      call this%us2M%initialize(nMstart,nMstop)
      call this%up2M%initialize(nMstart,nMstop)
      call this%enstM%initialize(nMstart,nMstop)

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

      call this%enstM%finalize()
      call this%up2M%finalize()
      call this%us2M%finalize()
      if ( l_2D_spectra ) then
         call this%enst%finalize()
         call this%up2%finalize()
         call this%us2%finalize()
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
      real(cp) :: us2R(n_r_max), up2R(n_r_max), enstR(n_r_max)
      real(cp) :: sd
      integer :: n_r, n_m, m

      this%n_calls = this%n_calls+1
      this%dt = time-this%timeLast

      !-- This is not cache-friendly but hopefully it's happening only
      !-- once in a while (otherwise we need (n_r, n_m) arrays
      do n_m=nMstart,nMstop
         m = idx2m(n_m)
         do n_r=1,n_r_max
            us2R(n_r) =cc2real(us_Mloc(n_m,n_r),m)
            up2R(n_r) =cc2real(up_Mloc(n_m,n_r),m)
            enstR(n_r)=cc2real(om_Mloc(n_m,n_r),m)
            us2R(n_r) =pi*us2R(n_r)*r(n_r)*height(n_r)
            up2R(n_r) =pi*up2R(n_r)*r(n_r)*height(n_r)
            enstR(n_r)=pi*enstR(n_r)*r(n_r)*height(n_r)
            if ( l_2D_spectra ) then
               if ( this%us2%l_SD ) then
                  call getMSD2(this%us2%mean(n_m,n_r), this%us2%SD(n_m,n_r), &
                       &       us2R(n_r),this%n_calls, this%dt, time)
                  call getMSD2(this%up2%mean(n_m,n_r), this%up2%SD(n_m,n_r), &
                       &       up2R(n_r),this%n_calls, this%dt, time)
                  call getMSD2(this%enst%mean(n_m,n_r), this%enst%SD(n_m,n_r), &
                       &       enstR(n_r), this%n_calls, this%dt, time)
               else
                  call getMSD2(this%us2%mean(n_m,n_r), sd, us2R(n_r),   &
                       &       this%n_calls, this%dt, time)
                  call getMSD2(this%up2%mean(n_m,n_r), sd, up2R(n_r),   &
                       &       this%n_calls, this%dt, time)
                  call getMSD2(this%enst%mean(n_m,n_r), sd, enstR(n_r), &
                       &       this%n_calls, this%dt, time)
               end if
            end if
         end do
         !--Radial integration
         us2_m(n_m) =rInt_R(us2R, r, rscheme)
         up2_m(n_m) =rInt_R(up2R, r, rscheme)
         enst_m(n_m)=rInt_R(enstR, r, rscheme)

         !-- Mean and SD of m-spectra
         call getMSD2(this%us2M%mean(n_m), this%us2M%SD(n_m), us2_m(n_m), &
              &       this%n_calls, this%dt, time)
         call getMSD2(this%up2M%mean(n_m), this%up2M%SD(n_m), up2_m(n_m), &
              &       this%n_calls, this%dt, time)
         call getMSD2(this%enstM%mean(n_m), this%enstM%SD(n_m), enst_m(n_m), &
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
      integer :: file_handle, n_m, m, n_p, n_r
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
      call MPI_GatherV(this%us2M%mean, nm_per_rank, MPI_DEF_REAL,  &
           &           us2_mean_global, recvcounts, displs,        &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%us2M%SD, nm_per_rank, MPI_DEF_REAL,    &
           &           us2_SD_global, recvcounts, displs,          &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%up2M%mean, nm_per_rank, MPI_DEF_REAL,  &
           &           up2_mean_global, recvcounts, displs,        &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%up2M%SD, nm_per_rank, MPI_DEF_REAL,    &
           &           up2_SD_global, recvcounts, displs,          &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%enstM%mean, nm_per_rank, MPI_DEF_REAL, &
           &           enst_mean_global, recvcounts, displs,       &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_GatherV(this%enstM%SD, nm_per_rank, MPI_DEF_REAL,   &
           &           enst_SD_global, recvcounts, displs,         &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)

      !-- Only rank==0 writes the spec_avg.TAG file
      if ( rank == 0 ) then
         open(newunit=file_handle, file='spec_avg.'//tag)
         do n_m=1,n_m_max
            m = idx2m(n_m)
            us2_SD_global(n_m) =sqrt(us2_SD_global(n_m)/this%timeLast)
            up2_SD_global(n_m) =sqrt(up2_SD_global(n_m)/this%timeLast)
            enst_SD_global(n_m)=sqrt(enst_SD_global(n_m)/this%timeLast)
            write(file_handle, '(I5, 6es16.8)') m,                                &
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
         if ( this%us2%l_SD ) then
            do n_r=1,n_r_max
               do n_m=nMstart,nMstop
                  this%us2%SD(n_m,n_r) =sqrt(this%us2%SD(n_m,n_r)/this%timeLast)
                  this%up2%SD(n_m,n_r) =sqrt(this%up2%SD(n_m,n_r)/this%timeLast)
                  this%enst%SD(n_m,n_r)=sqrt(this%enst%SD(n_m,n_r)/this%timeLast)
               end do
            end do
         end if

         call this%write_2D_spectra()
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
            write(file_handle, '(I5, 3es16.8)') m,       &
            &              round_off(us2_m_global(n_m)), &
            &              round_off(up2_m_global(n_m)), &
            &              round_off(enst_m_global(n_m))
         end do
         close(file_handle)

         this%ispec_counter = this%ispec_counter+1

      end if

   end subroutine write_spectra
!----------------------------------------------------------------------
   subroutine write_2D_spectra(this)

      class(spectra_type) :: this

      !-- Local variables:
      integer :: info, fh, version, header_size
      integer :: istat(MPI_STATUS_SIZE), filetype
      integer :: arr_size(2), arr_loc_size(2), arr_start(2)
      integer(kind=MPI_OFFSET_KIND) :: disp
      character(len=100) :: spec_file

      spec_file='2D_spec_avg.'//tag

      if ( this%us2%l_SD ) then
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
      call MPI_File_Write_all(fh, this%us2%mean, nm_per_rank*n_r_max,  & 
           &                  MPI_DEF_REAL, istat, ierr)
      if ( this%us2%l_SD ) then
         call MPI_File_Write_all(fh, this%us2%SD, nm_per_rank*n_r_max, & 
              &                  MPI_DEF_REAL, istat, ierr)
      end if
      call MPI_File_Write_all(fh, this%up2%mean, nm_per_rank*n_r_max,  &
           &                  MPI_DEF_REAL, istat, ierr)
      if ( this%up2%l_SD ) then
         call MPI_File_Write_all(fh, this%up2%SD, nm_per_rank*n_r_max, & 
              &                  MPI_DEF_REAL, istat, ierr)
      end if
      call MPI_File_Write_all(fh, this%enst%mean, nm_per_rank*n_r_max, &
           &                  MPI_DEF_REAL, istat, ierr)
      if ( this%enst%l_SD ) then
         call MPI_File_Write_all(fh, this%enst%SD, nm_per_rank*n_r_max,& 
              &                  MPI_DEF_REAL, istat, ierr)
      end if

      call MPI_Info_free(info, ierr)
      call MPI_File_close(fh, ierr)

   end subroutine write_2D_spectra
!----------------------------------------------------------------------
end module spectra
