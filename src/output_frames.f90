module output_frames
   !
   ! This module handles the writing of the snapshot. It uses MPI-IO
   !

   use parallel_mod
   use precision_mod
   use namelists, only:  ra, pr, ek, radratio, raxi, sc
   use horizontal, only: theta
   use blocking, only: nRstart, nRstop, nr_per_rank, nMstart, nMstop, &
       &               nm_per_rank, nRstart3D, nRstop3D, nR_per_rank_3D
   use truncation, only: n_r_max, n_m_max, m_max, minc, n_phi_max, idx2m
   use truncation_3D, only: n_r_max_3D, n_theta_max, n_phi_max_3D, minc_3D, &
       &                    m_max_3D, lm_max, l_max
   use radial_functions, only: r, tcond, xicond, tcond_3D, r_3D

   implicit none

   private

   public :: write_snapshot, open_snapshot_3D, close_snapshot_3D, &
   &         write_bulk_snapshot_3D

contains

   subroutine write_snapshot(filename, time, arr_Mloc)

      !-- Input variables
      character(len=*), intent(in) :: filename
      complex(cp),      intent(in) :: arr_Mloc(n_m_max,nMstart:nMstop)
      real(cp),         intent(in) :: time

      !-- Local variables
      integer :: version, fh, filetype, info
      integer :: istat(MPI_STATUS_SIZE)
      integer :: arr_size(2), arr_loc_size(2), arr_start(2)
      integer(kind=MPI_OFFSET_KIND) :: disp, offset

      version = 2

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

      call MPI_File_Open(MPI_COMM_WORLD, filename, ior(MPI_MODE_WRONLY, &
           &             MPI_MODE_CREATE), info, fh, ierr)

      !-- Only rank=0 writes the header of the file
      if ( rank == 0 ) then
         call MPI_File_Write(fh, version, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, time, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, ra, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, ek, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, pr, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, radratio, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, sc, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, raxi, 1, MPI_DEF_REAL, istat, ierr)

         call MPI_File_Write(fh, n_r_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, n_m_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, m_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, minc, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, n_phi_max, 1, MPI_INTEGER, istat, ierr)

         call MPI_File_Write(fh, r, n_r_max, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, tcond, n_r_max, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, xicond, n_r_max, MPI_DEF_REAL, istat, ierr)

         call MPI_File_Write(fh, idx2m, n_m_max, MPI_INTEGER, istat, ierr)
      end if

      if ( rank == 0 ) then
         !-- Rank 0 gets the displacement
         call MPI_File_get_position(fh, offset, ierr)
         call MPI_File_get_byte_offset(fh, offset, disp, ierr)
      end if
      !-- Broadcast the displacement
      call MPI_Bcast(disp, 1, MPI_OFFSET, 0, MPI_COMM_WORLD, ierr)

      arr_size(1) = n_m_max
      arr_size(2) = n_r_max
      arr_loc_size(1) = nm_per_rank
      arr_loc_size(2) = n_r_max
      arr_start(1) = nMstart-1
      arr_start(2) = 0
      call MPI_Type_Create_Subarray(2,arr_size, arr_loc_size, arr_start, &
           &                        MPI_ORDER_FORTRAN, MPI_DEF_COMPLEX,  &
           &                        filetype, ierr)
      call MPI_Type_Commit(filetype, ierr)

      !-- Set the view after the header
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype, "native", &
           &                 info, ierr)
  
      call MPI_File_Write_all(fh, arr_Mloc, nm_per_rank*n_r_max, &
           &                  MPI_DEF_COMPLEX, istat, ierr)

      call MPI_Type_free(filetype, ierr)
      call MPI_Info_free(info, ierr)
      call MPI_File_close(fh, ierr)

   end subroutine write_snapshot
!------------------------------------------------------------------------------
   subroutine open_snapshot_3D(filename, time, fh, info)

      !-- Input variables
      character(len=100), intent(in) :: filename
      real(cp),           intent(in) :: time

      !-- Output variables
      integer, intent(out) :: fh   ! file handler
      integer, intent(out) :: info ! file info handler

      !-- Local variables
      integer :: istat(MPI_STATUS_SIZE)
      integer :: version
      integer(kind=MPI_OFFSET_KIND) :: offset, disp

      version = 1

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
      call MPI_File_Open(MPI_COMM_WORLD, filename, ior(MPI_MODE_WRONLY, &
           &             MPI_MODE_CREATE), info, fh, ierr)

      !-- Only rank=0 writes the header of the graphic file
      if ( rank == 0 ) then
         call MPI_File_Write(fh, version, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, time, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, ra, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, ek, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, pr, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, radratio, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, sc, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, raxi, 1, MPI_DEF_REAL, istat, ierr)

         call MPI_File_Write(fh, n_r_max_3D, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, l_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, m_max_3D, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, lm_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, minc_3D, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, n_theta_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, n_phi_max_3D, 1, MPI_INTEGER, istat, ierr)

         call MPI_File_Write(fh, r_3D, n_r_max_3D, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, tcond_3D, n_r_max_3D, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, theta, n_theta_max, MPI_DEF_REAL, istat, ierr)
      end if

      if ( rank == 0 ) then
         !-- Rank 0 gets the displacement
         call MPI_File_get_position(fh, offset, ierr)
         call MPI_File_get_byte_offset(fh, offset, disp, ierr)
      end if
      !-- Broadcast the displacement
      call MPI_Bcast(disp, 1, MPI_OFFSET, 0, MPI_COMM_WORLD, ierr)

      disp = disp+(nRstart3D-1)*n_phi_max_3D*n_theta_max*4

      call MPI_File_Set_View(fh, disp, MPI_REAL4, MPI_REAL4, "native", &
           &                 info, ierr)

   end subroutine open_snapshot_3D
!------------------------------------------------------------------------------
   subroutine write_bulk_snapshot_3D(fh, arr)
      !
      ! This routine writes a 3-D file radius by radius (instead of dumping
      ! everything as in write_snapshot_3D)
      !

      !-- Input variables
      integer,  intent(in) :: fh   ! file handler
      real(cp), intent(in) :: arr(n_phi_max_3D,n_theta_max)

      !-- local variables
      integer :: n_phi, n_th, counter
      integer :: istat(MPI_STATUS_SIZE)
      integer(kind=MPI_OFFSET_KIND) :: offset
      real(outp) :: dummy(n_phi_max_3D, n_theta_max)

      !-- Copy local array to a single precision file
      do n_th=1,n_theta_max
         do n_phi=1,n_phi_max_3D
            dummy(n_phi, n_th)=real(arr(n_phi,n_th),outp)
         end do
      end do

      !-- Make sure everything is written by looping and counting
      counter = 0
      do while ( n_phi_max_3D*n_theta_max /= counter ) 
         offset =  -counter*4
         if (counter /= 0 ) call MPI_File_seek(fh, offset, MPI_SEEK_CUR, ierr)
         call MPI_File_Write(fh, dummy, n_phi_max_3D*n_theta_max, &
              &              MPI_REAL4, istat, ierr)
         call MPI_Get_count(istat, MPI_REAL4, counter, ierr)
      end do

   end subroutine write_bulk_snapshot_3D
!------------------------------------------------------------------------------
   subroutine close_snapshot_3D(fh, info)

      !-- Input variables
      integer, intent(inout) :: fh   ! file handler
      integer, intent(inout) :: info ! file info

      call MPI_Info_free(info, ierr)
      call MPI_File_close(fh, ierr)

   end subroutine close_snapshot_3D
!------------------------------------------------------------------------------
   subroutine write_snapshot_3D(filename, time, arr)

      !-- Input variables
      character(len=100), intent(in) :: filename
      real(cp),           intent(in) :: time
      real(cp), intent(in) :: arr(n_phi_max_3D,n_theta_max,nRstart3D:nRstop3D)

      !-- Local variables
      integer :: istat(MPI_STATUS_SIZE)
      integer :: info, fh, version, filetype
      integer :: arr_size(3), arr_loc_size(3), arr_start(3)
      integer(kind=MPI_OFFSET_KIND) :: disp, offset

      version = 1

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
      call MPI_File_Open(MPI_COMM_WORLD, filename, ior(MPI_MODE_WRONLY, &
           &             MPI_MODE_CREATE), info, fh, ierr)

      !-- Only rank=0 writes the header of the graphic file
      if ( rank == 0 ) then
         call MPI_File_Write(fh, version, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, time, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, ra, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, ek, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, pr, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, radratio, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, sc, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, raxi, 1, MPI_DEF_REAL, istat, ierr)

         call MPI_File_Write(fh, n_r_max_3D, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, l_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, m_max_3D, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, lm_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, minc_3D, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, n_theta_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, n_phi_max_3D, 1, MPI_INTEGER, istat, ierr)

         call MPI_File_Write(fh, r_3D, n_r_max_3D, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, tcond_3D, n_r_max_3D, MPI_DEF_REAL, istat, ierr)
      end if

      if ( rank == 0 ) then
         !-- Rank 0 gets the displacement
         call MPI_File_get_position(fh, offset, ierr)
         call MPI_File_get_byte_offset(fh, offset, disp, ierr)
      end if
      !-- Broadcast the displacement
      call MPI_Bcast(disp, 1, MPI_OFFSET, 0, MPI_COMM_WORLD, ierr)

      arr_size(1) = n_phi_max_3D
      arr_size(2) = n_theta_max
      arr_size(3) = n_r_max_3D
      arr_loc_size(1) = n_phi_max_3D
      arr_loc_size(2) = n_theta_max
      arr_loc_size(3) = nR_per_rank_3D
      arr_start(1) = 0
      arr_start(2) = 0
      arr_start(3) = nRstart3D-1
      call MPI_Type_Create_Subarray(3,arr_size, arr_loc_size, arr_start, &
           &                        MPI_ORDER_FORTRAN, MPI_REAL4,        &
           &                        filetype, ierr)
      call MPI_Type_Commit(filetype, ierr)

      !-- Set the view after the header
      call MPI_File_Set_View(fh, disp, MPI_REAL4, filetype, "native", &
           &                 info, ierr)

      call MPI_File_Write_all(fh, real(arr,kind=outp), n_phi_max_3D* &
           &                  n_theta_max*nR_per_rank_3D, MPI_REAL4, istat, ierr)

      call MPI_Type_free(filetype, ierr)
      call MPI_Info_free(info, ierr)
      call MPI_File_close(fh, ierr)

   end subroutine write_snapshot_3D
!------------------------------------------------------------------------------
end module output_frames
