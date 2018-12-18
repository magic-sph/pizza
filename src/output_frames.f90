module output_frames
   !
   ! This module handles the writing of the snapshot. It uses MPI-IO
   !

   use parallel_mod
   use precision_mod
   use namelists, only:  ra, pr, ek, radratio, raxi, sc
   use blocking, only: nRstart, nRstop, nr_per_rank, nMstart, nMstop, &
       &               nm_per_rank
   use truncation, only: n_r_max, n_m_max, m_max, minc, n_phi_max, idx2m
   use radial_functions, only: r, tcond, xicond

   implicit none

   private

   public :: write_snapshot_rloc, write_snapshot_mloc

contains

   subroutine write_snapshot_rloc(filename, time, arr_Rloc)

      !-- Input variables
      character(len=*), intent(in) :: filename
      complex(cp),      intent(in) :: arr_Rloc(nRstart:nRstop,n_r_max)
      real(cp),         intent(in) :: time

      !-- Local variables
      integer :: version, header_size, fh, filetype
      integer :: istat(MPI_STATUS_SIZE)
      integer(kind=MPI_OFFSET_KIND) :: disp

      version = 1

      header_size = SIZEOF_INTEGER+8+SIZEOF_DEF_REAL+8+6*SIZEOF_DEF_REAL+8+&
      &             5*SIZEOF_INTEGER+8+n_r_max*SIZEOF_DEF_REAL+8+n_r_max*  &
      &             SIZEOF_DEF_REAL+8+n_m_max*SIZEOF_INTEGER

      call MPI_File_Open(MPI_COMM_WORLD, filename, ior(MPI_MODE_WRONLY, &
           &             MPI_MODE_CREATE), MPI_INFO_NULL, fh, ierr)

      if ( rank == 0 ) then
         disp = 0
      else
         disp = 8+header_size
      end if
      call MPI_File_Set_View(fh, disp, MPI_BYTE, MPI_BYTE, "native", &
           &                 MPI_INFO_NULL, ierr)

      !-- Only rank=0 writes the header of the file
      if ( rank == 0 ) then
         call MPI_File_Write(fh, SIZEOF_INTEGER, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, version, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, SIZEOF_INTEGER, 1, MPI_INTEGER, istat, ierr)

         call MPI_File_Write(fh, SIZEOF_DEF_REAL, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, time, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, SIZEOF_DEF_REAL, 1, MPI_INTEGER, istat, ierr)

         call MPI_File_Write(fh, 6*SIZEOF_DEF_REAL, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, ra, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, ek, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, pr, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, radratio, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, sc, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, raxi, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, 6*SIZEOF_DEF_REAL, 1, MPI_INTEGER, istat, ierr)

         call MPI_File_Write(fh, 5*SIZEOF_INTEGER, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, n_r_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, n_m_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, m_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, minc, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, n_phi_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, 5*SIZEOF_INTEGER, 1, MPI_INTEGER, istat, ierr)

         call MPI_File_Write(fh, n_r_max*SIZEOF_DEF_REAL, 1, MPI_INTEGER, &
              &              istat, ierr)
         call MPI_File_Write(fh, r, n_r_max, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, n_r_max*SIZEOF_DEF_REAL, 1, MPI_INTEGER, &
              &              istat, ierr)

         call MPI_File_Write(fh, n_r_max*SIZEOF_DEF_REAL, 1, MPI_INTEGER, &
              &              istat, ierr)
         call MPI_File_Write(fh, tcond, n_r_max, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, n_r_max*SIZEOF_DEF_REAL, 1, MPI_INTEGER, &
              &              istat, ierr)

         call MPI_File_Write(fh, n_m_max*SIZEOF_INTEGER, 1, MPI_INTEGER, &
              &              istat, ierr)
         call MPI_File_Write(fh, idx2m, n_m_max, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, n_m_max*SIZEOF_INTEGER, 1, MPI_INTEGER, &
              &              istat, ierr)

         !-- Record marker for the field
         call MPI_File_Write(fh, n_m_max*n_r_max*SIZEOF_DEF_COMPLEX, 1, &
              &              MPI_INTEGER, istat, ierr)
      end if

      call MPI_Type_Vector(1,n_m_max*nR_per_rank, n_m_max*n_r_max,  &
           &               MPI_DEF_COMPLEX, filetype, ierr)
      call MPI_Type_Commit(filetype, ierr)

      !-- Set the view after the header
      disp = 8+header_size+4+rank*nR_per_rank*n_m_max*SIZEOF_DEF_COMPLEX
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype, "native", &
           &                  MPI_INFO_NULL, ierr)
  
      call MPI_File_Write(fh, arr_Rloc, nR_per_rank*n_m_max, MPI_DEF_COMPLEX, &
           &              istat, ierr)

      !-- Record marker
      disp = 8+header_size+4+n_m_max*n_r_max*SIZEOF_DEF_COMPLEX
      call MPI_File_Set_View(fh, disp, MPI_BYTE, MPI_BYTE,"native", &
           &                  MPI_INFO_NULL, ierr)
      if ( rank == n_procs-1 ) then
         call MPI_File_Write(fh, n_m_max*n_r_max*SIZEOF_DEF_COMPLEX,1, &
              &              MPI_INTEGER, istat, ierr)
      end if

      call MPI_File_close(fh, ierr)

   end subroutine write_snapshot_rloc
!------------------------------------------------------------------------------
   subroutine write_snapshot_mloc(filename, time, arr_Mloc)

      !-- Input variables
      character(len=*), intent(in) :: filename
      complex(cp),      intent(in) :: arr_Mloc(n_m_max,nMstart:nMstop)
      real(cp),         intent(in) :: time

      !-- Local variables
      integer :: version, header_size, fh, filetype, info
      integer :: istat(MPI_STATUS_SIZE)
      integer :: arr_size(2), arr_loc_size(2), arr_start(2)
      integer(kind=MPI_OFFSET_KIND) :: disp

      version = 2

      header_size = SIZEOF_INTEGER+7*SIZEOF_DEF_REAL+5*SIZEOF_INTEGER+ &
      &             3*n_r_max*SIZEOF_DEF_REAL+n_m_max*SIZEOF_INTEGER

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
      disp = header_size
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype, "native", &
           &                 info, ierr)
  
      call MPI_File_Write_all(fh, arr_Mloc, nm_per_rank*n_r_max, &
           &                  MPI_DEF_COMPLEX, istat, ierr)

      call MPI_Info_free(info, ierr)
      call MPI_File_close(fh, ierr)

   end subroutine write_snapshot_mloc
!------------------------------------------------------------------------------
end module output_frames
