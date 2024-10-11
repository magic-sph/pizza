module energy_transfer

   use precision_mod
   use parallel_mod
   use namelists, only: ra, pr, ek, raxi, sc, radratio, tag
   use fields, only: us_Rloc, up_Rloc
   use mem_alloc, only: bytes_allocated
   use constants, only: one, two, ci, pi
   use communications, only: m2r_single
   use blocking, only: nMstart, nMstop, nRstart, nRstop, getBlocks, load, &
       &               nR_per_rank, radial_balance
   use radial_functions, only: rscheme, m_R, or1, r
   use horizontal, only: phi
   use truncation, only: n_m_max, n_r_max, idx2m, n_phi_max, minc, m_max
   use radial_der, only: get_dr
   use integration, only: rInt_R
   use fourier, only: ifft

   implicit none

   private

   integer :: n_m_transfer,nMtranstart,nMtranstop,nmtran_per_rank
   integer :: n_calls ! Number of computation over the run
   type(load), allocatable :: m_tran_bal(:)
   integer, allocatable :: counts(:), disp(:), rtype(:), stype(:)
   real(cp), allocatable :: trans_dist(:,:)

   public :: calc_energy_transfer, initialize_transfer, finalize_transfer

contains

   subroutine initialize_transfer(m_max_transfer)

      !-- Input variable
      integer, intent(in) :: m_max_transfer ! Maximum m values

      !-- Local variables
      integer :: arr_size(3), arr_loc_size(3), arr_start(3)
      integer :: p, my_m_counts

      allocate ( m_tran_bal(0:n_procs-1) )
      n_m_transfer=m_max_transfer/minc+1

      call getBlocks(m_tran_bal, n_m_transfer, n_procs)
      nMtranstart = m_tran_bal(rank)%nStart
      nMtranstop = m_tran_bal(rank)%nStop
      nmtran_per_rank = m_tran_bal(rank)%n_per_rank

      allocate( counts(0:n_procs-1), disp(0:n_procs-1) )
      allocate( rtype(0:n_procs-1), stype(0:n_procs-1) )
      bytes_allocated = bytes_allocated + 4*n_procs*SIZEOF_INTEGER

      allocate( trans_dist(nMtranstart:nMtranstop,n_m_transfer) )
      bytes_allocated=bytes_allocated+n_m_transfer*nmtran_per_rank*SIZEOF_DEF_REAL
      trans_dist(:,:)=0.0_cp
      n_calls=0

      do p=0,n_procs-1
         my_m_counts = m_tran_bal(p)%n_per_rank

         counts(p)=1
         disp(p)  =0

         arr_size(1)=n_m_transfer
         arr_size(2)=n_m_transfer
         arr_size(3)=nR_per_rank
         arr_loc_size(1)=n_m_transfer
         arr_loc_size(2)=my_m_counts
         arr_loc_size(3)=nR_per_rank
         arr_start(1)=0
         arr_start(2)=m_tran_bal(p)%nStart-1
         arr_start(3)=0
         call MPI_Type_Create_Subarray(3, arr_size, arr_loc_size, arr_start, &
              &                        MPI_ORDER_FORTRAN, MPI_DEF_REAL,      &
              &                        stype(p), ierr)
         call MPI_Type_Commit(stype(p), ierr)

         arr_size(1)=n_m_transfer
         arr_size(2)=nmtran_per_rank
         arr_size(3)=n_r_max
         arr_loc_size(1)=n_m_transfer
         arr_loc_size(2)=nmtran_per_rank
         arr_loc_size(3)=radial_balance(p)%n_per_rank
         arr_start(1)=0
         arr_start(2)=0
         arr_start(3)=radial_balance(p)%nStart-1
         call MPI_Type_Create_Subarray(3, arr_size, arr_loc_size, arr_start, &
              &                        MPI_ORDER_FORTRAN, MPI_DEF_REAL,      &
              &                        rtype(p), ierr)
         call MPI_Type_Commit(rtype(p), ierr)
      end do

   end subroutine initialize_transfer
!-------------------------------------------------------------------------------
   subroutine finalize_transfer

      deallocate( counts, disp, rtype, stype, m_tran_bal, trans_dist )

   end subroutine finalize_transfer
!-------------------------------------------------------------------------------
   subroutine calc_energy_transfer(us_Mloc, up_Mloc, l_stop_time)

      !-- Input arrays
      complex(cp), intent(inout) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: up_Mloc(nMstart:nMstop,n_r_max)
      logical,     intent(in) :: l_stop_time

      !-- Local variables
      integer :: n_m, m, n_r, n_phi, n_m1, n_m2
      real(cp) :: fac, dusds, dupds, dusdp, dupdp, dphi
      real(cp) :: tmpR(n_r_max)
      real(cp) :: us_grid(n_phi_max), up_grid(n_phi_max)
      real(cp) :: tmp(n_m_transfer,n_m_transfer)
      real(cp) :: trans_Rloc(n_m_transfer,n_m_transfer,nRstart:nRstop)
      real(cp) :: trans_Mloc(n_m_transfer,nMtranstart:nMtranstop,n_r_max)
      real(cp) :: us_filt(n_m_transfer,n_phi_max,nRstart:nRstop)
      real(cp) :: up_filt(n_m_transfer,n_phi_max,nRstart:nRstop)
      real(cp) :: ugradus_filt(n_m_transfer,n_phi_max,nRstart:nRstop)
      real(cp) :: ugradup_filt(n_m_transfer,n_phi_max,nRstart:nRstop)
      complex(cp) :: dus_Mloc(nMstart:nMstop,n_r_max),dup_Mloc(nMstart:nMstop,n_r_max)
      complex(cp) :: dus_Rloc(n_m_max,nRstart:nRstop),dup_Rloc(n_m_max,nRstart:nRstop)

      !-- Get the radial derivatives
      call get_dr(us_Mloc,dus_Mloc,nMstart,nMstop,n_r_max,rscheme,nocopy=.false.,&
           &      l_dct_in=.true.)
      call get_dr(up_Mloc,dup_Mloc,nMstart,nMstop,n_r_max,rscheme,nocopy=.false.,&
           &      l_dct_in=.true.)

      !-- Transpose them in the r-distributed configuration
      call m2r_single%transp_m2r(dus_Mloc,dus_Rloc)
      call m2r_single%transp_m2r(dup_Mloc,dup_Rloc)

      !-- First construct help arrays
      do n_r=nRstart,nRstop
         call ifft(us_Rloc(:,n_r), us_grid, m_R(n_r))
         call ifft(up_Rloc(:,n_r), up_grid, m_R(n_r))

         do n_phi=1,n_phi_max
            do n_m=1,n_m_transfer
               m=idx2m(n_m)
               if ( m == 0 ) then
                  fac=one
               else
                  fac=two
               end if

               !-- Flow components filtered at degree m expressed on the grid
               us_filt(n_m,n_phi,n_r)=fac*real(exp(ci*m*phi(n_phi))* &
               &                               us_Rloc(n_m,n_r))
               up_filt(n_m,n_phi,n_r)=fac*real(exp(ci*m*phi(n_phi))* &
               &                               up_Rloc(n_m,n_r))

               !-- radial derivatives
               dusds=fac*real(exp(ci*m*phi(n_phi))*dus_Rloc(n_m,n_r))
               dupds=fac*real(exp(ci*m*phi(n_phi))*dup_Rloc(n_m,n_r))

               !-- phi derivatives
               dusdp=fac*real(ci*m*exp(ci*m*phi(n_phi))*us_Rloc(n_m,n_r))
               dupdp=fac*real(ci*m*exp(ci*m*phi(n_phi))*up_Rloc(n_m,n_r))

               !-- radial component of the advection term filtered at m
               ugradus_filt(n_m,n_phi,n_r)=us_grid(n_phi)*dusds           + &
               &                           up_grid(n_phi)*or1(n_r)*(dusdp - &
               &                           up_filt(n_m,n_phi,n_r))

               !-- phi component of the advection term
               ugradup_filt(n_m,n_phi,n_r)=us_grid(n_phi)*dupds           + &
               &                           up_grid(n_phi)*or1(n_r)*(dupdp + &
               &                           us_filt(n_m,n_phi,n_r))
            end do
         end do
      end do

      !-- Now compute the transfer function
      dphi=two*pi/n_phi_max
      do n_r=nRstart,nRstop
         tmp(:,:)=0.0_cp
         do n_phi=1,n_phi_max
            do n_m=1,n_m_transfer
                  tmp(:,n_m)=tmp(:,n_m)+                                      &
                  &          us_filt(n_m,n_phi,n_r)*ugradus_filt(:,n_phi,n_r)+&
                  &          up_filt(n_m,n_phi,n_r)*ugradup_filt(:,n_phi,n_r)
            end do
         end do
         trans_Rloc(:,:,n_r)=tmp(:,:)*dphi
      end do

      !-- Transpose the data to handle the radial integration
      call transp_Rloc_Mloc(trans_Rloc, trans_Mloc)

      !-- Compute the radial integration
      do n_m1=1,n_m_transfer
         do n_m2=nMtranstart,nMtranstop
            tmpR=trans_Mloc(n_m1,n_m2,:)
            !-- Sum the snapshots to stack them
            trans_dist(n_m2,n_m1)=trans_dist(n_m2,n_m1)+rInt_R(tmpR*r,r,rscheme)
         end do
      end do
      n_calls=n_calls+1

      if ( l_stop_time ) call write_output('energy_transfer.'//tag)

   end subroutine calc_energy_transfer
!-------------------------------------------------------------------------------
   subroutine write_output(filename)

      !-- Input variable
      character(len=*), intent(in) :: filename

      !-- Local variables
      integer :: version, header_size, fh, filetype, info
      integer :: istat(MPI_STATUS_SIZE)
      integer :: arr_size(2), arr_loc_size(2), arr_start(2)
      integer(kind=MPI_OFFSET_KIND) :: disp

      trans_dist(:,:)=trans_dist(:,:)/n_calls

      version=1
      header_size=SIZEOF_INTEGER+6*SIZEOF_DEF_REAL+6*SIZEOF_INTEGER+ &
      &            n_m_max*SIZEOF_INTEGER

      !-- MPI-IO setup
      call mpiio_setup(info)

      !-- Open file
      call MPI_File_Open(MPI_COMM_WORLD, filename, ior(MPI_MODE_WRONLY, &
           &             MPI_MODE_CREATE), info, fh, ierr)

      !-- Only rank=0 writes the header of the file
      if ( rank == 0 ) then
         call MPI_File_Write(fh, version, 1, MPI_INTEGER, istat, ierr)
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
         call MPI_File_Write(fh, n_m_transfer, 1, MPI_INTEGER, istat, ierr)

         call MPI_File_Write(fh, idx2m, n_m_max, MPI_INTEGER, istat, ierr)
      end if

      arr_size(1) = n_m_transfer
      arr_size(2) = n_m_transfer
      arr_loc_size(1) = nmtran_per_rank
      arr_loc_size(2) = n_m_transfer
      arr_start(1) = nMtranstart-1
      arr_start(2) = 0
      call MPI_Type_Create_Subarray(2, arr_size, arr_loc_size, arr_start, &
           &                        MPI_ORDER_FORTRAN, MPI_DEF_REAL,      &
           &                        filetype, ierr)
      call MPI_Type_Commit(filetype, ierr)

      !-- Set the view after the header
      disp = header_size
      call MPI_File_Set_View(fh, disp, MPI_DEF_REAL, filetype, "native", &
           &                 info, ierr)

      call MPI_File_Write_all(fh, trans_dist, nmtran_per_rank*n_m_transfer, &
           &                  MPI_DEF_REAL, istat, ierr)

      call MPI_Info_free(info, ierr)
      call MPI_File_close(fh, ierr)

   end subroutine write_output
!-------------------------------------------------------------------------------
   subroutine transp_Rloc_Mloc(arr_Rloc, arr_Mloc)

      real(cp), intent(in) :: arr_Rloc(n_m_transfer,n_m_max,nRstart:nRstop)
      real(cp), intent(out) :: arr_Mloc(n_m_transfer,nMtranstart:nMtranstop,n_r_max)

      call MPI_Alltoallw(arr_Rloc, counts, disp, stype, &
           &             arr_Mloc, counts, disp, rtype, &
           &             MPI_COMM_WORLD, ierr)

   end subroutine transp_Rloc_Mloc
!-------------------------------------------------------------------------------
end module energy_transfer
