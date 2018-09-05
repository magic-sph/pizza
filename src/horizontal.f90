module horizontal
   !
   ! This module implements the calculation of m-dependent array
   ! such as hyperdiffusion or heat flux patterns
   !

   use precision_mod
   use parallel_mod
   use mem_alloc
   use constants, only: one
   use blocking, only: nMstart, nMstop, m_balance, nm_per_rank
   use truncation, only: idx2m, n_m_max
   use namelists, only: hdif_m, hdif_exp, hdif_vel, hdif_temp, tag

   implicit none

   private

   real(cp), public, allocatable :: hdif_V(:)
   real(cp), public, allocatable :: hdif_T(:)

   public :: initialize_mfunctions, finalize_mfunctions, mfunctions

contains

   subroutine initialize_mfunctions

      allocate( hdif_V(nMstart:nMstop) )
      allocate( hdif_T(nMstart:nMstop) )
      bytes_allocated = bytes_allocated+2*(nMstop-nMstart+1)*SIZEOF_DEF_REAL

   end subroutine initialize_mfunctions
!--------------------------------------------------------------------------------
   subroutine finalize_mfunctions

      deallocate( hdif_V, hdif_T )

   end subroutine finalize_mfunctions
!--------------------------------------------------------------------------------
   subroutine mfunctions

      !-- Local variables
      integer :: n_m, m, n_p, file_handle
      integer :: displs(0:n_procs-1), recvcounts(0:n_procs-1)
      real(cp) :: eps
      real(cp) :: hdif_T_global(n_m_max), hdif_V_global(n_m_max)

      eps = 10.0_cp*epsilon(one)
      if ( abs(hdif_temp) > eps .or. abs(hdif_vel) > eps ) then

         do n_m=nMstart,nMstop
            m = idx2m(n_m)

            if ( m > hdif_m ) then ! This is the Nataf & Schaffer (2015) form
               hdif_T(n_m) = hdif_temp**(m-hdif_m)
               hdif_V(n_m) = hdif_vel**(m-hdif_m)
            else
               hdif_T(n_m) = one
               hdif_V(n_m) = one
            end if

         end do

         !-- Gather the profiles on rank 0 to write the profiles in hdif.TAG
         do n_p=0,n_procs-1
            recvcounts(n_p)=m_balance(n_p)%n_per_rank
         end do
         displs(0)=0
         do n_p=1,n_procs-1
            displs(n_p)=displs(n_p-1)+recvcounts(n_p-1)
         end do
         call MPI_GatherV(hdif_T, nm_per_rank, MPI_DEF_REAL,  &
              &           hdif_T_global, recvcounts, displs,  &
              &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
         call MPI_GatherV(hdif_V, nm_per_rank, MPI_DEF_REAL,  &
              &           hdif_V_global, recvcounts, displs,  &
              &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)


         !-- Now rank0 writes an output file
         if ( rank== 0 ) then
            open(newunit=file_handle, file='hdif.'//tag, status='new')

            do n_m=1,n_m_max
               m = idx2m(n_m)
               write(file_handle, '(I4, 2es16.8)') m, hdif_T_global(n_m), &
               &                                      hdif_V_global(n_m)
            end do
            close(file_handle)
         end if

      else

         hdif_V(:) = one
         hdif_T(:) = one

      end if

   end subroutine mfunctions
!--------------------------------------------------------------------------------
end module horizontal
