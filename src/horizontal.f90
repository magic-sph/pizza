module horizontal
   !
   ! This module implements the calculation of m-dependent array
   ! such as hyperdiffusion or heat flux patterns
   !

   use precision_mod
   use parallel_mod
   use mem_alloc
   use constants, only: zero, one
   use blocking, only: nMstart, nMstop, m_balance, nm_per_rank
   use truncation, only: idx2m, n_m_max
   use namelists, only: hdif_m, hdif_exp, hdif_vel, hdif_temp, tag, &
       &                t_bot,t_top

   implicit none

   private

   !-- Arrays for Inhomogeneous temperature B.Cs
   complex(cp), allocatable, public :: bott_Mloc(:)
   complex(cp), allocatable, public :: topt_Mloc(:)

   real(cp), public, allocatable :: hdif_V(:)
   real(cp), public, allocatable :: hdif_T(:)

   public :: initialize_mfunctions, finalize_mfunctions, mfunctions

contains

   subroutine initialize_mfunctions

      allocate( bott_Mloc(nMStart:nMstop) )
      allocate( topt_Mloc(nMStart:nMstop) )
      bytes_allocated = bytes_allocated + &
      &                 2*(nMstop-nMstart+1)*SIZEOF_DEF_COMPLEX

      topt_Mloc(:)=zero
      bott_Mloc(:)=zero

      allocate( hdif_V(nMstart:nMstop) )
      allocate( hdif_T(nMstart:nMstop) )
      bytes_allocated = bytes_allocated+2*(nMstop-nMstart+1)*SIZEOF_DEF_REAL

   end subroutine initialize_mfunctions
!--------------------------------------------------------------------------------
   subroutine finalize_mfunctions

      deallocate( bott_Mloc, topt_Mloc )
      deallocate( hdif_V, hdif_T )

   end subroutine finalize_mfunctions
!--------------------------------------------------------------------------------
   subroutine mfunctions

      !-- Local variables
      integer :: n, m_bot, m_top
      integer :: n_m, m, n_p, file_handle
      integer :: displs(0:n_procs-1), recvcounts(0:n_procs-1)
      real(cp) :: eps, tr_bot, ti_bot, tr_top, ti_top
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

      !-- Build matrices Inhomogeneous B.Cs
      !-- Inspired by MagIC (J.Wicht; T.Gastine; A.Barik; R.Raynaud; B.Putigny;
      !--                    R.Yadav; ; L.Duarte; B.dintrans; T.Schwaiger)!
      do n_m=nMstart,nMstop
         m = idx2m(n_m)
         bott_Mloc(n_m)=zero
         topt_Mloc(n_m)=zero
         do n=1,size(t_bot)/3
            m_bot =int(t_bot(3*n-2))
            tr_bot=t_bot(3*n-1)
            ti_bot=t_bot(3*n)
            m_top =int(t_top(3*n-2))
            tr_top=t_top(3*n-1)
            ti_top=t_top(3*n)
            if ( m_bot == m .and. &
                cmplx(tr_bot,ti_bot,kind=cp) /= zero ) then
               if ( m == 0 ) ti_bot=0.0_cp
               bott_Mloc(n_m)=cmplx(tr_bot,ti_bot,kind=cp)
            end if
            if ( m_top == m .and. &
                cmplx(tr_top,ti_top,kind=cp) /= zero ) then
               if ( m == 0 ) ti_top=0.0_cp
               topt_Mloc(n_m)=cmplx(tr_top,ti_top,kind=cp)
            end if
         end do
      end do

   end subroutine mfunctions
!--------------------------------------------------------------------------------
end module horizontal
