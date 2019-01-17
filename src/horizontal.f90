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
   use namelists, only: hdif_m, hdif_exp, hdif_vel, hdif_temp, tag,   &
       &                t_bot, t_top, xi_bot, xi_top, hdif_comp,      &
       &                l_heat, l_chem

   implicit none

   private

   !-- Arrays for Inhomogeneous temperature B.Cs
   complex(cp), allocatable, public :: bott_Mloc(:)
   complex(cp), allocatable, public :: topt_Mloc(:)
   complex(cp), allocatable, public :: botxi_Mloc(:)
   complex(cp), allocatable, public :: topxi_Mloc(:)

   real(cp), public, allocatable :: hdif_V(:)
   real(cp), public, allocatable :: hdif_T(:)
   real(cp), public, allocatable :: hdif_Xi(:)

   public :: initialize_mfunctions, finalize_mfunctions, mfunctions

contains

   subroutine initialize_mfunctions

      allocate( hdif_V(nMstart:nMstop) )
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_DEF_REAL
      if ( l_heat ) then
         allocate( bott_Mloc(nMStart:nMstop) )
         allocate( topt_Mloc(nMStart:nMstop) )
         allocate( hdif_T(nMstart:nMstop) )
         bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_DEF_REAL+&
         &                 2*(nMstop-nMstart+1)*SIZEOF_DEF_COMPLEX
      end if
      if ( l_chem ) then
         allocate( botxi_Mloc(nMStart:nMstop) )
         allocate( topxi_Mloc(nMStart:nMstop) )
         allocate( hdif_Xi(nMstart:nMstop) )
         bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_DEF_REAL+&
         &                 2*(nMstop-nMstart+1)*SIZEOF_DEF_COMPLEX
      end if

   end subroutine initialize_mfunctions
!--------------------------------------------------------------------------------
   subroutine finalize_mfunctions

      if ( l_heat ) deallocate( bott_Mloc, topt_Mloc, hdif_T )
      if ( l_chem ) deallocate( botxi_Mloc, topxi_Mloc, hdif_Xi )
      deallocate( hdif_V )

   end subroutine finalize_mfunctions
!--------------------------------------------------------------------------------
   subroutine mfunctions

      !-- Local variables
      integer :: n, m_bot, m_top
      integer :: n_m, m, n_p, file_handle
      integer :: displs(0:n_procs-1), recvcounts(0:n_procs-1)
      real(cp) :: eps, tr_bot, ti_bot, tr_top, ti_top
      real(cp) :: hdif_T_global(n_m_max), hdif_V_global(n_m_max)
      real(cp) :: hdif_Xi_global(n_m_max)

      eps = 10.0_cp*epsilon(one)
      if ( abs(hdif_temp) > eps .or. abs(hdif_vel) > eps .or. &
         & abs(hdif_comp) > eps ) then

         do n_m=nMstart,nMstop
            m = idx2m(n_m)

            if ( m > hdif_m ) then ! This is the Nataf & Schaffer (2015) form
               if ( l_heat ) hdif_T(n_m) = hdif_temp**(m-hdif_m)
               if ( l_chem ) hdif_Xi(n_m) = hdif_comp**(m-hdif_m)
               hdif_V(n_m) = hdif_vel**(m-hdif_m)
            else
               if ( l_heat ) hdif_T(n_m) = one
               if ( l_chem ) hdif_Xi(n_m) = one
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
         if ( l_heat ) then
            call MPI_GatherV(hdif_T, nm_per_rank, MPI_DEF_REAL,  &
                 &           hdif_T_global, recvcounts, displs,  &
                 &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
         end if
         if ( l_chem ) then
            call MPI_GatherV(hdif_Xi, nm_per_rank, MPI_DEF_REAL,  &
                 &           hdif_Xi_global, recvcounts, displs,  &
                 &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
         end if
         call MPI_GatherV(hdif_V, nm_per_rank, MPI_DEF_REAL,  &
              &           hdif_V_global, recvcounts, displs,  &
              &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)


         !-- Now rank0 writes an output file
         if ( rank== 0 ) then
            open(newunit=file_handle, file='hdif.'//tag, status='new')

            do n_m=1,n_m_max
               m = idx2m(n_m)
               write(file_handle, '(I5, 3es16.8)') m, hdif_T_global(n_m), &
               &                                      hdif_Xi_global(n_m),&
               &                                      hdif_V_global(n_m)
            end do
            close(file_handle)
         end if

      else

         hdif_V(:) = one
         if ( l_heat ) hdif_T(:) = one
         if ( l_chem ) hdif_Xi(:) = one

      end if

      !-- Build matrices Inhomogeneous B.Cs
      !-- Inspired by MagIC (J.Wicht; T.Gastine; A.Barik; R.Raynaud; B.Putigny;
      !--                    R.Yadav; ; L.Duarte; B.dintrans; T.Schwaiger)!
      if ( l_heat ) then
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
      end if

      if ( l_chem ) then
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            botxi_Mloc(n_m)=zero
            topxi_Mloc(n_m)=zero
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
                  botxi_Mloc(n_m)=cmplx(tr_bot,ti_bot,kind=cp)
               end if
               if ( m_top == m .and. &
                   cmplx(tr_top,ti_top,kind=cp) /= zero ) then
                  if ( m == 0 ) ti_top=0.0_cp
                  topxi_Mloc(n_m)=cmplx(tr_top,ti_top,kind=cp)
               end if
            end do
         end do

      end if

   end subroutine mfunctions
!--------------------------------------------------------------------------------
end module horizontal
