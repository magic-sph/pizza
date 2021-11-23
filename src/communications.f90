module communications

   use iso_fortran_env, only: output_unit
   use mpimod
   use precision_mod
   use blocking
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, n_m_max
   use namelists, only: mpi_transp, l_packed_transp, l_heat, l_chem
   use parallel_mod, only: n_procs, rank, ierr
   use mpi_alltoall_mod, only: type_mpiatoav, type_mpiatoaw
   use mpi_transp_mod, only: type_mpitransp
   use char_manip, only: capitalize

   implicit none

   private

   type, public :: help_transp
      integer,     allocatable :: rcounts(:)
      integer,     allocatable :: scounts(:)
      integer,     allocatable :: rdisp(:)
      integer,     allocatable :: sdisp(:)
      complex(cp), allocatable :: rbuff(:)
      complex(cp), allocatable :: sbuff(:)
      integer :: max_send, max_recv
   end type

   class(type_mpitransp), public, pointer :: r2m_fields, r2m_single
   class(type_mpitransp), public, pointer :: m2r_fields, m2r_single

   public :: initialize_communications, gather_from_Rloc,           &
   &         gather_from_mloc_to_rank0, scatter_from_rank0_to_mloc, &
   &         finalize_communications, reduce_radial_on_rank,        &
   &         my_reduce_mean, my_allreduce_maxloc,                   &
   &         scatter_from_rank0_to_Rloc

contains

   subroutine initialize_communications(n_log_file)

      !-- Input variable
      integer, intent(in) :: n_log_file

      !-- Local variables
      integer :: idx, n, n_out

      call capitalize(mpi_transp)
      if ( index(mpi_transp, 'AUTO') /= 0 ) then
         call find_faster_comm(idx,1,n_log_file)
      else if ( index(mpi_transp, 'ATOAV') /= 0 .or. index(mpi_transp, 'A2AV') /=0&
      &         .or. index(mpi_transp, 'ALLTOALLV') /= 0 .or. &
      &         index(mpi_transp, 'ALL2ALLV') /= 0 .or. &
      &         index(mpi_transp, 'ALL-TO-ALLV') /= 0 ) then
         idx = 1
      else if ( index(mpi_transp, 'ATOAW') /= 0 .or. index(mpi_transp, 'A2AW') /=0&
      &         .or. index(mpi_transp, 'ALLTOALLW') /= 0 .or. &
      &         index(mpi_transp, 'ALL2ALLW') /= 0 .or. &
      &         index(mpi_transp, 'ALL-TO-ALLW') /= 0 ) then
         idx = 2
      end if

      !-- Write choosen MPI communicator
      if ( rank == 0 ) then
         do n=1,2
            if ( n==1 ) then
               n_out = n_log_file
            else
               n_out = output_unit
            end if
            if ( idx == 1 ) then
               write(n_out,*) '! -> I choose alltoallv'
            else if ( idx == 2 ) then
               write(n_out,*) '! -> I choose alltoallw'
            end if
            !write(n_out,*)
         end do
      end if

      if ( idx == 1 ) then
         allocate(type_mpiatoav :: r2m_fields)
         allocate(type_mpiatoav :: r2m_single)
         allocate(type_mpiatoav :: m2r_fields)
         allocate(type_mpiatoav :: m2r_single)
      else if (idx == 2 ) then
         allocate(type_mpiatoaw :: r2m_fields)
         allocate(type_mpiatoaw :: r2m_single)
         allocate(type_mpiatoaw :: m2r_fields)
         allocate(type_mpiatoav :: m2r_single)
      end if

      if ( l_packed_transp ) then
         call r2m_single%create_comm(1)
         if ( l_heat ) then
            if ( l_chem ) then ! Heat and composition
               call m2r_fields%create_comm(5)
               call r2m_fields%create_comm(6)
            else ! Heat and no composition
               call m2r_fields%create_comm(4)
               call r2m_fields%create_comm(4)
            end if
         else
            if ( l_chem ) then ! No heat and composition
               call m2r_fields%create_comm(4)
               call r2m_fields%create_comm(4)
            else ! No heat, no composition
               call m2r_fields%create_comm(3)
               call r2m_fields%create_comm(2)
            end if
         end if
      else
         call m2r_single%create_comm(1)
         call r2m_single%create_comm(1)
      end if

   end subroutine initialize_communications
!------------------------------------------------------------------------------
   subroutine finalize_communications

      if ( l_packed_transp ) then
         call r2m_fields%destroy_comm()
         call m2r_fields%destroy_comm()
      else
         call m2r_single%destroy_comm()
      end if
      call r2m_single%destroy_comm()

   end subroutine finalize_communications
!------------------------------------------------------------------------------
   subroutine gather_from_mloc_to_rank0(arr_Mloc, arr_full)

      !-- Input variable:
      complex(cp), intent(in) :: arr_Mloc(nMstart:nMstop, n_r_max)

      !-- Output variable:
      complex(cp), intent(out) :: arr_full(n_m_max, n_r_max)

      !-- Local variables:
      complex(cp), allocatable :: rbuff(:), sbuff(:)
      integer, allocatable :: rcounts(:)
      integer, allocatable :: rdisp(:)
      integer ::p, ii, n_r, n_m

      allocate( rbuff(n_m_max*n_r_max), sbuff(nm_per_rank*n_r_max) )
      allocate ( rcounts(0:n_procs-1), rdisp(0:n_procs-1) )

      do p=0,n_procs-1
         rcounts(p)=n_r_max*m_balance(p)%n_per_rank
      end do

      rdisp(0)=0
      do p=1,n_procs-1
         rdisp(p)=rdisp(p-1)+rcounts(p-1)
      end do

      ii = 1
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            sbuff(ii)=arr_Mloc(n_m,n_r)
            ii = ii +1
         end do
      end do

      call MPI_Gatherv(sbuff, nm_per_rank*n_r_max, MPI_DEF_COMPLEX, &
           &           rbuff, rcounts, rdisp, MPI_DEF_COMPLEX, 0,   &
           &           MPI_COMM_WORLD, ierr)

      if ( rank == 0 ) then
         do p = 0, n_procs-1
            ii = rdisp(p)+1
            do n_r=1,n_r_max
               do n_m=m_balance(p)%nStart,m_balance(p)%nStop
                  arr_full(n_m,n_r)=rbuff(ii)
                  ii=ii+1
               end do
            end do
         end do
      end if

      deallocate( rbuff, sbuff, rcounts, rdisp)

   end subroutine gather_from_mloc_to_rank0
!------------------------------------------------------------------------------
   subroutine scatter_from_rank0_to_mloc(arr_full, arr_Mloc)

      !-- Input variable:
      complex(cp), intent(in) :: arr_full(n_m_max, n_r_max)

      !-- Output variable:
      complex(cp), intent(out) :: arr_Mloc(nMstart:nMstop, n_r_max)

      !-- Local variables:
      complex(cp), allocatable :: rbuff(:), sbuff(:)
      integer, allocatable :: scounts(:)
      integer, allocatable :: sdisp(:)
      integer :: p, ii, n_r, n_m

      allocate( rbuff(nm_per_rank*n_r_max), sbuff(n_m_max*n_r_max) )
      allocate ( scounts(0:n_procs-1), sdisp(0:n_procs-1) )

      do p=0,n_procs-1
         scounts(p)=n_r_max*m_balance(p)%n_per_rank
      end do

      sdisp(0)=0
      do p=1,n_procs-1
         sdisp(p)=sdisp(p-1)+scounts(p-1)
      end do

      if ( rank == 0 ) then
         do p = 0, n_procs-1
            ii = sdisp(p)+1
            do n_r=1,n_r_max
               do n_m=m_balance(p)%nStart,m_balance(p)%nStop
                  sbuff(ii) = arr_full(n_m,n_r)
                  ii=ii+1
               end do
            end do
         end do
      end if

      call MPI_Scatterv(sbuff, scounts, sdisp, MPI_DEF_COMPLEX,         &
           &            rbuff, nm_per_rank*n_r_max, MPI_DEF_COMPLEX, 0, &
           &            MPI_COMM_WORLD, ierr)

      ii = 1
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            arr_Mloc(n_m,n_r)=rbuff(ii)
            ii = ii +1
         end do
      end do

      deallocate ( scounts, sdisp, rbuff, sbuff )

   end subroutine scatter_from_rank0_to_mloc
!------------------------------------------------------------------------------
   subroutine gather_from_Rloc(arr_Rloc, arr_glob, irank)
      !
      ! This subroutine gather a r-distributed array on rank=irank
      !

      !-- Input variable
      real(cp), intent(in) :: arr_Rloc(nRstart:nRstop)
      integer,  intent(in) :: irank

      !-- Output variable
      real(cp), intent(out) :: arr_glob(1:n_r_max)

      !-- Local variables:
      integer :: p
      integer :: scount, rcounts(0:n_procs-1), rdisp(0:n_procs-1)

      scount = nRstop-nRstart+1
      do p=0,n_procs-1
         rcounts(p)=radial_balance(p)%n_per_rank
      end do
      rdisp(0)=0
      do p=1,n_procs-1
         rdisp(p)=rdisp(p-1)+rcounts(p-1)
      end do

      call MPI_GatherV(arr_Rloc, scount, MPI_DEF_REAL, arr_glob, rcounts, &
           &           rdisp, MPI_DEF_REAL, irank, MPI_COMM_WORLD, ierr)

   end subroutine gather_from_Rloc
!------------------------------------------------------------------------------
   subroutine scatter_from_rank0_to_Rloc(arr_full, arr_Rloc)

      !-- Input variable:
      complex(cp), intent(in) :: arr_full(n_m_max, n_r_max)

      !-- Output variable:
      complex(cp), intent(out) :: arr_Rloc(n_m_max, nRstart:nRstop)

      !-- Local variables:
      complex(cp), allocatable :: rbuff(:), sbuff(:)
      integer, allocatable :: scounts(:)
      integer, allocatable :: sdisp(:)
      integer :: p, ii, n_r, n_m

      allocate( rbuff(nR_per_rank*n_m_max), sbuff(n_m_max*n_r_max) )
      allocate ( scounts(0:n_procs-1), sdisp(0:n_procs-1) )

      do p=0,n_procs-1
         scounts(p)=n_m_max*radial_balance(p)%n_per_rank
      end do

      sdisp(0)=0
      do p=1,n_procs-1
         sdisp(p)=sdisp(p-1)+scounts(p-1)
      end do

      if ( rank == 0 ) then
         do p = 0, n_procs-1
            ii = sdisp(p)+1
            do n_r=radial_balance(p)%nStart,radial_balance(p)%nStop
               do n_m=1,n_m_max
                  sbuff(ii) = arr_full(n_m,n_r)
                  ii=ii+1
               end do
            end do
         end do
      end if

      call MPI_Scatterv(sbuff, scounts, sdisp, MPI_DEF_COMPLEX,         &
           &            rbuff, nR_per_rank*n_m_max, MPI_DEF_COMPLEX, 0, &
           &            MPI_COMM_WORLD, ierr)

      ii = 1
      do n_r=nRstart,nRstop
         do n_m=1,n_m_max
            arr_Rloc(n_m,n_r)=rbuff(ii)
            ii = ii +1
         end do
      end do

      deallocate ( scounts, sdisp, rbuff, sbuff )

   end subroutine scatter_from_rank0_to_Rloc
!------------------------------------------------------------------------------
   subroutine reduce_radial_on_rank(arr_dist, irank)

      !-- Input variable
      integer,  intent(in) :: irank

      !-- Output variable
      real(cp), intent(inout) :: arr_dist(n_r_max)

      !-- Local variable
      integer :: n_r
      real(cp) :: work(n_r_max)

      call MPI_Reduce(arr_dist, work, n_r_max, MPI_DEF_REAL, &
           &          MPI_SUM, irank, MPI_COMM_WORLD, ierr)

      if ( rank == irank ) then
         do n_r=1,n_r_max
            arr_dist(n_r) = work(n_r)
         end do
      end if

   end subroutine reduce_radial_on_rank
!------------------------------------------------------------------------------
   subroutine my_reduce_mean(scalar, irank)

      !-- Input variable
      integer, intent(in) :: irank

      !-- Output variable
      real(cp), intent(inout) :: scalar

      !-- Local variable
      real(cp) :: tmp

      call MPI_Reduce(scalar, tmp, 1, MPI_DEF_REAL, MPI_SUM, &
           &          irank, MPI_COMM_WORLD, ierr)

      if ( rank == irank ) scalar = tmp/real(n_procs,cp)

   end subroutine my_reduce_mean
!------------------------------------------------------------------------------
   function my_allreduce_maxloc(arr, nMstart) result(ind)
      !
      ! This function is the MPI version of the intrinsic Fortran 'maxloc'
      ! function
      !

      !-- Input variable:
      integer,  intent(in) :: nMstart
      real(cp), intent(in) :: arr(:)

      !-- Local variables:
      integer :: ind
      real(cp) :: idx(2), tmp(2)

      idx(2) = maxloc(arr,dim=1)
      idx(1) = maxval(arr)
      idx(2) = nMstart-1+idx(2)
      call MPI_AllReduce(idx, tmp, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, &
           &             MPI_COMM_WORLD, ierr)
      ind = int(tmp(2))

   end function my_allreduce_maxloc
!------------------------------------------------------------------------------
   subroutine find_faster_comm(idx,n_fields,n_log_file)
      !
      ! This subroutine tests two MPI transposition strategies and
      ! selects the fastest one.
      !

      !-- Input variables
      integer, intent(in) :: n_fields   ! number of fields
      integer, intent(in) :: n_log_file ! log.TAG file unit

      !-- Output variable:
      integer,  intent(out) :: idx

      !-- Local variables
      class(type_mpitransp), pointer :: m2r_test
      complex(cp) :: arr_Rloc(n_m_max,nRstart:nRstop,n_fields)
      complex(cp) :: arr_Mloc(nMstart:nMstop,n_r_max,n_fields)
      real(cp) :: tStart, tStop, tAlltoAllv, tAlltoAllw, minTime
      real(cp) :: rdm_real, rdm_imag, timers(2)
      integer :: n_f, n_r, n_m, n_t, n, n_out, ind(1)
      integer, parameter :: n_transp=50
      character(len=80) :: message

      !-- First fill an array with random numbers
      do n_f=1,n_fields
         do n_r=nRstart,nRstop
            do n_m=1,n_m_max
               call random_number(rdm_real)
               call random_number(rdm_imag)
               arr_Rloc(n_m,n_r,n_f)=cmplx(rdm_real,rdm_imag,kind=cp)
            end do
         end do
      end do

      !-- Try the all-to-allv strategy (50 back and forth transposes)
      allocate( type_mpiatoav :: m2r_test )
      call m2r_test%create_comm(n_fields)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      tStart = MPI_Wtime()
      do n_t=1,n_transp
         call m2r_test%transp_r2m(arr_Rloc, arr_Mloc)
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
         call m2r_test%transp_m2r(arr_Mloc, arr_Rloc)
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
      end do
      tStop = MPI_Wtime()
      tAlltoAllv = tStop-tStart
      call m2r_test%destroy_comm()
      deallocate( m2r_test)

      !-- Try the all-to-allw strategy (50 back and forth transposes)
      allocate( type_mpiatoaw :: m2r_test )
      call m2r_test%create_comm(n_fields)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      tStart = MPI_Wtime()
      do n_t=1,n_transp
         call m2r_test%transp_r2m(arr_Rloc, arr_Mloc)
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
         call m2r_test%transp_m2r(arr_Mloc, arr_Rloc)
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
      end do
      tStop = MPI_Wtime()
      tAlltoAllw = tStop-tStart
      call m2r_test%destroy_comm()
      deallocate( m2r_test)

      !-- Now determine the average over the ranks and send it to rank=0
      call MPI_Reduce(tAlltoAllv, timers(1), 1, MPI_DEF_REAL, MPI_SUM, 0, &
           &          MPI_COMM_WORLD, ierr)
      call MPI_Reduce(tAlltoAllw, timers(2), 1, MPI_DEF_REAL, MPI_SUM, 0, &
           &          MPI_COMM_WORLD, ierr)

      if ( rank == 0 ) then
         !-- Average over procs and number of transposes
         timers(:) = timers(:)/real(n_procs,cp)/real(n_transp,cp)

         !-- Determine the fastest
         ind = minloc(timers)
         minTime = minval(timers)
         idx = ind(1)

         do n=1,2
            if ( n==1 ) then
               n_out = n_log_file
            else
               n_out = output_unit
            end if
            write(n_out,*)
            if ( n_fields == 1 ) then
               write(message,'('' ! MPI transpose strategy for '', I1,'' field'')') &
               &     n_fields
            else
               write(message,'('' ! MPI transpose strategy for '', I1,'' fields'')') &
               &     n_fields
            end if
            write(n_out,'(A80)') message
            write(message,'('' ! alltoallv communicator          ='', &
            &               ES10.3, '' s'')') timers(1)
            write(n_out,'(A80)') message
            write(message,'('' ! alltoallw communicator          ='', &
            &               ES10.3, '' s'')') timers(2)
            write(n_out,'(A80)') message
         end do

      end if

      call MPI_Bcast(idx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(minTime,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)

   end subroutine find_faster_comm
!------------------------------------------------------------------------------
end module communications

