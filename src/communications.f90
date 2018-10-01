module communications

   use mpimod
   use precision_mod
   use blocking
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, n_m_max
   use parallel_mod, only: n_procs, rank, ierr

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

   type(help_transp), public :: r2m_fields
   type(help_transp), public :: m2r_fields

   public :: initialize_communications, transp_r2m, transp_m2r,     &
   &         gather_from_mloc_to_rank0, scatter_from_rank0_to_mloc, &
   &         finalize_communications, reduce_radial_on_rank,        &
   &         my_reduce_mean, my_allreduce_maxloc

contains

   subroutine initialize_communications()

      call create_r2m_type(r2m_fields)
      call create_m2r_type(m2r_fields)

   end subroutine initialize_communications
!------------------------------------------------------------------------------
   subroutine finalize_communications

      call destroy_communicator(m2r_fields)
      call destroy_communicator(r2m_fields)

   end subroutine finalize_communications
!------------------------------------------------------------------------------
   subroutine create_r2m_type(self)

      type(help_transp) :: self
      integer :: p

      allocate ( self%rcounts(0:n_procs-1), self%scounts(0:n_procs-1) )
      allocate ( self%rdisp(0:n_procs-1), self%sdisp(0:n_procs-1) )

      do p=0,n_procs-1
         self%scounts(p)=nR_per_rank*m_balance(p)%n_per_rank
         self%rcounts(p)=radial_balance(p)%n_per_rank*nm_per_rank
      end do

      self%rdisp(0)=0
      self%sdisp(0)=0
      do p=1,n_procs-1
         self%sdisp(p)=self%sdisp(p-1)+self%scounts(p-1)
         self%rdisp(p)=self%rdisp(p-1)+self%rcounts(p-1)
      end do

      self%max_send = sum(self%scounts)
      self%max_recv = sum(self%rcounts)

      bytes_allocated = bytes_allocated+4*n_procs*SIZEOF_INTEGER

      allocate( self%sbuff(1:self%max_send) )
      allocate( self%rbuff(1:self%max_recv) )

      bytes_allocated = bytes_allocated+(self%max_send+self%max_recv)*&
      &                 SIZEOF_DEF_COMPLEX

   end subroutine create_r2m_type
!------------------------------------------------------------------------------
   subroutine destroy_communicator(self)

      type(help_transp) :: self

      deallocate( self%rbuff, self%sbuff )
      deallocate( self%sdisp, self%rdisp )
      deallocate( self%scounts, self%rcounts )

   end subroutine destroy_communicator
!------------------------------------------------------------------------------
   subroutine create_m2r_type(self)

      type(help_transp) :: self
      integer :: p

      allocate ( self%rcounts(0:n_procs-1), self%scounts(0:n_procs-1) )
      allocate ( self%rdisp(0:n_procs-1), self%sdisp(0:n_procs-1) )

      do p=0,n_procs-1
         self%scounts(p)=radial_balance(p)%n_per_rank*nm_per_rank
         self%rcounts(p)=nR_per_rank*m_balance(p)%n_per_rank
      end do

      self%rdisp(0)=0
      self%sdisp(0)=0
      do p=1,n_procs-1
         self%sdisp(p)=self%sdisp(p-1)+self%scounts(p-1)
         self%rdisp(p)=self%rdisp(p-1)+self%rcounts(p-1)
      end do

      self%max_send = sum(self%scounts)
      self%max_recv = sum(self%rcounts)

      bytes_allocated = bytes_allocated+4*n_procs*SIZEOF_INTEGER

      allocate( self%sbuff(1:self%max_send) )
      allocate( self%rbuff(1:self%max_recv) )

   end subroutine create_m2r_type
!------------------------------------------------------------------------------
   subroutine transp_r2m(self, arr_Rloc, arr_Mloc)

      !-- Input variables
      type(help_transp), intent(inout) :: self
      complex(cp),       intent(in) :: arr_Rloc(n_m_max,nRstart:nRstop)

      !-- Output variable
      complex(cp), intent(out) :: arr_Mloc(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: p, ii, n_r, n_m

      do p = 0, n_procs-1
         ii = self%sdisp(p)+1
         do n_r=nRstart,nRstop
            do n_m=m_balance(p)%nStart,m_balance(p)%nStop
               self%sbuff(ii)=arr_Rloc(n_m,n_r)
               ii = ii +1
            end do
         end do
      end do

      call MPI_Alltoallv(self%sbuff, self%scounts, self%sdisp, MPI_DEF_COMPLEX, &
           &             self%rbuff, self%rcounts, self%rdisp, MPI_DEF_COMPLEX, &
           &             MPI_COMM_WORLD, ierr) 

      do p = 0, n_procs-1
         ii = self%rdisp(p)+1
         do n_r=radial_balance(p)%nStart,radial_balance(p)%nStop
            do n_m=nMstart,nMstop
               arr_Mloc(n_m,n_r)=self%rbuff(ii)
               ii=ii+1
            end do
         end do
      end do

   end subroutine transp_r2m
!------------------------------------------------------------------------------
   subroutine transp_m2r(self, arr_Mloc, arr_Rloc)

      !-- Input variables
      type(help_transp), intent(inout) :: self
      complex(cp),       intent(in) :: arr_Mloc(nMstart:nMstop,n_r_max)

      !-- Output variable
      complex(cp), intent(out) :: arr_Rloc(n_m_max,nRstart:nRstop)

      !-- Local variables
      integer :: p, ii, n_r, n_m

      do p = 0, n_procs-1
         ii = self%sdisp(p)+1
         do n_r=radial_balance(p)%nStart,radial_balance(p)%nStop
            do n_m=nMstart,nMstop
               self%sbuff(ii)=arr_Mloc(n_m,n_r)
               ii = ii+1
            end do
         end do
      end do

      call MPI_Alltoallv(self%sbuff, self%scounts, self%sdisp, MPI_DEF_COMPLEX, &
           &             self%rbuff, self%rcounts, self%rdisp, MPI_DEF_COMPLEX, &
           &             MPI_COMM_WORLD, ierr) 

      do p = 0, n_procs-1
         ii = self%rdisp(p)+1
         do n_r=nRstart,nRstop
            do n_m=m_balance(p)%nStart,m_balance(p)%nStop
               arr_Rloc(n_m,n_r)=self%rbuff(ii)
               ii=ii+1
            end do
         end do
      end do

   end subroutine transp_m2r
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

   end subroutine scatter_from_rank0_to_mloc
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
   function my_allreduce_maxloc(arr) result(ind)
      !
      ! This function is the MPI version of the intrinsic Fortran 'maxloc'
      ! function
      !

      !-- Input variable:
      real(cp), intent(in) :: arr(:)

      !-- Local variables:
      integer :: ind
      real(cp) :: idx(2), tmp(2)

      idx(2) = maxloc(arr,dim=1)
      idx(1) = maxval(arr)
      call MPI_AllReduce(idx, tmp, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, &
           &             MPI_COMM_WORLD, ierr)
      ind = int(tmp(2))

   end function my_allreduce_maxloc
!------------------------------------------------------------------------------
end module communications
