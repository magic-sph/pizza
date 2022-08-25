module communications

   use mpimod
   use precision_mod
   use blocking
   use blocking_lm, only: st_map, lo_map
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, n_m_max
   use truncation_3D, only: n_r_max_3D, lm_max, n_phi_max_3D
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
   type(help_transp), public :: lm2r_fields
   type(help_transp), public :: r2lm_fields
   type(help_transp), public :: r2all_fields

   public :: initialize_communications, transp_r2m, transp_m2r,     &
   &         gather_from_mloc_to_rank0, scatter_from_rank0_to_mloc, &
   &         finalize_communications, reduce_radial_on_rank,        &
   &         my_reduce_mean, my_allreduce_maxloc, transp_lm2r,      &
   &         transp_r2lm, allgather_from_Rloc,                      &
   &         allgather_from_Rloc_3D, scatter_from_rank0_to_lmloc,   &
   &         scatter_from_rank0_to_rloc, exchange_Nbound_from_Rloc_3D

contains

   subroutine initialize_communications(l_3D)

      logical, intent(in) :: l_3D

      call create_r2m_type(r2m_fields)
      call create_m2r_type(m2r_fields)
      call create_r2all_type(r2all_fields)

      if ( l_3D ) then
         call create_lm2r_type(lm2r_fields)
         call create_r2lm_type(r2lm_fields)
      end if

   end subroutine initialize_communications
!------------------------------------------------------------------------------
   subroutine finalize_communications(l_3D)

      logical, intent(in) :: l_3D

      if ( l_3D ) then
         call destroy_communicator(r2lm_fields)
         call destroy_communicator(lm2r_fields)
      end if

      call destroy_communicator(m2r_fields)
      call destroy_communicator(r2m_fields)
      call destroy_communicator(r2all_fields)

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
   subroutine create_lm2r_type(self)

      type(help_transp) :: self
      integer :: p

      allocate ( self%rcounts(0:n_procs-1), self%scounts(0:n_procs-1) )
      allocate ( self%rdisp(0:n_procs-1), self%sdisp(0:n_procs-1) )

      do p=0,n_procs-1
         self%scounts(p)=radial_balance_3D(p)%n_per_rank*nlm_per_rank
         self%rcounts(p)=nR_per_rank_3D*lm_balance(p)%n_per_rank
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

   end subroutine create_lm2r_type
!------------------------------------------------------------------------------
   subroutine create_r2lm_type(self)

      type(help_transp) :: self
      integer :: p

      allocate ( self%rcounts(0:n_procs-1), self%scounts(0:n_procs-1) )
      allocate ( self%rdisp(0:n_procs-1), self%sdisp(0:n_procs-1) )

      do p=0,n_procs-1
         self%scounts(p)=nR_per_rank_3D*lm_balance(p)%n_per_rank
         self%rcounts(p)=radial_balance_3D(p)%n_per_rank*nlm_per_rank
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

   end subroutine create_r2lm_type
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
   subroutine create_r2all_type(self)

      type(help_transp) :: self
      integer :: p

      allocate ( self%rcounts(0:n_procs-1), self%scounts(0:n_procs-1) )
      allocate ( self%rdisp(0:n_procs-1), self%sdisp(0:n_procs-1) )

      do p=0,n_procs-1
         self%scounts(p)=nR_per_rank*radial_balance(p)%n_per_rank
         self%rcounts(p)=n_r_max
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

   end subroutine create_r2all_type
!------------------------------------------------------------------------------
   subroutine transp_lm2r(self, arr_LMloc, arr_Rloc)

      !-- Input variables
      type(help_transp), intent(inout) :: self
      complex(cp),       intent(in) :: arr_LMloc(lmStart:lmStop, n_r_max_3D)

      !-- Output variables
      complex(cp), intent(out) :: arr_Rloc(lm_max, nRstart3D:nRstop3D)

      !-- Local variables
      integer :: p, ii, n_r, lm, l, m, lm_st

      do p = 0, n_procs-1
         ii = self%sdisp(p)+1
         do n_r=radial_balance_3D(p)%nStart,radial_balance_3D(p)%nStop
            do lm=lmStart,lmStop
               self%sbuff(ii)=arr_LMloc(lm,n_r)
               ii = ii+1
            end do
         end do
      end do

      call MPI_Alltoallv(self%sbuff, self%scounts, self%sdisp, MPI_DEF_COMPLEX, &
           &             self%rbuff, self%rcounts, self%rdisp, MPI_DEF_COMPLEX, &
           &             MPI_COMM_WORLD, ierr) 

      do p = 0, n_procs-1
         ii = self%rdisp(p)+1
         do n_r=nRstart3D,nRstop3D
            do lm=lm_balance(p)%nStart,lm_balance(p)%nStop
               !arr_Rloc(lm,n_r)=self%rbuff(ii)
               l = lo_map%lm2l(lm)
               m = lo_map%lm2m(lm)
               lm_st = st_map%lm2(l,m)
               arr_Rloc(lm_st,n_r)=self%rbuff(ii)
               ii=ii+1
            end do
         end do
      end do

   end subroutine transp_lm2r
!------------------------------------------------------------------------------
   subroutine transp_r2lm(self, arr_Rloc, arr_LMloc)

      !-- Input variables
      type(help_transp), intent(inout) :: self
      complex(cp),       intent(in) :: arr_Rloc(lm_max,nRstart3D:nRstop3D)

      !-- Output variable
      complex(cp), intent(out) :: arr_LMloc(lMstart:lMstop,n_r_max_3D)

      !-- Local variables
      integer :: p, ii, n_r, lm, l, m, lm_st

      do p = 0, n_procs-1
         ii = self%sdisp(p)+1
         do n_r=nRstart3D,nRstop3D
            do lm=lm_balance(p)%nStart,lm_balance(p)%nStop
               l = lo_map%lm2l(lm)
               m = lo_map%lm2m(lm)
               lm_st = st_map%lm2(l,m)
               self%sbuff(ii)=arr_Rloc(lm_st,n_r)
               ii = ii +1
            end do
         end do
      end do

      call MPI_Alltoallv(self%sbuff, self%scounts, self%sdisp, MPI_DEF_COMPLEX, &
           &             self%rbuff, self%rcounts, self%rdisp, MPI_DEF_COMPLEX, &
           &             MPI_COMM_WORLD, ierr) 

      do p = 0, n_procs-1
         ii = self%rdisp(p)+1
         do n_r=radial_balance_3D(p)%nStart,radial_balance_3D(p)%nStop
            do lm=lmStart,lmStop
               arr_LMloc(lm,n_r)=self%rbuff(ii)
               ii=ii+1
            end do
         end do
      end do

   end subroutine transp_r2lm
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
   subroutine allgather_from_Rloc(arr_Rloc, arr_full, len_arr)
      !
      ! This routine allgather the R-distributed array to a global array
      !

      !-- Input variable:
      integer,  intent(in) :: len_arr
      real(cp), intent(in) :: arr_Rloc(len_arr, nRstart:nRstop)

      !-- Output variable:
      real(cp), intent(out) :: arr_full(len_arr, n_r_max)

      !-- Local variables:
      real(cp), allocatable :: rbuff(:), sbuff(:)
      integer, allocatable :: rcounts(:)
      integer, allocatable :: rdisp(:)
      integer :: p, ii, n_r, n_m

      allocate( rbuff(len_arr*n_r_max), sbuff(len_arr*nR_per_rank) )
      allocate ( rcounts(0:n_procs-1), rdisp(0:n_procs-1) )

      do p=0,n_procs-1
         rcounts(p)=len_arr*radial_balance(p)%n_per_rank
      end do

      rdisp(0)=0
      do p=1,n_procs-1
         rdisp(p)=rdisp(p-1)+rcounts(p-1)
      end do

      ii = 1
      do n_r=nRstart,nRstop
         do n_m=1,len_arr
            sbuff(ii)=arr_Rloc(n_m,n_r)
            ii = ii +1
         end do
      end do

      call MPI_Allgatherv(sbuff, nR_per_rank*len_arr, MPI_DEF_REAL, &
           &              rbuff, rcounts, rdisp, MPI_DEF_REAL,      &
           &              MPI_COMM_WORLD, ierr)

      do p = 0, n_procs-1
         ii = rdisp(p)+1
         do n_r=radial_balance(p)%nStart,radial_balance(p)%nStop
            do n_m=1,len_arr
               arr_full(n_m,n_r)=rbuff(ii)
               ii=ii+1
            end do
         end do
      end do

      deallocate( rbuff, sbuff, rcounts, rdisp)

   end subroutine allgather_from_Rloc
!------------------------------------------------------------------------------
   subroutine allgather_from_Rloc_3D(arr_Rloc, arr_full, len_arr)
      !
      ! This routine allgather the R-distributed array to a global array
      !

      !-- Input variable:
      integer,  intent(in) :: len_arr
      real(cp), intent(in) :: arr_Rloc(len_arr, nRstart3D:nRstop3D)

      !-- Output variable:
      real(cp), intent(out) :: arr_full(len_arr, n_r_max_3D)

      !-- Local variables:
      real(cp), allocatable :: rbuff(:), sbuff(:)
      integer, allocatable :: rcounts(:)
      integer, allocatable :: rdisp(:)
      integer :: p, ii, n_r, n_m

      allocate( rbuff(len_arr*n_r_max_3D), sbuff(len_arr*nR_per_rank_3D) )
      allocate ( rcounts(0:n_procs-1), rdisp(0:n_procs-1) )

      do p=0,n_procs-1
         rcounts(p)=len_arr*radial_balance_3D(p)%n_per_rank
      end do

      rdisp(0)=0
      do p=1,n_procs-1
         rdisp(p)=rdisp(p-1)+rcounts(p-1)
      end do

      ii = 1
      do n_r=nRstart3D,nRstop3D
         do n_m=1,len_arr
            sbuff(ii)=arr_Rloc(n_m,n_r)
            ii = ii +1
         end do
      end do

      call MPI_Allgatherv(sbuff, nR_per_rank_3D*len_arr, MPI_DEF_REAL, &
           &              rbuff, rcounts, rdisp, MPI_DEF_REAL,         &
           &              MPI_COMM_WORLD, ierr)

      do p = 0, n_procs-1
         ii = rdisp(p)+1
         do n_r=radial_balance_3D(p)%nStart,radial_balance_3D(p)%nStop
            do n_m=1,len_arr
               arr_full(n_m,n_r)=rbuff(ii)
               ii=ii+1
            end do
         end do
      end do

      deallocate( rbuff, sbuff, rcounts, rdisp)

   end subroutine allgather_from_Rloc_3D
!------------------------------------------------------------------------------
   subroutine exchange_Nbound_from_Rloc_3D(arr_Rloc, arr_Rext, Nbound, len_arr)
      !
      ! This routine exchange Nbound edge points from an 3D-R-distributed array 
      ! to an extended 3D_Rloc array: rank N sends nRstart3D+Nbound pts to N-1
      ! and nRstop3D-Nbound pts to N+1 resp.; and receives nRstart3D+Nbound and
      ! nRstop3D-Nbound pts from N+1 and N-1 respectively
      !

      !-- Input variable:
      integer,  intent(in) :: Nbound, len_arr
      real(cp), intent(in) :: arr_Rloc(len_arr, nRstart3D:nRstop3D)

      !-- Output variable:
      real(cp), intent(out) :: arr_Rext(len_arr, nRstart3D-Nbound:nRstop3D+Nbound)

      !-- Local variables:
      real(cp), allocatable :: rbuff(:,:), sbuff(:,:)
      integer :: to, up, from, low
      integer :: n_m, mrank, status(MPI_STATUS_SIZE)

      allocate( rbuff(len_arr,Nbound), sbuff(len_arr,Nbound) )

      low=rank-1
      up =rank+1
      if ( low < 0 ) low=n_procs-1
      if ( up == n_procs ) up =0

      mrank = mod(rank,2)

      do n_m=1,len_arr
         arr_Rext(n_m,nRstart3D:nRstop3D) = arr_Rloc(n_m,:)
      end do

      !print*, rbuff(4,:), arr_Rext(4,:)
      !print*, 'checking NBound', Nbound
      !print*, arr_Rext(4,nRstart3D-Nbound:nRstart-1), " <- nRstart - ",  arr_Rext(4,nRstop3D-Nbound+1:nRstop3D), " - nRstop -> "
      !print*, arr_Rext(4,nRstop3D+1:nRstop3D+Nbound), " <- nRstop + ",  arr_Rext(4,nRstart3D:nRstart3D+Nbound-1), " + nRstart -> "

      !-- Send to up receive from low
      from=low
      to  =up

      if ( mrank == 0 .and. rank > 0 ) then
         call MPI_Recv(rbuff, Nbound*len_arr, MPI_DEF_REAL, &
              &        from, MPI_ANY_TAG,                   &
              &        MPI_COMM_WORLD, status, ierr)
         arr_Rext(:,nRstart3D-Nbound:nRstart-1) = rbuff(:,:)
      else if ( mrank > 0 .and. rank < n_procs-1 ) then
         sbuff(:,:) = arr_Rext(:,nRstop3D-Nbound+1:nRstop3D)
         call MPI_Send(sbuff, Nbound*len_arr, MPI_DEF_REAL, &
              &        to, rank,                            &
              &        MPI_COMM_WORLD, ierr)
      end if

      if ( mrank > 0 .and. rank > 0 ) then
         call MPI_Recv(rbuff, Nbound*len_arr, MPI_DEF_REAL, &
              &        from, MPI_ANY_TAG,                   &
              &        MPI_COMM_WORLD, status, ierr)
         arr_Rext(:,nRstart3D-Nbound:nRstart-1) = rbuff(:,:)
      else if ( mrank == 0 .and. rank < n_procs-1 ) then
         sbuff(:,:) = arr_Rext(:,nRstop3D-Nbound+1:nRstop3D)
         call MPI_Send(sbuff, Nbound*len_arr, MPI_DEF_REAL, &
              &        to, rank,                            &
              &        MPI_COMM_WORLD, ierr)
      end if

      !-- Send to low receive from up
      from=up
      to  =low

      if ( mrank == 0 .and. rank < n_procs-1 ) then
         call MPI_Recv(rbuff, Nbound*len_arr, MPI_DEF_REAL, &
              &        from, MPI_ANY_TAG,                   &
              &        MPI_COMM_WORLD, status, ierr)
         arr_Rext(:,nRstop3D+1:nRstop3D+Nbound) = rbuff(:,:)
      else if ( mrank > 0 .and. rank > 0 ) then
         sbuff(:,:) = arr_Rext(:,nRstart3D:nRstart3D+Nbound-1)
         call MPI_Send(sbuff, Nbound*len_arr, MPI_DEF_REAL, &
              &        to, rank,                            &
              &        MPI_COMM_WORLD, ierr)
      end if

      if ( mrank > 0 .and. rank < n_procs-1 ) then
         call MPI_Recv(rbuff, Nbound*len_arr, MPI_DEF_REAL, &
              &        from, MPI_ANY_TAG,                   &
              &        MPI_COMM_WORLD, status, ierr)
         arr_Rext(:,nRstop3D+1:nRstop3D+Nbound) = rbuff(:,:)
      else if ( mrank == 0 .and. rank > 0 ) then
         sbuff(:,:) = arr_Rext(:,nRstart3D:nRstart3D+Nbound-1)
         call MPI_Send(sbuff, Nbound*len_arr, MPI_DEF_REAL, &
              &        to, rank,                            &
              &        MPI_COMM_WORLD, ierr)
      end if

      deallocate( rbuff, sbuff )

   end subroutine exchange_Nbound_from_Rloc_3D
!------------------------------------------------------------------------------
   subroutine scatter_from_rank0_to_rloc(arr_full, arr_Rloc, len_arr)!, N_rank)
      !
      ! This routine scatter a global array from rank 0 to a radial distributed array
      !

      !-- Input variable:
      integer, intent(in) :: len_arr!N_rank,
      complex(cp), intent(in) :: arr_full(len_arr, n_r_max)

      !-- Output variable:
      complex(cp), intent(out) :: arr_Rloc(len_arr,nRstart:nRstop)

      !-- Local variables:
      complex(cp), allocatable :: rbuff(:), sbuff(:)
      integer, allocatable :: scounts(:)
      integer, allocatable :: sdisp(:)
      integer :: p, ii, n_r, n_m

      allocate( rbuff(len_arr*nR_per_rank), sbuff(len_arr*n_r_max) )
      allocate ( scounts(0:n_procs-1), sdisp(0:n_procs-1) )

      do p=0,n_procs-1
         scounts(p)=len_arr*radial_balance(p)%n_per_rank
      end do

      sdisp(0)=0
      do p=1,n_procs-1
         sdisp(p)=sdisp(p-1)+scounts(p-1)
      end do

      if ( rank == 0 ) then!N_rank ) then
         do p = 0, n_procs-1
            ii = sdisp(p)+1
            do n_r=radial_balance(p)%nStart,radial_balance(p)%nStop
               do n_m=1,len_arr
                  sbuff(ii) = arr_full(n_m,n_r)
                  ii=ii+1
               end do
            end do
         end do
      end if

      call MPI_Scatterv(sbuff, scounts, sdisp, MPI_DEF_COMPLEX,         &
           &            rbuff, len_arr*nR_per_rank, MPI_DEF_COMPLEX, 0, &
           &            MPI_COMM_WORLD, ierr)

      ii = 1
      do n_r=nRstart,nRstop
         do n_m=1,len_arr
            arr_Rloc(n_m,n_r)=rbuff(ii)
            ii = ii +1
         end do
      end do

      deallocate( rbuff, sbuff, scounts, sdisp)

   end subroutine scatter_from_rank0_to_rloc
!------------------------------------------------------------------------------
   subroutine scatter_from_rank0_to_mloc(arr_full, arr_Mloc)
      !
      ! This routine scatter a global array from rank 0 to a m-distributed array
      !

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

      deallocate( rbuff, sbuff, scounts, sdisp)

   end subroutine scatter_from_rank0_to_mloc
!------------------------------------------------------------------------------
   subroutine scatter_from_rank0_to_lmloc(arr_full, arr_LMloc)
      !
      ! This routine scatter a global array from rank 0 to a lm-distributed array
      !

      !-- Input variable:
      complex(cp), intent(in) :: arr_full(lm_max, n_r_max_3D)

      !-- Output variable:
      complex(cp), intent(out) :: arr_LMloc(lmStart:lmStop, n_r_max_3D)

      !-- Local variables:
      complex(cp), allocatable :: rbuff(:), sbuff(:)
      integer, allocatable :: scounts(:)
      integer, allocatable :: sdisp(:)
      integer :: p, ii, n_r, lm, l, m, lm_st

      allocate( rbuff(nlm_per_rank*n_r_max_3D), sbuff(lm_max*n_r_max_3D) )
      allocate ( scounts(0:n_procs-1), sdisp(0:n_procs-1) )

      do p=0,n_procs-1
         scounts(p)=n_r_max_3D*lm_balance(p)%n_per_rank
      end do

      sdisp(0)=0
      do p=1,n_procs-1
         sdisp(p)=sdisp(p-1)+scounts(p-1)
      end do

      if ( rank == 0 ) then
         do p = 0, n_procs-1
            ii = sdisp(p)+1
            do n_r=1,n_r_max_3D
               do lm=lm_balance(p)%nStart,lm_balance(p)%nStop
                  l = lo_map%lm2l(lm)
                  m = lo_map%lm2m(lm)
                  lm_st = st_map%lm2(l,m)
                  sbuff(ii) = arr_full(lm_st,n_r)
                  ii=ii+1
               end do
            end do
         end do
      end if

      call MPI_Scatterv(sbuff, scounts, sdisp, MPI_DEF_COMPLEX,             &
           &            rbuff, nlm_per_rank*n_r_max_3D, MPI_DEF_COMPLEX, 0, &
           &            MPI_COMM_WORLD, ierr)

      ii = 1
      do n_r=1,n_r_max_3D
         do lm=lmStart,lmStop
            arr_LMloc(lm,n_r)=rbuff(ii)
            ii = ii +1
         end do
      end do

      deallocate( rbuff, sbuff, scounts, sdisp)

   end subroutine scatter_from_rank0_to_lmloc
!------------------------------------------------------------------------------
   subroutine reduce_radial_on_rank(arr_dist, len_arr_r, irank)

      !-- Input variable
      integer,  intent(in) :: len_arr_r
      integer,  intent(in) :: irank

      !-- Output variable
      real(cp), intent(inout) :: arr_dist(len_arr_r)

      !-- Local variable
      integer :: n_r
      real(cp) :: work(len_arr_r)

      call MPI_Reduce(arr_dist, work, len_arr_r, MPI_DEF_REAL, &
           &          MPI_SUM, irank, MPI_COMM_WORLD, ierr)

      if ( rank == irank ) then
         do n_r=1,len_arr_r
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
end module communications
