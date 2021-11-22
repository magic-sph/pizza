#define M_LOOP(mStart, mStop, action) do n_m=mStart,mStop; action; end do;
module parallel_solvers
   !
   ! This module contains the routines that are used to solve linear banded problems
   ! with R-distributed arrays.
   !

   use precision_mod
   use parallel_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: one, zero
   use truncation, only: n_m_max, n_r_max

   implicit none

   private

   type, public :: type_tri_par   
      integer :: nRMin
      integer :: nRMax
      integer :: mMin 
      integer :: mMax
      real(cp), allocatable :: low(:,:)
      real(cp), allocatable :: diag(:,:)
      real(cp), allocatable :: up(:,:)
   contains
      procedure :: initialize => initialize_3
      procedure :: finalize => finalize_3
      procedure :: prepare_mat => prepare_mat_3
      procedure :: solver_up => solver_up_3
      procedure :: solver_dn => solver_dn_3
      procedure :: solver_single ! Used for one single right hand side
      procedure :: solver_finish => solver_finish_3
   end type

   type, public :: type_penta_par
      integer :: nRMin
      integer :: nRMax
      integer :: mMin
      integer :: mMax
      complex(cp), allocatable :: low2(:,:)
      complex(cp), allocatable :: low1(:,:)
      complex(cp), allocatable :: diag(:,:)
      complex(cp), allocatable :: up1(:,:)
      complex(cp), allocatable :: up2(:,:)
   contains
      procedure :: initialize => initialize_5
      procedure :: finalize => finalize_5
      procedure :: prepare_mat => prepare_mat_5
      procedure :: solver_up => solver_up_5
      procedure :: solver_dn => solver_dn_5
      procedure :: solver_finish => solver_finish_5
   end type

contains

   subroutine initialize_3(this, nRstart, nRstop, mMin, mMax)
      !
      ! Memory allocation of a parallel tridiagonal matrix
      !

      class(type_tri_par) :: this
      integer, intent(in) :: nRstart
      integer, intent(in) :: nRstop
      integer, intent(in) :: mMin
      integer, intent(in) :: mMax

      this%nRMin=nRstart
      this%nRMax=nRstop
      this%mMin=mMin
      this%mMax=mMax
      allocate( this%low(mMin:mMax,nRstart-1:nRstop+1) )
      allocate( this%diag(mMin:mMax,nRstart-1:nRstop+1) )
      allocate( this%up(mMin:mMax,nRstart-1:nRstop+1) )
      bytes_allocated = bytes_allocated+3*(mMax-mMin+1)*(nRstop-nRstart+3)*SIZEOF_DEF_REAL

      !-- Fill an identity matrix by default
      this%low(:,:) =0.0_cp
      this%diag(:,:)=one
      this%up(:,:)  =0.0_cp

   end subroutine initialize_3
!-------------------------------------------------------------------------------------
   subroutine initialize_5(this, nRstart, nRstop, mMin, mMax)
      !
      ! Memory allocation of a parallel tridiagonal matrix
      !

      class(type_penta_par) :: this
      integer, intent(in) :: nRstart
      integer, intent(in) :: nRstop
      integer, intent(in) :: mMin
      integer, intent(in) :: mMax

      this%nRMin=nRstart
      this%nRMax=nRstop
      this%mMin=mMin
      this%mMax=mMax
      allocate( this%low1(mMin:mMax,nRstart-2:nRstop+2) )
      allocate( this%low2(mMin:mMax,nRstart-2:nRstop+2) )
      allocate( this%diag(mMin:mMax,nRstart-2:nRstop+2) )
      allocate( this%up1(mMin:mMax,nRstart-2:nRstop+2) )
      allocate( this%up2(mMin:mMax,nRstart-2:nRstop+2) )
      bytes_allocated = bytes_allocated+5*(mMax-mMin+1)*(nRstop-nRstart+5)*SIZEOF_DEF_REAL

      !-- Fill an identity matrix by default
      this%low2(:,:)=zero
      this%low1(:,:)=zero
      this%diag(:,:)=one
      this%up1(:,:) =zero
      this%up2(:,:) =zero

   end subroutine initialize_5
!-------------------------------------------------------------------------------------
   subroutine finalize_3(this)
      !
      ! Memory deallocation of a parallel tridiagonal solver
      !
      class(type_tri_par) :: this

      deallocate(this%low, this%diag, this%up)

   end subroutine finalize_3
!-------------------------------------------------------------------------------------
   subroutine finalize_5(this)
      !
      ! Memory deallocation of a parallel pentadiagonal solver
      !
      class(type_penta_par) :: this

      deallocate(this%low1, this%low2, this%diag, this%up1, this%up2)

   end subroutine finalize_5
!-------------------------------------------------------------------------------------
   subroutine prepare_mat_3(this)
      !
      ! LU factorisation of a tridiagonal matrix: the diagonal is overwritten
      !
      class(type_tri_par) :: this

      !-- Local variables
      real(cp) :: p
      integer :: m, nR, n_ms, tag

      tag = 97437
      n_ms = this%mMax-this%mMin+1

      !-- Set 'out-of-bound' values to zero for safety
      do m=this%mMin, this%mMax
         if ( this%nRMin == 1 ) this%low(m,this%nRMin)=0.0_cp
         if ( this%nRMax == n_r_max ) this%up(m,this%nRMax) =0.0_cp
      end do

      if ( rank == 0 ) then ! Rank 0 starts the flow
         do nR=this%nRMin+1,this%nRMax
            do m=this%mMin,this%mMax
               p=this%diag(m,nR)-this%low(m,nR)*this%up(m,nR-1)* &
               & this%diag(m,nR-1)
               this%diag(m,nR)=one/p
            end do
         end do
      else
         call MPI_Recv(this%diag(this%mMin:this%mMax,this%nRMin-1), n_ms, MPI_DEF_REAL, &
              &        rank-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         call MPI_Recv(this%up(this%mMin:this%mMax,this%nRMin-1), n_ms, MPI_DEF_REAL, &
              &        rank-1, tag+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         do nR=this%nRMin,this%nRMax
            do m=this%mMin,this%mMax
               p=this%diag(m,nR)-this%low(m,nR)*this%up(m,nR-1)* &
               & this%diag(m,nR-1)
               this%diag(m,nR)=one/p
            end do
         end do
      end if
      
      if ( rank < n_procs-1 ) then ! If this is not the last rank
         call MPI_Send(this%diag(this%mMin:this%mMax,this%nRMax), n_ms, MPI_DEF_REAL, &
              &        rank+1, tag, MPI_COMM_WORLD, ierr)
         call MPI_Send(this%up(this%mMin:this%mMax,this%nRMax), n_ms, MPI_DEF_REAL, &
              &        rank+1, tag+1, MPI_COMM_WORLD, ierr)
      end if

   end subroutine prepare_mat_3
!-------------------------------------------------------------------------------------
   subroutine prepare_mat_5(this)
      !
      ! LU factorisation of a pentadiagonal matrix
      !
      class(type_penta_par) :: this

      !-- Local variables
      integer :: nR, m, n_ms, tag

      n_ms = this%mMax-this%mMin+1
      tag = 69843

      !-- Set 'out-of-bound' values to zero for safety
      do m=this%mMin,this%mMax
         if ( this%nRMin == 1 ) then
            this%low1(m,this%nRMin)  =0.0_cp
            this%low2(m,this%nRMin)  =0.0_cp
            this%low2(m,this%nRMin+1)=0.0_cp
         end if
         if ( this%nRMax == n_r_max ) then
            this%up1(m,this%nRMax)  =0.0_cp
            this%up2(m,this%nRMax)  =0.0_cp
            this%up2(m,this%nRMax-1)=0.0_cp
         end if
      end do

      !-- Now proper LU factorisation
      if ( rank == 0 ) then
         nR=2
         do m=this%mMin,this%mMax
            this%up1(m,nR)=this%up1(m,nR)-this%low1(m,nR)*this%up2(m,nR-1)/ &
            &              this%diag(m,nR-1)
            this%diag(m,nR)=this%diag(m,nR)-this%low1(m,nR)*this%up1(m,nR-1)/ &
            &               this%diag(m,nR-1)
         end do

         do nR=this%nRmin+2,this%nRMax
            do m=this%mMin,this%mMax
               this%low1(m,nR)=this%low1(m,nR)-this%low2(m,nR)*this%up1(m,nR-2)/ &
               &               this%diag(m,nR-2)
               this%up1(m,nR)=this%up1(m,nR)-this%low1(m,nR)*this%up2(m,nR-1)/ &
               &              this%diag(m,nR-1)
               this%diag(m,nR)=this%diag(m,nR)-this%low1(m,nR)*this%up1(m,nR-1)/ &
               &               this%diag(m,nR-1)-this%low2(m,nR)*this%up2(m,nR-2)/ &
               &               this%diag(m,nR-2)
             end do
         end do
      else
         call MPI_Recv(this%diag(this%mMin:this%mMax,this%nRMin-2), n_ms, MPI_DEF_COMPLEX, &
              &        rank-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         call MPI_Recv(this%diag(this%mMin:this%mMax,this%nRMin-1), n_ms, MPI_DEF_COMPLEX, &
              &        rank-1, tag+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         call MPI_Recv(this%up1(this%mMin:this%mMax,this%nRMin-2), n_ms, MPI_DEF_COMPLEX, &
              &        rank-1, tag+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         call MPI_Recv(this%up1(this%mMin:this%mMax,this%nRMin-1), n_ms, MPI_DEF_COMPLEX, &
              &        rank-1, tag+2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         call MPI_Recv(this%up2(this%mMin:this%mMax,this%nRMin-2), n_ms, MPI_DEF_COMPLEX, &
              &        rank-1, tag+3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         call MPI_Recv(this%up2(this%mMin:this%mMax,this%nRMin-1), n_ms, MPI_DEF_COMPLEX, &
              &        rank-1, tag+4, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

         do nR=this%nRmin,this%nRMax
            do m=this%mMin,this%mMax
               this%low1(m,nR)=this%low1(m,nR)-this%low2(m,nR)*this%up1(m,nR-2)/ &
               &               this%diag(m,nR-2)
               this%up1(m,nR)=this%up1(m,nR)-this%low1(m,nR)*this%up2(m,nR-1)/ &
               &              this%diag(m,nR-1)
               this%diag(m,nR)=this%diag(m,nR)-this%low1(m,nR)*this%up1(m,nR-1)/ &
               &               this%diag(m,nR-1)-this%low2(m,nR)*this%up2(m,nR-2)/ &
               &               this%diag(m,nR-2)
             end do
         end do
      end if

      if ( rank < n_procs-1 ) then ! If not last rank -> send
         call MPI_Send(this%diag(this%mMin:this%mMax,this%nRMax-1), n_ms, MPI_DEF_COMPLEX, &
              &        rank+1, tag, MPI_COMM_WORLD, ierr)
         call MPI_Send(this%diag(this%mMin:this%mMax,this%nRMax), n_ms, MPI_DEF_COMPLEX, &
              &        rank+1, tag+1, MPI_COMM_WORLD, ierr)
         call MPI_Send(this%up1(this%mMin:this%mMax,this%nRMax-1), n_ms, MPI_DEF_COMPLEX, &
              &        rank+1, tag+1, MPI_COMM_WORLD, ierr)
         call MPI_Send(this%up1(this%mMin:this%mMax,this%nRMax), n_ms, MPI_DEF_COMPLEX, &
              &        rank+1, tag+2, MPI_COMM_WORLD, ierr)
         call MPI_Send(this%up2(this%mMin:this%mMax,this%nRMax-1), n_ms, MPI_DEF_COMPLEX, &
              &        rank+1, tag+3, MPI_COMM_WORLD, ierr)
         call MPI_Send(this%up2(this%mMin:this%mMax,this%nRMax), n_ms, MPI_DEF_COMPLEX, &
              &        rank+1, tag+4, MPI_COMM_WORLD, ierr)
      end if

      do nR=this%nRmin,this%nRMax
         do m=this%mMin,this%mMax
            this%diag(m,nR)=one/this%diag(m,nR)
            this%up1(m,nR) =this%up1(m,nR)*this%diag(m,nR)
            this%up2(m,nR) =this%up2(m,nR)*this%diag(m,nR)
            this%low1(m,nR)=this%low1(m,nR)*this%diag(m,nR)
            this%low2(m,nR)=this%low2(m,nR)*this%diag(m,nR)
          end do
      enddo

   end subroutine prepare_mat_5
!-------------------------------------------------------------------------------------
   subroutine solver_single(this, x, nRstart, nRstop)
      !
      ! This routine is used to solve a solve one single linear system that does not
      ! depend on lm. This is for instance used for z10 when l_z10mat is required.
      !
      class(type_tri_par) :: this

      !-- Input variables
      integer, intent(in) :: nRstart, nRstop

      !-- Output variables
      complex(cp), intent(inout) :: x(nRstart-1:nRstop+1)

      !-- Local variables
      integer :: nR0,nR
      integer :: tag, n_r_cmb, n_r_icb

      n_r_cmb=1
      n_r_icb=n_r_max
      tag = 53976

      nR0 = nRstart
      if ( nRstart > n_r_cmb ) then ! Not the first block
         nR0 = nR0-1
         call MPI_Recv(x(nR0), 1, MPI_DEF_COMPLEX, rank-1, tag, MPI_COMM_WORLD, &
              &        MPI_STATUS_IGNORE, ierr)
      else ! Lower boundary: x -> x - low * x(i-1)
         x(nR0)=x(nR0)-this%low(1,nR0)*x(nR0-1)
      end if

      do nR=nR0+1,nRstop
         x(nR)=x(nR)-this%diag(1,nR-1)*this%low(1,nR)*x(nR-1)
      end do

      if ( nRstop < n_r_icb ) then ! Not the last block
         call MPI_SSend(x(nRstop), 1, MPI_DEF_COMPLEX, rank+1, tag, MPI_COMM_WORLD, ierr)
      end if

      tag = tag+1
      if ( nRstop < n_r_icb ) then ! This is not the last chunk
         call MPI_Recv(x(nRstop+1), 1, MPI_DEF_COMPLEX, rank+1, &
              &        tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if

      do nR=nRstop,nRstart,-1
         x(nR)=(x(nR)-this%up(1,nR)*x(nR+1))*this%diag(1,nR)
      end do

      if ( nRstart > n_r_cmb ) then
         call MPI_SSend(x(nRstart), 1, MPI_DEF_COMPLEX, rank-1, tag, MPI_COMM_WORLD, ierr)
      end if

      tag = tag+1
      if ( nRstart /= n_r_cmb ) then
         if ( nRstart == n_r_cmb+1 ) then ! send down
            call MPI_Ssend(x(nRstart), 1, MPI_DEF_COMPLEX, rank-1, tag, &
                 &         MPI_COMM_WORLD, ierr)
         end if

         if ( nRstop == n_r_cmb ) then ! send down
            call MPI_Recv(x(nRstop+1), 1, MPI_DEF_COMPLEX, rank+1, tag, &
                 &        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         end if
      end if

      tag = tag+1
      if ( nRstart > n_r_cmb ) then ! This is not the first block
         call MPI_Recv(x(nRstart-1), 1, MPI_DEF_COMPLEX, rank-1, &
              &         tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if

      if ( nRstop < n_r_icb .and. nRstop >= n_r_cmb ) then ! This is not the last block
         call MPI_Ssend(x(nRstop), 1, MPI_DEF_COMPLEX, rank+1, tag, MPI_COMM_WORLD, ierr)
      end if

   end subroutine solver_single
!-------------------------------------------------------------------------------------
   subroutine solver_up_3(this, x, mStart, mStop, nRstart, nRstop, tag, array_req, &
              &           req, ms_block, n_m_block)
      !
      ! First part of the parallel tridiag solver: forward substitution
      !

      class(type_tri_par) :: this

      !-- Input variables
      integer, intent(in) :: mStart ! Starting m (OMP thread dependent)
      integer, intent(in) :: mStop  ! Stopping m (OMP thread dependent)
      integer, intent(in) :: nRstart ! Starting nR
      integer, intent(in) :: nRstop  ! Stopping nR
      integer, intent(in) :: ms_block ! Starting block-index of m
      integer, intent(in) :: n_m_block ! Size of the block
      integer, intent(in) :: tag

      !-- Output variables
      integer,     intent(inout) :: array_req(:)
      complex(cp), intent(inout) :: x(1:n_m_max, nRstart-1:nRstop+1)
      integer,     intent(inout) :: req

      !-- Local variables
      integer :: nR, nR0, n_m, mb, mu

      mb = ms_block
      mu = mb+n_m_block-1

      nR0 = nRstart
      if ( nRstart > 1 ) then ! Not the first block
         nR0 = nR0-1
         !$omp master
         call MPI_Recv(x(mb:mu,nR0), mu-mb+1, MPI_DEF_COMPLEX, rank-1, tag, &
              &        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         !$omp end master
         !$omp barrier
      else ! Lower boundary: x -> x - low * x(i-1)
         M_LOOP(mStart, mStop, x(n_m,nR0)=x(n_m,nR0)-this%low(n_m,nR0)*x(n_m,nR0-1))
         !$omp barrier
      end if

      do nR=nR0+1,nRstop
         M_LOOP(mStart, mStop, x(n_m,nR)=x(n_m,nR)-this%diag(n_m,nR-1)*this%low(n_m,nR)*x(n_m,nR-1))
      end do
      !$omp barrier

      if ( nRstop < n_r_max ) then ! Not the last block
         !$omp barrier
         !$omp master
         call MPI_ISend(x(mb:mu,nRstop), mu-mb+1, MPI_DEF_COMPLEX, rank+1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         req = req+1
         !call MPI_SSend(x(mb:mu,nRstop), mu-mb+1, MPI_DEF_COMPLEX, rank+1, &
         !     &         tag, MPI_COMM_WORLD, ierr)
         !$omp end master
      end if

   end subroutine solver_up_3
!-------------------------------------------------------------------------------------
   subroutine solver_up_5(this, x, mStart, mStop, nRstart, nRstop, tag, array_req, &
              &           req, ms_block, n_m_block)
      !
      ! First part of the parallel pentadiag solver: forward substitution
      !
      class(type_penta_par) :: this

      !-- Input variables
      integer, intent(in) :: mStart ! Starting m (OMP thread dependent)
      integer, intent(in) :: mStop  ! Stopping m (OMP thread dependent)
      integer, intent(in) :: nRstart ! Starting nR
      integer, intent(in) :: nRstop  ! Stopping nR
      integer, intent(in) :: ms_block ! Starting block-index of m
      integer, intent(in) :: n_m_block ! Size of the block
      integer, intent(in) :: tag

      !-- Output variables
      integer,     intent(inout) :: array_req(:)
      complex(cp), intent(inout) :: x(1:n_m_max, nRstart-2:nRstop+2)
      integer,     intent(inout) :: req

      !-- Local variables
      integer :: nR, n_m, recv(2), lb, lu, mb, mu

      mb = ms_block
      mu = mb+n_m_block-1
      lb=mStart
      lu=mStop

      if ( nRstart > 1 ) then ! This is not the first block
         !$omp master
         call MPI_IRecv(x(mb:mu,nRstart-2), mu-mb+1, MPI_DEF_COMPLEX, rank-1,&
              &        tag, MPI_COMM_WORLD, recv(1), ierr)
         !-- Non blocking receive
         call MPI_IRecv(x(mb:mu,nRstart-1), mu-mb+1, MPI_DEF_COMPLEX, rank-1,&
              &        tag+1, MPI_COMM_WORLD, recv(2), ierr)
         !-- Non-blocking receive
         call MPI_Waitall(2, recv, MPI_STATUSES_IGNORE, ierr) ! wait to receive
         !$omp end master
         !$omp barrier
      end if

      do nR=nRstart,nRstop
         M_LOOP(lb,lu,x(n_m,nR)=this%diag(n_m,nR)*x(n_m,nR)-this%low1(n_m,nR)*x(n_m,nR-1)-this%low2(n_m,nR)*x(n_m,nR-2))
      end do
      !$omp barrier

      if ( nRstop < n_r_max ) then ! This is not the last block
         !$omp barrier
         !$omp master
         call MPI_ISend(x(mb:mu,nRstop-1), mu-mb+1, MPI_DEF_COMPLEX, rank+1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         call MPI_ISend(x(mb:mu,nRstop), mu-mb+1, MPI_DEF_COMPLEX, rank+1, &
              &         tag+1, MPI_COMM_WORLD, array_req(req+1), ierr)
         req = req+2
         !$omp end master
      end if

   end subroutine solver_up_5
!-------------------------------------------------------------------------------------
   subroutine solver_dn_3(this, x, mStart, mStop, nRstart, nRstop, tag, array_req, &
              &           req, ms_block, n_m_block)
      !
      ! Second part of the parallel tridiag solver: backward substitution
      !

      class(type_tri_par) :: this

      !-- Input variables
      integer, intent(in) :: mStart ! Starting m (OMP thread dependent)
      integer, intent(in) :: mStop  ! Stopping m (OMP thread dependent)
      integer, intent(in) :: nRstart ! Starting nR
      integer, intent(in) :: nRstop  ! Stopping nR
      integer, intent(in) :: ms_block ! Starting block-index of m
      integer, intent(in) :: n_m_block ! Size of the block
      integer, intent(in) :: tag

      !-- Output variables
      integer,     intent(inout) :: array_req(:)
      complex(cp), intent(inout) :: x(1:n_m_max, nRstart-1:nRstop+1)
      integer,     intent(inout) :: req

      !-- Local variables
      integer :: nR, n_m, mb, mu

      mb = ms_block
      mu = mb+n_m_block-1

      if ( nRstop < n_r_max ) then ! This is not the last chunk
         !$omp master
         call MPI_Recv(x(mb:mu,nRstop+1), mu-mb+1, MPI_DEF_COMPLEX, rank+1, &
              &        tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         !$omp end master
         !$omp barrier
      end if

      do nR=nRstop,nRstart,-1
         M_LOOP(mStart,mStop,x(n_m,nR)=(x(n_m,nR)-this%up(n_m,nR)*x(n_m,nR+1))*this%diag(n_m,nR))
      end do

      if ( nRstart > 1 ) then
         !$omp barrier
         !$omp master
         call MPI_ISend(x(mb:mu,nRstart), mu-mb+1, MPI_DEF_COMPLEX, rank-1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         req=req+1
         !call MPI_SSend(x(mb:mu,nRstart), mu-mb+1, MPI_DEF_COMPLEX, rank-1, &
         !     &         tag, MPI_COMM_WORLD, ierr)
         !$omp end master
      end if

   end subroutine solver_dn_3
!-------------------------------------------------------------------------------------
   subroutine solver_dn_5(this, x, mStart, mStop, nRstart, nRstop, tag, array_req, &
              &           req, ms_block, n_m_block)
      !
      ! Second part of the parallel pentadiagonal solver: backward substitution
      !

      class(type_penta_par) :: this

      !-- Input variables
      integer, intent(in) :: mStart ! Starting m (OMP thread dependent)
      integer, intent(in) :: mStop  ! Stopping m (OMP thread dependent)
      integer, intent(in) :: nRstart ! Starting nR
      integer, intent(in) :: nRstop  ! Stopping nR
      integer, intent(in) :: ms_block ! Starting block-index of m
      integer, intent(in) :: n_m_block ! Size of the block
      integer, intent(in) :: tag

      !-- Output variables
      integer,     intent(inout) :: array_req(:)
      complex(cp), intent(inout) :: x(n_m_max, nRstart-2:nRstop+2)
      integer,     intent(inout) :: req

      !-- Local variables
      integer :: nR, n_m, rcv(2), mb, mu

      mb = ms_block
      mu = mb+n_m_block-1

      if ( nRstop < n_r_max ) then ! This is not the last rank
         !$omp master
         call MPI_IRecv(x(mb:mu,nRstop+1), mu-mb+1, MPI_DEF_COMPLEX, rank+1, &
              &        tag, MPI_COMM_WORLD, rcv(1), ierr)
         call MPI_IRecv(x(mb:mu,nRstop+2), mu-mb+1, MPI_DEF_COMPLEX, rank+1, &
              &        tag+1, MPI_COMM_WORLD, rcv(2), ierr)
         call MPI_Waitall(2, rcv, MPI_STATUSES_IGNORE, ierr)
         !$omp end master
         !$omp barrier
      end if

      do nR=nRstop,nRstart,-1
         M_LOOP(mStart,mStop,x(n_m,nR)=x(n_m,nR)-this%up1(n_m,nR)*x(n_m,nR+1)-this%up2(n_m,nR)*x(n_m,nR+2))
      end do

      !$omp barrier
      !$omp master
      if ( nRstart > 1 ) then
         call MPI_Isend(x(mb:mu,nRstart), mu-mb+1, MPI_DEF_COMPLEX, rank-1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         call MPI_Isend(x(mb:mu,nRstart+1), mu-mb+1, MPI_DEF_COMPLEX, rank-1, &
              &         tag+1, MPI_COMM_WORLD, array_req(req+1), ierr)
         req=req+2
      end if
      !$omp end master

   end subroutine solver_dn_5
!-------------------------------------------------------------------------------------
   subroutine solver_finish_3(this, x, ms_block, n_m_block, nRstart, nRstop, tag, &
              &               array_req, req)
      !
      ! Last part of the parallel tridiag solver: halo synchronisation
      !

      class(type_tri_par) :: this

      !-- Input variables
      integer, intent(in) :: nRstart  ! Starting index in radius
      integer, intent(in) :: nRstop   ! Stopping index in radius
      integer, intent(in) :: ms_block ! Starting block-index of m
      integer, intent(in) :: n_m_block ! Size of the block
      integer, intent(in) :: tag

      !-- Output variables
      integer,     intent(inout) :: array_req(:)
      complex(cp), intent(inout) :: x(1:n_m_max, nRstart-1:nRstop+1)
      integer,     intent(inout) :: req
      
      !-- Local variables
      integer :: mb, mu

      mb = ms_block
      mu = mb+n_m_block-1

      !$omp master
      if ( nRstart /= 1 ) then
         if ( nRstart == 2 ) then ! send down
            call MPI_Isend(x(mb:mu,nRstart), mu-mb+1, MPI_DEF_COMPLEX, &
                 &         rank-1, tag, MPI_COMM_WORLD, array_req(req), ierr)
            req = req+1
         end if

         if ( nRstop == 1 ) then ! send down
            call MPI_Irecv(x(mb:mu,nRstop+1), mu-mb+1, MPI_DEF_COMPLEX, &
                 &         rank+1, tag, MPI_COMM_WORLD, array_req(req), ierr)
            req = req+1
         end if
      end if

      if ( nRstart > 1 ) then ! This is not the first block
         call MPI_Irecv(x(mb:mu,nRstart-1), mu-mb+1, MPI_DEF_COMPLEX, rank-1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         req = req+1
      end if

      if ( nRstop < n_r_max .and. nRstop >= 1 ) then ! This is not the last block
         !call MPI_Ssend(x(mb:mu,nRstop), mu-mb+1, MPI_DEF_COMPLEX, rank+1, &
         !     &         tag, MPI_COMM_WORLD, ierr)
         call MPI_Isend(x(mb:mu,nRstop), mu-mb+1, MPI_DEF_COMPLEX, rank+1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         req = req+1
      end if
      !$omp end master

   end subroutine solver_finish_3
!-------------------------------------------------------------------------------------
   subroutine solver_finish_5(this, x, ms_block, n_m_block, nRstart, nRstop, tag, &
              &               array_req, req)
      !
      ! Last part of the parallel pentadiag solver: halo synchronisation
      !

      class(type_penta_par) :: this

      !-- Input variables
      integer, intent(in) :: nRstart  ! Starting index in radius
      integer, intent(in) :: nRstop   ! Stopping index in radius
      integer, intent(in) :: ms_block ! Starting block-index of m
      integer, intent(in) :: n_m_block ! Size of the block
      integer, intent(in) :: tag

      !-- Output variables
      integer,     intent(inout) :: array_req(:) ! MPI requests
      complex(cp), intent(inout) :: x(n_m_max, nRstart-2:nRstop+2) ! Solution
      integer,     intent(inout) :: req

      !-- Local variables
      integer :: istart, istop, mb, mu
      logical :: l_send_dn, l_send_up

      mb = ms_block
      mu = mb+n_m_block-1
      l_send_dn = .false.
      l_send_up = .false.
      istart = 1
      istop = n_r_max

      !$omp master
      if ( nRstart == 2 ) then ! Send down
         call MPI_Isend(x(mb:mu,nRstart), mu-mb+1, MPI_DEF_COMPLEX, rank-1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         call MPI_Isend(x(mb:mu,nRstart+1), mu-mb+1, MPI_DEF_COMPLEX, rank-1, &
              &         tag+1, MPI_COMM_WORLD, array_req(req+1), ierr)
         req=req+2
      end if
      if ( nRStop==1 ) then ! Receive up
         call MPI_Irecv(x(mb:mu,nRstop+1), mu-mb+1, MPI_DEF_COMPLEX, rank+1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         call MPI_Irecv(x(mb:mu,nRstop+2), mu-mb+1, MPI_DEF_COMPLEX, rank+1, &
              &         tag+1, MPI_COMM_WORLD, array_req(req+1), ierr)
         req=req+2
      end if

      if ( nRstart > 1 ) then ! This is not the first block
         l_send_dn=.true.
         istart=nRstart
      end if
      if ( (nRstop < n_r_max) .and. (nRstop>=istart) ) then
         l_send_up=.true.
         istop=nRstop
      end if

      if ( l_send_up ) then
         !-- Update halo
         call MPI_Isend(x(mb:mu,istop), mu-mb+1, MPI_DEF_COMPLEX, rank+1, &
              &         tag, MPI_COMM_WORLD, array_req(req), ierr)
         req=req+1
      end if
      if ( l_send_dn ) then
         call MPI_Irecv(x(mb:mu,istart-2), mu-mb+1, MPI_DEF_COMPLEX, rank-1, &
              &         tag+1, MPI_COMM_WORLD, array_req(req), ierr)
         req=req+1
         if ( (istop > istart) .or. (.not. l_send_up) ) then
            call MPI_Irecv(x(mb:mu,istart-1), mu-mb+1, MPI_DEF_COMPLEX, &
                 &         rank-1, tag, MPI_COMM_WORLD, array_req(req), ierr)
            req=req+1
         else
            call MPI_Recv(x(mb:mu,istart-1), mu-mb+1, MPI_DEF_COMPLEX, &
                 &        rank-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         end if
      end if
      if ( l_send_up ) then
         call MPI_Isend(x(mb:mu,istop-1), mu-mb+1, MPI_DEF_COMPLEX, rank+1, &
              &         tag+1, MPI_COMM_WORLD, array_req(req), ierr)
         req=req+1
      end if
      !$omp end master

   end subroutine solver_finish_5
!-------------------------------------------------------------------------------------
end module parallel_solvers
