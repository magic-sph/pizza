module mpi_transp_mod
   !
   ! This is an abstract class that will be used to define MPI transposers
   !

   use precision_mod
   use truncation, only: n_m_max, n_r_max
   use blocking, only: nRstart, nRstop, nMstart,nMstop

   implicit none

   private

   type, abstract, public :: type_mpitransp
      integer :: n_fields
   contains
      procedure(create_if), deferred :: create_comm
      procedure(destroy_if), deferred :: destroy_comm
      procedure(transp_m2r_if), deferred :: transp_m2r
      procedure(transp_r2m_if), deferred :: transp_r2m
   end type type_mpitransp

   interface 

      subroutine create_if(this,n_fields)
         import
         class(type_mpitransp) :: this
         integer, intent(in) :: n_fields
      end subroutine create_if

      subroutine destroy_if(this)
         import
         class(type_mpitransp) :: this
      end subroutine destroy_if

      subroutine transp_m2r_if(this, arr_Mloc, arr_Rloc)
         import
         class(type_mpitransp) :: this
         complex(cp), intent(in) :: arr_Mloc(nMstart:nMstop,n_r_max,*)
         complex(cp), intent(out) :: arr_Rloc(n_m_max,nRstart:nRstop,*)
      end subroutine transp_m2r_if

      subroutine transp_r2m_if(this, arr_Rloc, arr_Mloc)
         import
         class(type_mpitransp) :: this
         complex(cp), intent(in) :: arr_Rloc(n_m_max,nRstart:nRstop,*)
         complex(cp), intent(out) :: arr_Mloc(nMstart:nMstop,1:n_r_max,*)
      end subroutine transp_r2m_if

   end interface

end module mpi_transp_mod


!----------------------------------------------------------------------------------
module  mpi_alltoall_mod
   !
   ! This module contains the implementation of all-to-all global communicators
   !

   use precision_mod
   use parallel_mod
   use mem_alloc
   use truncation, only: n_m_max, n_r_max
   use blocking
   use mpi_transp_mod, only: type_mpitransp

   implicit none

   private

   type, public, extends(type_mpitransp) :: type_mpiatoav
      integer, allocatable :: rcounts(:)
      integer, allocatable :: scounts(:)
      integer, allocatable :: rdisp(:)
      integer, allocatable :: sdisp(:)
      integer :: max_send, max_recv
   contains
      procedure :: create_comm => create_comm_alltoallv
      procedure :: destroy_comm => destroy_comm_alltoallv
      procedure :: transp_m2r => transp_m2r_alltoallv
      procedure :: transp_r2m => transp_r2m_alltoallv
   end type type_mpiatoav

   type, public, extends(type_mpitransp) :: type_mpiatoaw
      integer, allocatable :: rtype(:)
      integer, allocatable :: stype(:)
      integer, allocatable :: disp(:)
      integer, allocatable :: counts(:)
   contains
      procedure :: create_comm => create_comm_alltoallw
      procedure :: destroy_comm => destroy_comm_alltoallw
      procedure :: transp_m2r => transp_m2r_alltoallw
      procedure :: transp_r2m => transp_r2m_alltoallw
   end type type_mpiatoaw

contains

   subroutine create_comm_alltoallv(this,n_fields)

      class(type_mpiatoav) :: this

      !-- Input variable
      integer, intent(in) :: n_fields

      !-- Local variables
      integer :: p

      this%n_fields=n_fields

      allocate ( this%rcounts(0:n_procs-1), this%scounts(0:n_procs-1) )
      allocate ( this%rdisp(0:n_procs-1), this%sdisp(0:n_procs-1) )

      do p=0,n_procs-1
         this%scounts(p)=nR_per_rank*m_balance(p)%n_per_rank*this%n_fields
         this%rcounts(p)=radial_balance(p)%n_per_rank*nm_per_rank*this%n_fields
      end do

      this%rdisp(0)=0
      this%sdisp(0)=0
      do p=1,n_procs-1
         this%sdisp(p)=this%sdisp(p-1)+this%scounts(p-1)
         this%rdisp(p)=this%rdisp(p-1)+this%rcounts(p-1)
      end do

      this%max_send = sum(this%scounts)
      this%max_recv = sum(this%rcounts)

      bytes_allocated = bytes_allocated+4*n_procs*SIZEOF_INTEGER

   end subroutine create_comm_alltoallv
!----------------------------------------------------------------------------------
   subroutine create_comm_alltoallw(this,n_fields)

      class(type_mpiatoaw) :: this

      !-- Input variable
      integer, intent(in) :: n_fields

      !-- Local variables
      integer :: arr_size(3), arr_loc_size(3), arr_start(3)
      integer :: p, my_m_counts

      this%n_fields = n_fields

      allocate ( this%counts(0:n_procs-1), this%disp(0:n_procs-1) )
      allocate ( this%rtype(0:n_procs-1), this%stype(0:n_procs-1) )
      bytes_allocated = bytes_allocated+4*n_procs*SIZEOF_INTEGER

      do p=0,n_procs-1
         my_m_counts = m_balance(p)%n_per_rank

         this%counts(p)=1
         this%disp(p)  =0

         arr_size(1)=n_m_max
         arr_size(2)=nR_per_rank
         arr_size(3)=this%n_fields
         arr_loc_size(1)=my_m_counts
         arr_loc_size(2)=nR_per_rank
         arr_loc_size(3)=this%n_fields
         arr_start(1)=m_balance(p)%nStart-1
         arr_start(2)=0
         arr_start(3)=0
         call MPI_Type_Create_Subarray(3, arr_size, arr_loc_size, arr_start, &
              &                        MPI_ORDER_FORTRAN, MPI_DEF_COMPLEX,   &
              &                        this%stype(p), ierr)
         call MPI_Type_Commit(this%stype(p), ierr)

         arr_size(1)=nm_per_rank
         arr_size(2)=n_r_max
         arr_size(3)=this%n_fields
         arr_loc_size(1)=nm_per_rank
         arr_loc_size(2)=radial_balance(p)%n_per_rank
         arr_loc_size(3)=this%n_fields
         arr_start(1)=0
         arr_start(2)=radial_balance(p)%nStart-1
         arr_start(3)=0
         call MPI_Type_Create_Subarray(3, arr_size, arr_loc_size, arr_start, &
              &                        MPI_ORDER_FORTRAN, MPI_DEF_COMPLEX,   &
              &                        this%rtype(p), ierr)
         call MPI_Type_Commit(this%rtype(p), ierr)
      end do

   end subroutine create_comm_alltoallw
!----------------------------------------------------------------------------------
   subroutine destroy_comm_alltoallv(this)

      class(type_mpiatoav) :: this

      deallocate( this%sdisp, this%rdisp, this%scounts, this%rcounts )

   end subroutine destroy_comm_alltoallv
!----------------------------------------------------------------------------------
   subroutine destroy_comm_alltoallw(this)

      class(type_mpiatoaw) :: this

      !-- Local variables
      integer :: p

      do p = 0, n_procs-1
         call MPI_Type_Free(this%rtype(p), ierr)
      end do
      deallocate( this%counts, this%disp, this%rtype, this%stype )

   end subroutine destroy_comm_alltoallw
!----------------------------------------------------------------------------------
   subroutine transp_m2r_alltoallv(this, arr_Mloc, arr_Rloc)

      class(type_mpiatoav) :: this

      !-- Input variables
      complex(cp), intent(in) :: arr_Mloc(nMstart:nMstop,1:n_r_max,*)

      !-- Output variable
      complex(cp), intent(out) :: arr_Rloc(1:n_m_max,nRstart:nRstop,*)

      !-- Local variables
      complex(cp) :: sbuff(1:this%max_send),rbuff(1:this%max_recv)
      integer :: p, ii, n_r, n_m, n_f

      do p = 0, n_procs-1
         ii = this%rdisp(p)+1
         do n_f=1,this%n_fields
            do n_r=radial_balance(p)%nStart,radial_balance(p)%nStop
               do n_m=nMstart,nMstop
                  rbuff(ii)=arr_Mloc(n_m,n_r,n_f)
                  ii = ii+1
               end do
            end do
         end do
      end do

      call MPI_Alltoallv(rbuff, this%rcounts, this%rdisp, MPI_DEF_COMPLEX, &
           &             sbuff, this%scounts, this%sdisp, MPI_DEF_COMPLEX, &
           &             MPI_COMM_WORLD, ierr) 

      do p = 0, n_procs-1
         ii = this%sdisp(p)+1
         do n_f=1,this%n_fields
            do n_r=nRstart,nRstop
               do n_m=m_balance(p)%nStart,m_balance(p)%nStop
                  arr_Rloc(n_m,n_r,n_f)=sbuff(ii)
                  ii=ii+1
               end do
            end do
         end do
      end do

   end subroutine transp_m2r_alltoallv
!----------------------------------------------------------------------------------
   subroutine transp_m2r_alltoallw(this, arr_Mloc, arr_Rloc)

      class(type_mpiatoaw) :: this
      complex(cp), intent(in) :: arr_Mloc(nMstart:nMstop,n_r_max,*)
      complex(cp), intent(out) :: arr_Rloc(n_m_max,nRstart:nRstop,*)

      call MPI_Alltoallw(arr_Mloc, this%counts, this%disp, this%rtype, &
           &             arr_Rloc, this%counts, this%disp, this%stype, &
           &             MPI_COMM_WORLD, ierr)

   end subroutine transp_m2r_alltoallw
!----------------------------------------------------------------------------------
   subroutine transp_r2m_alltoallv(this, arr_Rloc, arr_Mloc)

      !-- Input variables
      class(type_mpiatoav) :: this
      complex(cp), intent(in) :: arr_Rloc(n_m_max,nRstart:nRstop,*)

      !-- Output variable
      complex(cp), intent(out) :: arr_Mloc(nMstart:nMstop,n_r_max,*)

      !-- Local variables
      complex(cp) :: sbuff(1:this%max_send),rbuff(1:this%max_recv)
      integer :: p, ii, n_r, n_m,n_f

      do p = 0, n_procs-1
         ii = this%sdisp(p)+1
         do n_f=1,this%n_fields
            do n_r=nRstart,nRstop
               do n_m=m_balance(p)%nStart,m_balance(p)%nStop
                  sbuff(ii)=arr_Rloc(n_m,n_r,n_f)
                  ii = ii +1
               end do
            end do
         end do
      end do

      call MPI_Alltoallv(sbuff, this%scounts, this%sdisp, MPI_DEF_COMPLEX, &
           &             rbuff, this%rcounts, this%rdisp, MPI_DEF_COMPLEX, &
           &             MPI_COMM_WORLD, ierr) 

      do p = 0, n_procs-1
         ii = this%rdisp(p)+1
         do n_f=1,this%n_fields
            do n_r=radial_balance(p)%nStart,radial_balance(p)%nStop
               do n_m=nMstart,nMstop
                  arr_Mloc(n_m,n_r,n_f)=rbuff(ii)
                  ii=ii+1
               end do
            end do
         end do
      end do

   end subroutine transp_r2m_alltoallv
!----------------------------------------------------------------------------------
   subroutine transp_r2m_alltoallw(this, arr_Rloc, arr_Mloc)

      class(type_mpiatoaw) :: this
      complex(cp), intent(in) :: arr_Rloc(n_m_max,nRstart:nRstop,*)
      complex(cp), intent(out) :: arr_Mloc(nMstart:nMstop,n_r_max,*)

      call MPI_Alltoallw(arr_Rloc, this%counts, this%disp, this%stype, &
           &             arr_Mloc, this%counts, this%disp, this%rtype, &
           &             MPI_COMM_WORLD, ierr)

   end subroutine transp_r2m_alltoallw
!----------------------------------------------------------------------------------
end module mpi_alltoall_mod
