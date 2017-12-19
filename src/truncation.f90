module truncation

   use precision_mod
   use mem_alloc, only: bytes_allocated

   implicit none

   private

   integer, public :: n_r_max
   integer, public :: n_cheb_max
   integer, public :: m_max
   integer, public :: minc
   integer, public :: n_m_max
   integer, public :: n_phi_max
   integer, public :: n_phi_tot
   integer, allocatable, public :: idx2m(:)
   integer, allocatable, public :: m2idx(:)

   public :: initialize_truncation, finalize_truncation

contains

   subroutine initialize_truncation

      integer :: n_m, m

      n_m_max = m_max/minc+1

      n_phi_tot = 3*m_max!+1
      n_phi_max = n_phi_tot/minc

      allocate( idx2m(n_m_max) )
      allocate( m2idx(0:m_max) )

      bytes_allocated = bytes_allocated+(n_m_max+m_max+1)*SIZEOF_INTEGER

      do n_m=1,n_m_max
         idx2m(n_m)=(n_m-1)*minc
      end do


      m2idx(:)=-1
      n_m = 1
      do m=0,m_max,minc
         m2idx(m)=n_m
         n_m = n_m+1
      end do

   end subroutine initialize_truncation
!------------------------------------------------------------------------------
   subroutine finalize_truncation

      deallocate( m2idx, idx2m)

   end subroutine finalize_truncation
!------------------------------------------------------------------------------
end module truncation
