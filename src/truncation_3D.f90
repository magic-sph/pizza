module truncation_3D

   use precision_mod
   use mem_alloc, only: bytes_allocated

   implicit none

   private

   integer, public :: minc_3D
   integer, public :: n_r_max_3D
   integer, public :: n_cheb_max_3D
   !integer, public :: n_z_max
   integer, public :: l_max         ! max degree of Plms
   integer, public :: n_theta_max
   integer, public :: m_max_3D
   integer, public :: n_m_max_3D
   integer, public :: n_phi_tot_3D
   integer, public :: n_phi_max_3D
   integer, public :: lm_max        ! number of l/m combinations
   integer, public :: lmP_max       ! number of l/m combination if l runs to l_max+1
   integer, allocatable, public :: idx2m3D(:)
   integer, allocatable, public :: m3D_2idx(:)

   public :: initialize_truncation_3D, finalize_truncation_3D, &
   &         write_truncation_3D_info

contains

   subroutine initialize_truncation_3D

      !-- Local variables
      integer :: n_m, m

      ! number of points for z-grid in the 3D interpolation
      !n_z_max = 2*n_r_max_3D! 2xn_s -> In truncation.f90

      ! number of theta grid-points & phi 3D grid
      n_phi_max_3D=n_phi_tot_3D/minc_3D
      n_theta_max =n_phi_tot_3D/2

      l_max   =(2*n_theta_max)/3
      m_max_3D=(l_max/minc_3D)*minc_3D

      ! max number of ms (different oders)
      n_m_max_3D=m_max_3D/minc_3D+1

      ! number of l/m combinations
      lm_max=m_max_3D*(l_max+1)/minc_3D - m_max_3D*(m_max_3D-minc_3D)/ &
      &      (2*minc_3D)+(l_max+1-m_max_3D)
      ! number of l/m combination if l runs to l_max+1
      lmP_max=lm_max+n_m_max_3D

      !-- retrieve indices from m and vice-versa functions
      allocate( idx2m3D(n_m_max_3D) )
      allocate( m3D_2idx(0:m_max_3D) )

      bytes_allocated = bytes_allocated+(n_m_max_3D+m_max_3D+1)*SIZEOF_INTEGER

      do n_m=1,n_m_max_3D
         idx2m3D(n_m)=(n_m-1)*minc_3D
      end do

      m3D_2idx(:)=-1
      n_m = 1
      do m=0,m_max_3D,minc_3D
         m3D_2idx(m)=n_m
         n_m = n_m+1
      end do

   end subroutine initialize_truncation_3D
!------------------------------------------------------------------------------
   subroutine finalize_truncation_3D

      deallocate( m3D_2idx, idx2m3D)

   end subroutine finalize_truncation_3D
!------------------------------------------------------------------------------
   subroutine write_truncation_3D_info(n_out)

      integer, intent(in) :: n_out

      write(n_out,*)
      write(n_out,*) '! 3-D Grid parameters:'
      write(n_out,'(''  n_r_max_3D   ='',i6, &
           &   '' = number of radial grid points'')') n_r_max_3D
      write(n_out,'(''  n_cheb_max_3D='',i6)') n_cheb_max_3D
      write(n_out,'(''  n_phi_max_3D ='',i6, &
           &   '' = no of longitude grid points'')') n_phi_max_3D
      write(n_out,'(''  n_theta_max  ='',i6, &
           &   '' = no of latitude grid points'')') n_theta_max
      write(n_out,'(''  l_max        ='',i6, '' = max degree of Plm'')') l_max
      write(n_out,'(''  m_max_3D     ='',i6, '' = max oder of Plm'')') m_max_3D
      write(n_out,'(''  lm_max       ='',i6, '' = no of l/m combinations'')') lm_max
      write(n_out,'(''  minc_3D      ='',i6, '' = longitude symmetry wave no'')') minc_3D

   end subroutine write_truncation_3D_info
!------------------------------------------------------------------------------
end module truncation_3D
