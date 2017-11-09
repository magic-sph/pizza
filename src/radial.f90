module radial_functions
   !
   ! This module calculates some arrays that depend on the radial direction
   ! only. This is called only once at the initialisation of the code.
   !

   use truncation, only: n_r_max, n_cheb_max, m_max
   use radial_der, only: get_dr
   use constants, only: one, two, three, pi, half
   use namelists, only: tag, alph1, alph2, l_newmap, radratio, g0, g1, g2, &
       &                l_non_rot, ek, l_ek_pump, l_temp_3D, tcond_fac,    &
       &                r_cmb, r_icb, l_cheb_coll
   use mem_alloc, only: bytes_allocated
   use radial_scheme, only: type_rscheme
   use chebyshev, only: type_cheb
   use parallel_mod
   use precision_mod

   implicit none

   private
 
   !-- arrays depending on r:
   real(cp), public, allocatable :: r(:)         ! radii
   real(cp), public, allocatable :: or1(:)       ! :math:`1/r`
   real(cp), public, allocatable :: or2(:)       ! :math:`1/r^2`
   real(cp), public, allocatable :: beta(:)      ! 1/h dh/ds
   real(cp), public, allocatable :: dbeta(:)     ! Radial gradient of beta
   real(cp), public, allocatable :: height(:)    ! Spherical shell height
   real(cp), public, allocatable :: ekpump(:)    ! Ekman pumping function
   real(cp), public, allocatable :: oheight(:)   ! 1/h

   real(cp), public, allocatable :: delxr2(:) ! Auxiliary arrays containing effective Courant grid intervals
   real(cp), public, allocatable :: delxh2(:) ! Auxiliary arrays containing effective Courant grid intervals

   real(cp), public, allocatable :: tcond(:) ! Conducting temperature
   real(cp), public, allocatable :: dtcond(:) ! Conducting temperature


   real(cp), public :: alpha1   ! Input parameter for non-linear map to define degree of spacing (0.0:2.0)
   real(cp), public :: alpha2   ! Input parameter for non-linear map to define central point of different spacing (-1.0:1.0)
 
   !-- arrays for buoyancy, depend on Ra and Pr:
   real(cp), public, allocatable :: rgrav(:)     ! Buoyancy term `dtemp0/Di`

   !-- Radial scheme
   class(type_rscheme), public, pointer :: rscheme
 
   public :: initialize_radial_functions, radial, finalize_radial_functions

contains

   subroutine initialize_radial_functions
      !
      ! Initial memory allocation
      !

      integer :: n_in, n_in_2

      ! allocate the arrays
      allocate( r(n_r_max), or1(n_r_max), or2(n_r_max) )
      allocate( beta(n_r_max), dbeta(n_r_max), height(n_r_max) )
      allocate( rgrav(n_r_max), ekpump(n_r_max), oheight(n_r_max) )
      allocate( delxr2(n_r_max), delxh2(n_r_max) )
      allocate( tcond(n_r_max), dtcond(n_r_max))
      bytes_allocated = bytes_allocated+13*n_r_max*SIZEOF_DEF_REAL

      allocate ( type_cheb :: rscheme )

      n_in = n_cheb_max
      if ( l_newmap ) then
         n_in_2 = 1
      else
         n_in_2 = 0
      end if

      call rscheme%initialize(n_r_max,n_in,n_in_2,l_cheb_coll)

   end subroutine initialize_radial_functions
!------------------------------------------------------------------------------
   subroutine finalize_radial_functions

      call rscheme%finalize()

      deallocate( tcond, dtcond )
      deallocate( delxr2, delxh2 )
      deallocate( rgrav )
      deallocate( beta, dbeta, height, ekpump, oheight )
      deallocate( r, or1, or2 )

   end subroutine finalize_radial_functions
!------------------------------------------------------------------------------
   subroutine radial
      !
      !  Calculates everything needed for radial functions, transforms etc.
      !

      integer :: n_r, file_handle
      character(len=100) :: file_name
      real(cp) :: ratio1, ratio2, delmin, c1
      real(cp) :: ek_pump_fac, h

      if ( l_ek_pump .and. (.not. l_non_rot)) then
         ek_pump_fac = one
      else
         ek_pump_fac = 0.0_cp
      end if

      ratio1=alph1
      ratio2=alph2

      call rscheme%get_grid(n_r_max, r_icb, r_cmb, ratio1, ratio2, r)
      call rscheme%get_der_mat(n_r_max, l_cheb_coll)

      if ( rank == 0 ) then
         file_name = 'radius.'//tag
         open(newunit=file_handle, file=file_name, status='unknown')
         do n_r=1,n_r_max
            write(file_handle,'(I4, ES16.8)') n_r, r(n_r)
         end do
         close(file_handle)
      end if

      rgrav(:)=g0+g1*r(:)/r_cmb+g2*(r_cmb/r)**2

      or1(:)=one/r(:)      ! 1/r
      or2(:)=or1*or1(:)    ! 1/r**2

      !-- arrays for Courant conditions
      c1=(two*pi/(three*real(m_max,cp)))**2
      delxh2(1)      =c1*r_cmb**2
      delxh2(n_r_max)=c1*r_icb**2
      delxr2(1)      =(r(1)-r(2))**2
      delxr2(n_r_max)=(r(n_r_max-1)-r(n_r_max))**2
      do n_r=2,n_r_max-1
         delxh2(n_r)=c1*r(n_r)**2
         delmin=min((r(n_r-1)-r(n_r)),(r(n_r)-r(n_r+1)))
         delxr2(n_r)=delmin*delmin
      end do

      if ( l_temp_3D ) then
         tcond(1) = 0.0_cp
         do n_r=2,n_r_max
            h = sqrt(r_cmb**2-r(n_r)**2)
            tcond(n_r) = r_icb*(r_cmb/h*asinh(h/r(n_r))-one)
         end do
         call get_dr(tcond, dtcond, n_r_max, rscheme)
      else 
         tcond(:) = tcond_fac*(log(r(:))/log(radratio)-log(r_cmb)/log(radratio))
         dtcond(:)= tcond_fac*or1(:)/log(radratio)
      end if

      if ( l_non_rot ) then
         beta(:)   =0.0_cp
         dbeta(:)  =0.0_cp
         ekpump(:) =0.0_cp
         oheight(:)=1.0_cp
         height(:) =1.0_cp
      else ! Calculate beta only when this is rotating !
         height(1) = 0.0_cp
         beta(1)   = 0.0_cp
         dbeta(1)  = 0.0_cp
         ekpump(1) = ek_pump_fac*0.0e0_cp
         oheight(1)= 0.0e0_cp
         do n_r=2,n_r_max
            height(n_r) = two * sqrt(r_cmb**2-r(n_r)**2)
            oheight(n_r)= half/sqrt(r_cmb**2-r(n_r)**2)
            beta(n_r)   = -r(n_r)/(r_cmb**2-r(n_r)**2)
            dbeta(n_r)  = -(r_cmb**2+r(n_r)**2)/(r_cmb**2-r(n_r)**2)**2
            ekpump(n_r) = half*ek_pump_fac*sqrt(ek*r_cmb)/ &
            &             (r_cmb**2-r(n_r)**2)**(3.0_cp/4.0_cp)
         end do
      end if

   end subroutine radial
!------------------------------------------------------------------------------
end module radial_functions
