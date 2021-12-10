module radial_functions
   !
   ! This module calculates some arrays that depend on the radial direction
   ! only. This is called only once at the initialisation of the code.
   !

   use iso_fortran_env, only: output_unit
   use truncation, only: n_r_max, n_cheb_max, m_max
   use radial_der, only: get_dr
   use constants, only: one, two, three, pi, half, third, tiny_number
   use namelists, only: tag, alph1, alph2, l_newmap, radratio, g0, g1, g2, &
       &                l_non_rot, ek, l_ek_pump, l_temp_3D, tcond_fac,    &
       &                r_cmb, r_icb, l_cheb_coll, beta_shift, xicond_fac, &
       &                ktopt, kbott, t_bot, t_top, l_heat, l_chem, xi_bot,&
       &                xi_top, l_xi_3D, ktopxi, kbotxi, h_temp, h_xi,     &
       &                l_finite_diff, fd_stretch, fd_ratio, container,    &
       &                beta_fac
   use mem_alloc, only: bytes_allocated
   use radial_scheme, only: type_rscheme
   use chebyshev, only: type_cheb
   use finite_differences, only: type_fd
   use parallel_mod
   use precision_mod

   implicit none

   private
 
   !-- arrays depending on r:
   real(cp), public, allocatable :: r(:)         ! radii
   real(cp), public, allocatable :: or1(:)       ! :math:`1/r`
   real(cp), public, allocatable :: or2(:)       ! :math:`1/r^2`
   real(cp), public, allocatable :: or3(:)       ! :math:`1/r^3`
   real(cp), public, allocatable :: beta(:)      ! 1/h dh/ds
   real(cp), public, allocatable :: dbeta(:)     ! Radial gradient of beta
   real(cp), public, allocatable :: d2beta(:)    ! Radial gradient of beta
   real(cp), public, allocatable :: d3beta(:)    ! Radial gradient of beta
   real(cp), public, allocatable :: height(:)    ! Spherical shell height
   real(cp), public, allocatable :: ekpump(:)    ! Ekman pumping function
   real(cp), public, allocatable :: ekp_up(:)    ! Ekman pumping function (in front of up)
   real(cp), public, allocatable :: ekp_us(:)    ! Ekman pumping function (in front of us)
   real(cp), public, allocatable :: ekp_dusdp(:) ! Ekman pumping function (in front of dus/dphi)
   real(cp), public, allocatable :: oheight(:)   ! 1/h

   real(cp), public, allocatable :: delxr2(:) ! Auxiliary arrays containing effective Courant grid intervals
   real(cp), public, allocatable :: delxh2(:) ! Auxiliary arrays containing effective Courant grid intervals

   real(cp), public, allocatable :: tcond(:) ! Conducting temperature
   real(cp), public, allocatable :: dtcond(:)! radial derivative of conducting temperature
   real(cp), public, allocatable :: xicond(:) ! Conducting chemical composition
   real(cp), public, allocatable :: dxicond(:)! radial derivative


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
      allocate( r(n_r_max), or1(n_r_max), or2(n_r_max), or3(n_r_max) )
      allocate( beta(n_r_max), dbeta(n_r_max), height(n_r_max) )
      allocate( d2beta(n_r_max), d3beta(n_r_max) )
      allocate( rgrav(n_r_max), ekpump(n_r_max), oheight(n_r_max) )
      allocate( delxr2(n_r_max), delxh2(n_r_max), ekp_up(n_r_max) )
      allocate( tcond(n_r_max), dtcond(n_r_max), ekp_us(n_r_max) )
      allocate( xicond(n_r_max), dxicond(n_r_max), ekp_dusdp(n_r_max) )
      bytes_allocated = bytes_allocated+21*n_r_max*SIZEOF_DEF_REAL

      if ( .not. l_finite_diff ) then
         allocate ( type_cheb :: rscheme )
         n_in = n_cheb_max
         if ( l_newmap ) then
            n_in_2 = 1
         else
            n_in_2 = 0
         end if
         call rscheme%initialize(n_r_max,n_in,n_in_2,l_cheb_coll)
      else
         allocate ( type_fd :: rscheme )
         call rscheme%initialize(n_r_max,2,2,.true.)
      end if


   end subroutine initialize_radial_functions
!------------------------------------------------------------------------------
   subroutine finalize_radial_functions

      call rscheme%finalize()

      deallocate( tcond, dtcond, xicond, dxicond )
      deallocate( delxr2, delxh2, rgrav )
      deallocate( beta, dbeta, height, ekpump, oheight, d2beta, d3beta )
      deallocate( r, or1, or2, or3, ekp_up, ekp_us, ekp_dusdp )

   end subroutine finalize_radial_functions
!------------------------------------------------------------------------------
   subroutine radial
      !
      !  Calculates everything needed for radial functions, transforms etc.
      !

      integer :: n_r, file_handle
      character(len=100) :: file_name
      real(cp) :: ratio1, ratio2, delmin, c1
      real(cp) :: ek_pump_fac, fac, Lin, Lout

      if ( l_ek_pump .and. (.not. l_non_rot)) then
         ek_pump_fac = one
      else
         ek_pump_fac = 0.0_cp
      end if

      if ( l_finite_diff ) then
         ratio1=fd_stretch
         ratio2=fd_ratio
      else
         ratio1=alph1
         ratio2=alph2
      end if

      call rscheme%get_grid(n_r_max, r_icb, r_cmb, ratio1, ratio2, r)
      if ( .not. l_finite_diff ) then
         call rscheme%get_der_mat(n_r_max, l_cheb_coll)
      end if

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
      or3(:)=or2*or1(:)    ! 1/r**3

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

      !-- Calculate conducting temperature
      if ( l_heat ) call get_conducting_state(tcond,dtcond,l_temp_3D,tcond_fac,&
                         &                    ktopt,kbott,t_top,t_bot, h_temp)

      !-- Calculate chemical composition
      if ( l_chem ) call get_conducting_state(xicond,dxicond,l_xi_3D,xicond_fac,&
                         &                    ktopxi,kbotxi,xi_top,xi_bot, h_xi)

      if ( l_non_rot ) then
         beta(:)     =0.0_cp
         dbeta(:)    =0.0_cp
         ekpump(:)   =0.0_cp
         ekp_up(:)   =0.0_cp
         ekp_us(:)   =0.0_cp
         ekp_dusdp(:)=0.0_cp
         oheight(:)  =1.0_cp
         height(:)   =1.0_cp
         d2beta(:)   =0.0_cp
         d3beta(:)   =0.0_cp
      else ! Calculate beta only when this is rotating !
         if ( index(container, 'SPHERE') == 1 ) then
            if ( abs(beta_shift) <= 10.0_cp*epsilon(beta_shift) ) then
               height(1) = 0.0_cp
               beta(1)   = 0.0_cp
               dbeta(1)  = 0.0_cp
               d2beta(1) = 0.0_cp
               d3beta(1) = 0.0_cp
               ekpump(1) = ek_pump_fac*0.0e0_cp
               oheight(1)= 0.0e0_cp
               do n_r=2,n_r_max
                  height(n_r) = two * sqrt(r_cmb**2-r(n_r)**2)
                  oheight(n_r)= half/sqrt(r_cmb**2-r(n_r)**2)
                  beta(n_r)   = -r(n_r)/(r_cmb**2-r(n_r)**2)
                  dbeta(n_r)  = -(r_cmb**2+r(n_r)**2)/(r_cmb**2-r(n_r)**2)**2
                  d2beta(n_r) = -128.0_cp*r(n_r)*(three*r_cmb**2+r(n_r)**2)/height(n_r)**6
                  d3beta(n_r) = -1536.0_cp*(r_cmb**4+r(n_r)**4+6.0_cp*r(n_r)**2*r_cmb**2)/&
                  &             height(n_r)**8
                  ekpump(n_r) = half*ek_pump_fac*sqrt(ek*r_cmb)/ &
                  &             (r_cmb**2-r(n_r)**2)**0.75_cp
               end do
            else
               do n_r=1,n_r_max
                  height(n_r) = two * sqrt((r_cmb+beta_shift)**2-r(n_r)**2)
                  oheight(n_r)= half/sqrt((r_cmb+beta_shift)**2-r(n_r)**2)
                  beta(n_r)   = -r(n_r)/((r_cmb+beta_shift)**2-r(n_r)**2)
                  dbeta(n_r)  = -((r_cmb+beta_shift)**2+r(n_r)**2)/   &
                  &              ((r_cmb+beta_shift)**2-r(n_r)**2)**2
                  d2beta(n_r) = -two*r(n_r)*(three*(r_cmb+beta_shift)**2+r(n_r)**2)/ &
                  &              height(n_r)**6
                  d3beta(n_r) = -6.0_cp*((r_cmb+beta_shift)**4+r(n_r)**4+6.0_cp* &
                  &              r(n_r)**2*(r_cmb+beta_shift)**2)/height(n_r)**8
                  ekpump(n_r) = half*ek_pump_fac*sqrt(ek*r_cmb)/      &
                  &             ((r_cmb+beta_shift)**2-r(n_r)**2)**0.75_cp
               end do
            end if
            ekp_up(:)   =-half*beta(:)
            ekp_us(:)   =-5.0_cp*r_cmb*beta(:)*oheight(:)
            ekp_dusdp(:)=beta(:)
         else if ( index(container, 'EXP') == 1 ) then
            fac=beta_fac/r_cmb ! container with L=exp(-fac * r)
            height(:) =exp(-fac*r(:))
            beta(:)   =-fac
            dbeta(:)  =0.0_cp
            d2beta(:) =0.0_cp
            d3beta(:) =0.0_cp
            ! n.ez = cos(theta); tan(theta)=dh/ds; n.ez=1/sqrt((dh/ds)^2+1)
            ! Here: n.ez=1/sqrt(fac^2*exp(-2*gamma*r)+1)
            ! Ekpump = E^{-1/2} / h / sqrt(n.ez)
            ekpump(:)   =half *ek_pump_fac * sqrt(ek) / height(:) * &
            &            (one+fac**2*height(:)**2)**0.25_cp
            ekp_up(:)   =-half*fac**3*height(:)**2/(one+fac**2*height(:)**2)
            ekp_us(:)   =half*(two*fac-fac**3*height(:)**2)/sqrt(one+fac**2*height(:)**2)
            ekp_dusdp(:)=-fac**2*height**2*or1(:)
            height(:)   =two*exp(-fac*r(:)) ! twice for total height afterwards
            oheight(:)  =one/height(:)
         else if ( index(container, 'BUSSE') == 1 ) then
            Lin =one ! Normalised at one in full disks
            Lout=beta_fac
            fac =(Lin-Lout)/(r_cmb-r_icb)/r_cmb
            height(:)   =-fac*r(:)+Lin
            beta(:)     =-fac/height
            dbeta(:)    =-fac**2/height(:)**2
            d2beta(:)   =two*fac**3/height**3
            d3beta(:)   =-6.0_cp*fac**4/height(:)**4
            ekpump(:)   =half*ek_pump_fac*sqrt(ek) / height(:) *(one+fac**2)**0.25_cp
            ekp_up(:)   =0.0_cp
            ekp_us(:)   =fac*(1+fac**2)**0.75_cp/height(:)
            ekp_dusdp(:)=-fac**2*or1(:)
            height(:) =two*(fac*r(:)+Lin) ! twice for total height afterwards
            oheight(:)=one/height(:)
         end if
      end if

   end subroutine radial
!------------------------------------------------------------------------------
   subroutine get_conducting_state(tcond, dtcond, l_3D, cond_fac, ktop, kbot, &
              &                    t_top, t_bot, epsc0)
      !
      !  Calculates the conducting state of the temperature equation
      !

      !-- Input variables
      logical,  intent(in) :: l_3D
      integer,  intent(in) :: ktop
      integer,  intent(in) :: kbot
      real(cp), intent(in) :: cond_fac
      real(cp), intent(in) :: t_bot(:)
      real(cp), intent(in) :: t_top(:)
      real(cp), intent(inout) :: epsc0

      !-- Output variables
      real(cp), intent(inout) ::  tcond(n_r_max)
      real(cp), intent(inout) :: dtcond(n_r_max)

      !-- Local variables
      integer :: n, n_r, m_bot, m_top
      real(cp) :: h(n_r_max)     ! :math:`sqrt(r_cmb^2-r^2)`
      real(cp) :: tr_bot, tr_top
      real(cp) :: f_bot, f_top

      f_bot=0.0_cp
      f_top=0.0_cp
      do n=1,size(t_bot)/3
         m_bot =int(t_bot(3*n-2))
         tr_bot=t_bot(3*n-1)
         m_top =int(t_top(3*n-2))
         tr_top=t_top(3*n-1)
         if ( m_bot == 0 .and. abs(tr_bot) > 10.0_cp*epsilon(one) ) then
            f_bot=real(tr_bot,kind=cp)
         end if
         if ( m_top == 0 .and. abs(tr_top) > 10.0_cp*epsilon(one) ) then
            f_top=real(tr_top,kind=cp)
         end if
      end do

      if ( rank == 0 ) then
         if ( ktopt == 1 ) then
            write(output_unit,*) '! Top value', f_top
         else
            write(output_unit,*) '! Top flux', f_top
         end if
         if ( kbott == 1 ) then
            write(output_unit,*) '! Bot value', f_bot
         else
            write(output_unit,*) '! Bot flux', f_bot
         end if
      end if

      !-- Conductive Temperature profile 2D and 3D projected onto QG
      if ( l_3D ) then
         !3D-tcond profiles -- Warning! In 3D, tcond(r_cmb) induces division by 0 (in h)
         tcond(1) = f_top!0.0_cp
         h = half*height(:)!sqrt(r_cmb**2-r(:)**2)
         !-- 3D heat sources to compensate for top/bottom fluxes
         epsc0 = three*(r_icb**2*f_bot - r_cmb**2*f_top)/ &
         &             (r_cmb**3 - r_icb**3)
         if ( epsc0<10.0_cp*epsilon(one) ) epsc0=0.0_cp
         if ( ktop==1 .and. kbot==1 ) then
            do n_r=2,n_r_max
               tcond(n_r) = (r_icb/(r_icb-r_cmb))*(one-r_cmb*asinh(h(n_r)/r(n_r))/h(n_r))
            end do
         elseif ( ktop==2 .and. kbot==1 ) then
            do n_r=2,n_r_max
               tcond(n_r) = one + f_top*r_cmb**2*(one/r_icb - asinh(h(n_r)/r(n_r))/h(n_r)) &
               &          + (epsc0/two)*(r_icb**2 - r(n_r)**2 - third*h(n_r)**2 +          &
               &             two*r_cmb**3*(one/r_icb-asinh(h(n_r)/r(n_r))/h(n_r)))
            end do
         elseif ( ktop==1 .and. kbot==2 ) then
            do n_r=2,n_r_max
               tcond(n_r) = f_bot*r_icb**2*(one/r_cmb - asinh(h(n_r)/r(n_r))/h(n_r)) &
               &          + (epsc0/two)*(r_cmb**2 - r(n_r)**2 - third*h(n_r)**2 +    &
               &             two*r_icb**3*(one/r_cmb-asinh(h(n_r)/r(n_r))/h(n_r)))
            enddo
         else
            do n_r=2,n_r_max
               tcond(n_r) = f_top*r_cmb - f_top*r_cmb**2*(asinh(h(n_r)/r(n_r))/h(n_r)) &
               &          - (epsc0/two)*(r(n_r)**2 + third*h(n_r)**2 +                 &
               &             two*r_cmb**3*(asinh(h(n_r)/r(n_r))/h(n_r)))
            enddo 
         endif
         call get_dr(tcond, dtcond, n_r_max, rscheme)
      else
         !2D-tcond profiles -- Warning! In 2D, need a geometric factor
         !-- 2D heat sources to compensate for top/bottom fluxes
         if ( ktop==1 .and. kbot==1 ) then
            tcond(:) = cond_fac*(log(r(:)/r_cmb)/log(radratio))
            dtcond(:)= cond_fac*(or1(:)/log(radratio))
         elseif ( ktop==2 .and. kbot==1 ) then
            tcond(:) =cond_fac*(f_bot + f_top*r_cmb*log(r(:)/r_icb) + 0.25_cp*epsc0*  &
            &                   (r_icb**2 - r(:)**2 + two*r_cmb**2*log(r(:)/ r_icb)))
            dtcond(:)=cond_fac*(f_top*r_cmb*or1(:)+half*epsc0*(r_cmb**2*or1(:) - r(:)))
         elseif( ktop==1 .and. kbot==2 ) then
            tcond(:) =cond_fac*(f_bot*r_icb*log(r(:)/r_cmb) + 0.25_cp*epsc0* &
            &                    (r_cmb**2 - r(:)**2 + two*r_icb**2*log(r(:)/ r_cmb))+&
            &                    f_top)
            dtcond(:)=cond_fac*(f_bot*r_icb*or1(:)+half*epsc0*(r_icb**2*or1(:) - r(:)))
         else
            epsc0 = -two*(r_cmb*f_top-r_icb*f_bot)/(r_cmb**2-r_icb**2)
            if ( abs(epsc0)<tiny_number ) epsc0=0.0_cp
            tcond(:) = cond_fac*( half*epsc0*( r_cmb**2*log(r(:)/r_cmb)-half*(r(:)**2- &
            &                     r_cmb**2))+r_cmb*f_top*log(r(:)/r_cmb))
            dtcond(:)=cond_fac*(f_top*r_cmb*or1(:)+half*epsc0*(r_cmb**2*or1(:)-r(:)))
         endif
         if (rank == 0) write(output_unit,*) &
         &         '! Warning: Sources introduced to balance surface heat flux'
         if (rank == 0) write(output_unit,'(''!      epsc0*pr='',ES16.6)') epsc0
      end if

#ifdef WITH_DEBUG
      block
         integer :: file_handle
         if ( rank == 0 ) then
            open(newunit=file_handle, file='cond')
            do n_r=1,n_r_max
               write(file_handle, '(3ES16.8)') r(n_r), tcond(n_r), dtcond(n_r)
            end do
            close(file_handle)
         end if
      end block
#endif

   end subroutine get_conducting_state
!------------------------------------------------------------------------------
end module radial_functions
