module radial_functions
   !
   ! This module calculates some arrays that depend on the radial direction
   ! only. This is called only once at the initialisation of the code.
   !

   use iso_fortran_env, only: output_unit
   use blocking, only: nMstart, nMstop, lmStart, lmStop
   use truncation, only: n_r_max, n_cheb_max, m_max
   use truncation_3D, only: n_r_max_3D, n_cheb_max_3D, l_max
   use radial_der, only: get_dr
   use constants, only: one, two, three, pi, half, third, sq4pi
   use namelists, only: tag, alph1, alph2, l_newmap, radratio, g0, g1, g2,   &
       &                l_non_rot, ek, l_ek_pump, pr, l_tcond_3D, tcond_fac, &
       &                r_cmb, r_icb, l_cheb_coll, beta_shift, xicond_fac,   &
       &                ktopt, kbott, t_bot, t_top, l_heat, l_chem, xi_bot,  &
       &                xi_top, l_xi_3D, ktopxi, kbotxi, l_3D, l_heat_3D,    &
       &                epsrc0, epsrc0xi
   use mem_alloc, only: bytes_allocated
   use radial_scheme, only: type_rscheme
   use chebyshev, only: type_cheb
   use useful, only: abortRun
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
   real(cp), public, allocatable :: dtcond(:)! radial derivative of conducting temperature
   real(cp), public, allocatable :: xicond(:) ! Conducting chemical composition
   real(cp), public, allocatable :: dxicond(:)! radial derivative


   real(cp), public :: alpha1   ! Input parameter for non-linear map to define degree of spacing (0.0:2.0)
   real(cp), public :: alpha2   ! Input parameter for non-linear map to define central point of different spacing (-1.0:1.0)
 
   !-- arrays for buoyancy, depend on Ra and Pr:
   real(cp), public, allocatable :: rgrav(:)     ! Buoyancy term `dtemp0/Di`

   !-- arrays depending on spherical r:
   real(cp), public, allocatable :: r_3D(:)         ! spherical radii
   real(cp), public, allocatable :: or1_3D(:)       ! :math:`1/r_3D`
   real(cp), public, allocatable :: or2_3D(:)       ! :math:`1/r_3D^2`
   real(cp), public, allocatable :: rgrav_3D(:)
   real(cp), public, allocatable :: tcond_3D(:)
   real(cp), public, allocatable :: dtcond_3D(:)

   real(cp), public, allocatable :: delxr2_3D(:) ! Auxiliary arrays containing effective Courant 3D grid intervals
   real(cp), public, allocatable :: delxh2_3D(:) ! Auxiliary arrays containing effective Courant 3D grid intervals

   !-- Radial scheme
   class(type_rscheme), public, pointer :: rscheme, rscheme_3D
 
   public :: initialize_radial_functions, radial, finalize_radial_functions, &
   &         radial_3D

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
      allocate( xicond(n_r_max), dxicond(n_r_max))
      tcond(:)=0.0_cp ; dtcond(:)=0.0_cp; xicond(:)=0.0_cp; dxicond(:)=0.0_cp
      bytes_allocated = bytes_allocated+15*n_r_max*SIZEOF_DEF_REAL

      allocate ( type_cheb :: rscheme )

      n_in = n_cheb_max
      if ( l_newmap ) then
         n_in_2 = 1
      else
         n_in_2 = 0
      end if

      call rscheme%initialize(nMstart,nMstop,n_r_max,n_in,n_in_2,l_cheb_coll)

      if ( l_3D ) then
         allocate( r_3D(n_r_max_3D), or1_3D(n_r_max_3D), or2_3D(n_r_max_3D) )
         allocate( rgrav_3D(n_r_max_3D) )
         allocate( delxr2_3D(n_r_max_3D), delxh2_3D(n_r_max_3D) )
         bytes_allocated = bytes_allocated+6*n_r_max_3D*SIZEOF_DEF_REAL
         allocate( tcond_3D(n_r_max_3D),dtcond_3D(n_r_max_3D) )
         bytes_allocated = bytes_allocated+2*n_r_max_3D*SIZEOF_DEF_REAL

         allocate ( type_cheb :: rscheme_3D )
         call rscheme_3D%initialize(lmStart,lmStop,n_r_max_3D,n_cheb_max_3D,0,.true.)
      end if

   end subroutine initialize_radial_functions
!------------------------------------------------------------------------------
   subroutine finalize_radial_functions

      call rscheme%finalize()

      deallocate( tcond, dtcond, xicond, dxicond )
      deallocate( delxr2, delxh2, rgrav )
      deallocate( beta, dbeta, height, ekpump, oheight )
      deallocate( r, or1, or2 )

      if ( l_3D ) then
         deallocate( tcond_3D, dtcond_3D )
         deallocate( delxr2_3D, delxh2_3D )
         deallocate( r_3D, or1_3D, or2_3D, rgrav_3D )
         call rscheme_3D%finalize()
      end if

   end subroutine finalize_radial_functions
!------------------------------------------------------------------------------
   subroutine radial
      !
      !  Calculates everything needed for radial functions, transforms etc.
      !

      integer :: n_r, file_handle
      character(len=100) :: file_name
      real(cp) :: ratio1, ratio2, delmin, c1
      real(cp) :: ek_pump_fac

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

      if ( l_non_rot ) then
         beta(:)   =0.0_cp
         dbeta(:)  =0.0_cp
         ekpump(:) =0.0_cp
         oheight(:)=1.0_cp
         height(:) =1.0_cp
      else ! Calculate beta only when this is rotating !
         if ( abs(beta_shift) <= 10.0_cp*epsilon(beta_shift) ) then
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
         else
            do n_r=1,n_r_max
               height(n_r) = two * sqrt((r_cmb+beta_shift)**2-r(n_r)**2)
               oheight(n_r)= half/sqrt((r_cmb+beta_shift)**2-r(n_r)**2)
               beta(n_r)   = -r(n_r)/((r_cmb+beta_shift)**2-r(n_r)**2)
               dbeta(n_r)  = -((r_cmb+beta_shift)**2+r(n_r)**2)/   &
               &              ((r_cmb+beta_shift)**2-r(n_r)**2)**2
               ekpump(n_r) = half*ek_pump_fac*sqrt(ek*r_cmb)/      &
               &             ((r_cmb+beta_shift)**2-r(n_r)**2)**(3.0_cp/4.0_cp)
            end do
         end if
      end if

      !-- Calculate conducting temperature
      if ( l_heat ) call get_conducting_state(tcond,dtcond,l_tcond_3D,tcond_fac,&
                         &                    ktopt,kbott,t_top,t_bot,epsrc0)

      !-- Calculate chemical composition
      if ( l_chem ) call get_conducting_state(xicond,dxicond,l_xi_3D,xicond_fac,&
                         &                    ktopxi,kbotxi,xi_top,xi_bot,epsrc0xi)

   end subroutine radial
!------------------------------------------------------------------------------
   subroutine radial_3D
      !
      !  Calculates everything needed for 3D radial functions, transforms etc.
      !

      integer :: n, n_r, file_handle
      character(len=100) :: file_name
      real(cp) :: ratio1, ratio2, delmin
      real(cp) :: tr_bot, tr_top, m_bot, l_bot, m_top, l_top
      real(cp) :: bott_bc, topt_bc
      real(cp) :: tmpsc0

      ratio1=alph1 ! To be changed to alph1_3D at some point
      ratio2=alph2

      call rscheme_3D%get_grid(n_r_max_3D, r_icb, r_cmb, ratio1, ratio2, r_3D)
      call rscheme_3D%get_der_mat(n_r_max_3D, .true.)

      if ( rank == 0 ) then
         file_name = 'radius_3D.'//tag
         open(newunit=file_handle, file=file_name, status='unknown')
         do n_r=1,n_r_max_3D
            write(file_handle,'(I4, ES16.8)') n_r, r_3D(n_r)
         end do
         close(file_handle)
      end if


      rgrav_3D(:)=g0+g1*r_3D(:)/r_cmb+g2*(r_cmb/r_3D)**2

      or1_3D(:)=one/r_3D(:)      ! 1/r_3D
      or2_3D(:)=or1_3D*or1_3D(:) ! 1/r_3D**2

      !-- arrays for Courant and Alfven conditions
      delxh2_3D(1)         =r_cmb**2/real(l_max*(l_max+1),kind=cp)
      delxh2_3D(n_r_max_3D)=r_icb**2/real(l_max*(l_max+1),kind=cp)
      delxr2_3D(1)         =(r_3D(1)-r_3D(2))**2
      delxr2_3D(n_r_max_3D)=(r_3D(n_r_max_3D-1)-r_3D(n_r_max_3D))**2
      do n_r=2,n_r_max_3D-1
         delxh2_3D(n_r)=r_3D(n_r)**2/real(l_max*(l_max+1),kind=cp)
         delmin=min((r_3D(n_r-1)-r_3D(n_r)),(r_3D(n_r)-r_3D(n_r+1)))
         delxr2_3D(n_r)=delmin*delmin
      end do

      !-- Get Conducting state 3D
      if ( l_heat_3D ) then
         bott_bc=0.0_cp
         topt_bc=0.0_cp
         do n=1,size(t_bot)/4
            l_bot =int(t_bot(4*n-3))
            m_bot =int(t_bot(4*n-2))
            tr_bot=t_bot(4*n-1)
            l_top =int(t_top(4*n-3))
            m_top =int(t_top(4*n-2))
            tr_top=t_top(4*n-1)
            if ( l_bot == 0 .and. m_bot == 0 .and. abs(tr_bot) > 10.0_cp*epsilon(one) ) then
               bott_bc=real(tr_bot,kind=cp)
            end if
            if ( l_top == 0 .and. m_top == 0 .and. abs(tr_top) > 10.0_cp*epsilon(one) ) then
               topt_bc=real(tr_top,kind=cp)
            end if
         end do

         !-- Actual 3D Conducting state
         if ( kbott == 1 .and. ktopt == 1 ) then
            !-- Tcond = DeltaT(=1)*re*ri(1/re - 1/r)/d(=1) + Te
            if ( (bott_bc-topt_bc) /= one ) then
               !-- Because T_3D is written in Full Temperature:
               !-- When fixed temperature contrast, DeltaT has to be 1
               !-- If Te specified and =/= 0, then Ti should be modified accordingly!
               if ( rank == 0 ) write(output_unit,*) &
               &                '! Warning: Top/Bot Temperature Boundary conditions corrected:'
               if ( rank == 0 ) write(output_unit,*) &
               &                '! Delta T should be 1 in fixed temperature contrast!'
               if ( bott_bc /= 0.0_cp ) then
                  topt_bc = bott_bc - one
               else if ( topt_bc /= 0.0_cp ) then
                  bott_bc = topt_bc + one
               else
                  !-- Default is given by T_cmb (most likely everything at 0)
                  bott_bc = topt_bc + one
               end if
            end if
            tcond_3D(:) = r_icb*r_cmb*or1_3D(:)-r_icb + topt_bc
            dtcond_3D(:)=-r_icb*r_cmb*or2_3D(:)

         else if ( ktopt==2 .and. kbott==1 ) then
            tcond_3D(:) = bott_bc + topt_bc*r_cmb**2*(one/r_icb - or1_3D(:))
            if ( epsrc0 /= 0.0_cp ) then
               tcond_3D(:) = tcond_3D(:) + epsrc0/6.0_cp* &
               &            (r_icb**2 - r_3D(n_r)**2 + two*r_cmb**3*(one/r_icb-or1_3D(:)))
            end if
            call get_dr(tcond_3D, dtcond_3D, n_r_max_3D, rscheme_3D)

         else if ( ktopt==1 .and. kbott==2 ) then
            tcond_3D(:) = topt_bc + bott_bc*r_icb**2*(one/r_cmb - or1_3D(:))
            if ( epsrc0 /= 0.0_cp ) then
               tcond_3D(:) = tcond_3D(:) + epsrc0/6.0_cp* &
               &            (r_cmb**2 - r_3D(n_r)**2 + two*r_icb**3*(one/r_cmb-or1_3D(:)))
            end if
            call get_dr(tcond_3D, dtcond_3D, n_r_max_3D, rscheme_3D)

         else if ( ktopt==2 .and. kbott==2 ) then !testing this
            !-- When top flux and bot flux conditions are used,
            !-- 3D heat sources have to compensate for top/bottom fluxes
            tmpsc0 = three*(r_icb**2*bott_bc - r_cmb**2*topt_bc)/ &
            &              (r_cmb**3 - r_icb**3)/pr
            if ( tmpsc0 /= epsrc0 ) then
               if ( rank == 0 ) write(output_unit,*) '! Warning: Sources introduced/corrected to balance heat fluxes'
               if ( rank == 0 ) write(output_unit,'(''!      epsrc0 ='',ES16.6)') tmpsc0
               epsrc0 = tmpsc0
            end if
            tcond_3D(:) = - topt_bc*r_cmb**2*or1_3D(:) - epsrc0/6.0_cp* &
            &               (r_3D(n_r)**2 + two*r_cmb**3*or1_3D(:))
            call get_dr(tcond_3D, dtcond_3D, n_r_max_3D, rscheme_3D)

         else
            if ( rank == 0 ) write(output_unit,*) 'Only kbott x ktopt = {1,2}x{1,2} are supported !'
            call abortRun('3D Tcond with BCs =/= {1,2} not implemented')
         end if
      else
         tcond_3D(:) = 0.0_cp
         dtcond_3D(:) = 0.0_cp
      end if

      if ( rank == 0 ) write(output_unit,*) '! 3D Temperature Boundary Condition at the CMB', topt_bc
      if ( rank == 0 ) write(output_unit,*) '! 3D Temperature Boundary Condition at the ICB', bott_bc

      tmpsc0 = epsrc0*pr*4.0_cp*third*pi*(r_cmb**3-r_icb**3)
      if ( rank == 0 ) write(output_unit,*) '! Total volumetric buoyancy Sources = ', tmpsc0

   end subroutine radial_3D
!------------------------------------------------------------------------------
   subroutine get_conducting_state(tcond, dtcond, l_tcond_3D, cond_fac, ktop, &
              &                    kbot, t_top, t_bot, epsc0)
      !
      !  Calculates the conducting state of the temperature equation
      !

      !-- Input variables
      logical,  intent(in) :: l_tcond_3D
      integer,  intent(in) :: ktop
      integer,  intent(in) :: kbot
      real(cp), intent(in) :: cond_fac
      real(cp), intent(in) :: t_bot(:)
      real(cp), intent(in) :: t_top(:)

      !-- Output variables
      real(cp), intent(inout) ::  tcond(n_r_max)
      real(cp), intent(inout) :: dtcond(n_r_max)
      real(cp), intent(inout) :: epsc0

      !-- Local variables
      integer :: n, n_r, m_bot, m_top!, n_out
      real(cp) :: h(n_r_max)     ! :math:`sqrt(r_cmb^2-r^2)`
      real(cp) :: tr_bot, tr_top
      real(cp) :: bott_bc, topt_bc
      real(cp) :: tmpsc0

      bott_bc=0.0_cp
      topt_bc=0.0_cp
      do n=1,size(t_bot)/3
         m_bot =int(t_bot(3*n-2))
         tr_bot=t_bot(3*n-1)
         m_top =int(t_top(3*n-2))
         tr_top=t_top(3*n-1)
         if ( m_bot == 0 .and. abs(tr_bot) > 10.0_cp*epsilon(one) ) then
            bott_bc=real(tr_bot,kind=cp)
         end if
         if ( m_top == 0 .and. abs(tr_top) > 10.0_cp*epsilon(one) ) then
            topt_bc=real(tr_top,kind=cp)
         end if

         if ( ktop==1 .and. kbot==1 ) then
            if ( (bott_bc-topt_bc) /= one ) then
               !-- When fixed temperature contrast, DeltaT has to be 1
               !-- If Te specified and =/= 0, then Ti should be modified accordingly!
               if ( rank == 0 ) write(output_unit,*) &
               &                '! Warning: Top/Bot Temperature Boundary conditions corrected:' 
               if ( rank == 0 ) write(output_unit,*) &
               &                '! Delta T should be 1 in fixed temperature contrast!'
               if ( bott_bc /= 0.0_cp ) then
                  topt_bc = bott_bc - one
               else if ( topt_bc /= 0.0_cp ) then
                  bott_bc = topt_bc + one
               else
                  !-- Default is given by T_cmb (most likely everything at 0)
                  bott_bc = topt_bc + one
               end if
            end if
         end if
      end do

      if ( abs(epsc0)<10.0_cp*epsilon(one) ) epsc0=0.0_cp

      !-- Conductive Temperature profile 2D and 3D projected onto QG
      if ( l_tcond_3D ) then

         !-- Warning! if l_tcond_3D, tcond(r_cmb) induces division by 0 (in h)
         h = half*height(:)!sqrt(r_cmb**2-r(:)**2)

         if ( ktop==1 .and. kbot==1 ) then
            !-- Tcond = <[DeltaT(=1) * re*ri/(ri-re)*(1/re - 1/r) + Te]>_z
            do n_r=2,n_r_max
               tcond(n_r) = (r_icb/(r_icb-r_cmb))*(one-r_cmb*asinh(h(n_r)/r(n_r))/h(n_r)) + topt_bc
            end do
            tcond(1) = topt_bc !--> Tcond(r_cmb) induces division by 0 (in h)

         else if ( ktop==2 .and. kbot==1 ) then
            do n_r=2,n_r_max
               tcond(n_r) = bott_bc + topt_bc*r_cmb**2*(one/r_icb - asinh(h(n_r)/r(n_r))/h(n_r))
            end do
            tcond(1) = bott_bc + topt_bc*r_cmb**2*(one/r_icb - or1(1)) !--> Approximated by tcond_3D(1)
            if ( epsc0 /= 0.0_cp ) then
               do n_r=2,n_r_max
                  tcond(n_r) = tcond(n_r) + epsc0/6.0_cp*(r_icb**2 - r(n_r)**2 - third*h(n_r)**2 + &
                  &                         two*r_cmb**3*(one/r_icb-asinh(h(n_r)/r(n_r))/h(n_r)))
               end do
               tcond(1) = tcond(1) + epsc0/6.0_cp*(r_icb**2 - r(1)**2 + &
               &                     two*r_cmb**3*(one/r_icb-or1(1)))
            end if

         else if ( ktop==1 .and. kbot==2 ) then
            do n_r=2,n_r_max
               tcond(n_r) = topt_bc + bott_bc*r_icb**2*(one/r_cmb - asinh(h(n_r)/r(n_r))/h(n_r))
            enddo
            tcond(1) = topt_bc + bott_bc*r_icb**2*(one/r_cmb - or1(1)) !--> Approximated by tcond_3D(1)
            if ( epsc0 /= 0.0_cp ) then
               do n_r=2,n_r_max
                  tcond(n_r) = tcond(n_r) + epsc0/6.0_cp*(r_cmb**2 - r(n_r)**2 - third*h(n_r)**2 + &
                  &                         two*r_icb**3*(one/r_cmb-asinh(h(n_r)/r(n_r))/h(n_r)))
               end do
               tcond(1) = tcond(1) + epsc0/6.0_cp*(r_cmb**2 - r(1)**2 + &
               &                     two*r_icb**3*(one/r_cmb-or1(1)))
            end if

         else if ( ktop==2 .and. kbot==2 ) then
            !-- When top flux and bot flux conditions are used,
            !-- 3D heat sources have to compensate for top/bottom fluxes
            tmpsc0 = three*(r_icb**2*bott_bc - r_cmb**2*topt_bc)/ &
            &              (r_cmb**3 - r_icb**3)/pr
            if ( tmpsc0 /= epsc0 ) then
               if ( rank == 0 ) write(output_unit,*) '! Warning: Sources introduced/corrected to balance heat fluxes'
               if ( rank == 0 ) write(output_unit,'(''!      epsrc0 ='',ES16.6)') tmpsc0
               epsc0 = tmpsc0
            end if
            do n_r=2,n_r_max
               tcond(n_r) = - topt_bc*r_cmb**2*(asinh(h(n_r)/r(n_r))/h(n_r)) &
               &            - epsc0/6.0_cp*(r(n_r)**2 + third*h(n_r)**2 +    &
               &              two*r_cmb**3*(asinh(h(n_r)/r(n_r))/h(n_r)))
            enddo
            tcond(1) = - topt_bc*r_cmb**2*or1(1)     &
            &          - epsc0/6.0_cp*(r(1)**2 + two*r_cmb**3*or1(1)) !--> Approximated by tcond_3D(1)

         else
            if ( rank == 0 ) write(output_unit,*) 'Only kbott x ktopt = {1,2}x{1,2} are supported !'
            call abortRun('Tcond with BCs =/= {1,2} not implemented')
         endif
         call get_dr(tcond, dtcond, n_r_max, rscheme)

      else

         !2D-tcond profiles -- QG-2D profiles may recquire a geometric factor
         if ( ktop==1 .and. kbot==1 ) then
            !-- Tcond = DeltaT(=1) * one/(ln(ri)-ln(re))*(ln(r) - ln(re)) + Te
            tcond(:) = cond_fac*(log(r(:)/r_cmb)/log(radratio) + topt_bc)
            dtcond(:)= cond_fac*(or1(:)/log(radratio))

         else if ( ktop==2 .and. kbot==1 ) then
            tcond(:) = cond_fac*(bott_bc + topt_bc*r_cmb*log(r(:)/r_icb))
            dtcond(:)= cond_fac*(topt_bc*r_cmb*or1(:))
            if ( epsc0 /= 0.0_cp ) then
               tcond(:) = tcond(:) + cond_fac*(0.25_cp*epsc0*(r_icb**2 - r(:)**2 + &
               &                               two*r_cmb**2*log(r(:)/r_icb)))
               dtcond(:)= dtcond(:)+ cond_fac*(half*epsc0*(r_cmb**2*or1(:) - r(:)))
            end if

         else if( ktop==1 .and. kbot==2 ) then
            tcond(:) = cond_fac*(topt_bc + bott_bc*r_icb*log(r(:)/r_cmb))
            dtcond(:)= cond_fac*(bott_bc*r_icb*or1(:))
            if ( epsc0 /= 0.0_cp ) then
               tcond(:) = tcond(:) + cond_fac*(0.25_cp*epsc0*(r_cmb**2 - r(:)**2 + &
               &                               two*r_icb**2*log(r(:)/r_cmb)))
               dtcond(:)= dtcond(:)+ cond_fac*(half*epsc0*(r_icb**2*or1(:) - r(:)))
            end if

         else if( ktop==2 .and. kbot==2 ) then
            !-- When top flux and bot flux conditions are used,
            !-- 2D heat sources have to compensate for top/bottom fluxes
            tmpsc0 = two*(r_icb*bott_bc - r_cmb*topt_bc)/ &
            &            (r_cmb**2 - r_icb**2)/pr
            if ( tmpsc0 /= epsc0 ) then
               if ( rank == 0 ) write(output_unit,*) '! Warning: Sources introduced/corrected to balance heat fluxes'
               if ( rank == 0 ) write(output_unit,'(''!      epsrc0 ='',ES16.6)') tmpsc0
               epsc0 = tmpsc0
            end if
            tcond(:) = cond_fac*(topt_bc*r_cmb*log(r(:)) + 0.25_cp*epsc0* &
            &                   (two*r_cmb**2*log(r(:)) - r(:)**2))
            dtcond(:)= cond_fac*(topt_bc*r_cmb*or1(:) + half*epsc0*(r_cmb**2*or1(:) - r(:)))

         else
            if ( rank == 0 ) write(output_unit,*) 'Only kbott x ktopt = {1,2}x{1,2} are supported !'
            call abortRun('Tcond with BCs =/= {1,2} not implemented')
         endif

      end if

      if ( rank == 0 ) write(output_unit,*) '! Temperature Boundary Condition at the CMB', topt_bc
      if ( rank == 0 ) write(output_unit,*) '! Temperature Boundary Condition at the ICB', bott_bc

      tmpsc0 = epsc0*pr*4.0_cp*third*pi*(r_cmb**3-r_icb**3)
      if ( rank == 0 ) write(output_unit,*) '! Total volumetric buoyancy Sources = ', tmpsc0

      !call logWrite('! Sources introduced to balance surface heat flux!')
      !do n=1,2
      !   if ( n == 1 ) n_out=output_unit
      !   if ( n == 2 ) n_out=n_log_file
      !   write(n_out,*) ''
      !   write(n_out, '(''!      Top Flux='',ES16.6)') topt_bc
      !   write(n_out, '(''!      Bot Flux='',ES16.6)') bott_bc
      !   write(n_out,*) ' ! Warning: Sources introduced to balance surface heat flux'
      !   write(n_out, '(''!      epsrc0*pr='',ES16.6)') epsc0
      !   write(n_out,*) ''
      !end do
      !call logWrite(message)

   end subroutine get_conducting_state
!------------------------------------------------------------------------------
end module radial_functions
