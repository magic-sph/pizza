module horizontal
   !
   ! This module implements the calculation of m-dependent array
   ! such as hyperdiffusion or heat flux patterns
   !

   use precision_mod
   use parallel_mod
   use mem_alloc
   use constants, only: zero, one, pi, two, half
   use blocking, only: nMstart, nMstop, m_balance, nm_per_rank, lmStart, lmStop
   use truncation, only: idx2m, n_m_max
   use truncation_3D, only: n_theta_max, lm_max, l_max, minc_3D, m_max_3D
   use blocking_lm, only: lmP2l, lmP2lm, lm2l, lm2m
   use namelists, only: hdif_m, hdif_exp, hdif_vel, hdif_temp, tag,   &
       &                t_bot, t_top, xi_bot, xi_top, hdif_comp,      &
       &                l_heat, l_chem, l_3D, l_mag_3D, hdif_mag

   implicit none

   private

   !-- Arrays for Inhomogeneous temperature B.Cs
   complex(cp), allocatable, public :: bott_Mloc(:)
   complex(cp), allocatable, public :: topt_Mloc(:)
   complex(cp), allocatable, public :: botxi_Mloc(:)
   complex(cp), allocatable, public :: topxi_Mloc(:)

   real(cp), public, allocatable :: hdif_V(:)
   real(cp), public, allocatable :: hdif_T(:)
   real(cp), public, allocatable :: hdif_Xi(:)
   real(cp), public, allocatable :: hdif_B(:)   !-- WARNING:: Totally wrong -> to be modified!!

   !-- Arrays depending on theta (colatitude):
   real(cp), public, allocatable :: theta(:)
   real(cp), public, allocatable :: sint(:)
   real(cp), public, allocatable :: cost(:)
   real(cp), public, allocatable :: osint1(:)
   real(cp), public, allocatable :: osint2(:)

   !-- Arrays depending on l and m:
   complex(cp), public, allocatable :: dPhi(:)
   complex(cp), public, allocatable :: dPhi0(:)
   real(cp), public, allocatable :: dLh(:)
   real(cp), public, allocatable :: dTheta1S(:),dTheta1A(:)
   real(cp), public, allocatable :: dTheta2S(:),dTheta2A(:)
   real(cp), public, allocatable :: dTheta3S(:),dTheta3A(:)
   real(cp), public, allocatable :: dTheta4S(:),dTheta4A(:)
   real(cp), public, allocatable :: D_m(:),D_l(:),D_lP1(:)


   public :: initialize_mfunctions, finalize_mfunctions, mfunctions, &
   &         spherical_functions

contains

   subroutine initialize_mfunctions

      allocate( hdif_V(nMstart:nMstop) )
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_DEF_REAL
      if ( l_heat ) then
         allocate( bott_Mloc(nMStart:nMstop) )
         allocate( topt_Mloc(nMStart:nMstop) )
         allocate( hdif_T(nMstart:nMstop) )
         bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_DEF_REAL+&
         &                 2*(nMstop-nMstart+1)*SIZEOF_DEF_COMPLEX
      end if
      if ( l_chem ) then
         allocate( botxi_Mloc(nMStart:nMstop) )
         allocate( topxi_Mloc(nMStart:nMstop) )
         allocate( hdif_Xi(nMstart:nMstop) )
         bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_DEF_REAL+&
         &                 2*(nMstop-nMstart+1)*SIZEOF_DEF_COMPLEX
      end if
      if ( l_mag_3D ) then
         allocate( hdif_B(lmStart:lmStop) )
         bytes_allocated = bytes_allocated+(lmStop-lmStart+1)*SIZEOF_DEF_REAL
      end if

      if ( l_3D ) then

         !-- add theta functions
         allocate( theta(n_theta_max), sint(n_theta_max), cost(n_theta_max) )
         allocate( osint1(n_theta_max), osint2(n_theta_max) )
         bytes_allocated = bytes_allocated+5*n_theta_max*SIZEOF_DEF_REAL

         !-- Arrays depending on l and m:
         allocate( dPhi(lm_max) )
         allocate( dPhi0(lm_max) )
         allocate( dLh(lm_max) )
         allocate( dTheta1S(lm_max),dTheta1A(lm_max) )
         allocate( dTheta2S(lm_max),dTheta2A(lm_max) )
         allocate( dTheta3S(lm_max),dTheta3A(lm_max) )
         allocate( dTheta4S(lm_max),dTheta4A(lm_max) )
         allocate( D_m(lm_max),D_l(lm_max),D_lP1(lm_max) )
         bytes_allocated = bytes_allocated+14*lm_max*SIZEOF_DEF_REAL

      end if

   end subroutine initialize_mfunctions
!--------------------------------------------------------------------------------
   subroutine finalize_mfunctions

      if ( l_heat ) deallocate( bott_Mloc, topt_Mloc, hdif_T )
      if ( l_chem ) deallocate( botxi_Mloc, topxi_Mloc, hdif_Xi )
      deallocate( hdif_V )
      if ( l_mag_3D ) deallocate( hdif_B )

      if ( l_3D ) then
         !-- theta functions
         deallocate( sint, cost, osint1, osint2 )
         !-- lm functions
         deallocate( dPhi, dPhi0, dLh, dTheta1S, dTheta1A )
         deallocate( dTheta2S, dTheta2A, dTheta3S, dTheta3A, dTheta4S, dTheta4A )
         deallocate( D_m, D_l, D_lP1 )
      end if

   end subroutine finalize_mfunctions
!--------------------------------------------------------------------------------
   subroutine mfunctions

      !-- Local variables
      integer :: n, m_bot, m_top
      integer :: n_m, m, n_p, file_handle
      integer :: displs(0:n_procs-1), recvcounts(0:n_procs-1)
      real(cp) :: eps, tr_bot, ti_bot, tr_top, ti_top
      real(cp) :: hdif_T_global(n_m_max), hdif_V_global(n_m_max)
      real(cp) :: hdif_Xi_global(n_m_max)

      eps = 10.0_cp*epsilon(one)
      if ( abs(hdif_temp) > eps .or. abs(hdif_vel) > eps .or. &
         & abs(hdif_comp) > eps ) then

         do n_m=nMstart,nMstop
            m = idx2m(n_m)

            if ( m > hdif_m ) then ! This is the Nataf & Schaffer (2015) form
               if ( l_heat ) hdif_T(n_m) = hdif_temp**(m-hdif_m)
               if ( l_chem ) hdif_Xi(n_m) = hdif_comp**(m-hdif_m)
               hdif_V(n_m) = hdif_vel**(m-hdif_m)
            else
               if ( l_heat ) hdif_T(n_m) = one
               if ( l_chem ) hdif_Xi(n_m) = one
               hdif_V(n_m) = one
            end if

         end do

         !-- Gather the profiles on rank 0 to write the profiles in hdif.TAG
         do n_p=0,n_procs-1
            recvcounts(n_p)=m_balance(n_p)%n_per_rank
         end do
         displs(0)=0
         do n_p=1,n_procs-1
            displs(n_p)=displs(n_p-1)+recvcounts(n_p-1)
         end do
         if ( l_heat ) then
            call MPI_GatherV(hdif_T, nm_per_rank, MPI_DEF_REAL,  &
                 &           hdif_T_global, recvcounts, displs,  &
                 &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
         end if
         if ( l_chem ) then
            call MPI_GatherV(hdif_Xi, nm_per_rank, MPI_DEF_REAL,  &
                 &           hdif_Xi_global, recvcounts, displs,  &
                 &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
         end if
         call MPI_GatherV(hdif_V, nm_per_rank, MPI_DEF_REAL,  &
              &           hdif_V_global, recvcounts, displs,  &
              &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)

         !-- Now rank0 writes an output file
         if ( rank== 0 ) then
            open(newunit=file_handle, file='hdif.'//tag, status='new')

            do n_m=1,n_m_max
               m = idx2m(n_m)
               write(file_handle, '(I5, 3es16.8)') m, hdif_T_global(n_m), &
               &                                      hdif_Xi_global(n_m),&
               &                                      hdif_V_global(n_m)
            end do
            close(file_handle)
         end if

      else

         hdif_V(:) = one
         if ( l_heat ) hdif_T(:) = one
         if ( l_chem ) hdif_Xi(:) = one

      end if

      !-- Build matrices Inhomogeneous B.Cs
      !-- Inspired by MagIC (J.Wicht; T.Gastine; A.Barik; R.Raynaud; B.Putigny;
      !--                    R.Yadav; ; L.Duarte; B.dintrans; T.Schwaiger)!
      if ( l_heat ) then
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            bott_Mloc(n_m)=zero
            topt_Mloc(n_m)=zero
            do n=1,size(t_bot)/3
               m_bot =int(t_bot(3*n-2))
               tr_bot=t_bot(3*n-1)
               ti_bot=t_bot(3*n)
               m_top =int(t_top(3*n-2))
               tr_top=t_top(3*n-1)
               ti_top=t_top(3*n)
               if ( m_bot == m .and. &
                   cmplx(tr_bot,ti_bot,kind=cp) /= zero ) then
                  if ( m == 0 ) ti_bot=0.0_cp
                  bott_Mloc(n_m)=cmplx(tr_bot,ti_bot,kind=cp)
               end if
               if ( m_top == m .and. &
                   cmplx(tr_top,ti_top,kind=cp) /= zero ) then
                  if ( m == 0 ) ti_top=0.0_cp
                  topt_Mloc(n_m)=cmplx(tr_top,ti_top,kind=cp)
               end if
            end do
         end do
      end if

      if ( l_chem ) then
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            botxi_Mloc(n_m)=zero
            topxi_Mloc(n_m)=zero
            do n=1,size(t_bot)/3
               m_bot =int(t_bot(3*n-2))
               tr_bot=t_bot(3*n-1)
               ti_bot=t_bot(3*n)
               m_top =int(t_top(3*n-2))
               tr_top=t_top(3*n-1)
               ti_top=t_top(3*n)
               if ( m_bot == m .and. &
                   cmplx(tr_bot,ti_bot,kind=cp) /= zero ) then
                  if ( m == 0 ) ti_bot=0.0_cp
                  botxi_Mloc(n_m)=cmplx(tr_bot,ti_bot,kind=cp)
               end if
               if ( m_top == m .and. &
                   cmplx(tr_top,ti_top,kind=cp) /= zero ) then
                  if ( m == 0 ) ti_top=0.0_cp
                  topxi_Mloc(n_m)=cmplx(tr_top,ti_top,kind=cp)
               end if
            end do
         end do

      end if

   end subroutine mfunctions
!------------------------------------------------------------------------------
   subroutine spherical_functions

      !-- Local variables
      integer :: l, m
      integer :: n_theta, lm, n_p, file_handle
      integer :: displs(0:n_procs-1), recvcounts(0:n_procs-1)
      real(cp) :: eps, hdif_B_global(lm_max)
      real(cp) :: clm(0:l_max+1,0:l_max+1), colat
      real(cp) :: theta_ord(n_theta_max), gauss(n_theta_max)

      !-- Calculate grid points and weights for the
      call gauleg(-one,one,theta_ord,gauss,n_theta_max)

      !-- cos(theta) and some basic theta functions:
      do n_theta=1,n_theta_max  ! Loop over colat in NHS
         colat =theta_ord(n_theta)
         theta(n_theta) =colat
         sint(n_theta)  =sin(colat)
         cost(n_theta)  =cos(colat)
         osint1(n_theta)=one/sin(colat)
         osint2(n_theta)=one/(sin(colat)*sin(colat))
      end do

      !-- Build arrays depending on degree l and order m
      do m=0,m_max_3D,minc_3D  ! Build auxiliary array clm
         do l=m,l_max+1
            clm(l,m)=sqrt( real((l+m)*(l-m),cp) / real((2*l-1)*(2*l+1),cp) )
         end do
      end do

      do lm=1,lm_max
         l=lm2l(lm)
         m=lm2m(lm)

         !-- Help arrays:
         D_l(lm)  =real(l,cp)
         D_lP1(lm)=real(l+1,cp)
         D_m(lm)  =real(m,cp)

         !---- Operators for derivatives:
         !-- Phi derivate:
         dPhi(lm)=cmplx(0.0_cp,real(m,cp),cp)
         if ( l < l_max ) then
            dPhi0(lm)    =cmplx(0.0_cp,real(m,cp),cp)
         else
            dPhi0(lm)    =zero
         end if
         !-- Negative horizontal Laplacian *r^2
         dLh(lm)     =real(l*(l+1),cp)                 ! = qll1
         !-- Operator ( 1/sin(theta) * d/d theta * sin(theta)**2 )
         dTheta1S(lm)=real(l+1,cp)        *clm(l,m)    ! = qcl1
         dTheta1A(lm)=real(l,cp)          *clm(l+1,m)  ! = qcl
         !-- Operator ( sin(theta) * d/d theta )
         dTheta2S(lm)=real(l-1,cp)        *clm(l,m)    ! = qclm1
         dTheta2A(lm)=real(l+2,cp)        *clm(l+1,m)  ! = qcl2
         !-- Operator ( sin(theta) * d/d theta + cos(theta) dLh )
         dTheta3S(lm)=real((l-1)*(l+1),cp)*clm(l,m)    ! = q0l1lm1(lm)
         dTheta3A(lm)=real(l*(l+2),cp)    *clm(l+1,m)  ! = q0cll2(lm)
         !-- Operator ( 1/sin(theta) * d/d theta * sin(theta)**2 ) * dLh
         dTheta4S(lm)=dTheta1S(lm)*real((l-1)*l,cp)
         dTheta4A(lm)=dTheta1A(lm)*real((l+1)*(l+2),cp)
      enddo

      eps = 10.0_cp*epsilon(one)
      if ( abs(hdif_mag) > eps ) then

         do lm=lmStart,lmStop
            l=lm2l(lm)
            m=lm2l(lm)

            if ( l > hdif_m ) then ! This is the Nataf & Schaffer (2015) form
               if ( l_mag_3D ) hdif_B(lm) = hdif_mag**(l-hdif_m)
            else
               if ( l_mag_3D ) hdif_B(lm) = one
            end if

         end do

         !-- Gather the profiles on rank 0 to write the profiles in hdif.TAG
         do n_p=0,n_procs-1
            recvcounts(n_p)=m_balance(n_p)%n_per_rank
         end do
         displs(0)=0
         do n_p=1,n_procs-1
            displs(n_p)=displs(n_p-1)+recvcounts(n_p-1)
         end do
         if ( l_mag_3D ) then
            call MPI_GatherV(hdif_B, nm_per_rank, MPI_DEF_REAL,  &
                 &           hdif_B_global, recvcounts, displs,  &
                 &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)
         end if

         !-- Now rank0 writes an output file
         if ( rank== 0 ) then
            open(newunit=file_handle, file='hdif.'//tag, status='new')

            do lm=1,lm_max
               l=lm2l(lm)
               m=lm2l(lm)
               write(file_handle, '(I5, 3es16.8)') l, m, hdif_B_global(lm)
            end do
            close(file_handle)
         end if

      else

         if ( l_mag_3D ) hdif_B(:) = one

      end if

   end subroutine spherical_functions
!--------------------------------------------------------------------------------
   subroutine gauleg(sinThMin,sinThMax,theta_ord,gauss,n_th_max)
      !
      ! Subroutine is based on a NR code.
      ! Calculates N zeros of legendre polynomial P(l=N) in
      ! the interval [sinThMin,sinThMax].
      ! Zeros are returned in radiants theta_ord(i)
      ! The respective weights for Gauss-integration are given in gauss(i).
      !

      !-- Input variables:
      real(cp), intent(in) :: sinThMin ! lower bound in radiants
      real(cp), intent(in) :: sinThMax ! upper bound in radiants
      integer,  intent(in) :: n_th_max ! desired maximum degree
    
      !-- Output variables:
      real(cp), intent(out) :: theta_ord(n_th_max) ! zeros cos(theta)
      real(cp), intent(out) :: gauss(n_th_max)     ! associated Gauss-Legendre weights
    
      !-- Local variables:
      integer :: m,i,j
      real(cp) :: sinThMean,sinThDiff,p1,p2,p3,pp,z,z1
      real(cp), parameter :: eps = 10.0_cp*epsilon(one)
    
      ! use symmetry
      m=(n_th_max+1)/2
    
      !-- Map on symmetric interval:
      sinThMean=half*(sinThMax+sinThMin)
      sinThDiff=half*(sinThMax-sinThMin)
    
      do i=1,m
         !----- Initial guess for zeros:
         z  = cos( pi*( (real(i,cp)-0.25_cp)/(real(n_th_max,cp)+half)) )
         z1 = z+10.0_cp*eps
     
         do while( abs(z-z1) > eps)
            !----- Use recurrence to calulate P(l=n_th_max,z=cos(theta))
            p1=one
            p2=0.0_cp
            ! do loop over degree !
            do j=1,n_th_max
               p3=p2
               p2=p1
               p1=( real(2*j-1,cp)*z*p2-real(j-1,cp)*p3 )/real(j,cp)
            end do
      
            !----- Newton method to refine zero: pp is derivative !
            pp=real(n_th_max,cp)*(z*p1-p2)/(z*z-one)
            z1=z
            z=z1-p1/pp
         end do
     
         !----- Another zero found
         theta_ord(i)           =acos(sinThMean+sinThDiff*z)
         theta_ord(n_th_max+1-i)=acos(sinThMean-sinThDiff*z)
         gauss(i)               =two*sinThDiff/((one-z*z)*pp*pp)
         gauss(n_th_max+1-i)    =gauss(i)
    
      end do
     
   end subroutine gauleg
!------------------------------------------------------------------------------
end module horizontal
