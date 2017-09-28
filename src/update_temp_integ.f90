module update_temp_integ

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: one, zero, four, ci, three, half
   use pre_calculations, only: opr
   use namelists, only: kbott, ktopt, tadvz_fac, ra
   use radial_functions, only: rscheme, or1, or2, dtcond, tcond, beta, &
       &                       rgrav, r
   use blocking, only: nMstart, nMstop
   use truncation, only: n_r_max, idx2m
   use radial_der, only: get_ddr, get_dr
   use fields, only: work_Mloc
   use algebra, only: prepare_bordered_mat, solve_bordered_mat
   use useful, only: abortRun, roll
   use time_schemes, only: type_tscheme
   use chebsparselib, only: intcheb2, intcheb1rmult1, rmult2, intcheb2rmult2

   implicit none
   
   private

   integer, parameter :: klA4=4
   integer, parameter :: kuA4=4
   integer, parameter :: klB=4
   integer, parameter :: kuB=4
   integer, parameter :: n_boundaries=2
   integer, parameter :: n_bands_Amat=2*klA4+kuA4+1
   integer, parameter :: n_bands_Bmat = klB+kuB+1

   integer :: lenA4
   real(cp), allocatable :: A1mat(:,:)
   real(cp), allocatable :: A2mat(:,:)
   real(cp), allocatable :: Bmat(:,:)
   real(cp), allocatable :: A3mat(:,:,:)
   real(cp), allocatable :: A4mat(:,:,:)
   logical,  allocatable :: lTmat(:)
   integer,  allocatable :: pivotA4(:, :)
   integer,  allocatable :: pivotA1(:)
   complex(cp), allocatable :: rhs(:)

   public :: update_temp_int, initialize_temp_integ, finalize_temp_integ, &
   &         get_temp_rhs_imp

contains

   subroutine initialize_temp_integ

      lenA4 = n_r_max-n_boundaries

      allocate( A1mat(n_boundaries, n_boundaries) )
      allocate( A2mat(n_boundaries, lenA4) )
      allocate( A3mat(lenA4, n_boundaries,nMstart:nMstop) )
      allocate( A4mat(n_bands_Amat,lenA4,nMstart:nMstop) )
      allocate( Bmat(n_bands_Bmat,n_r_max) )

      A1mat(:,:)  =0.0_cp
      A2mat(:,:)  =0.0_cp
      A3mat(:,:,:)=0.0_cp
      A4mat(:,:,:)=0.0_cp
      Bmat(:,:)   =0.0_cp

      call get_rhs_mat(Bmat)

      allocate( lTmat(nMstart:nMstop) )
      lTmat(:)=.false.

      allocate( pivotA4(n_r_max, nMstart:nMstop) )
      allocate( pivotA1(n_r_max) )
      allocate( rhs(n_r_max) )
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_COMPLEX

   end subroutine initialize_temp_integ
!------------------------------------------------------------------------------
   subroutine finalize_temp_integ

      deallocate( rhs )
      deallocate( lTmat)

   end subroutine finalize_temp_integ
!------------------------------------------------------------------------------
   subroutine update_temp_int(us_Mloc, temp_Mloc, dtemp_Mloc, dVsT_Mloc,    &
              &           dtemp_exp_Mloc, dtemp_imp_Mloc, buo_imp_Mloc, &
              &           tscheme, lMat, l_roll_imp, l_log_next)

      !-- Input variables
      type(type_tscheme), intent(in) :: tscheme
      logical,            intent(in) :: lMat
      logical,            intent(in) :: l_roll_imp
      logical,            intent(in) :: l_log_next
      complex(cp),        intent(in) :: us_Mloc(nMstart:nMstop, n_r_max)

      !-- Output variables
      complex(cp), intent(out) :: temp_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(out) :: dtemp_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(inout) :: dtemp_exp_Mloc(nMstart:nMstop,n_r_max,tscheme%norder_exp)
      complex(cp), intent(inout) :: dtemp_imp_Mloc(nMstart:nMstop,n_r_max,tscheme%norder_imp-1)
      complex(cp), intent(inout) :: buo_imp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: dVsT_Mloc(nMstart:nMstop, n_r_max)

      !-- Local variables
      real(cp) :: rhsr(n_r_max), rhsi(n_r_max), tmp(n_r_max)
      integer :: n_r, n_m, n_r_out, m, n_o

      if ( lMat ) lTMat(:)=.false.

      !-- Assemble first buoyancy part from T^{n}
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m /= 0 ) then
               buo_imp_Mloc(n_m,n_r)=-tscheme%wimp_lin(2)*rgrav(n_r)*or1(n_r) &
               &                      *ra*opr*ci*real(m,cp)*temp_Mloc(n_m,n_r)
            end if
         end do
      end do

      !-- Finish calculation of advection
      call get_dr( dVsT_Mloc, work_Mloc, nMstart, nMstop, n_r_max, &
           &       rscheme, nocopy=.true. )

      !-- Finish calculation of the explicit part for current time step
      do n_r=1,n_r_max
         do n_m=nMstart, nMstop
            dtemp_exp_Mloc(n_m,n_r,1)=dtemp_exp_Mloc(n_m,n_r,1)   &
            &                     -or1(n_r)*work_Mloc(n_m,n_r)    &
            &                     -us_Mloc(n_m,n_r)*(dtcond(n_r)- &
            &                     tadvz_fac*beta(n_r)*tcond(n_r))
         end do
      end do

      !-- Calculation of the implicit part
      call get_temp_rhs_imp(temp_Mloc, dtemp_Mloc, tscheme%wimp_lin(2), &
           &                dtemp_imp_Mloc(:,:,1))


      do n_m=nMstart, nMstop

         m = idx2m(n_m)
         
         if ( .not. lTmat(n_m) ) then
            call get_lhs_mats(tscheme, m, A1mat(:,:), A2mat(:,:), &
                 &           A3mat(:,:,n_m), A4mat(:,:,n_m),     &
                 &           pivotA1(:), pivotA4(:,n_m))
            lTmat(n_m)=.true.
         end if

         do n_r=2,n_r_max-1
            do n_o=1,tscheme%norder_imp-1
               if ( n_o == 1 ) then
                  rhs(n_r)=tscheme%wimp(n_o+1)*dtemp_imp_Mloc(n_m,n_r,n_o)
               else
                  rhs(n_r)=rhs(n_r)+tscheme%wimp(n_o+1)*dtemp_imp_Mloc(n_m,n_r,n_o)
               end if
            end do
            do n_o=1,tscheme%norder_exp
               rhs(n_r)=rhs(n_r)+tscheme%wexp(n_o)*dtemp_exp_Mloc(n_m,n_r,n_o)
            end do
         end do

         do n_r=1,n_r_max
            rhsr(n_r)=real(rhs(n_r))
            rhsi(n_r)=aimag(rhs(n_r))
         end do
         call rscheme%costf1(rhsr,n_r_max)
         call rscheme%costf1(rhsi,n_r_max)
         tmp(:)=rhsr(:)
         call dgbmv('N', n_r_max, n_r_max, klB, kuB, one, Bmat, n_bands_Bmat, &
              &     tmp, 1, 0.0_cp, rhsr, 1)
         tmp(:)=rhsi(:)
         call dgbmv('N', n_r_max, n_r_max, klB, kuB, one, Bmat, n_bands_Bmat, &
              &     tmp, 1, 0.0_cp, rhsi, 1)
         rhsr(1)=0.0_cp
         rhsi(2)=0.0_cp

         call solve_bordered_mat(A1mat,A2mat,A3mat(:,:,n_m),A4mat(:,:,n_m), &
              &                  n_boundaries,lenA4, klA4, kuA4, pivotA1,   &
              &                  pivotA4(:,n_m), rhsr, n_r_max)

         call solve_bordered_mat(A1mat,A2mat,A3mat(:,:,n_m),A4mat(:,:,n_m), &
              &                  n_boundaries,lenA4, klA4, kuA4, pivotA1,   &
              &                  pivotA4(:,n_m), rhsi, n_r_max)

         do n_r=1,n_r_max
            rhs(n_r)=cmplx(rhsr(n_r), rhsi(n_r), kind=cp)
         end do

         do n_r_out=1,rscheme%n_max
            temp_Mloc(n_m, n_r_out)=rhs(n_r_out)
         end do

      end do

      !-- set cheb modes > rscheme%n_max to zero (dealiazing)
      do n_r_out=rscheme%n_max+1,n_r_max
         do n_m=nMstart,nMstop
            temp_Mloc(n_m,n_r_out)=zero
         end do
      end do

      !-- Bring temperature back to physical space
      call rscheme%costf1(temp_Mloc, nMstart, nMstop, n_r_max)

      !print*, 'T[n+1]', sum(abs(temp_Mloc))

      !-- Assemble second buoyancy part from T^{n+1}
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m /= 0 ) then
               buo_imp_Mloc(n_m,n_r)=            buo_imp_Mloc(n_m,n_r)-&
               &               tscheme%wimp_lin(1)*rgrav(n_r)*or1(n_r) &
               &               *ra*opr*ci*real(m,cp)*temp_Mloc(n_m,n_r)
            end if
         end do
      end do

      !-- Roll the arrays before filling again the first block
      call roll(dtemp_exp_Mloc, nMstart, nMstop, n_r_max, tscheme%norder_exp)
      if ( l_roll_imp ) then
         call roll(dtemp_imp_Mloc, nMstart, nMstop, n_r_max, tscheme%norder_imp-1)
      end if

      !-- In case log is needed on the next iteration, recalculate dT/dr
      if ( l_log_next ) then
         call get_dr(temp_Mloc, dtemp_Mloc, nMstart, nMstop, n_r_max, rscheme)
      end if

   end subroutine update_temp_int
!------------------------------------------------------------------------------
   subroutine get_temp_rhs_imp(temp_Mloc, dtemp_Mloc, wimp, dtemp_imp_Mloc_last)

      !-- Input variables
      complex(cp), intent(in) :: temp_Mloc(nMstart:nMstop,n_r_max)
      real(cp),    intent(in) :: wimp

      !-- Output variable
      complex(cp), intent(out) :: dtemp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(out) :: dtemp_imp_Mloc_last(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_m, m
      real(cp) :: dm2

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            dtemp_imp_Mloc_last(n_m,n_r)=temp_Mloc(n_m,n_r)
         end do
      end do

      if ( wimp /= 0.0_cp ) then

         call get_ddr(temp_Mloc, dtemp_Mloc, work_Mloc, nMstart, nMstop, &
              &       n_r_max, rscheme)
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               dm2 = real(m,cp)*real(m,cp)
               dtemp_imp_Mloc_last(n_m,n_r)=dtemp_imp_Mloc_last(n_m,n_r) &
               &                          +wimp*opr*( work_Mloc(n_m,n_r) &
               &                         +or1(n_r)*  dtemp_Mloc(n_m,n_r) &
               &                     -dm2*or2(n_r)*   temp_Mloc(n_m,n_r) )
            end do
         end do

      end if

   end subroutine get_temp_rhs_imp
!------------------------------------------------------------------------------
   subroutine get_lhs_mats(tscheme, m, A1, A2, A3, A4, pivA1, &
              &            pivA4)

      !-- Input variables
      type(type_tscheme), intent(in) :: tscheme        ! time step
      integer,            intent(in) :: m

      !-- Output variables
      real(cp), intent(out) :: A1(n_boundaries,n_boundaries)
      real(cp), intent(out) :: A2(n_boundaries,lenA4)
      real(cp), intent(out) :: A3(lenA4,n_boundaries)
      real(cp), intent(out) :: A4(n_bands_Amat,lenA4)
      integer,  intent(out) :: pivA1(n_boundaries)
      integer,  intent(out) :: pivA4(lenA4)

      !-- Local variables
      real(cp) :: stencilA4(klA4+kuA4+1)
      integer :: n_r, info, i_r, n_band
      real(cp) :: dm2, a, b

      dm2 = real(m,cp)*real(m,cp)

      a = half*(r(1)-r(n_r_max))
      b = half*(r(1)+r(n_r_max))

      !-- Fill A4 banded-block
      do n_r=1,lenA4
         i_r = n_r+n_boundaries

         !-- Define the equations
         stencilA4 = intcheb2rmult2(a,b,i_r-1,klA4+kuA4+1)-     &
         &                       tscheme%wimp_lin(1)*opr*(      &
         &                           rmult2(a,b,klA4+kuA4+1)    &
         &      -three*intcheb1rmult1(a,b,i_r-1,klA4+kuA4+1)    &
         &            +(one-dm2)*intcheb2(a,i_r-1,klA4+kuA4+1) )

         !-- Roll the array for band storage
         do n_band=1,klA4+kuA4+1
            if ( n_r+kuA4+1-n_band <= lenA4 .and. n_r+kuA4+1-n_band >= 1 ) then
               A4(klA4+n_band,n_r+kuA4+1-n_band) = rscheme%rnorm*stencilA4(n_band)
            end if
         end do
      end do

      !-- Fill A3
      do n_r=1,n_boundaries
         i_r = n_r+n_boundaries
         stencilA4 = intcheb2rmult2(a,b,i_r-1,klA4+kuA4+1)-     &
         &                       tscheme%wimp_lin(1)*opr*(      &
         &                           rmult2(a,b,klA4+kuA4+1)    &
         &      -three*intcheb1rmult1(a,b,i_r-1,klA4+kuA4+1)    &
         &            +(one-dm2)*intcheb2(a,i_r-1,klA4+kuA4+1) )

         !-- Only the lower bands can contribute to the matrix A3
         do n_band=1,klA4
            if ( n_r <= n_band .and. n_r+n_boundaries-n_band >= 1 ) then
               A3(n_r,n_r+n_boundaries-n_band) = rscheme%rnorm*stencilA4(kuA4+1+n_band)
            end if
         end do
      end do

      !-- Add the tau Lines for boundary conditions into A1 and A2
      do n_r=1,n_r_max
         if ( n_r <= n_boundaries) then
            if ( ktopt == 1 ) then
               A1(1,n_r)=rscheme%rnorm*rscheme%rMat(1,n_r)
            else
               A1(1,n_r)=rscheme%rnorm*rscheme%drMat(1,n_r)
            end if
            if ( kbott == 1 ) then
               A1(2,n_r)=rscheme%rnorm*rscheme%rMat(n_r_max,n_r)
            else
               A1(2,n_r)=rscheme%rnorm*rscheme%drMat(n_r_max,n_r)
            end if
         else
            if ( ktopt == 1 ) then
               A2(1,n_r-n_boundaries)=rscheme%rnorm*rscheme%rMat(1,n_r)
            else
               A2(1,n_r-n_boundaries)=rscheme%rnorm*rscheme%drMat(1,n_r)
            end if
            if ( kbott == 1 ) then
               A2(2,n_r-n_boundaries)=rscheme%rnorm*rscheme%rMat(n_r_max,n_r)
            else
               A2(2,n_r-n_boundaries)=rscheme%rnorm*rscheme%drMat(n_r_max,n_r)
            end if
         end if
      end do

      !-- Cheb factor for boundary conditions
      do n_r=1,n_boundaries
         A1(n_r,1)    =rscheme%boundary_fac*A1(n_r,1)
         A2(n_r,lenA4)=rscheme%boundary_fac*A2(1,lenA4)
      end do


      !-- LU factorisation
      call prepare_bordered_mat(A1,A2,A3,A4,n_boundaries,lenA4, &
           &                    klA4,kuA4, pivA1, pivA4)

   end subroutine get_lhs_mats
!------------------------------------------------------------------------------
   subroutine get_rhs_mat(Bmat)

      !-- Output variable
      real(cp), intent(out) :: Bmat(n_bands_Bmat,n_r_max)

      !-- Local variables
      real(cp) :: stencilB(klB+kuB+1)
      real(cp) :: a, b
      integer :: n_band, n_r, i_r

      a = half*(r(1)-r(n_r_max))
      b = half*(r(1)+r(n_r_max))


      !-- Fill right-hand side matrix
      do n_r=1,lenA4
         i_r = n_r+n_boundaries

         !-- Define right-hand side equations
         stencilB = intcheb2rmult2(a,b,i_r-1,klB+kuB+1)

         !-- Roll array for band storage
         do n_band=1,klB+kuB+1
            if ( i_r+kuB+1-n_band <= n_r_max .and. i_r+kuB+1-n_band >= 1 ) then
               Bmat(n_band,i_r+kuB+1-n_band) = rscheme%rnorm*stencilB(n_band)
            end if
         end do
      end do

   end subroutine get_rhs_mat
!------------------------------------------------------------------------------
end module update_temp_integ
