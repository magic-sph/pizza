module update_temp_integ

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: zero, one, ci, half
   use namelists, only: kbott, ktopt, tadvz_fac, BuoFac, r_cmb, r_icb, &
       &                TdiffFac, l_non_rot
   use radial_functions, only: rscheme, or1, or2, dtcond, tcond, rgrav, r
   use blocking, only: nMstart, nMstop
   use truncation, only: n_r_max, idx2m, n_cheb_max
   use radial_der, only: get_dr
   use fields, only: work_Mloc
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use matrix_types, only: type_bandmat_real, type_bordmat_real
   use chebsparselib, only: intcheb2rmult2, intcheb2rmult2lapl

   implicit none
   
   private

   integer, parameter :: n_boundaries=2

   logical,  allocatable :: lTmat(:)
   complex(cp), allocatable :: rhs(:)

   type(type_bordmat_real), allocatable :: LHS_mat(:)
   type(type_bandmat_real) :: RHSE_mat
   type(type_bandmat_real), allocatable :: RHSI_mat(:)
   real(cp), allocatable :: tempfac(:,:) ! Preconditon matrix

   public :: update_temp_int, initialize_temp_integ, finalize_temp_integ, &
   &         get_temp_rhs_imp_int

contains

   subroutine initialize_temp_integ

      integer :: n_m, m

      call RHSE_mat%initialize(4, 4, n_cheb_max)

      allocate( RHSI_mat(nMstart:nMstop) )
      allocate( LHS_mat(nMstart:nMstop) )

      do n_m=nMstart,nMstop
         call RHSI_mat(n_m)%initialize(4, 4, n_cheb_max)
         call LHS_mat(n_m)%initialize(4, 4, n_boundaries, n_cheb_max)
      end do

      call get_rhs_exp_mat(RHSE_mat)
      do n_m=nMstart,nMstop
         m = idx2m(n_m)
         call get_rhs_imp_mat(RHSI_mat(n_m), m)
      end do

      allocate( lTmat(nMstart:nMstop) )
      lTmat(:)=.false.
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_LOGICAL

      allocate( rhs(n_cheb_max) )
      bytes_allocated = bytes_allocated+n_cheb_max*SIZEOF_DEF_COMPLEX

      allocate( tempfac(n_cheb_max, nMstart:nMstop) )
      bytes_allocated = bytes_allocated + n_cheb_max*(nMstop-nMstart+1)* &
      &                 SIZEOF_DEF_REAL

   end subroutine initialize_temp_integ
!------------------------------------------------------------------------------
   subroutine finalize_temp_integ

      !-- Local variables
      integer :: n_m

      call RHSE_mat%finalize()
      do n_m=nMstart,nMstop
         call RHSI_mat(n_m)%finalize()
         call LHS_mat(n_m)%finalize()
      end do
      deallocate( tempfac, LHS_mat, RHSI_mat, rhs, lTmat )

   end subroutine finalize_temp_integ
!------------------------------------------------------------------------------
   subroutine update_temp_int(psi_Mloc, temp_Mloc, dtemp_Mloc, dVsT_Mloc,    &
              &               dtemp_exp_Mloc, temp_old_Mloc, dtemp_imp_Mloc, &
              &               buo_imp_Mloc, tscheme, lMat, l_log_next)

      !-- Input variables
      type(type_tscheme), intent(in) :: tscheme
      logical,            intent(in) :: lMat
      logical,            intent(in) :: l_log_next
      complex(cp),        intent(in) :: psi_Mloc(nMstart:nMstop, n_r_max)

      !-- Output variables
      complex(cp), intent(out) :: temp_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(out) :: dtemp_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(inout) :: dtemp_exp_Mloc(nMstart:nMstop,n_r_max,tscheme%norder_exp)
      complex(cp), intent(inout) :: temp_old_Mloc(nMstart:nMstop,n_r_max,tscheme%norder_imp-1)
      complex(cp), intent(inout) :: dtemp_imp_Mloc(nMstart:nMstop,n_r_max,tscheme%norder_imp_lin-1)
      complex(cp), intent(inout) :: buo_imp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: dVsT_Mloc(nMstart:nMstop, n_r_max)

      !-- Local variables
      real(cp) :: h2
      integer :: n_r, n_m, m, n_cheb

      if ( lMat ) lTMat(:)=.false.

      !-- Assemble first buoyancy part from T^{n}
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m /= 0 ) then
               buo_imp_Mloc(n_m,n_r)=-tscheme%wimp_lin(2)*rgrav(n_r)*or1(n_r) &
               &                      *BuoFac*ci*real(m,cp)*temp_Mloc(n_m,n_r)
            end if
         end do
      end do

      !-- Finish calculation of advection
      call get_dr( dVsT_Mloc, work_Mloc, nMstart, nMstop, n_r_max, &
           &       rscheme, nocopy=.true. )

      !-- Finish calculation of the explicit part for current time step
      if ( l_non_rot ) then
         do n_r=1,n_r_max
            do n_m=nMstart, nMstop
               m = idx2m(n_m)
               dtemp_exp_Mloc(n_m,n_r,1)=dtemp_exp_Mloc(n_m,n_r,1)   &
               &                     -or1(n_r)*work_Mloc(n_m,n_r)    &
               &             -ci*real(m,cp)*or1(n_r)*dtcond(n_r)*    &
               &                                psi_Mloc(n_m,n_r)
            end do
         end do
      else ! this is rotating
         do n_r=1,n_r_max
            h2 = r_cmb*r_cmb-r(n_r)*r(n_r)
            do n_m=nMstart, nMstop
               m = idx2m(n_m)
               dtemp_exp_Mloc(n_m,n_r,1)=dtemp_exp_Mloc(n_m,n_r,1)   &
               &                     -or1(n_r)*work_Mloc(n_m,n_r)    &
               &          -ci*real(m,cp)*(h2*or1(n_r)*dtcond(n_r)+   &
               &                            tadvz_fac* tcond(n_r))*  &
               &                                psi_Mloc(n_m,n_r)
            end do
         end do
      end if

      !-- Transform the explicit part to chebyshev space
      call rscheme%costf1(dtemp_exp_Mloc(:,:,1), nMstart, nMstop, n_r_max)

      !-- Matrix-vector multiplication by the operator \int\int r^2 .
      do n_m=nMstart,nMstop
         do n_cheb=1,n_cheb_max
            rhs(n_cheb)=dtemp_exp_Mloc(n_m,n_cheb,1)
         end do

         call RHSE_mat%mat_vec_mul(rhs)

         rhs(1)=zero
         rhs(2)=zero
         do n_cheb=1,n_cheb_max
            dtemp_exp_Mloc(n_m,n_cheb,1)=rhs(n_cheb)
         end do
         !-- Pad with zeros
         do n_cheb=n_cheb_max+1,n_r_max
            dtemp_exp_Mloc(n_m,n_cheb,1)=zero
         end do
      end do

      !-- Calculation of the implicit part
      call get_temp_rhs_imp_int(temp_Mloc, temp_old_Mloc(:,:,1), &
           &                    dtemp_imp_Mloc(:,:,1), tscheme%l_calc_lin_rhs)

      !-- Now assemble the right hand side and store it in work_Mloc
      call tscheme%set_imex_rhs(work_Mloc, dtemp_imp_Mloc, dtemp_exp_Mloc, &
           &                    temp_old_Mloc, nMstart, nMstop, n_r_max)

      do n_m=nMstart, nMstop

         m = idx2m(n_m)
         
         if ( .not. lTmat(n_m) ) then
            call get_lhs_mat( tscheme, LHS_mat(n_m), tempfac(:,n_m), m )
            lTmat(n_m)=.true.
         end if

         do n_cheb=1,n_cheb_max
            rhs(n_cheb)=work_Mloc(n_m,n_cheb)
         end do

         !-- Multiply rhs by precond matrix
         do n_cheb=1,n_cheb_max
            rhs(n_cheb)=rhs(n_cheb)*tempfac(n_cheb,n_m)
         end do

         call LHS_mat(n_m)%solve(rhs, n_cheb_max)

         do n_cheb=1,n_cheb_max
            temp_Mloc(n_m, n_cheb)=rhs(n_cheb)
         end do

      end do

      !-- set cheb modes > n_cheb_max to zero (dealiazing)
      do n_cheb=n_cheb_max+1,n_r_max
         do n_m=nMstart,nMstop
            temp_Mloc(n_m,n_cheb)=zero
         end do
      end do

      !-- Bring temperature back to physical space
      call rscheme%costf1(temp_Mloc, nMstart, nMstop, n_r_max)

      !-- Finish assembling buoyancy 
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m /= 0 ) then
               buo_imp_Mloc(n_m,n_r)=            buo_imp_Mloc(n_m,n_r)-&
               &               tscheme%wimp_lin(1)*rgrav(n_r)*or1(n_r) &
               &               *BuoFac*ci*real(m,cp)*temp_Mloc(n_m,n_r)
            end if
         end do
      end do

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dtemp_imp_Mloc, dtemp_exp_Mloc, temp_old_Mloc, &
           &                   nMstart, nMstop, n_r_max)

      !-- In case log is needed on the next iteration, recalculate dT/dr
      !-- This is needed to estimate the heat fluxes
      if ( l_log_next ) then
         call get_dr(temp_Mloc, dtemp_Mloc, nMstart, nMstop, n_r_max, rscheme)
      end if

   end subroutine update_temp_int
!------------------------------------------------------------------------------
   subroutine get_temp_rhs_imp_int(temp_Mloc, temp_old, dtemp_imp_Mloc_last,  &
              &                    l_calc_rhs_lin)

      !-- Input variables
      complex(cp), intent(in) :: temp_Mloc(nMstart:nMstop,n_r_max)
      logical,     intent(in) :: l_calc_rhs_lin

      !-- Output variable
      complex(cp), intent(out) :: temp_old(nMstart:nMstop,n_r_max)
      complex(cp), intent(out) :: dtemp_imp_Mloc_last(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_m, n_cheb

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            temp_old(n_m,n_r)=temp_Mloc(n_m,n_r)
         end do
      end do

      !-- Transform the implicit part to chebyshev space
      call rscheme%costf1(temp_old, nMstart, nMstop, n_r_max)

      !-- Matrix-vector multiplication by the operator \int\int r^2 .
      do n_m=nMstart,nMstop

         do n_cheb=1,n_cheb_max
            rhs(n_cheb)= temp_old(n_m,n_cheb)
         end do

         call RHSE_mat%mat_vec_mul(rhs)

         rhs(1)=zero
         rhs(2)=zero

         do n_cheb=1,n_cheb_max
            temp_old(n_m,n_cheb)=rhs(n_cheb)
         end do

         !-- Pad with zeros
         do n_cheb=n_cheb_max+1,n_r_max
            temp_old(n_m,n_cheb)=zero
         end do

      end do

      if ( l_calc_rhs_lin ) then

         !-- Copy temp_Mloc into work_Mloc
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               work_Mloc(n_m,n_r)=temp_Mloc(n_m,n_r)
            end do
         end do

         !-- Transform work_Mloc to Cheb space
         call rscheme%costf1(work_Mloc, nMstart, nMstop, n_r_max)

         !-- Matrix-vector multiplication by the LHS operator
         do n_m=nMstart,nMstop
            do n_cheb=1,n_cheb_max
               rhs(n_cheb)=work_Mloc(n_m,n_cheb)
            end do
            call RHSI_mat(n_m)%mat_vec_mul(rhs)
            rhs(1)=zero
            rhs(2)=zero
            do n_cheb=1,n_cheb_max
               work_Mloc(n_m,n_cheb)=rhs(n_cheb)
            end do
            !-- Pad with zeros
            do n_cheb=n_cheb_max+1,n_r_max
               work_Mloc(n_m,n_cheb)=zero
            end do
         end do

         !-- Finally assemble the right hand side
         do n_cheb=1,n_cheb_max
            do n_m=nMstart,nMstop
               dtemp_imp_Mloc_last(n_m,n_cheb)=TdiffFac*work_Mloc(n_m,n_cheb) 
            end do
         end do

         !-- Pad with zeros
         do n_cheb=n_cheb_max+1,n_r_max
            do n_m=nMstart,nMstop
               dtemp_imp_Mloc_last(n_m,n_cheb)=zero
            end do
         end do

      end if

   end subroutine get_temp_rhs_imp_int
!------------------------------------------------------------------------------
   subroutine get_lhs_mat(tscheme, A_mat, tempMat_fac, m)

      !-- Input variables
      type(type_tscheme), intent(in) :: tscheme    ! time step
      integer,            intent(in) :: m          ! Azimuthal wavenumber

      !-- Output variables
      type(type_bordmat_real), intent(inout) :: A_mat
      real(cp),                intent(inout) :: tempMat_fac(A_mat%nlines)

      !-- Local variables
      real(cp) :: stencilA4(A_mat%nbands)
      integer :: n_r, i_r, n_band, n_b
      real(cp) :: a, b

      !-- We have to fill A3 with zeros again otherwise on the next iteration
      !-- with a different dt there might be some issues with spurious values
      do n_r=1,A_mat%nlines_band
         do n_b=1,A_mat%ntau
            A_mat%A3(n_r,n_b)=0.0_cp
         end do
      end do

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      !-- Fill A4 banded-block
      do n_r=1,A_mat%nlines_band
         i_r = n_r+n_boundaries

         !-- Define the equations
         stencilA4 = intcheb2rmult2(a,b,i_r-1,A_mat%nbands)-     &
         &                       tscheme%wimp_lin(1)*TdiffFac*   &
         &     intcheb2rmult2lapl(a,b,m,i_r-1,A_mat%nbands)  

         !-- Roll the array for band storage
         do n_band=1,A_mat%nbands
            if ( n_r+A_mat%ku+1-n_band <= A_mat%nlines_band .and. n_r+A_mat%ku+1-n_band >= 1 ) then
               A_mat%A4(A_mat%kl+n_band,n_r+A_mat%ku+1-n_band) = rscheme%rnorm*stencilA4(n_band)
            end if
         end do

         tempMat_fac(n_r+A_mat%ntau)=one/rscheme%rnorm/maxval(abs(stencilA4))
      end do

      !-- Fill A3
      do n_r=1,A_mat%nlines_band
         i_r = n_r+A_mat%ntau
         stencilA4 = intcheb2rmult2(a,b,i_r-1,A_mat%nbands)-    &
         &                       tscheme%wimp_lin(1)*TdiffFac*  &
         &       intcheb2rmult2lapl(a,b,m,i_r-1,A_mat%nbands)  

         !-- Only the lower bands can contribute to the matrix A3
         do n_band=1,A_mat%kl
            if ( n_r <= n_band .and. n_r+A_mat%ntau-n_band >= 1 ) then
               A_mat%A3(n_r,n_r+A_mat%ntau-n_band) = rscheme%rnorm*stencilA4(A_mat%ku+1+n_band)
            end if
         end do
      end do

      !-- Add the tau Lines for boundary conditions into A1 and A2
      do n_r=1,A_mat%nlines
         if ( n_r <= A_mat%ntau ) then
            if ( ktopt == 1 ) then
               A_mat%A1(1,n_r)=rscheme%rnorm*rscheme%rMat(1,n_r)
            else
               A_mat%A1(1,n_r)=rscheme%rnorm*rscheme%drMat(1,n_r)
            end if
            if ( kbott == 1 ) then
               A_mat%A1(2,n_r)=rscheme%rnorm*rscheme%rMat(2,n_r)
            else
               A_mat%A1(2,n_r)=rscheme%rnorm*rscheme%drMat(2,n_r)
            end if
         else
            if ( ktopt == 1 ) then
               A_mat%A2(1,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%rMat(1,n_r)
            else
               A_mat%A2(1,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%drMat(1,n_r)
            end if
            if ( kbott == 1 ) then
               A_mat%A2(2,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%rMat(2,n_r)
            else
               A_mat%A2(2,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%drMat(2,n_r)
            end if
         end if
      end do

      !-- Cheb factor for boundary conditions
      do n_r=1,A_mat%ntau
         A_mat%A1(n_r,1)                =rscheme%boundary_fac*A_mat%A1(n_r,1)
         A_mat%A2(n_r,A_mat%nlines_band)=rscheme%boundary_fac*A_mat%A2(n_r,A_mat%nlines_band)
      end do

      !-- Continue to assemble precond matrix
      do n_r=1,A_mat%ntau
         tempMat_fac(n_r)=one/maxval(abs(A_mat%A2(n_r,:)))
      end do

      !-- Multiply the lines of the matrix by precond
      do n_r=1,A_mat%ntau
         A_mat%A1(n_r,:) = A_mat%A1(n_r,:)*tempMat_fac(n_r)
         A_mat%A2(n_r,:) = A_mat%A2(n_r,:)*tempMat_fac(n_r)
      end do

      do n_r=A_mat%ntau+1,A_mat%nlines
         A_mat%A3(n_r-A_mat%ntau,:) = A_mat%A3(n_r-A_mat%ntau,:)* tempMat_fac(n_r)
      end do

      do n_r=1,A_mat%nlines_band
         i_r = n_r+A_mat%ntau

         do n_band=1,A_mat%nbands
            if ( n_r+A_mat%ku+1-n_band <= A_mat%nlines_band .and. &
            &                         n_r+A_mat%ku+1-n_band >= 1 ) then
               A_mat%A4(A_mat%kl+n_band,n_r+A_mat%ku+1-n_band) =  &
               & A_mat%A4(A_mat%kl+n_band,n_r+A_mat%ku+1-n_band)* &
               &                      tempMat_fac(n_r+A_mat%ntau)
            end if
         end do

      end do

      !-- LU factorisation
      call A_mat%prepare_LU()

   end subroutine get_lhs_mat
!------------------------------------------------------------------------------
   subroutine get_rhs_exp_mat(B_mat)

      !-- Output variable
      type(type_bandmat_real), intent(inout) :: B_mat

      !-- Local variables
      real(cp) :: stencilB(B_mat%nbands)
      real(cp) :: a, b
      integer :: n_band, n_r, i_r

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      !-- Fill right-hand side matrix
      do n_r=1,B_mat%nlines
         i_r = n_r+n_boundaries

         !-- Define right-hand side equations
         stencilB = intcheb2rmult2(a,b,i_r-1,B_mat%nbands)

         !-- Roll array for band storage
         do n_band=1,B_mat%nbands
            if ( i_r+B_mat%ku+1-n_band <= B_mat%nlines .and. i_r+B_mat%ku+1-n_band >= 1 ) then
               B_mat%dat(n_band,i_r+B_mat%ku+1-n_band) = rscheme%rnorm*stencilB(n_band)
            end if
         end do
      end do

   end subroutine get_rhs_exp_mat
!------------------------------------------------------------------------------
   subroutine get_rhs_imp_mat(Cmat, m)

      !-- Input variable
      integer, intent(in) :: m

      !-- Output variable
      type(type_bandmat_real), intent(inout) :: Cmat

      !-- Local variables
      real(cp) :: stencilC(Cmat%nbands)
      real(cp) :: a, b
      integer :: n_band, n_r, i_r

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      !-- Fill right-hand side matrix
      do n_r=1,Cmat%nlines
         i_r = n_r+n_boundaries

         !-- Define right-hand side equations
         stencilC = intcheb2rmult2lapl(a,b,m,i_r-1,Cmat%nbands)

         !-- Roll array for band storage
         do n_band=1,Cmat%nbands
            if ( i_r+Cmat%ku+1-n_band <= Cmat%nlines .and. i_r+Cmat%ku+1-n_band >= 1 ) then
               Cmat%dat(n_band,i_r+Cmat%ku+1-n_band) = rscheme%rnorm*stencilC(n_band)
            end if
         end do
      end do

   end subroutine get_rhs_imp_mat
!------------------------------------------------------------------------------
end module update_temp_integ
