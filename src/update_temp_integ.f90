module update_temp_integ

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: zero, ci, half
   use pre_calculations, only: opr
   use namelists, only: kbott, ktopt, tadvz_fac, ra, r_cmb, r_icb, l_non_rot
   use radial_functions, only: rscheme, or1, r, dtcond, tcond, beta, &
       &                       rgrav
   use blocking, only: nMstart, nMstop
   use truncation, only: n_r_max, idx2m
   use radial_der, only: get_dr
   use fields, only: work_Mloc
   use useful, only: abortRun, roll
   use time_schemes, only: type_tscheme
   use matrix_types, only: type_bandmat_real, type_bordmat_real
   use chebsparselib, only: intcheb2rmult2hmult2, intcheb2rmult2hmult2lapl, &
       &                    intcheb2rmult1, intcheb2rmult2, intcheb2rmult2lapl

   implicit none
   
   private

   integer, parameter :: n_boundaries=2 ! Number of boundary conditions

   logical,  allocatable :: lTmat(:)
   complex(cp), allocatable :: rhs(:)

   type(type_bordmat_real), allocatable :: LHS_mat(:)
   type(type_bandmat_real) :: RHSE_mat, RHSI_mat
   type(type_bandmat_real), allocatable :: RHSIL_mat(:)

   public :: update_temp_int, initialize_temp_integ, finalize_temp_integ, &
   &         get_temp_rhs_imp_int

contains

   subroutine initialize_temp_integ

      integer :: n_m, m

      call RHSE_mat%initialize(3, 3, n_r_max)
      if ( l_non_rot ) then
         call RHSI_mat%initialize(4, 4, n_r_max)
      else
         call RHSI_mat%initialize(8, 8, n_r_max)
      end if

      allocate( RHSIL_mat(nMstart:nMstop) )
      allocate( LHS_mat(nMstart:nMstop) )

      if ( l_non_rot ) then
         do n_m=nMstart,nMstop
            call RHSIL_mat(n_m)%initialize(2, 2, n_r_max)
            call LHS_mat(n_m)%initialize(4, 4, n_boundaries, n_r_max)
         end do
      else
         do n_m=nMstart,nMstop
            call RHSIL_mat(n_m)%initialize(4, 4, n_r_max)
            call LHS_mat(n_m)%initialize(8, 8, n_boundaries, n_r_max)
         end do
      end if

      call get_rhs_exp_mat(RHSE_mat)
      call get_rhs_imp_mat(RHSI_mat)
      do n_m=nMstart,nMstop
         m = idx2m(n_m)
         call get_rhs_imp_lin_mat(RHSIL_mat(n_m), m)
      end do

      allocate( lTmat(nMstart:nMstop) )
      lTmat(:)=.false.
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_LOGICAL

      allocate( rhs(n_r_max) )
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_COMPLEX

   end subroutine initialize_temp_integ
!------------------------------------------------------------------------------
   subroutine finalize_temp_integ

      !-- Local variables
      integer :: n_m

      call RHSE_mat%finalize()
      call RHSI_mat%finalize()
      do n_m=nMstart,nMstop
         call RHSIL_mat(n_m)%finalize()
         call LHS_mat(n_m)%finalize()
      end do
      deallocate( LHS_mat, RHSIL_mat, rhs, lTmat )

   end subroutine finalize_temp_integ
!------------------------------------------------------------------------------
   subroutine update_temp_int(psi_Mloc, temp_Mloc, dtemp_Mloc, dVsT_Mloc,&
              &           dtemp_exp_Mloc, dtemp_imp_Mloc, buo_imp_Mloc,  &
              &           tscheme, lMat, l_roll_imp, l_log_next)

      !-- Input variables
      type(type_tscheme), intent(in) :: tscheme
      logical,            intent(in) :: lMat
      logical,            intent(in) :: l_roll_imp
      logical,            intent(in) :: l_log_next
      complex(cp),        intent(in) :: psi_Mloc(nMstart:nMstop, n_r_max)

      !-- Output variables
      complex(cp), intent(out) :: temp_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(out) :: dtemp_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(inout) :: dtemp_exp_Mloc(nMstart:nMstop,n_r_max,tscheme%norder_exp)
      complex(cp), intent(inout) :: dtemp_imp_Mloc(nMstart:nMstop,n_r_max,tscheme%norder_imp-1)
      complex(cp), intent(inout) :: buo_imp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: dVsT_Mloc(nMstart:nMstop, n_r_max)

      !-- Local variables
      integer :: n_r, n_m, n_r_out, m, n_o

      if ( lMat ) lTMat(:)=.false.

      !-- Assemble first buoyancy part from T^{n}
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m /= 0 ) then
               if ( l_non_rot ) then
                  buo_imp_Mloc(n_m,n_r)=-tscheme%wimp_lin(2)*rgrav(n_r)*r(n_r)**3* &
                  &                      ra*opr*ci*real(m,cp)*temp_Mloc(n_m,n_r)

               else
                  buo_imp_Mloc(n_m,n_r)=-tscheme%wimp_lin(2)*rgrav(n_r)*r(n_r)**3* &
                  &                      (r_cmb**2-r(n_r)**2)**3*ra*opr*           &
                  &                      ci*real(m,cp)*temp_Mloc(n_m,n_r)
               end if
            end if
         end do
      end do

      !-- Finish calculation of advection
      call get_dr( dVsT_Mloc, work_Mloc, nMstart, nMstop, n_r_max, &
           &       rscheme, nocopy=.true. )

      !-- Finish calculation of the explicit part for current time step
      do n_r=1,n_r_max
         do n_m=nMstart, nMstop
            m = idx2m(n_m)
            if ( l_non_rot ) then
               dtemp_exp_Mloc(n_m,n_r,1)=dtemp_exp_Mloc(n_m,n_r,1)   &
               &                             -work_Mloc(n_m,n_r)     &
               &    -ci*real(m,cp)*dtcond(n_r)*psi_Mloc(n_m,n_r)
            else
               dtemp_exp_Mloc(n_m,n_r,1)=dtemp_exp_Mloc(n_m,n_r,1)   &
               &         -(r_cmb**2-r(n_r)**2)*work_Mloc(n_m,n_r)    &
               &                 -ci*real(m,cp)*psi_Mloc(n_m,n_r)*(  &
               &               (r_cmb**2-r(n_r)**2)*dtcond(n_r)+     &
               &                    tadvz_fac*r(n_r)*tcond(n_r) )
            end if
         end do
      end do

      !-- Transform the explicit part to chebyshev space
      call rscheme%costf1(dtemp_exp_Mloc(:,:,1), nMstart, nMstop, n_r_max)

      !-- Matrix-vector multiplication by the operator \int\int r .
      do n_m=nMstart,nMstop
         do n_r=1,n_r_max
            rhs(n_r)=dtemp_exp_Mloc(n_m,n_r,1)
         end do
         call RHSE_mat%mat_vec_mul(rhs)
         rhs(1)=zero
         rhs(2)=zero
         do n_r=1,n_r_max
            dtemp_exp_Mloc(n_m,n_r,1)=rhs(n_r)
         end do
      end do

      !-- Calculation of the implicit part
      call get_temp_rhs_imp_int(temp_Mloc, tscheme%wimp_lin(2), &
           &                    dtemp_imp_Mloc(:,:,1))

      do n_m=nMstart, nMstop

         m = idx2m(n_m)
         
         if ( .not. lTmat(n_m) ) then
            call get_lhs_mat( tscheme, m, LHS_mat(n_m) )
            lTmat(n_m)=.true.
         end if

         do n_r=1,n_r_max
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

         call LHS_mat(n_m)%solve(rhs, n_r_max)

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

      !-- Assemble second buoyancy part from T^{n+1}
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m /= 0 ) then
               if ( l_non_rot ) then
                  buo_imp_Mloc(n_m,n_r)=            buo_imp_Mloc(n_m,n_r)-      &
                  &               tscheme%wimp_lin(1)*rgrav(n_r)*r(n_r)**3*     &
                  &               ra*opr*ci*real(m,cp)*temp_Mloc(n_m,n_r)
               else
                  buo_imp_Mloc(n_m,n_r)=            buo_imp_Mloc(n_m,n_r)-      &
                  &               tscheme%wimp_lin(1)*rgrav(n_r)*r(n_r)**3*     &
                  &               (r_cmb**2-r(n_r)**2)**3*ra*opr*ci*real(m,cp)* &
                  &                                    temp_Mloc(n_m,n_r)
               end if
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
   subroutine get_temp_rhs_imp_int(temp_Mloc, wimp, dtemp_imp_Mloc_last)

      !-- Input variables
      complex(cp), intent(in) :: temp_Mloc(nMstart:nMstop,n_r_max)
      real(cp),    intent(in) :: wimp

      !-- Output variable
      complex(cp), intent(out) :: dtemp_imp_Mloc_last(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_m

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            dtemp_imp_Mloc_last(n_m,n_r)=temp_Mloc(n_m,n_r)
         end do
      end do

      !-- Transform the implicit part to chebyshev space
      call rscheme%costf1(dtemp_imp_Mloc_last, nMstart, nMstop, n_r_max)

      !-- Matrix-vector multiplication by the operator \int\int r^2 .
      do n_m=nMstart,nMstop
         do n_r=1,n_r_max
            rhs(n_r)= dtemp_imp_Mloc_last(n_m,n_r)
         end do
         call RHSI_mat%mat_vec_mul(rhs)
         rhs(1)=zero
         rhs(2)=zero
         do n_r=1,n_r_max
            dtemp_imp_Mloc_last(n_m,n_r)=rhs(n_r)
         end do
      end do

      if ( wimp /= 0.0_cp ) then

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
            do n_r=1,n_r_max
               rhs(n_r)=work_Mloc(n_m,n_r)
            end do
            call RHSIL_mat(n_m)%mat_vec_mul(rhs)
            rhs(1)=zero
            rhs(2)=zero
            do n_r=1,n_r_max
               work_Mloc(n_m,n_r)=rhs(n_r)
            end do
         end do

         !-- Finally assemble the right hand side
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               dtemp_imp_Mloc_last(n_m,n_r)=dtemp_imp_Mloc_last(n_m,n_r) &
               &                            +wimp*opr*work_Mloc(n_m,n_r) 
            end do
         end do

      end if

   end subroutine get_temp_rhs_imp_int
!------------------------------------------------------------------------------
   subroutine get_lhs_mat(tscheme, m, A_mat)

      !-- Input variables
      type(type_tscheme), intent(in) :: tscheme    ! time step
      integer,            intent(in) :: m          ! Azimuthal wavenumber

      !-- Output variables
      type(type_bordmat_real), intent(inout) :: A_mat

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
         if ( l_non_rot ) then
            stencilA4 = intcheb2rmult2(a,b,i_r-1,A_mat%nbands)-    &
            &                       tscheme%wimp_lin(1)*opr*       &
            &     intcheb2rmult2lapl(a,b,m,i_r-1,A_mat%nbands)  
         else
            stencilA4 = intcheb2rmult2hmult2(a,b,r_cmb,i_r-1,A_mat%nbands)- &
            &                       tscheme%wimp_lin(1)*opr*                &
            &     intcheb2rmult2hmult2lapl(a,b,r_cmb,m,i_r-1,A_mat%nbands)  
         end if

         !-- Roll the array for band storage
         do n_band=1,A_mat%nbands
            if ( n_r+A_mat%ku+1-n_band <= A_mat%nlines_band .and. n_r+A_mat%ku+1-n_band >= 1 ) then
               A_mat%A4(A_mat%kl+n_band,n_r+A_mat%ku+1-n_band) = rscheme%rnorm*stencilA4(n_band)
            end if
         end do
      end do

      !-- Fill A3
      do n_r=1,A_mat%nlines_band
         i_r = n_r+A_mat%ntau
         if ( l_non_rot ) then
            stencilA4 = intcheb2rmult2(a,b,i_r-1,A_mat%nbands)-              &
            &                                   tscheme%wimp_lin(1)*opr*     &
            &       intcheb2rmult2lapl(a,b,m,i_r-1,A_mat%nbands)  
         else
            stencilA4 = intcheb2rmult2hmult2(a,b,r_cmb,i_r-1,A_mat%nbands)-    &
            &                                   tscheme%wimp_lin(1)*opr*       &
            &       intcheb2rmult2hmult2lapl(a,b,r_cmb,m,i_r-1,A_mat%nbands)  
         end if

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
         stencilB = intcheb2rmult1(a,b,i_r-1,B_mat%nbands)

         !-- Roll array for band storage
         do n_band=1,B_mat%nbands
            if ( i_r+B_mat%ku+1-n_band <= B_mat%nlines .and. i_r+B_mat%ku+1-n_band >= 1 ) then
               B_mat%dat(n_band,i_r+B_mat%ku+1-n_band) = rscheme%rnorm*stencilB(n_band)
            end if
         end do
      end do

   end subroutine get_rhs_exp_mat
!------------------------------------------------------------------------------
   subroutine get_rhs_imp_lin_mat(Cmat, m)

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
         if ( l_non_rot ) then
            stencilC = intcheb2rmult2lapl(a,b,m,i_r-1,Cmat%nbands)
         else
            stencilC = intcheb2rmult2hmult2lapl(a,b,r_cmb,m,i_r-1,Cmat%nbands)
         end if

         !-- Roll array for band storage
         do n_band=1,Cmat%nbands
            if ( i_r+Cmat%ku+1-n_band <= Cmat%nlines .and. i_r+Cmat%ku+1-n_band >= 1 ) then
               Cmat%dat(n_band,i_r+Cmat%ku+1-n_band) = rscheme%rnorm*stencilC(n_band)
            end if
         end do
      end do

   end subroutine get_rhs_imp_lin_mat
!------------------------------------------------------------------------------
   subroutine get_rhs_imp_mat(B_mat)

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
         if ( l_non_rot ) then
            stencilB = intcheb2rmult2(a,b,i_r-1,B_mat%nbands)
         else
            stencilB = intcheb2rmult2hmult2(a,b,r_cmb,i_r-1,B_mat%nbands)
         end if

         !-- Roll array for band storage
         do n_band=1,B_mat%nbands
            if ( i_r+B_mat%ku+1-n_band <= B_mat%nlines .and. i_r+B_mat%ku+1-n_band >= 1 ) then
               B_mat%dat(n_band,i_r+B_mat%ku+1-n_band) = rscheme%rnorm*stencilB(n_band)
            end if
         end do
      end do

   end subroutine get_rhs_imp_mat
!------------------------------------------------------------------------------
end module update_temp_integ
