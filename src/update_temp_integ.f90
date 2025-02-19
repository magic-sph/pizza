module update_temp_integ

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: zero, one, ci, half
   use namelists, only: kbott, ktopt, tadvz_fac, r_cmb, r_icb, &
       &                TdiffFac, l_non_rot, l_galerkin, l_full_disk
   use horizontal, only: hdif_T, bott_Mloc, topt_Mloc
   use radial_functions, only: rscheme, or2, dtcond, tcond, r
   use blocking, only: nMstart, nMstop
   use truncation, only: n_r_max, idx2m, n_cheb_max, m2idx
   use radial_der, only: get_dr
   use fields, only: work_Mloc
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use band_matrix, only: type_bandmat_real, band_band_product
   use bordered_matrix, only: type_bordmat_real
   use galerkin
   use chebsparselib, only: intcheb2rmult1, intcheb2rmult2lapl, intcheb2rmult2

   implicit none

   private

   integer, parameter :: n_boundaries=2 ! Number of BCs for this equation
   integer, parameter :: klA=4
   integer, parameter :: kuA=4

   logical,  allocatable :: lTmat(:), lAssmat(:)
   complex(cp), allocatable :: rhs(:)

   type(type_bordmat_real), allocatable :: LHS_mat_tau(:)
   type(type_bandmat_real), allocatable :: LHS_mat_gal(:)
   type(type_bordmat_real), allocatable :: Ass_mat_tau(:)
   type(type_bandmat_real), allocatable :: Ass_mat_gal(:)
   type(type_bandmat_real) :: RHSE_mat(2), gal_sten, gal_sten_m0
   type(type_bandmat_real), allocatable :: RHSI_mat(:)
   real(cp), allocatable :: tempfac(:,:) ! Precondition matrix
   real(cp), allocatable :: assfac(:,:)

   public :: update_temp_int, initialize_temp_integ, finalize_temp_integ, &
   &         get_temp_rhs_imp_int, finish_exp_temp_int, assemble_temp_int

contains

   subroutine initialize_temp_integ(tscheme)

      !-- Input variable
      class(type_tscheme), intent(in) :: tscheme

      !-- Local variables
      integer :: n_m, m

      if ( l_galerkin ) then
         !-- Define Galerkin basis and stencils
         if ( l_full_disk ) then
            if ( ktopt == 1 ) then
               call get_galerkin_stencil(gal_sten, n_r_max, 1)
               call get_galerkin_stencil(gal_sten_m0, n_r_max, 3)
            else
               call get_galerkin_stencil(gal_sten, n_r_max, 4)
               call get_galerkin_stencil(gal_sten_m0, n_r_max, 2)
            end if
         else
            if ( ktopt == 1 .and. kbott == 1 ) then
               call get_galerkin_stencil(gal_sten, n_r_max, 1)
            else if ( ktopt /= 1 .and. kbott /= 1 ) then
               call get_galerkin_stencil(gal_sten, n_r_max, 2)
            else if ( ktopt == 1 .and. kbott /= 1 ) then
               call get_galerkin_stencil(gal_sten, n_r_max, 3)
            else if ( ktopt /= 1 .and. kbott == 1 ) then
               call get_galerkin_stencil(gal_sten, n_r_max, 4)
            end if
         end if
      end if

      call RHSE_mat(1)%initialize(3, 3, n_r_max)
      call RHSE_mat(2)%initialize(4, 4, n_r_max)

      allocate( RHSI_mat(nMstart:nMstop) )
      do n_m=nMstart,nMstop
         call RHSI_mat(n_m)%initialize(2, 2, n_r_max)
      end do

      if ( l_galerkin ) then
         allocate( LHS_mat_gal(nMstart:nMstop) )
         if ( tscheme%l_assembly ) allocate( Ass_mat_gal(nMstart:nMstop) )
         !do n_m=nMstart,nMstop
         !   call LHS_mat_gal(n_m)%initialize(klA, kuA, n_r_max, l_lhs=.true.)
         !end do
      else
         allocate( LHS_mat_tau(nMstart:nMstop) )
         do n_m=nMstart,nMstop
            call LHS_mat_tau(n_m)%initialize(klA, kuA, n_boundaries, n_r_max)
         end do
         if ( tscheme%l_assembly ) then
            allocate( Ass_mat_tau(nMstart:nMstop) )
            do n_m=nMstart,nMstop
               call Ass_mat_tau(n_m)%initialize(klA, kuA, n_boundaries, n_r_max)
            end do
         end if
      end if

      call get_rhs_exp_mat(RHSE_mat(1),1)
      call get_rhs_exp_mat(RHSE_mat(2),2)
      do n_m=nMstart,nMstop
         m = idx2m(n_m)
         call get_rhs_imp_mat(RHSI_mat(n_m), m)
      end do

      allocate( lTmat(nMstart:nMstop) )
      lTmat(:)=.false.
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_LOGICAL

      if ( tscheme%l_assembly ) then
         allocate( lAssmat(nMstart:nMstop) )
         lAssmat(:)=.false.
         bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_LOGICAL
         allocate( assfac(n_r_max, nMstart:nMstop) )
         bytes_allocated = bytes_allocated + n_r_max*(nMstop-nMstart+1)* &
         &                 SIZEOF_DEF_REAL
      end if

      allocate( rhs(n_r_max) )
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_COMPLEX

      allocate( tempfac(n_r_max, nMstart:nMstop) )
      bytes_allocated = bytes_allocated + n_r_max*(nMstop-nMstart+1)* &
      &                 SIZEOF_DEF_REAL

   end subroutine initialize_temp_integ
!------------------------------------------------------------------------------
   subroutine finalize_temp_integ(tscheme)

      !-- Input variable
      class(type_tscheme), intent(in) :: tscheme

      !-- Local variables
      integer :: n_m

      call RHSE_mat(1)%finalize()
      call RHSE_mat(2)%finalize()
      do n_m=nMstart,nMstop
         call RHSI_mat(n_m)%finalize()
      end do
      deallocate( tempfac, RHSI_mat, rhs, lTmat )
      if ( tscheme%l_assembly ) deallocate(assfac, lAssmat)

      if ( l_galerkin ) then
         do n_m=nMstart,nMstop
            call LHS_mat_gal(n_m)%finalize()
         end do
         deallocate( LHS_mat_gal )
         call destroy_galerkin_stencil(gal_sten)
         if ( l_full_disk ) call destroy_galerkin_stencil(gal_sten_m0)
         if ( tscheme%l_assembly ) then
            do n_m=nMstart,nMstop
               call Ass_mat_gal(n_m)%finalize()
            end do
            deallocate( Ass_mat_gal )
         end if
      else
         do n_m=nMstart,nMstop
            call LHS_mat_tau(n_m)%finalize()
         end do
         deallocate( LHS_mat_tau )
         if ( tscheme%l_assembly ) then
            do n_m=nMstart,nMstop
               call Ass_mat_tau(n_m)%finalize()
            end do
            deallocate( Ass_mat_tau )
         end if
      end if

   end subroutine finalize_temp_integ
!------------------------------------------------------------------------------
   subroutine update_temp_int(temp_hat_Mloc, temp_Mloc, dtemp_Mloc, dTdt,  &
              &               tscheme, lMat, l_log_next)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lMat
      logical,             intent(in) :: l_log_next

      !-- Output variables
      complex(cp),       intent(out) :: temp_hat_Mloc(nMstart:nMstop, n_r_max)
      complex(cp),       intent(out) :: temp_Mloc(nMstart:nMstop, n_r_max)
      complex(cp),       intent(out) :: dtemp_Mloc(nMstart:nMstop, n_r_max)
      type(type_tarray), intent(inout) :: dTdt

      !-- Local variables
      integer :: n_m, m, n_cheb

      if ( lMat ) lTMat(:)=.false.

      !-- Now assemble the right hand side and store it in work_Mloc
      call tscheme%set_imex_rhs(work_Mloc, dTdt)

      do n_m=nMstart, nMstop

         m = idx2m(n_m)

         if ( .not. lTmat(n_m) ) then
            if ( l_galerkin ) then
               call get_lhs_mat_gal(tscheme%wimp_lin(1), LHS_mat_gal(n_m), &
                    &               tempfac(:,n_m), m)
            else
               call get_lhs_mat_tau(tscheme%wimp_lin(1), LHS_mat_tau(n_m), &
                    &               tempfac(:,n_m), m)
            end if
            lTmat(n_m)=.true.
         end if

         do n_cheb=1,n_r_max
            rhs(n_cheb)=work_Mloc(n_m,n_cheb)
         end do
         !-- Inhomogeneous heat B.Cs (if not zero or not galerkin)
         if ( .not. l_galerkin ) then
            rhs(1)=topt_Mloc(n_m)
            rhs(2)=bott_Mloc(n_m)
         endif

         !-- Multiply rhs by precond matrix
         do n_cheb=1,n_r_max
            rhs(n_cheb)=rhs(n_cheb)*tempfac(n_cheb,n_m)
         end do

         if ( l_galerkin ) then
            call LHS_mat_gal(n_m)%solve(rhs(1+n_boundaries:n_r_max), &
                 &                      n_r_max-n_boundaries)
            !-- Put the two first zero lines at the end
            rhs = cshift(rhs, n_boundaries)

            !-- Transform from Galerkin space to Chebyshev space
            if ( l_full_disk ) then
               if ( m == 0 ) then
                  call galerkin2cheb(gal_sten_m0, rhs)
               else
                  call galerkin2cheb(gal_sten, rhs)
               end if
            else
               call galerkin2cheb(gal_sten, rhs)
            end if
         else
            call LHS_mat_tau(n_m)%solve(rhs, n_r_max)
         end if

         do n_cheb=1,n_r_max
            temp_hat_Mloc(n_m, n_cheb)=rhs(n_cheb)
         end do

      end do

      !-- set cheb modes > n_cheb_max to zero (dealiazing)
      !do n_cheb=n_cheb_max+1,n_r_max
      !   do n_m=nMstart,nMstop
      !      temp_Mloc(n_m,n_cheb)=zero
      !   end do
      !end do

      !-- Copy temp_hat into temp_Mloc
      do n_cheb=1,n_r_max
         do n_m=nMstart,nMstop
            temp_Mloc(n_m,n_cheb)=temp_hat_Mloc(n_m,n_cheb)
         end do
      end do

      !-- Bring temperature back to physical space
      call rscheme%costf1(temp_Mloc, nMstart, nMstop, n_r_max)

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dTdt)

      !-- Compute implicit state
      if ( tscheme%istage == tscheme%nstages ) then
         call get_temp_rhs_imp_int(temp_hat_Mloc, dTdt, 1, tscheme%l_imp_calc_rhs(1))
      else
         call get_temp_rhs_imp_int(temp_hat_Mloc, dTdt, tscheme%istage+1,    &
              &                    tscheme%l_imp_calc_rhs(tscheme%istage+1))
      end if

      !-- In case log is needed on the next iteration, recalculate dT/dr
      !-- This is needed to estimate the heat fluxes
      if ( l_log_next ) then
         call get_dr(temp_Mloc, dtemp_Mloc, nMstart, nMstop, n_r_max, rscheme)
      end if

   end subroutine update_temp_int
!------------------------------------------------------------------------------
   subroutine assemble_temp_int(temp_hat_Mloc, temp_Mloc, dtemp_Mloc, dTdt,  &
              &                 tscheme, l_log_next)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: l_log_next

      !-- Output variables
      complex(cp),       intent(out) :: temp_hat_Mloc(nMstart:nMstop, n_r_max)
      complex(cp),       intent(out) :: temp_Mloc(nMstart:nMstop, n_r_max)
      complex(cp),       intent(out) :: dtemp_Mloc(nMstart:nMstop, n_r_max)
      type(type_tarray), intent(inout) :: dTdt

      !-- Local variables
      integer :: n_m, m, n_cheb

      !-- Now assemble the right hand side and store it in work_Mloc
      call tscheme%assemble_imex(work_Mloc, dTdt)

      do n_m=nMstart, nMstop

         m = idx2m(n_m)

         if ( .not. lAssmat(n_m) ) then
            if ( l_galerkin ) then
               call get_lhs_mat_gal(0.0_cp, Ass_mat_gal(n_m), assfac(:,n_m), m)
            else
               call get_lhs_mat_tau(0.0_cp, Ass_mat_tau(n_m), assfac(:,n_m), m)
            end if
            lAssmat(n_m)=.true.
         end if

         do n_cheb=1,n_r_max
            rhs(n_cheb)=work_Mloc(n_m,n_cheb)
         end do
         !-- Inhomogeneous heat B.Cs (if not zero or not galerkin)
         if ( .not. l_galerkin ) then
            rhs(1)=topt_Mloc(n_m)
            rhs(2)=bott_Mloc(n_m)
         endif

         !-- Multiply rhs by precond matrix
         do n_cheb=1,n_r_max
            rhs(n_cheb)=rhs(n_cheb)*assfac(n_cheb,n_m)
         end do

         if ( l_galerkin ) then
            call Ass_mat_gal(n_m)%solve(rhs(1+n_boundaries:n_r_max), &
                 &                      n_r_max-n_boundaries)
            !-- Put the two first zero lines at the end
            rhs = cshift(rhs, n_boundaries)

            !-- Transform from Galerkin space to Chebyshev space
            if ( l_full_disk ) then
               if ( m == 0 ) then
                  call galerkin2cheb(gal_sten_m0, rhs)
               else
                  call galerkin2cheb(gal_sten, rhs)
               end if
            else
               call galerkin2cheb(gal_sten, rhs)
            end if
         else
            call Ass_mat_tau(n_m)%solve(rhs, n_r_max)
         end if

         do n_cheb=1,n_r_max
            temp_hat_Mloc(n_m, n_cheb)=rhs(n_cheb)
         end do

      end do

      !-- Copy temp_hat into temp_Mloc
      do n_cheb=1,n_r_max
         do n_m=nMstart,nMstop
            temp_Mloc(n_m,n_cheb)=temp_hat_Mloc(n_m,n_cheb)
         end do
      end do

      !-- Bring temperature back to physical space
      call rscheme%costf1(temp_Mloc, nMstart, nMstop, n_r_max)

      !-- Compute implicit stage for first state of next iteration if required
      call get_temp_rhs_imp_int(temp_hat_Mloc, dTdt, 1, tscheme%l_imp_calc_rhs(1))

      !-- In case log is needed on the next iteration, recalculate dT/dr
      !-- This is needed to estimate the heat fluxes
      if ( l_log_next ) then
         call get_dr(temp_Mloc, dtemp_Mloc, nMstart, nMstop, n_r_max, rscheme)
      end if

   end subroutine assemble_temp_int
!------------------------------------------------------------------------------
   subroutine finish_exp_temp_int(psi_Mloc, dVsT_Mloc, dtemp_exp_last)

      !-- Input variables
      complex(cp), intent(in) :: psi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: dVsT_Mloc(nMstart:nMstop,n_r_max)

      !-- Output variables
      complex(cp), intent(inout) :: dtemp_exp_last(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_m, m, n_cheb
      real(cp) :: h2

      !-- Finish calculation of advection
      call rscheme%costf1(dVsT_Mloc, nMstart, nMstop, n_r_max)
      do n_cheb=n_cheb_max+1,n_r_max
         do n_m=nMstart,nMstop
            dVsT_Mloc(n_m,n_cheb)=zero
         end do
      end do
      call get_dr( dVsT_Mloc, work_Mloc, nMstart, nMstop, n_r_max, &
           &       rscheme, nocopy=.true., l_dct_in=.false. )

      !-- Finish calculation of the explicit part for current time step
      if ( l_non_rot ) then
         do n_r=1,n_r_max
            do n_m=nMstart, nMstop
               m = idx2m(n_m)
               dtemp_exp_last(n_m,n_r)= r(n_r)* dtemp_exp_last(n_m,n_r)   &
               &           -ci*real(m,cp)*dtcond(n_r)*psi_Mloc(n_m,n_r)
            end do
         end do
      else ! this is rotating
         do n_r=1,n_r_max
            h2 = r_cmb*r_cmb-r(n_r)*r(n_r)
            do n_m=nMstart, nMstop
               m = idx2m(n_m)
               dtemp_exp_last(n_m,n_r)=   r(n_r)*   dtemp_exp_last(n_m,n_r)   &
               &         -ci*real(m,cp)*(h2*dtcond(n_r)+                      &
               &           tadvz_fac*r(n_r)* tcond(n_r))* psi_Mloc(n_m,n_r)
            end do
         end do
      end if

      !-- Add the advection term
      do n_r=1,n_r_max
         dtemp_exp_last(:,n_r)=dtemp_exp_last(:,n_r)-work_Mloc(:,n_r)
      end do

      !-- Transform the explicit part to Chebyshev space
      call rscheme%costf1(dtemp_exp_last, nMstart, nMstop, n_r_max)

      do n_cheb=n_cheb_max+1,n_r_max
         do n_m=nMstart,nMstop
            dtemp_exp_last(n_m,n_cheb)=zero
         end do
      end do

      !-- Matrix-vector multiplication by the operator \int\int r^2 .
      do n_m=nMstart,nMstop
         do n_cheb=1,n_r_max
            rhs(n_cheb)=dtemp_exp_last(n_m,n_cheb)
         end do

         call RHSE_mat(1)%mat_vec_mul(rhs)

         rhs(1)=zero
         rhs(2)=zero
         do n_cheb=1,n_r_max
            dtemp_exp_last(n_m,n_cheb)=rhs(n_cheb)
         end do
      end do

   end subroutine finish_exp_temp_int
!------------------------------------------------------------------------------
   subroutine get_temp_rhs_imp_int(temp_hat_Mloc, dTdt, istage, l_calc_lin)

      !-- Input variables
      complex(cp),       intent(in) :: temp_hat_Mloc(nMstart:nMstop,n_r_max)
      logical,           intent(in) :: l_calc_lin
      integer,           intent(in) :: istage

      !-- Output variable
      type(type_tarray), intent(inout) :: dTdt

      !-- Local variables
      integer :: n_m, n_cheb

      !-- Matrix-vector multiplication by the operator \int\int r^2 .
      if ( istage == 1 ) then
         do n_m=nMstart,nMstop

            do n_cheb=1,n_r_max
               rhs(n_cheb)= temp_hat_Mloc(n_m,n_cheb)
            end do

            call RHSE_mat(2)%mat_vec_mul(rhs)

            rhs(1)=zero
            rhs(2)=zero

            do n_cheb=1,n_r_max
               dTdt%old(n_m,n_cheb,istage)=rhs(n_cheb)
            end do

         end do
      end if

      if ( l_calc_lin ) then

         !-- Matrix-vector multiplication by the LHS operator
         do n_m=nMstart,nMstop
            do n_cheb=1,n_r_max
               rhs(n_cheb)=temp_hat_Mloc(n_m,n_cheb)
            end do
            call RHSI_mat(n_m)%mat_vec_mul(rhs)
            rhs(1)=zero
            rhs(2)=zero
            do n_cheb=1,n_r_max
               work_Mloc(n_m,n_cheb)=rhs(n_cheb)
            end do
         end do

         !-- Finally assemble the right hand side
         do n_cheb=1,n_r_max
            do n_m=nMstart,nMstop
               dTdt%impl(n_m,n_cheb,istage)=TdiffFac*hdif_T(n_m)*&
               &                            work_Mloc(n_m,n_cheb)
            end do
         end do

      end if

   end subroutine get_temp_rhs_imp_int
!------------------------------------------------------------------------------
   subroutine get_lhs_mat_gal(wimp, Cmat, tempMat_fac, m)

      !-- Input variables
      real(cp), intent(in) :: wimp ! weight of time integrator
      integer,  intent(in) :: m    ! Azimuthal wavenumber

      !-- Output variables
      type(type_bandmat_real), intent(out) :: Cmat
      real(cp),                intent(inout) :: tempMat_fac(:)

      !-- Local variables
      type(type_bandmat_real) :: Amat
      real(cp), allocatable :: stencilA(:)
      integer :: n_r, i_r, n_band, n_m
      real(cp) :: a, b

      n_m = m2idx(m)
      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      call Amat%initialize(klA, kuA, n_r_max)
      allocate( stencilA(Amat%nbands) )

      !-- Fill A matrix
      do n_r=1,Amat%nlines
         i_r = n_r+n_boundaries

         !-- Define the equations
         stencilA = intcheb2rmult2(a,b,i_r-1,Amat%nbands)-      &
         &          wimp*TdiffFac*hdif_T(n_m)*                  &
         &          intcheb2rmult2lapl(a,b,m,i_r-1,Amat%nbands)

         !-- Roll the array for band storage
         do n_band=1,Amat%nbands
            if ( i_r+Amat%ku+1-n_band <= Amat%nlines .and. i_r+Amat%ku+1-n_band >= 1 ) then
               Amat%dat(n_band,i_r+Amat%ku+1-n_band) = rscheme%rnorm*stencilA(n_band)
            end if
         end do

      end do

      !-- Multiply by the Galerkin matrix and allocate Cmat
      if ( l_full_disk ) then
         if ( m == 0 ) then
            call band_band_product(Amat, gal_sten_m0, Cmat, l_lhs=.true.)
         else
            call band_band_product(Amat, gal_sten, Cmat, l_lhs=.true.)
         end if
      else
         call band_band_product(Amat, gal_sten, Cmat, l_lhs=.true.)
      end if

      !-- Remove first blank rows (careful remapping of kl, ku and room for LU factorisation)
      call Cmat%remove_leading_blank_rows(n_boundaries)

      !-- Truncate the last N lines
      call Cmat%remove_last_rows(n_boundaries)

      !-- Quick preconditionner
      deallocate( stencilA )
      allocate( stencilA(Cmat%nbands) )
      stencilA(:)=0.0_cp
      do n_r=1,Cmat%nlines
        do n_band=1,Cmat%nbands
           if ( n_r+Cmat%ku+1-n_band <= Cmat%nlines .and. &
           &                         n_r+Cmat%ku+1-n_band >= 1 ) then

               stencilA(n_band)=Cmat%dat(Cmat%kl+n_band,n_r+Cmat%ku+1-n_band)
           end if
         end do
         tempMat_fac(n_r+n_boundaries)=one/maxval(abs(stencilA))
      end do
      do n_r=1,Cmat%nlines
        do n_band=1,Cmat%nbands
           if ( n_r+Cmat%ku+1-n_band <= Cmat%nlines .and. &
           &                         n_r+Cmat%ku+1-n_band >= 1 ) then
              Cmat%dat(Cmat%kl+n_band,n_r+Cmat%ku+1-n_band) =  &
              & Cmat%dat(Cmat%kl+n_band,n_r+Cmat%ku+1-n_band)* &
              &                      tempMat_fac(n_r+n_boundaries)
          end if
       end do
      end do
      tempMat_fac(1:n_boundaries)=one

      !-- LU factorisation
      call Cmat%prepare_LU()

      !-- Remove temporary arrays
      deallocate( stencilA )
      call Amat%finalize()

   end subroutine get_lhs_mat_gal
!------------------------------------------------------------------------------
   subroutine get_lhs_mat_tau(wimp, A_mat, tempMat_fac, m)

      !-- Input variables
      real(cp), intent(in) :: wimp ! weight of time integrator
      integer,  intent(in) :: m          ! Azimuthal wavenumber

      !-- Output variables
      type(type_bordmat_real), intent(inout) :: A_mat
      real(cp),                intent(inout) :: tempMat_fac(A_mat%nlines)

      !-- Local variables
      real(cp) :: stencilA4(A_mat%nbands)
      integer :: n_r, i_r, n_band, n_b, n_m
      real(cp) :: a, b

      n_m = m2idx(m)

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
         &           wimp*TdiffFac*hdif_T(n_m)*                  &
         &           intcheb2rmult2lapl(a,b,m,i_r-1,A_mat%nbands)

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
         &           wimp*TdiffFac*hdif_T(n_m)*                 &
         &           intcheb2rmult2lapl(a,b,m,i_r-1,A_mat%nbands)

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
            if ( l_full_disk ) then
               if ( m == 0 ) then
                  A_mat%A1(2,n_r)=rscheme%rnorm*rscheme%drMat(2,n_r)
               else
                  A_mat%A1(2,n_r)=rscheme%rnorm*rscheme%rMat(2,n_r)
               end if
            else
               if ( kbott == 1 ) then
                  A_mat%A1(2,n_r)=rscheme%rnorm*rscheme%rMat(2,n_r)
               else
                  A_mat%A1(2,n_r)=rscheme%rnorm*rscheme%drMat(2,n_r)
               end if
            end if
         else
            if ( ktopt == 1 ) then
               A_mat%A2(1,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%rMat(1,n_r)
            else
               A_mat%A2(1,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%drMat(1,n_r)
            end if
            if ( l_full_disk ) then
               if ( m == 0 ) then
                  A_mat%A2(2,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%drMat(2,n_r)
               else
                  A_mat%A2(2,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%rMat(2,n_r)
               end if
            else
               if ( kbott == 1 ) then
                  A_mat%A2(2,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%rMat(2,n_r)
               else
                  A_mat%A2(2,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%drMat(2,n_r)
               end if
            end if
         end if
      end do

      !-- Cheb factor for boundary conditions
      do n_r=1,A_mat%ntau
         A_mat%A1(n_r,1)                =rscheme%boundary_fac*A_mat%A1(n_r,1)
         if ( n_r_max == A_mat%nlines_band ) then
            A_mat%A2(n_r,A_mat%nlines_band)=rscheme%boundary_fac*A_mat%A2(n_r,A_mat%nlines_band)
         end if
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

   end subroutine get_lhs_mat_tau
!------------------------------------------------------------------------------
   subroutine get_rhs_exp_mat(B_mat, i_type)

      !-- Input variable
      integer, intent(in) :: i_type

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
         if ( i_type == 1 ) then
            stencilB = intcheb2rmult1(a,b,i_r-1,B_mat%nbands)
         else if ( i_type == 2 ) then
            stencilB = intcheb2rmult2(a,b,i_r-1,B_mat%nbands)
         end if

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
