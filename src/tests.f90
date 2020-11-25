module tests
   !
   ! This module contains several testing subroutines
   !
   use iso_fortran_env, only: output_unit
   use precision_mod
   use parallel_mod
   use galerkin
   use constants, only: one, half, two, pi
   use namelists, only: l_newmap, radratio, l_rerror_fix, rerror_fac
   use chebyshev, only: type_cheb
   use radial_functions, only: rscheme
   use radial_der, only: initialize_der_arrays, finalize_der_arrays, get_ddr
   use radial_scheme, only: type_rscheme
   use algebra, only: prepare_full_mat, solve_full_mat
   use useful, only: abortRun
   use chebsparselib, only: intcheb2rmult1, rmult1, intcheb1, intcheb4, &
       &                    intcheb2, eye, intcheb4rmult4,  &
       &                    intcheb4rmult4laplrot2, intcheb4rmult4lapl2
   use integration, only: rInt_R
   use band_matrix, only: type_bandmat_real, band_band_product
   use bordered_matrix, only: type_bordmat_real

   implicit none

   private

   public :: solve_laplacian, solve_biharmo, test_radial_der, test_i4, &
   &         test_i4r4laplro

contains

   subroutine solve_laplacian(nMstart, nMstop)

      !-- Input variables
      integer, intent(in) :: nMstart
      integer, intent(in) :: nMstop

      !-- Local variables
      integer, allocatable :: nrs(:)
      real(cp), allocatable :: r_loc(:), sol(:), sol_theo(:)
      integer :: file_handle, n_r_max_loc, n_in, n_r
      real(cp) :: r_cmb, r_icb, eps, alph1, alph2, c1, c2
      real(cp) :: errColloc, errInteg, bbot, btop, epsh, errGal
      real(cp) :: tStart, tStop, tColloc, tInteg, tGal
      real(cp) :: timeLuColl, timeLuInt, timeSolveColl, timeSolveInt
      real(cp) :: timeLuGal, timeSolveGal

      eps = epsilon(1.0_cp)
      if ( l_newmap ) then
         n_in = 1
      else
         n_in = 0
      end if

      if ( rank == 0 ) then
         allocate ( type_cheb :: rscheme )
         r_cmb=one/(one-radratio)
         r_icb=r_cmb-one
         nrs = [10, 12, 16, 20, 24, 32, 36, 40, 48, 56, 64, 80, 96, 128,   &
         &      144, 160, 180, 200, 220, 256, 320, 400, 512, 640, 768, 1024, &
         &      1536,  2048, 3072, 4096]
         if ( l_newmap ) then
            open( newunit=file_handle, file='error_laplacian_map')
         else
            open( newunit=file_handle, file='error_laplacian')
         end if
         !btop = -0.25_cp
         !bbot =  0.75_cp
         btop = 0._cp
         bbot = 0._cp
         epsh =  3.0_cp

         do n_r=1,size(nrs)
            n_r_max_loc = nrs(n_r)
            allocate( r_loc(n_r_max_loc), sol(n_r_max_loc) )
            allocate( sol_theo(n_r_max_loc) )
            ! allocate( df_num(n_r_max_loc), d2f_num(n_r_max_loc) )
            ! allocate( df_theo(n_r_max_loc), d2f_theo(n_r_max_loc) )
            call rscheme%initialize(n_r_max_loc,n_r_max_loc, n_in, &
                 &                  l_cheb_coll=.true.)
            if ( l_newmap ) then
               alph1=1.0_cp/cosh(abs(log(eps))/(n_r_max_loc-1))
               alph2=0.0_cp
            else
               alph1=0.0_cp
               alph2=0.0_cp
            end if
            call rscheme%get_grid(n_r_max_loc, r_icb, r_cmb, alph1, alph2, r_loc)
            call initialize_der_arrays(n_r_max_loc, nMstart, nMstop, l_rerror_fix, &
                 &                     rerror_fac)
            call rscheme%get_der_mat(n_r_max_loc, l_cheb_coll=.true.)

            c1 = (btop-bbot-0.25_cp*epsh*(r_cmb**2-r_icb**2))/log(r_cmb/r_icb)
            c2 = btop-0.25_cp*epsh*r_cmb**2-c1*log(r_cmb)
            sol_theo(:)=0.25_cp*epsh*r_loc(:)**2+c1*log(r_loc)+c2

            tStart = MPI_Wtime()
            call solve_laplacian_colloc(n_r_max_loc, r_loc, rscheme, sol, &
                 &                      timeLuColl, timeSolveColl, bbot,  &
                 &                      btop, epsh)
            tStop = MPI_Wtime()
            tColloc = tStop-tStart
            errColloc = maxval(abs(sol(:)-sol_theo(:)))

            tStart = MPI_Wtime()
            call solve_laplacian_integ_tau(n_r_max_loc, r_loc, rscheme, sol,    & 
                 &                         timeLuInt, timeSolveInt, bbot, btop, &
                 &                         epsh)
            tStop = MPI_Wtime()
            tInteg = tStop-tStart
            errInteg = maxval(abs(sol(:)-sol_theo(:)))

            tStart = MPI_Wtime()
            call solve_laplacian_integ_galerkin(n_r_max_loc, r_loc, rscheme, sol, &
                &                              timeLuGal, timeSolveGal, epsh)
            tStop = MPI_Wtime()
            tGal = tStop-tStart
            errGal = maxval(abs(sol(:)-sol_theo(:)))

            !-- Store it in a file
            write(file_handle, '(i5,12es14.6)') n_r_max_loc, errColloc,  &
            &                                   tColloc, timeLuColl,     &
            &                                   timeSolveColl, errInteg, &
            &                                   tInteg, timeLuInt,       &
            &                                   timeSolveInt, errGal,    &
            &                                   tGal, timeLuGal, timeSolveGal
            write(output_unit, '(i5,6es14.6)') n_r_max_loc, errColloc, tColloc,    &
            &                                  errInteg, tInteg, errGal, tGal


            call finalize_der_arrays()
            call rscheme%finalize()
            deallocate( r_loc, sol, sol_theo )
         end do

         close(file_handle)

      end if

   end subroutine solve_laplacian
!------------------------------------------------------------------------------
   subroutine test_i4r4laplro(nMstart, nMstop)

      integer, intent(in) :: nMstart
      integer, intent(in) :: nMstop

      !-- Local variables
      integer, allocatable :: nrs(:)
      real(cp), allocatable :: r_loc(:), rhs1(:), rhs2(:), or2(:), rhs3(:)
      integer :: file_handle, n_r_max_loc, n_in, n_r
      real(cp) :: r_cmb, r_icb, eps, alph1, alph2, err
      real(cp) :: ro
      type(type_bordmat_real) :: A_mat
      type(type_bandmat_real) :: B_mat

      ro = 4.15_cp

      eps = epsilon(1.0_cp)
      if ( l_newmap ) then
         n_in = 1
      else
         n_in = 0
      end if

      if ( rank == 0 ) then
         allocate ( type_cheb :: rscheme )
         r_cmb=one/(one-radratio)
         r_icb=radratio/(one-radratio)
         nrs = [16, 20, 24, 32, 36, 40, 48, 56, 64, 80, 96, 128,            &
         &      144, 160, 180, 200, 220, 256, 320, 400, 512, 640, 768, 1024,&
         &      2048, 4096]
         open( newunit=file_handle, file='sol', form='unformatted')

         do n_r=1,size(nrs)
            n_r_max_loc = nrs(n_r)
            allocate( r_loc(n_r_max_loc), rhs1(n_r_max_loc), rhs2(n_r_max_loc) )
            allocate( or2(n_r_max_loc), rhs3(n_r_max_loc) )

            call rscheme%initialize(n_r_max_loc,n_r_max_loc, n_in, &
                 &                  l_cheb_coll=.true.)
            alph1=0.0_cp
            alph2=0.0_cp
            call rscheme%get_grid(n_r_max_loc, r_icb, r_cmb, alph1, alph2, r_loc)
            or2(:)=one/r_loc(:)/r_loc(:)
            call initialize_der_arrays(n_r_max_loc, nMstart, nMstop, l_rerror_fix, &
                 &                     rerror_fac)
            call rscheme%get_der_mat(n_r_max_loc, l_cheb_coll=.true.)

            call A_mat%initialize(12, 12, 4, n_r_max_loc)
            ! call A_mat%initialize(4, 4, 4, n_r_max_loc)
            call B_mat%initialize(1, 1, n_r_max_loc)

            call get_lhs_mat(A_mat, m=5, n_r_max=n_r_max_loc, &
                 &           r_icb=r_icb, r_cmb=r_cmb)
            call get_rhs_mat(B_mat, r_icb, r_cmb)

            !-- Define a RHS
            rhs1(:) = cos(r_loc(:))

            !-- Bring it to Chebyshev space
            call rscheme%costf1(rhs1, n_r_max_loc)
            call B_mat%mat_vec_mul(rhs1)

            rhs1(1)=0.0_cp
            rhs1(2)=0.0_cp
            rhs1(3)=0.0_cp
            rhs1(4)=0.0_cp

            call A_mat%solve(rhs1, n_r_max_loc)

            call rscheme%costf1(rhs1, n_r_max_loc)
            write(file_handle) nrs(n_r)
            write(file_handle) r_loc
            write(file_handle) rhs1

            err = sum(abs(rhs1))/n_r_max_loc
            write(output_unit, '(i5,5es20.12)') n_r_max_loc, err, rhs1(1), &
            &                                   rhs1(n_r_max_loc/2)

            call A_mat%finalize()
            call B_mat%finalize()
            call finalize_der_arrays()

            call rscheme%finalize()
            deallocate( r_loc, rhs1, rhs2, rhs3, or2 )
         end do

         close(file_handle)

      end if

   end subroutine test_i4r4laplro
!------------------------------------------------------------------------------
   subroutine solve_biharmo(nMstart, nMstop)

      !-- Input variables
      integer, intent(in) :: nMstart
      integer, intent(in) :: nMstop

      !-- Local variables
      integer, allocatable :: nrs(:)
      real(cp), allocatable :: r_loc(:), sol(:), sol_theo(:)
      integer :: file_handle, n_r_max_loc, n_in, n_r
      real(cp) :: r_cmb, r_icb, eps, alph1, alph2
      real(cp) :: errColloc, errInteg, errGal
      real(cp) :: tStart, tStop, tColloc, tInteg, tGal
      real(cp) :: timeLuColl, timeLuInt, timeSolveColl, timeSolveInt
      real(cp) :: timeLuGal, timeSolveGal

      eps = epsilon(1.0_cp)
      if ( l_newmap ) then
         n_in = 1
      else
         n_in = 0
      end if

      if ( rank == 0 ) then
         allocate ( type_cheb :: rscheme )
         r_cmb=two*pi
         r_icb=0.0_cp
         nrs = [12, 16, 20, 24, 32, 36, 40, 48, 56, 64, 80, 96, 128,        &
         &      144, 160, 180, 200, 220, 256, 320, 400, 512, 640, 768, 1024,&
         &      2048, 4096]
         if ( l_newmap ) then
            open( newunit=file_handle, file='error_biharmo_map')
         else
            open( newunit=file_handle, file='error_biharmo')
         end if

         do n_r=1,size(nrs)
            n_r_max_loc = nrs(n_r)
            allocate( r_loc(n_r_max_loc), sol(n_r_max_loc) )
            allocate( sol_theo(n_r_max_loc) )
            ! allocate( df_num(n_r_max_loc), d2f_num(n_r_max_loc) )
            ! allocate( df_theo(n_r_max_loc), d2f_theo(n_r_max_loc) )
            call rscheme%initialize(n_r_max_loc,n_r_max_loc, n_in, &
                 &                  l_cheb_coll=.true.)
            if ( l_newmap ) then
               alph1=1.0_cp/cosh(abs(log(eps))/(n_r_max_loc-1))
               alph2=0.0_cp
            else
               alph1=0.0_cp
               alph2=0.0_cp
            end if
            call rscheme%get_grid(n_r_max_loc, r_icb, r_cmb, alph1, alph2, r_loc)
            call initialize_der_arrays(n_r_max_loc, nMstart, nMstop, l_rerror_fix, &
                 &                     rerror_fac)
            call rscheme%get_der_mat(n_r_max_loc, l_cheb_coll=.true.)

            sol_theo(:)=one/8.0_cp*(two*(r_loc(:)-pi)*sin(r_loc)+ &
            &           (two*pi-r_loc(:))*r_loc(:)*cos(r_loc))

            tStart = MPI_Wtime()
            call solve_biharmo_colloc(n_r_max_loc, r_loc, rscheme, sol, &
                 &                      timeLuColl, timeSolveColl)
            tStop = MPI_Wtime()
            tColloc = tStop-tStart
            errColloc = maxval(abs(sol(:)-sol_theo(:)))

            tStart = MPI_Wtime()
            call solve_biharmo_integ(n_r_max_loc, r_loc, rscheme, sol, &
                 &                  timeLuInt, timeSolveInt)
            tStop = MPI_Wtime()
            tInteg = tStop-tStart
            errInteg = maxval(abs(sol(:)-sol_theo(:)))

            tStart = MPI_Wtime()
            call solve_biharmo_galerkin(n_r_max_loc, r_loc, rscheme, sol, &
                 &                      timeLuGal, timeSolveGal)
            tStop = MPI_Wtime()
            tGal = tStop-tStart
            errGal = maxval(abs(sol(:)-sol_theo(:)))

            !-- Store it in a file
            write(file_handle, '(i5,12es14.6)') n_r_max_loc, errColloc,  &
            &                                   tColloc, timeLuColl,     &
            &                                   timeSolveColl, errInteg, &
            &                                   tInteg, timeLuInt,       &
            &                                   timeSolveInt, errGal,    &
            &                                   tGal, timeLuGal, timeSolveGal
            write(output_unit, '(i5,6es14.6)') n_r_max_loc, errColloc, tColloc,  &
            &                                  errInteg, tInteg, errGal, tGal


            call finalize_der_arrays()
            call rscheme%finalize()
            deallocate( r_loc, sol, sol_theo )
            ! deallocate( r_loc, f, df_theo, d2f_theo, df_num, d2f_num )
         end do

         close(file_handle)

      end if

   end subroutine solve_biharmo
!------------------------------------------------------------------------------
   subroutine solve_laplacian_integ_galerkin(n_r_max, r, rscheme, rhs, timeLu, &
              &                              timeSolve, epsh)

      !-- Input variables
      integer,             intent(in) :: n_r_max
      real(cp),            intent(in) :: r(n_r_max)
      class(type_rscheme), intent(in) :: rscheme
      real(cp),            intent(in) :: epsh

      !-- Output variable
      real(cp),            intent(out) :: rhs(n_r_max)
      real(cp),            intent(out) :: timeLu
      real(cp),            intent(out) :: timeSolve

      !-- Local variables
      type(type_bandmat_real) :: gal_sten, Amat, Cmat, Bmat
      real(cp), allocatable :: stencilA(:), stencilB(:)
      real(cp) :: a, b, tStart, tStop
      integer :: i_r, klA, kuA, klB, kuB
      integer :: n_band, n_r, n_boundaries, lenA

      klA = 1
      kuA = 1
      klB = 3
      kuB = 3
      n_boundaries = 2
      lenA = n_r_max-n_boundaries

      call Amat%initialize(klA, kuA, n_r_max)
      call Bmat%initialize(klB, kuB, n_r_max)

      allocate ( stencilA(klA+kuA+1), stencilB(klB+kuB+1) )

      call get_galerkin_stencil(gal_sten, n_r_max, 1)

      a = half*(r(1)-r(n_r_max))
      b = half*(r(1)+r(n_r_max))

      !-- Fill right-hand side matrix
      do n_r=1,lenA
         i_r = n_r+n_boundaries

         !-- Define the equations 
         stencilA = rmult1(a,b,i_r-1,klA+kuA+1)-intcheb1(a,i_r-1,klA+kuA+1)
         stencilB = intcheb2rmult1(a,b,i_r-1,klB+kuB+1)

         !-- Roll the arrays for band storage
         do n_band=1,klB+kuB+1
            if ( i_r+kuB+1-n_band <= n_r_max .and. i_r+kuB+1-n_band >= 1 ) then
               Bmat%dat(n_band,i_r+kuB+1-n_band) = rscheme%rnorm*stencilB(n_band)
            end if
         end do
         !-- Roll the arrays for band storage
         do n_band=1,klA+kuA+1
            if ( i_r+kuA+1-n_band <= n_r_max .and. i_r+kuA+1-n_band >= 1 ) then
               Amat%dat(n_band,i_r+kuA+1-n_band) = rscheme%rnorm*stencilA(n_band)
            end if
         end do
      end do

      !-- Multiply by the Galerkin matrix
      call band_band_product(Amat, gal_sten, Cmat, l_lhs=.true.)

      !-- Assemble right-hand side
      rhs(:)= epsh
      call rscheme%costf1(rhs, n_r_max)
      call Bmat%mat_vec_mul(rhs)

      !-- Remove first blank rows (careful remapping of kl, ku and room for LU factorisation)
      call Cmat%remove_leading_blank_rows(n_boundaries)

      !-- Truncate the last N lines
      call Cmat%remove_last_rows(n_boundaries)

      !-- LU factorisation
      tStart = MPI_Wtime()
      call Cmat%prepare_LU()
      tStop = MPI_Wtime()
      timeLu= tStop-tStart

      !-- Solve
      tStart = MPI_Wtime()
      call Cmat%solve(rhs(1+n_boundaries:n_r_max),n_r_max-n_boundaries)
      tStop = MPI_Wtime()
      timeSolve= tStop-tStart

      !-- Put the two first zeroes at the end
      rhs = cshift(rhs, n_boundaries)

      !-- Transform from Galerkin space to Chebyshev space
      call galerkin2cheb(gal_sten, rhs)

      !-- Final DCT to bring the solution back to physical space
      call rscheme%costf1(rhs, n_r_max)

      call destroy_galerkin_stencil(gal_sten)
      call Amat%finalize()
      call Cmat%finalize()
      call Bmat%finalize()
      deallocate ( stencilA, stencilB )

   end subroutine solve_laplacian_integ_galerkin
!------------------------------------------------------------------------------
   subroutine solve_laplacian_integ_tau(n_r_max, r, rscheme, rhs, timeLu, &
              &                         timeSolve, bbot, btop, epsh)

      !-- Input variables
      integer,             intent(in) :: n_r_max
      real(cp),            intent(in) :: r(n_r_max)
      class(type_rscheme), intent(in) :: rscheme
      real(cp),            intent(in) :: bbot
      real(cp),            intent(in) :: btop
      real(cp),            intent(in) :: epsh

      !-- Output variable
      real(cp),            intent(out) :: rhs(n_r_max)
      real(cp),            intent(out) :: timeLu
      real(cp),            intent(out) :: timeSolve

      !-- Local variables
      type(type_bandmat_real) :: Bmat
      type(type_bordmat_real) :: Amat
      real(cp), allocatable :: stencilA(:), stencilB(:)
      real(cp) :: a, b, tStart, tStop
      integer :: i_r, klA, kuA, klB, kuB
      integer :: n_band, n_r, n_boundaries, lenA

      klA = 1
      kuA = 1
      klB  = 3
      kuB  = 3
      n_boundaries = 2
      lenA = n_r_max-n_boundaries

      call Amat%initialize(klA, kuA, n_boundaries, n_r_max)
      call Bmat%initialize(klB, kuB, n_r_max)

      allocate ( stencilA(klA+kuA+1), stencilB(klB+kuB+1) )
      
      a = half*(r(1)-r(n_r_max))
      b = half*(r(1)+r(n_r_max))

      !-- Fill A4 banded-block
      do n_r=1,lenA
         i_r = n_r+n_boundaries
         
         !-- Define the equations
         stencilA = rmult1(a,b,i_r-1,klA+kuA+1)-intcheb1(a,i_r-1,klA+kuA+1)

         !-- Roll the array for band storage
         do n_band=1,klA+kuA+1
            if ( n_r+kuA+1-n_band <= lenA .and. n_r+kuA+1-n_band >= 1 ) then
               Amat%A4(klA+n_band,n_r+kuA+1-n_band) = rscheme%rnorm*stencilA(n_band)
            end if
         end do
      end do

      !-- Fill A3
      do n_r=1,lenA
         i_r = n_r+n_boundaries
         stencilA = rmult1(a,b,i_r-1,klA+kuA+1)-intcheb1(a,i_r-1,klA+kuA+1)

         !-- Only the lower bands can contribute to the matrix A3
         do n_band=1,klA
            if ( n_r <= n_band .and. n_r+n_boundaries-n_band >= 1 ) then
               Amat%A3(n_r,n_r+n_boundaries-n_band) = rscheme%rnorm*stencilA(kuA+1+n_band)
            end if
         end do
      end do

      !-- Fill right-hand side matrix
      do n_r=1,lenA
         i_r = n_r+n_boundaries

         !-- Define the equations 
         stencilB = intcheb2rmult1(a,b,i_r-1,klB+kuB+1)

         !-- Roll the arrays for band storage
         do n_band=1,klB+kuB+1
            if ( i_r+kuB+1-n_band <= n_r_max .and. i_r+kuB+1-n_band >= 1 ) then
               Bmat%dat(n_band,i_r+kuB+1-n_band) = rscheme%rnorm*stencilB(n_band)
            end if
         end do
      end do

      !-- Add the tau Lines for boundary conditions into A1 and A2
      do n_r=1,n_r_max
         if ( n_r <= n_boundaries) then
            Amat%A1(1,n_r)=rscheme%rnorm*rscheme%rMat(1,n_r)
            Amat%A1(2,n_r)=rscheme%rnorm*rscheme%rMat(n_r_max,n_r)
         else
            Amat%A2(1,n_r-n_boundaries)=rscheme%rnorm*rscheme%rMat(1,n_r)
            Amat%A2(2,n_r-n_boundaries)=rscheme%rnorm*rscheme%rMat(n_r_max,n_r)
         end if
      end do

      !-- Cheb factor for boundary conditions
      do n_r=1,n_boundaries
         Amat%A1(n_r,1)               =rscheme%boundary_fac*Amat%A1(n_r,1)
         Amat%A2(n_r,Amat%nlines_band)=rscheme%boundary_fac* &
         &                                  Amat%A2(n_r,Amat%nlines_band)
      end do

      !-- Assemble right-hand side
      rhs(:)=epsh
      call rscheme%costf1(rhs, n_r_max)
      call Bmat%mat_vec_mul(rhs)
      !-- Boundary conditions
      rhs(1) = btop
      rhs(2) = bbot

      !-- LU factorisation
      tStart = MPI_Wtime()
      call Amat%prepare_LU()
      tStop = MPI_Wtime()
      timeLu= tStop-tStart

      !-- Solve
      tStart = MPI_Wtime()
      call Amat%solve(rhs, n_r_max)
      tStop = MPI_Wtime()
      timeSolve= tStop-tStart

      !-- Final DCT to bring the solution back to physical space
      call rscheme%costf1(rhs, n_r_max)

      call Amat%finalize()
      call Bmat%finalize()
      deallocate ( stencilA, stencilB )

   end subroutine solve_laplacian_integ_tau
!------------------------------------------------------------------------------
   subroutine solve_laplacian_colloc(n_r_max, r, rscheme, rhs, timeLu, &
              &                      timeSolve, bbot, btop, epsh)

      !-- Input variables
      integer,             intent(in) :: n_r_max
      real(cp),            intent(in) :: r(n_r_max)
      class(type_rscheme), intent(in) :: rscheme
      real(cp),            intent(in) :: bbot
      real(cp),            intent(in) :: btop
      real(cp),            intent(in) :: epsh

      !-- Output variable
      real(cp),            intent(out) :: rhs(n_r_max)
      real(cp),            intent(out) :: timeLu
      real(cp),            intent(out) :: timeSolve

      !-- Local variables
      real(cp), allocatable :: mat(:,:)
      real(cp) :: tStart, tStop
      integer, allocatable :: pivot(:)
      integer :: nR_out, nR, info

      allocate( mat(n_r_max, n_r_max), pivot(n_r_max) )

      !-- Boundary conditions =  fixed values
      do nR_out=1,rscheme%n_max
         mat(1,nR_out)=rscheme%rnorm*rscheme%rMat(1,nR_out)
         mat(n_r_max,nR_out)=rscheme%rnorm*rscheme%rMat(n_r_max,nR_out)
      end do
      rhs(:)       = epsh
      rhs(1)       = btop
      rhs(n_r_max) = bbot

      !----- Other points:
      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            mat(nR,nR_out)= rscheme%rnorm * ( rscheme%d2rMat(nR,nR_out) + &
            &             1.0_cp/r(nR)*        rscheme%drMat(nR,nR_out) )
         end do
      end do

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         mat(nR,1)      =rscheme%boundary_fac*mat(nR,1)
         mat(nR,n_r_max)=rscheme%boundary_fac*mat(nR,n_r_max)
      end do

      !----- LU decomposition:
      tStart = MPI_Wtime()
      call prepare_full_mat(mat,n_r_max,n_r_max,pivot,info)
      tStop = MPI_Wtime()
      timeLu = tStop-tStart
      if ( info /= 0 ) then
         call abortRun('Singular matrix mat!')
      end if

      tStart = MPI_Wtime()
      call solve_full_mat(mat, n_r_max, n_r_max, pivot, rhs)
      tStop = MPI_Wtime()
      timeSolve = tStop-tStart
      call rscheme%costf1(rhs, n_r_max)

      deallocate( mat, pivot )

   end subroutine solve_laplacian_colloc
!------------------------------------------------------------------------------
   subroutine solve_biharmo_galerkin(n_r_max, r, rscheme, rhs, timeLu, timeSolve)

      !-- Input variables
      integer,             intent(in) :: n_r_max
      real(cp),            intent(in) :: r(n_r_max)
      class(type_rscheme), intent(in) :: rscheme

      !-- Output variable
      real(cp),            intent(out) :: rhs(n_r_max)
      real(cp),            intent(out) :: timeLu
      real(cp),            intent(out) :: timeSolve

      !-- Local variables
      type(type_bandmat_real) :: Amat, Bmat, Cmat, gal_sten
      real(cp), allocatable :: stencilB(:), stencilA(:)
      real(cp) :: a, b, tStart, tStop
      integer :: i_r, klA, kuA, klB, kuB, n_band, n_r
      integer :: n_boundaries, lenA

      klA = 4
      kuA = 4
      klB  = 4
      kuB  = 4
      n_boundaries = 4
      lenA = n_r_max-n_boundaries

      call Amat%initialize(klA, kuA, n_r_max)
      call Bmat%initialize(klB, kuB, n_r_max)
      allocate ( stencilB(klB+kuB+1), stencilA(klA+kuA+1) )

      call get_galerkin_stencil(gal_sten, n_r_max, 2)
      
      a = half*(r(1)-r(n_r_max))
      b = half*(r(1)+r(n_r_max))

      !-- Fill matrices
      do n_r=1,lenA
         i_r = n_r+n_boundaries

         !-- Define right-hand side equations
         stencilB = intcheb4(a,i_r-1,klB+kuB+1)

         stencilA = eye(klA+kuA+1) + two*intcheb2(a,i_r-1,klA+kuA+1)+ &
         &           intcheb4(a,i_r-1,klA+kuA+1)

         !-- Roll array for band storage
         do n_band=1,klB+kuB+1
            if ( i_r+kuB+1-n_band <= n_r_max .and. i_r+kuB+1-n_band >= 1 ) then
               Bmat%dat(n_band,i_r+kuB+1-n_band) = rscheme%rnorm*stencilB(n_band)
            end if
         end do

         !-- Roll the array for band storage
         do n_band=1,klA+kuA+1
            if ( i_r+kuA+1-n_band <= n_r_max .and. i_r+kuA+1-n_band >= 1 ) then
               Amat%dat(n_band,i_r+kuA+1-n_band) = rscheme%rnorm*stencilA(n_band)
            end if
         end do
      end do

      !-- Multiply by the Galerkin matrix
      call band_band_product(Amat, gal_sten, Cmat, l_lhs=.true.)

      !-- Assemble right-hand side
      rhs(:) = cos(r)
      call rscheme%costf1(rhs, n_r_max)
      call Bmat%mat_vec_mul(rhs)

      !-- Remove first blank rows (careful remapping of kl, ku and room for LU factorisation)
      call Cmat%remove_leading_blank_rows(n_boundaries)

      !-- Truncate the last N lines
      call Cmat%remove_last_rows(n_boundaries)

      !-- LU factorisation of A matrix
      tStart = MPI_Wtime()
      call Cmat%prepare_LU()
      tStop = MPI_Wtime()
      timeLu= tStop-tStart

      !-- Solve
      tStart = MPI_Wtime()
      call Cmat%solve(rhs(1+n_boundaries:n_r_max),n_r_max-n_boundaries)
      tStop = MPI_Wtime()
      timeSolve= tStop-tStart

      !-- Put the two first zeroes at the end
      rhs = cshift(rhs, n_boundaries)

      !-- Transform from Galerkin space to Chebyshev space
      call galerkin2cheb(gal_sten, rhs)

      !-- Final DCT to bring back the data to physical space
      call rscheme%costf1(rhs, n_r_max)

      call destroy_galerkin_stencil(gal_sten)
      call Amat%finalize()
      call Bmat%finalize()
      call Cmat%finalize()
      deallocate ( stencilA, stencilB )

   end subroutine solve_biharmo_galerkin
!------------------------------------------------------------------------------
   subroutine solve_biharmo_integ(n_r_max, r, rscheme, rhs, timeLu, timeSolve)

      !-- Input variables
      integer,             intent(in) :: n_r_max
      real(cp),            intent(in) :: r(n_r_max)
      class(type_rscheme), intent(in) :: rscheme

      !-- Output variable
      real(cp),            intent(out) :: rhs(n_r_max)
      real(cp),            intent(out) :: timeLu
      real(cp),            intent(out) :: timeSolve

      !-- Local variables
      type(type_bordmat_real) :: Amat
      type(type_bandmat_real) :: Bmat
      real(cp), allocatable :: stencilB(:), stencilA(:)
      real(cp) :: a, b, tStart, tStop
      integer :: i_r, klA, kuA, klB, kuB, n_band, n_r
      integer :: n_boundaries, lenA

      klA = 4
      kuA = 4
      klB = 4
      kuB = 4
      n_boundaries = 4
      lenA = n_r_max-n_boundaries

      call Amat%initialize(klA,kuA,n_boundaries,n_r_max)
      call Bmat%initialize(klB,kuB,n_r_max)
      allocate ( stencilB(klB+kuB+1), stencilA(klA+kuA+1) )
      
      a = half*(r(1)-r(n_r_max))
      b = half*(r(1)+r(n_r_max))

      !-- Fill A4 banded-block
      do n_r=1,lenA
         i_r = n_r+n_boundaries

         !-- Define the equations
         stencilA = eye(klA+kuA+1) + two*intcheb2(a,i_r-1,klA+kuA+1)+ &
         &          intcheb4(a,i_r-1,klA+kuA+1)

         !-- Roll the array for band storage
         do n_band=1,klA+kuA+1
            if ( n_r+kuA+1-n_band <= lenA .and. n_r+kuA+1-n_band >= 1 ) then
               Amat%A4(klA+n_band,n_r+kuA+1-n_band) = rscheme%rnorm*stencilA(n_band)
            end if
         end do
      end do

      !-- Fill A3
      do n_r=1,lenA
         i_r = n_r+n_boundaries
         stencilA = eye(klA+kuA+1) + two*intcheb2(a,i_r-1,klA+kuA+1)+ &
         &          intcheb4(a,i_r-1,klA+kuA+1)

         !-- Only the lower bands can contribute to the matrix A3
         do n_band=1,klA
            if ( n_r <= n_band .and. n_r+n_boundaries-n_band >= 1 ) then
               Amat%A3(n_r,n_r+n_boundaries-n_band) = rscheme%rnorm*stencilA(kuA+1+n_band)
            end if
         end do
      end do

      !-- Fill right-hand side matrix
      do n_r=1,lenA
         i_r = n_r+n_boundaries

         !-- Define right-hand side equations
         stencilB = intcheb4(a,i_r-1,klB+kuB+1)

         !-- Roll array for band storage
         do n_band=1,klB+kuB+1
            if ( i_r+kuB+1-n_band <= n_r_max .and. i_r+kuB+1-n_band >= 1 ) then
               Bmat%dat(n_band,i_r+kuB+1-n_band) = rscheme%rnorm*stencilB(n_band)
            end if
         end do
      end do

      !-- Add the tau Lines for boundary conditions into A1 and A2
      do n_r=1,n_r_max
         if ( n_r <= n_boundaries) then
            Amat%A1(1,n_r)=rscheme%rnorm*rscheme%rMat(1,n_r)
            Amat%A1(2,n_r)=rscheme%rnorm*rscheme%rMat(n_r_max,n_r)
            Amat%A1(3,n_r)=rscheme%rnorm*rscheme%drMat(1,n_r)
            Amat%A1(4,n_r)=rscheme%rnorm*rscheme%drMat(n_r_max,n_r)
         else
            Amat%A2(1,n_r-n_boundaries)=rscheme%rnorm*rscheme%rMat(1,n_r)
            Amat%A2(2,n_r-n_boundaries)=rscheme%rnorm*rscheme%rMat(n_r_max,n_r)
            Amat%A2(3,n_r-n_boundaries)=rscheme%rnorm*rscheme%drMat(1,n_r)
            Amat%A2(4,n_r-n_boundaries)=rscheme%rnorm*rscheme%drMat(n_r_max,n_r)
         end if
      end do

      !-- Cheb factors for boundary conditions
      do n_r=1,n_boundaries
         Amat%A1(n_r,1)    =rscheme%boundary_fac*Amat%A1(n_r,1)
         Amat%A2(n_r,lenA)=rscheme%boundary_fac*Amat%A2(n_r,lenA)
      end do

      !-- Assemble right-hand side
      rhs(:) = cos(r)
      call rscheme%costf1(rhs, n_r_max)
      call Bmat%mat_vec_mul(rhs)
      !-- Boundary conditions
      rhs(1) = 0.0_cp
      rhs(2) = 0.0_cp
      rhs(3) = 0.0_cp
      rhs(4) = 0.0_cp

      !-- LU factorisation of A matrix
      tStart = MPI_Wtime()
      call Amat%prepare_LU()
      tStop = MPI_Wtime()
      timeLu= tStop-tStart

      !-- Solve
      tStart = MPI_Wtime()
      call Amat%solve(rhs, n_r_max)
      tStop = MPI_Wtime()
      timeSolve= tStop-tStart

      !-- Final DCT to bring back the data to physical space
      call rscheme%costf1(rhs, n_r_max)

      call Amat%finalize()
      call Bmat%finalize()
      deallocate ( stencilA, stencilB )

   end subroutine solve_biharmo_integ
!------------------------------------------------------------------------------
   subroutine solve_biharmo_colloc(n_r_max, r, rscheme, rhs, timeLu, &
              &                    timeSolve)

      !-- Input variables
      integer,             intent(in) :: n_r_max
      real(cp),            intent(in) :: r(n_r_max)
      class(type_rscheme), intent(in) :: rscheme

      !-- Output variable
      real(cp),            intent(out) :: rhs(n_r_max)
      real(cp),            intent(out) :: timeLu
      real(cp),            intent(out) :: timeSolve

      !-- Local variables
      real(cp), allocatable :: mat(:,:), tmp(:)
      real(cp) :: tStart, tStop
      integer, allocatable :: pivot(:)
      integer :: nR_out, nR, info, nR_out_psi, nR_psi

      allocate( mat(2*n_r_max, 2*n_r_max), tmp(2*n_r_max), pivot(2*n_r_max) )
      mat(:,:)=0.0_cp
      tmp(:)=0.0_cp

      !-- Boundary conditions =  fixed values
      do nR_out=1,rscheme%n_max
         mat(1,nR_out)      =rscheme%rnorm*rscheme%rMat(1,nR_out)
         mat(n_r_max,nR_out)=rscheme%rnorm*rscheme%rMat(n_r_max,nR_out)
         mat(n_r_max+1,nR_out)=rscheme%rnorm*rscheme%drMat(1,nR_out)
         mat(2*n_r_max,nR_out)=rscheme%rnorm*rscheme%drMat(n_r_max,nR_out)
      end do
      tmp(1:n_r_max) = cos(r)
      tmp(1)         = 0.0_cp
      tmp(n_r_max)   = 0.0_cp
      tmp(n_r_max+1) = 0.0_cp
      tmp(2*n_r_max) = 0.0_cp

      !----- Other points:
      do nR_out=1,n_r_max
         nR_out_psi=nR_out+n_r_max
         do nR=2,n_r_max-1
            nR_psi=nR+n_r_max
            mat(nR,nR_out)= rscheme%rnorm * rscheme%rMat(nR,nR_out) 
            mat(nR,nR_out_psi)= rscheme%rnorm *( rscheme%d2rMat(nR,nR_out)+&
            &                                  two*rscheme%rMat(nR,nR_out) )

            mat(nR_psi,nR_out) = rscheme%rnorm * rscheme%d2rMat(nR,nR_out)
            mat(nR_psi,nR_out_psi) = -rscheme%rnorm * rscheme%rMat(nR,nR_out)
         end do
      end do

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         nR_psi = nR+n_r_max
         mat(nR,1)            =rscheme%boundary_fac*mat(nR,1)
         mat(nR,n_r_max)      =rscheme%boundary_fac*mat(nR,n_r_max)
         mat(nR,n_r_max+1)    =rscheme%boundary_fac*mat(nR,n_r_max+1)
         mat(nR,2*n_r_max)    =rscheme%boundary_fac*mat(nR,2*n_r_max)
         mat(nR_psi,1)        =rscheme%boundary_fac*mat(nR_psi,1)
         mat(nR_psi,n_r_max)  =rscheme%boundary_fac*mat(nR_psi,n_r_max)
         mat(nR_psi,n_r_max+1)=rscheme%boundary_fac*mat(nR_psi,n_r_max+1)
         mat(nR_psi,2*n_r_max)=rscheme%boundary_fac*mat(nR_psi,2*n_r_max)
      end do


      !----- LU decomposition:
      tStart = MPI_Wtime()
      call prepare_full_mat(mat,2*n_r_max,2*n_r_max,pivot,info)
      tStop = MPI_Wtime()
      timeLu = tStop-tStart
      if ( info /= 0 ) then
         call abortRun('Singular matrix mat!')
      end if

      tStart = MPI_Wtime()
      call solve_full_mat(mat, 2*n_r_max, 2*n_r_max, pivot, tmp)
      tStop = MPI_Wtime()
      timeSolve = tStop-tStart
      rhs(:)=tmp(1:n_r_max)
      call rscheme%costf1(rhs, n_r_max)

      deallocate( mat, pivot )

   end subroutine solve_biharmo_colloc
!------------------------------------------------------------------------------
   subroutine get_rhs_mat(D_mat, r_icb, r_cmb)
      !
      ! This corresponds to the matrix that goes in front of the explicit terms
      !

      real(cp),           intent(in) :: r_icb
      real(cp),           intent(in) :: r_cmb

      !-- Output variable
      type(type_bandmat_real), intent(inout) :: D_mat

      !-- Local variables
      real(cp) :: stencilD(D_mat%nbands)
      real(cp) :: a, b
      integer :: n_band, n_r, i_r, n_bounds

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      n_bounds = 4

      !-- Fill right-hand side matrix
      do n_r=1,D_mat%nlines
         i_r = n_r+n_bounds

         !-- Define right-hand side equations
         stencilD = eye(D_mat%nbands)

         !-- Roll array for band storage
         do n_band=1,D_mat%nbands
            if ( i_r+D_mat%ku+1-n_band <= D_mat%nlines .and. i_r+D_mat%ku+1-n_band >= 1 ) then
               D_mat%dat(n_band,i_r+D_mat%ku+1-n_band) = rscheme%rnorm*stencilD(n_band)
            end if
         end do
      end do

   end subroutine get_rhs_mat
!------------------------------------------------------------------------------
   subroutine get_lhs_mat(A_mat, m, n_r_max, r_icb, r_cmb)

      !-- Input variables
      integer,            intent(in) :: m          ! Azimuthal wavenumber
      integer,            intent(in) :: n_r_max
      real(cp),           intent(in) :: r_icb
      real(cp),           intent(in) :: r_cmb

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
         i_r = n_r+A_mat%ntau

         !-- Define the equations
         stencilA4 = intcheb4rmult4laplrot2(a,b,m,i_r-1,A_mat%nbands)  
         ! stencilA4 = intcheb4rmult4lapl2(a,b,m,i_r-1,A_mat%nbands)  

         print*, n_r, maxval(abs(stencilA4)), minval(abs(stencilA4))

         !-- Roll the array for band storage
         do n_band=1,A_mat%nbands
            if ( n_r+A_mat%ku+1-n_band <= A_mat%nlines_band .and. n_r+A_mat%ku+1-n_band >= 1 ) then
               A_mat%A4(A_mat%kl+n_band,n_r+A_mat%ku+1-n_band) = rscheme%rnorm* &
               &             stencilA4(n_band)
            end if
         end do
      end do

      !print*, 'min/max', minval(abs(A_mat%A4)), maxval(abs(A_mat%A4))

      !-- Fill A3
      do n_r=1,A_mat%nlines_band
         i_r = n_r+A_mat%ntau

         stencilA4 = intcheb4rmult4laplrot2(a,b,m,i_r-1,A_mat%nbands)  
         ! stencilA4 = intcheb4rmult4lapl2(a,b,m,i_r-1,A_mat%nbands)  

         !-- Only the lower bands can contribute to the matrix A3
         do n_band=1,A_mat%kl
            if ( n_r <= n_band .and. n_r+A_mat%ntau-n_band >= 1 ) then
               A_mat%A3(n_r,n_r+A_mat%ntau-n_band) = rscheme%rnorm*    &
               &      stencilA4(A_mat%ku+1+n_band)
            end if
         end do
      end do

      !-- Add the tau Lines for boundary conditions into A1 and A2
      do n_r=1,A_mat%nlines
         if ( n_r <= A_mat%ntau ) then
            A_mat%A1(1,n_r)=rscheme%rnorm*rscheme%rMat(1,n_r)
            A_mat%A1(2,n_r)=rscheme%rnorm*rscheme%rMat(n_r_max,n_r)
            A_mat%A1(3,n_r)=rscheme%rnorm*rscheme%drMat(1,n_r)
            A_mat%A1(4,n_r)=rscheme%rnorm*rscheme%drMat(n_r_max,n_r)
         else
            A_mat%A2(1,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%rMat(1,n_r)
            A_mat%A2(2,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%rMat(n_r_max,n_r)
            A_mat%A2(3,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%drMat(1,n_r)
            A_mat%A2(4,n_r-A_mat%ntau)=rscheme%rnorm*rscheme%drMat(n_r_max,n_r)
         end if
      end do

      !-- Cheb factor for boundary conditions
      do n_r=1,A_mat%ntau
         A_mat%A1(n_r,1)                =rscheme%boundary_fac*A_mat%A1(n_r,1)
         A_mat%A2(n_r,A_mat%nlines_band)=rscheme%boundary_fac*A_mat%A2(n_r,A_mat%nlines_band)
         print*, n_r, maxval(abs(A_mat%A2(n_r, :))), minval(abs(A_mat%A2(n_r,:)))
      end do

      ! do n_r=1,A_mat%ntau
         ! print*, A_mat%A1(n_r,:)
      ! end do
      stop
      !-- LU factorisation
      call A_mat%prepare_LU()

   end subroutine get_lhs_mat
!------------------------------------------------------------------------------
   subroutine test_i4()
      !
      ! This subroutine is used to compare the implementation of 4th order
      ! recurrence relation: i4r4() vs i4(). This is a simple dgemv comparison
      ! of I4*f with I4R4*f/r^4
      !

      !-- Local variables
      integer, allocatable :: nrs(:)
      real(cp), allocatable :: r_loc(:), rhs1(:), rhs2(:), or2(:), rhs3(:)
      integer :: file_handle, n_r_max_loc, n_in, n_r
      real(cp) :: r_cmb, r_icb, eps, alph1, alph2, err
      real(cp) :: ro
      type(type_bandmat_real) :: I4_mat, I4R4_mat, I4R4H6_mat

      ro = 4.15_cp

      eps = epsilon(1.0_cp)
      if ( l_newmap ) then
         n_in = 1
      else
         n_in = 0
      end if

      if ( rank == 0 ) then
         allocate ( type_cheb :: rscheme )
         r_cmb=one/(one-radratio)
         r_icb=radratio/(one-radratio)
         nrs = [16, 20, 24, 32, 36, 40, 48, 56, 64, 80, 96, 128,            &
         &      144, 160, 180, 200, 220, 256, 320, 400, 512, 640, 768, 1024,&
         &      2048, 4096, 8192, 16384]
         open( newunit=file_handle, file='error_i4')

         do n_r=1,size(nrs)
            n_r_max_loc = nrs(n_r)
            allocate( r_loc(n_r_max_loc), rhs1(n_r_max_loc), rhs2(n_r_max_loc) )
            allocate( or2(n_r_max_loc), rhs3(n_r_max_loc) )

            call rscheme%initialize(n_r_max_loc,n_r_max_loc, n_in, &
                 &                  l_cheb_coll=.true.)
            alph1=0.0_cp
            alph2=0.0_cp
            call rscheme%get_grid(n_r_max_loc, r_icb, r_cmb, alph1, alph2, r_loc)
            or2(:)=one/r_loc(:)/r_loc(:)

            call I4_mat%initialize(4, 4, n_r_max_loc)
            call I4R4_mat%initialize(8, 8, n_r_max_loc)
            call I4R4H6_mat%initialize(14, 14, n_r_max_loc)

            call fill_i4_mat(I4_mat, r_cmb, r_icb, 4)
            call fill_i4r4_mat(I4R4_mat, r_cmb, r_icb, 4)
            call fill_i4r4h6_mat(I4R4H6_mat, r_cmb, r_icb, 4)

            !-- Define a RHS
            rhs1(:) = sqrt(r_loc)
            rhs2(:) = sqrt(r_loc)*or2(:)*or2(:)
            rhs3(:) = sqrt(r_loc)/(ro**2-r_loc(:)**2)**3

            !-- Bring it to Chebyshev space
            call rscheme%costf1(rhs1, n_r_max_loc)
            call rscheme%costf1(rhs2, n_r_max_loc)
            call rscheme%costf1(rhs3, n_r_max_loc)

            !-- DGEMM
            !call I4_mat%mat_vec_mul(rhs1)
            !call I4R4_mat%mat_vec_mul(rhs2)
            call I4R4_mat%mat_vec_mul(rhs1)
            call I4R4H6_mat%mat_vec_mul(rhs3)
            rhs1(1)=0.0_cp
            rhs1(2)=0.0_cp
            rhs1(3)=0.0_cp
            rhs1(4)=0.0_cp

            rhs3(1)=0.0_cp
            rhs3(2)=0.0_cp
            rhs3(3)=0.0_cp
            rhs3(4)=0.0_cp

            err = maxval(abs(rhs1-rhs3))
            write(file_handle, '(i5,5es20.12)') n_r_max_loc, err, rhs1(5), &
            &                                   rhs3(5), rhs1(9), rhs3(9)
            write(output_unit, '(i5,5es20.12)') n_r_max_loc, err, rhs1(5), rhs3(5), &
            &                                   rhs1(9), rhs3(9)



            call I4_mat%finalize()
            call I4R4_mat%finalize()
            call I4R4H6_mat%finalize()

            call rscheme%finalize()
            deallocate( r_loc, rhs1, rhs2, rhs3, or2 )
         end do

         close(file_handle)

      end if

   end subroutine test_i4
!------------------------------------------------------------------------------
   subroutine fill_i4_mat(D_mat, r_cmb, r_icb, n_boundaries)

      !-- Input variables
      real(cp), intent(in) :: r_cmb
      real(cp), intent(in) :: r_icb
      integer,  intent(in) :: n_boundaries

      !-- Output variable
      type(type_bandmat_real), intent(inout) :: D_mat

      !-- Local variables
      real(cp) :: stencilD(D_mat%nbands)
      real(cp) :: a, b
      integer :: n_band, n_r, i_r

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      !-- Fill right-hand side matrix
      do n_r=1,D_mat%nlines
         i_r = n_r+n_boundaries

         !-- Define right-hand side equations
         stencilD = intcheb4(a,i_r-1,D_mat%nbands)

         !-- Roll array for band storage
         do n_band=1,D_mat%nbands
            if ( i_r+D_mat%ku+1-n_band <= D_mat%nlines .and. i_r+D_mat%ku+1-n_band >= 1 ) then
               D_mat%dat(n_band,i_r+D_mat%ku+1-n_band) = rscheme%rnorm*stencilD(n_band)
            end if
         end do
      end do

   end subroutine fill_i4_mat
!------------------------------------------------------------------------------
   subroutine fill_i4r4_mat(D_mat, r_cmb, r_icb, n_boundaries)

      !-- Input variables
      real(cp), intent(in) :: r_cmb
      real(cp), intent(in) :: r_icb
      integer,  intent(in) :: n_boundaries

      !-- Output variable
      type(type_bandmat_real), intent(inout) :: D_mat

      !-- Local variables
      real(cp) :: stencilD(D_mat%nbands)
      real(cp) :: a, b
      integer :: n_band, n_r, i_r

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      !-- Fill right-hand side matrix
      do n_r=1,D_mat%nlines
         i_r = n_r+n_boundaries

         !-- Define right-hand side equations
         stencilD = intcheb4rmult4(a,b,i_r-1,D_mat%nbands)

         !-- Roll array for band storage
         do n_band=1,D_mat%nbands
            if ( i_r+D_mat%ku+1-n_band <= D_mat%nlines .and. i_r+D_mat%ku+1-n_band >= 1 ) then
               D_mat%dat(n_band,i_r+D_mat%ku+1-n_band) = rscheme%rnorm*stencilD(n_band)
            end if
         end do
      end do

   end subroutine fill_i4r4_mat
!------------------------------------------------------------------------------
   subroutine fill_i4r4h6_mat(D_mat, r_cmb, r_icb, n_boundaries)

      !-- Input variables
      real(cp), intent(in) :: r_cmb
      real(cp), intent(in) :: r_icb
      integer,  intent(in) :: n_boundaries

      !-- Output variable
      type(type_bandmat_real), intent(inout) :: D_mat

      !-- Local variables
      real(cp) :: stencilD(D_mat%nbands)
      real(cp) :: a, b
      integer :: n_band, n_r, i_r

      a = half*(r_cmb-r_icb)
      b = half*(r_cmb+r_icb)

      !-- Fill right-hand side matrix
      do n_r=1,D_mat%nlines
         i_r = n_r+n_boundaries

         !-- Define right-hand side equations
         stencilD = intcheb4rmult4(a,b,i_r-1,D_mat%nbands)

         !-- Roll array for band storage
         do n_band=1,D_mat%nbands
            if ( i_r+D_mat%ku+1-n_band <= D_mat%nlines .and. i_r+D_mat%ku+1-n_band >= 1 ) then
               D_mat%dat(n_band,i_r+D_mat%ku+1-n_band) = rscheme%rnorm*stencilD(n_band)
            end if
         end do
      end do

   end subroutine fill_i4r4h6_mat
!------------------------------------------------------------------------------
   subroutine test_radial_der(nMstart, nMstop)

      !-- Input variables
      integer, intent(in) :: nMstart
      integer, intent(in) :: nMstop

      !-- Local variables
      integer, allocatable :: nrs(:)
      integer :: n_r, n_r_max_loc, n_in
      real(cp), allocatable :: r_loc(:), f(:)
      real(cp), allocatable :: df_num(:), d2f_num(:)
      real(cp), allocatable :: df_num_matmul(:)
      real(cp), allocatable :: df_theo(:), d2f_theo(:)
      real(cp), allocatable :: dcheb(:,:), d2cheb(:,:)
      real(cp) :: r_cmb, r_icb, err_d1, err_d2
      real(cp) :: err_d1_matmul, err_d2_matmul
      real(cp) :: alph1, alph2
      integer :: file_handle
      real(cp) :: eps

      eps = epsilon(1.0_cp)

      if ( l_newmap ) then
         n_in = 1
      else
         n_in = 0
      end if

      if ( rank == 0 ) then
         allocate ( type_cheb :: rscheme )
         r_cmb=one/(one-radratio)
         r_icb=r_cmb-one
         nrs = [6, 8, 12, 16, 20, 24, 32, 36, 40, 48, 56, 64, 80, 96, 128,   &
         &      144, 160, 180, 200, 220, 256, 320, 400, 512, 640, 768, 1024, &
         &      1536,  2048]
         if ( l_newmap ) then
            open( newunit=file_handle, file='error_drMat_map')
         else
            open( newunit=file_handle, file='error_drMat')
         end if

         do n_r=1,size(nrs)
            n_r_max_loc = nrs(n_r)
            allocate( r_loc(n_r_max_loc), f(n_r_max_loc) )
            allocate( df_num(n_r_max_loc), d2f_num(n_r_max_loc) )
            allocate( dcheb(n_r_max_loc,n_r_max_loc) )
            allocate( d2cheb(n_r_max_loc,n_r_max_loc) )
            allocate( df_num_matmul(n_r_max_loc) )
            allocate( df_theo(n_r_max_loc), d2f_theo(n_r_max_loc) )
            call rscheme%initialize(n_r_max_loc, n_r_max_loc, n_in, &
                 &                  l_cheb_coll=.true.)
            if ( l_newmap ) then
               alph1=one/cosh(abs(log(eps))/(n_r_max_loc-1))
               alph2=0.0_cp
            else
               alph1=0.0_cp
               alph2=0.0_cp
            end if
            call rscheme%get_grid(n_r_max_loc, r_icb, r_cmb, alph1, alph2, r_loc)
            call initialize_der_arrays(n_r_max_loc, nMstart, nMstop, l_rerror_fix, &
                 &                     rerror_fac)
            call rscheme%get_der_mat(n_r_max_loc, l_cheb_coll=.true.)

            !-- Define a function and its analytical derivatives
            f(:) = sqrt(r_loc(:))
            df_theo(:) = 0.5_cp/sqrt(r_loc(:))
            d2f_theo(:) = -0.25_cp/r_loc(:)**(1.5_cp)

            !amp = 50.0_cp
            !f(:) = tanh(amp*(r_loc(:)-r_icb))*tanh(amp*(r_cmb-r_loc(:)))
            !df_theo(:)= amp*((1.0_cp-(tanh(amp*(r_loc(:)-r_icb)))**2.0_cp)*   &
            !&           tanh(amp*(r_cmb-r_loc(:)))-tanh(amp*(r_loc(:)-r_icb))*&
            !&           (1.0_cp-(tanh(amp*(r_cmb-r_loc(:))))**2.0_cp))
            !d2f_theo(:)=amp*amp*(-2.0_cp*tanh(amp*(r_loc(:)-r_icb))*    &
            !&           (1.0_cp-(tanh(amp*(r_loc(:)-r_icb)))**2.0_cp)*  &
            !&           tanh(amp*(r_cmb-r_loc(:)))-                     &
            !&           2.0_cp*tanh(amp*(r_loc(:)-r_icb))*              &
            !&           tanh(amp*(r_cmb-r_loc(:)))*                     &
            !&           (1.0_cp-(tanh(amp*(r_cmb-r_loc(:))))**2.0_cp))

            !-- Get the derivatives numerically
            call get_ddr(f, df_num, d2f_num, n_r_max_loc, rscheme)
            err_d1 = maxval(abs(df_num(:)-df_theo(:)))
            err_d2 = maxval(abs(d2f_num(:)-d2f_theo(:)))

            dcheb(:,:) = rscheme%rnorm*rscheme%drMat(:,:)
            dcheb(:,1) = rscheme%boundary_fac*dcheb(:,1)
            dcheb(:,n_r_max_loc) = rscheme%boundary_fac*dcheb(:,n_r_max_loc)
            d2cheb(:,:) = rscheme%rnorm*rscheme%d2rMat(:,:)
            d2cheb(:,1) = rscheme%boundary_fac*d2cheb(:,1)
            d2cheb(:,n_r_max_loc) = rscheme%boundary_fac*d2cheb(:,n_r_max_loc)
            call rscheme%costf1(f,n_r_max_loc)
            call dgemv('N', n_r_max_loc, n_r_max_loc, 1.0_cp, dcheb,  &
                 &      n_r_max_loc, f, 1, 0.0_cp, df_num_matmul,&
                 &      1)
            err_d1_matmul = maxval(abs(df_num_matmul(:)-df_theo(:)))
            call dgemv('N', n_r_max_loc, n_r_max_loc, 1.0_cp, d2cheb,  &
                 &      n_r_max_loc, f, 1, 0.0_cp, df_num_matmul,&
                 &      1)
            err_d2_matmul = maxval(abs(df_num_matmul(:)-d2f_theo(:)))

            !-- Store it in a file
            write(file_handle, '(i5,4es20.12)') n_r_max_loc, err_d1, err_d1_matmul, &
            &                                   err_d2, err_d2_matmul
            write(output_unit, '(i5,4es20.12)') n_r_max_loc, err_d1, err_d1_matmul, &
            &                                   err_d2, err_d2_matmul

            call finalize_der_arrays()
            call rscheme%finalize()
            deallocate( r_loc, f, df_theo, d2f_theo, df_num, d2f_num )
            deallocate( df_num_matmul, dcheb, d2cheb )
         end do

         close(file_handle)
      end if

   end subroutine test_radial_der
!------------------------------------------------------------------------------
end module tests
