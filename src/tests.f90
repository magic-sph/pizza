module tests
   !
   ! This module contains several testing subroutines
   !
   use precision_mod
   use parallel_mod
   use constants, only: one, half, two, pi
   use namelists, only: l_newmap, radratio
   use chebyshev, only: type_cheb
   use radial_functions, only: rscheme
   use radial_der, only: initialize_der_arrays, finalize_der_arrays, get_ddr
   use radial_scheme, only: type_rscheme
   use algebra, only: sgefa, rgesl, prepare_bordered_mat, solve_bordered_mat
   use useful, only: abortRun
   use chebsparselib, only: intcheb2rmult1, rmult1, intcheb1, intcheb4, &
       &                    intcheb2, eye

   implicit none

   private

   public :: solve_laplacian, solve_biharmo, test_radial_der

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
      real(cp) :: errColloc, errInteg
      real(cp) :: tStart, tStop, tColloc, tInteg
      real(cp) :: timeLuColl, timeLuInt, timeSolveColl, timeSolveInt

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

         do n_r=1,size(nrs)
            n_r_max_loc = nrs(n_r)
            allocate( r_loc(n_r_max_loc), sol(n_r_max_loc) )
            allocate( sol_theo(n_r_max_loc) )
            ! allocate( df_num(n_r_max_loc), d2f_num(n_r_max_loc) )
            ! allocate( df_theo(n_r_max_loc), d2f_theo(n_r_max_loc) )
            call rscheme%initialize(n_r_max_loc,n_r_max_loc, n_in)
            if ( l_newmap ) then
               alph1=1.0_cp/cosh(abs(log(eps))/(n_r_max_loc-1))
               alph2=0.0_cp
            else
               alph1=0.0_cp
               alph2=0.0_cp
            end if
            call rscheme%get_grid(n_r_max_loc, r_icb, r_cmb, alph1, alph2, r_loc)
            call initialize_der_arrays(n_r_max_loc, nMstart, nMstop)
            call rscheme%get_der_mat(n_r_max_loc)

            c1 = (-1.0_cp-0.75_cp*(r_cmb**2-r_icb**2))/log(r_cmb/r_icb)
            c2 = -0.25_cp-0.75_cp*r_cmb**2-c1*log(r_cmb)
            sol_theo(:)=0.75_cp*r_loc(:)**2+c1*log(r_loc)+c2

            tStart = MPI_Wtime()
            call solve_laplacian_colloc(n_r_max_loc, r_loc, rscheme, sol, &
                 &                      timeLuColl, timeSolveColl)
            tStop = MPI_Wtime()
            tColloc = tStop-tStart
            errColloc = maxval(abs(sol(:)-sol_theo(:)))

            tStart = MPI_Wtime()
            call solve_laplacian_integ(n_r_max_loc, r_loc, rscheme, sol, &
                 &                     timeLuInt, timeSolveInt)
            tStop = MPI_Wtime()
            tInteg = tStop-tStart
            errInteg = maxval(abs(sol(:)-sol_theo(:)))

            !-- Store it in a file
            write(file_handle, '(i5,8es14.6)') n_r_max_loc, errColloc,   &
            &                                   tColloc, timeLuColl,     &
            &                                   timeSolveColl, errInteg, &
            &                                   tInteg, timeLuInt, timeSolveInt
            write(6, '(i5,8es14.6)') n_r_max_loc, errColloc, tColloc,     &
            &                        timeLuColl, timeSolveColl, errInteg, &
            &                        tInteg, timeLuInt, timeSolveInt


            call finalize_der_arrays()
            call rscheme%finalize()
            deallocate( r_loc, sol, sol_theo )
            ! deallocate( r_loc, f, df_theo, d2f_theo, df_num, d2f_num )
         end do

         close(file_handle)

      end if

   end subroutine solve_laplacian
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
      real(cp) :: errColloc, errInteg
      real(cp) :: tStart, tStop, tColloc, tInteg
      real(cp) :: timeLuColl, timeLuInt, timeSolveColl, timeSolveInt

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
            call rscheme%initialize(n_r_max_loc,n_r_max_loc, n_in)
            if ( l_newmap ) then
               alph1=1.0_cp/cosh(abs(log(eps))/(n_r_max_loc-1))
               alph2=0.0_cp
            else
               alph1=0.0_cp
               alph2=0.0_cp
            end if
            call rscheme%get_grid(n_r_max_loc, r_icb, r_cmb, alph1, alph2, r_loc)
            call initialize_der_arrays(n_r_max_loc, nMstart, nMstop)
            call rscheme%get_der_mat(n_r_max_loc)

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
                 &                     timeLuInt, timeSolveInt)
            tStop = MPI_Wtime()
            tInteg = tStop-tStart
            errInteg = maxval(abs(sol(:)-sol_theo(:)))

            !-- Store it in a file
            write(file_handle, '(i5,8es14.6)') n_r_max_loc, errColloc,   &
            &                                   tColloc, timeLuColl,     &
            &                                   timeSolveColl, errInteg, &
            &                                   tInteg, timeLuInt, timeSolveInt
            write(6, '(i5,8es14.6)') n_r_max_loc, errColloc, tColloc,     &
            &                        timeLuColl, timeSolveColl, errInteg, &
            &                        tInteg, timeLuInt, timeSolveInt


            call finalize_der_arrays()
            call rscheme%finalize()
            deallocate( r_loc, sol, sol_theo )
            ! deallocate( r_loc, f, df_theo, d2f_theo, df_num, d2f_num )
         end do

         close(file_handle)

      end if

   end subroutine solve_biharmo
!------------------------------------------------------------------------------
   subroutine solve_laplacian_integ(n_r_max, r, rscheme, rhs, timeLu, timeSolve)

      !-- Input variables
      integer,             intent(in) :: n_r_max
      real(cp),            intent(in) :: r(n_r_max)
      class(type_rscheme), intent(in) :: rscheme

      !-- Output variable
      real(cp),            intent(out) :: rhs(n_r_max)
      real(cp),            intent(out) :: timeLu
      real(cp),            intent(out) :: timeSolve

      !-- Local variables
      real(cp), allocatable :: tmp(:), stencilA4(:), stencilB(:)
      real(cp), allocatable :: Bmat(:,:), A4mat(:,:)
      real(cp), allocatable :: A1mat(:,:), A2mat(:,:), A3mat(:,:)
      integer, allocatable :: pivotA1(:), pivotA4(:)
      real(cp) :: a, b, tStart, tStop
      integer :: i_r, klA4, kuA4, klB, kuB, nStart
      integer :: n_band, n_r, n_bands_Bmat, n_bands_Amat, n_boundaries, lenA4

      klA4 = 1
      kuA4 = 1
      klB  = 3
      kuB  = 3
      n_bands_Bmat = klB+kuB+1
      !-- Factor 2 in front of klA4 is needed for LU factorisation
      n_bands_Amat = 2*klA4+kuA4+1
      n_boundaries = 2
      lenA4 = n_r_max-n_boundaries

      allocate ( A4mat(n_bands_Amat, lenA4) )
      allocate ( stencilA4(klA4+kuA4+1), stencilB(klB+kuB+1) )
      allocate ( A1mat(n_boundaries,n_boundaries), A2mat(n_boundaries,lenA4) )
      allocate ( A3mat(lenA4, n_boundaries) )
      allocate ( Bmat(n_bands_Bmat, n_r_max) )
      allocate ( tmp(n_r_max), pivotA1(n_boundaries), pivotA4(lenA4) )
      
      !-- Set the matrices to zero
      A1mat(:,:)=0.0_cp
      A2mat(:,:)=0.0_cp
      A3mat(:,:)=0.0_cp
      A4mat(:,:)=0.0_cp
      Bmat(:,:) =0.0_cp
      
      a = half*(r(1)-r(n_r_max))
      b = half*(r(1)+r(n_r_max))

      !-- Fill A4 banded-block
      do n_r=1,lenA4
         i_r = n_r+n_boundaries
         
         !-- Define the equations
         stencilA4 = rmult1(a,b,i_r-1,klA4+kuA4+1)-intcheb1(a,i_r-1,klA4+kuA4+1)

         !-- Roll the array for band storage
         do n_band=1,klA4+kuA4+1
            if ( n_r+kuA4+1-n_band <= lenA4 .and. n_r+kuA4+1-n_band >= 1 ) then
               A4mat(klA4+n_band,n_r+kuA4+1-n_band) = rscheme%rnorm*stencilA4(n_band)
            end if
         end do
      end do

      !-- Fill A3
      do n_r=1,lenA4
         i_r = n_r+n_boundaries
         stencilA4 = rmult1(a,b,i_r-1,klA4+kuA4+1)-intcheb1(a,i_r-1,klA4+kuA4+1)

         !-- Only the lower bands can contribute to the matrix A3
         do n_band=1,klA4
            if ( n_r <= n_band .and. n_r+n_boundaries-n_band >= 1 ) then
               A3mat(n_r,n_r+n_boundaries-n_band) = rscheme%rnorm*stencilA4(kuA4+1+n_band)
            end if
         end do
      end do

      !-- Fill right-hand side matrix
      do n_r=1,lenA4
         i_r = n_r+n_boundaries

         !-- Define the equations 
         stencilB = intcheb2rmult1(a,b,i_r-1,klB+kuB+1)

         !-- Roll the arrays for band storage
         do n_band=1,klB+kuB+1
            if ( i_r+kuB+1-n_band <= n_r_max .and. i_r+kuB+1-n_band >= 1 ) then
               Bmat(n_band,i_r+kuB+1-n_band) = rscheme%rnorm*stencilB(n_band)
            end if
         end do
      end do

      !-- Add the tau Lines for boundary conditions into A1 and A2
      do n_r=1,n_r_max
         if ( n_r <= n_boundaries) then
            A1mat(1,n_r)=rscheme%rnorm*rscheme%rMat(1,n_r)
            A1mat(2,n_r)=rscheme%rnorm*rscheme%rMat(n_r_max,n_r)
         else
            A2mat(1,n_r-n_boundaries)=rscheme%rnorm*rscheme%rMat(1,n_r)
            A2mat(2,n_r-n_boundaries)=rscheme%rnorm*rscheme%rMat(n_r_max,n_r)
         end if
      end do

      !-- Cheb factor for boundary conditions
      do n_r=1,n_boundaries
         A1mat(n_r,1)    =rscheme%boundary_fac*A1mat(n_r,1)
         A2mat(n_r,lenA4)=rscheme%boundary_fac*A2mat(n_r,lenA4)
      end do

      !-- Assemble right-hand side
      tmp(:)       = 3.0_cp
      call rscheme%costf1(tmp, n_r_max)
      nStart = n_boundaries+1
      call dgbmv('N', n_r_max, n_r_max, klB, kuB, one, Bmat, n_bands_Bmat, &
           &     tmp, 1, 0.0_cp, rhs, 1)
      !-- Boundary conditions
      rhs(1) = -0.25_cp
      rhs(2) = 0.75_cp

      !-- LU factorisation
      tStart = MPI_Wtime()
      call prepare_bordered_mat(A1mat,A2mat,A3mat,A4mat,n_boundaries,lenA4, &
           &                    klA4,kuA4, pivotA1, pivotA4)
      tStop = MPI_Wtime()
      timeLu= tStop-tStart

      !-- Solve
      tStart = MPI_Wtime()
      call solve_bordered_mat(A1mat,A2mat,A3mat,A4mat,n_boundaries,lenA4, &
           &                  klA4, kuA4, pivotA1, pivotA4, rhs, n_r_max)
      tStop = MPI_Wtime()
      timeSolve= tStop-tStart

      !-- Final DCT to bring the solution back to physical space
      call rscheme%costf1(rhs, n_r_max)

      deallocate ( A4mat, A1mat, A2mat, A3mat, Bmat, tmp, pivotA4 )
      deallocate ( stencilA4, stencilB )

   end subroutine solve_laplacian_integ
!------------------------------------------------------------------------------
   subroutine solve_laplacian_colloc(n_r_max, r, rscheme, rhs, timeLu, &
              &                      timeSolve)

      !-- Input variables
      integer,             intent(in) :: n_r_max
      real(cp),            intent(in) :: r(n_r_max)
      class(type_rscheme), intent(in) :: rscheme

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
      rhs(:)       = 3.0_cp
      rhs(1)       = -0.25_cp
      rhs(n_r_max) = 0.75_cp

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
      call sgefa(mat,n_r_max,n_r_max,pivot,info)
      tStop = MPI_Wtime()
      timeLu = tStop-tStart
      if ( info /= 0 ) then
         call abortRun('Singular matrix mat!')
      end if

      tStart = MPI_Wtime()
      call rgesl(mat, n_r_max, n_r_max, pivot, rhs)
      tStop = MPI_Wtime()
      timeSolve = tStop-tStart
      call rscheme%costf1(rhs, n_r_max)

      deallocate( mat, pivot )

   end subroutine solve_laplacian_colloc
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
      real(cp), allocatable :: tmp(:), stencilB(:), stencilA4(:)
      real(cp), allocatable :: Bmat(:,:), A4mat(:,:)
      real(cp), allocatable :: A1mat(:,:), A2mat(:,:), A3mat(:,:)
      integer, allocatable :: pivotA1(:), pivotA4(:)
      real(cp) :: a, b, tStart, tStop
      integer :: i_r, klA4, kuA4, klB, kuB, nStart, n_band, n_r
      integer :: n_bands_Bmat, n_bands_Amat, n_boundaries, lenA4

      klA4 = 4
      kuA4 = 4
      klB  = 4
      kuB  = 4
      n_boundaries = 4
      n_bands_Bmat = klB+kuB+1
      !-- Factor 2 in front of klA4 is needed for LU factorisation
      n_bands_Amat = 2*klA4+kuA4+1
      lenA4 = n_r_max-n_boundaries

      allocate ( stencilB(klB+kuB+1), stencilA4(klA4+kuA4+1) )
      allocate ( A4mat(n_bands_Amat, lenA4) )
      allocate ( A1mat(n_boundaries,n_boundaries), A2mat(n_boundaries,lenA4) )
      allocate ( A3mat(lenA4, n_boundaries) )
      allocate ( Bmat(n_bands_Bmat, n_r_max) )
      allocate ( tmp(n_r_max), pivotA1(n_boundaries), pivotA4(lenA4) )
      
      !-- Set the matrices to zero
      A1mat(:,:)=0.0_cp
      A2mat(:,:)=0.0_cp
      A3mat(:,:)=0.0_cp
      A4mat(:,:)=0.0_cp
      Bmat(:,:) =0.0_cp
      
      a = half*(r(1)-r(n_r_max))
      b = half*(r(1)+r(n_r_max))

      !-- Fill A4 banded-block
      do n_r=1,lenA4
         i_r = n_r+n_boundaries

         !-- Define the equations
         stencilA4 = eye(klA4+kuA4+1) + two*intcheb2(a,i_r-1,klA4+kuA4+1)+ &
         &           intcheb4(a,i_r-1,klA4+kuA4+1)

         !-- Roll the array for band storage
         do n_band=1,klA4+kuA4+1
            if ( n_r+kuA4+1-n_band <= lenA4 .and. n_r+kuA4+1-n_band >= 1 ) then
               A4mat(klA4+n_band,n_r+kuA4+1-n_band) = rscheme%rnorm*stencilA4(n_band)
            end if
         end do
      end do

      !-- Fill A3
      do n_r=1,lenA4
         i_r = n_r+n_boundaries
         stencilA4 = eye(klA4+kuA4+1) + two*intcheb2(a,i_r-1,klA4+kuA4+1)+ &
         &           intcheb4(a,i_r-1,klA4+kuA4+1)

         !-- Only the lower bands can contribute to the matrix A3
         do n_band=1,klA4
            if ( n_r <= n_band .and. n_r+n_boundaries-n_band >= 1 ) then
               A3mat(n_r,n_r+n_boundaries-n_band) = rscheme%rnorm*stencilA4(kuA4+1+n_band)
            end if
         end do
      end do

      !-- Fill right-hand side matrix
      do n_r=1,lenA4
         i_r = n_r+n_boundaries

         !-- Define right-hand side equations
         stencilB = intcheb4(a,i_r-1,klB+kuB+1)

         !-- Roll array for band storage
         do n_band=1,klB+kuB+1
            if ( i_r+kuB+1-n_band <= n_r_max .and. i_r+kuB+1-n_band >= 1 ) then
               Bmat(n_band,i_r+kuB+1-n_band) = rscheme%rnorm*stencilB(n_band)
            end if
         end do
      end do

      !-- Add the tau Lines for boundary conditions into A1 and A2
      do n_r=1,n_r_max
         if ( n_r <= n_boundaries) then
            A1mat(1,n_r)=rscheme%rnorm*rscheme%rMat(1,n_r)
            A1mat(2,n_r)=rscheme%rnorm*rscheme%rMat(n_r_max,n_r)
            A1mat(3,n_r)=rscheme%rnorm*rscheme%drMat(1,n_r)
            A1mat(4,n_r)=rscheme%rnorm*rscheme%drMat(n_r_max,n_r)
         else
            A2mat(1,n_r-n_boundaries)=rscheme%rnorm*rscheme%rMat(1,n_r)
            A2mat(2,n_r-n_boundaries)=rscheme%rnorm*rscheme%rMat(n_r_max,n_r)
            A2mat(3,n_r-n_boundaries)=rscheme%rnorm*rscheme%drMat(1,n_r)
            A2mat(4,n_r-n_boundaries)=rscheme%rnorm*rscheme%drMat(n_r_max,n_r)
         end if
      end do

      !-- Cheb factors for boundary conditions
      do n_r=1,n_boundaries
         A1mat(n_r,1)    =rscheme%boundary_fac*A1mat(n_r,1)
         A2mat(n_r,lenA4)=rscheme%boundary_fac*A2mat(n_r,lenA4)
      end do

      !-- Assemble right-hand side
      tmp(:) = cos(r)
      call rscheme%costf1(tmp, n_r_max)
      nStart = n_boundaries+1
      call dgbmv('N', n_r_max, n_r_max, klB, kuB, one, Bmat, n_bands_Bmat, &
           &     tmp, 1, 0.0_cp, rhs, 1)
      !-- Boundary conditions
      rhs(1) = 0.0_cp
      rhs(2) = 0.0_cp
      rhs(3) = 0.0_cp
      rhs(4) = 0.0_cp

      !-- LU factorisation of A matrix
      tStart = MPI_Wtime()
      call prepare_bordered_mat(A1mat,A2mat,A3mat,A4mat,n_boundaries,lenA4, &
           &                    klA4,kuA4, pivotA1, pivotA4)
      tStop = MPI_Wtime()
      timeLu= tStop-tStart

      !-- Solve
      tStart = MPI_Wtime()
      call solve_bordered_mat(A1mat,A2mat,A3mat,A4mat,n_boundaries,lenA4, &
           &                  klA4, kuA4, pivotA1, pivotA4, rhs, n_r_max)
      tStop = MPI_Wtime()
      timeSolve= tStop-tStart

      !-- Final DCT to bring back the data to physical space
      call rscheme%costf1(rhs, n_r_max)

      deallocate ( A4mat, A1mat, A2mat, A3mat, Bmat, tmp, pivotA4 )
      deallocate ( stencilA4, stencilB )

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
      call sgefa(mat,2*n_r_max,2*n_r_max,pivot,info)
      tStop = MPI_Wtime()
      timeLu = tStop-tStart
      if ( info /= 0 ) then
         call abortRun('Singular matrix mat!')
      end if

      tStart = MPI_Wtime()
      call rgesl(mat, 2*n_r_max, 2*n_r_max, pivot, tmp)
      tStop = MPI_Wtime()
      timeSolve = tStop-tStart
      rhs(:)=tmp(1:n_r_max)
      call rscheme%costf1(rhs, n_r_max)

      deallocate( mat, pivot )

   end subroutine solve_biharmo_colloc
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
            call rscheme%initialize(n_r_max_loc,n_r_max_loc, n_in)
            if ( l_newmap ) then
               alph1=one/cosh(abs(log(eps))/(n_r_max_loc-1))
               alph2=0.0_cp
            else
               alph1=0.0_cp
               alph2=0.0_cp
            end if
            call rscheme%get_grid(n_r_max_loc, r_icb, r_cmb, alph1, alph2, r_loc)
            call initialize_der_arrays(n_r_max_loc, nMstart, nMstop)
            call rscheme%get_der_mat(n_r_max_loc)

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
            write(file_handle, '(i5,4es20.12)') n_r_max_loc, err_d1, err_d1_matmul,&
            &                                   err_d2, err_d2_matmul
            write(6, '(i5,4es20.12)') n_r_max_loc, err_d1, err_d1_matmul, err_d2, &
            &                         err_d2_matmul

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
