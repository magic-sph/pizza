module tests
   !
   ! This module contains several testing subroutines
   !
   use precision_mod
   use parallel_mod
   use constants, only: one, half
   use namelists, only: l_newmap, radratio
   use chebyshev, only: type_cheb
   use radial_functions, only: rscheme
   use radial_der, only: initialize_der_arrays, finalize_der_arrays, get_ddr
   use radial_scheme, only: type_rscheme
   use algebra, only: sgefa, rgesl, prepare_bordered_mat, solve_bordered_mat
   use useful, only: abortRun

   implicit none

   private

   public :: solve_laplacian, test_radial_der

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
         &      1536,  2048]
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
            call solve_laplacian_colloc(n_r_max_loc, r_loc, rscheme, sol)
            tStop = MPI_Wtime()
            tColloc = tStop-tStart
            errColloc = maxval(abs(sol(:)-sol_theo(:)))

            tStart = MPI_Wtime()
            call solve_laplacian_integ(n_r_max_loc, r_loc, rscheme, sol)
            tStop = MPI_Wtime()
            tInteg = tStop-tStart
            errInteg = maxval(abs(sol(:)-sol_theo(:)))

            !-- Store it in a file
            write(file_handle, '(i5,5es20.12)') n_r_max_loc, errColloc, &
            &                                   tColloc, errInteg, tInteg
            write(6, '(i5,5es20.12)') n_r_max_loc, errColloc, tColloc,  &
                                      errInteg, tInteg


            call finalize_der_arrays()
            call rscheme%finalize()
            deallocate( r_loc, sol, sol_theo )
            ! deallocate( r_loc, f, df_theo, d2f_theo, df_num, d2f_num )
         end do

         close(file_handle)

      end if

   end subroutine solve_laplacian
!------------------------------------------------------------------------------
   subroutine solve_laplacian_integ(n_r_max, r, rscheme, rhs)

      !-- Input variables
      integer,             intent(in) :: n_r_max
      real(cp),            intent(in) :: r(n_r_max)
      class(type_rscheme), intent(in) :: rscheme

      !-- Output variable
      real(cp),            intent(out) :: rhs(n_r_max)

      !-- Local variables
      real(cp), allocatable :: tmp(:)
      real(cp), allocatable :: Bmat(:,:), A4mat(:,:)
      real(cp), allocatable :: A1mat(:,:), A2mat(:,:), A3mat(:,:)
      integer :: pivotA1(2)
      integer, allocatable :: pivotA4(:)
      real(cp) :: a, b
      integer :: n_band, n_r, n_bands_Bmat, n_bands_Amat

      n_bands_Bmat = 7
      !n_bands_Amat = 3
      n_bands_Amat = 4

      allocate ( A4mat(n_bands_Amat, n_r_max-2) )
      allocate ( A1mat(2,2), A2mat(2,n_r_max-2), A3mat(n_r_max-2, 2) )
      allocate ( Bmat(n_bands_Bmat, n_r_max) )
      allocate ( tmp(n_r_max), pivotA4(n_r_max) )
      
      !-- Set the matrices to zero
      A1mat(:,:)=0.0_cp
      A2mat(:,:)=0.0_cp
      A3mat(:,:)=0.0_cp
      A4mat(:,:)=0.0_cp
      Bmat(:,:)=0.0_cp
      
      a = half*(r(1)-r(n_r_max))
      b = half*(r(1)+r(n_r_max))

      do n_r=1,n_r_max
         if ( n_r > 2 .and. n_r < n_r_max-2 ) Bmat(1,n_r+3)= a**3/8.0_cp/n_r/(n_r-1)
         if ( n_r > 2 .and. n_r < n_r_max-1 ) Bmat(2,n_r+2)= a**2*b/4.0_cp/n_r/(n_r-1)
         if ( n_r > 2 .and. n_r < n_r_max ) then
            Bmat(3,n_r+1)= -a**3/8.0_cp/(n_r-1)/(n_r-2)
            A4mat(2,n_r-1)= half*a*(one+one/(n_r-1))
         end if
         if ( n_r > 2 ) then
            Bmat(4,n_r)= -a**2*b/2.0_cp/(n_r)/(n_r-2)
            A4mat(3,n_r-2)= b
         end if
         if ( n_r > 2 ) Bmat(5,n_r-1)= -a**3/8.0_cp/(n_r)/(n_r-1)
         if ( n_r > 3 ) A4mat(4,n_r-3)= half*a*(one-one/(n_r-1))
         if ( n_r > 2 ) Bmat(6,n_r-2)= a**2*b/4.0_cp/(n_r-1)/(n_r-2)
         if ( n_r > 3 ) Bmat(7,n_r-3)= a**3/8.0_cp/(n_r-1)/(n_r-2)
      end do
      A3mat(1,2)=0.25_cp*rscheme%rnorm*a


      do n_r=1,n_r_max
         do n_band=1,n_bands_Bmat
            Bmat(n_band,n_r)=rscheme%rnorm*Bmat(n_band,n_r)
         end do
      end do
      do n_r=1,n_r_max-2
         do n_band=1,n_bands_Amat
            A4mat(n_band,n_r)=rscheme%rnorm*A4mat(n_band,n_r)
         end do
      end do

      !-- Add the tau Lines for boundary conditions
      do n_r=1,n_r_max
         if ( n_r <= 2) then
            A1mat(1,n_r)=rscheme%rnorm*rscheme%rMat(1,n_r)
            A1mat(2,n_r)=rscheme%rnorm*rscheme%rMat(n_r_max,n_r)
         else
            A2mat(1,n_r-2)=rscheme%rnorm*rscheme%rMat(1,n_r)
            A2mat(2,n_r-2)=rscheme%rnorm*rscheme%rMat(n_r_max,n_r)
         end if
      end do
      A1mat(1,1)=rscheme%boundary_fac*A1mat(1,1)
      A1mat(2,1)=rscheme%boundary_fac*A1mat(2,1)
      A2mat(1,n_r_max-2)=rscheme%boundary_fac*A2mat(1,n_r_max-2)
      A2mat(2,n_r_max-2)=rscheme%boundary_fac*A2mat(2,n_r_max-2)

      !-- Assemble right-hand side
      tmp(:)       = 3.0_cp
      call rscheme%costf1(tmp, n_r_max)
      call dgbmv('N', n_r_max, n_r_max, 3, 3, one, Bmat, n_bands_Bmat, tmp, &
           &     1, 0.0_cp, rhs, 1)
      !-- Boundary conditions
      rhs(1) = -0.25_cp
      rhs(2) = 0.75_cp

      call prepare_bordered_mat(A1mat,A2mat,A3mat,A4mat,2,n_r_max-2,1,1, &
           &                    pivotA1, pivotA4)

      call solve_bordered_mat(A1mat,A2mat,A3mat,A4mat,2,n_r_max-2,1,1,  &
           &                  pivotA1, pivotA4, rhs, n_r_max)

      call rscheme%costf1(rhs, n_r_max)

      deallocate ( A4mat, A1mat, A2mat, A3mat, Bmat, tmp, pivotA4 )

   end subroutine solve_laplacian_integ
!------------------------------------------------------------------------------
   subroutine solve_laplacian_colloc(n_r_max, r, rscheme, rhs)

      !-- Input variables
      integer,             intent(in) :: n_r_max
      real(cp),            intent(in) :: r(n_r_max)
      class(type_rscheme), intent(in) :: rscheme

      !-- Output variable
      real(cp),            intent(out) :: rhs(n_r_max)

      !-- Local variables
      real(cp), allocatable :: mat(:,:)
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
      call sgefa(mat,n_r_max,n_r_max,pivot,info)
      if ( info /= 0 ) then
         call abortRun('Singular matrix mat!')
      end if

      call rgesl(mat, n_r_max, n_r_max, pivot, rhs)
      call rscheme%costf1(rhs, n_r_max)

      deallocate( mat, pivot )

   end subroutine solve_laplacian_colloc
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
