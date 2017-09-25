module tests
   !
   ! This module contains several testing subroutines
   !

   use precision_mod
   use parallel_mod
   use constants, only: one
   use namelists, only: l_newmap, radratio
   use chebyshev, only: type_cheb
   use radial_functions, only: rscheme
   use radial_der, only: initialize_der_arrays, finalize_der_arrays, get_ddr
   use radial_scheme, only: type_rscheme
   use algebra, only: sgefa, rgesl
   use useful, only: abortRun

   public :: solve_laplacian_collocation, test_radial_der

contains

   subroutine solve_laplacian_collocation(nMstart, nMstop)

      !-- Input variables
      integer, intent(in) :: nMstart
      integer, intent(in) :: nMstop

      !-- Local variables
      integer, allocatable :: nrs(:)
      real(cp), allocatable :: r_loc(:), sol(:), sol_theo(:)
      integer :: file_handle, n_r_max_loc, n_in
      real(cp) :: r_cmb, r_icb, eps, alph1, alph2, c1, c2, err

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

            call assemble_and_solve_lap(n_r_max_loc, r_loc, rscheme, sol)

            err = maxval(abs(sol(:)-sol_theo(:)))

            !-- Store it in a file
            write(file_handle, '(i5,2es20.12)') n_r_max_loc, err
            write(6, '(i5,2es20.12)') n_r_max_loc, err


            call finalize_der_arrays()
            call rscheme%finalize()
            deallocate( r_loc, sol, sol_theo )
            ! deallocate( r_loc, f, df_theo, d2f_theo, df_num, d2f_num )
         end do

         close(file_handle)

      end if

   end subroutine solve_laplacian_collocation
!------------------------------------------------------------------------------
   subroutine assemble_and_solve_lap(n_r_max, r, rscheme, rhs)

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

      !if ( rscheme%n_max < n_r_max ) then ! fill with zeros !
      !   do nR_out=rscheme%n_max+1,n_r_max
      !      mat(1,nR_out)      =0.0_cp
      !      mat(n_r_max,nR_out)=0.0_cp
      !   end do
      !end if

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

   end subroutine assemble_and_solve_lap
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
