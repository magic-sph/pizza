module outputs_phase
   !
   ! This module handles the computation and the writing of phase-field
   ! related diagnostics
   !

   use precision_mod
   use parallel_mod
   use fields, only: dtemp_Rloc
   use fourier, only: ifft
   use communications, only: gather_from_Rloc
   use constants, only: half, one, two, pi
   use radial_functions, only: r, rscheme, tcond, dtcond, m_R
   use namelists, only: r_cmb, tag, ra, pr, raxi, sc, radratio, stef, &
       &                tmelt, ek
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, n_phi_max, minc
   use blocking, only: nRstart, nRstop, nR_per_rank, load, getBlocks, &
       &               radial_balance
   use integration, only: rInt_R

   implicit none

   private

   integer :: n_phase_file, n_rmelt_file
   integer :: nPstart, nPstop
   type(load), allocatable :: phi_balance(:)
   real(cp), allocatable :: rmelt_loc(:), dt_rmelt_loc(:)
   real(cp), allocatable :: ekinS(:), ekinL(:), volS(:)
   real(cp), allocatable :: phase_phys_Rloc(:,:), temp_phys_Rloc(:,:)
   real(cp), allocatable :: phase_phys_Ploc(:,:), temp_phys_Ploc(:,:)
   real(cp), allocatable :: dtemp_phys_Rloc(:,:), dtemp_phys_Ploc(:,:)

   public :: initialize_out_phase, finalize_out_phase, calc_out_phase, &
   &         write_out_phase

contains

   subroutine initialize_out_phase()
      !
      ! Allocate the arrays needed to handle phase field outputs and
      ! open the output files.
      !

      real(cp) :: fac
      integer, parameter :: version=1
      integer :: n_phi

      allocate( ekinS(nRstart:nRstop), ekinL(nRstart:nRstop) )
      allocate( volS(nRstart:nRstop) )
      ekinS(:)=0.0_cp
      ekinL(:)=0.0_cp
      volS(:) =0.0_cp
      bytes_allocated=bytes_allocated+3*nR_per_rank*SIZEOF_DEF_REAL

      allocate( phase_phys_Rloc(n_phi_max,nRstart:nRstop) )
      allocate( temp_phys_Rloc(n_phi_max,nRstart:nRstop) )
      allocate( dtemp_phys_Rloc(n_phi_max,nRstart:nRstop) )
      phase_phys_Rloc(:,:)=0.0_cp
      temp_phys_Rloc(:,:) =0.0_cp
      dtemp_phys_Rloc(:,:)=0.0_cp
      bytes_allocated=bytes_allocated+3*n_phi_max*nR_per_rank*SIZEOF_DEF_REAL

      !-- Distribute phi over the ranks
      allocate(phi_balance(0:n_procs-1))
      call getBlocks(phi_balance, n_phi_max, n_procs)
      nPstart = phi_balance(rank)%nStart
      nPstop = phi_balance(rank)%nStop

      allocate( phase_phys_Ploc(nPstart:nPstop,n_r_max) )
      allocate( temp_phys_Ploc(nPstart:nPstop,n_r_max) )
      allocate( dtemp_phys_Ploc(nPstart:nPstop,n_r_max) )
      bytes_allocated=bytes_allocated+3*n_r_max*(nPstop-nPstart+1)*SIZEOF_DEF_REAL

      allocate( rmelt_loc(nPstart:nPstop), dt_rmelt_loc(nPstart:nPstop) )
      rmelt_loc(:)   =0.0_cp
      dt_rmelt_loc(:)=0.0_cp
      bytes_allocated=bytes_allocated+2*(nPstop-nPstart+1)*SIZEOF_DEF_REAL

      !-- Now save the melting line into a binary file
      if ( rank == 0 ) then
         open(newunit=n_phase_file, file='phase.'//tag, status='new')
         open(newunit=n_rmelt_file, file='rmelt.'//tag, status='new', &
         &    form='unformatted', access='stream')

         write(n_rmelt_file) version
         write(n_rmelt_file) n_phi_max
         fac=two*pi/real(n_phi_max*minc,cp)
         write(n_rmelt_file) ((n_phi-1)*fac, n_phi=1,n_phi_max)
         write(n_rmelt_file) ra, ek, pr, radratio, raxi, sc, stef, tmelt
      end if

   end subroutine initialize_out_phase
!---------------------------------------------------------------------------------------
   subroutine finalize_out_phase()
      !
      ! Close the files related to phase field diagnostics and deallocates
      ! the corresponding arrays.
      !

      if ( rank == 0 ) then
         close(n_rmelt_file)
         close(n_phase_file)
      end if

      deallocate( rmelt_loc, dtemp_phys_Rloc, dtemp_phys_Ploc )
      deallocate( dt_rmelt_loc, volS, ekinS, ekinL )
      deallocate( phase_phys_Rloc, temp_phys_Rloc, phase_phys_Ploc, temp_phys_Ploc)
      deallocate( phi_balance )

   end subroutine finalize_out_phase
!---------------------------------------------------------------------------------------
   subroutine calc_out_phase(us, up, phase, temp, nR)
      !
      ! This routine copies the arrays on the physical grid into global 2D arrays.
      ! It also computes the kinetic energy content in the solid and liquid phases.
      !

      !-- Input variables
      integer :: nR
      real(cp), intent(in) :: us(:)
      real(cp), intent(in) :: up(:)
      real(cp), intent(in) :: phase(:)
      real(cp), intent(in) :: temp(:)

      !-- Local variables:
      integer :: n_phi
      real(cp) :: ekin, dphi

      !-- Copy into 2-D arrays
      phase_phys_Rloc(:,nR)= phase(:)
      temp_phys_Rloc(:,nR) = temp(:)+tcond(nR)

      !-- Get the temperature gradient on the grid and store it:
      call ifft(dtemp_Rloc(:,nR), dtemp_phys_Rloc(:,nR), m_R(nR))
      dtemp_phys_Rloc(:,nR)=dtemp_phys_Rloc(:,nR)+dtcond(nR)

      !-- Measure the kinetic energy content in both phases
      ekinS(nR)=0.0_cp
      ekinL(nR)=0.0_cp
      volS(nR) =0.0_cp
      dphi=two*pi/real(n_phi_max,cp)
      do n_phi=1,n_phi_max
         ekin=half*(us(n_phi)**2+up(n_phi)**2)
         if (phase(n_phi) >= half) then
            ekinS(nR)=ekinS(nR)+dphi*ekin*r(nR)
            volS(nR) =volS(nR) +dphi*r(nR)
         else
            ekinL(nR)=ekinL(nR)+dphi*ekin*r(nR)
         end if
      end do

   end subroutine calc_out_phase
!---------------------------------------------------------------------------------------
   subroutine write_out_phase(time)
      !
      ! This routine handles the final MPI gather/reduction of diagnostics and
      ! the writing of the output files: phase.TAG and rmelt.TAG
      !

      !-- Input variable
      real(cp), intent(in) :: time

      !-- Local variables
      integer :: n_p
      real(cp) :: rmelt_mean_loc, tmelt_mean_loc, ekL, ekS, vS
      real(cp) :: tmelt_loc, norm, rmelt_max_loc, rmelt_min_loc
      real(cp) :: rmelt_max, rmelt_min, rmelt_mean, tmelt_mean
      real(cp) :: phase_max_loc, phase_min_loc, phase_max, phase_min
      real(cp) :: ekinS_global(n_r_max), ekinL_global(n_r_max), volS_global(n_r_max)
      real(cp) :: rmelt_global(n_phi_max), dt_rmelt_global(n_phi_max)

      !-- MPI transposes
      call transp_R2Phi(temp_phys_Rloc, temp_phys_Ploc)
      call transp_R2Phi(dtemp_phys_Rloc, dtemp_phys_Ploc)
      call transp_R2Phi(phase_phys_Rloc, phase_phys_Ploc)

      !-- Compute melt radius and melt temperature along rmelt
      rmelt_mean_loc=0.0_cp
      tmelt_mean_loc=0.0_cp
      norm=one/n_phi_max
      do n_p=nPstart,nPstop
         call get_rmelt_tmelt(phase_phys_Ploc(n_p,:), temp_phys_Ploc(n_p,:), &
              &               dtemp_phys_Ploc(n_p,:), rmelt_loc(n_p),        &
              &               tmelt_loc, dt_rmelt_loc(n_p))
         rmelt_mean_loc=rmelt_mean_loc+norm*rmelt_loc(n_p)
         tmelt_mean_loc=tmelt_mean_loc+norm*tmelt_loc
      end do

      rmelt_max_loc=maxval(rmelt_loc)
      rmelt_min_loc=minval(rmelt_loc)

      !-- MPI reductions on rank 0
      call MPI_Reduce(rmelt_mean_loc, rmelt_mean, 1, MPI_DEF_REAL, MPI_SUM, &
           &          0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(tmelt_mean_loc, tmelt_mean, 1, MPI_DEF_REAL, MPI_SUM, &
           &          0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(rmelt_min_loc, rmelt_min, 1, MPI_DEF_REAL, MPI_MIN, &
           &          0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(rmelt_max_loc, rmelt_max, 1, MPI_DEF_REAL, MPI_MAX, &
           &          0, MPI_COMM_WORLD, ierr)

      !-- MPI gather r-distributed arrays on rank=0
      call gather_from_Rloc(ekinS,ekinS_global,0)
      call gather_from_Rloc(ekinL,ekinL_global,0)
      call gather_from_Rloc(volS,volS_global,0)

      !-- Phase field min/max
      phase_max_loc=maxval(phase_phys_Rloc)
      phase_min_loc=minval(phase_phys_RLoc)
      call MPI_Reduce(phase_max_loc, phase_max, 1, MPI_DEF_REAL, MPI_MAX, &
           &          0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(phase_min_loc, phase_min, 1, MPI_DEF_REAL, MPI_MIN, &
           &          0, MPI_COMM_WORLD, ierr)

      !-- MPI gather phi-distributed array on rank=0
      call gather_from_Ploc(rmelt_loc,rmelt_global)
      call gather_from_Ploc(dt_rmelt_loc,dt_rmelt_global)

      !-- Write outputs
      if ( rank == 0 ) then
         ekL = rInt_R(ekinL_global, r, rscheme)
         ekS = rInt_R(ekinS_global, r, rscheme)
         vS = rInt_R(volS_global, r, rscheme)

         write(n_phase_file, '(ES20.12,7ES16.8,2ES13.5)') time, rmelt_mean, &
         &     tmelt_mean, rmelt_min, rmelt_max, vS, ekS, ekL, phase_min,   &
         &     phase_max

         write(n_rmelt_file) time
         write(n_rmelt_file) rmelt_global
         write(n_rmelt_file) dt_rmelt_global
      end if

   end subroutine write_out_phase
!---------------------------------------------------------------------------------------
   subroutine get_rmelt_tmelt(phase, temp, dtemp, rphase, tphase, dtphase)
      !
      ! This subroutine determines the melting point by approximating it by
      ! the radius where phi=0.5. It returns the radius and the temperature.
      ! It computes a 4th order Lagrangian interpolation between the four closest
      ! radii.
      !

      !-- Input variables
      real(cp), intent(in) :: phase(:) ! Phase field
      real(cp), intent(in) :: temp(:)  ! Temperature
      real(cp), intent(in) :: dtemp(:)  ! Radial derivative of temperature

      !-- Output variables
      real(cp), intent(out) :: rphase ! Radius of the melting point
      real(cp), intent(out) :: tphase ! Temperature along the interface
      real(cp), intent(out) :: dtphase ! Temperature gradient along the interface

      !-- Local variables
      integer :: n_r, n_r_phase, n_r_start, n_r_stop
      real(cp) :: x(4), y(4)

      !-- Determine the radial level where \phi=0.5
      n_r_phase=1
      do n_r=2,n_r_max
         if ( phase(n_r) < half .and. phase(n_r-1) > half ) then
            n_r_phase=n_r
         end if
      end do
      !-- 4th order Lagrangian interpolation of melting point
      if ( n_r_phase == 1 ) then
         rphase =r_cmb
         tphase =temp(1)
         dtphase=dtemp(1)
      else
         if ( n_r_phase == 2 ) then
            n_r_start=n_r_phase-1
            n_r_stop =n_r_phase+2
         else if ( n_r_phase == n_r_max ) then
            n_r_start=n_r_phase-3
            n_r_stop =n_r_phase
         else
            n_r_start=n_r_phase-2
            n_r_stop =n_r_phase+1
         end if
         x(:)=phase(n_r_start:n_r_stop)
         y(:)=r(n_r_start:n_r_stop)
         rphase=lagrange_interp(x,half,y)
         x(:)=r(n_r_start:n_r_stop)
         y(:)=temp(n_r_start:n_r_stop)
         tphase=lagrange_interp(x,rphase,y)
         y(:)=dtemp(n_r_start:n_r_stop)
         dtphase=lagrange_interp(x,rphase,y)
      end if

   end subroutine get_rmelt_tmelt
!---------------------------------------------------------------------------------------
   subroutine transp_R2Phi(arr_Rloc, arr_Ploc)
      !
      ! This subroutine is used to compute a MPI transpose between a R-distributed
      ! array and a Phi-distributed array
      !

      !-- Input array
      real(cp), intent(in) :: arr_Rloc(n_phi_max,nRstart:nRstop)

      !-- Output array
      real(cp), intent(out) :: arr_Ploc(nPstart:nPstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_p
      integer, allocatable :: rcounts(:), scounts(:), rdisp(:), sdisp(:)
      real(cp), allocatable :: sbuff(:), rbuff(:)
      integer :: p, ii, my_phi_counts

      !-- Set displacements vectors and buffer sizes
      allocate( rcounts(0:n_procs-1), scounts(0:n_procs-1) )
      allocate( rdisp(0:n_procs-1), sdisp(0:n_procs-1) )
      do p=0,n_procs-1
         my_phi_counts=phi_balance(p)%n_per_rank
         scounts(p)=nR_per_rank*my_phi_counts
         rcounts(p)=radial_balance(p)%n_per_rank*(nPStop-nPStart+1)
      end do

      rdisp(0)=0
      sdisp(0)=0
      do p=1,n_procs-1
         sdisp(p)=sdisp(p-1)+scounts(p-1)
         rdisp(p)=rdisp(p-1)+rcounts(p-1)
      end do
      allocate( sbuff(sum(scounts)), rbuff(sum(rcounts)) )
      sbuff(:)=0.0_cp
      rbuff(:)=0.0_cp

      !-- Prepare buffer
      do p=0,n_procs-1
         ii=sdisp(p)+1
         do n_r=nRstart,nRstop
            do n_p=phi_balance(p)%nStart,phi_balance(p)%nStop
               sbuff(ii)=arr_Rloc(n_p,n_r)
               ii=ii+1
            end do
         end do
      end do

      !-- All to all
      call MPI_Alltoallv(sbuff, scounts, sdisp, MPI_DEF_REAL, &
           &             rbuff, rcounts, rdisp, MPI_DEF_REAL, &
           &             MPI_COMM_WORLD, ierr)

      !-- Reassemble array
      do p=0,n_procs-1
         ii=rdisp(p)+1
         do n_r=radial_balance(p)%nStart,radial_balance(p)%nStop
            do n_p=nPstart,nPstop
               arr_Ploc(n_p,n_r)=rbuff(ii)
               ii=ii+1
            end do
         end do
      end do

      !-- Clear memory from temporary arrays
      deallocate( rcounts, scounts, rdisp, sdisp, rbuff, sbuff )

   end subroutine transp_R2Phi
!---------------------------------------------------------------------------------------
   real(cp) function lagrange_interp(xp, x, yp)
      !
      ! This function performs a Lagrange interpolation around the point
      ! x. The order depends on the size of the input arrays x and y
      !

      !-- Input variables
      real(cp), intent(in) :: xp(:) ! Grid points where the quantity is known
      real(cp), intent(in) :: x ! Point where the quantity is interpolated
      real(cp), intent(in) :: yp(:) ! value

      !-- Local variables
      real(cp) :: lag_i ! Lagrange polynomial of order n
      integer :: i, j, n

      n = size(xp) ! Degree of the Lagrange interpolant
      lagrange_interp=0.0_cp
      do i=1,n
         lag_i=one
         do j=1,n
            if (i /= j) lag_i = lag_i*(x - xp(j))/(xp(i) - xp(j))
         end do
         lagrange_interp=lagrange_interp+lag_i*yp(i)
      end do

   end function lagrange_interp
!----------------------------------------------------------------------------
   subroutine gather_from_Ploc(arr_Ploc, arr_full)
      !
      ! This subroutine gathers a phi-distributed array on rank 0.
      !

      !-- Input variable
      real(cp), intent(in) :: arr_Ploc(nPstart:nPstop)

      !-- Output variable
      real(cp), intent(out) :: arr_full(n_phi_max)

      !-- Local variables
      integer :: rcounts(0:n_procs-1), rdisp(0:n_procs-1), scount
      integer :: p

      scount = (nPstop-nPstart+1)
      do p=0,n_procs-1
         rcounts(p)=phi_balance(p)%n_per_rank
      end do
      rdisp(0)=0
      do p=1,n_procs-1
         rdisp(p)=rdisp(p-1)+rcounts(p-1)
      end do

      call MPI_GatherV(arr_Ploc, scount, MPI_DEF_REAL, arr_full, rcounts, rdisp, &
           &           MPI_DEF_REAL, 0, MPI_COMM_WORLD, ierr)

   end subroutine gather_from_Ploc
!----------------------------------------------------------------------------
end module outputs_phase
