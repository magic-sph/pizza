module pre_calculations

   use iso_fortran_env, only: output_unit
   use truncation, only: n_r_max, n_m_max, n_phi_max, minc, m_max, &
   &                     n_cheb_max
   use constants, only: one, four, third, pi, vol_oc, surf, two, vol_otc
   use outputs, only: n_log_file
   use radial_functions, only: r, rscheme, radial, height
   use namelists, only: l_newmap, tag, dtMax, dtMin, pr, l_non_rot, ek, &
       &                r_cmb, r_icb, container
   use horizontal, only: mfunctions
   use integration, only: rInt_R
   use precision_mod
   use parallel_mod

   implicit none

   private
   
   public :: preCalc

contains

   subroutine preCalc

      !-- Local variables
      integer :: file_handle, n_r
      character(len=100) :: file_name

      dtMin = dtMax/1e6_cp

      !-- Set r-dependent arrays
      call radial()

      if ( ( l_newmap ) .and. (rank == 0) ) then
         file_name='rNM.'//tag
         open(newunit=file_handle, file=file_name, status='unknown')
         do n_r=1,n_r_max
            write(file_handle,'(I5,3ES16.8)') n_r, r(n_r),     &
            &                                rscheme%drx(n_r), &
            &                                rscheme%ddrx(n_r)
         end do
         close(file_handle)
      end if

      !-- Set m-dependent arrays
      call mfunctions()

      !-- Compute some constants
      vol_oc =four*third*pi*(r_cmb**3-r_icb**3)
      if ( index(container, 'SPHERE') == 1 ) then
         vol_otc=four*third*pi*(r_cmb**2-r_icb**2)**(1.5_cp)
      else if ( index(container, 'EXP') == 1 ) then
         vol_otc=two * pi * rInt_R(height(:)*r(:), r(:), rscheme)
      else if ( index(container, 'BUSSE') == 1 ) then
         vol_otc=two * pi * rInt_R(height(:)*r(:), r(:), rscheme)
      end if

      surf   =pi*(r_cmb**2-r_icb**2)

      !-- Write some informations
      if ( rank == 0 ) then
         call write_info(output_unit)
         call write_info(n_log_file)
      end if

   end subroutine preCalc
!---------------------------------------------------------------------------
   subroutine write_info(n_out)

      integer, intent(in) :: n_out

      write(n_out,*) ''
      write(n_out, '('' ! Spherical shell volume  :'',es14.6)') vol_oc
      write(n_out, '('' ! Volume outside tan. cyl.:'',es14.6)') vol_otc
      write(n_out, '('' ! Annulus surface         :'',es14.6)') surf

      write(n_out,*) ''
      write(n_out,*) '! MPI ranks:'
      write(n_out,'(''  n_procs      ='',i6)') n_procs

      write(n_out,*) ''
      write(n_out,*) '! Grid parameters:'
      write(n_out,'(''  n_r_max      ='',i6, &
           &   '' = number of radial grid points'')') n_r_max
      write(n_out,'(''  n_cheb_max   ='',i6)') n_cheb_max
      write(n_out,'(''  n_phi_max    ='',i6, &
           &   '' = no of azimuthal grid points'')') n_phi_max
      write(n_out,'(''  m_max        ='',i6, '' = max order'')') m_max
      write(n_out,'(''  n_m_max      ='',i6, '' = number of m s'')') n_m_max
      write(n_out,'(''  minc         ='',i6, '' = longitude symmetry wave no'')') minc
      write(n_out,*) ''

   end subroutine write_info
!---------------------------------------------------------------------------

end module pre_calculations
