module checkpoints

   use parallel_mod
   use precision_mod
   use constants, only: zero, two
   use char_manip, only: dble2str
   use blocking, only: nMstart,nMstop
   use communications, only: gather_from_mloc_to_rank0, &
       &                     scatter_from_rank0_to_mloc
   use truncation, only: n_r_max, m_max, minc, n_m_max, idx2m
   use namelists, only: ra,raxi,pr,sc,ek,radratio,alph1,alph2,tag, l_AB1, &
       &                start_file, scale_u, scale_t, l_heat, l_chem
   use radial_scheme, only: type_rscheme
   use radial_functions, only: rscheme, r
   use chebyshev, only: type_cheb
   use useful, only: abortRun, polynomial_interpolation

   implicit none

   private

   public :: write_checkpoint, read_checkpoint

   class(type_rscheme), pointer :: rscheme_old

contains

   subroutine write_checkpoint(time, dt, dtNew, n_time_step, n_log_file, &
              &                l_stop_time, t_Mloc, us_Mloc, up_Mloc,    &
              &                dtempdtLast_Mloc, dpsidtLast_Mloc)


      !-- Input variables
      real(cp),    intent(in) :: time
      real(cp),    intent(in) :: dt
      real(cp),    intent(in) :: dtNew
      integer,     intent(in) :: n_time_step
      integer,     intent(in) :: n_log_file
      logical,     intent(in) :: l_stop_time
      complex(cp), intent(in) :: t_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: dtempdtLast_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: dpsidtLast_Mloc(nMstart:nMstop,n_r_max)

      !-- Local variables
      complex(cp), allocatable :: work(:,:)
      integer :: n_rst_file, version
      character(len=100) :: rst_file, string

      if ( l_stop_time ) then
         rst_file="checkpoint_end."//tag
      else
         call dble2str(time,string)
         rst_file='checkpoint_t='//trim(string)//'.'//tag
      end if

      version = 1

      if ( rank == 0 ) then
         open(newunit=n_rst_file, file=rst_file, status='unknown', &
         &    form='unformatted')

         !-- Write the header of the file
         write(n_rst_file) version
         write(n_rst_file) time,dt,dtNew
         write(n_rst_file) ra,pr,raxi,sc,ek,radratio
         write(n_rst_file) n_r_max,m_max,minc

         !-- Store radius and scheme version (FD or CHEB)
         if ( rscheme%version == 'cheb' ) then
            write(n_rst_file) rscheme%version, rscheme%n_max, &
            &                 rscheme%order_boundary, alph1, alph2
         ! else
            ! write(n_rst_file) rscheme%version, rscheme%order, &
            ! &                 rscheme%order_boundary, fd_stretch, fd_ratio
         end if

         write(n_rst_file) r

         !-- Write the number of following fields
         write(n_rst_file) l_heat, l_chem

      end if

      !-- Memory allocation of global arrays to write outputs
      ! allocate( work(n_m_max,n_r_max) )
      if ( rank == 0 ) then
         allocate( work(n_m_max,n_r_max) )
      else
         allocate( work(1,1) )
      end if

      !-- Gather fields on rank 0 and write

      !-- us and uphi
      call gather_from_mloc_to_rank0(us_Mloc, work)
      if ( rank == 0 ) write(n_rst_file) work
      call gather_from_mloc_to_rank0(up_Mloc, work)
      if ( rank == 0 ) write(n_rst_file) work
      call gather_from_mloc_to_rank0(dpsidtLast_Mloc, work)
      if ( rank == 0 ) write(n_rst_file) work

      !-- Temperature
      if ( l_heat ) then
         call gather_from_mloc_to_rank0(t_Mloc, work)
         if ( rank == 0 ) write(n_rst_file) work
         call gather_from_mloc_to_rank0(dtempdtLast_Mloc, work)
         if ( rank == 0 ) write(n_rst_file) work
      end if

      deallocate( work )

      !-- Close checkpoint file and display a message in the log file
      if ( rank == 0 ) then

         close(n_rst_file)

         write(*,'(/,1P,A,/,A,ES20.10,/,A,I15,/,A,A)')&
         &    " ! Storing checkpoint file:",          &
         &    "             at time=",time,           &
         &    "            step no.=",n_time_step,    &
         &    "           into file=",rst_file

         write(n_log_file,'(/,1P,A,/,A,ES20.10,/,A,I15,/,A,A)') &
         &    " ! Storing checkpoint file:",                    &
         &    "             at time=",time,                     &
         &    "            step no.=",n_time_step,              &
         &    "           into file=",rst_file

      end if

   end subroutine write_checkpoint
!------------------------------------------------------------------------------
   subroutine read_checkpoint(us_Mloc, up_Mloc, dpsidtLast_Mloc, temp_Mloc, &
              &               dtempdtLast_Mloc, time, dt_old, dt_new)

      !-- Output variables
      complex(cp), intent(out) :: us_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(out) :: up_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(out) :: temp_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(out) :: dpsidtLast_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(out) :: dtempdtLast_Mloc(nMstart:nMstop, n_r_max)
      real(cp),    intent(out) :: time
      real(cp),    intent(out) :: dt_old
      real(cp),    intent(out) :: dt_new

      !-- Local variables
      logical :: startfile_does_exist, l_heat_old, l_chem_old
      integer :: n_start_file, version
      integer,     allocatable :: m2idx_old(:)
      real(cp),    allocatable :: r_old(:)
      complex(cp), allocatable :: work(:,:), work_old(:,:)
      real(cp) :: ra_old, raxi_old, sc_old, pr_old, radratio_old, ek_old
      integer :: n_r_max_old, m_max_old, minc_old, n_m_max_old
      character(len=72) :: rscheme_version_old
      real(cp) :: ratio1, ratio2
      integer :: n_in, n_in_2, m, n_m, n_r_max_max, m_max_max

      if ( rank == 0 ) then
         inquire(file=start_file, exist=startfile_does_exist)

         if ( startfile_does_exist ) then
            open(newunit=n_start_file, file=start_file, status='old', &
            &    form='unformatted')
         else
            call abortRun('! The restart file does not exist !')
         end if

         read(n_start_file) version
         read(n_start_file) time, dt_old, dt_new
         read(n_start_file) ra_old,pr_old,raxi_old,sc_old,ek_old,radratio_old
         read(n_start_file) n_r_max_old,m_max_old,minc_old

         n_r_max_max = max(n_r_max,n_r_max_old)
         m_max_max = max(m_max,m_max_old)

         !---- Compare parameters:
         if ( ra /= ra_old ) &
            write(*,'(/,'' ! New Rayleigh number (old/new):'',2ES16.6)') ra_old,ra
         if ( ek /= ek_old ) &
            write(*,'(/,'' ! New Ekman number (old/new):'',2ES16.6)') ek_old,ek
         if ( pr /= pr_old ) &
            write(*,'(/,'' ! New Prandtl number (old/new):'',2ES16.6)') pr_old,pr
         if ( raxi /= raxi_old ) &
            write(*,'(/,'' ! New composition-based Rayleigh number (old/new):'',2ES16.6)') raxi_old,raxi
         if ( sc /= sc_old ) &
            write(*,'(/,'' ! New Schmidt number (old/new):'',2ES16.6)') sc_old,sc
         if ( radratio /= radratio_old )                                    &
            write(*,'(/,'' ! New mag aspect ratio (old/new):'',2ES16.6)') &
            radratio_old,radratio

         if ( m_max_old /= m_max ) &
            write(*,*) '! New m_max (old,new)    :',m_max_old,m_max
         if ( minc_old /= minc ) &
            write(*,*) '! New minc (old,new)     :',minc_old,minc
         if ( n_r_max_old /= n_r_max ) &
            write(*,*) '! New n_r_max (old,new)  :',n_r_max_old,n_r_max

         read(n_start_file) rscheme_version_old, n_in, n_in_2, ratio1, ratio2
         if ( rscheme_version_old == 'cheb' ) then
            allocate ( type_cheb :: rscheme_old )
         ! else
            ! allocate ( type_fd :: rscheme_old )
         end if

         call rscheme_old%initialize(n_r_max_old, n_in, n_in_2, &
              &                      no_work_array=.true.)

         ! call rscheme_old%get_grid(n_r_max_old, r_icb_old, r_cmb_old, ratio1, &
              ! &                       ratio2, r_old)

         if ( rscheme%version /= rscheme_old%version ) &
            & write(*,'(/,'' ! New radial scheme (old/new):'',2A4)') &
            & rscheme_old%version, rscheme%version

         allocate( r_old(n_r_max_old) )
         read(n_start_file) r_old
         read(n_start_file) l_heat_old, l_chem_old

         n_m_max_old=m_max_old/minc_old+1
         allocate( m2idx_old(0:m_max_max) )

         m2idx_old(:)=-1
         n_m = 1
         do m=0,m_max_old,minc_old
            m2idx_old(m)=n_m
            n_m = n_m+1
         end do

         allocate( work_old(n_m_max_old, n_r_max_old) )
         allocate(     work(n_m_max, n_r_max) )

      else
         allocate( r_old(1), work_old(1,1), work(1,1), m2idx_old(1) )

      end if

      call MPI_Bcast(time,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(dt_old,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(dt_new,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ek_old,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(l_heat_old,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(l_chem_old,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

      if ( ek_old /= ek ) then ! If Ekman is different, let's use AB1
         l_AB1 = .true.
      end if

      !-- Read fields with rank0 and scatter them

      !-- us
      if ( rank == 0 ) then
         read( n_start_file ) work_old
         call map_field(work_old, work, r_old, m2idx_old, scale_u, &
              &         n_m_max_old, n_r_max_old, n_r_max_max,.false.)
      end if
      call scatter_from_rank0_to_mloc(work, us_Mloc)

      !-- uphi
      if ( rank == 0 ) then
         read( n_start_file ) work_old
         call map_field(work_old, work, r_old, m2idx_old, scale_u, &
              &         n_m_max_old, n_r_max_old, n_r_max_max,.false.)
      end if
      call scatter_from_rank0_to_mloc(work, up_Mloc)

      !-- dpsidtLast
      if ( rank == 0 ) then
         read( n_start_file ) work_old
         call map_field(work_old, work, r_old, m2idx_old, scale_u, &
              &         n_m_max_old, n_r_max_old, n_r_max_max,.true.)
      end if
      call scatter_from_rank0_to_mloc(work, dpsidtLast_Mloc)

      if ( l_heat_old ) then
         !-- Temperature
         if ( rank == 0 ) then
            read( n_start_file ) work_old
            call map_field(work_old, work, r_old, m2idx_old, scale_t, &
                 &         n_m_max_old, n_r_max_old, n_r_max_max,.false.)
         end if
         call scatter_from_rank0_to_mloc(work, temp_Mloc)

         !-- dTdtLast
         if ( rank == 0 ) then
            read( n_start_file ) work_old
            call map_field(work_old, work, r_old, m2idx_old, scale_t, &
                 &         n_m_max_old, n_r_max_old, n_r_max_max,.true.)
         end if
         call scatter_from_rank0_to_mloc(work, dtempdtLast_Mloc)
      end if

      if ( rank == 0 ) then
         call rscheme_old%finalize(no_work_array=.true.)
      end if
      deallocate( r_old, work_old, work, m2idx_old )

   end subroutine read_checkpoint
!------------------------------------------------------------------------------
   subroutine map_field(field_old, field_new, r_old, m2idx_old, scale_field, &
              &         n_m_max_old, n_r_max_old, n_r_max_max, lBc)

      !-- Input variables:
      integer,     intent(in) :: n_r_max_old
      integer,     intent(in) :: n_m_max_old
      integer,     intent(in) :: n_r_max_max
      complex(cp), intent(in) :: field_old(n_m_max_old, n_r_max_old)
      real(cp),    intent(in) :: r_old(n_r_max_old)
      integer,     intent(in) :: m2idx_old(0:)
      real(cp),    intent(in) :: scale_field
      logical,     intent(in) :: lBc

      !-- Output variable:
      complex(cp), intent(out) :: field_new(n_m_max, n_r_max)

      !-- Local variables:
      complex(cp) :: radial_data(n_r_max_max)
      integer :: n_m, m, n_m_old, n_r

      do n_m=1,n_m_max

         m = idx2m(n_m)
         n_m_old = m2idx_old(m)

         if ( n_m_old > 0 ) then

            if ( n_r_max /= n_r_max_old .or.                               &
            &    rscheme%order_boundary /= rscheme_old%order_boundary .or. &
            &    rscheme%version /= rscheme_old%version ) then
               do n_r=1,n_r_max_old
                  radial_data(n_r)=field_old(n_m_old,n_r)
               end do
               call map_field_r(radial_data, r_old, n_r_max_old,n_r_max_max,lBc)
               do n_r=1,n_r_max
                  field_new(n_m, n_r) = scale_field*radial_data(n_r)
               end do
            else
               do n_r=1,n_r_max
                  field_new(n_m, n_r) = scale_field*field_old(n_m_old,n_r)
               end do
            end if
         else
            do n_r=1,n_r_max
               field_new(n_m,n_r)=zero
            end do
         end if

      end do

   end subroutine map_field
!------------------------------------------------------------------------------
   subroutine map_field_r(radial_data,r_old,n_r_max_old,n_r_max_max,lBc)

      !-- Input variables
      integer,  intent(in) :: n_r_max_old
      integer,  intent(in) :: n_r_max_max
      real(cp), intent(in) :: r_old(n_r_max_old)
      logical,  intent(in) :: lBc

      !-- Output data
      complex(cp), intent(inout) :: radial_data(:)

      !-- Local variables
      real(cp) :: radial_data_real(n_r_max_max),radial_data_imag(n_r_max_max)
      integer :: n_r, n_r_old, n_r_index_start
      real(cp) :: xold(4)
      complex(cp) :: yold(4)
      complex(cp), allocatable :: work(:)
      real(cp) :: cheb_norm_old,cheb_fac

      !-- Since beta is singular on the first grid point
      !-- it is not possible to take a dct of the dpsidtLast array
      !-- Hence we here overwrite the boundary value by a simple
      !-- linear interpolation. This is harmless for the rest since
      !-- the boundary value is useless in the time advance
      if ( lBc ) then
         radial_data(1)=two*radial_data(2)-radial_data(3)
      end if

      !-- If **both** the old and the new schemes are Chebyshev, we can
      !-- use costf to get the new data
      !-- Both have also to use the same mapping (order_boundary is a proxy of l_map
      if ( rscheme%version == 'cheb' .and. rscheme_old%version == 'cheb' &
      &   .and. rscheme%order_boundary == rscheme_old%order_boundary  ) then

         do n_r=1,n_r_max_old
            radial_data_real(n_r) = real(radial_data(n_r))
            radial_data_imag(n_r) = aimag(radial_data(n_r))
         end do
         call rscheme_old%costf1(radial_data_real,n_r_max_old)
         call rscheme_old%costf1(radial_data_imag,n_r_max_old)

         !----- Fill up cheb polynomial with zeros:
         if ( n_r_max>n_r_max_old ) then
            n_r_index_start=n_r_max_old+1
            do n_r=n_r_index_start,n_r_max
               radial_data_real(n_r)=0.0_cp
               radial_data_imag(n_r)=0.0_cp
            end do
         end if

         !----- Now transform to new radial grid points:
         call rscheme%costf1(radial_data_real,n_r_max)
         call rscheme%costf1(radial_data_imag,n_r_max)
         !----- Rescale :
         cheb_norm_old=sqrt(two/real(n_r_max_old-1,kind=cp))
         cheb_fac=cheb_norm_old/rscheme%rnorm
         do n_r=1,n_r_max
            radial_data(n_r)=cheb_fac*cmplx(radial_data_real(n_r), &
            &                               radial_data_imag(n_r),cp)
         end do

      !-- If either the old grid or the new grid is FD, we use a 
      !-- polynomial interpolation
      else

         allocate( work(n_r_max) )

         !-- Interpolate data and store into a work array
         do n_r=1,n_r_max

            n_r_old=minloc(abs(r_old-r(n_r)),1)
            if ( n_r_old < 3 ) n_r_old=3
            if ( n_r_old == n_r_max_old ) n_r_old=n_r_max_old-1

            xold(1)=r_old(n_r_old-2)
            xold(2)=r_old(n_r_old-1)
            xold(3)=r_old(n_r_old)
            xold(4)=r_old(n_r_old+1)

            yold(1)=radial_data(n_r_old-2)
            yold(2)=radial_data(n_r_old-1)
            yold(3)=radial_data(n_r_old)
            yold(4)=radial_data(n_r_old+1)

            call polynomial_interpolation(xold, yold, r(n_r), work(n_r))

         end do

         !-- Copy interpolated data
         do n_r=1,n_r_max
            radial_data(n_r)=work(n_r)
         end do

         deallocate( work )

      end if

   end subroutine map_field_r
!------------------------------------------------------------------------------
end module checkpoints
