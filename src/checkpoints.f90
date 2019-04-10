module checkpoints

   use parallel_mod
   use precision_mod
   use communications, only: lm2r_fields, transp_lm2r
   use constants, only: zero, two
   use char_manip, only: dble2str
   use fields, only: work_3D_Rloc
   use blocking, only: nMstart, nMstop, nm_per_rank, nRstart3D, nRstop3D, &
       &               nR_per_rank_3D, lmStart, lmStop
   use blocking_lm, only: st_map
   use communications, only: gather_from_mloc_to_rank0, lm2r_fields,  &
       &                     scatter_from_rank0_to_mloc, transp_lm2r, &
       &                     scatter_from_rank0_to_lmloc
   use truncation, only: n_r_max, m_max, minc, n_m_max, idx2m
   use truncation_3D, only: lm_max, l_max, m_max_3D, n_phi_tot_3D, n_theta_max, &
       &                    n_r_max_3D, minc_3D
   use namelists, only: ra,raxi,pr,sc,ek,radratio,alph1,alph2,tag, l_AB1, &
       &                start_file, scale_u, scale_t, l_heat, l_chem,     &
       &                l_bridge_step, l_cheb_coll, scale_xi, l_3D,       &
       &                l_heat_3D, l_mag
   use radial_scheme, only: type_rscheme
   use radial_functions, only: rscheme, r, rscheme_3D, r_3D
   use chebyshev, only: type_cheb
   use useful, only: abortRun, polynomial_interpolation
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray

   implicit none

   private

   public :: read_checkpoint, write_checkpoint

contains

   subroutine write_checkpoint(time, tscheme, n_time_step, n_log_file,     &
              &                l_stop_time, t_Mloc, xi_Mloc, us_Mloc,      &
              &                up_Mloc, dTdt, dxidt, dpsidt, temp_3D_Rloc, &
              &                dTdt_3D)
      !
      ! This subroutine writes the checkpoint files using MPI-IO. For the sake
      ! of simplicity we do not include the record marker. Classical Fortran can
      ! hence only read it when using access='stream'
      !


      !-- Input variables
      real(cp),            intent(in) :: time
      class(type_tscheme), intent(in) :: tscheme
      integer,             intent(in) :: n_time_step
      integer,             intent(in) :: n_log_file
      logical,             intent(in) :: l_stop_time
      complex(cp),         intent(in) :: t_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(in) :: xi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(in) :: temp_3D_Rloc(lm_max,nRstart3D:nRstop3D)
      type(type_tarray),   intent(in) :: dTdt
      type(type_tarray),   intent(in) :: dxidt
      type(type_tarray),   intent(in) :: dpsidt
      type(type_tarray),   intent(in) :: dTdt_3D

      !-- Local variables
      integer :: info, fh, filetype, n_o, lmtype
      integer :: version
      integer :: istat(MPI_STATUS_SIZE)
      integer :: arr_size(2), arr_loc_size(2), arr_start(2)
      integer(kind=MPI_OFFSET_KIND) :: disp, offset
      character(len=100) :: rst_file, string

      if ( l_stop_time ) then
         rst_file="checkpoint_end."//tag
      else
         call dble2str(time,string)
         rst_file='checkpoint_t='//trim(string)//'.'//tag
      end if

      version = 6

      call MPI_Info_create(info, ierr)

      !-- Enable collective buffering
      call MPI_Info_set(info, "romio_cb_write", "automatic", ierr)
      call MPI_Info_set(info, "romio_cb_read", "automatic", ierr)

      !-- Disable data sieving (let the filesystem handles it)
      call MPI_Info_set(info, "romio_ds_write", "disable", ierr)
      call MPI_Info_set(info, "romio_ds_read", "disable", ierr)

      !-- Set the stripping unit to 4M
      call MPI_Info_set(info, "stripping_unit", "4194304", ierr)

      !-- Set the buffer size to 4M
      call MPI_Info_set(info,"cb_buffer_size","4194304", ierr)

      !-- Open file
      call MPI_File_Open(MPI_COMM_WORLD, rst_file, ior(MPI_MODE_WRONLY, &
           &             MPI_MODE_CREATE), info, fh, ierr)

      !-- Only rank=0 writes the header of the file
      if ( rank == 0 ) then
         call MPI_File_Write(fh, version, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, time, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, l_cheb_coll, 1, MPI_LOGICAL, istat, ierr)
         call MPI_File_Write(fh, tscheme%family, len(tscheme%family),  &
              &              MPI_CHARACTER, istat, ierr)
         call MPI_File_Write(fh, tscheme%norder_exp, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, tscheme%norder_imp_lin, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, tscheme%norder_imp, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, tscheme%dt, size(tscheme%dt), MPI_DEF_REAL, &
              &              istat, ierr)
         call MPI_File_Write(fh, ra, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, pr, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, raxi, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, sc, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, ek, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, radratio, 1, MPI_DEF_REAL, istat, ierr)

         call MPI_File_Write(fh, n_r_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, m_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, minc, 1, MPI_INTEGER, istat, ierr)

         call MPI_File_Write(fh, rscheme%version, len(rscheme%version),  &
              &              MPI_CHARACTER, istat, ierr)
         call MPI_File_Write(fh, rscheme%n_max, 1, MPI_INTEGER, istat, ierr)
         call MPI_File_Write(fh, rscheme%order_boundary, 1, MPI_INTEGER, istat, &
              &              ierr)
         call MPI_File_Write(fh, alph1, 1, MPI_DEF_REAL, istat, ierr)
         call MPI_File_Write(fh, alph2, 1, MPI_DEF_REAL, istat, ierr)

         call MPI_File_Write(fh, r, n_r_max, MPI_DEF_REAL, istat, ierr)

         call MPI_File_Write(fh, l_heat, 1, MPI_LOGICAL, istat, ierr)
         call MPI_File_Write(fh, l_chem, 1, MPI_LOGICAL, istat, ierr)

         call MPI_File_Write(fh, l_3D, 1, MPI_LOGICAL, istat, ierr)

         if ( l_3D ) then
            call MPI_File_Write(fh, n_r_max_3D, 1, MPI_INTEGER, istat, ierr)
            call MPI_File_Write(fh, n_theta_max, 1, MPI_INTEGER, istat, ierr)
            call MPI_File_Write(fh, n_phi_tot_3D, 1, MPI_INTEGER, istat, ierr)
            call MPI_File_Write(fh, l_max, 1, MPI_INTEGER, istat, ierr)
            call MPI_File_Write(fh, m_max_3D, 1, MPI_INTEGER, istat, ierr)
            call MPI_File_Write(fh, minc_3D, 1, MPI_INTEGER, istat, ierr)
            call MPI_File_Write(fh, lm_max, 1, MPI_INTEGER, istat, ierr)

            call MPI_File_Write(fh, rscheme_3D%version, len(rscheme_3D%version), &
                 &              MPI_CHARACTER, istat, ierr)
            call MPI_File_Write(fh, rscheme_3D%n_max, 1, MPI_INTEGER, istat, ierr)
            call MPI_File_Write(fh, rscheme_3D%order_boundary, 1, MPI_INTEGER, &
                 &              istat, ierr)
            !-- Mapping for the 3-D grid: disabled for now -- !
            call MPI_File_Write(fh, 1.0_cp, 1, MPI_DEF_REAL, istat, ierr)
            call MPI_File_Write(fh, 0.0_cp, 1, MPI_DEF_REAL, istat, ierr)
            call MPI_File_Write(fh, r_3D, n_r_max_3D, MPI_DEF_REAL, istat, ierr)
            call MPI_File_Write(fh, l_heat_3D, 1, MPI_LOGICAL, istat, ierr)
            call MPI_File_Write(fh, l_mag, 1, MPI_LOGICAL, istat, ierr)
         end if

      end if

      arr_size(1) = n_m_max
      arr_size(2) = n_r_max
      arr_loc_size(1) = nm_per_rank
      arr_loc_size(2) = n_r_max
      arr_start(1) = nMstart-1
      arr_start(2) = 0
      call MPI_Type_Create_Subarray(2, arr_size, arr_loc_size, arr_start, &
           &                        MPI_ORDER_FORTRAN, MPI_DEF_COMPLEX,   &
           &                        filetype, ierr)
      call MPI_Type_Commit(filetype, ierr)

      if ( rank == 0 ) then
         !-- Rank 0 gets the displacement
         call MPI_File_get_position(fh, offset, ierr)
         call MPI_File_get_byte_offset(fh, offset, disp, ierr)
      end if
      !-- Broadcast the displacement
      call MPI_Bcast(disp, 1, MPI_OFFSET, 0, MPI_COMM_WORLD, ierr)

      !-- Set the view after the header
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype, "native", &
           &                 info, ierr)

      !-- Now finally write the fields
      call MPI_File_Write_all(fh, us_Mloc, nm_per_rank*n_r_max, &
           &                  MPI_DEF_COMPLEX, istat, ierr)
      disp = disp+n_r_max*n_m_max*SIZEOF_DEF_COMPLEX
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype, "native", &
           &                 info, ierr)
      call MPI_File_Write_all(fh, up_Mloc, nm_per_rank*n_r_max, &
           &                  MPI_DEF_COMPLEX, istat, ierr)
      disp = disp+n_r_max*n_m_max*SIZEOF_DEF_COMPLEX
      call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype, "native", &
           &                 info, ierr)

      if ( tscheme%family == 'MULTISTEP' ) then
         do n_o=2,tscheme%norder_exp
            call MPI_File_Write_all(fh, dpsidt%expl(:,:,n_o), nm_per_rank*n_r_max, &
                 &                  MPI_DEF_COMPLEX, istat, ierr)
            disp = disp+n_r_max*n_m_max*SIZEOF_DEF_COMPLEX
            call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype, "native", &
                 &                 info, ierr)
         end do
         do n_o=2,tscheme%norder_imp_lin-1
            call MPI_File_Write_all(fh, dpsidt%impl(:,:,n_o), nm_per_rank*n_r_max, &
                 &                  MPI_DEF_COMPLEX, istat, ierr)
            disp = disp+n_r_max*n_m_max*SIZEOF_DEF_COMPLEX
            call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype, "native", &
                 &                 info, ierr)
         end do
         do n_o=2,tscheme%norder_imp-1
            call MPI_File_Write_all(fh, dpsidt%old(:,:,n_o), nm_per_rank*n_r_max, &
                 &                  MPI_DEF_COMPLEX, istat, ierr)
            disp = disp+n_r_max*n_m_max*SIZEOF_DEF_COMPLEX
            call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype, "native", &
                 &                 info, ierr)
         end do
      end if

      !-- Temperature
      if ( l_heat ) then
         call MPI_File_Write_all(fh, t_Mloc, nm_per_rank*n_r_max, MPI_DEF_COMPLEX, &
              &                  istat, ierr)
         disp = disp+n_r_max*n_m_max*SIZEOF_DEF_COMPLEX
         call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype, "native", &
              &                 info, ierr)

         if ( tscheme%family == 'MULTISTEP' ) then
            do n_o=2,tscheme%norder_exp
               call MPI_File_Write_all(fh, dTdt%expl(:,:,n_o), nm_per_rank*n_r_max, &
                    &                  MPI_DEF_COMPLEX, istat, ierr)
               disp = disp+n_r_max*n_m_max*SIZEOF_DEF_COMPLEX
               call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype, &
                    &                 "native", info, ierr)
            end do
            do n_o=2,tscheme%norder_imp_lin-1
               call MPI_File_Write_all(fh, dTdt%impl(:,:,n_o), nm_per_rank*n_r_max, &
                    &                  MPI_DEF_COMPLEX, istat, ierr)
               disp = disp+n_r_max*n_m_max*SIZEOF_DEF_COMPLEX
               call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype, &
                    &                 "native", info, ierr)
            end do
            do n_o=2,tscheme%norder_imp-1
               call MPI_File_Write_all(fh, dTdt%old(:,:,n_o), nm_per_rank*n_r_max, &
                    &                  MPI_DEF_COMPLEX, istat, ierr)
               disp = disp+n_r_max*n_m_max*SIZEOF_DEF_COMPLEX
               call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype, &
                    &                 "native", info, ierr)
            end do
         end if

      end if

      !-- Chemical composition
      if ( l_chem ) then
         call MPI_File_Write_all(fh, xi_Mloc, nm_per_rank*n_r_max, MPI_DEF_COMPLEX, &
              &                  istat, ierr)
         disp = disp+n_r_max*n_m_max*SIZEOF_DEF_COMPLEX
         call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype, "native", &
              &                 info, ierr)

         if ( tscheme%family == 'MULTISTEP' ) then
            do n_o=2,tscheme%norder_exp
               call MPI_File_Write_all(fh, dxidt%expl(:,:,n_o), nm_per_rank*n_r_max, &
                    &                  MPI_DEF_COMPLEX, istat, ierr)
               disp = disp+n_r_max*n_m_max*SIZEOF_DEF_COMPLEX
               call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype, &
                    &                 "native", info, ierr)
            end do
            do n_o=2,tscheme%norder_imp_lin-1
               call MPI_File_Write_all(fh, dxidt%impl(:,:,n_o), nm_per_rank*n_r_max, &
                    &                  MPI_DEF_COMPLEX, istat, ierr)
               disp = disp+n_r_max*n_m_max*SIZEOF_DEF_COMPLEX
               call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype,  &
                    &                 "native", info, ierr)
            end do
            do n_o=2,tscheme%norder_imp-1
               call MPI_File_Write_all(fh, dxidt%old(:,:,n_o), nm_per_rank*n_r_max, &
                    &                  MPI_DEF_COMPLEX, istat, ierr)
               disp = disp+n_r_max*n_m_max*SIZEOF_DEF_COMPLEX
               call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, filetype, &
                    &                 "native", info, ierr)
            end do
         end if

      end if

      !-- Write 3-D temperature
      if ( l_heat_3D ) then
         arr_size(1) = lm_max
         arr_size(2) = n_r_max_3D
         arr_loc_size(1) = lm_max
         arr_loc_size(2) = nR_per_rank_3D
         arr_start(1) = 0
         arr_start(2) = nRStart3D-1
         call MPI_Type_Create_Subarray(2,arr_size, arr_loc_size, arr_start, &
              &                        MPI_ORDER_FORTRAN, MPI_DEF_COMPLEX,  &
              &                        lmtype, ierr)
         call MPI_Type_Commit(lmtype, ierr)

         !-- Set the view after the QG fields
         call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, &
              &                 lmtype, "native", MPI_INFO_NULL, ierr)
         call MPI_File_Write_all(fh, temp_3D_Rloc, lm_max*nR_per_rank_3D,  &
              &                  MPI_DEF_COMPLEX, istat, ierr)

         if ( tscheme%family == 'MULTISTEP' ) then
            do n_o=2,tscheme%norder_exp
               call transp_lm2r(lm2r_fields, dTdt_3D%expl(:,:,n_o), &
                    &           work_3D_Rloc)
               disp = disp+n_r_max_3D*lm_max*SIZEOF_DEF_COMPLEX
               call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, lmtype, "native", &
                 &                  MPI_INFO_NULL, ierr)
               call MPI_File_Write_all(fh, work_3D_Rloc, lm_max*nR_per_rank_3D, &
                    &              MPI_DEF_COMPLEX, istat, ierr)
            end do
            do n_o=2,tscheme%norder_imp_lin-1
               call transp_lm2r(lm2r_fields, dTdt_3D%impl(:,:,n_o), &
                    &           work_3D_Rloc)
               disp = disp+n_r_max_3D*lm_max*SIZEOF_DEF_COMPLEX
               call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, lmtype, "native", &
                    &                  MPI_INFO_NULL, ierr)
               call MPI_File_Write_all(fh, work_3D_Rloc, lm_max*nR_per_rank_3D, &
                    &                  MPI_DEF_COMPLEX, istat, ierr)
            end do
            do n_o=2,tscheme%norder_imp-1
               call transp_lm2r(lm2r_fields, dTdt_3D%old(:,:,n_o), &
                    &           work_3D_Rloc)
               disp = disp+n_r_max_3D*lm_max*SIZEOF_DEF_COMPLEX
               call MPI_File_Set_View(fh, disp, MPI_DEF_COMPLEX, lmtype, "native", &
                    &                  MPI_INFO_NULL, ierr)
               call MPI_File_Write_all(fh, work_3D_Rloc, lm_max*nR_per_rank_3D, &
                    &                  MPI_DEF_COMPLEX, istat, ierr)
            end do
         end if

         call MPI_Type_Free(lmtype, ierr)
      end if

      call MPI_Type_Free(filetype, ierr)
      call MPI_Info_free(info, ierr)
      call MPI_File_close(fh, ierr)

      !-- Close checkpoint file and display a message in the log file
      if ( rank == 0 ) then

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
   subroutine read_checkpoint(us_Mloc, up_Mloc, temp_Mloc, xi_Mloc, dpsidt, &
              &               dTdt, dxidt, temp_3D_LMloc, dTdt_3D, time, tscheme)

      !-- Output variables
      class(type_tscheme), intent(inout) :: tscheme
      complex(cp),         intent(out) :: us_Mloc(nMstart:nMstop, n_r_max)
      complex(cp),         intent(out) :: up_Mloc(nMstart:nMstop, n_r_max)
      complex(cp),         intent(out) :: temp_Mloc(nMstart:nMstop, n_r_max)
      complex(cp),         intent(out) :: xi_Mloc(nMstart:nMstop, n_r_max)
      complex(cp),         intent(out) :: temp_3D_LMloc(lmStart:lmStop,n_r_max_3D)
      type(type_tarray),   intent(inout) :: dpsidt
      type(type_tarray),   intent(inout) :: dTdt
      type(type_tarray),   intent(inout) :: dxidt
      type(type_tarray),   intent(inout) :: dTdt_3D
      real(cp),            intent(out) :: time

      !-- Local variables
      class(type_rscheme), pointer :: rscheme_old
      class(type_rscheme), pointer :: rscheme_3D_old
      logical :: startfile_does_exist, l_heat_old, l_chem_old
      integer :: n_start_file, version
      integer,     allocatable :: m2idx_old(:)
      real(cp),    allocatable :: r_old(:), r_3D_old(:)
      complex(cp), allocatable :: work(:,:), work_old(:,:)
      real(cp) :: ra_old, raxi_old, sc_old, pr_old, radratio_old, ek_old
      integer :: n_r_max_old, m_max_old, minc_old, n_m_max_old
      character(len=72) :: rscheme_version_old
      character(len=10) :: tscheme_family_old
      real(cp) :: ratio1, ratio2
      integer :: n_in, n_in_2, m, n_m, n_r_max_max, m_max_max
      integer :: norder_imp_lin_old, norder_exp_old, n_o, norder_imp_old
      logical :: l_coll_old, l_3D_old, l_heat_3D_old, l_mag_old
      real(cp), allocatable :: dt_array_old(:)
      integer, allocatable :: lm2lmo(:)
      integer :: n_r_max_3D_old, n_theta_max_old, n_phi_tot_old
      integer :: n_r_max_3D_max
      integer :: l_max_old, m_max_3D_old, minc_3D_old, lm_max_old

      if ( rank == 0 ) then
         inquire(file=start_file, exist=startfile_does_exist)

         if ( startfile_does_exist ) then
            open(newunit=n_start_file, file=start_file, status='old', &
            &    form='unformatted', access='stream')
         else
            call abortRun('! The restart file does not exist !')
         end if

         read(n_start_file) version
         if ( version == 1 ) then ! This was CN/AB2 in the initial version
            allocate( dt_array_old(max(2,tscheme%norder_exp)) )
            dt_array_old(:)=0.0_cp
            read(n_start_file) time, dt_array_old(2), dt_array_old(1)
            norder_imp_lin_old = 2
            norder_imp_old = 2
            norder_exp_old = 2
            l_coll_old = .true.
            tscheme_family_old = 'MULTISTEP'
         else if ( version == 2 ) then ! This was without l_cheb_coll
            read(n_start_file) time
            read(n_start_file) norder_exp_old, norder_imp_old
            norder_imp_lin_old = norder_imp_old
            allocate( dt_array_old(max(norder_exp_old,tscheme%norder_exp) ) )
            dt_array_old(:)=0.0_cp
            read(n_start_file) dt_array_old(1:norder_exp_old)
            l_coll_old = .true.
            tscheme_family_old = 'MULTISTEP'
         else if ( version == 3 ) then
            read(n_start_file) time
            read(n_start_file) l_coll_old
            read(n_start_file) norder_exp_old, norder_imp_old
            norder_imp_lin_old=norder_imp_old
            allocate( dt_array_old(max(norder_exp_old,tscheme%norder_exp) ) )
            dt_array_old(:)=0.0_cp
            read(n_start_file) dt_array_old(1:norder_exp_old)
            tscheme_family_old = 'MULTISTEP'
         else if ( version == 4 ) then
            read(n_start_file) time
            read(n_start_file) l_coll_old
            read(n_start_file) norder_exp_old, norder_imp_lin_old, norder_imp_old
            allocate( dt_array_old(max(norder_exp_old,tscheme%norder_exp) ) )
            dt_array_old(:)=0.0_cp
            read(n_start_file) dt_array_old(1:norder_exp_old)
            tscheme_family_old = 'MULTISTEP'
         else if ( version == 5 .or. version == 6 ) then
            read(n_start_file) time
            read(n_start_file) l_coll_old
            read(n_start_file) tscheme_family_old
            read(n_start_file) norder_exp_old, norder_imp_lin_old, norder_imp_old
            if ( tscheme_family_old == 'MULTISTEP' ) then
               allocate( dt_array_old(max(norder_exp_old,tscheme%norder_exp) ) )
               dt_array_old(:)=0.0_cp
               read(n_start_file) dt_array_old(1:norder_exp_old)
            else if ( tscheme_family_old == 'DIRK' ) then
               allocate( dt_array_old(max(1,size(tscheme%dt))) )
               dt_array_old(:)=0.0_cp
               read(n_start_file) dt_array_old(1)
            end if
         end if
         if ( tscheme_family_old == 'MULTISTEP' ) then
            dt_array_old(norder_exp_old:size(tscheme%dt))=dt_array_old(norder_exp_old)
         else if ( tscheme_family_old == 'DIRK' ) then
            dt_array_old(1:size(tscheme%dt))=dt_array_old(1)
         end if
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

         call rscheme_old%initialize(nMstart, nMstop, n_r_max_old, n_in, &
              &                      n_in_2,l_cheb_coll=.true.,          &
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

         l_heat_3D_old = .false.
         l_mag_old = .false.
         if ( version == 6 ) then
            read(n_start_file) l_3D_old
         else
            l_3D_old = .false.
         end if

         if ( l_3D_old ) then
            read(n_start_file) n_r_max_3D_old, n_theta_max_old, n_phi_tot_old
            n_r_max_3D_max = max(n_r_max_3D_old, n_r_max_3D)
            read(n_start_file) l_max_old, m_max_3D_old, minc_3D_old, lm_max_old
            read(n_start_file) rscheme_version_old, n_in, n_in_2, ratio1, ratio2
            if ( rscheme_version_old == 'cheb' ) then
               allocate ( type_cheb :: rscheme_3D_old )
            end if

            call rscheme_3D_old%initialize(lmStart, lmStop, n_r_max_3D_old, n_in, &
                 &                         n_in_2,l_cheb_coll=.true.,             &
                 &                         no_work_array=.true.)

            if ( rscheme_3D%version /= rscheme_3D_old%version ) &
               & write(*,'(/,'' ! New (3D) radial scheme (old/new):'',2A4)') &
               & rscheme_3D_old%version, rscheme_3D%version

            allocate( r_3D_old(n_r_max_3D_old) )
            read(n_start_file) r_3D_old
            read(n_start_file) l_heat_3D_old, l_mag_old
         end if

         if ( l_3D_old .and. l_3D ) then
            allocate( lm2lmo(lm_max) )
            call get_lm2lmo(lm2lmo, lm_max, l_max_old, minc_3D_old)
         else
            allocate( lm2lmo(1) )
         end if

      else
         allocate( r_old(1), work_old(1,1), work(1,1), m2idx_old(1) )
         allocate( r_3D_old(1), lm2lmo(1) )
      end if

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_Bcast(time,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(version,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(norder_exp_old,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(norder_imp_old,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(norder_imp_lin_old,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(tscheme_family_old,len(tscheme_family_old),MPI_CHARACTER,0, &
           &         MPI_COMM_WORLD,ierr)
      if ( tscheme_family_old == 'MULTISTEP' ) then
         if ( rank /= 0 ) allocate( dt_array_old(max(norder_exp_old,tscheme%norder_exp)) )
      else if ( tscheme_family_old == 'DIRK' ) then
         if ( rank /= 0 ) allocate( dt_array_old(max(1,size(tscheme%dt))) )
      end if
      call MPI_Bcast(dt_array_old,size(dt_array_old),MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ek_old,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(l_heat_old,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(l_chem_old,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(l_coll_old,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(l_heat_old,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(l_3D_old,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(l_heat_3D_old,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(l_mag_old,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

      !-- Fill the time step array
      do n_o=1,size(tscheme%dt)
         !-- If the new scheme has higher order one fill the missing dt values
         !-- with the oldest
         if ( n_o > size(dt_array_old) ) then
            tscheme%dt(n_o)=dt_array_old(size(dt_array_old))
         else
            tscheme%dt(n_o)=dt_array_old(n_o)
         end if
      end do

      !-- If new Ekman and old Ekman differ, we use AB1 for the first time step.
      if ( ek_old /= ek ) then 
         l_AB1 = .true.
      end if

      !-- If old and new schemes differ in precision, one has to use a bridging step
      if ( tscheme%family == 'MULTISTEP' ) then
         if ( tscheme_family_old == 'DIRK' ) then
            l_bridge_step = .true.
         else
            if ( tscheme%norder_imp > norder_imp_old .or. tscheme%norder_imp_lin > norder_imp_lin_old ) then
               l_bridge_step = .true.
            else
               l_bridge_step = .false.
            end if
         end if
      else if ( tscheme%family == 'DIRK' ) then
         l_bridge_step = .false.
      end if

      !-- If old and new radial schemes differ, one has to use AB1 for the first
      !-- time step and a bridging CN step for the implicit part.
      if ( tscheme%family == 'MULTISTEP' ) then
         if ( l_coll_old .neqv. l_cheb_coll ) then
            l_AB1 = .true.
            l_bridge_step = .true.
         end if
      end if

      !-- If version of the old file is < 4 we have to bridge the steps as well
      if ( tscheme%family == 'MULTISTEP' ) then
         if ( version < 4 .and. norder_imp_old > 2 ) then
            l_AB1 = .true.
            l_bridge_step = .true.
         end if
      end if

      !-- Read fields with rank0 and scatter them

      !-- us
      if ( rank == 0 ) then
         read( n_start_file ) work_old
         call map_field(work_old, work, rscheme_old, r_old, m2idx_old, scale_u, &
              &         n_m_max_old, n_r_max_old, n_r_max_max, lBc=.false.,     &
              &         l_phys_space=.true.)
      end if
      call scatter_from_rank0_to_mloc(work, us_Mloc)

      !-- uphi
      if ( rank == 0 ) then
         read( n_start_file ) work_old
         call map_field(work_old, work, rscheme_old, r_old, m2idx_old, scale_u, &
              &         n_m_max_old, n_r_max_old, n_r_max_max, lBc=.false.,     &
              &         l_phys_space=.true.)
      end if
      call scatter_from_rank0_to_mloc(work, up_Mloc)

      if ( tscheme_family_old == 'MULTISTEP' ) then
         !-- Explicit time step
         do n_o=2,norder_exp_old
            if ( rank == 0 ) then
               read( n_start_file ) work_old
               call map_field(work_old, work, rscheme_old, r_old, m2idx_old,  &
                    &         scale_u, n_m_max_old, n_r_max_old, n_r_max_max, &
                    &         lBc=.true., l_phys_space=l_coll_old)
            end if
            if ( n_o <= tscheme%norder_exp ) then
               call scatter_from_rank0_to_mloc(work, dpsidt%expl(:,:,n_o))
            end if
         end do

         !-- Implicit time step
         do n_o=2,norder_imp_lin_old-1
            if ( rank == 0 ) then
               read( n_start_file ) work_old
               call map_field(work_old, work, rscheme_old, r_old, m2idx_old,  &
                    &         scale_u, n_m_max_old, n_r_max_old, n_r_max_max, &
                    &         lBc=.true., l_phys_space=l_coll_old)
            end if
            if ( n_o <= tscheme%norder_imp_lin-1 ) then
               call scatter_from_rank0_to_mloc(work, dpsidt%impl(:,:,n_o))
            end if
         end do
         if ( version > 3 ) then
            do n_o=2,norder_imp_old-1
               if ( rank == 0 ) then
                  read( n_start_file ) work_old
                  call map_field(work_old, work, rscheme_old, r_old, m2idx_old,  &
                       &         scale_u, n_m_max_old, n_r_max_old, n_r_max_max, &
                       &         lBc=.false., l_phys_space=l_coll_old)
               end if
               if ( n_o <= tscheme%norder_imp-1 ) then
                  call scatter_from_rank0_to_mloc(work, dpsidt%old(:,:,n_o))
               end if
            end do
         end if
      end if

      if ( l_heat_old ) then
         !-- Temperature
         if ( rank == 0 ) then
            read( n_start_file ) work_old
            call map_field(work_old, work, rscheme_old, r_old, m2idx_old,  &
                 &         scale_t, n_m_max_old, n_r_max_old, n_r_max_max, &
                 &         lBc=.false., l_phys_space=.true.)
         end if
         call scatter_from_rank0_to_mloc(work, temp_Mloc)

         if ( tscheme_family_old == 'MULTISTEP' ) then
            !-- Explicit time step
            do n_o=2,norder_exp_old
               if ( rank == 0 ) then
                  read( n_start_file ) work_old
                  call map_field(work_old, work, rscheme_old, r_old, m2idx_old,  &
                       &         scale_t, n_m_max_old, n_r_max_old, n_r_max_max, &
                       &         lBc=.true., l_phys_space=l_coll_old)
               end if
               if ( n_o <= tscheme%norder_exp ) then
                  call scatter_from_rank0_to_mloc(work, dTdt%expl(:,:,n_o))
               end if
            end do

            !-- Implicit time step
            do n_o=2,norder_imp_lin_old-1
               if ( rank == 0 ) then
                  read( n_start_file ) work_old
                  call map_field(work_old, work, rscheme_old, r_old, m2idx_old,  &
                       &         scale_t, n_m_max_old, n_r_max_old, n_r_max_max, &
                       &         lBc=.true., l_phys_space=l_coll_old)
               end if
               if ( n_o <= tscheme%norder_imp_lin-1 ) then
                  call scatter_from_rank0_to_mloc(work, dTdt%impl(:,:,n_o))
               end if
            end do
            if ( version > 3 ) then
               do n_o=2,norder_imp_old-1
                  if ( rank == 0 ) then
                     read( n_start_file ) work_old
                     call map_field(work_old, work, rscheme_old, r_old, m2idx_old,&
                          &         scale_t, n_m_max_old, n_r_max_old,            &
                          &         n_r_max_max, lBc=.true.,                      &
                          &         l_phys_space=l_coll_old)
                  end if
                  if ( n_o <= tscheme%norder_imp-1 ) then
                     call scatter_from_rank0_to_mloc(work, dTdt%old(:,:,n_o))
                  end if
               end do
            end if
         end if

      end if

      if ( l_chem_old ) then
         !-- Chemical composition
         if ( rank == 0 ) then
            read( n_start_file ) work_old
            call map_field(work_old, work, rscheme_old, r_old, m2idx_old,   &
                 &         scale_xi, n_m_max_old, n_r_max_old, n_r_max_max, &
                 &         lBc=.false., l_phys_space=.true.)
         end if
         call scatter_from_rank0_to_mloc(work, xi_Mloc)

         if ( tscheme_family_old == 'MULTISTEP' ) then
            !-- Explicit time step
            do n_o=2,norder_exp_old
               if ( rank == 0 ) then
                  read( n_start_file ) work_old
                  call map_field(work_old, work, rscheme_old, r_old, m2idx_old,   &
                       &         scale_xi, n_m_max_old, n_r_max_old, n_r_max_max, &
                       &         lBc=.true., l_phys_space=l_coll_old)
               end if
               if ( n_o <= tscheme%norder_exp ) then
                  call scatter_from_rank0_to_mloc(work, dxidt%expl(:,:,n_o))
               end if
            end do

            !-- Implicit time step
            do n_o=2,norder_imp_lin_old-1
               if ( rank == 0 ) then
                  read( n_start_file ) work_old
                  call map_field(work_old, work, rscheme_old, r_old, m2idx_old,   &
                       &         scale_xi, n_m_max_old, n_r_max_old, n_r_max_max, &
                       &         lBc=.true., l_phys_space=l_coll_old)
               end if
               if ( n_o <= tscheme%norder_imp_lin-1 ) then
                  call scatter_from_rank0_to_mloc(work, dxidt%impl(:,:,n_o))
               end if
            end do
            if ( version > 3 ) then
               do n_o=2,norder_imp_old-1
                  if ( rank == 0 ) then
                     read( n_start_file ) work_old
                     call map_field(work_old, work, rscheme_old, r_old, m2idx_old,&
                          &         scale_xi, n_m_max_old, n_r_max_old,           &
                          &         n_r_max_max, lBc=.false.,                     &
                          &         l_phys_space=l_coll_old)
                  end if
                  if ( n_o <= tscheme%norder_imp-1 ) then
                     call scatter_from_rank0_to_mloc(work, dxidt%old(:,:,n_o))
                  end if
               end do
            end if
         end if

      end if

      !-- Reading 3-D fields
      if ( l_3D ) then
         if ( l_heat_3D_old ) then
            if ( rank == 0 ) then
               deallocate( work_old, work )
               allocate( work_old(lm_max_old,n_r_max_3D_old) )
               allocate( work(lm_max,n_r_max_3D) )
               read( n_start_file ) work_old
               call map_field_3D(work_old, work, rscheme_3D_old, r_3D_old,      &
                    &            lm2lmo,  scale_t, lm_max_old, n_r_max_3D_old,  &
                    &            n_r_max_3D_max, lBc=.false., l_phys_space=.true.)
            end if
            call scatter_from_rank0_to_lmloc(work, temp_3D_LMloc)

            if ( tscheme_family_old == 'MULTISTEP' ) then
               !-- Explicit time step
               do n_o=2,norder_exp_old
                  if ( rank == 0 ) then
                     read( n_start_file ) work_old
                     call map_field_3D(work_old, work, rscheme_3D_old, r_3D_old,  &
                          &            lm2lmo,  scale_t, lm_max_old,              &
                          &            n_r_max_3D_old, n_r_max_3D_max, lBc=.true.,&
                          &            l_phys_space=.true.)
                  end if
                  if ( n_o <= tscheme%norder_exp ) then
                     call scatter_from_rank0_to_lmloc(work, dTdt_3D%expl(:,:,n_o))
                  end if
               end do

               !-- Implicit time step
               do n_o=2,norder_imp_lin_old-1
                  if ( rank == 0 ) then
                     read( n_start_file ) work_old
                     call map_field_3D(work_old, work, rscheme_3D_old, r_3D_old,  &
                          &            lm2lmo,  scale_t, lm_max_old,              &
                          &            n_r_max_3D_old, n_r_max_3D_max, lBc=.true.,&
                          &            l_phys_space=.true.)
                  end if
                  if ( n_o <= tscheme%norder_imp_lin-1 ) then
                     call scatter_from_rank0_to_lmloc(work, dTdt_3D%impl(:,:,n_o))
                  end if
               end do

               do n_o=2,norder_imp_old-1
                  if ( rank == 0 ) then
                     read( n_start_file ) work_old
                     call map_field_3D(work_old, work, rscheme_3D_old, r_3D_old, &
                          &            lm2lmo,  scale_t, lm_max_old,             &
                          &            n_r_max_3D_old, n_r_max_3D_max,           &
                          &            lBc=.false., l_phys_space=.true.)
                  end if
                  if ( n_o <= tscheme%norder_imp-1 ) then
                     call scatter_from_rank0_to_lmloc(work, dTdt_3D%old(:,:,n_o))
                  end if
               end do
            end if
         else
            call abortRun('2-D to 3-D Not implemented yet')
         end if
      end if

      if ( rank == 0 ) then
         call rscheme_old%finalize(no_work_array=.true.)
         if ( l_3D ) call rscheme_3D_old%finalize(no_work_array=.true.)
      end if
      deallocate( r_old, work_old, work, m2idx_old, lm2lmo )


   end subroutine read_checkpoint
!------------------------------------------------------------------------------
   subroutine get_lm2lmo(lm2lmo, lm_max, l_max_old, minc_old)

      integer, intent(in) :: lm_max
      integer, intent(in) :: l_max_old
      integer, intent(in) :: minc_old

      !-- Output variable
      integer, intent(out) :: lm2lmo(lm_max)

      integer :: lm, l, m, lo, mo, lmo

      do lm=1,lm_max
         l = st_map%lm2l(lm)
         m = st_map%lm2m(lm)
         lm2lmo(lm) = -1 ! -1 means there's no data
         lmo = 0
         do mo=0,l_max_old,minc_old
            do lo=mo,l_max_old
               lmo=lmo+1
               if ( lo==l .and. mo==m ) then
                  lm2lmo(lm)=lmo ! data found in checkpoint
                  cycle
               end if
            end do
         end do
      end do

   end subroutine get_lm2lmo
!------------------------------------------------------------------------------
   subroutine map_field_3D(field_old, field_new, rscheme_old, r_3D_old,     &
              &            lm2lmo, scale_field, lm_max_old, n_r_max_3D_old, &
              &            n_r_max_3D_max, lBc, l_phys_space)

      !-- Input variables:
      class(type_rscheme), intent(in) :: rscheme_old
      integer,             intent(in) :: n_r_max_3D_old
      integer,             intent(in) :: lm_max_old
      complex(cp),         intent(in) :: field_old(lm_max_old,n_r_max_3D_old)
      real(cp),            intent(in) :: scale_field
      real(cp),            intent(in) :: r_3D_old(n_r_max_3D_old)
      integer,             intent(in) :: lm2lmo(lm_max)
      integer,             intent(in) :: n_r_max_3D_max
      logical,             intent(in) :: lBc
      logical,             intent(in) :: l_phys_space

      !-- Output variable:
      complex(cp), intent(out) :: field_new(lm_max, n_r_max_3D)

      !-- Local variables:
      integer :: lm, lmo, n_r
      complex(cp) :: radial_data(n_r_max_3D_max)

      do lm=1,lm_max

         lmo = lm2lmo(lm)

         if ( lmo > 0 ) then
            if ( n_r_max_3D /= n_r_max_3D_old .or.         &
            &    rscheme_3D%order_boundary /= rscheme_old%order_boundary&
            &    .or. rscheme_3D%version /= rscheme_old%version ) then

               do n_r=1,n_r_max_3D_old  ! copy on help arrays
                  radial_data(n_r)=field_old(lmo,n_r)
               end do
               call map_field_r(radial_data, rscheme_3D, rscheme_old, r_3D_old,  &
                    &           n_r_max_3D, n_r_max_3D_old, n_r_max_3D_max, lBc, &
                    &           l_phys_space)
               do n_r=1,n_r_max_3D
                  field_new(lm,n_r)=scale_field*radial_data(n_r)
               end do
            else
               do n_r=1,n_r_max_3D
                  field_new(lm,n_r)=scale_field*field_old(lmo,n_r)
               end do
            end if
         else
            do n_r=1,n_r_max_3D
               field_new(lm,n_r)=zero
            end do
         end if

      end do

   end subroutine map_field_3D
!------------------------------------------------------------------------------
   subroutine map_field(field_old, field_new, rscheme_old, r_old, m2idx_old, &
              &         scale_field, n_m_max_old, n_r_max_old, n_r_max_max,  &
              &         lBc, l_phys_space)

      !-- Input variables:
      class(type_rscheme), intent(in) :: rscheme_old
      integer,             intent(in) :: n_r_max_old
      integer,             intent(in) :: n_m_max_old
      integer,             intent(in) :: n_r_max_max
      complex(cp),         intent(in) :: field_old(n_m_max_old, n_r_max_old)
      real(cp),            intent(in) :: r_old(n_r_max_old)
      integer,             intent(in) :: m2idx_old(0:)
      real(cp),            intent(in) :: scale_field
      logical,             intent(in) :: lBc
      logical,             intent(in) :: l_phys_space

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
               call map_field_r(radial_data, rscheme, rscheme_old, r_old, &
                    &           n_r_max, n_r_max_old, n_r_max_max, lBc,   &
                    &           l_phys_space)
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
   subroutine map_field_r(radial_data, rscheme_new, rscheme_old, r_old, &
              &           n_r_max_new, n_r_max_old, n_r_max_max,lBc, l_phys_space)

      !-- Input variables
      class(type_rscheme), intent(in) :: rscheme_new
      class(type_rscheme), intent(in) :: rscheme_old
      integer,             intent(in) :: n_r_max_new
      integer,             intent(in) :: n_r_max_old
      integer,             intent(in) :: n_r_max_max
      real(cp), intent(in) :: r_old(n_r_max_old)
      logical,  intent(in) :: lBc
      logical,  intent(in) :: l_phys_space

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
      if ( lBc .and. l_phys_space ) then
         radial_data(1)=two*radial_data(2)-radial_data(3)
      end if

      !-- If **both** the old and the new schemes are Chebyshev, we can
      !-- use costf to get the new data
      !-- Both have also to use the same mapping (order_boundary is a proxy of l_map
      if ( rscheme%version == 'cheb' .and. rscheme_old%version == 'cheb' &
      &   .and. rscheme%order_boundary == rscheme_old%order_boundary  ) then

         do n_r=1,n_r_max_old
            radial_data_real(n_r) =  real(radial_data(n_r))
            radial_data_imag(n_r) = aimag(radial_data(n_r))
         end do
         if ( l_phys_space ) then
            call rscheme_old%costf1(radial_data_real,n_r_max_old)
            call rscheme_old%costf1(radial_data_imag,n_r_max_old)
         end if

         !----- Fill up cheb polynomial with zeros:
         if ( n_r_max_new>n_r_max_old ) then
            n_r_index_start=n_r_max_old+1
            do n_r=n_r_index_start,n_r_max_new
               radial_data_real(n_r)=0.0_cp
               radial_data_imag(n_r)=0.0_cp
            end do
         end if

         !----- Now transform to new radial grid points:
         if ( l_phys_space ) then
            call rscheme_new%costf1(radial_data_real,n_r_max_new)
            call rscheme_new%costf1(radial_data_imag,n_r_max_new)
         end if
         !----- Rescale :
         cheb_norm_old=sqrt(two/real(n_r_max_old-1,kind=cp))
         cheb_fac=cheb_norm_old/rscheme_new%rnorm
         do n_r=1,n_r_max_new
            radial_data(n_r)=cheb_fac*cmplx(radial_data_real(n_r), &
            &                               radial_data_imag(n_r),cp)
         end do

      !-- If either the old grid or the new grid is FD, we use a 
      !-- polynomial interpolation
      else

         allocate( work(n_r_max_new) )

         !-- Interpolate data and store into a work array
         do n_r=1,n_r_max_new

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
         do n_r=1,n_r_max_new
            radial_data(n_r)=work(n_r)
         end do

         deallocate( work )

      end if

   end subroutine map_field_r
!------------------------------------------------------------------------------
end module checkpoints
