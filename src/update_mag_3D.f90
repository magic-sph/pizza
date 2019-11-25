module update_mag_3D_mod

   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation_3D, only: n_r_max_3D, lm_max, l_max
   use radial_functions, only: r_3D, or1_3D, or2_3D, rscheme_3D
   use namelists, only: kbotb, ktopb, BdiffFac
   use blocking_lm, only: st_map, lo_map, lo_sub_map, chunksize
   use blocking, only: lmStart, lmStop
   !use horizontal, only: hdif_B   !-- WARNING:: Totally wrong -> to be modified!!
   use init_fields, only: bpeaktop, bpeakbot
   use parallel_mod, only: rank
   use algebra, only: prepare_full_mat, solve_full_mat
   use radial_der, only: get_ddr, get_dr
   use fields, only:  work_b_LMloc, work_j_LMloc ! reduce number of  arrays
   use constants, only: zero, one, two, sq4pi
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   
   implicit none

   private

   !-- Local work arrays:
   complex(cp), allocatable :: rhs1(:,:,:),rhs2(:,:,:)
   real(cp), allocatable :: bMat(:,:,:)
   real(cp), allocatable :: jMat(:,:,:)
   integer, allocatable :: bPivot(:,:)
   integer, allocatable :: jPivot(:,:)
#ifdef WITH_PRECOND_BJ
   real(cp), allocatable :: bMat_fac(:,:)
   real(cp), allocatable :: jMat_fac(:,:)
#endif
   logical, public, allocatable :: lBmat(:)


   integer :: maxThreads

   public :: initialize_update_mag_3D, update_mag_3D, finalize_update_mag_3D, &
   &         finish_exp_mag_3D, get_mag_3D_rhs_imp

contains

   subroutine initialize_update_mag_3D

      allocate( bMat(n_r_max_3D,n_r_max_3D,l_max) )
      allocate( jMat(n_r_max_3D,n_r_max_3D,l_max) )
      bytes_allocated = bytes_allocated+(2*n_r_max_3D*n_r_max_3D*l_max) &
      &                 *SIZEOF_DEF_REAL
      allocate( bPivot(n_r_max_3D,l_max) )
      allocate( jPivot(n_r_max_3D,l_max) )
      bytes_allocated = bytes_allocated+2*n_r_max_3D*l_max*SIZEOF_INTEGER

#ifdef WITH_PRECOND_BJ
      allocate(bMat_fac(n_r_max_3D,l_max))
      allocate(jMat_fac(n_r_max_3D,l_max))
      bytes_allocated = bytes_allocated+2*n_r_max_3D*l_max*SIZEOF_DEF_REAL
#endif
      allocate( lBmat(0:l_max) )
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif

      allocate(rhs1(n_r_max_3D,lo_sub_map%sizeLMB2max,0:maxThreads-1))
      allocate(rhs2(n_r_max_3D,lo_sub_map%sizeLMB2max,0:maxThreads-1))
      bytes_allocated = bytes_allocated+2*n_r_max_3D*lo_sub_map%sizeLMB2max*&
      &                 maxThreads*SIZEOF_DEF_COMPLEX


   end subroutine initialize_update_mag_3D
!-----------------------------------------------------------------------------
   subroutine finalize_update_mag_3D

      deallocate( bMat, jMat, bPivot, jPivot, lBmat )

#ifdef WITH_PRECOND_BJ
      deallocate( bMat_fac, jMat_fac )
#endif
      deallocate( rhs1, rhs2 )

   end subroutine finalize_update_mag_3D
!-----------------------------------------------------------------------------
   subroutine update_mag_3D(b_3D, db_3D, ddB_3D, dBdt_3D, aj_3D, dj_3D, &
              &             djdt_3D, tscheme, lMat, nLMB)
      !                                                                   
      !  Calculated update of magnetic field potential and the time       
      !  stepping arrays dbdtLast, ...                                    
      !                                                                   
      !  updates the magnetic field potentials b, aj and
      !  their derivatives,
      !  adds explicit part to time derivatives of b and j
      !

      !-- Input of variables:
      logical,             intent(in) :: lMat
      class(type_tscheme), intent(in) :: tscheme
      integer,             intent(in) :: nLMB

      !-- Input/output of scalar potentials and time stepping arrays:
      complex(cp), intent(inout) :: b_3D(lmStart:lmStop,n_r_max_3D)
      complex(cp), intent(inout) :: aj_3D(lmStart:lmStop,n_r_max_3D)
      complex(cp), intent(out) :: db_3D(lmStart:lmStop,n_r_max_3D)
      complex(cp), intent(out) :: ddb_3D(lmStart:lmStop,n_r_max_3D)
      complex(cp), intent(out) :: dj_3D(lmStart:lmStop,n_r_max_3D)

      type(type_tarray), intent(inout) :: dBdt_3D
      type(type_tarray), intent(inout) :: djdt_3D

      !-- Local variables:
      integer :: l1,m1               ! degree and order
      integer :: lmStart_00          ! max(2,lmStart)
      integer :: lm1,lm,lmB          ! position of (l,m) in array
      integer :: nLMB2
      integer :: n_r_out             ! No of cheb polynome (degree+1)
      integer :: n_r                 ! No of radial grid point

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: threadid
      integer :: iChunk,nChunks,size_of_last_chunk,lmB0
      
      if ( lMat ) lBmat(:)=.false.

      nLMBs2(1:) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m
      lmStart_00 = max(2,lmStart)

      !-- Calculation of the implicit part
      call get_mag_3D_rhs_imp(b_3D, db_3D, ddb_3D,               &
           &                 dBdt_3D%old(:,:,tscheme%istage),    &
           &                 dBdt_3D%impl(:,:,tscheme%istage),   &
           &                 aj_3D, dj_3D,                       &
           &                 djdt_3D%old(:,:,tscheme%istage),    &
           &                 djdt_3D%impl(:,:,tscheme%istage),   &
           &                  tscheme%l_imp_calc_rhs(tscheme%istage))

      !-- Now assemble the right hand side and store it in work_b_LMloc and work_j_LMloc
      call tscheme%set_imex_rhs(work_b_LMloc, dBdt_3D, lmStart, lmStop, n_r_max_3D)
      call tscheme%set_imex_rhs(work_j_LMloc, djdt_3D, lmStart, lmStop, n_r_max_3D)

      ! This is a loop over all l values which should be treated on
      ! the actual MPI rank
      !$OMP PARALLEL default(shared)
      !$OMP SINGLE
      do nLMB2=1,nLMBs2(nLMB)
         !$OMP TASK default(shared) &
         !$OMP firstprivate(nLMB2) &
         !$OMP private(lmB,lm,lm1,l1,m1,nR,iChunk,nChunks,size_of_last_chunk) &
         !$OMP private(dbdt_ic,djdt_ic,fac,bpeaktop,ff,cimp,aimp,threadid)

         ! determine the number of chunks of m
         ! total number for l1 is sizeLMB2(nLMB2,nLMB)
         ! chunksize is given
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk = chunksize + (sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         l1=lm22l(1,nLMB2,nLMB)
         if ( l1 > 0 ) then
            if ( .not. lBmat(l1) ) then
#ifdef WITH_PRECOND_BJ
               call get_bMat(tscheme,l1,one,      & !-> 'one' to fix for hdif: to be modified
                    &        bMat(:,:,l1),bPivot(:,l1), bMat_fac(:,l1),&
                    &        jMat(:,:,l1),jPivot(:,l1), jMat_fac(:,l1))
#else
               call get_bMat(tscheme,l1,one, & !-> 'one' to fix for hdif: to be modified
                    &        bMat(:,:,l1),bPivot(:,l1),           &
                    &        jMat(:,:,l1),jPivot(:,l1) )
#endif
               lBmat(l1)=.true.
            end if
         end if

         do iChunk=1,nChunks
            !$OMP TASK if (nChunks>1) default(shared) &
            !$OMP firstprivate(iChunk) &
            !$OMP private(lmB0,lmB,lm,lm1,m1,n_r) &
            !$OMP private(threadid)
#ifdef WITHOMP
            threadid = omp_get_thread_num()
#else
            threadid = 0
#endif
            lmB0=(iChunk-1)*chunksize
            lmB=lmB0

            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               !do lm=1,sizeLMB2(nLMB2,nLMB)
               lm1=lm22lm(lm,nLMB2,nLMB)
               !l1 =lm22l(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( l1 > 0 ) then
                  lmB=lmB+1
                  !-------- Magnetic boundary conditions, outer core:
                  !         Note: the CMB condition is not correct if we assume free slip
                  !         and a conducting mantle

                  rhs1(1,lmB,threadid)         =zero
                  rhs1(n_r_max_3D,lmB,threadid)=zero
                  if ( kbotb == 2 ) rhs1(n_r_max_3D-1,lmB,threadid)=zero

                  rhs2(1,lmB,threadid)         =zero
                  rhs2(n_r_max_3D,lmB,threadid)=zero

                  do n_r=2,n_r_max_3D-1
                     rhs1(n_r,lmB,threadid)=work_b_LMloc(lm1,n_r)
                     rhs2(n_r,lmB,threadid)=work_j_LMloc(lm1,n_r)
                  end do
               end if ! l>0
            end do ! loop over lm in block

            if ( lmB > lmB0 ) then
#ifdef WITH_PRECOND_BJ
               do lm=lmB0+1,lmB
                  do n_r=1,n_r_max_3D
                     rhs1(n_r,lm,threadid)=rhs1(n_r,lm,threadid)*bMat_fac(n_r,l1)
                     rhs2(n_r,lm,threadid)=rhs2(n_r,lm,threadid)*jMat_fac(n_r,l1)
                  end do
               end do
#endif
               call solve_full_mat(bMat(:,:,l1),n_r_max_3D,n_r_max_3D, &
                    &         bPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid),&
                    &         lmB-lmB0)
               call solve_full_mat(jMat(:,:,l1),n_r_max_3D,n_r_max_3D, &
                    &         jPivot(:,l1),rhs2(:,lmB0+1:lmB,threadid),&
                    &         lmB-lmB0)
            end if

            !----- Update magnetic field in cheb space:
            !PERFON('upB_set')
            lmB=lmB0
            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( l1 > 0 ) then
                  lmB=lmB+1
                  if ( m1 > 0 ) then
                     do n_r_out=1,rscheme_3D%n_max  ! outer core
                        b_3D(lm1,n_r_out) =rhs1(n_r_out,lmB,threadid)
                        aj_3D(lm1,n_r_out)=rhs2(n_r_out,lmB,threadid)
                     end do
                  else
                     do n_r_out=1,rscheme_3D%n_max   ! outer core
                        b_3D(lm1,n_r_out)= cmplx(real(rhs1(n_r_out,lmB,threadid)), &
                        &                        0.0_cp,kind=cp)
                        aj_3D(lm1,n_r_out)= cmplx(real(rhs2(n_r_out,lmB,threadid)),&
                        &                         0.0_cp,kind=cp)
                     end do
                  end if
               end if
            end do
            !PERFOFF
            !$OMP END TASK
         end do
         !$OMP END TASK
      end do      ! end of loop over lm blocks
      !$OMP END SINGLE
      !$OMP END PARALLEL

      !-- Set cheb modes > rscheme_3D%n_max to zero (dealiazing)
      do n_r_out=rscheme_3D%n_max+1,n_r_max_3D
         do lm1=lmStart_00,lmStop
            b_3D(lm1,n_r_out) =zero
            aj_3D(lm1,n_r_out)=zero
         end do
      end do

      !-- Bring magnetic vectors back to physical space
      call get_ddr(b_3D, db_3D, ddb_3D, lmStart, lmStop, &
           &       n_r_max_3D, rscheme_3D, l_dct=.false.)
      call rscheme_3D%costf1(b_3D, lmStart, lmStop, n_r_max_3D)
      call get_dr(aj_3D, dj_3D, lmStart, lmStop,          &
           &      n_r_max_3D, rscheme_3D, l_dct_in=.false.)
      call rscheme_3D%costf1(aj_3D, lmStart, lmStop, n_r_max_3D)

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dBdt_3D, lmStart, lmStop, n_r_max_3D)
      call tscheme%rotate_imex(djdt_3D, lmStart, lmStop, n_r_max_3D)

   end subroutine update_mag_3D
!-------------------------------------------------------------------------------
   subroutine finish_exp_mag_3D(dVxBh_LMloc, dj_exp_last)

      !-- Input variables
      complex(cp), intent(inout) :: dVxBh_LMloc(lmStart:lmStop,n_r_max_3D)

      !-- Output variables
      complex(cp), intent(inout) :: dj_exp_last(lmStart:lmStop,n_r_max_3D)

      !-- Local variables
      integer :: n_r, lm

      call get_dr( dVxBh_LMloc,work_j_LMloc, lmStart, lmStop, &
           &       n_r_max_3D, rscheme_3D, nocopy=.true. )

      do n_r=1,n_r_max_3D
         do lm=max(2,lmStart),lmStop
            dj_exp_last(lm,n_r)=              dj_exp_last(lm,n_r) &
            &                   +or2_3D(n_r)*work_j_LMloc(lm,n_r)
         end do
      end do

   end subroutine finish_exp_mag_3D
!-------------------------------------------------------------------------------
   subroutine get_mag_3D_rhs_imp(b_3D, db_3D, ddb_3D, B_last, dB_imp_last,  &
              &                  aj_3D,dj_3D, aj_last, dj_imp_last,         &
              &                  l_calc_lin_rhs)

      !-- Input variables
      complex(cp), intent(inout) :: b_3D(lmStart:lmStop,n_r_max_3D)
      complex(cp), intent(inout) :: aj_3D(lmStart:lmStop,n_r_max_3D)
      logical,     intent(in) :: l_calc_lin_rhs

      !-- Output variable
      complex(cp), intent(inout) :: db_3D(lmStart:lmStop,n_r_max_3D)
      complex(cp), intent(inout) :: ddb_3D(lmStart:lmStop,n_r_max_3D)
      complex(cp), intent(inout) :: dj_3D(lmStart:lmStop,n_r_max_3D)

      complex(cp), intent(out) :: B_Last(lmStart:lmStop,n_r_max_3D)
      complex(cp), intent(out) :: aj_Last(lmStart:lmStop,n_r_max_3D)
      complex(cp), intent(out) :: dB_imp_last(lmStart:lmStop,n_r_max_3D)
      complex(cp), intent(out) :: dj_imp_last(lmStart:lmStop,n_r_max_3D)

      !-- Local variables
      real(cp) :: dL
      integer :: n_r, lm, l1, lmStart_00
      integer, pointer :: lm2l(:)

      lm2l(1:lm_max) => lo_map%lm2l
      lmStart_00 = max(2,lmStart)

      do n_r=1,n_r_max_3D
         do lm=lmStart_00,lmStop
            l1 = lm2l(lm)
            dL = real(l1*(l1+1),cp)
            B_last(lm,n_r) = dL*or2_3D(n_r)* b_3D(lm,n_r)
            aj_last(lm,n_r)= dL*or2_3D(n_r)*aj_3D(lm,n_r)
         end do
      end do

      if ( l_calc_lin_rhs ) then
         call get_ddr(b_3D, db_3D, ddb_3D, lmStart, lmStop,  &
              &       n_r_max_3D, rscheme_3D)
         call get_ddr(aj_3D,dj_3D,work_j_LMloc,lmStart,lmStop,&
              &       n_r_max_3D, rscheme_3D)

         !-- Calculate explicit time step part:
         do n_r=1,n_r_max_3D
            do lm=lmStart_00,lmStop
               l1 = lm2l(lm)
               dL = real(l1*(l1+1),cp)
               dB_imp_last(lm,n_r)=BdiffFac* (      dL * or2_3D(n_r) * (   &
               &             ddb_3D(lm,n_r)-dL*or2_3D(n_r)*b_3D(lm,n_r) ) )
               dj_imp_last(lm,n_r)=BdiffFac* (      dL * or2_3D(n_r) * (   &
               &       work_j_LMloc(lm,n_r)-dL*or2_3D(n_r)*aj_3D(lm,n_r) ) )
            end do
         end do
      end if

   end subroutine get_mag_3D_rhs_imp
!-----------------------------------------------------------------------------
#ifdef WITH_PRECOND_BJ
   subroutine get_bMat(tscheme,l,hdif,bMat,bPivot,bMat_fac,jMat,jPivot,jMat_fac)
#else
   subroutine get_bMat(tscheme,l,hdif,bMat,bPivot,jMat,jPivot)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrices 
      !  bmat(i,j) and ajmat for the dynamo equations.                    
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme        ! time step
      integer,  intent(in) :: l
      real(cp), intent(in) :: hdif
    
      !-- Output variables:
      real(cp), intent(out) :: bMat(n_r_max_3D,n_r_max_3D)
      integer,  intent(out) :: bPivot(n_r_max_3D)
      real(cp), intent(out) :: jMat(n_r_max_3D,n_r_max_3D)
      integer,  intent(out) :: jPivot(n_r_max_3D)
#ifdef WITH_PRECOND_BJ
      real(cp), intent(out) :: bMat_fac(n_r_max_3D), jMat_fac(n_r_max_3D)
#endif
 
      !-- local variables:
      integer :: info, n_r, n_r_out
      real(cp) :: dLh
 
      dLh=real(l*(l+1),kind=cp)

      !-- matricies depend on degree l but not on order m,
      !   we thus have to construct bmat and ajmat for each l:

      !----- boundary conditions for outer core field:
      do n_r_out=1,rscheme_3D%n_max
         if ( ktopb == 1 .or. ktopb == 3 ) then
            !-------- at CMB (n_r=1):
            !         the internal poloidal field should fit a potential
            !         field (matrix bmat) and the toroidal field has to
            !         vanish (matrix ajmat).
            bMat(1,n_r_out)=rscheme_3D%rnorm*(                   &
            &                      rscheme_3D%drMat(1,n_r_out) + &
            &  real(l,cp)*or1_3D(1)*rscheme_3D%rMat(1,n_r_out)  )

            jMat(1,n_r_out)=rscheme_3D%rnorm*rscheme_3D%rMat(1,n_r_out)
         else if ( ktopb == 2 ) then
            call abortRun('! Boundary condition ktopb=2 not defined!')
         else if ( ktopb == 3 ) then
            call abortRun('! Boundary condition ktopb=3 not defined!')
         else if ( ktopb == 4 ) then
            !----- pseudo vacuum condition, field has only
            !      a radial component, horizontal components
            !      vanish when aj and db are zero:
            bMat(1,n_r_out)=rscheme_3D%rnorm*rscheme_3D%drMat(1,n_r_out)
            jMat(1,n_r_out)=rscheme_3D%rnorm* rscheme_3D%rMat(1,n_r_out)
         end if
         !-------- at IC (n_r=n_r_max_3D):
         if ( kbotb == 1 ) then
            !----------- insulating IC, field has to fit a potential field:
            bMat(n_r_max_3D,n_r_out)=rscheme_3D%rnorm*(                            &
            &                               rscheme_3D%drMat(n_r_max_3D,n_r_out) - &
            & real(l+1,cp)*or1_3D(n_r_max_3D)*rscheme_3D%rMat(n_r_max_3D,n_r_out) )
            jMat(n_r_max_3D,n_r_out)=rscheme_3D%rnorm*rscheme_3D%rMat(n_r_max_3D,n_r_out)
         else if ( kbotb == 2 ) then
            !----------- perfect conducting IC
            bMat(n_r_max_3D-1,n_r_out)=rscheme_3D%rnorm*rscheme_3D%d2rMat(n_r_max_3D,n_r_out)
            jMat(n_r_max_3D,n_r_out)  =rscheme_3D%rnorm* rscheme_3D%drMat(n_r_max_3D,n_r_out)
         else if ( kbotb == 3 ) then
            call abortRun('! Boundary condition kbotb=3 not defined!')
         else if ( kbotb == 4 ) then
            !----- Pseudovacuum conduction at lower boundary:
            bMat(n_r_max_3D,n_r_out)=rscheme_3D%rnorm*rscheme_3D%drMat(n_r_max_3D,n_r_out)
            jMat(n_r_max_3D,n_r_out)= rscheme_3D%rnorm*rscheme_3D%rMat(n_r_max_3D,n_r_out)
         end if
      end do ! loop over cheb modes !
      !----- fill up with zeros:
      if ( rscheme_3D%n_max < n_r_max_3D ) then ! fill with zeros !
         do n_r_out=rscheme_3D%n_max+1,n_r_max_3D
            bMat(1,n_r_out)=0.0_cp
            jMat(1,n_r_out)=0.0_cp
            if ( kbotb == 1 ) then
               bMat(n_r_max_3D,n_r_out)  =0.0_cp
               jMat(n_r_max_3D,n_r_out)  =0.0_cp
            else if ( kbotb == 2 ) then
               bMat(n_r_max_3D-1,n_r_out)=0.0_cp
               jMat(n_r_max_3D,n_r_out)  =0.0_cp
            else if ( kbotb == 3 ) then
               call abortRun('! Boundary condition kbotb=3 not defined!')
            else if ( kbotb == 4 ) then
               bMat(n_r_max_3D,n_r_out)  =0.0_cp
               jMat(n_r_max_3D,n_r_out)  =0.0_cp
            end if
         end do
      end if
    
      !----- Other points:
      do n_r_out=1,n_r_max_3D
         do n_r=2,n_r_max_3D-1
            bMat(n_r,n_r_out)=rscheme_3D%rnorm*(                    &
            &        dLh*or2_3D(n_r)*rscheme_3D%rMat(n_r,n_r_out) - &
            & BdiffFac*tscheme%wimp_lin(1)*hdif*dLh*or2_3D(n_r) * ( &
            &                      rscheme_3D%d2rMat(n_r,n_r_out) - &
            &        dLh*or2_3D(n_r)*rscheme_3D%rMat(n_r,n_r_out) ) )
    
            jMat(n_r,n_r_out)=rscheme_3D%rnorm*(                    &
            &        dLh*or2_3D(n_r)*rscheme_3D%rMat(n_r,n_r_out) - &
            & BdiffFac*tscheme%wimp_lin(1)*hdif*dLh*or2_3D(n_r) * ( &
            &                      rscheme_3D%d2rMat(n_r,n_r_out) - &
            &        dLh*or2_3D(n_r)*rscheme_3D%rMat(n_r,n_r_out) ) )
         end do
      end do
    
      !----- normalization for highest and lowest Cheb mode:
      do n_r=1,n_r_max_3D
         bMat(n_r,1)         =rscheme_3D%boundary_fac*bMat(n_r,1)
         bMat(n_r,n_r_max_3D)=rscheme_3D%boundary_fac*bMat(n_r,n_r_max_3D)
         jMat(n_r,1)         =rscheme_3D%boundary_fac*jMat(n_r,1)
         jMat(n_r,n_r_max_3D)=rscheme_3D%boundary_fac*jMat(n_r,n_r_max_3D)
      end do
    
#ifdef WITH_PRECOND_BJ
      ! compute the linesum of each line
      do n_r=1,n_r_max_3D
         bMat_fac(n_r)=one/maxval(abs(bMat(n_r,:)))
         jMat_fac(n_r)=one/maxval(abs(jMat(n_r,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do n_r=1,n_r_max_3D
         bMat(n_r,:) = bMat(n_r,:)*bMat_fac(n_r)
         jMat(n_r,:) = jMat(n_r,:)*jMat_fac(n_r)
      end do
#endif

      !----- LU decomposition:
      call prepare_full_mat(bMat,n_r_max_3D,n_r_max_3D,bPivot,info)
      if ( info /= 0 ) then
         call abortRun('Singular matrix BMat in get_bmat!')
      end if

      !----- LU decomposition:
      call prepare_full_mat(jMat,n_r_max_3D,n_r_max_3D,jPivot,info)
      if ( info /= 0 ) then
         call abortRun('! Singular matrix ajmat in get_bmat!')
      end if

   end subroutine get_bMat
!-----------------------------------------------------------------------------
end module update_mag_3D_mod
