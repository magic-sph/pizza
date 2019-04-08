module update_temp_3D_mod

   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation_3D, only: n_r_max_3D, lm_max, l_max
   use radial_functions, only: or1_3D, or2_3D, rscheme_3D
   use namelists, only: kbott, ktopt, TdiffFac
   use blocking_lm, only: nLMBs, st_map, lo_map, lo_sub_map, lmStartB, lmStopB, &
       &                  llm, ulm, chunksize
   use horizontal, only: dLh
   use parallel_mod, only: rank
   use algebra, only: prepare_full_mat, solve_full_mat
   use radial_der, only: get_ddr, get_dr
   use fields, only:  work_LMloc
   use constants, only: zero, one, two, sq4pi
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray

   implicit none

   private

   !-- Local variables
   complex(cp), allocatable :: rhs1(:,:,:)
   real(cp), allocatable :: T0Mat(:,:)     ! for l=m=0  
   real(cp), allocatable :: TMat(:,:,:)
   integer, allocatable :: T0Pivot(:)
   integer, allocatable :: TPivot(:,:)
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: TMat_fac(:,:)
   real(cp), allocatable :: T0Mat_fac(:)
#endif
   logical, allocatable :: lTmat(:)

   integer :: maxThreads

   public :: initialize_update_temp_3D, update_temp_3D, finalize_update_temp_3D, &
   &         finish_exp_temp_3D, get_temp_3D_rhs_imp

contains

   subroutine initialize_update_temp_3D

      allocate( T0Mat(n_r_max_3D,n_r_max_3D) )      ! for l=m=0  
      allocate( TMat(n_r_max_3D,n_r_max_3D,l_max) )
      bytes_allocated = bytes_allocated+(n_r_max_3D*n_r_max_3D*(1+l_max))* &
      &                 SIZEOF_DEF_REAL
      allocate( T0Pivot(n_r_max_3D) )
      allocate( TPivot(n_r_max_3D,l_max) )
      bytes_allocated = bytes_allocated+(n_r_max_3D+n_r_max_3D*l_max)*SIZEOF_INTEGER
#ifdef WITH_PRECOND_S
      allocate(TMat_fac(n_r_max_3D,l_max))
      allocate(T0Mat_fac(n_r_max_3D))
      bytes_allocated = bytes_allocated+n_r_max_3D*SIZEOF_DEF_REAL
      bytes_allocated = bytes_allocated+n_r_max_3D*l_max*SIZEOF_DEF_REAL
#endif
      allocate( lTmat(0:l_max) )
      lTmat(:)=.false.
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif
      allocate( rhs1(n_r_max_3D,lo_sub_map%sizeLMB2max,0:maxThreads-1) )
      bytes_allocated = bytes_allocated + n_r_max_3D*lo_sub_map%sizeLMB2max*&
      &                 maxThreads*SIZEOF_DEF_COMPLEX

   end subroutine initialize_update_temp_3D
!------------------------------------------------------------------------------
   subroutine finalize_update_temp_3D

      deallocate( T0Mat, TMat, T0Pivot, TPivot, lTmat )
#ifdef WITH_PRECOND_S
      deallocate( TMat_fac )
      deallocate( T0Mat_fac )
#endif
      deallocate( rhs1 )
  
   end subroutine finalize_update_temp_3D
!------------------------------------------------------------------------------
   subroutine update_temp_3D(temp_3D, dtemp_3D, dTdt_3D, tscheme, lMat, nLMB)
      !
      !  updates the temperature field T and its radial derivatives
      !  adds explicit part to time derivatives of s
      !

      !-- Input of variables:
      logical,             intent(in) :: lMat
      class(type_tscheme), intent(in) :: tscheme
      integer,             intent(in) :: nLMB

      !-- Input/output of scalar fields:
      complex(cp),       intent(inout) :: temp_3D(llm:ulm,n_r_max_3D)
      complex(cp),       intent(inout) :: dtemp_3D(llm:ulm,n_r_max_3D)
      type(type_tarray), intent(inout) :: dTdt_3D


      !-- Local variables:
      integer :: l1,m1              ! degree and order
      integer :: lm1,lmB,lm         ! position of (l,m) in array
      integer :: nLMB2
      integer :: nR                 ! counts radial grid points
      integer :: n_r_out             ! counts cheb modes
      real(cp) ::  rhs(n_r_max_3D) ! real RHS for l=m=0

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: threadid
      integer :: iChunk,nChunks,size_of_last_chunk,lmB0

      if ( lMat ) lTmat(:)=.false.

      nLMBs2(1:nLMBs) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      !-- Calculation of the implicit part
      call get_temp_3D_rhs_imp(temp_3D, dtemp_3D,                 &
           &                   dTdt_3D%old(:,:,tscheme%istage),   &
           &                   dTdt_3D%impl(:,:,tscheme%istage),  &
           &                   tscheme%l_imp_calc_rhs(tscheme%istage))

      !-- Now assemble the right hand side and store it in work_LMloc
      call tscheme%set_imex_rhs(work_LMloc, dTdt_3D, llm, ulm, n_r_max_3D)

      ! one subblock is linked to one l value and needs therefore once the matrix
      !$OMP PARALLEL default(shared)
      !$OMP SINGLE
      do nLMB2=1,nLMBs2(nLMB)
         ! this inner loop is in principle over the m values which belong to the
         ! l value
         !$OMP TASK default(shared) &
         !$OMP firstprivate(nLMB2) &
         !$OMP private(lm,lm1,l1,m1,lmB,threadid) &
         !$OMP private(nChunks,size_of_last_chunk,iChunk)
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk = chunksize + (sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         ! This task treats one l given by l1
         l1=lm22l(1,nLMB2,nLMB)
         !write(*,"(3(A,I3),A)") "Launching task for nLMB2=",nLMB2," (l=",l1,") and scheduling ",nChunks," subtasks."

         if ( l1 == 0 ) then
            if ( .not. lTmat(l1) ) then
#ifdef WITH_PRECOND_S
               call get_T0Mat(tscheme,T0Mat,T0Pivot,T0Mat_fac)
#else
               call get_T0Mat(tscheme,T0Mat,T0Pivot)
#endif
               lTmat(l1)=.true.
            end if
         else
            if ( .not. lTmat(l1) ) then
#ifdef WITH_PRECOND_S
               call get_TMat(tscheme,l1, TMat(:,:,l1),TPivot(:,l1),TMat_fac(:,l1))
#else
               call get_TMat(tscheme,l1, TMat(:,:,l1),TPivot(:,l1))
#endif
               lTmat(l1)=.true.
            end if
          end if

         do iChunk=1,nChunks
            !$OMP TASK default(shared) &
            !$OMP firstprivate(iChunk) &
            !$OMP private(lmB0,lmB,lm,lm1,m1,nR,n_r_out) &
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

               if ( l1 == 0 ) then
                  rhs(1)         =0.0_cp*sq4pi
                  rhs(n_r_max_3D)=1.0_cp*sq4pi ! To be fixed later
                  do nR=2,n_r_max_3D-1
                     rhs(nR)=real(work_LMloc(lm1,nR))
                  end do

#ifdef WITH_PRECOND_S
                  rhs(:) = T0Mat_fac(:)*rhs(:)
#endif

                  call solve_full_mat(T0Mat,n_r_max_3D,n_r_max_3D,T0Pivot,rhs)

               else ! l1  /=  0
                  lmB=lmB+1

                  rhs1(1,lmB,threadid)         =zero
                  rhs1(n_r_max_3D,lmB,threadid)=zero
#ifdef WITH_PRECOND_S
                  rhs1(1,lmB,threadid)         =TMat_fac(1,l1)*rhs1(1,lmB,threadid)
                  rhs1(n_r_max_3D,lmB,threadid)=TMat_fac(1,l1)*rhs1(n_r_max_3D,lmB,threadid)
#endif
                  do nR=2,n_r_max_3D-1
                     rhs1(nR,lmB,threadid)=work_LMloc(lm1,nR)
#ifdef WITH_PRECOND_S
                     rhs1(nR,lmB,threadid)=TMat_fac(nR,l1)*rhs1(nR,lmB,threadid)
#endif
                  end do
               end if
            end do
            !PERFOFF

            if ( lmB  >  lmB0 ) then
               call solve_full_mat(TMat(:,:,l1),n_r_max_3D,n_r_max_3D,      &
                    &              TPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid),&
                    &              lmB-lmB0)
            end if

            lmB=lmB0
            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( l1 == 0 ) then
                  do n_r_out=1,rscheme_3D%n_max
                     temp_3D(lm1,n_r_out)=rhs(n_r_out)
                  end do
               else
                  lmB=lmB+1
                  if ( m1 > 0 ) then
                     do n_r_out=1,rscheme_3D%n_max
                        temp_3D(lm1,n_r_out)=rhs1(n_r_out,lmB,threadid)
                     end do
                  else
                     do n_r_out=1,rscheme_3D%n_max
                        temp_3D(lm1,n_r_out)= cmplx(real(rhs1(n_r_out,lmB,threadid)), &
                        &                          0.0_cp,kind=cp)
                     end do
                  end if
               end if
            end do
            !PERFOFF
            !$OMP END TASK
         end do
         !$OMP END TASK
      end do     ! loop over lm blocks
      !$OMP END SINGLE
      !$OMP END PARALLEL

      do n_r_out=rscheme_3D%n_max+1,n_r_max_3D
         do lm1=llm,ulm
            temp_3D(lm1,n_r_out)=zero
         end do
      end do

      !-- Bring temperature back to physical space
      call rscheme_3D%costf1(temp_3D, llm, ulm, n_r_max_3D)

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dTdt_3D, llm, ulm, n_r_max_3D)

   end subroutine update_temp_3D
!-------------------------------------------------------------------------------
   subroutine finish_exp_temp_3D(dVrT_LMloc, dtemp_exp_last)

      !-- Input variables
      complex(cp), intent(inout) :: dVrT_LMloc(llm:ulm,n_r_max_3D)

      !-- Output variables
      complex(cp), intent(inout) :: dtemp_exp_last(llm:ulm,n_r_max_3D)

      !-- Local variables
      integer :: n_r, lm

      call get_dr( dVrT_LMloc,work_LMloc, llm, ulm, &
           &       n_r_max_3D, rscheme_3D, nocopy=.true. )

      do n_r=1,n_r_max_3D
         do lm=llm,ulm
            dtemp_exp_last(lm,n_r)=         dtemp_exp_last(lm,n_r) &
            &                      -or2_3D(n_r)*work_LMloc(lm,n_r)
         end do
      end do

   end subroutine finish_exp_temp_3D
!-------------------------------------------------------------------------------
   subroutine get_temp_3D_rhs_imp(temp_3D, dtemp_3D, temp_last, &
              &                   dtemp_imp_last, l_calc_lin_rhs)

      !-- Input variables
      complex(cp), intent(in) :: temp_3D(llm:ulm,n_r_max_3D)
      logical,     intent(in) :: l_calc_lin_rhs

      !-- Output variable
      complex(cp), intent(out) :: dtemp_3D(llm:ulm,n_r_max_3D)
      complex(cp), intent(out) :: temp_last(llm:ulm,n_r_max_3D)
      complex(cp), intent(out) :: dtemp_imp_last(llm:ulm,n_r_max_3D)

      !-- Local variables
      integer :: n_r, lm
      integer, pointer :: lm2l(:),lm2m(:)

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      do n_r=1,n_r_max_3D
         do lm=llm,ulm
            temp_last(lm,n_r)=temp_3D(lm,n_r)
         end do
      end do

      if ( l_calc_lin_rhs ) then
         call get_ddr(temp_3D, dtemp_3D, work_LMloc, llm, ulm, &
              &       n_r_max_3D, rscheme_3D)

         !-- Calculate explicit time step part:
         do n_r=1,n_r_max_3D
            do lm=llm,ulm
               dtemp_imp_last(lm,n_r)=TdiffFac* (                           & 
               &                                         work_LMloc(lm,n_r) &
               &        + two*or1_3D(n_r)                * dtemp_3D(lm,n_r) &
               &        - dLh(st_map%lm2(lm2l(lm),lm2m(lm))) * or2_3D(n_r)  &
               &                                         *  temp_3D(lm,n_r) )
            end do
         end do

      end if

   end subroutine get_temp_3D_rhs_imp
!-------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_T0Mat(tscheme,tMat,tPivot,tMat_fac)
#else
   subroutine get_T0Mat(tscheme,tMat,tPivot)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrix   
      !  TMat0                                                            
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step

      !-- Output variables
      real(cp), intent(out) :: tMat(n_r_max_3D,n_r_max_3D)
      integer,  intent(out) :: tPivot(n_r_max_3D)
#ifdef WITH_PRECOND_S
      real(cp), intent(out) :: tMat_fac(n_r_max_3D)
#endif

      !-- Local variables:
      integer :: info,nR_out,nR

      !----- Boundary condition:
      do nR_out=1,rscheme_3D%n_max
    
         if ( ktopt == 1 ) then
            !--------- Constant temperature at CMB:
            tMat(1,nR_out)=rscheme_3D%rnorm*rscheme_3D%rMat(1,nR_out)
         else
            !--------- Constant flux at CMB:
            tMat(1,nR_out)=rscheme_3D%rnorm*rscheme_3D%drMat(1,nR_out)
         end if
         if ( kbott == 1 ) then
            !--------- Constant temperature at ICB:
            tMat(n_r_max_3D,nR_out)=rscheme_3D%rnorm* &
            &                    rscheme_3D%rMat(n_r_max_3D,nR_out)
         else
            !--------- Constant flux at ICB:
            tMat(n_r_max_3D,nR_out)=rscheme_3D%rnorm* &
            &                    rscheme_3D%drMat(n_r_max_3D,nR_out)
         end if
      end do
      if ( rscheme_3D%n_max < n_r_max_3D ) then ! fill with zeros !
         do nR_out=rscheme_3D%n_max+1,n_r_max_3D
            tMat(1,nR_out)         =0.0_cp
            tMat(n_r_max_3D,nR_out)=0.0_cp
         end do
      end if
    
      do nR_out=1,n_r_max_3D
         do nR=2,n_r_max_3D-1
            tMat(nR,nR_out)= rscheme_3D%rnorm * (                         &
            &                                rscheme_3D%rMat(nR,nR_out) - & 
            &     tscheme%wimp_lin(1)*TdiffFac*(                          &
            &                              rscheme_3D%d2rMat(nR,nR_out) + &
            &        two*or1_3D(nR)*        rscheme_3D%drMat(nR,nR_out) ) )
         end do
      end do
    
      !----- Factors for highest and lowest cheb mode:
      do nR=1,n_r_max_3D
         tMat(nR,1)         =rscheme_3D%boundary_fac*tMat(nR,1)
         tMat(nR,n_r_max_3D)=rscheme_3D%boundary_fac*tMat(nR,n_r_max_3D)
      end do
    
#ifdef WITH_PRECOND_S
      ! compute the linesum of each line
      do nR=1,n_r_max_3D
         tMat_fac(nR)=one/maxval(abs(tMat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max_3D
         tMat(nR,:) = tMat(nR,:)*tMat_fac(nR)
      end do
#endif
    
      !---- LU decomposition:
      call prepare_full_mat(tMat,n_r_max_3D,n_r_max_3D,tPivot,info)
      if ( info /= 0 ) then
         call abortRun('! Singular matrix TMat0!')
      end if

   end subroutine get_T0Mat
!-----------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_Tmat(tscheme,l,tMat,tPivot,tMat_fac)
#else
   subroutine get_Tmat(dt,l,tMat,tPivot)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matricies
      !  tMat(i,j) and s0mat for the temperature equation.
      !
      
      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step
      integer,  intent(in) :: l

      !-- Output variables
      real(cp), intent(out) :: tMat(n_r_max_3D,n_r_max_3D)
      integer,  intent(out) :: tPivot(n_r_max_3D)
#ifdef WITH_PRECOND_S
      real(cp),intent(out) :: tMat_fac(n_r_max_3D)
#endif

      !-- Local variables:
      integer :: info,nR_out,nR
      real(cp) :: dLh

      dLh=real(l*(l+1),kind=cp)

      !----- Boundary coditions:
      do nR_out=1,rscheme_3D%n_max
         if ( ktopt == 1 ) then
            tMat(1,nR_out)=rscheme_3D%rnorm*rscheme_3D%rMat(1,nR_out)
         else
            tMat(1,nR_out)=rscheme_3D%rnorm*rscheme_3D%drMat(1,nR_out)
         end if
         if ( kbott == 1 ) then
            tMat(n_r_max_3D,nR_out)=rscheme_3D%rnorm* &
            &                    rscheme_3D%rMat(n_r_max_3D,nR_out)
         else
            tMat(n_r_max_3D,nR_out)=rscheme_3D%rnorm* &
            &                    rscheme_3D%drMat(n_r_max_3D,nR_out)
         end if
      end do
      if ( rscheme_3D%n_max < n_r_max_3D ) then ! fill with zeros !
         do nR_out=rscheme_3D%n_max+1,n_r_max_3D
            tMat(1,nR_out)         =0.0_cp
            tMat(n_r_max_3D,nR_out)=0.0_cp
         end do
      end if

      !----- Other points:
      do nR_out=1,n_r_max_3D
         do nR=2,n_r_max_3D-1
            tMat(nR,nR_out)= rscheme_3D%rnorm * (                        &
            &                               rscheme_3D%rMat(nR,nR_out) - &
            &     tscheme%wimp_lin(1)*TdiffFac*(                         &
            &                             rscheme_3D%d2rMat(nR,nR_out) + &
            &              two*or1_3D(nR)* rscheme_3D%drMat(nR,nR_out) - &
            &      dLh*or2_3D(nR)*          rscheme_3D%rMat(nR,nR_out) ) )
         end do
      end do

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max_3D
         tMat(nR,1)         =rscheme_3D%boundary_fac*tMat(nR,1)
         tMat(nR,n_r_max_3D)=rscheme_3D%boundary_fac*tMat(nR,n_r_max_3D)
      end do

#ifdef WITH_PRECOND_S
      ! compute the linesum of each line
      do nR=1,n_r_max_3D
         tMat_fac(nR)=one/maxval(abs(tMat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max_3D
         tMat(nR,:) = tMat(nR,:)*tMat_fac(nR)
      end do
#endif

      !----- LU decomposition:
      call prepare_full_mat(tMat,n_r_max_3D,n_r_max_3D,tPivot,info)
      if ( info /= 0 ) then
         call abortRun('Singular matrix tMat!')
      end if
            
   end subroutine get_Tmat
!-----------------------------------------------------------------------------
end module update_temp_3D_mod
