module blocking_lm
   !
   !  Module containing blocking information
   !

   use precision_mod
   use mem_alloc, only: memWrite, bytes_allocated
   use parallel_mod, only: rank, n_procs
   use truncation_3D, only: m_max_3D, lmP_max, lm_max, l_max, n_m_max_3D, &
       &                    n_theta_max, minc_3D, n_r_max_3D
   use lm_mapping, only: mappings, allocate_mappings, deallocate_mappings, &
       &                 allocate_subblocks_mappings,                      &
       &                 deallocate_subblocks_mappings,                    &
       &                 subblocks_mappings
   use useful, only: logWrite, abortRun
   use constants, only: one
 
   implicit none
 
   private
 
   !------------------------------------------------------------------------
   !  Dividing loops over all spherical harmonics into blocks that
   !  contain approx. nChunk harmonics . The number of blocks, nLMBs,   
   !  is a multiple of nThreadsUse (the number of processors used).
   !  Parameter nChunk controls the number (and size) of blocks.
   !  When nThreadUse is large, the size of the blocks may be 
   !  considerably smaller than the chosen nChunk,
   !  since nLMBs must be a multiple of nThreadsUse!
 
   integer, public, pointer :: lm2(:,:),lm2l(:),lm2m(:)
   integer, public, pointer :: lm2mc(:),l2lmAS(:)
   integer, public, pointer :: lm2lmS(:),lm2lmA(:)
 
   integer, public, pointer :: lmP2(:,:),lmP2l(:)
   integer, public, pointer :: lmP2lmPS(:),lmP2lmPA(:)
 
   integer, public, pointer :: lm2lmP(:),lmP2lm(:)
   
   integer :: nLMBs_per_rank, rank_with_l1m0
 
   
   type(mappings), public, target :: st_map
   type(mappings), public, target :: lo_map
 
   integer, public :: nLMBs,sizeLMB
 
   integer, public, allocatable :: lmStartB(:),lmStopB(:)
 
   integer, public, pointer :: nLMBs2(:),sizeLMB2(:,:)
   integer, public, pointer :: lm22lm(:,:,:)
   integer, public, pointer :: lm22l(:,:,:)
   integer, public, pointer :: lm22m(:,:,:)
 
   type(subblocks_mappings), public, target :: st_sub_map, lo_sub_map

   !------------------------------------------------------------------------
   !  Following divides loops over points in theta-direction (index ic) into
   !  blocks. Enhances performance by trying to decrease memory access
   !  but is not relevant for SMP parallel processing.
   !  It is tried to divide the theta loop into parts whose data
   !  fit into cache.
   !
 
   integer, public :: n_theta_blocks
   integer, public :: llm, ulm
   integer, parameter, public :: chunksize=16
   integer, public :: lm_per_rank, lm_on_last_rank
 
   public :: initialize_blocking, finalize_blocking, get_lm_blocks

contains

   subroutine initialize_blocking(n_log_file)

      !-- Input variables:
      integer,            intent(in) :: n_log_file

      !-- Local variables:
      real(cp) :: load
      integer :: iLoad
      integer :: n
      integer(lip) :: local_bytes_used
      integer :: LMB_with_l1m0,l1m0,irank
      integer :: n_out

      local_bytes_used = bytes_allocated
      call allocate_mappings(st_map,l_max,lm_max,lmP_max)
      call allocate_mappings(lo_map,l_max,lm_max,lmP_max)

      nLMBs = n_procs
      nLMBs_per_rank = nLMBs/n_procs
      sizeLMB=(lm_max-1)/nLMBs+1

      if ( nLMBs*sizeLMB > lm_max ) then
         load=real(lm_max-(nLMBs-1)*sizeLMB,cp)/sizeLMB
         iLoad=int(load)
         if ( iLoad >= 1 ) then
            nLMBs=nLMBs-iLoad
         end if
      end if
      allocate( lmStartB(nLMBs),lmStopB(nLMBs) )
      bytes_allocated = bytes_allocated+2*nLMBS*SIZEOF_INTEGER

      !--- Calculate lm and ml blocking:
      do n=1,nLMBs
         lmStartB(n)=(n-1)*sizeLMB+1
         lmStopB(n) =min(n*sizeLMB,lm_max)
         if ( lmStopB(n) == lm_max ) exit
      end do

      call get_standard_lm_blocking(st_map,minc_3D)
      call get_snake_lm_blocking(lo_map,minc_3D)

      do n=1,nLMBs
         if ( lmStopB(n) == lm_max ) exit
      end do

      !--- Get the block (rank+1) with the l1m0 mode
      l1m0 = lo_map%lm2(1,0)
      do n=1,nLMBs
         if ( (l1m0 >= lmStartB(n)) .and. (l1m0 <= lmStopB(n)) ) then
            LMB_with_l1m0=n
            exit
         end if
      end do

      ! which rank does have the LMB with LMB_with_l1m0?
      do irank=0,n_procs-1
         if ((LMB_with_l1m0-1 >= irank*nLMBs_per_rank).and.&
              &(LMB_with_l1m0-1 <= (irank+1)*nLMBs_per_rank-1)) then
            rank_with_l1m0 = irank
         end if
      end do
      
      ! set the standard ordering as default
      lm2(0:,0:) => st_map%lm2
      lm2l(1:lm_max) => st_map%lm2l
      lm2m(1:lm_max) => st_map%lm2m
      lm2mc(1:lm_max)=> st_map%lm2mc
      l2lmAS(0:l_max)=> st_map%l2lmAS
      lm2lmS(1:lm_max) => st_map%lm2lmS
      lm2lmA(1:lm_max) => st_map%lm2lmA
      lmP2(0:,0:) => st_map%lmP2
      lmP2l(1:lmP_max) => st_map%lmP2l
      lmP2lmPS(1:lmP_max) => st_map%lmP2lmPS
      lmP2lmPA(1:lmP_max) => st_map%lmP2lmPA
      lm2lmP(1:lm_max) => st_map%lm2lmP
      lmP2lm(1:lmP_max) => st_map%lmP2lm

      call allocate_subblocks_mappings(st_sub_map,st_map,nLMBs,l_max,&
           &                           lmStartB,lmStopB)
      call allocate_subblocks_mappings(lo_sub_map,lo_map,nLMBs,l_max,&
           &                           lmStartB,lmStopB)

      !--- Getting lm sub-blocks:
      !PRINT*," ---------------- Making the standard subblocks -------------- "
      call get_subblocks(st_map, st_sub_map) !nLMBs,lm2l,lm2m, nLMBs2,sizeLMB2,lm22lm,lm22l,lm22m)
      !PRINT*," ---------------- Making the lorder subblocks ---------------- "
      call get_subblocks(lo_map, lo_sub_map)
      !PRINT*," ---------------- Making the snake order subblocks ----------- "

      ! default mapping
      nLMBs2(1:nLMBs) => st_sub_map%nLMBs2
      sizeLMB2(1:,1:) => st_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => st_sub_map%lm22lm
      lm22l(1:,1:,1:) => st_sub_map%lm22l
      lm22m(1:,1:,1:) => st_sub_map%lm22m

      !-- Calculate blocking parameters for blocking loops over theta:
      n_theta_blocks = 1
      if ( rank == 0 ) then
         do n=1,2
            if ( n==1 ) n_out=6
            if ( n==2 ) n_out=n_log_file
            
            write(n_out,*) '!-- Blocking information:'
            write(n_out,*)
            write(n_out,*) '!    Number of LM-blocks:',nLMBs
            write(n_out,*) '!    Size   of LM-blocks:',sizeLMB
            write(n_out,*)
            write(n_out,*) '! Number of theta blocks:',n_theta_blocks
            write(n_out,*) '!   size of theta blocks:',n_theta_max
         end do
      end if

      local_bytes_used = bytes_allocated-local_bytes_used
      call memWrite('blocking_lm.f90', local_bytes_used)

   end subroutine initialize_blocking
!------------------------------------------------------------------------
   subroutine finalize_blocking

      call deallocate_mappings(st_map)
      call deallocate_mappings(lo_map)

      deallocate( lmStartB,lmStopB )

      call deallocate_subblocks_mappings(st_sub_map)
      call deallocate_subblocks_mappings(lo_sub_map)

   end subroutine finalize_blocking
!------------------------------------------------------------------------
   subroutine get_lm_blocks
    
      !-- Local variables:
      integer :: nLMB_start,nLMB_end

      ! set the local lower and upper index for lm
      ! we have nLMBs LM blocks which are distributed over the ranks
      ! with nLMBs_per_rank blocks per rank (computed in m_blocking.F90)
      
      nLMB_start = 1+rank*nLMBs_per_rank
      nLMB_end   = min((rank+1)*nLMBs_per_rank,nLMBs)
      llm = lmStartB(nLMB_start)
      ulm = lmStopB(nLMB_end)
      lm_per_rank=nLMBs_per_rank*sizeLMB
      lm_on_last_rank=lmStopB (min(n_procs*nLMBs_per_rank,nLMBs))- &
                      lmStartB(1+(n_procs-1)*nLMBs_per_rank)+1

   end subroutine get_lm_blocks
!---------------------------------------------------------------------------
   subroutine get_subblocks(map,sub_map) 

      !-- Input variables:
      type(mappings),           intent(in) :: map
      type(subblocks_mappings), intent(inout) :: sub_map
      
      !-- Local variables:
      integer :: number_of_blocks
      logical :: lAfter,lStop
      integer :: size
      integer :: nB2,n,n2,n3
      integer :: help,help1(lm_max),help2(lm_max),help3(lm_max)
      integer :: lm,l,m
      integer :: check(0:l_max,0:l_max)

      number_of_blocks=sub_map%nLMBs
      
      check = 0
      lStop=.false.
      size=0
      nB2=0
      do n=1,nLMBs
         sub_map%nLMBs2(n)=1
         lm=lmStartB(n)
         !------ Start first sub-block:
         sub_map%sizeLMB2(1,n) =1
         sub_map%lm22lm(1,1,n) =lm
         sub_map%lm22l(1,1,n)  =map%lm2l(lm)
         sub_map%lm22m(1,1,n)  =map%lm2m(lm)
         do lm=lmStartB(n)+1,lmStopB(n)
            do n2=1,sub_map%nLMBs2(n)
               if ( sub_map%lm22l(1,n2,n) == map%lm2l(lm) ) then
                  !------ Add to old block
                  sub_map%sizeLMB2(n2,n)=sub_map%sizeLMB2(n2,n)+1
                  lAfter=.false.
                  exit
               else
                  lAfter=.true.
               end if
            end do
            if ( lAfter ) then
               !------ Start new l-block:
               n2 = sub_map%nLMBs2(n)+1
               sub_map%nLMBs2(n)     =n2
               sub_map%sizeLMB2(n2,n)=1
            end if
            sub_map%lm22lm(sub_map%sizeLMB2(n2,n),n2,n)=lm
            sub_map%lm22l(sub_map%sizeLMB2(n2,n),n2,n) =map%lm2l(lm)
            sub_map%lm22m(sub_map%sizeLMB2(n2,n),n2,n) =map%lm2m(lm)
         end do

         !------ Resort:
         if ( sub_map%nLMBs2(n) > 1 ) then
            do n2=1,sub_map%nLMBs2(n)
               do n3=n2+1,sub_map%nLMBs2(n)
                  if  ( sub_map%lm22m(1,n2,n) > sub_map%lm22m(1,n3,n) ) then
                     help=sub_map%sizeLMB2(n2,n)
                     do lm=1,help
                        help1(lm)=sub_map%lm22l(lm,n2,n)
                        help2(lm)=sub_map%lm22m(lm,n2,n)
                        help3(lm)=sub_map%lm22lm(lm,n2,n)
                     end do
                     sub_map%sizeLMB2(n2,n)=sub_map%sizeLMB2(n3,n)
                     do lm=1,sub_map%sizeLMB2(n2,n)
                        sub_map%lm22l(lm,n2,n) =sub_map%lm22l(lm,n3,n)
                        sub_map%lm22m(lm,n2,n) =sub_map%lm22m(lm,n3,n)
                        sub_map%lm22lm(lm,n2,n)=sub_map%lm22lm(lm,n3,n)
                     end do
                     sub_map%sizeLMB2(n3,n)=help
                     do lm=1,help
                        sub_map%lm22l(lm,n3,n) =help1(lm)
                        sub_map%lm22m(lm,n3,n) =help2(lm)
                        sub_map%lm22lm(lm,n3,n)=help3(lm)
                     end do
                  end if
               end do
            end do
         end if

         nB2=nB2+sub_map%nLMBs2(n)
         do n2=1,sub_map%nLMBs2(n)
            if ( sub_map%sizeLMB2(n2,n) > sub_map%sizeLMB2max ) then
               lStop=.true.
               size=max(size,sub_map%sizeLMB2(n2,n))
            end if
            do n3=1,sub_map%sizeLMB2(n2,n)
               l=sub_map%lm22l(n3,n2,n)
               m=sub_map%lm22m(n3,n2,n)
               check(l,m)=check(l,m)+1
            end do
         end do

      end do

      if ( lStop ) then
         write(*,*) '! Increase sizeLMB2max in m_blocking.F90!'
         write(*,*) '! to at least:',size
         call abortRun('Stop run in blocking')
      end if

      do m=0,m_max_3D,minc_3D
         do l=m,l_max
            if ( check(l,m) == 0 ) then
               write(*,*) 'Warning, forgotten l,m:',l,m,map%lm2(l,m)
               call abortRun('Stop run in blocking')
            else if ( check(l,m) > 1 ) then
               write(*,*) 'Warning, too much l,m:',l,m,check(l,m)
               call abortRun('Stop run in blocking')
            end if
         end do
      end do
   end subroutine get_subblocks
!------------------------------------------------------------------------
   subroutine get_standard_lm_blocking(map,minc)

      type(mappings), intent(inout) :: map
      integer,        intent(in) :: minc
      
      ! Local variables
      integer :: m,l,lm,lmP,mc
      
      do m=0,map%l_max
         do l=m,map%l_max
            map%lm2(l,m)  =-1
            map%lmP2(l,m) =-1
         end do
         l=map%l_max+1
         map%lmP2(l,m)=-1
      end do

      lm =0
      lmP=0
      mc =0
      do m=0,map%m_max,minc
         mc=mc+1
         !m2mc(m)=mc
         do l=m,map%l_max
            lm         =lm+1
            map%lm2l(lm)   =l
            map%lm2m(lm)   =m
            map%lm2mc(lm)  =mc
            map%lm2(l,m)   =lm
            if ( m == 0 ) map%l2lmAS(l)=lm
            lmP        =lmP+1
            map%lmP2l(lmP) = l
            map%lmP2m(lmP) = m
            map%lmP2(l,m)  =lmP
            !if ( m == 0 ) l2lmPAS(l)=lmP
            map%lm2lmP(lm) =lmP
            map%lmP2lm(lmP)=lm
         end do
         l=map%l_max+1    ! Extra l for lmP
         lmP=lmP+1
         map%lmP2l(lmP) =l
         map%lmP2m(lmP) = m
         map%lmP2(l,m)  =lmP
         !if ( m == 0 ) l2lmPAS(l)=lmP
         map%lmP2lm(lmP)=-1
      end do
      if ( lm /= map%lm_max ) then
         write(*,"(2(A,I6))") 'Wrong lm=',lm," != map%lm_max = ",map%lm_max
         call abortRun('Stop run in blocking')
      end if
      if ( lmP /= map%lmP_max ) then
         write(*,*) 'Wrong lmP!'
         call abortRun('Stop run in blocking')
      end if
      do lm=1,map%lm_max
         l=map%lm2l(lm)
         m=map%lm2m(lm)
         if ( l > 0 .and. l > m ) then
            map%lm2lmS(lm)=map%lm2(l-1,m)
         else
            map%lm2lmS(lm)=-1
         end if
         if ( l < map%l_max ) then
            map%lm2lmA(lm)=map%lm2(l+1,m)
         else
            map%lm2lmA(lm)=-1
         end if
      end do
      do lmP=1,map%lmP_max
         l=map%lmP2l(lmP)
         m=map%lmP2m(lmP)
         if ( l > 0 .and. l > m ) then
            map%lmP2lmPS(lmP)=map%lmP2(l-1,m)
         else
            map%lmP2lmPS(lmP)=-1
         end if
         if ( l < map%l_max+1 ) then
            map%lmP2lmPA(lmP)=map%lmP2(l+1,m)
         else
            map%lmP2lmPA(lmP)=-1
         end if
      end do

   end subroutine get_standard_lm_blocking
!------------------------------------------------------------------------
   subroutine get_snake_lm_blocking(map,minc)

      type(mappings), intent(inout) :: map
      integer,        intent(in) :: minc
      
      ! Local variables
      integer :: l,proc,lm,m,i_l,lmP,mc
      logical :: Ascending
      integer :: l_list(0:n_procs-1,map%l_max+1)
      integer :: l_counter(0:n_procs-1)
      integer :: temp_l_counter,l0proc,pc,src_proc,temp
      integer :: temp_l_list(map%l_max+1)

      do m=0,map%l_max
         do l=m,map%l_max
            map%lm2(l,m)  =-1
            map%lmP2(l,m) =-1
         end do
         l=map%l_max+1
         map%lmP2(l,m)=-1
      end do

      ! First we loop over all l values and jump for each
      ! new l value to the next process in a snake like fashion.
      proc=0
      Ascending=.true.
      l_counter=1
      do l=map%l_max,0,-1
         ! this l block is distributed to the actual proc
         l_list(proc,l_counter(proc))=l
         l_counter(proc) = l_counter(proc)+1
         if (l == 0) l0proc=proc
         ! now determine on which proc to put the next l value
         if (Ascending) then
            if (proc < n_procs-1) then
               proc=proc+1
            else if (proc == n_procs-1) then
               Ascending=.false.
            end if
         else
            if (proc > 0) then
               proc=proc-1
            else if (proc == 0) then
               Ascending=.true.
            end if
         end if
      end do

      ! Now distribution is as equal as possible. We rotate the distribution
      ! now to have the l0proc as first process. 
      if (l0proc /= 0) then
         temp_l_list=l_list(0,:)
         temp_l_counter=l_counter(0)
         pc = 0
         do while (.true.)
            src_proc=modulo(l0proc+pc,n_procs)
            if (src_proc /= 0) then
               l_list(pc,:)=l_list(src_proc,:)
               l_counter(pc)=l_counter(src_proc)
            else
               l_list(pc,:)=temp_l_list
               l_counter(pc)=temp_l_counter
               exit
            end if
            ! now we can overwrite src_proc
            pc=src_proc
         end do
         l0proc=0
      end if

      ! Last step in preparation is to put the l=0 on process 0
      ! as the first l in the list
      do i_l=1,l_counter(0)-1
         if (l_list(0,i_l) == 0) then
            temp=l_list(0,1)
            l_list(0,1)=l_list(0,i_l)
            l_list(0,i_l)=temp
            exit
         end if
      end do

      lm=1
      lmP=1
      do proc=0,n_procs-1
         lmStartB(proc+1)=lm
         do i_l=1,l_counter(proc)-1
            l=l_list(proc,i_l)
            mc = 0
            do m=0,l,minc
               mc = mc+1
               map%lm2(l,m)=lm
               map%lm2l(lm)=l
               map%lm2m(lm)=m
               map%lm2mc(lm)=mc

               map%lmP2(l,m)=lmP
               map%lmP2l(lmP)=l
               map%lmP2m(lmP)= m
               map%lm2lmP(lm)=lmP
               map%lmP2lm(lmP)=lm

               lm = lm+1
               lmP = lmP+1
            end do
         end do
         lmStopB(proc+1)=lm-1
      end do
      
      if ( lm-1 /= map%lm_max ) then
         write(*,"(2(A,I6))") 'get_snake_lm_blocking: Wrong lm-1 = ',lm-1,&
              & " != map%lm_max = ",map%lm_max
         call abortRun('Stop run in blocking')
      end if

      l=map%l_max+1    ! Extra l for lmP
      mc =0
      do m=0,map%m_max,minc
         mc=mc+1

         map%lmP2l(lmP) =l
         map%lmP2m(lmP) = m
         map%lmP2(l,m)  =lmP
         map%lmP2lm(lmP)=-1
         lmP=lmP+1
      end do

      if ( lmP-1 /= map%lmP_max ) then
         call abortRun('Wrong lmP!')
      end if

      do lm=1,map%lm_max
         l=map%lm2l(lm)
         m=map%lm2m(lm)
         if ( l > 0 .and. l > m ) then
            map%lm2lmS(lm)=map%lm2(l-1,m)
         else
            map%lm2lmS(lm)=-1
         end if
         if ( l < map%l_max ) then
            map%lm2lmA(lm)=map%lm2(l+1,m)
         else
            map%lm2lmA(lm)=-1
         end if
      end do
      do lmP=1,map%lmP_max
         l=map%lmP2l(lmP)
         m=map%lmP2m(lmP)
         if ( l > 0 .and. l > m ) then
            map%lmP2lmPS(lmP)=map%lmP2(l-1,m)
         else
            map%lmP2lmPS(lmP)=-1
         end if
         if ( l < map%l_max+1 ) then
            map%lmP2lmPA(lmP)=map%lmP2(l+1,m)
         else
            map%lmP2lmPA(lmP)=-1
         end if
      end do

   end subroutine get_snake_lm_blocking
!------------------------------------------------------------------------
end module blocking_lm
