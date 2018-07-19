module radial_der
   !
   ! Radial derivatives functions
   !

   use constants, only: zero, one, three
   use parallel_mod, only: n_threads
   use precision_mod
   use mem_alloc
   use radial_scheme, only: type_rscheme

   implicit none

   private

   interface get_dcheb
      module procedure get_dcheb_real_1d
      module procedure get_dcheb_complex_2d
   end interface get_dcheb

   interface get_ddcheb
      module procedure get_ddcheb_real_1d
      module procedure get_ddcheb_complex_2d
   end interface get_ddcheb

   interface get_dr
      module procedure get_dr_real_1d
      module procedure get_dr_complex_2d
   end interface get_dr

   interface get_ddr
      module procedure get_ddr_real_1d
      module procedure get_ddr_complex_2d
   end interface get_ddr

   public :: get_ddr, get_dcheb, get_dr, initialize_der_arrays, &
   &         finalize_der_arrays

   complex(cp), allocatable :: work(:,:)
   real(cp), allocatable :: work_1d_real(:)
   real(cp) :: thr

contains

!------------------------------------------------------------------------------
   subroutine initialize_der_arrays(n_r_max, nMstart, nMstop, l_rerror, rerror_fac)
      !
      ! Allocate work arrays to compute derivatives
      !

      integer,  intent(in) :: n_r_max
      integer,  intent(in) :: nMstart
      integer,  intent(in) :: nMstop
      logical,  intent(in) :: l_rerror
      real(cp), intent(in) :: rerror_fac

      if ( l_rerror ) then
         thr = rerror_fac * epsilon(1.0_cp)
      else
         thr = 0.0_cp
      end if

      allocate( work_1d_real(n_r_max) )
      allocate( work(nMstart:nMstop,n_r_max) )
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL+&
      &                 n_r_max*(nMstop-nMstart+1)*SIZEOF_DEF_COMPLEX

   end subroutine initialize_der_arrays
!------------------------------------------------------------------------------
   subroutine finalize_der_arrays
      !
      ! Deallocate work arrays
      !
      deallocate( work_1d_real, work )

   end subroutine finalize_der_arrays
!------------------------------------------------------------------------------
   subroutine get_dcheb_complex_2d(f,df,nMstart,nMstop,n_r_max,n_cheb_max)
      !
      !  Returns chebychev coeffitiens of first derivative df and second  
      !  derivative ddf for a function whose cheb-coeff. are given as     
      !  columns in array f(nMstart:nMstop,n_r_max).                             
      !

      !-- Input variables:
      integer,     intent(in) :: nMstart  ! No of function to start with
      integer,     intent(in) :: nMstop   ! No of function to stop with
      integer,     intent(in) :: n_r_max    ! second dimension of f,df,ddf
      integer,     intent(in) :: n_cheb_max ! Number of cheb modes
      complex(cp), intent(in) :: f(nMstart:nMstop,n_r_max)

      !-- Output variables:
      complex(cp), intent(out) ::  df(nMstart:nMstop,n_r_max)

      !-- Local variables:
      integer :: n_m,n_cheb
      real(cp) :: fac_cheb
      real(cp) :: threshold(nMstart:nMstop)

      do n_m=nMstart, nMstop
         threshold(n_m) = maxval(abs(f(n_m,:)))*thr
      end do

      !-- initialize derivatives:
      do n_cheb=n_cheb_max,n_r_max
         do n_m=nMstart,nMstop
            df(n_m,n_cheb)=zero
         end do
      end do
      n_cheb  =n_cheb_max-1
      if ( n_r_max == n_cheb_max ) then
         fac_cheb=real(n_cheb,kind=cp)
      else
         fac_cheb=real(2*n_cheb,kind=cp)
      end if
      do n_m=nMstart,nMstop
         if ( abs(df(n_m,n_cheb)) >= threshold(n_m) ) then
            df(n_m,n_cheb)=fac_cheb*f(n_m,n_cheb+1)
         else
            df(n_m,n_cheb)=zero
         end if
      end do

      !----- Recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=real(2*n_cheb,kind=cp)
         do n_m=nMstart,nMstop
            if ( abs(f(n_m,n_cheb+1)) >= threshold(n_m) ) then
               df(n_m,n_cheb)=df(n_m,n_cheb+2) + fac_cheb*f(n_m,n_cheb+1)
            else
               df(n_m,n_cheb)=df(n_m,n_cheb+2)
            end if
         end do
      end do

   end subroutine get_dcheb_complex_2d
!------------------------------------------------------------------------------
   subroutine get_dcheb_real_1d(f,df,n_r_max,n_cheb_max)

      !-- Input variables:
      integer,  intent(in) :: n_r_max    ! second dimension of f,df,ddf
      integer,  intent(in) :: n_cheb_max ! Number of cheb modes
      real(cp), intent(in) :: f(n_r_max)

      !-- Output variables:
      real(cp), intent(out) ::  df(n_r_max)

      !-- Local variables:
      integer :: n_cheb
      real(cp) :: fac_cheb, threshold

      threshold = maxval(abs(f))*thr

      !-- initialize derivatives:
      do n_cheb=n_cheb_max,n_r_max
         df(n_cheb)=0.0_cp
      end do
      n_cheb  =n_cheb_max-1
      if ( n_r_max == n_cheb_max ) then
         fac_cheb=real(n_cheb,kind=cp)
      else
         fac_cheb=real(2*n_cheb,kind=cp)
      end if
      if ( abs(f(n_cheb+1)) >= threshold ) then
         df(n_cheb)=fac_cheb*f(n_cheb+1)
      else
         df(n_cheb)=0.0_cp
      end if

      !----- Recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=real(2*n_cheb,kind=cp)
         if ( abs(f(n_cheb+1)) >= threshold ) then
            df(n_cheb)=df(n_cheb+2) + fac_cheb*f(n_cheb+1)
         else
            df(n_cheb)=0.0_cp
         end if
      end do

   end subroutine get_dcheb_real_1d
!------------------------------------------------------------------------------
   subroutine get_ddcheb_complex_2d(f,df,ddf,nMstart,nMstop,n_r_max,n_cheb_max)
      !
      !  Returns chebychev coefficents of first derivative df and second  
      !  derivative ddf for a function whose cheb-coeff. are given as     
      !  columns in array f(n_c_tot,n_r_max).                             
      !
    
      !-- Input variables:
      integer,     intent(in) :: nMstart  ! No of column to start with
      integer,     intent(in) :: nMstop   ! No of column to stop with
      integer,     intent(in) :: n_r_max    ! second dimension of f,df,ddf
      integer,     intent(in) :: n_cheb_max ! Number of cheb modes
      complex(cp), intent(in) :: f(nMstart:nMstop,n_r_max)
    
      !-- Output variables:
      complex(cp), intent(out) ::  df(nMstart:nMstop,n_r_max)
      complex(cp), intent(out) ::  ddf(nMstart:nMstop,n_r_max)
    
      !-- local variables:
      integer :: n_m,n_cheb
      real(cp) :: fac_cheb
      real(cp) :: threshold(nMstart:nMstop)

      do n_m=nMstart, nMstop
         threshold(n_m) = maxval(abs(f(n_m,:)))*thr
      end do

      !----- initialize derivatives:
      do n_cheb=n_cheb_max,n_r_max
         do n_m=nMstart,nMstop
            df(n_m,n_cheb) =zero
            ddf(n_m,n_cheb)=zero
         end do
      end do
      n_cheb=n_cheb_max-1
      if ( n_cheb_max == n_r_max ) then
         fac_cheb=real(n_cheb,kind=cp)
      else
         fac_cheb=real(2*n_cheb,kind=cp)
      end if
      do n_m=nMstart,nMstop
         if ( abs(df(n_m,n_cheb)) >= threshold(n_m) ) then
            df(n_m,n_cheb)=fac_cheb*f(n_m,n_cheb+1)
         else
            df(n_m,n_cheb)=zero
         end if
         ddf(n_m,n_cheb)=zero
      end do
    
      !----- recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=real(2*n_cheb,kind=cp)
         do n_m=nMstart,nMstop
            if ( abs(f(n_m,n_cheb+1)) >= threshold(n_m) ) then
               df(n_m,n_cheb) = df(n_m,n_cheb+2) + fac_cheb* f(n_m,n_cheb+1)
            else
               df(n_m,n_cheb) = df(n_m,n_cheb+2)
            end if
            if ( abs(df(n_m,n_cheb+1)) >= threshold(n_m) ) then
               ddf(n_m,n_cheb)=ddf(n_m,n_cheb+2) + fac_cheb*df(n_m,n_cheb+1)
            else
               ddf(n_m,n_cheb)=ddf(n_m,n_cheb+2)
            end if
         end do
      end do

   end subroutine get_ddcheb_complex_2d
!------------------------------------------------------------------------------
   subroutine get_ddcheb_real_1d(f,df,ddf,n_r_max,n_cheb_max)
      !
      !  Returns chebychev coefficents of first derivative df and second  
      !  derivative ddf for a function whose cheb-coeff. are given as     
      !  columns in array f(n_c_tot,n_r_max).                             
      !
    
      !-- Input variables:
      integer,  intent(in) :: n_r_max    ! second dimension of f,df,ddf
      integer,  intent(in) :: n_cheb_max ! Number of cheb modes
      real(cp), intent(in) :: f(n_r_max)
    
      !-- Output variables:
      real(cp), intent(out) ::  df(n_r_max)
      real(cp), intent(out) ::  ddf(n_r_max)
    
      !-- local variables:
      integer :: n_cheb
      real(cp) :: fac_cheb, threshold

      threshold = maxval(abs(f))*thr
    
      !----- initialize derivatives:
      do n_cheb=n_cheb_max,n_r_max
         df(n_cheb) =0.0_cp
         ddf(n_cheb)=0.0_cp
      end do
      n_cheb=n_cheb_max-1
      if ( n_cheb_max == n_r_max ) then
         fac_cheb=real(n_cheb,kind=cp)
      else
         fac_cheb=real(2*n_cheb,kind=cp)
      end if
      if ( abs(f(n_cheb+1)) >= threshold ) then
         df(n_cheb) =fac_cheb*f(n_cheb+1)
      else
         df(n_cheb) =0.0_cp
      end if
      ddf(n_cheb)=0.0_cp
    
      !----- recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=real(2*n_cheb,kind=cp)
         if ( abs(f(n_cheb+1)) >= threshold ) then
            df(n_cheb) = df(n_cheb+2) + fac_cheb* f(n_cheb+1)
         else
            df(n_cheb) = df(n_cheb+2)
         end if

         if ( abs(df(n_cheb+1)) >= threshold ) then
            ddf(n_cheb)=ddf(n_cheb+2) + fac_cheb*df(n_cheb+1)
         else
            ddf(n_cheb)=ddf(n_cheb+2)
         end if
      end do

   end subroutine get_ddcheb_real_1d
!------------------------------------------------------------------------------
   subroutine get_dr_complex_2d(f,df,nMstart,nMstop,n_r_max,r_scheme,nocopy,  &
              &                 l_dct_in,l_dct_out)
      !
      !  Returns first radial derivative df of the input function f.      
      !  Array f(nMstart:nMstop,*) may contain several functions numbered by     
      !  the first index. The subroutine calculates the derivaties of     
      !  the functions f(nMstart,*) to f(nMstop) by transforming      
      !  to a Chebychev representation using n_r_max radial grid points . 
      !
    
      !-- Input variables:
      integer,             intent(in) :: n_r_max  ! number of radial grid points
      integer,             intent(in) :: nMstart  ! first function to be treated
      integer,             intent(in) :: nMstop   ! last function to be treated
      complex(cp),         intent(inout) :: f(nMstart:nMstop,n_r_max)
      class(type_rscheme), intent(in) :: r_scheme
      logical, optional,   intent(in) :: nocopy
      logical, optional,   intent(in) :: l_dct_in
      logical, optional,   intent(in) :: l_dct_out
    
      !-- Output variables:
      complex(cp), intent(out) :: df(nMstart:nMstop,n_r_max) ! first derivative of f
    
      !-- Local:
      integer :: iThread, start_m, stop_m, per_thread, all_ms
      integer :: n_r,n_f,od
      logical :: copy_array, l_dct_in_loc, l_dct_out_loc

      if ( present(l_dct_in) ) then
         l_dct_in_loc=l_dct_in
      else
         l_dct_in_loc=.true.
      end if

      if ( present(l_dct_out) ) then
         l_dct_out_loc=l_dct_out
      else
         l_dct_out_loc=.true.
      end if
    
      if ( r_scheme%version == 'cheb' ) then

         if ( present(nocopy) ) then
            copy_array=.false.
         else
            copy_array=.true.
         end if
    
         if ( copy_array )  then
            !$omp parallel do default(shared) &
            !$omp private(n_r,n_f)
            do n_r=1,n_r_max
               do n_f=nMstart,nMstop
                  work(n_f,n_r)=f(n_f,n_r)
               end do
            end do
            !$omp end parallel do
       
            !-- Transform f to cheb space:
            if ( l_dct_in_loc ) call r_scheme%costf1(work,nMstart,nMstop,n_r_max)
          
            !-- Get derivatives:
            !call get_dcheb(work,df,nMstart,nMstop,n_r_max,r_scheme%n_max)
            !$omp parallel default(shared) &
            !$omp private(iThread, start_m, stop_m)
            all_ms = nMstop-nMstart+1
            per_thread = all_ms/n_threads
            !$omp do
            do iThread=0,n_threads-1
               start_m=nMstart+iThread*per_thread
               stop_m = start_m+per_thread-1
               if ( iThread == n_threads-1 ) stop_m = nMstop
               call get_dcheb(work(start_m:stop_m,:),df(start_m:stop_m,:), &
                    &              start_m,stop_m,n_r_max,n_r_max)
            end do
            !$omp end do
            !$omp end parallel
          
            !-- Transform back:
            if ( l_dct_out_loc ) call r_scheme%costf1(df,nMstart,nMstop,n_r_max)

         else

            !-- Transform f to cheb space:
            if ( l_dct_in_loc ) call r_scheme%costf1(f,nMstart,nMstop,n_r_max)
          
            !-- Get derivatives:
            !call get_dcheb(f,df,nMstart,nMstop,n_r_max,r_scheme%n_max)
            !$omp parallel default(shared) &
            !$omp private(iThread, start_m, stop_m)
            all_ms = nMstop-nMstart+1
            per_thread = all_ms/n_threads
            !$omp do
            do iThread=0,n_threads-1
               start_m=nMstart+iThread*per_thread
               stop_m = start_m+per_thread-1
               if ( iThread == n_threads-1 ) stop_m = nMstop
               call get_dcheb(f(start_m:stop_m,:),df(start_m:stop_m,:), start_m, &
                    &           stop_m,n_r_max,n_r_max)
            end do
            !$omp end do
            !$omp end parallel
          
            !-- Transform back:
            if ( l_dct_in_loc ) call r_scheme%costf1(f,nMstart,nMstop,n_r_max)
            if ( l_dct_out_loc ) call r_scheme%costf1(df,nMstart,nMstop,n_r_max)

         end if
       
         !-- New map:
         !$omp parallel do default(shared) &
         !$omp private(n_r,n_f)
         do n_r=1,n_r_max
            do n_f=nMstart,nMstop
               df(n_f,n_r)=r_scheme%drx(n_r)*df(n_f,n_r)
            end do
         end do
         !$omp end parallel do

      else

         !-- Initialise to zero:
         do n_r=1,n_r_max
            do n_f=nMstart,nMstop
               df(n_f,n_r) =zero
            end do
         end do

         !-- Bulk points for 1st derivative
         do od=0,r_scheme%order
            do n_r=1+r_scheme%order/2,n_r_max-r_scheme%order/2
               do n_f=nMstart,nMstop
                  df(n_f,n_r)=df(n_f,n_r)+r_scheme%dr(n_r,od)*f(n_f,n_r-r_scheme%order/2+od)
               end do
            end do
         end do

         !-- Boundary points for 1st derivative
         do od=0,r_scheme%order_boundary
            do n_r=1,r_scheme%order/2
               do n_f=nMstart,nMstop
                  df(n_f,n_r) = df(n_f,n_r)+r_scheme%dr_top(n_r,od) * f(n_f,od+1)
               end do
            end do
            do n_r=1,r_scheme%order/2
               do n_f=nMstart,nMstop
                  df(n_f,n_r_max-n_r+1) = df(n_f,n_r_max-n_r+1)+               &
                  &                       r_scheme%dr_bot(n_r,od)*f(n_f,n_r_max-od)
               end do
            end do
         end do

      end if

   end subroutine get_dr_complex_2d
!------------------------------------------------------------------------------
   subroutine get_dr_real_1d(f,df,n_r_max,r_scheme,l_dct)
    
      !-- Input variables:
      integer,             intent(in) :: n_r_max  ! number of radial grid points
      real(cp),            intent(inout) :: f(n_r_max)
      class(type_rscheme), intent(in) :: r_scheme
      logical, optional,   intent(in) :: l_dct
    
      !-- Output variables:
      real(cp), intent(out) :: df(n_r_max) ! first derivative of f
    
      !-- Local:
      integer :: n_r,od
      logical :: l_dct_loc

      if ( present(l_dct) ) then
         l_dct_loc=l_dct
      else
         l_dct_loc=.true.
      end if
    
      if ( r_scheme%version == 'cheb' ) then

         do n_r=1,n_r_max
            work_1d_real(n_r)=f(n_r)
         end do
       
         !-- Transform f to cheb space:
         if ( l_dct_loc ) call r_scheme%costf1(work_1d_real,n_r_max)
          
         !-- Get derivatives:
         !call get_dcheb(work_1d_real,df,n_r_max,r_scheme%n_max)
         call get_dcheb(work_1d_real,df,n_r_max,n_r_max)
          
         !-- Transform back:
         call r_scheme%costf1(df,n_r_max)

         !-- New map:
         do n_r=1,n_r_max
            df(n_r)=r_scheme%drx(n_r)*df(n_r)
         end do

      else

         !-- Initialise to zero:
         do n_r=1,n_r_max
            df(n_r)=0.0_cp
         end do

         !-- Bulk points for 1st derivative
         do od=0,r_scheme%order
            do n_r=1+r_scheme%order/2,n_r_max-r_scheme%order/2
               df(n_r)=df(n_r)+r_scheme%dr(n_r,od)*f(n_r-r_scheme%order/2+od)
            end do
         end do

         !-- Boundary points for 1st derivative
         do od=0,r_scheme%order_boundary
            do n_r=1,r_scheme%order/2
               df(n_r) = df(n_r)+r_scheme%dr_top(n_r,od) * f(od+1)
            end do
            do n_r=1,r_scheme%order/2
               df(n_r_max-n_r+1) = df(n_r_max-n_r+1)+               &
               &                   r_scheme%dr_bot(n_r,od)*f(n_r_max-od)
            end do
         end do

      end if

   end subroutine get_dr_real_1d
!------------------------------------------------------------------------------
   subroutine get_ddr_complex_2d(f,df,ddf,nMstart,nMstop,n_r_max,r_scheme,l_dct)
      !
      !  Returns first radial derivative df and second radial             
      !  derivative ddf of the input function f.                          
      !  Array f(nMstart:nMstop,*) may contain several functions numbered by     
      !  the first index. The subroutine calculates the derivatives of    
      !  the functions f(nMstart,*) to f(nMstop) by transforming      
      !  to a Chebychev representation using n_r_max radial grid points.  
      !
    
      !-- Input variables:
      integer,             intent(in) :: n_r_max  ! number of radial grid points
      integer,             intent(in) :: nMstart  ! first function to be treated
      integer,             intent(in) :: nMstop   ! last function to be treated
      complex(cp),         intent(in) :: f(nMstart:nMstop,n_r_max)
      class(type_rscheme), intent(in) :: r_scheme
      logical, optional,   intent(in) :: l_dct
    
      !-- Output variables:
      complex(cp), intent(out) :: df(nMstart:nMstop,n_r_max) ! first derivative of f
      complex(cp), intent(out) :: ddf(nMstart:nMstop,n_r_max)! second derivative of f
    
      !-- Local variables:
      logical :: l_dct_loc
      integer :: iThread, start_m, stop_m, per_thread, all_ms
      integer :: n_r,n_f,od

      if ( present(l_dct) ) then
         l_dct_loc=l_dct
      else
         l_dct_loc=.true.
      end if

      if ( r_scheme%version == 'cheb' ) then
    
         !-- Copy input functions:
         !$omp parallel do default(shared) &
         !$omp private(n_r,n_f)
         do n_r=1,n_r_max
            do n_f=nMstart,nMstop
               work(n_f,n_r)=f(n_f,n_r)
            end do
         end do
         !$omp end parallel do
    
         !-- Transform f to cheb space:
         if ( l_dct_loc ) call r_scheme%costf1(work,nMstart,nMstop,n_r_max)
    
         !-- Get derivatives:
         !call get_ddcheb(work,df,ddf,nMstart,nMstop,n_r_max,r_scheme%n_max)
         !$omp parallel default(shared) &
         !$omp private(iThread, start_m, stop_m)
         all_ms = nMstop-nMstart+1
         per_thread = all_ms/n_threads
         !$omp do
         do iThread=0,n_threads-1
            start_m=nMstart+iThread*per_thread
            stop_m = start_m+per_thread-1
            if ( iThread == n_threads-1 ) stop_m = nMstop
            call get_ddcheb(work(start_m:stop_m,:),df(start_m:stop_m,:),  &
                 &          ddf(start_m:stop_m,:),start_m,stop_m,n_r_max,n_r_max)
         end do
         !$omp end do
         !$omp end parallel
    
         !-- Transform back:
         call r_scheme%costf1(df,nMstart,nMstop,n_r_max)
         call r_scheme%costf1(ddf,nMstart,nMstop,n_r_max)
    
         !-- New map:
         !$omp parallel do default(shared) &
         !$omp private(n_r,n_f)
         do n_r=1,n_r_max
            do n_f=nMstart,nMstop
               ddf(n_f,n_r)=r_scheme%ddrx(n_r)*df(n_f,n_r)+&
               &            r_scheme%drx(n_r)*r_scheme%drx(n_r)*ddf(n_f,n_r)
               df(n_f,n_r) =r_scheme%drx(n_r)*df(n_f,n_r)
            end do
         end do
         !$omp end parallel do

      else

         !-- Initialise to zero:
         do n_r=1,n_r_max
            do n_f=nMstart,nMstop
               df(n_f,n_r) =zero
               ddf(n_f,n_r)=zero
            end do
         end do

         !-- Bulk points for 1st and 2nd derivatives
         do od=0,r_scheme%order
            do n_r=1+r_scheme%order/2,n_r_max-r_scheme%order/2
               do n_f=nMstart,nMstop
                  df(n_f,n_r)  = df(n_f,n_r) + r_scheme%dr(n_r,od) * f(n_f,n_r-r_scheme%order/2+od)
                  ddf(n_f,n_r) = ddf(n_f,n_r)+r_scheme%ddr(n_r,od) * f(n_f,n_r-r_scheme%order/2+od)
               end do
            end do
         end do

         !-- Boundary points for 1st derivative
         do od=0,r_scheme%order_boundary
            do n_r=1,r_scheme%order/2
               do n_f=nMstart,nMstop
                  df(n_f,n_r) = df(n_f,n_r)+r_scheme%dr_top(n_r,od) * f(n_f,od+1)
               end do
            end do
            do n_r=1,r_scheme%order/2
               do n_f=nMstart,nMstop
                  df(n_f,n_r_max-n_r+1) = df(n_f,n_r_max-n_r+1)+               &
                  &                       r_scheme%dr_bot(n_r,od)*f(n_f,n_r_max-od)
               end do
            end do
         end do

         !-- Boundary points for 2nd derivative
         do od=0,r_scheme%order_boundary+1
            do n_r=1,r_scheme%order/2
               do n_f=nMstart,nMstop
                  ddf(n_f,n_r) = ddf(n_f,n_r)+r_scheme%ddr_top(n_r,od) * f(n_f,od+1)
               end do
            end do
            do n_r=1,r_scheme%order/2
               do n_f=nMstart,nMstop
                  ddf(n_f,n_r_max-n_r+1) = ddf(n_f,n_r_max-n_r+1)+               &
                  &                       r_scheme%ddr_bot(n_r,od)*f(n_f,n_r_max-od)
               end do
            end do
         end do

      end if

   end subroutine get_ddr_complex_2d
!------------------------------------------------------------------------------
   subroutine get_ddr_real_1d(f,df,ddf,n_r_max,r_scheme)
    
      !-- Input variables:
      integer,             intent(in) :: n_r_max  ! number of radial grid points
      real(cp),            intent(in) :: f(n_r_max)
      class(type_rscheme), intent(in) :: r_scheme
    
      !-- Output variables:
      real(cp), intent(out) :: df(n_r_max) ! first derivative of f
      real(cp), intent(out) :: ddf(n_r_max)! second derivative of f
    
      !-- Local variables:
      integer :: n_r,od

      if ( r_scheme%version == 'cheb' ) then
    
         !-- Copy input functions:
         do n_r=1,n_r_max
            work_1d_real(n_r)=f(n_r)
         end do
    
         !-- Transform f to cheb space:
         call r_scheme%costf1(work_1d_real,n_r_max)
    
         !-- Get derivatives:
         !call get_ddcheb(work_1d_real,df,ddf,n_r_max,r_scheme%n_max)
         call get_ddcheb(work_1d_real,df,ddf,n_r_max,n_r_max)
    
         !-- Transform back:
         call r_scheme%costf1(df,n_r_max)
         call r_scheme%costf1(ddf,n_r_max)
    
         !-- New map:
         do n_r=1,n_r_max
            ddf(n_r)=r_scheme%ddrx(n_r)*df(n_r)+&
            &         r_scheme%drx(n_r)*r_scheme%drx(n_r)*ddf(n_r)
            df(n_r) =r_scheme%drx(n_r)*df(n_r)
         end do

      else

         !-- Initialise to zero:
         do n_r=1,n_r_max
            df(n_r) =0.0_cp
            ddf(n_r)=0.0_cp
         end do

         !-- Bulk points for 1st and 2nd derivatives
         do od=0,r_scheme%order
            do n_r=1+r_scheme%order/2,n_r_max-r_scheme%order/2
               df(n_r)  = df(n_r) + r_scheme%dr(n_r,od) * f(n_r-r_scheme%order/2+od)
               ddf(n_r) = ddf(n_r)+r_scheme%ddr(n_r,od) * f(n_r-r_scheme%order/2+od)
            end do
         end do

         !-- Boundary points for 1st derivative
         do od=0,r_scheme%order_boundary
            do n_r=1,r_scheme%order/2
               df(n_r) = df(n_r)+r_scheme%dr_top(n_r,od) * f(od+1)
            end do
            do n_r=1,r_scheme%order/2
               df(n_r_max-n_r+1) = df(n_r_max-n_r+1)+               &
               &                   r_scheme%dr_bot(n_r,od)*f(n_r_max-od)
            end do
         end do

         !-- Boundary points for 2nd derivative
         do od=0,r_scheme%order_boundary+1
            do n_r=1,r_scheme%order/2
               ddf(n_r) = ddf(n_r)+r_scheme%ddr_top(n_r,od) * f(od+1)
            end do
            do n_r=1,r_scheme%order/2
               ddf(n_r_max-n_r+1) = ddf(n_r_max-n_r+1)+               &
               &                    r_scheme%ddr_bot(n_r,od)*f(n_r_max-od)
            end do
         end do

      end if

   end subroutine get_ddr_real_1d
!------------------------------------------------------------------------------
end module radial_der
