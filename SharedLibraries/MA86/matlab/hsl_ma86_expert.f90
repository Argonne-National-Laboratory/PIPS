!
! This routine offers direct access to hsl_ma86 functionality
!     [handle, info] = ma86_expert('factor', A[, control, P])
!     [x, info] = ma86_expert('solve', handle, b[, control])
!     [x, info, handle] = ma86_expert('backslash', A, b[, control, P])
!     ma86_expert('destroy', handle)
! P is an optional ordering.
!

! Structure of code:
! The module ma86_matlab_main provides a series of functions, one for each
!    possible 'action' in real or complex data.
! The module ma86_handles provides a series of functions for dealing with
!    integer handles associated on the Fortran side with keep and order
! The mexFunction unwraps handles and calls the corrected subroutine
!    for each 'action'.

module ma86_matlab_main
!$ use omp_lib
   use hsl_matlab
   use hsl_mc68_double
   use hsl_mc69_double, only: mc69_cscl_clean_real => mc69_cscl_clean, &
                              HSL_MATRIX_REAL_SYM_INDEF, &
                              HSL_MATRIX_CPLX_SYM, &
                              HSL_MATRIX_CPLX_HERM_INDEF
   use hsl_mc69_double_complex, only: mc69_cscl_clean_complex => mc69_cscl_clean
   use hsl_ma86_double, only: ma86_rkeep => ma86_keep,         &
                              ma86_rcontrol => ma86_control,   &
                              ma86_rinfo => ma86_info,         &
                              ma86_analyse,                    &
                              ma86_factor,                     &
                              ma86_factor_solve,               &
                              ma86_solve,                      &
                              ma86_finalise,                   &
                              ma86_get_n__
   use hsl_ma86_double_complex, only: ma86_ckeep => ma86_keep, &
                              ma86_ccontrol => ma86_control,   &
                              ma86_cinfo => ma86_info,         &
                              ma86_analyse,                    &
                              ma86_factor,                     &
                              ma86_factor_solve,               &
                              ma86_solve,                      &
                              ma86_finalise,                   &
                              ma86_get_n__
   implicit none

   integer, parameter :: wp = kind(0d0)
   integer, parameter :: long = selected_int_kind(18)

   interface ma86_matlab_finalise
      module procedure ma86_matlab_finalise_real
      module procedure ma86_matlab_finalise_complex
   end interface ma86_matlab_finalise
   interface copy_control_in
      module procedure copy_control_in_real
      module procedure copy_control_in_complex
   end interface copy_control_in
   interface find_scaling
      module procedure find_scaling_real
      module procedure find_scaling_complex
   end interface find_scaling

contains
   subroutine copy_control_in_real(pm, control, num_threads, scalopt)
      integer(mwp_) :: pm
      type(ma86_rcontrol), intent(inout) :: control
      integer, intent(out) :: num_threads
      integer, intent(out) :: scalopt

      integer(mwp_) :: pc
      integer(int4_) :: fnum
      character(80) :: fname
      character(200) :: warnmsg

      num_threads = 0
      scalopt = 1 ! mc77_1

      do fnum = 1, MATLAB_get_no_fields(pm)
         fname = MATLAB_get_field_name_by_no(pm, fnum)
         select case(trim(fname))
         case("nemin")
            call MATLAB_get_value(pm, fname, pc, control%nemin)
         case("nb")
            call MATLAB_get_value(pm, fname, pc, control%nb)
         case("small")
            call MATLAB_get_value(pm, fname, pc, control%small)
         case("static")
            call MATLAB_get_value(pm, fname, pc, control%static)
         case("u")
            call MATLAB_get_value(pm, fname, pc, control%u)
         case("umin")
            call MATLAB_get_value(pm, fname, pc, control%umin)
         case("num_threads")
            call MATLAB_get_value(pm, fname, pc, num_threads)
         case("scaling")
            call MATLAB_get_value(pm, fname, pc, scalopt)
         case default
            write(warnmsg, "(3a)") "Ignored unrecognised control parameter '", &
               trim(fname), "'"
            call MATLAB_warning(warnmsg)
         end select
      end do
   end subroutine copy_control_in_real

   subroutine copy_control_in_complex(pm, control, num_threads, matrix_type, &
         scalopt)
      integer(mwp_) :: pm
      type(ma86_ccontrol), intent(inout) :: control
      integer, intent(out) :: num_threads
      integer, intent(out) :: matrix_type
      integer, intent(out) :: scalopt

      integer(mwp_) :: pc
      integer(int4_) :: fnum
      character(80) :: fname
      character(200) :: warnmsg
      logical :: herm

      num_threads = 0
      matrix_type = HSL_MATRIX_CPLX_SYM

      do fnum = 1, MATLAB_get_no_fields(pm)
         fname = MATLAB_get_field_name_by_no(pm, fnum)
         select case(trim(fname))
         case("nemin")
            call MATLAB_get_value(pm, fname, pc, control%nemin)
         case("nb")
            call MATLAB_get_value(pm, fname, pc, control%nb)
         case("small")
            call MATLAB_get_value(pm, fname, pc, control%small)
         case("static")
            call MATLAB_get_value(pm, fname, pc, control%static)
         case("u")
            call MATLAB_get_value(pm, fname, pc, control%u)
         case("umin")
            call MATLAB_get_value(pm, fname, pc, control%umin)
         case("num_threads")
            call MATLAB_get_value(pm, fname, pc, num_threads)
         case("hermitian")
            call MATLAB_get_value(pm, fname, pc, herm)
            if(herm) matrix_type = HSL_MATRIX_CPLX_HERM_INDEF
         case("scaling")
            call MATLAB_get_value(pm, fname, pc, scalopt)
         case default
            write(warnmsg, "(3a)") "Ignored unrecognised control parameter '", &
               trim(fname), "'"
            call MATLAB_warning(warnmsg)
         end select
      end do

   end subroutine copy_control_in_complex

   subroutine unknownError(routine, flag)
      character(len=*) :: routine
      integer :: flag
      character(len=200) :: errmsg

      write(errmsg, "(3a,i3)") "Error return from ", routine, ". flag = ", flag
      call MATLAB_error(errmsg)
   end subroutine unknownError

   ! [handle, info] = ma86_expert('factor', A[, P])
   ! handle and 'factor' already dealt with
   subroutine ma86_matlab_analyse_factor_real(nlhs_in, plhs, nrhs_in, prhs, &
         keep, order, scaling)
      integer*4 :: nlhs_in, nrhs_in
      integer(mwp_) :: plhs(*), prhs(*)
      type(ma86_rkeep), intent(inout) :: keep
      integer, dimension(:), allocatable, intent(inout) :: order
      real(wp), dimension(:), allocatable, intent(inout) :: scaling

      integer, parameter :: A_in = 1, &
                            control_in = 2, &
                            P_in = 3
      integer, parameter :: info_out = 1

      integer(mws_) :: mwm, mwn, mwnrhs, nz, temp
      integer(mwp_) :: ml_ptr, ml_row, ml_val
      integer :: m, n, nrhs, i
      integer, dimension(:), allocatable :: ptr, row, invp
      real(wp), dimension(:), allocatable :: val
      real(wp), dimension(:,:), allocatable :: rhs
      integer :: flag
      integer :: st
      integer :: num_threads, orig_num_threads, scalopt

      type(mc68_control) :: control68
      type(mc68_info) :: info68

      type(ma86_rcontrol) :: control
      type(ma86_rinfo) :: info

      character(len=200) :: errmsg

      real(wp) :: atime, ftime, otime
      integer :: t_start, t_stop, t_rate
      character(len=10) :: orderused

      ! Check number of arguments
      if(nrhs_in.lt.1 .or. nrhs_in .gt.3) &
         call MATLAB_error("Wrong number of input arguments")
      if(nlhs_in.ne.0 .and. nlhs_in.ne.1) &
         call MATLAB_error("Wrong number of output arguments")

      ! Check matrix is sparse
      call matlab_to_fortran(prhs(A_in), mwm, mwn, ptr, row, val, 'A')
      if(mwm.ne.mwn) call MATLAB_error("The matrix must be square")
      n = mwn

      ! Convert to HSL standard form
      call mc69_cscl_clean_real(HSL_MATRIX_REAL_SYM_INDEF, n, n, ptr, row, &
         flag, val=val)
      select case(flag)
      case(0,1,4,5)
         ! Sucess (expect oor entries as we are throwing away upr triangle!)
      case(-12)
         ! Hermitian but non-zero imaginary component on diagonal
         call MATLAB_context_error("hsl_ma86:ExpectsHermitian", &
            "HSL_MA86 Error: Hermitian factorization requested, but matrix has &
            &non-zero imaginary component on diagonal.")
      case default
         call unknownError("mc69_cscl_clean", flag)
      end select

      ! Read control in
      num_threads = 0
      if(nrhs_in.ge.control_in) then
         ! User supplied control structure
         if(MATLAB_is_structure(prhs(control_in))) then
            call copy_control_in(prhs(control_in), control, num_threads, &
               scalopt)
         endif
      endif

      ! Set number of threads
      orig_num_threads = 0
      if(num_threads.ne.0) then
         ! User specified a specific number of threads
         if(num_threads.le.0) &
            call MATLAB_error("control.num_threads must be positive")
!$       orig_num_threads = omp_get_max_threads()
         if(orig_num_threads.eq.0 .and. num_threads.gt.1) then
            call MATLAB_warning( &
               "Parallel run requested, but compiled without OpenMP. Running in&
               & serial")
         endif
!$       if(num_threads .gt. orig_num_threads) &
!$          call MATLAB_warning( &
!$             "control.num_threads exceeds omp_get_max_threads(). Proceeding &
!$             &with control.num_threads anyway.")
!$       call omp_set_num_threads(num_threads)
      endif

      ! Obtain ordering
      call system_clock(t_start)
      if(nrhs_in.ge.P_in) then
         ! User supplied (note: input as double, inverse of order)
         call matlab_to_fortran(prhs(P_in), invp, temp, 'P')
         if(mwn .ne. temp) then
            write(errmsg, "(4(a,i5))") &
               "Dimensions of A and P inconsistent: A=", &
               n, "x", n, ", size(P)=", temp
            call MATLAB_error(errmsg)
         endif
         allocate(order(n), stat=st)
         if(st.ne.0) goto 100
         do i = 1, n
            order(invp(i)) = i
         end do
         orderused = 'user'
      else
         ! call mc68
         allocate(order(n), stat=st)
         if(st.ne.0) goto 100
         ! disable printing
         control68%lp = 0
         control68%wp = 0
         control68%mp = 0
         ! 3 = MeTiS ordering
         call mc68_order(3, n, ptr, row, order, control68, info68)
         select case(info68%flag)
         case(0)
            ! OK. do nothing.
            orderused = 'MeTiS'
         case(-5)
            ! No MeTiS available.
            ! 1 = AMD ordering
            call mc68_order(1, n, ptr, row, order, control68, info68)
            if(info68%flag.ne.0) then
!$             if(orig_num_threads.ne.0) &
!$                call omp_set_num_threads(orig_num_threads)
               call unknownError("mc68_order", info68%flag)
            endif
            orderused = 'AMD'
         case default
!$          if(orig_num_threads .ne. 0) &
!$             call omp_set_num_threads(orig_num_threads)
            call unknownError("mc68_order", info68%flag)
         end select
      endif
      call system_clock(t_stop, t_rate)
      otime = real(t_stop-t_start)/t_rate

      ! Call analyse
      call system_clock(t_start)
      call ma86_analyse(n, ptr, row, order, keep, control, info)
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0)
         ! Success. Do nothing
      case(-2)
         ! Problem with order
         if(allocated(invp)) then
            write(errmsg, "(a, 10i4)") "Problem with P. P = ", &
               invp(1:min(n,10))
         else
            write(errmsg, "(a, 10i4)") "Problem with order. order = ", &
               order(1:min(n,10))
         endif
         call MATLAB_error(errmsg)
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma86_analyse", info%flag)
      end select
      atime = real(t_stop-t_start)/t_rate

      ! Obtain scaling
      deallocate(scaling, stat=st)
      allocate(scaling(n), stat=st)
      if(st.ne.0) goto 100
      call find_scaling(scalopt, n, ptr, row, val, scaling)

      ! Call factor
      call system_clock(t_start)
      if(allocated(scaling)) then
         call ma86_factor(n, ptr, row, val, order, keep, control, info, &
            scale=scaling)
      else
         call ma86_factor(n, ptr, row, val, order, keep, control, info)
      endif
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0:1)
         ! Success. Do nothing.
      case(-5)
         ! IEEE infinities in factorization
         call MATLAB_context_error("hsl_ma86:Infs", &
            "HSL_MA86 Error: IEEE infinities encountered during factorization. &
            &Try larger control.u or control.small.")
      case(-7)
         ! Bad value of static.
         call MATLAB_context_error("hsl_ma86:BadStatic", &
            "HSL_MA86 Error: control.static < abs(control.small) and &
            &control.static ~= 0.0.")
      case(2:3)
         ! Singular matrix
         call MATLAB_warning("HSL_MA86: Matrix was found to be singular.")
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma86_analyse", info%flag)
      end select
      ftime = real(t_stop-t_start)/t_rate

      ! Setup info structure if required
      if(nlhs_in.ge.info_out) then
         call copy_info_out(plhs(info_out), info, otime, atime, ftime, &
            orderused)
      endif

!$    if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)

      return
      100 continue
      call MATLAB_error("Insufficient memory")

   contains
      subroutine copy_info_out(ml_info, info, otime, atime, ftime, order)
         integer(mwp_), intent(inout) :: ml_info
         type(ma86_rinfo), intent(in) :: info
         real(wp), intent(in) :: otime, atime, ftime
         character(len=*) :: order

         integer(mwp_) :: comp
         integer*4 :: szf

         character(len=50), dimension(11) :: fields = (/ &
            "matrix_rank   ", &
            "num_delay     ", &
            "num_factor    ", &
            "num_flops     ", &
            "num_neg       ", &
            "num_perturbed ", &
            "num_two       ", &
            "order         ", &
            "order_time    ", &
            "analyse_time  ", &
            "factor_time   " /)

         ml_info = MATLAB_create_structure(fields)

         call MATLAB_set_field(ml_info, 'matrix_rank', info%matrix_rank)
         call MATLAB_set_field(ml_info, 'num_delay', info%num_delay)
         call MATLAB_set_field(ml_info, 'num_factor', info%num_factor)
         call MATLAB_set_field(ml_info, 'num_flops', info%num_flops)
         call MATLAB_set_field(ml_info, 'num_neg', info%num_neg)
         call MATLAB_set_field(ml_info, 'num_perturbed', info%num_perturbed)
         call MATLAB_set_field(ml_info, 'num_two', info%num_two)
         call MATLAB_set_field(ml_info, 'order', order)
         call MATLAB_set_field(ml_info, 'order_time', otime)
         call MATLAB_set_field(ml_info, 'analyse_time', atime)
         call MATLAB_set_field(ml_info, 'factor_time', ftime)
      end subroutine copy_info_out

   end subroutine ma86_matlab_analyse_factor_real

   subroutine find_scaling_real(scalopt, n, ptr, row, val, scaling)
      integer, intent(in) :: scalopt
      integer, intent(in) :: n
      integer, dimension(n+1), intent(in) :: ptr
      integer, dimension(ptr(n+1)-1), intent(in) :: row
      real(wp), dimension(ptr(n+1)-1), intent(in) :: val
      real(wp), allocatable, dimension(:), intent(out) :: scaling

      integer :: icntl(10), info(10)
      real(wp) :: cntl(10), rinfo(10)
      integer, dimension(:), allocatable :: iw
      real(wp), dimension(:), allocatable :: dw

      integer :: liw, ldw, st
      character(len=80) :: errstr

      deallocate(scaling, stat=st)
      if(scalopt.eq.0) return ! No scaling
      ! Otherwise assume mc77 1-norm scaling

      ! Allocate space
      allocate(scaling(n), stat=st)
      if(st.ne.0) goto 100

      liw = n
      ldw = ptr(n+1)-1 + 2*n ! note: mc77 user doc lies when it says m*(m+1)/2
      allocate(iw(liw), dw(ldw), stat=st)
      if(st.ne.0) goto 100

      call mc77id(icntl, cntl)
      icntl(4) = -1 ! Disable checking
      icntl(6) = -1 ! Symmetric matrix

      call mc77ad(1, n, n, ptr(n+1)-1, ptr, row, val, iw, liw, dw, &
         ldw, icntl, cntl, info, rinfo)
      if(info(1).lt.0) then
         write(errstr, "(a,i4)") "MC77 returned error", info(1)
         call MATLAB_error(errstr)
      endif
      scaling(1:n) = dw(1:n)

      return
      100 continue
      call MATLAB_error("Insufficient memory")
   end subroutine find_scaling_real

   ! [handle, info] = ma86_expert('factor', A[, P])
   ! handle and 'factor' already dealt with
   subroutine ma86_matlab_analyse_factor_complex(nlhs_in, plhs, nrhs_in, prhs, &
         keep, order, scaling)
      integer*4 :: nlhs_in, nrhs_in
      integer(mwp_) :: plhs(*), prhs(*)
      type(ma86_ckeep), intent(inout) :: keep
      integer, dimension(:), allocatable, intent(inout) :: order
      real(wp), dimension(:), allocatable, intent(inout) :: scaling

      integer, parameter :: A_in = 1, &
                            control_in = 2, &
                            P_in = 3
      integer, parameter :: info_out = 1

      integer(mws_) :: mwm, mwn, mwnrhs, nz
      integer(mwp_) :: ml_ptr, ml_row, ml_val
      integer :: m, n, nrhs, i
      integer, dimension(:), allocatable :: ptr, row, invp
      complex(wp), dimension(:), allocatable :: val
      integer :: flag
      integer :: st
      integer(mws_) :: temp
      integer :: num_threads, orig_num_threads, scalopt
      character(len=10) :: orderused

      type(mc68_control) :: control68
      type(mc68_info) :: info68

      type(ma86_ccontrol) :: control
      type(ma86_cinfo) :: info

      character(len=200) :: errmsg

      real(wp) :: atime, ftime, otime
      integer :: t_start, t_stop, t_rate

      integer :: matrix_type

      ! Check number of arguments
      if(nrhs_in.lt.1 .or. nrhs_in .gt.3) &
         call MATLAB_error("Wrong number of input arguments")
      if(nlhs_in.ne.0 .and. nlhs_in.ne.1) &
         call MATLAB_error("Wrong number of output arguments")
         
      ! Copy matrix in
      call matlab_to_fortran(prhs(A_in), mwm, mwn, ptr, row, val, 'A')
      if(mwm.ne.mwn) call MATLAB_error("The matrix must be square")
      n = mwn

      ! Copy control in
      num_threads = 0
      matrix_type = HSL_MATRIX_CPLX_SYM
      if(nrhs_in.ge.control_in) then
         ! User supplied control structure
         if(MATLAB_is_structure(prhs(control_in))) then
            call copy_control_in(prhs(control_in), control, num_threads, &
               matrix_type, scalopt)
         endif
      endif

      ! Convert to HSL standard form
      call mc69_cscl_clean_complex(matrix_type, n, n, ptr, row, flag, val=val)
      select case(flag)
      case(0,1,4,5)
         ! Sucess (expect oor entries as we are throwing away upr triangle!)
      case(-12)
         ! Hermitian but non-zero imaginary component on diagonal
         call MATLAB_context_error("hsl_ma86:ExpectsHermitian", &
            "HSL_MA86 Error: Hermitian factorization requested, but matrix has &
            &non-zero imaginary component on diagonal.")
      case default
         call unknownError("mc69_cscl_clean", flag)
      end select

      ! Set number of threads
      orig_num_threads = 0
      if(num_threads.ne.0) then
         ! User specified a specific number of threads
         if(num_threads.le.0) &
            call MATLAB_error("control.num_threads must be positive")
!$       orig_num_threads = omp_get_max_threads()
         if(orig_num_threads.eq.0 .and. num_threads.gt.1) then
            call MATLAB_warning( &
               "Parallel run requested, but compiled without OpenMP. Running in&
               & serial")
         endif
!$       if(num_threads .gt. orig_num_threads) &
!$          call MATLAB_warning( &
!$             "control.num_threads exceeds omp_get_max_threads(). Proceeding &
!$             &with control.num_threads anyway.")
!$       call omp_set_num_threads(num_threads)
      endif

      ! Obtain ordering
      call system_clock(t_start)
      if(nrhs_in.ge.P_in) then
         ! User supplied (note: input as double, inverse of order)
         call matlab_to_fortran(prhs(P_in), invp, temp, 'P')
         if(mwn .ne. temp) then
            write(errmsg, "(4(a,i5))") &
               "Dimensions of A and P inconsistent: A=", &
               n, "x", n, ", size(P)=", temp
            call MATLAB_error(errmsg)
         endif
         allocate(order(n), stat=st)
         if(st.ne.0) goto 100
         do i = 1, n
            order(invp(i)) = i
         end do
         orderused = 'user'
      else
         ! call mc68
         allocate(order(n), stat=st)
         if(st.ne.0) goto 100
         ! disable printing
         control68%lp = 0
         control68%wp = 0
         control68%mp = 0
         ! 3 = MeTiS ordering
         call mc68_order(3, n, ptr, row, order, control68, info68)
         select case(info68%flag)
         case(0)
            ! OK. do nothing.
            orderused = 'MeTiS'
         case(-5)
            ! No MeTiS available.
            ! 1 = AMD ordering
            call mc68_order(1, n, ptr, row, order, control68, info68)
            if(info68%flag.ne.0) then
!$             if(orig_num_threads.ne.0) &
!$                call omp_set_num_threads(orig_num_threads)
               call unknownError("mc68_order", info68%flag)
            endif
            orderused = 'AMD'
         case default
!$          if(orig_num_threads .ne. 0) &
!$             call omp_set_num_threads(orig_num_threads)
            call unknownError("mc68_order", info68%flag)
         end select
      endif
      call system_clock(t_stop, t_rate)
      otime = real(t_stop-t_start)/t_rate

      ! Call analyse
      call system_clock(t_start)
      call ma86_analyse(n, ptr, row, order, keep, control, info)
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0)
         ! Success. Do nothing
      case(-2)
         ! Problem with order
         if(allocated(invp)) then
            write(errmsg, "(a, 10i4)") "Problem with P. P = ", &
               invp(1:min(n,10))
         else
            write(errmsg, "(a, 10i4)") "Problem with order. order = ", &
               order(1:min(n,10))
         endif
         call MATLAB_error(errmsg)
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma86_analyse", info%flag)
      end select
      atime = real(t_stop-t_start)/t_rate

      ! Obtain scaling
      deallocate(scaling, stat=st)
      allocate(scaling(n), stat=st)
      if(st.ne.0) goto 100
      call find_scaling(scalopt, n, ptr, row, val, scaling)

      ! Call factor
      call system_clock(t_start)
      if(allocated(scaling)) then
         call ma86_factor(matrix_type, n, ptr, row, val, order, keep, control, &
            info, scale=scaling)
      else
         call ma86_factor(matrix_type, n, ptr, row, val, order, keep, control, &
            info)
      endif
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0:1)
         ! Success. Do nothing.
      case(-5)
         ! IEEE infinities in factorization
         call MATLAB_context_error("hsl_ma86:Infs", &
            "HSL_MA86 Error: IEEE infinities encountered during factorization. &
            &Try larger control.u or control.small.")
      case(-7)
         ! Bad value of static.
         call MATLAB_context_error("hsl_ma86:BadStatic", &
            "HSL_MA86 Error: control.static < abs(control.small) and &
            &control.static ~= 0.0.")
      case(2:3)
         ! Singular matrix
         call MATLAB_warning("HSL_MA86: Matrix was found to be singular.")
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma86_analyse", info%flag)
      end select
      ftime = real(t_stop-t_start)/t_rate

      ! Setup info structure if required
      if(nlhs_in.ge.info_out) then
         call copy_info_out(plhs(info_out), info, otime, atime, ftime, &
            orderused)
      endif

!$    if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)

      return
      100 continue
      call MATLAB_error("Insufficient memory")

   contains
      subroutine copy_info_out(ml_info, info, otime, atime, ftime, order)
         integer(mwp_), intent(inout) :: ml_info
         type(ma86_cinfo), intent(in) :: info
         real(wp), intent(in) :: otime, atime, ftime
         character(len=*) :: order

         integer(mwp_) :: comp
         integer*4 :: szf

         character(len=50), dimension(11) :: fields = (/ &
            "matrix_rank   ", &
            "num_delay     ", &
            "num_factor    ", &
            "num_flops     ", &
            "num_neg       ", &
            "num_perturbed ", &
            "num_two       ", &
            "order         ", &
            "order_time    ", &
            "analyse_time  ", &
            "factor_time   " /)

         ml_info = MATLAB_create_structure(fields)

         call MATLAB_set_field(ml_info, 'matrix_rank', info%matrix_rank)
         call MATLAB_set_field(ml_info, 'num_delay', info%num_delay)
         call MATLAB_set_field(ml_info, 'num_factor', info%num_factor)
         call MATLAB_set_field(ml_info, 'num_flops', info%num_flops)
         call MATLAB_set_field(ml_info, 'num_neg', info%num_neg)
         call MATLAB_set_field(ml_info, 'num_perturbed', info%num_perturbed)
         call MATLAB_set_field(ml_info, 'num_two', info%num_two)
         call MATLAB_set_field(ml_info, 'order', order)
         call MATLAB_set_field(ml_info, 'order_time', otime)
         call MATLAB_set_field(ml_info, 'analyse_time', atime)
         call MATLAB_set_field(ml_info, 'factor_time', ftime)
      end subroutine copy_info_out

   end subroutine ma86_matlab_analyse_factor_complex

   subroutine find_scaling_complex(scalopt, n, ptr, row, val, scaling)
      integer, intent(in) :: scalopt
      integer, intent(in) :: n
      integer, dimension(n+1), intent(in) :: ptr
      integer, dimension(ptr(n+1)-1), intent(in) :: row
      complex(wp), dimension(ptr(n+1)-1), intent(in) :: val
      real(wp), dimension(:), allocatable, intent(out) :: scaling

      integer :: icntl(10), info(10)
      real(wp) :: cntl(10), rinfo(10)
      integer, dimension(:), allocatable :: iw
      real(wp), dimension(:), allocatable :: dw

      integer :: liw, ldw, st
      character(len=80) :: errstr

      ! Deallocate scaling to be safe then check if we need to scale at all
      deallocate(scaling, stat=st)
      if(scalopt.eq.0) return ! No scaling
      ! Otherwise assume mc77 1-norm scaling

      ! Allocate space
      allocate(scaling(n), stat=st)
      if(st.ne.0) goto 100

      liw = n
      ldw = ptr(n+1)-1 + 2*n ! note: mc77 user doc lies when it says m*(m+1)/2
      allocate(iw(liw), dw(ldw), stat=st)
      if(st.ne.0) goto 100

      call mc77iz(icntl, cntl)
      icntl(4) = -1 ! Disable checking
      icntl(6) = -1 ! Symmetric matrix

      call mc77az(1, n, n, ptr(n+1)-1, ptr, row, val, iw, liw, dw, &
         ldw, icntl, cntl, info, rinfo)
      if(info(1).lt.0) then
         write(errstr, "(a,i4)") "MC77 returned error", info(1)
         call MATLAB_error(errstr)
      endif
      scaling(1:n) = dw(1:n)

      return
      100 continue
      call MATLAB_error("Insufficient memory")
   end subroutine find_scaling_complex

   ! [x, info] = ma86_expert('solve', handle, b)
   ! 'solve' and handle already dealt with
   subroutine ma86_matlab_solve_real(nlhs_in, plhs, nrhs_in, prhs, keep, order,&
         scaling)
      integer*4 :: nlhs_in, nrhs_in
      integer(mwp_) :: plhs(*), prhs(*)
      type(ma86_rkeep), intent(inout) :: keep
      integer, dimension(ma86_get_n__(keep)), intent(in) :: order
      real(wp), dimension(:), allocatable, intent(inout) :: scaling

      integer, parameter :: b_in = 1, &
                            control_in = 2
      integer, parameter :: x_out = 1, &
                            info_out = 2

      integer(mws_) :: mwn, mwnrhs
      integer(mwp_) :: ml_ptr, ml_row, ml_val
      real(wp), dimension(:,:), allocatable :: rhs
      integer :: flag
      integer :: st
      integer :: temp

      type(ma86_rcontrol) :: control
      type(ma86_rinfo) :: info

      character(len=200) :: errmsg
      integer :: n, nrhs
      integer :: num_threads, orig_num_threads, scalopt

      real(wp) :: stime
      integer :: t_start, t_stop, t_rate

      ! Check number of arguments
      if(nrhs_in.lt.1) call MATLAB_error("Insufficient input arguments")
      if(nrhs_in.gt.2) call MATLAB_error("Too many input arguments")
      if(nlhs_in.lt.1) call MATLAB_error("Insufficient output arguments")
      if(nlhs_in.gt.2) call MATLAB_error("Too many output arguments")

      n = ma86_get_n__(keep)

      ! Get rhs
      call matlab_to_fortran(prhs(b_in), rhs, mwn, mwnrhs, 'b')
      if(mwn .ne. n) then
         write(errmsg, "(4(a,i5))") "Dimensions of A and b inconsistent: A=", &
            n, "x", n, ", b=", mwn, "x", mwnrhs
         call MATLAB_error(errmsg)
      endif
      nrhs = mwnrhs

      ! Copy control in
      num_threads = 0
      if(nrhs_in.ge.control_in) then
         ! User supplied control structure
         if(MATLAB_is_structure(prhs(control_in))) then
            call copy_control_in(prhs(control_in), control, num_threads, &
               scalopt)
         endif
      endif

      ! Set number of threads
      orig_num_threads = 0
      if(num_threads.ne.0) then
         ! User specified a specific number of threads
         if(num_threads.le.0) &
            call MATLAB_error("control.num_threads must be positive")
!$          orig_num_threads = omp_get_max_threads()
         if(orig_num_threads.eq.0 .and. num_threads.gt.1) then
            call MATLAB_warning( &
               "Parallel run requested, but compiled without OpenMP. Running in&
               & serial")
         endif
!$       if(num_threads .gt. orig_num_threads) &
!$          call MATLAB_warning( &
!$             "control.num_threads exceeds omp_get_max_threads(). Proceeding &
!$             &with control.num_threads anyway.")
!$       call omp_set_num_threads(num_threads)
      endif
      
      ! Call solve
      call system_clock(t_start)
      if(allocated(scaling)) then
         call ma86_solve(nrhs, n, rhs, order, keep, control, info, &
            scale=scaling)
      else
         call ma86_solve(nrhs, n, rhs, order, keep, control, info)
      endif
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0)
         ! Success. Do nothing.
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma86_analyse", info%flag)
      end select
      stime = real(t_stop-t_start)/t_rate

      ! Setup info structure if required
      if(nlhs_in.ge.info_out) then
         call copy_info_out(plhs(info_out), info, stime)
      endif

      ! Copy solution out
      plhs(x_out) = fortran_to_matlab(rhs)

!$    if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)

      return
      100 continue
      call MATLAB_error("Insufficient memory")

   contains
      subroutine copy_info_out(ml_info, info, stime)
         integer(mwp_), intent(inout) :: ml_info
         type(ma86_rinfo), intent(in) :: info
         real(wp), intent(in) :: stime

         integer(mwp_) :: comp
         integer*4 :: szf

         character(len=50), dimension(1) :: fields = (/ &
            "solve_time    " /)

         ml_info = MATLAB_create_structure(fields)

         call MATLAB_set_field(ml_info, 'solve_time', stime)
      end subroutine copy_info_out

   end subroutine ma86_matlab_solve_real

   ! [x, info] = ma86_expert('solve', handle, b)
   ! 'solve' and handle already dealt with
   subroutine ma86_matlab_solve_complex(nlhs_in, plhs, nrhs_in, prhs, keep, &
         order, scaling)
      integer*4 :: nlhs_in, nrhs_in
      integer(mwp_) :: plhs(*), prhs(*)
      type(ma86_ckeep), intent(inout) :: keep
      integer, dimension(ma86_get_n__(keep)), intent(in) :: order
      real(wp), dimension(:), allocatable, intent(inout) :: scaling

      integer, parameter :: b_in = 1, &
                            control_in = 2
      integer, parameter :: x_out = 1, &
                            info_out = 2

      integer(mws_) :: mwn, mwnrhs
      integer(mwp_) :: ml_ptr, ml_row, ml_val
      complex(wp), dimension(:,:), allocatable :: rhs
      integer :: flag
      integer :: st
      integer :: temp

      type(ma86_ccontrol) :: control
      type(ma86_cinfo) :: info

      character(len=200) :: errmsg
      integer :: n, nrhs
      integer :: num_threads, orig_num_threads, scalopt

      real(wp) :: stime
      integer :: t_start, t_stop, t_rate
      integer :: matrix_type

      ! Check number of arguments
      if(nrhs_in.lt.1) call MATLAB_error("Insufficient input arguments")
      if(nrhs_in.gt.2) call MATLAB_error("Too many input arguments")
      if(nlhs_in.lt.1) call MATLAB_error("Insufficient output arguments")
      if(nlhs_in.gt.2) call MATLAB_error("Too many output arguments")

      n = ma86_get_n__(keep)

      ! Obtain rhs vector
      call matlab_to_fortran(prhs(b_in), rhs, mwn, mwnrhs, 'b')
      if(mwn .ne. n) then
         write(errmsg, "(4(a,i5))") "Dimensions of A and b inconsistent: A=", &
           n, "x", n, ", b=", mwn, "x", mwnrhs
         call MATLAB_error(errmsg)
      endif
      nrhs = mwnrhs

      ! Copy control in
      num_threads = 0
      matrix_type = HSL_MATRIX_CPLX_SYM
      if(nrhs_in.ge.control_in) then
         ! User supplied control structure
         if(MATLAB_is_structure(prhs(control_in))) then
            call copy_control_in(prhs(control_in), control, num_threads, &
               matrix_type, scalopt)
         endif
      endif

      ! Set number of threads
      orig_num_threads = 0
      if(num_threads.ne.0) then
         ! User specified a specific number of threads
         if(num_threads.le.0) &
            call MATLAB_error("control.num_threads must be positive")
!$          orig_num_threads = omp_get_max_threads()
         if(orig_num_threads.eq.0 .and. num_threads.gt.1) then
            call MATLAB_warning( &
               "Parallel run requested, but compiled without OpenMP. Running in&
               & serial")
         endif
!$       if(num_threads .gt. orig_num_threads) &
!$          call MATLAB_warning( &
!$             "control.num_threads exceeds omp_get_max_threads(). Proceeding &
!$             &with control.num_threads anyway.")
!$       call omp_set_num_threads(num_threads)
      endif

      ! Call solve
      call system_clock(t_start)
      if(allocated(scaling)) then
         call ma86_solve(nrhs, n, rhs, order, keep, control, info, &
            scale=scaling)
      else
         call ma86_solve(nrhs, n, rhs, order, keep, control, info)
      endif
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0)
         ! Success. Do nothing.
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma86_analyse", info%flag)
      end select
      stime = real(t_stop-t_start)/t_rate

      ! Setup info structure if required
      if(nlhs_in.ge.info_out) then
         call copy_info_out(plhs(info_out), info, stime)
      endif

      ! Copy solution out
      plhs(x_out) = fortran_to_matlab(rhs)

      ! Restore number of threads
!$    if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)

      return
      100 continue
      call MATLAB_error("Insufficient memory")
   contains
      subroutine copy_info_out(ml_info, info, stime)
         integer(mwp_), intent(inout) :: ml_info
         type(ma86_cinfo), intent(in) :: info
         real(wp), intent(in) :: stime

         integer(mwp_) :: comp
         integer*4 :: szf

         character(len=50), dimension(1) :: fields = (/ &
            "solve_time    " /)

         ml_info = MATLAB_create_structure(fields)

         call MATLAB_set_field(ml_info, 'solve_time', stime)
      end subroutine copy_info_out


   end subroutine ma86_matlab_solve_complex

   ! [x, info, handle] = ma86_expert('backslash', A, b[, control, P])
   ! handle and 'backslash' already dealt with
   subroutine ma86_matlab_backslash_real(nlhs_in, plhs, nrhs_in, prhs, &
         keep, order, scaling)
      integer*4 :: nlhs_in, nrhs_in
      integer(mwp_) :: plhs(*), prhs(*)
      type(ma86_rkeep), intent(inout) :: keep
      integer, dimension(:), allocatable, intent(inout) :: order
      real(wp), dimension(:), allocatable, intent(inout) :: scaling

      integer, parameter :: A_in = 1, &
                            b_in = 2, &
                            control_in = 3, &
                            P_in = 4
      integer, parameter :: x_out = 1, &
                            info_out = 2

      integer(mws_) :: mwm, mwn, mwnrhs, nz, temp
      integer(mwp_) :: ml_ptr, ml_row, ml_val
      integer :: m, n, nrhs, i
      integer, dimension(:), allocatable :: ptr, row, invp
      real(wp), dimension(:), allocatable :: val
      real(wp), dimension(:,:), allocatable :: rhs
      integer :: flag
      integer :: st
      integer :: num_threads, orig_num_threads, scalopt

      type(mc68_control) :: control68
      type(mc68_info) :: info68

      type(ma86_rcontrol) :: control
      type(ma86_rinfo) :: info

      character(len=200) :: errmsg

      real(wp) :: atime, ftime, otime
      integer :: t_start, t_stop, t_rate
      character(len=10) :: orderused

      ! Check number of arguments
      if(nrhs_in.lt.2) call MATLAB_error("Insufficient input arguments")
      if(nrhs_in.gt.4) call MATLAB_error("Too many input arguments")
      if(nlhs_in.lt.1) call MATLAB_error("Insufficient output arguments")
      if(nlhs_in.gt.3) call MATLAB_error("Too many output arguments")

      ! Get sparse matrix
      call matlab_to_fortran(prhs(A_in), mwm, mwn, ptr, row, val, 'A')
      if(mwm.ne.mwn) call MATLAB_error("The matrix must be square")
      n = mwn

      ! Convert to HSL standard form
      call mc69_cscl_clean_real(HSL_MATRIX_REAL_SYM_INDEF, n, n, ptr, row, &
         flag, val=val)
      select case(flag)
      case(0,1,4,5)
         ! Sucess (expect oor entries as we are throwing away upr triangle!)
      case(-12)
         ! Hermitian but non-zero imaginary component on diagonal
         call MATLAB_context_error("hsl_ma86:ExpectsHermitian", &
            "HSL_MA86 Error: Hermitian factorization requested, but matrix has &
            &non-zero imaginary component on diagonal.")
      case default
         call unknownError("mc69_cscl_clean", flag)
      end select

      ! Obtain rhs vector
      call matlab_to_fortran(prhs(b_in), rhs, temp, mwnrhs, 'b')
      if(mwn .ne. temp) then
         write(errmsg, "(4(a,i5))") "Dimensions of A and b inconsistent: A=", &
            mwn, "x", mwn, ", b=", temp, "x", mwnrhs
         call MATLAB_error(errmsg)
      endif
      nrhs = mwnrhs

      ! Read control in
      num_threads = 0
      if(nrhs_in.ge.control_in) then
         ! User supplied control structure
         if(MATLAB_is_structure(prhs(control_in))) then
            call copy_control_in(prhs(control_in), control, num_threads, &
               scalopt)
         endif
      endif

      ! Set number of threads
      orig_num_threads = 0
      if(num_threads.ne.0) then
         ! User specified a specific number of threads
         if(num_threads.le.0) &
            call MATLAB_error("control.num_threads must be positive")
!$       orig_num_threads = omp_get_max_threads()
         if(orig_num_threads.eq.0 .and. num_threads.gt.1) then
            call MATLAB_warning( &
               "Parallel run requested, but compiled without OpenMP. Running in&
               & serial")
         endif
!$       if(num_threads .gt. orig_num_threads) &
!$          call MATLAB_warning( &
!$             "control.num_threads exceeds omp_get_max_threads(). Proceeding &
!$             &with control.num_threads anyway.")
!$       call omp_set_num_threads(num_threads)
      endif

      ! Obtain ordering
      call system_clock(t_start)
      if(nrhs_in.ge.P_in) then
         ! User supplied (note: input as double, inverse of order)
         call matlab_to_fortran(prhs(P_in), invp, temp, 'P')
         if(mwn .ne. temp) then
            write(errmsg, "(4(a,i5))") &
               "Dimensions of A and P inconsistent: A=", &
               n, "x", n, ", size(P)=", temp
            call MATLAB_error(errmsg)
         endif
         allocate(order(n), stat=st)
         if(st.ne.0) goto 100
         do i = 1, n
            order(invp(i)) = i
         end do
         orderused = "user"
      else
         ! call mc68
         allocate(order(n), stat=st)
         if(st.ne.0) goto 100
         ! disable printing
         control68%lp = 0
         control68%wp = 0
         control68%mp = 0
         ! 3 = MeTiS ordering
         call mc68_order(3, n, ptr, row, order, control68, info68)
         select case(info68%flag)
         case(0)
            ! OK. do nothing.
            orderused = 'MeTiS'
         case(-5)
            ! No MeTiS available.
            ! 1 = AMD ordering
            call mc68_order(1, n, ptr, row, order, control68, info68)
            if(info68%flag.ne.0) then
!$             if(orig_num_threads.ne.0) &
!$                call omp_set_num_threads(orig_num_threads)
               call unknownError("mc68_order", info68%flag)
            endif
            orderused = 'AMD'
         case default
!$          if(orig_num_threads .ne. 0) &
!$             call omp_set_num_threads(orig_num_threads)
            call unknownError("mc68_order", info68%flag)
         end select
      endif
      call system_clock(t_stop, t_rate)
      otime = real(t_stop-t_start)/t_rate

      ! Call analyse
      call system_clock(t_start)
      call ma86_analyse(n, ptr, row, order, keep, control, info)
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0)
         ! Success. Do nothing
      case(-2)
         ! Problem with order
         if(allocated(invp)) then
            write(errmsg, "(a, 10i4)") "Problem with P. P = ", &
               invp(1:min(n,10))
         else
            write(errmsg, "(a, 10i4)") "Problem with order. order = ", &
               order(1:min(n,10))
         endif
         call MATLAB_error(errmsg)
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma86_analyse", info%flag)
      end select
      atime = real(t_stop-t_start)/t_rate

      ! Obtain scaling
      call find_scaling(scalopt, n, ptr, row, val, scaling)

      ! Call factor_solve
      call system_clock(t_start)
      if(allocated(scaling)) then
         call ma86_factor_solve(n, ptr, row, val, order, keep, control, info, &
            nrhs, n, rhs, scale=scaling)
      else
         call ma86_factor_solve(n, ptr, row, val, order, keep, control, info, &
            nrhs, n, rhs)
      endif
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0:1)
         ! Success. Do nothing.
      case(-5)
         ! IEEE infinities in factorization
         call MATLAB_context_error("hsl_ma86:Infs", &
            "HSL_MA86 Error: IEEE infinities encountered during factorization. &
            &Try larger control.u or control.small.")
      case(-7)
         ! Bad value of static.
         call MATLAB_context_error("hsl_ma86:BadStatic", &
            "HSL_MA86 Error: control.static < abs(control.small) and &
            &control.static ~= 0.0.")
      case(2:3)
         ! Singular matrix
         call MATLAB_warning("HSL_MA86: Matrix was found to be singular.")
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma86_analyse", info%flag)
      end select
      ftime = real(t_stop-t_start)/t_rate

      ! Setup info structure if required
      if(nlhs_in.ge.info_out) then
         call copy_info_out(plhs(info_out), info, otime, atime, ftime, &
            orderused)
      endif

      ! Copy solution out
      plhs(x_out) = fortran_to_matlab(rhs)

      ! Restore number of threads
!$    if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)

      return
      100 continue
      call MATLAB_error("Insufficient memory")

   contains
      subroutine copy_info_out(ml_info, info, otime, atime, ftime, order)
         integer(mwp_), intent(inout) :: ml_info
         type(ma86_rinfo), intent(in) :: info
         real(wp), intent(in) :: otime, atime, ftime
         character(len=*) :: order

         integer(mwp_) :: comp
         integer*4 :: szf

         character(len=50), dimension(11) :: fields = (/ &
            "matrix_rank      ", &
            "num_delay        ", &
            "num_factor       ", &
            "num_flops        ", &
            "num_neg          ", &
            "num_perturbed    ", &
            "num_two          ", &
            "order            ", &
            "order_time       ", &
            "analyse_time     ", &
            "factor_solve_time" /)

         ml_info = MATLAB_create_structure(fields)

         call MATLAB_set_field(ml_info, 'matrix_rank', info%matrix_rank)
         call MATLAB_set_field(ml_info, 'num_delay', info%num_delay)
         call MATLAB_set_field(ml_info, 'num_factor', info%num_factor)
         call MATLAB_set_field(ml_info, 'num_flops', info%num_flops)
         call MATLAB_set_field(ml_info, 'num_neg', info%num_neg)
         call MATLAB_set_field(ml_info, 'num_perturbed', info%num_perturbed)
         call MATLAB_set_field(ml_info, 'num_two', info%num_two)
         call MATLAB_set_field(ml_info, 'order', order)
         call MATLAB_set_field(ml_info, 'order_time', otime)
         call MATLAB_set_field(ml_info, 'analyse_time', atime)
         call MATLAB_set_field(ml_info, 'factor_solve_time', ftime)
      end subroutine copy_info_out

   end subroutine ma86_matlab_backslash_real

   ! [x, info, handle] = ma86_expert('backslash', A, b[, control, P])
   ! handle and 'backslash' already dealt with
   subroutine ma86_matlab_backslash_complex(nlhs_in, plhs, nrhs_in, prhs, &
         keep, order, scaling)
      integer*4 :: nlhs_in, nrhs_in
      integer(mwp_) :: plhs(*), prhs(*)
      type(ma86_ckeep), intent(inout) :: keep
      integer, dimension(:), allocatable, intent(inout) :: order
      real(wp), dimension(:), allocatable, intent(inout) :: scaling

      integer, parameter :: A_in = 1, &
                            b_in = 2, &
                            control_in = 3, &
                            P_in = 4
      integer, parameter :: x_out = 1, &
                            info_out = 2

      integer(mws_) :: mwm, mwn, mwnrhs, nz, temp
      integer(mwp_) :: ml_ptr, ml_row, ml_val
      integer :: m, n, nrhs, i
      integer, dimension(:), allocatable :: ptr, row, invp
      complex(wp), dimension(:), allocatable :: val
      complex(wp), dimension(:,:), allocatable :: rhs
      integer :: flag
      integer :: st
      integer :: num_threads, orig_num_threads, scalopt

      type(mc68_control) :: control68
      type(mc68_info) :: info68

      type(ma86_ccontrol) :: control
      type(ma86_cinfo) :: info

      character(len=200) :: errmsg

      real(wp) :: atime, ftime, otime
      integer :: t_start, t_stop, t_rate
      character(len=10) :: orderused
      integer :: matrix_type

      ! Check number of arguments
      if(nrhs_in.lt.2) call MATLAB_error("Insufficient input arguments")
      if(nrhs_in.gt.4) call MATLAB_error("Too many input arguments")
      if(nlhs_in.lt.1) call MATLAB_error("Insufficient output arguments")
      if(nlhs_in.gt.3) call MATLAB_error("Too many output arguments")

      ! Get sparse matrix
      call matlab_to_fortran(prhs(A_in), mwm, mwn, ptr, row, val, 'A')
      if(mwm.ne.mwn) call MATLAB_error("The matrix must be square")
      n = mwn

      ! Read control in
      num_threads = 0
      matrix_type = HSL_MATRIX_CPLX_SYM
      if(nrhs_in.ge.control_in) then
         ! User supplied control structure
         if(MATLAB_is_structure(prhs(control_in))) then
            call copy_control_in(prhs(control_in), control, num_threads, &
               matrix_type, scalopt)
         endif
      endif

      ! Convert to HSL standard form
      call mc69_cscl_clean_complex(matrix_type, n, n, ptr, row, &
         flag, val=val)
      select case(flag)
      case(0,1,4,5)
         ! Sucess (expect oor entries as we are throwing away upr triangle!)
      case(-12)
         ! Hermitian but non-zero imaginary component on diagonal
         call MATLAB_context_error("hsl_ma86:ExpectsHermitian", &
            "HSL_MA86 Error: Hermitian factorization requested, but matrix has &
            &non-zero imaginary component on diagonal.")
      case default
         call unknownError("mc69_cscl_clean", flag)
      end select

      ! Obtain rhs vector
      call matlab_to_fortran(prhs(b_in), rhs, temp, mwnrhs, 'b')
      if(mwn .ne. temp) then
         write(errmsg, "(4(a,i5))") "Dimensions of A and b inconsistent: A=", &
            mwn, "x", mwn, ", b=", temp, "x", mwnrhs
         call MATLAB_error(errmsg)
      endif
      nrhs = mwnrhs

      ! Set number of threads
      orig_num_threads = 0
      if(num_threads.ne.0) then
         ! User specified a specific number of threads
         if(num_threads.le.0) &
            call MATLAB_error("control.num_threads must be positive")
!$       orig_num_threads = omp_get_max_threads()
         if(orig_num_threads.eq.0 .and. num_threads.gt.1) then
            call MATLAB_warning( &
               "Parallel run requested, but compiled without OpenMP. Running in&
               & serial")
         endif
!$       if(num_threads .gt. orig_num_threads) &
!$          call MATLAB_warning( &
!$             "control.num_threads exceeds omp_get_max_threads(). Proceeding &
!$             &with control.num_threads anyway.")
!$       call omp_set_num_threads(num_threads)
      endif

      ! Obtain ordering
      call system_clock(t_start)
      if(nrhs_in.ge.P_in) then
         ! User supplied (note: input as double, inverse of order)
         call matlab_to_fortran(prhs(P_in), invp, temp, 'P')
         if(mwn .ne. temp) then
            write(errmsg, "(4(a,i5))") &
               "Dimensions of A and P inconsistent: A=", &
               n, "x", n, ", size(P)=", temp
            call MATLAB_error(errmsg)
         endif
         allocate(order(n), stat=st)
         if(st.ne.0) goto 100
         do i = 1, n
            order(invp(i)) = i
         end do
         orderused = 'user'
      else
         ! call mc68
         allocate(order(n), stat=st)
         if(st.ne.0) goto 100
         ! disable printing
         control68%lp = 0
         control68%wp = 0
         control68%mp = 0
         ! 3 = MeTiS ordering
         call mc68_order(3, n, ptr, row, order, control68, info68)
         select case(info68%flag)
         case(0)
            ! OK. do nothing.
            orderused = 'MeTiS'
         case(-5)
            ! No MeTiS available.
            ! 1 = AMD ordering
            call mc68_order(1, n, ptr, row, order, control68, info68)
            if(info68%flag.ne.0) then
!$             if(orig_num_threads.ne.0) &
!$                call omp_set_num_threads(orig_num_threads)
               call unknownError("mc68_order", info68%flag)
            endif
            orderused = 'AMD'
         case default
!$          if(orig_num_threads .ne. 0) &
!$             call omp_set_num_threads(orig_num_threads)
            call unknownError("mc68_order", info68%flag)
         end select
      endif
      call system_clock(t_stop, t_rate)
      otime = real(t_stop-t_start)/t_rate

      ! Call analyse
      call system_clock(t_start)
      call ma86_analyse(n, ptr, row, order, keep, control, info)
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0)
         ! Success. Do nothing
      case(-2)
         ! Problem with order
         if(allocated(invp)) then
            write(errmsg, "(a, 10i4)") "Problem with P. P = ", &
               invp(1:min(n,10))
         else
            write(errmsg, "(a, 10i4)") "Problem with order. order = ", &
               order(1:min(n,10))
         endif
         call MATLAB_error(errmsg)
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma86_analyse", info%flag)
      end select
      atime = real(t_stop-t_start)/t_rate

      ! Obtain scaling
      deallocate(scaling, stat=st)
      allocate(scaling(n), stat=st)
      if(st.ne.0) goto 100
      call find_scaling(scalopt, n, ptr, row, val, scaling)

      ! Call factor_solve
      call system_clock(t_start)
      if(allocated(scaling)) then
         call ma86_factor_solve(matrix_type, n, ptr, row, val, order, keep, &
            control, info, nrhs, n, rhs, scale=scaling)
      else
         call ma86_factor_solve(matrix_type, n, ptr, row, val, order, keep, &
            control, info, nrhs, n, rhs)
      endif
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0:1)
         ! Success. Do nothing.
      case(-5)
         ! IEEE infinities in factorization
         call MATLAB_context_error("hsl_ma86:Infs", &
            "HSL_MA86 Error: IEEE infinities encountered during factorization. &
            &Try larger control.u or control.small.")
      case(-7)
         ! Bad value of static.
         call MATLAB_context_error("hsl_ma86:BadStatic", &
            "HSL_MA86 Error: control.static < abs(control.small) and &
            &control.static ~= 0.0.")
      case(2:3)
         ! Singular matrix
         call MATLAB_warning("HSL_MA86: Matrix was found to be singular.")
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma86_analyse", info%flag)
      end select
      ftime = real(t_stop-t_start)/t_rate

      ! Setup info structure if required
      if(nlhs_in.ge.info_out) then
         call copy_info_out(plhs(info_out), info, otime, atime, ftime, &
            orderused)
      endif

      ! Copy solution out
      plhs(x_out) = fortran_to_matlab(rhs)

      ! Restore number of threads
!$    if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)

      return
      100 continue
      call MATLAB_error("Insufficient memory")

   contains
      subroutine copy_info_out(ml_info, info, otime, atime, ftime, order)
         integer(mwp_), intent(inout) :: ml_info
         type(ma86_cinfo), intent(in) :: info
         real(wp), intent(in) :: otime, atime, ftime
         character(len=*) :: order

         integer(mwp_) :: comp
         integer*4 :: szf

         character(len=50), dimension(11) :: fields = (/ &
            "matrix_rank      ", &
            "num_delay        ", &
            "num_factor       ", &
            "num_flops        ", &
            "num_neg          ", &
            "num_perturbed    ", &
            "num_two          ", &
            "order            ", &
            "order_time       ", &
            "analyse_time     ", &
            "factor_solve_time" /)

         ml_info = MATLAB_create_structure(fields)

         call MATLAB_set_field(ml_info, 'matrix_rank', info%matrix_rank)
         call MATLAB_set_field(ml_info, 'num_delay', info%num_delay)
         call MATLAB_set_field(ml_info, 'num_factor', info%num_factor)
         call MATLAB_set_field(ml_info, 'num_flops', info%num_flops)
         call MATLAB_set_field(ml_info, 'num_neg', info%num_neg)
         call MATLAB_set_field(ml_info, 'num_perturbed', info%num_perturbed)
         call MATLAB_set_field(ml_info, 'num_two', info%num_two)
         call MATLAB_set_field(ml_info, 'order', order)
         call MATLAB_set_field(ml_info, 'order_time', otime)
         call MATLAB_set_field(ml_info, 'analyse_time', atime)
         call MATLAB_set_field(ml_info, 'factor_solve_time', ftime)
      end subroutine copy_info_out

   end subroutine ma86_matlab_backslash_complex

   ! ma86_expert('destroy', handle)
   ! arguments 'destroy' and handle already removed
   subroutine ma86_matlab_finalise_real(keep)
      type(ma86_rkeep), intent(inout) :: keep

      type(ma86_rcontrol) :: control

      call ma86_finalise(keep, control)
   end subroutine ma86_matlab_finalise_real

   ! ma86_expert('destroy', handle)
   ! arguments 'destroy' and handle already removed
   subroutine ma86_matlab_finalise_complex(keep)
      type(ma86_ckeep), intent(inout) :: keep

      type(ma86_ccontrol) :: control

      call ma86_finalise(keep, control)
   end subroutine ma86_matlab_finalise_complex

end module ma86_matlab_main

! This module looks after a SAVEd set of variables mapping integer handles
! to Fortran keep and order variables
module ma86_handles
   use hsl_ma86_double, only: ma86_rkeep => ma86_keep
   use hsl_ma86_double_complex, only: ma86_ckeep => ma86_keep
   use ma86_matlab_main
   implicit none

   ! Data associated with the handle
   ! Considered to be empty if both associated(rkeep) and associated(ckeep)
   ! are .false.
   type ma86_hdl
      type(ma86_rkeep), pointer :: rkeep => null()
      type(ma86_ckeep), pointer :: ckeep => null()
      integer, dimension(:), allocatable :: order
      real(wp), dimension(:), allocatable :: scaling
   end type ma86_hdl

   ! How many handles initally and how much increase once exhausted
   integer, parameter :: initial_handles = 5
   double precision, parameter :: multiplier = 2.0

   ! SAVEd data
   integer, save :: next_handle = 1
   integer, save :: total_handles = 0
   type(ma86_hdl), dimension(:), allocatable, save :: handles

contains

   integer function ma86_new_handle(complexFlag)
      logical, intent(in) :: complexFlag

      type(ma86_hdl), dimension(:), allocatable :: temp
      integer :: i

      ! Do we need to expand the number of available handles?
      if (next_handle .gt. total_handles) then
         if(total_handles.ne.0) then
            ! Need to expand existing handle selection
            allocate(temp(total_handles))
            do i = 1, total_handles
               temp(i)%rkeep => handles(i)%rkeep
               temp(i)%ckeep => handles(i)%ckeep
               if(allocated(handles(i)%order)) then
                  allocate(temp(i)%order(size(handles(i)%order)))
                  temp(i)%order(:) = handles(i)%order(:)
               endif
            end do
            deallocate(handles)
            total_handles = max(int(multiplier*total_handles), total_handles)
            allocate(handles(total_handles))
            do i = 1, size(temp)
               handles(i)%rkeep => temp(i)%rkeep
               handles(i)%ckeep => temp(i)%ckeep
               if(allocated(temp(i)%order)) then
                  allocate(handles(i)%order(size(temp(i)%order)))
                  handles(i)%order(:) = temp(i)%order(:)
               endif
               if(allocated(temp(i)%scaling)) then
                  allocate(handles(i)%scaling(size(temp(i)%scaling)))
                  handles(i)%scaling(:) = temp(i)%scaling(:)
               endif
            end do
            deallocate(temp)
         else
            ! First call since module loaded
            total_handles = initial_handles
            allocate(handles(total_handles))

            ! Register clean function
            call mexAtExit(cleanup_all_handles)
         endif
      endif

      ma86_new_handle = next_handle
      if(complexFlag) then
         allocate(handles(next_handle)%ckeep)
      else
         allocate(handles(next_handle)%rkeep)
      endif
      next_handle = next_handle + 1
   end function ma86_new_handle

   ! This routine is called at unload of this module from MATLAB.
   ! It shuold cleanup all SAVEd data
   subroutine cleanup_all_handles()
      integer :: i, st

      do i = 1, next_handle-1
         call cleanup_handle(i)
      end do
   end subroutine cleanup_all_handles

   ! Destroy the data associated with a handle.
   ! Recover all free pointers at end of handle list.
   subroutine cleanup_handle(handle)
      integer, intent(in) :: handle

      integer :: current
      integer :: st

      if(associated(handles(handle)%rkeep)) then
         call ma86_matlab_finalise(handles(handle)%rkeep)
         deallocate(handles(handle)%rkeep)
         nullify(handles(handle)%rkeep)
      endif
      if(associated(handles(handle)%ckeep)) then
         call ma86_matlab_finalise(handles(handle)%ckeep)
         deallocate(handles(handle)%ckeep)
         nullify(handles(handle)%ckeep)
      endif
      deallocate(handles(handle)%order, stat=st)
      deallocate(handles(handle)%scaling, stat=st)

      do current = handle, 1, -1
         if(current.ne.next_handle-1) exit
         if(.not.associated(handles(current)%rkeep) .and. .not. &
               associated(handles(current)%ckeep)) then
            ! Current "last" element is unallocated, make it next available
            ! element.
            next_handle = next_handle - 1
         endif
      end do
   end subroutine cleanup_handle

end module ma86_handles


! Gateway routine
! Strip first argument and converts any handles, then calls relevant routine
subroutine mexFunction(nlhs_in, plhs, nrhs_in, prhs)
   use hsl_matlab
   use hsl_mc68_double
   use hsl_mc69_double
   use hsl_ma86_double
   use ma86_handles
   use ma86_matlab_main
   implicit none

   integer*4 :: nlhs_in, nrhs_in, tlhs, trhs
   integer(mwp_) :: plhs(*), prhs(*)
   integer(mws_) :: mwstemp

   character(len=200) :: act
   integer :: handle
   integer, dimension(1) :: thandle
   integer :: st
   integer :: n

   character(len=200) :: errmsg

   if(nrhs_in.lt.1) call MATLAB_error("Insufficient input arguments")

   if(.not. MATLAB_is_character(prhs(1))) &
      call MATLAB_error("First argument must be string")
   mwstemp = len(act)
   call MATLAB_get_string(prhs(1), act, mwstemp)

   select case(trim(act))
   case('factor')
      ! [handle, info] = ma86_expert('factor', A[, P])
      ! At least handle is required for output

      ! Setup handle and store in first output argument
      if(nlhs_in.lt.1) call MATLAB_error("Insufficient output arguments")
      handle = ma86_new_handle(MATLAB_is_complex(prhs(2)))
      plhs(1) = fortran_to_matlab(handle)

      ! Call work routine
      if(associated(handles(handle)%rkeep)) then
         call ma86_matlab_analyse_factor_real(nlhs_in-1_int4_, plhs(2), &
            nrhs_in-1_int4_, prhs(2), handles(handle)%rkeep, &
            handles(handle)%order, handles(handle)%scaling)
      else
         call ma86_matlab_analyse_factor_complex(nlhs_in-1_int4_, plhs(2), &
            nrhs_in-1_int4_, prhs(2), handles(handle)%ckeep, &
            handles(handle)%order, handles(handle)%scaling)
      endif

   case('backslash')
      ! [x, info, handle] = ma86_expert('backslash', A, b[, control, P])
      ! At least soln is required for output

      ! Setup handle and store in third output argument, if present
      handle = ma86_new_handle(MATLAB_is_complex(prhs(2)))
      if(nlhs_in.ge.3) plhs(3) = fortran_to_matlab(handle)

      ! Call work routine
      if(associated(handles(handle)%rkeep)) then
         call ma86_matlab_backslash_real(nlhs_in, plhs, nrhs_in-1_int4_, &
            prhs(2), handles(handle)%rkeep, handles(handle)%order, &
            handles(handle)%scaling)
      else
         call ma86_matlab_backslash_complex(nlhs_in, plhs, nrhs_in-1_int4_, &
            prhs(2), handles(handle)%ckeep, handles(handle)%order, &
            handles(handle)%scaling)
      endif

      ! Destroy handle if it is not being returned to the user
      if(nlhs_in.lt.3) call cleanup_handle(handle)

   case('solve')
      ! [x, info] = ma86_expert('solve', handle, b)

      ! Obtain and check handle
      if(nrhs_in.lt.2) call MATLAB_error("Insufficient input arguments")
      call matlab_to_fortran(prhs(2), handle, 'handle')
      if(handle.lt.1 .or. handle.gt.total_handles) &
         call MATLAB_error("Invalid handle")

      ! Call worker routine based on allocated parts of handle
      if(associated(handles(handle)%rkeep)) then
         call ma86_matlab_solve_real(nlhs_in, plhs(1), nrhs_in-2_int4_, &
            prhs(3), handles(handle)%rkeep, handles(handle)%order, &
            handles(handle)%scaling)
      elseif(associated(handles(handle)%ckeep)) then
         call ma86_matlab_solve_complex(nlhs_in, plhs(1), nrhs_in-2_int4_, &
            prhs(3), handles(handle)%ckeep, handles(handle)%order, &
            handles(handle)%scaling)
      else
         ! No parts allocated, probably a destroyed or unallocated handle
         call MATLAB_error("Invalid handle")
      endif
      
   case('destroy')
      ! ma86_expert('destroy', handle)

      ! Obtain and check handle
      if(nrhs_in.lt.2) call MATLAB_error("Insufficient input arguments")
      call matlab_to_fortran(prhs(2), handle, 'handle')
      if(handle.lt.1 .or. handle.gt.total_handles) &
         call MATLAB_error("Invalid handle")

      ! Destroy anything that is there
      call cleanup_handle(handle)

   case default
      write(errmsg, "(3a)") "Unrecognised action: '", trim(act), "'"
      call MATLAB_error(errmsg)
   end select

   return
   100 continue
   call MATLAB_error("Insufficient memory")
end subroutine mexFunction

