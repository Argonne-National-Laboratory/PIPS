program hsl_ma86ss
   use hsl_ma86_single
   use hsl_mc69_single
   implicit none

   integer, parameter :: wp = kind(0.0)

   type (ma86_keep)    :: keep
   type (ma86_control) :: control
   type (ma86_info)    :: info

   integer :: i, n, flag, more
   integer, dimension(:), allocatable  :: ptr, row, order
   real(wp), dimension(:), allocatable :: val, x


   ! Read the lower triangle of the matrix
   read(*,*) n
   allocate(ptr(n+1));        read(*,*) ptr(:)
   allocate(row(ptr(n+1)-1)); read(*,*) row(:)
   allocate(val(ptr(n+1)-1)); read(*,*) val(:)
   ! Read the right hand side
   allocate(x(n)); read(*,*) x(:)

   ! Use the input order
   allocate(order(n))
   do i = 1,n
      order(i) = i
   end do

   ! Uncomment the following lines to enable checking (performance overhead)
   !call mc69_verify(6, HSL_MATRIX_REAL_SYM_INDEF, n, n, ptr, row, flag, more)
   !if(flag.ne.0) then
   !   write(*,*) "Matrix not in HSL standard format. flag, more = ", flag, more
   !   stop
   !endif

   ! Analyse
   call ma86_analyse(n, ptr, row, order, keep, control, info)
   if(info%flag.lt.0) then
      write(*,*) "Failure during analyse with info%flag = ", info%flag
      stop
   endif

   ! Factor
   call ma86_factor(n, ptr, row, val, order, keep, control, info)
   if(info%flag.lt.0) then
      write(*,*) "Failure during factor with info%flag = ", info%flag
      stop
   endif

   ! Solve
   call ma86_solve(x, order, keep, control, info)
   if(info%flag.lt.0) then
      write(*,*) "Failure during solve with info%flag = ", info%flag
      stop
   endif

   write(*,'(a)') ' Computed solution:'
   write(*,'(8f10.3)') x(1:n)

   ! Finalize
   call ma86_finalise(keep, control)

end program hsl_ma86ss
