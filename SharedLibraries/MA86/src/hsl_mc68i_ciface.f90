module hsl_mc68_integer_ciface
   use iso_c_binding
   use hsl_mc68_integer, only: &
      f_mc68_control       => mc68_control,     &
      f_mc68_info          => mc68_info,        &
      f_mc68_order         => mc68_order
   implicit none

   type, bind(C) :: mc68_control
      integer(C_INT) :: f_array_in ! 0 is false, otherwise true
      integer(C_INT) :: f_array_out! 0 is false, otherwise true
      integer(C_INT) :: min_l_workspace
      integer(C_INT) :: lp
      integer(C_INT) :: wp
      integer(C_INT) :: mp
      integer(C_INT) :: nemin
      integer(C_INT) :: print_level
      integer(C_INT) :: row_full_thresh
      integer(C_INT) :: row_search
   end type mc68_control

   type, bind(C) :: mc68_info
     integer(C_INT) :: flag
     integer(C_INT) :: iostat
     integer(C_INT) :: stat
     integer(C_INT) :: out_range
     integer(C_INT) :: duplicate
     integer(C_INT) :: n_compressions
     integer(C_INT) :: n_zero_eigs
     integer(C_INT) :: l_workspace
     integer(C_INT) :: zb01_info
     integer(C_INT) :: n_dense_rows
   end type mc68_info
contains
   subroutine copy_control_in(ccontrol, fcontrol, f_array_in, f_array_out, &
         min_l_workspace)
      type(mc68_control), intent(in) :: ccontrol
      type(f_mc68_control), intent(out) :: fcontrol
      logical, intent(out) :: f_array_in
      logical, intent(out) :: f_array_out
      integer, intent(out) :: min_l_workspace

      f_array_in                 = (ccontrol%f_array_in .ne. 0)
      f_array_out                = (ccontrol%f_array_out .ne. 0)
      min_l_workspace            = ccontrol%min_l_workspace
      fcontrol%lp                = ccontrol%lp
      fcontrol%wp                = ccontrol%wp
      fcontrol%mp                = ccontrol%mp
      fcontrol%nemin             = ccontrol%nemin
      fcontrol%print_level       = ccontrol%print_level
      fcontrol%row_full_thresh   = ccontrol%row_full_thresh
      fcontrol%row_search        = ccontrol%row_search
   end subroutine copy_control_in

   subroutine copy_info_out(finfo, cinfo)
      type(f_mc68_info), intent(in) :: finfo
      type(mc68_info), intent(out) :: cinfo

      cinfo%flag              = finfo%flag
      cinfo%iostat            = finfo%iostat
      cinfo%stat              = finfo%stat
      cinfo%out_range         = finfo%out_range
      cinfo%duplicate         = finfo%duplicate
      cinfo%n_compressions    = finfo%n_compressions
      cinfo%n_zero_eigs       = finfo%n_zero_eigs
      cinfo%l_workspace       = finfo%l_workspace
      cinfo%zb01_info         = finfo%zb01_info
      cinfo%n_dense_rows      = finfo%n_dense_rows
   end subroutine copy_info_out
end module hsl_mc68_integer_ciface

subroutine mc68_default_control(ccontrol) bind(C)
   use hsl_mc68_integer_ciface
   implicit none

   type(mc68_control), intent(out) :: ccontrol

   type(f_mc68_control) :: fcontrol

   ccontrol%f_array_in        = 0 ! C array indexing for input
   ccontrol%f_array_out       = 0 ! C array indexing for output
   ccontrol%min_l_workspace   = 0 ! Equivalent to not present in Fortran
   ccontrol%lp                = fcontrol%lp
   ccontrol%wp                = fcontrol%wp
   ccontrol%mp                = fcontrol%mp
   ccontrol%nemin             = fcontrol%nemin
   ccontrol%print_level       = fcontrol%print_level
   ccontrol%row_full_thresh   = fcontrol%row_full_thresh
   ccontrol%row_search        = fcontrol%row_search
end subroutine mc68_default_control

subroutine mc68_order(ord, n, cptr, crow, cperm, ccontrol, cinfo) bind(C)
   use hsl_mc68_integer_ciface
   implicit none

   integer(C_INT), value, intent(in) :: ord
   integer(C_INT), value, intent(in) :: n
   type(C_PTR), value, intent(in) :: cptr
   type(C_PTR), value, intent(in) :: crow
   type(C_PTR), value, intent(in) :: cperm
   type(mc68_control), intent(in) :: ccontrol
   type(mc68_info), intent(out) :: cinfo

   integer, dimension(:), pointer :: fptr
   integer, dimension(:), allocatable, target :: fptr_alloc
   integer, dimension(:), pointer :: frow
   integer, dimension(:), allocatable, target :: frow_alloc
   integer, dimension(:), pointer :: fperm
   type(f_mc68_control) :: fcontrol
   type(f_mc68_info) :: finfo
   logical :: f_array_in, f_array_out
   integer :: min_l_workspace

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_array_in, f_array_out, &
      min_l_workspace)
   call C_F_POINTER(cptr, fptr, shape = (/ n+1 /))
   if(.not.f_array_in) then
      allocate(fptr_alloc(n+1))
      fptr_alloc(:) = fptr(:) + 1
      fptr => fptr_alloc
   endif
   call C_F_POINTER(crow, frow, shape = (/ fptr(n+1)-1 /))
   if(.not.f_array_in) then
      allocate(frow_alloc(fptr(n+1)-1))
      frow_alloc(:) = frow(:) + 1
      frow => frow_alloc
   endif
   call C_F_POINTER(cperm, fperm, shape = (/ n /))

   ! Call the Fortran routine
   if(min_l_workspace.le.0) then
      call f_mc68_order(ord, n, fptr, frow, fperm, fcontrol, finfo)
   else
      call f_mc68_order(ord, n, fptr, frow, fperm, fcontrol, finfo, &
         min_l_workspace=min_l_workspace)
   endif

   ! Adjust perm for C indexing
   if(.not.f_array_out) then
      ! Note: 2x2 pivoting are distinguished by negative values. However in
      ! C zero is a valid index and we cannot distinguish between +/- 0.
      ! Hence we do not offer 2x2 pivoting information if C indexing is
      ! requested.
      fperm(:) = abs(fperm(:)) - 1
   endif

   ! Copy information out to C structure
   call copy_info_out(finfo, cinfo)
end subroutine mc68_order
