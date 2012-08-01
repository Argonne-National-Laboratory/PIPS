!
! COPYRIGHT (c) 2011 Science and Technology Facilities Council
! Original date 25 Feburary 2011
!
! Written by: Jonathan Hogg and Jennifer Scott
!

module hsl_ma86d_ciface
   use iso_c_binding
   use hsl_ma86_double, only :                     &
      f_ma86_keep          => ma86_keep,           &
      f_ma86_control       => ma86_control,        &
      f_ma86_info          => ma86_info,           &
      f_ma86_analyse       => ma86_analyse,        &
      f_ma86_factor        => ma86_factor,         &
      f_ma86_factor_solve  => ma86_factor_solve,   &
      f_ma86_solve         => ma86_solve,          &
      f_ma86_finalise      => ma86_finalise,       &
      f_ma86_get_n__       => ma86_get_n__
   implicit none

   ! Data type for user controls
   type, bind(C) :: ma86_control
      ! C/Fortran interface related controls
      integer(C_INT) :: f_arrays ! 0 is false, otherwise is true
      ! Printing controls
      integer(C_INT) :: diagnostics_level
      integer(C_INT) :: unit_diagnostics
      integer(C_INT) :: unit_error
      integer(C_INT) :: unit_warning
      ! Controls used by ma86_analyse
      integer(C_INT) :: nemin
      integer(C_INT) :: nb
      ! Controls used by ma86_factor and ma86_factor_solve
      integer(C_INT) :: action ! 0 is false, otherwise is true
      integer(C_INT) :: nbi
      integer(C_INT) :: pool_size
      real(C_DOUBLE) :: small
      real(C_DOUBLE) :: static
      real(C_DOUBLE) :: u
      real(C_DOUBLE) :: umin
      integer(C_INT) :: scaling
   end type ma86_control

   !*************************************************

   ! data type for returning information to user.
   type, bind(C) :: ma86_info 
      real(C_DOUBLE)  :: detlog
      integer(C_INT)  :: detsign
      integer(C_INT)  :: flag
      integer(C_INT)  :: matrix_rank
      integer(C_INT)  :: maxdepth
      integer(C_INT)  :: num_delay
      integer(C_LONG) :: num_factor
      integer(C_LONG) :: num_flops
      integer(C_INT)  :: num_neg
      integer(C_INT)  :: num_nodes
      integer(C_INT)  :: num_nothresh
      integer(C_INT)  :: num_perturbed
      integer(C_INT)  :: num_two
      integer(C_INT)  :: pool_size
      integer(C_INT)  :: stat
      real(C_DOUBLE)  :: usmall
   end type ma86_info
contains
   subroutine copy_control_in(ccontrol, fcontrol, f_arrays)
      type(ma86_control), intent(in) :: ccontrol
      type(f_ma86_control), intent(out) :: fcontrol
      logical, intent(out) :: f_arrays

      f_arrays                   = (ccontrol%f_arrays .ne. 0)
      fcontrol%diagnostics_level = ccontrol%diagnostics_level
      fcontrol%unit_diagnostics  = ccontrol%unit_diagnostics
      fcontrol%unit_error        = ccontrol%unit_error
      fcontrol%unit_warning      = ccontrol%unit_warning
      fcontrol%nemin             = ccontrol%nemin
      fcontrol%nb                = ccontrol%nb
      fcontrol%action            = (ccontrol%action .ne. 0)
      fcontrol%nbi               = ccontrol%nbi
      fcontrol%pool_size         = ccontrol%pool_size
      fcontrol%small             = ccontrol%small
      fcontrol%static            = ccontrol%static
      fcontrol%u                 = ccontrol%u
      fcontrol%umin              = ccontrol%umin
      fcontrol%scaling           = ccontrol%scaling
   end subroutine copy_control_in

   subroutine copy_info_out(finfo, cinfo)
      type(f_ma86_info), intent(in) :: finfo
      type(ma86_info), intent(out) :: cinfo

      cinfo%detlog = finfo%detlog
      cinfo%detsign = finfo%detsign
      cinfo%flag = finfo%flag
      cinfo%matrix_rank = finfo%matrix_rank
      cinfo%maxdepth = finfo%maxdepth
      cinfo%num_delay = finfo%num_delay
      cinfo%num_factor = finfo%num_factor
      cinfo%num_flops = finfo%num_flops
      cinfo%num_neg = finfo%num_neg
      cinfo%num_nodes = finfo%num_nodes
      cinfo%num_nothresh = finfo%num_nothresh
      cinfo%num_perturbed = finfo%num_perturbed
      cinfo%num_two = finfo%num_two
      cinfo%pool_size = finfo%pool_size
      cinfo%stat = finfo%stat
      cinfo%usmall = finfo%usmall
   end subroutine copy_info_out
end module hsl_ma86d_ciface

subroutine ma86_default_control_d(ccontrol) bind(C)
   use hsl_ma86d_ciface
   implicit none

   type(ma86_control), intent(out) :: ccontrol

   type(f_ma86_control) :: fcontrol

   ccontrol%f_arrays          = 0 ! false
   ccontrol%diagnostics_level = fcontrol%diagnostics_level
   ccontrol%unit_diagnostics  = fcontrol%unit_diagnostics
   ccontrol%unit_error        = fcontrol%unit_error
   ccontrol%unit_warning      = fcontrol%unit_warning
   ccontrol%nemin             = fcontrol%nemin
   ccontrol%nb                = fcontrol%nb
   if( fcontrol%action ) then
      ccontrol%action = 1 ! true
   else
      ccontrol%action = 0 ! false
   endif
   ccontrol%nbi               = fcontrol%nbi
   ccontrol%pool_size         = fcontrol%pool_size
   ccontrol%small             = fcontrol%small
   ccontrol%static            = fcontrol%static
   ccontrol%u                 = fcontrol%u
   ccontrol%umin              = fcontrol%umin
   ccontrol%scaling           = fcontrol%scaling
end subroutine ma86_default_control_d

subroutine ma86_analyse_d(n, cptr, crow, corder, ckeep, ccontrol, cinfo) bind(C)
   use hsl_ma86d_ciface
   implicit none

   integer(C_INT), value, intent(in) :: n
   type(C_PTR), value, intent(in) :: cptr
   type(C_PTR), value, intent(in) :: crow
   type(C_PTR), value, intent(in) :: corder
   type(C_PTR), intent(inout) :: ckeep
   type(ma86_control), intent(in) :: ccontrol
   type(ma86_info), intent(inout) :: cinfo

   integer(C_INT), dimension(:), pointer :: fptr
   integer, dimension(:), allocatable, target :: fptr_alloc
   integer(C_INT), dimension(:), pointer :: frow
   integer, dimension(:), allocatable, target :: frow_alloc
   integer(C_INT), dimension(:), pointer :: forder
   integer, dimension(:), allocatable, target :: forder_alloc
   type(f_ma86_keep), pointer :: fkeep
   type(f_ma86_control) :: fcontrol
   type(f_ma86_info) :: finfo
   logical :: f_arrays

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   call C_F_POINTER(cptr, fptr, shape = (/ n+1 /) )
   if(.not.f_arrays) then
      allocate(fptr_alloc(n+1))
      fptr_alloc(:) = fptr(:) + 1
      fptr => fptr_alloc
   endif
   call C_F_POINTER(crow, frow, shape = (/ fptr(n+1)-1 /) )
   if(.not.f_arrays) then
      allocate(frow_alloc(fptr_alloc(n+1)-1))
      frow_alloc(:) = frow(:) + 1
      frow => frow_alloc
   endif
   call C_F_POINTER(corder, forder, shape = (/ n /) )
   if(.not.f_arrays) then
      allocate(forder_alloc(n))
      forder_alloc(:) = forder(:) + 1
      forder => forder_alloc
   endif

   ! Allocate space to store keep and arrange a C pointer to it
   allocate(fkeep)
   ckeep = c_loc(fkeep)

   ! Call the Fortran routine
   call f_ma86_analyse(n, fptr, frow, forder, fkeep, fcontrol, finfo)

   ! Copy information out to C structure
   call copy_info_out(finfo, cinfo)

   ! Copy order out if using C indexing
   if(.not.f_arrays) then
      call C_F_POINTER(corder, forder, shape = (/ n /) )
      forder(:) = forder_alloc(:) - 1
   endif
end subroutine ma86_analyse_d

subroutine ma86_factor_d(n, cptr, crow, cval, corder, ckeep, ccontrol, cinfo, &
      cscale) bind(C)
   use hsl_ma86d_ciface
   implicit none

   integer(C_INT), value, intent(in) :: n
   type(C_PTR), value, intent(in) :: cptr
   type(C_PTR), value, intent(in) :: crow
   type(C_PTR), value, intent(in) :: cval
   type(C_PTR), value, intent(in) :: corder
   type(C_PTR), intent(inout) :: ckeep
   type(ma86_control), intent(in) :: ccontrol
   type(ma86_info), intent(inout) :: cinfo
   type(C_PTR), value, intent(in) :: cscale

   integer(C_INT), dimension(:), pointer :: fptr
   integer, dimension(:), allocatable, target :: fptr_alloc
   integer(C_INT), dimension(:), pointer :: frow
   integer, dimension(:), allocatable, target :: frow_alloc
   real(C_DOUBLE), dimension(:), pointer :: fval
   integer(C_INT), dimension(:), pointer :: forder
   integer, dimension(:), allocatable, target :: forder_alloc
   type(f_ma86_keep), pointer :: fkeep
   type(f_ma86_control) :: fcontrol
   type(f_ma86_info) :: finfo
   logical :: f_arrays
   real(C_DOUBLE), dimension(:), pointer :: fscale

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   call C_F_POINTER(cptr, fptr, shape = (/ n+1 /) )
   if(.not.f_arrays) then
      allocate(fptr_alloc(n+1))
      fptr_alloc(:) = fptr(:) + 1
      fptr => fptr_alloc
   endif
   call C_F_POINTER(crow, frow, shape = (/ fptr(n+1)-1 /) )
   if(.not.f_arrays) then
      allocate(frow_alloc(fptr_alloc(n+1)-1))
      frow_alloc(:) = frow(:) + 1
      frow => frow_alloc
   endif
   call C_F_POINTER(cval, fval, shape = (/ fptr(n+1) /) )
   call C_F_POINTER(corder, forder, shape = (/ n /) )
   if(.not.f_arrays) then
      allocate(forder_alloc(n))
      forder_alloc(:) = forder(:) + 1
      forder => forder_alloc
   endif
   call C_F_POINTER(ckeep, fkeep)
   if(C_ASSOCIATED(cscale)) then
      call C_F_POINTER(cscale, fscale, shape = (/ n /) )
   else
      nullify(fscale)
   endif

   ! Call the Fortran routine
   if(associated(fscale)) then
      call f_ma86_factor(n, fptr, frow, fval, forder, fkeep, fcontrol, finfo, &
         scale=fscale)
   else
      call f_ma86_factor(n, fptr, frow, fval, forder, fkeep, fcontrol, finfo)
   endif

   ! Copy information out to C structure
   call copy_info_out(finfo, cinfo)
end subroutine ma86_factor_d

subroutine ma86_factor_solve_d(n, cptr, crow, cval, corder, ckeep, ccontrol, &
      cinfo, nrhs, ldx, cx, cscale) bind(C)
   use hsl_ma86d_ciface
   implicit none

   integer(C_INT), value, intent(in) :: n
   type(C_PTR), value, intent(in) :: cptr
   type(C_PTR), value, intent(in) :: crow
   type(C_PTR), value, intent(in) :: cval
   type(C_PTR), value, intent(in) :: corder
   type(C_PTR), intent(inout) :: ckeep
   type(ma86_control), intent(in) :: ccontrol
   type(ma86_info), intent(inout) :: cinfo
   integer(C_INT), value, intent(in) :: nrhs
   integer(C_INT), value, intent(in) :: ldx
   type(C_PTR), value, intent(in) :: cx
   type(C_PTR), value, intent(in) :: cscale

   integer(C_INT), dimension(:), pointer :: fptr
   integer, dimension(:), allocatable, target :: fptr_alloc
   integer(C_INT), dimension(:), pointer :: frow
   integer, dimension(:), allocatable, target :: frow_alloc
   real(C_DOUBLE), dimension(:), pointer :: fval
   integer(C_INT), dimension(:), pointer :: forder
   integer, dimension(:), allocatable, target :: forder_alloc
   real(C_DOUBLE), dimension(:,:), pointer :: fx
   type(f_ma86_keep), pointer :: fkeep
   type(f_ma86_control) :: fcontrol
   type(f_ma86_info) :: finfo
   logical :: f_arrays
   real(C_DOUBLE), dimension(:), pointer :: fscale

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   call C_F_POINTER(cptr, fptr, shape = (/ n+1 /) )
   if(.not.f_arrays) then
      allocate(fptr_alloc(n+1))
      fptr_alloc(:) = fptr(:) + 1
      fptr => fptr_alloc
   endif
   call C_F_POINTER(crow, frow, shape = (/ fptr(n+1)-1 /) )
   if(.not.f_arrays) then
      allocate(frow_alloc(fptr_alloc(n+1)-1))
      frow_alloc(:) = frow(:) + 1
      frow => frow_alloc
   endif
   call C_F_POINTER(cval, fval, shape = (/ fptr(n+1)-1 /) )
   call C_F_POINTER(corder, forder, shape = (/ n /) )
   if(.not.f_arrays) then
      allocate(forder_alloc(n))
      forder_alloc(:) = forder(:) + 1
      forder => forder_alloc
   endif
   call C_F_POINTER(ckeep, fkeep)
   call C_F_POINTER(cx, fx, shape = (/ ldx, nrhs /))
   if(C_ASSOCIATED(cscale)) then
      call C_F_POINTER(cscale, fscale, shape = (/ n /) )
   else
      nullify(fscale)
   endif

   ! Call the Fortran routine
   if(associated(fscale)) then
      call f_ma86_factor_solve(n, fptr, frow, fval, forder, fkeep, &
         fcontrol, finfo, nrhs, ldx, fx, scale=fscale)
   else
      call f_ma86_factor_solve(n, fptr, frow, fval, forder, fkeep, &
         fcontrol, finfo, nrhs, ldx, fx)
   endif

   ! Copy information out to C structure
   call copy_info_out(finfo, cinfo)
end subroutine ma86_factor_solve_d

subroutine ma86_solve_d(job, nrhs, ldx, cx, corder, ckeep, ccontrol, cinfo, &
      cscale) bind(C)
   use hsl_ma86d_ciface
   implicit none

   integer(C_INT), value, intent(in) :: job
   integer(C_INT), value, intent(in) :: nrhs
   integer(C_INT), value, intent(in) :: ldx
   type(C_PTR), value, intent(in) :: cx
   type(C_PTR), value, intent(in) :: corder
   type(C_PTR), intent(inout) :: ckeep
   type(ma86_control), intent(in) :: ccontrol
   type(ma86_info), intent(inout) :: cinfo
   type(C_PTR), value, intent(in) :: cscale ! deprecated, ignored

   integer :: n
   real(C_DOUBLE), dimension(:,:), pointer :: fx
   integer(C_INT), dimension(:), pointer :: forder
   integer, dimension(:), allocatable, target :: forder_alloc
   type(f_ma86_keep), pointer :: fkeep
   type(f_ma86_control) :: fcontrol
   type(f_ma86_info) :: finfo
   logical :: f_arrays

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   call C_F_POINTER(ckeep, fkeep)
   call C_F_POINTER(cx, fx, shape = (/ nrhs, ldx /) )
   n = f_ma86_get_n__(fkeep)
   call C_F_POINTER(corder, forder, shape = (/ n /) )
   if(.not.f_arrays) then
      allocate(forder_alloc(n))
      forder_alloc(:) = forder(:) + 1
      forder => forder_alloc
   endif
   ! Note: scale is ignored - it has been removed from Fortran version

   ! Call the Fortran routine
   call f_ma86_solve(nrhs, ldx, fx, forder, fkeep, fcontrol, finfo, job=job)

   ! Copy information out to C structure
   call copy_info_out(finfo, cinfo)
end subroutine ma86_solve_d

subroutine ma86_finalise_d(ckeep, ccontrol) bind(C)
   use hsl_ma86d_ciface
   implicit none

   type(C_PTR), intent(inout) :: ckeep
   type(ma86_control), intent(in) :: ccontrol

   type(f_ma86_keep), pointer :: fkeep
   type(f_ma86_control) :: fcontrol
   logical :: f_arrays

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   call C_F_POINTER(ckeep, fkeep)

   ! Call the Fortran routine
   call f_ma86_finalise(fkeep, fcontrol)

   ! Free memory
   deallocate(fkeep)
   ckeep = C_NULL_PTR
end subroutine ma86_finalise_d
