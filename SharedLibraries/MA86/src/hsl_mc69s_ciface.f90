module hsl_mc69_single_ciface
   use hsl_mc69_single, only:                      &
      f_mc69_verify        => mc69_verify,         &
      f_mc69_print         => mc69_print,          &
      f_mc69_cscl_clean    => mc69_cscl_clean,     &
      f_mc69_cscl_convert  => mc69_cscl_convert,   &
      f_mc69_cscu_convert  => mc69_cscu_convert,   &
      f_mc69_csclu_convert => mc69_csclu_convert,  &
      f_mc69_csrl_convert  => mc69_csrl_convert,   &
      f_mc69_csru_convert  => mc69_csru_convert,   &
      f_mc69_csrlu_convert => mc69_csrlu_convert,  &
      f_mc69_coord_convert => mc69_coord_convert,  &
      f_mc69_set_values    => mc69_set_values
   use iso_c_binding
   implicit none

   integer, parameter :: MC69_ERROR_C_MAP_TOO_SHORT = -1001
   integer, parameter :: MC69_ERROR_C_ROW_TOO_SHORT = -1002
end module

integer(C_INT) function mc69_verify_s(unit, matrix_type, findex, m, n, cptr, &
      crow, cval, more) bind(C)
   use hsl_mc69_single_ciface
   implicit none

   integer(C_INT), value, intent(in) :: unit
   integer(C_INT), value, intent(in) :: matrix_type
   integer(C_INT), value, intent(in) :: findex
   integer(C_INT), value, intent(in) :: m
   integer(C_INT), value, intent(in) :: n
   type(C_PTR), value, intent(in) :: cptr
   type(C_PTR), value, intent(in) :: crow
   type(C_PTR), value, intent(in) :: cval
   integer(C_INT), intent(out) :: more

   integer :: ne
   integer, dimension(:), pointer :: fptr
   integer, dimension(:), allocatable, target :: fptr1
   integer, dimension(:), pointer :: frow
   integer, dimension(:), allocatable, target :: frow1
   real(C_FLOAT), dimension(:), pointer :: fval

   ! Copy data in and associate pointers correctly
   call C_F_POINTER(cptr, fptr, shape = (/ n+1 /))
   if(findex.eq.0) then
      allocate(fptr1(n+1))
      fptr1(:) = fptr(:) + 1
      fptr => fptr1
   endif
   ne = fptr(n+1)-1
   call C_F_POINTER(crow, frow, shape = (/ ne /))
   if(findex.eq.0) then
      allocate(frow1(ne))
      frow1(:) = frow(:) + 1
      frow => frow1
   endif
   if(C_ASSOCIATED(cval)) then
      call C_F_POINTER(cval, fval, shape = (/ ne /))
   else
      nullify(fval)
   endif

   ! Call the Fortran routine
   if(associated(fval)) then
      call f_mc69_verify(unit, matrix_type, m, n, fptr, frow, &
         mc69_verify_s, more, val=fval)
   else
      call f_mc69_verify(unit, matrix_type, m, n, fptr, frow, &
         mc69_verify_s, more)
   endif

   ! Set more correctly wrt findex
   if(findex.eq.0) more = more - 1
end function mc69_verify_s

subroutine mc69_print_s(unit, lines, matrix_type, findex, m, n, cptr, crow, &
      cval) bind(C)
   use hsl_mc69_single_ciface
   implicit none

   integer(C_INT), value, intent(in) :: unit
   integer(C_INT), value, intent(in) :: lines
   integer(C_INT), value, intent(in) :: matrix_type
   integer(C_INT), value, intent(in) :: findex
   integer(C_INT), value, intent(in) :: m
   integer(C_INT), value, intent(in) :: n
   type(C_PTR), value, intent(in) :: cptr
   type(C_PTR), value, intent(in) :: crow
   type(C_PTR), value, intent(in) :: cval

   integer :: ne
   integer, dimension(:), pointer :: fptr
   integer, dimension(:), allocatable, target :: fptr1
   integer, dimension(:), pointer :: frow
   integer, dimension(:), allocatable, target :: frow1
   real(C_FLOAT), dimension(:), pointer :: fval

   ! Copy data in and associate pointers correctly
   call C_F_POINTER(cptr, fptr, shape = (/ n+1 /))
   if(findex.eq.0) then
      allocate(fptr1(n+1))
      fptr1(:) = fptr(:) + 1
      fptr => fptr1
   endif
   ne = fptr(n+1)-1
   call C_F_POINTER(crow, frow, shape = (/ ne /))
   if(findex.eq.0) then
      allocate(frow1(ne))
      frow1(:) = frow(:) + 1
      frow => frow1
   endif
   if(C_ASSOCIATED(cval)) then
      call C_F_POINTER(cval, fval, shape = (/ ne /))
   else
      nullify(fval)
   endif

   ! Call the Fortran routine
   if(associated(fval)) then
      call f_mc69_print(unit, lines, matrix_type, m, n, fptr, frow, val=fval, &
         cbase=(findex.eq.0))
   else
      call f_mc69_print(unit, lines, matrix_type, m, n, fptr, frow, &
         cbase=(findex.eq.0))
   endif
end subroutine mc69_print_s

subroutine mc69_set_values_s(matrix_type, lmap, cmap, cval_in, ne, cval_out) &
      bind(C)
   use hsl_mc69_single_ciface
   implicit none

   integer(C_INT), value, intent(in) :: matrix_type
   integer(C_INT), value, intent(in) :: lmap
   type(C_PTR), value, intent(in) :: cmap ! Always 1-indexed
   type(C_PTR), value, intent(in) :: cval_in
   integer(C_INT), value, intent(in) :: ne
   type(C_PTR), value, intent(in) :: cval_out

   integer :: i
   integer :: maxin
   integer(C_INT), dimension(:), pointer :: fmap
   real(C_FLOAT), dimension(:), pointer :: fval_in
   real(C_FLOAT), dimension(:), pointer :: fval_out

   ! Copy data in and associate pointers correctly
   call C_F_POINTER(cmap, fmap, shape = (/ lmap /))
   maxin = 0
   do i = 1, lmap
      maxin = max(maxin, abs(fmap(i)))
   end do
   call C_F_POINTER(cval_in, fval_in, shape = (/ maxin /))
   call C_F_POINTER(cval_out, fval_out, shape = (/ ne /))

   ! Call the Fortran routine
   call f_mc69_set_values(matrix_type, lmap, fmap, fval_in, ne, fval_out)

end subroutine mc69_set_values_s

integer(C_INT) function mc69_cscl_clean_s(unit, matrix_type, findex, m, n, &
      cptr, crow, cval, noor, ndup, clmap, cmap) bind(C)
   use hsl_mc69_single_ciface
   implicit none

   integer(C_INT), value, intent(in) :: unit
   integer(C_INT), value, intent(in) :: matrix_type
   integer(C_INT), value, intent(in) :: findex
   integer(C_INT), value, intent(in) :: m
   integer(C_INT), value, intent(in) :: n
   type(C_PTR), value, intent(in) :: cptr
   type(C_PTR), value, intent(in) :: crow
   type(C_PTR), value, intent(in) :: cval
   integer(C_INT), intent(out) :: noor
   integer(C_INT), intent(out) :: ndup
   type(C_PTR), value, intent(in) :: clmap
   type(C_PTR), value, intent(in) :: cmap ! Always 1-based regardless of findex

   integer :: ne
   integer(C_INT), dimension(:), pointer :: fptr
   integer, dimension(:), allocatable, target :: fptr1
   integer(C_INT), dimension(:), pointer :: frow
   integer, dimension(:), allocatable, target :: frow1
   real(C_FLOAT), dimension(:), pointer :: fval
   integer(C_INT), pointer :: flmap
   integer(C_INT), dimension(:), pointer :: fmap
   integer, dimension(:), allocatable :: fmap_alloc

   ! Short cut matters if only one of map or lmap is present
   if(C_ASSOCIATED(cmap) .neqv. C_ASSOCIATED(clmap)) then
      mc69_cscl_clean_s = -16
      return
   endif

   ! Copy data in and associate pointers correctly
   call C_F_POINTER(cptr, fptr, shape = (/ n+1 /))
   if(findex.eq.0) then
      allocate(fptr1(n+1))
      fptr1(:) = fptr(:) + 1
      fptr => fptr1
   endif
   ne = fptr(n+1)-1
   call C_F_POINTER(crow, frow, shape = (/ ne /))
   if(findex.eq.0) then
      allocate(frow1(ne))
      frow1(:) = frow(:) + 1
      frow => frow1
   endif
   if(C_ASSOCIATED(cval)) then
      call C_F_POINTER(cval, fval, shape = (/ ne /))
   else
      nullify(fval)
   endif
   if(C_ASSOCIATED(clmap)) then
      ! Note: If clmap is present, then so is cmap by check above
      call C_F_POINTER(clmap, flmap)
      call C_F_POINTER(cmap, fmap, shape = (/ flmap /))
   else
      nullify(flmap)
      nullify(fmap)
   endif

   ! Call the Fortran routine
   if(associated(fval)) then
      if(associated(flmap)) then
         call f_mc69_cscl_clean(matrix_type, m, n, fptr, frow, &
            mc69_cscl_clean_s, val=fval, lmap=flmap, map=fmap_alloc, lp=unit, &
            noor=noor, ndup=ndup)
      else
         call f_mc69_cscl_clean(matrix_type, m, n, fptr, frow, &
            mc69_cscl_clean_s, val=fval, lp=unit, noor=noor, ndup=ndup)
      endif
   else
      if(associated(flmap)) then
         call f_mc69_cscl_clean(matrix_type, m, n, fptr, frow, &
            mc69_cscl_clean_s, lmap=flmap, map=fmap_alloc, lp=unit, noor=noor, &
            ndup=ndup)
      else
         call f_mc69_cscl_clean(matrix_type, m, n, fptr, frow, &
            mc69_cscl_clean_s, lp=unit, noor=noor, ndup=ndup)
      endif
   endif

   ! Return without copying out data if there was an error
   if(mc69_cscl_clean_s .lt. 0) return

   ! Unpack the data back to C, note that map is always 1-based
   if(findex.eq.0) then
      call C_F_POINTER(cptr, fptr, shape = (/ n+1 /))
      fptr(:) = fptr1(:)-1
      call C_F_POINTER(crow, frow, shape = (/ ne /))
      frow(:) = frow1(:)-1
   endif
   if(associated(flmap)) then
      if(flmap .gt. size(fmap)) then
         ! map array allocated by C is too small
         mc69_cscl_clean_s = MC69_ERROR_C_MAP_TOO_SHORT
         return
      endif
      fmap(1:flmap) = fmap_alloc(1:flmap)
   endif
end function mc69_cscl_clean_s

integer(C_INT) function mc69_cscl_convert_s(unit, matrix_type, findex, m, n, &
      cptr_in, crow_in, cval_in, cptr_out, lrow, crow_out, cval_out, noor, &
      ndup, clmap, cmap) bind(C) result(flag)
   use hsl_mc69_single_ciface
   implicit none

   integer(C_INT), value, intent(in) :: unit
   integer(C_INT), value, intent(in) :: matrix_type
   integer(C_INT), value, intent(in) :: findex
   integer(C_INT), value, intent(in) :: m
   integer(C_INT), value, intent(in) :: n
   type(C_PTR), value, intent(in) :: cptr_in
   type(C_PTR), value, intent(in) :: crow_in
   type(C_PTR), value, intent(in) :: cval_in
   type(C_PTR), value, intent(in) :: cptr_out
   integer(C_INT), value, intent(in) :: lrow
   type(C_PTR), value, intent(in) :: crow_out
   type(C_PTR), value, intent(in) :: cval_out
   integer(C_INT), intent(out) :: noor
   integer(C_INT), intent(out) :: ndup
   type(C_PTR), value, intent(in) :: clmap
   type(C_PTR), value, intent(in) :: cmap

   integer :: ne
   integer(C_INT), dimension(:), pointer :: fptr_in
   integer, dimension(:), allocatable, target :: fptr_in_alloc
   integer(C_INT), dimension(:), pointer :: fptr_out
   integer(C_INT), dimension(:), pointer :: frow_in
   integer, dimension(:), allocatable, target :: frow_in_alloc
   integer(C_INT), dimension(:), pointer :: frow_out
   integer, dimension(:), allocatable :: frow_out_alloc
   real(C_FLOAT), dimension(:), pointer :: fval_in
   real(C_FLOAT), dimension(:), pointer :: fval_out
   real(C_FLOAT), dimension(:), allocatable :: fval_out_alloc
   integer(C_INT), pointer :: flmap
   integer(C_INT), dimension(:), pointer :: fmap
   integer, dimension(:), allocatable :: fmap_alloc

   ! Deal with case where only one of val_in or val_out is present
   if(C_ASSOCIATED(cval_in) .neqv. C_ASSOCIATED(cval_out)) then
      flag = -15
      return
   endif
   ! Deal with case where only one of map or lmap is present
   if(C_ASSOCIATED(cmap) .neqv. C_ASSOCIATED(clmap)) then
      flag = -16
      return
   endif

   ! Copy data in and associate pointers correctly
   call C_F_POINTER(cptr_in, fptr_in, shape = (/ n+1 /))
   if(findex.eq.0) then
      allocate(fptr_in_alloc(n+1))
      fptr_in_alloc(:) = fptr_in(:) + 1
      fptr_in => fptr_in_alloc
   endif
   ne = fptr_in(n+1)-1
   call C_F_POINTER(crow_in, frow_in, shape = (/ ne /))
   if(findex.eq.0) then
      allocate(frow_in_alloc(ne))
      frow_in_alloc(:) = frow_in(:) + 1
      frow_in => frow_in_alloc
   endif
   if(C_ASSOCIATED(cval_in)) then
      call C_F_POINTER(cval_in, fval_in, shape = (/ ne /))
      call C_F_POINTER(cval_out, fval_out, shape = (/ lrow /))
   else
      nullify(fval_in)
      nullify(fval_out)
   endif
   if(C_ASSOCIATED(clmap)) then
      ! Note: If clmap is present, then so is cmap by check above
      call C_F_POINTER(clmap, flmap)
      call C_F_POINTER(cmap, fmap, shape = (/ flmap /))
   else
      nullify(flmap)
      nullify(fmap)
   endif
   call C_F_POINTER(cptr_out, fptr_out, shape = (/ n+1 /))
   call C_F_POINTER(crow_out, frow_out, shape = (/ lrow /))

   ! Call the Fortran routine
   if(associated(fval_in)) then
      if(associated(flmap)) then
         call f_mc69_cscl_convert(matrix_type, m, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, val_in=fval_in, &
            val_out=fval_out_alloc, lmap=flmap, map=fmap_alloc, lp=unit, &
            noor=noor, ndup=ndup)
      else
         call f_mc69_cscl_convert(matrix_type, m, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, val_in=fval_in, &
            val_out=fval_out_alloc, lp=unit, noor=noor, ndup=ndup)
      endif
   else
      if(associated(flmap)) then
         call f_mc69_cscl_convert(matrix_type, m, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, lmap=flmap, map=fmap_alloc, &
            lp=unit, noor=noor, ndup=ndup)
      else
         call f_mc69_cscl_convert(matrix_type, m, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, lp=unit, noor=noor, &
            ndup=ndup)
      endif
   endif

   ! Return without copying out data if there was an error
   if(flag .lt. 0) return

   ! Unpack the data back to C, adjusting index findex if necessary
   ne = fptr_out(n+1)-1
   if(findex.eq.0) fptr_out(:) = fptr_out(:) - 1
   if(ne .gt. lrow) then
      ! row_out and val_out are too small
      flag = MC69_ERROR_C_ROW_TOO_SHORT
      return
   endif
   if(findex.eq.0) then
      frow_out(1:ne) = frow_out_alloc(1:ne) - 1
   else
      frow_out(1:ne) = frow_out_alloc(1:ne)
   endif
   if(associated(fval_out)) then
      fval_out(1:ne) = fval_out_alloc(1:ne)
   endif
   ! Note: map is always 1-based
   if(associated(flmap)) then
      if(flmap .gt. size(fmap)) then
         ! map array allocated by C is too small
         flag = MC69_ERROR_C_MAP_TOO_SHORT
         return
      endif
      fmap(1:flmap) = fmap_alloc(1:flmap)
   endif
end function mc69_cscl_convert_s

integer(C_INT) function mc69_cscu_convert_s(unit, matrix_type, findex, n, &
      cptr_in, crow_in, cval_in, cptr_out, lrow, crow_out, cval_out, noor, &
      ndup, clmap, cmap) bind(C) result(flag)
   use hsl_mc69_single_ciface
   implicit none

   integer(C_INT), value, intent(in) :: unit
   integer(C_INT), value, intent(in) :: matrix_type
   integer(C_INT), value, intent(in) :: findex
   integer(C_INT), value, intent(in) :: n
   type(C_PTR), value, intent(in) :: cptr_in
   type(C_PTR), value, intent(in) :: crow_in
   type(C_PTR), value, intent(in) :: cval_in
   type(C_PTR), value, intent(in) :: cptr_out
   integer(C_INT), value, intent(in) :: lrow
   type(C_PTR), value, intent(in) :: crow_out
   type(C_PTR), value, intent(in) :: cval_out
   integer(C_INT), intent(out) :: noor
   integer(C_INT), intent(out) :: ndup
   type(C_PTR), value, intent(in) :: clmap
   type(C_PTR), value, intent(in) :: cmap

   integer :: ne
   integer(C_INT), dimension(:), pointer :: fptr_in
   integer, dimension(:), allocatable, target :: fptr_in_alloc
   integer(C_INT), dimension(:), pointer :: fptr_out
   integer(C_INT), dimension(:), pointer :: frow_in
   integer, dimension(:), allocatable, target :: frow_in_alloc
   integer(C_INT), dimension(:), pointer :: frow_out
   integer, dimension(:), allocatable :: frow_out_alloc
   real(C_FLOAT), dimension(:), pointer :: fval_in
   real(C_FLOAT), dimension(:), pointer :: fval_out
   real(C_FLOAT), dimension(:), allocatable :: fval_out_alloc
   integer(C_INT), pointer :: flmap
   integer(C_INT), dimension(:), pointer :: fmap
   integer, dimension(:), allocatable :: fmap_alloc

   ! Deal with case where only one of val_in or val_out is present
   if(C_ASSOCIATED(cval_in) .neqv. C_ASSOCIATED(cval_out)) then
      flag = -15
      return
   endif
   ! Deal with case where only one of map or lmap is present
   if(C_ASSOCIATED(cmap) .neqv. C_ASSOCIATED(clmap)) then
      flag = -16
      return
   endif

   ! Copy data in and associate pointers correctly
   call C_F_POINTER(cptr_in, fptr_in, shape = (/ n+1 /))
   if(findex.eq.0) then
      allocate(fptr_in_alloc(n+1))
      fptr_in_alloc(:) = fptr_in(:) + 1
      fptr_in => fptr_in_alloc
   endif
   ne = fptr_in(n+1)-1
   call C_F_POINTER(crow_in, frow_in, shape = (/ ne /))
   if(findex.eq.0) then
      allocate(frow_in_alloc(ne))
      frow_in_alloc(:) = frow_in(:) + 1
      frow_in => frow_in_alloc
   endif
   if(C_ASSOCIATED(cval_in)) then
      call C_F_POINTER(cval_in, fval_in, shape = (/ ne /))
      call C_F_POINTER(cval_out, fval_out, shape = (/ lrow /))
   else
      nullify(fval_in)
      nullify(fval_out)
   endif
   if(C_ASSOCIATED(clmap)) then
      ! Note: If clmap is present, then so is cmap by check above
      call C_F_POINTER(clmap, flmap)
      call C_F_POINTER(cmap, fmap, shape = (/ flmap /))
   else
      nullify(flmap)
      nullify(fmap)
   endif
   call C_F_POINTER(cptr_out, fptr_out, shape = (/ n+1 /))
   call C_F_POINTER(crow_out, frow_out, shape = (/ lrow /))

   ! Call the Fortran routine
   if(associated(fval_in)) then
      if(associated(flmap)) then
         call f_mc69_cscu_convert(matrix_type, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, val_in=fval_in, &
            val_out=fval_out_alloc, lmap=flmap, map=fmap_alloc, lp=unit, &
            noor=noor, ndup=ndup)
      else
         call f_mc69_cscu_convert(matrix_type, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, val_in=fval_in, &
            val_out=fval_out_alloc, lp=unit, noor=noor, ndup=ndup)
      endif
   else
      if(associated(flmap)) then
         call f_mc69_cscu_convert(matrix_type, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, lmap=flmap, map=fmap_alloc, &
            lp=unit, noor=noor, ndup=ndup)
      else
         call f_mc69_cscu_convert(matrix_type, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, lp=unit, noor=noor, &
            ndup=ndup)
      endif
   endif

   ! Return without copying out data if there was an error
   if(flag .lt. 0) return

   ! Unpack the data back to C, adjusting index findex if necessary
   ne = fptr_out(n+1)-1
   if(findex.eq.0) fptr_out(:) = fptr_out(:) - 1
   if(ne .gt. lrow) then
      ! row_out and val_out are too small
      flag = MC69_ERROR_C_ROW_TOO_SHORT
      return
   endif
   if(findex.eq.0) then
      frow_out(1:ne) = frow_out_alloc(1:ne) - 1
   else
      frow_out(1:ne) = frow_out_alloc(1:ne)
   endif
   if(associated(fval_out)) then
      fval_out(1:ne) = fval_out_alloc(1:ne)
   endif
   ! Note: map is always 1-based
   if(associated(flmap)) then
      if(flmap .gt. size(fmap)) then
         ! map array allocated by C is too small
         flag = MC69_ERROR_C_MAP_TOO_SHORT
         return
      endif
      fmap(1:flmap) = fmap_alloc(1:flmap)
   endif
end function mc69_cscu_convert_s

integer(C_INT) function mc69_csclu_convert_s(unit, matrix_type, findex, n, &
      cptr_in, crow_in, cval_in, cptr_out, lrow, crow_out, cval_out, noor, &
      ndup, clmap, cmap) bind(C) result(flag)
   use hsl_mc69_single_ciface
   implicit none

   integer(C_INT), value, intent(in) :: unit
   integer(C_INT), value, intent(in) :: matrix_type
   integer(C_INT), value, intent(in) :: findex
   integer(C_INT), value, intent(in) :: n
   type(C_PTR), value, intent(in) :: cptr_in
   type(C_PTR), value, intent(in) :: crow_in
   type(C_PTR), value, intent(in) :: cval_in
   type(C_PTR), value, intent(in) :: cptr_out
   integer(C_INT), value, intent(in) :: lrow
   type(C_PTR), value, intent(in) :: crow_out
   type(C_PTR), value, intent(in) :: cval_out
   integer(C_INT), intent(out) :: noor
   integer(C_INT), intent(out) :: ndup
   type(C_PTR), value, intent(in) :: clmap
   type(C_PTR), value, intent(in) :: cmap

   integer :: ne
   integer(C_INT), dimension(:), pointer :: fptr_in
   integer, dimension(:), allocatable, target :: fptr_in_alloc
   integer(C_INT), dimension(:), pointer :: fptr_out
   integer(C_INT), dimension(:), pointer :: frow_in
   integer, dimension(:), allocatable, target :: frow_in_alloc
   integer(C_INT), dimension(:), pointer :: frow_out
   integer, dimension(:), allocatable :: frow_out_alloc
   real(C_FLOAT), dimension(:), pointer :: fval_in
   real(C_FLOAT), dimension(:), pointer :: fval_out
   real(C_FLOAT), dimension(:), allocatable :: fval_out_alloc
   integer(C_INT), pointer :: flmap
   integer(C_INT), dimension(:), pointer :: fmap
   integer, dimension(:), allocatable :: fmap_alloc

   ! Deal with case where only one of val_in or val_out is present
   if(C_ASSOCIATED(cval_in) .neqv. C_ASSOCIATED(cval_out)) then
      flag = -15
      return
   endif
   ! Deal with case where only one of map or lmap is present
   if(C_ASSOCIATED(cmap) .neqv. C_ASSOCIATED(clmap)) then
      flag = -16
      return
   endif

   ! Copy data in and associate pointers correctly
   call C_F_POINTER(cptr_in, fptr_in, shape = (/ n+1 /))
   if(findex.eq.0) then
      allocate(fptr_in_alloc(n+1))
      fptr_in_alloc(:) = fptr_in(:) + 1
      fptr_in => fptr_in_alloc
   endif
   ne = fptr_in(n+1)-1
   call C_F_POINTER(crow_in, frow_in, shape = (/ ne /))
   if(findex.eq.0) then
      allocate(frow_in_alloc(ne))
      frow_in_alloc(:) = frow_in(:) + 1
      frow_in => frow_in_alloc
   endif
   if(C_ASSOCIATED(cval_in)) then
      call C_F_POINTER(cval_in, fval_in, shape = (/ ne /))
      call C_F_POINTER(cval_out, fval_out, shape = (/ lrow /))
   else
      nullify(fval_in)
      nullify(fval_out)
   endif
   if(C_ASSOCIATED(clmap)) then
      ! Note: If clmap is present, then so is cmap by check above
      call C_F_POINTER(clmap, flmap)
      call C_F_POINTER(cmap, fmap, shape = (/ flmap /))
   else
      nullify(flmap)
      nullify(fmap)
   endif
   call C_F_POINTER(cptr_out, fptr_out, shape = (/ n+1 /))
   call C_F_POINTER(crow_out, frow_out, shape = (/ lrow /))

   ! Call the Fortran routine
   if(associated(fval_in)) then
      if(associated(flmap)) then
         call f_mc69_csclu_convert(matrix_type, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, val_in=fval_in, &
            val_out=fval_out_alloc, lmap=flmap, map=fmap_alloc, lp=unit, &
            noor=noor, ndup=ndup)
      else
         call f_mc69_csclu_convert(matrix_type, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, val_in=fval_in, &
            val_out=fval_out_alloc, lp=unit, noor=noor, ndup=ndup)
      endif
   else
      if(associated(flmap)) then
         call f_mc69_csclu_convert(matrix_type, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, lmap=flmap, map=fmap_alloc, &
            lp=unit, noor=noor, ndup=ndup)
      else
         call f_mc69_csclu_convert(matrix_type, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, lp=unit, noor=noor, &
            ndup=ndup)
      endif
   endif

   ! Return without copying out data if there was an error
   if(flag .lt. 0) return

   ! Unpack the data back to C, adjusting index findex if necessary
   ne = fptr_out(n+1)-1
   if(findex.eq.0) fptr_out(:) = fptr_out(:) - 1
   if(ne .gt. lrow) then
      ! row_out and val_out are too small
      flag = MC69_ERROR_C_ROW_TOO_SHORT
      return
   endif
   if(findex.eq.0) then
      frow_out(1:ne) = frow_out_alloc(1:ne) - 1
   else
      frow_out(1:ne) = frow_out_alloc(1:ne)
   endif
   if(associated(fval_out)) then
      fval_out(1:ne) = fval_out_alloc(1:ne)
   endif
   ! Note: map is always 1-based
   if(associated(flmap)) then
      if(flmap .gt. size(fmap)) then
         ! map array allocated by C is too small
         flag = MC69_ERROR_C_MAP_TOO_SHORT
         return
      endif
      fmap(1:flmap) = fmap_alloc(1:flmap)
   endif
end function mc69_csclu_convert_s

integer(C_INT) function mc69_csrl_convert_s(unit, matrix_type, findex, m, n, &
      cptr_in, crow_in, cval_in, cptr_out, lrow, crow_out, cval_out, noor, &
      ndup, clmap, cmap) bind(C) result(flag)
   use hsl_mc69_single_ciface
   implicit none

   integer(C_INT), value, intent(in) :: unit
   integer(C_INT), value, intent(in) :: matrix_type
   integer(C_INT), value, intent(in) :: findex
   integer(C_INT), value, intent(in) :: m
   integer(C_INT), value, intent(in) :: n
   type(C_PTR), value, intent(in) :: cptr_in
   type(C_PTR), value, intent(in) :: crow_in
   type(C_PTR), value, intent(in) :: cval_in
   type(C_PTR), value, intent(in) :: cptr_out
   integer(C_INT), value, intent(in) :: lrow
   type(C_PTR), value, intent(in) :: crow_out
   type(C_PTR), value, intent(in) :: cval_out
   integer(C_INT), intent(out) :: noor
   integer(C_INT), intent(out) :: ndup
   type(C_PTR), value, intent(in) :: clmap
   type(C_PTR), value, intent(in) :: cmap

   integer :: ne
   integer(C_INT), dimension(:), pointer :: fptr_in
   integer, dimension(:), allocatable, target :: fptr_in_alloc
   integer(C_INT), dimension(:), pointer :: fptr_out
   integer(C_INT), dimension(:), pointer :: frow_in
   integer, dimension(:), allocatable, target :: frow_in_alloc
   integer(C_INT), dimension(:), pointer :: frow_out
   integer, dimension(:), allocatable :: frow_out_alloc
   real(C_FLOAT), dimension(:), pointer :: fval_in
   real(C_FLOAT), dimension(:), pointer :: fval_out
   real(C_FLOAT), dimension(:), allocatable :: fval_out_alloc
   integer(C_INT), pointer :: flmap
   integer(C_INT), dimension(:), pointer :: fmap
   integer, dimension(:), allocatable :: fmap_alloc

   ! Deal with case where only one of val_in or val_out is present
   if(C_ASSOCIATED(cval_in) .neqv. C_ASSOCIATED(cval_out)) then
      flag = -15
      return
   endif
   ! Deal with case where only one of map or lmap is present
   if(C_ASSOCIATED(cmap) .neqv. C_ASSOCIATED(clmap)) then
      flag = -16
      return
   endif

   ! Copy data in and associate pointers correctly
   call C_F_POINTER(cptr_in, fptr_in, shape = (/ m+1 /))
   if(findex.eq.0) then
      allocate(fptr_in_alloc(m+1))
      fptr_in_alloc(:) = fptr_in(:) + 1
      fptr_in => fptr_in_alloc
   endif
   ne = fptr_in(m+1)-1
   call C_F_POINTER(crow_in, frow_in, shape = (/ ne /))
   if(findex.eq.0) then
      allocate(frow_in_alloc(ne))
      frow_in_alloc(:) = frow_in(:) + 1
      frow_in => frow_in_alloc
   endif
   if(C_ASSOCIATED(cval_in)) then
      call C_F_POINTER(cval_in, fval_in, shape = (/ ne /))
      call C_F_POINTER(cval_out, fval_out, shape = (/ lrow /))
   else
      nullify(fval_in)
      nullify(fval_out)
   endif
   if(C_ASSOCIATED(clmap)) then
      ! Note: If clmap is present, then so is cmap by check above
      call C_F_POINTER(clmap, flmap)
      call C_F_POINTER(cmap, fmap, shape = (/ flmap /))
   else
      nullify(flmap)
      nullify(fmap)
   endif
   call C_F_POINTER(cptr_out, fptr_out, shape = (/ n+1 /))
   call C_F_POINTER(crow_out, frow_out, shape = (/ lrow /))

   ! Call the Fortran routine
   if(associated(fval_in)) then
      if(associated(flmap)) then
         call f_mc69_csrl_convert(matrix_type, m, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, val_in=fval_in, &
            val_out=fval_out_alloc, lmap=flmap, map=fmap_alloc, lp=unit, &
            noor=noor, ndup=ndup)
      else
         call f_mc69_csrl_convert(matrix_type, m, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, val_in=fval_in, &
            val_out=fval_out_alloc, lp=unit, noor=noor, ndup=ndup)
      endif
   else
      if(associated(flmap)) then
         call f_mc69_csrl_convert(matrix_type, m, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, lmap=flmap, map=fmap_alloc, &
            lp=unit, noor=noor, ndup=ndup)
      else
         call f_mc69_csrl_convert(matrix_type, m, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, lp=unit, noor=noor, &
            ndup=ndup)
      endif
   endif

   ! Return without copying out data if there was an error
   if(flag .lt. 0) return

   ! Unpack the data back to C, adjusting index findex if necessary
   ne = fptr_out(n+1)-1
   if(findex.eq.0) fptr_out(:) = fptr_out(:) - 1
   if(ne .gt. lrow) then
      ! row_out and val_out are too small
      flag = MC69_ERROR_C_ROW_TOO_SHORT
      return
   endif
   if(findex.eq.0) then
      frow_out(1:ne) = frow_out_alloc(1:ne) - 1
   else
      frow_out(1:ne) = frow_out_alloc(1:ne)
   endif
   if(associated(fval_out)) then
      fval_out(1:ne) = fval_out_alloc(1:ne)
   endif
   ! Note: map is always 1-based
   if(associated(flmap)) then
      if(flmap .gt. size(fmap)) then
         ! map array allocated by C is too small
         flag = MC69_ERROR_C_MAP_TOO_SHORT
         return
      endif
      fmap(1:flmap) = fmap_alloc(1:flmap)
   endif
end function mc69_csrl_convert_s

integer(C_INT) function mc69_csru_convert_s(unit, matrix_type, findex, n, &
      cptr_in, crow_in, cval_in, cptr_out, lrow, crow_out, cval_out, noor, &
      ndup, clmap, cmap) bind(C) result(flag)
   use hsl_mc69_single_ciface
   implicit none

   integer(C_INT), value, intent(in) :: unit
   integer(C_INT), value, intent(in) :: matrix_type
   integer(C_INT), value, intent(in) :: findex
   integer(C_INT), value, intent(in) :: n
   type(C_PTR), value, intent(in) :: cptr_in
   type(C_PTR), value, intent(in) :: crow_in
   type(C_PTR), value, intent(in) :: cval_in
   type(C_PTR), value, intent(in) :: cptr_out
   integer(C_INT), value, intent(in) :: lrow
   type(C_PTR), value, intent(in) :: crow_out
   type(C_PTR), value, intent(in) :: cval_out
   integer(C_INT), intent(out) :: noor
   integer(C_INT), intent(out) :: ndup
   type(C_PTR), value, intent(in) :: clmap
   type(C_PTR), value, intent(in) :: cmap

   integer :: ne
   integer(C_INT), dimension(:), pointer :: fptr_in
   integer, dimension(:), allocatable, target :: fptr_in_alloc
   integer(C_INT), dimension(:), pointer :: fptr_out
   integer(C_INT), dimension(:), pointer :: frow_in
   integer, dimension(:), allocatable, target :: frow_in_alloc
   integer(C_INT), dimension(:), pointer :: frow_out
   integer, dimension(:), allocatable :: frow_out_alloc
   real(C_FLOAT), dimension(:), pointer :: fval_in
   real(C_FLOAT), dimension(:), pointer :: fval_out
   real(C_FLOAT), dimension(:), allocatable :: fval_out_alloc
   integer(C_INT), pointer :: flmap
   integer(C_INT), dimension(:), pointer :: fmap
   integer, dimension(:), allocatable :: fmap_alloc

   ! Deal with case where only one of val_in or val_out is present
   if(C_ASSOCIATED(cval_in) .neqv. C_ASSOCIATED(cval_out)) then
      flag = -15
      return
   endif
   ! Deal with case where only one of map or lmap is present
   if(C_ASSOCIATED(cmap) .neqv. C_ASSOCIATED(clmap)) then
      flag = -16
      return
   endif

   ! Copy data in and associate pointers correctly
   call C_F_POINTER(cptr_in, fptr_in, shape = (/ n+1 /))
   if(findex.eq.0) then
      allocate(fptr_in_alloc(n+1))
      fptr_in_alloc(:) = fptr_in(:) + 1
      fptr_in => fptr_in_alloc
   endif
   ne = fptr_in(n+1)-1
   call C_F_POINTER(crow_in, frow_in, shape = (/ ne /))
   if(findex.eq.0) then
      allocate(frow_in_alloc(ne))
      frow_in_alloc(:) = frow_in(:) + 1
      frow_in => frow_in_alloc
   endif
   if(C_ASSOCIATED(cval_in)) then
      call C_F_POINTER(cval_in, fval_in, shape = (/ ne /))
      call C_F_POINTER(cval_out, fval_out, shape = (/ lrow /))
   else
      nullify(fval_in)
      nullify(fval_out)
   endif
   if(C_ASSOCIATED(clmap)) then
      ! Note: If clmap is present, then so is cmap by check above
      call C_F_POINTER(clmap, flmap)
      call C_F_POINTER(cmap, fmap, shape = (/ flmap /))
   else
      nullify(flmap)
      nullify(fmap)
   endif
   call C_F_POINTER(cptr_out, fptr_out, shape = (/ n+1 /))
   call C_F_POINTER(crow_out, frow_out, shape = (/ lrow /))

   ! Call the Fortran routine
   if(associated(fval_in)) then
      if(associated(flmap)) then
         call f_mc69_csru_convert(matrix_type, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, val_in=fval_in, &
            val_out=fval_out_alloc, lmap=flmap, map=fmap_alloc, lp=unit, &
            noor=noor, ndup=ndup)
      else
         call f_mc69_csru_convert(matrix_type, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, val_in=fval_in, &
            val_out=fval_out_alloc, lp=unit, noor=noor, ndup=ndup)
      endif
   else
      if(associated(flmap)) then
         call f_mc69_csru_convert(matrix_type, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, lmap=flmap, map=fmap_alloc, &
            lp=unit, noor=noor, ndup=ndup)
      else
         call f_mc69_csru_convert(matrix_type, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, lp=unit, noor=noor, &
            ndup=ndup)
      endif
   endif

   ! Return without copying out data if there was an error
   if(flag .lt. 0) return

   ! Unpack the data back to C, adjusting index findex if necessary
   ne = fptr_out(n+1)-1
   if(findex.eq.0) fptr_out(:) = fptr_out(:) - 1
   if(ne .gt. lrow) then
      ! row_out and val_out are too small
      flag = MC69_ERROR_C_ROW_TOO_SHORT
      return
   endif
   if(findex.eq.0) then
      frow_out(1:ne) = frow_out_alloc(1:ne) - 1
   else
      frow_out(1:ne) = frow_out_alloc(1:ne)
   endif
   if(associated(fval_out)) then
      fval_out(1:ne) = fval_out_alloc(1:ne)
   endif
   ! Note: map is always 1-based
   if(associated(flmap)) then
      if(flmap .gt. size(fmap)) then
         ! map array allocated by C is too small
         flag = MC69_ERROR_C_MAP_TOO_SHORT
         return
      endif
      fmap(1:flmap) = fmap_alloc(1:flmap)
   endif
end function mc69_csru_convert_s

integer(C_INT) function mc69_csrlu_convert_s(unit, matrix_type, findex, n, &
      cptr_in, crow_in, cval_in, cptr_out, lrow, crow_out, cval_out, noor, &
      ndup, clmap, cmap) bind(C) result(flag)
   use hsl_mc69_single_ciface
   implicit none

   integer(C_INT), value, intent(in) :: unit
   integer(C_INT), value, intent(in) :: matrix_type
   integer(C_INT), value, intent(in) :: findex
   integer(C_INT), value, intent(in) :: n
   type(C_PTR), value, intent(in) :: cptr_in
   type(C_PTR), value, intent(in) :: crow_in
   type(C_PTR), value, intent(in) :: cval_in
   type(C_PTR), value, intent(in) :: cptr_out
   integer(C_INT), value, intent(in) :: lrow
   type(C_PTR), value, intent(in) :: crow_out
   type(C_PTR), value, intent(in) :: cval_out
   integer(C_INT), intent(out) :: noor
   integer(C_INT), intent(out) :: ndup
   type(C_PTR), value, intent(in) :: clmap
   type(C_PTR), value, intent(in) :: cmap

   integer :: ne
   integer(C_INT), dimension(:), pointer :: fptr_in
   integer, dimension(:), allocatable, target :: fptr_in_alloc
   integer(C_INT), dimension(:), pointer :: fptr_out
   integer(C_INT), dimension(:), pointer :: frow_in
   integer, dimension(:), allocatable, target :: frow_in_alloc
   integer(C_INT), dimension(:), pointer :: frow_out
   integer, dimension(:), allocatable :: frow_out_alloc
   real(C_FLOAT), dimension(:), pointer :: fval_in
   real(C_FLOAT), dimension(:), pointer :: fval_out
   real(C_FLOAT), dimension(:), allocatable :: fval_out_alloc
   integer(C_INT), pointer :: flmap
   integer(C_INT), dimension(:), pointer :: fmap
   integer, dimension(:), allocatable :: fmap_alloc

   ! Deal with case where only one of val_in or val_out is present
   if(C_ASSOCIATED(cval_in) .neqv. C_ASSOCIATED(cval_out)) then
      flag = -15
      return
   endif
   ! Deal with case where only one of map or lmap is present
   if(C_ASSOCIATED(cmap) .neqv. C_ASSOCIATED(clmap)) then
      flag = -16
      return
   endif

   ! Copy data in and associate pointers correctly
   call C_F_POINTER(cptr_in, fptr_in, shape = (/ n+1 /))
   if(findex.eq.0) then
      allocate(fptr_in_alloc(n+1))
      fptr_in_alloc(:) = fptr_in(:) + 1
      fptr_in => fptr_in_alloc
   endif
   ne = fptr_in(n+1)-1
   call C_F_POINTER(crow_in, frow_in, shape = (/ ne /))
   if(findex.eq.0) then
      allocate(frow_in_alloc(ne))
      frow_in_alloc(:) = frow_in(:) + 1
      frow_in => frow_in_alloc
   endif
   if(C_ASSOCIATED(cval_in)) then
      call C_F_POINTER(cval_in, fval_in, shape = (/ ne /))
      call C_F_POINTER(cval_out, fval_out, shape = (/ lrow /))
   else
      nullify(fval_in)
      nullify(fval_out)
   endif
   if(C_ASSOCIATED(clmap)) then
      ! Note: If clmap is present, then so is cmap by check above
      call C_F_POINTER(clmap, flmap)
      call C_F_POINTER(cmap, fmap, shape = (/ flmap /))
   else
      nullify(flmap)
      nullify(fmap)
   endif
   call C_F_POINTER(cptr_out, fptr_out, shape = (/ n+1 /))
   call C_F_POINTER(crow_out, frow_out, shape = (/ lrow /))

   ! Call the Fortran routine
   if(associated(fval_in)) then
      if(associated(flmap)) then
         call f_mc69_csrlu_convert(matrix_type, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, val_in=fval_in, &
            val_out=fval_out_alloc, lmap=flmap, map=fmap_alloc, lp=unit, &
            noor=noor, ndup=ndup)
      else
         call f_mc69_csrlu_convert(matrix_type, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, val_in=fval_in, &
            val_out=fval_out_alloc, lp=unit, noor=noor, ndup=ndup)
      endif
   else
      if(associated(flmap)) then
         call f_mc69_csrlu_convert(matrix_type, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, lmap=flmap, map=fmap_alloc, &
            lp=unit, noor=noor, ndup=ndup)
      else
         call f_mc69_csrlu_convert(matrix_type, n, fptr_in, frow_in, &
            fptr_out, frow_out_alloc, flag, lp=unit, noor=noor, &
            ndup=ndup)
      endif
   endif

   ! Return without copying out data if there was an error
   if(flag .lt. 0) return

   ! Unpack the data back to C, adjusting index findex if necessary
   ne = fptr_out(n+1)-1
   if(findex.eq.0) fptr_out(:) = fptr_out(:) - 1
   if(ne .gt. lrow) then
      ! row_out and val_out are too small
      flag = MC69_ERROR_C_ROW_TOO_SHORT
      return
   endif
   if(findex.eq.0) then
      frow_out(1:ne) = frow_out_alloc(1:ne) - 1
   else
      frow_out(1:ne) = frow_out_alloc(1:ne)
   endif
   if(associated(fval_out)) then
      fval_out(1:ne) = fval_out_alloc(1:ne)
   endif
   ! Note: map is always 1-based
   if(associated(flmap)) then
      if(flmap .gt. size(fmap)) then
         ! map array allocated by C is too small
         flag = MC69_ERROR_C_MAP_TOO_SHORT
         return
      endif
      fmap(1:flmap) = fmap_alloc(1:flmap)
   endif
end function mc69_csrlu_convert_s

integer(C_INT) function mc69_coord_convert_s(unit, matrix_type, findex, m, n,  &
      ne, crow_in, ccol_in, cval_in, cptr_out, lrow, crow_out, cval_out, noor, &
      ndup, clmap, cmap) bind(C) result(flag)
   use hsl_mc69_single_ciface
   implicit none

   integer(C_INT), value, intent(in) :: unit
   integer(C_INT), value, intent(in) :: matrix_type
   integer(C_INT), value, intent(in) :: findex
   integer(C_INT), value, intent(in) :: m
   integer(C_INT), value, intent(in) :: n
   integer(C_INT), value, intent(in) :: ne
   type(C_PTR), value, intent(in) :: crow_in
   type(C_PTR), value, intent(in) :: ccol_in
   type(C_PTR), value, intent(in) :: cval_in
   type(C_PTR), value, intent(in) :: cptr_out
   integer(C_INT), value, intent(in) :: lrow
   type(C_PTR), value, intent(in) :: crow_out
   type(C_PTR), value, intent(in) :: cval_out
   integer(C_INT), intent(out) :: noor
   integer(C_INT), intent(out) :: ndup
   type(C_PTR), value, intent(in) :: clmap
   type(C_PTR), value, intent(in) :: cmap

   integer :: ne_out
   integer(C_INT), dimension(:), pointer :: frow_in
   integer, dimension(:), allocatable, target :: frow_in_alloc
   integer(C_INT), dimension(:), pointer :: fptr_out
   integer(C_INT), dimension(:), pointer :: fcol_in
   integer, dimension(:), allocatable, target :: fcol_in_alloc
   integer(C_INT), dimension(:), pointer :: frow_out
   integer, dimension(:), allocatable :: frow_out_alloc
   real(C_FLOAT), dimension(:), pointer :: fval_in
   real(C_FLOAT), dimension(:), pointer :: fval_out
   real(C_FLOAT), dimension(:), allocatable :: fval_out_alloc
   integer(C_INT), pointer :: flmap
   integer(C_INT), dimension(:), pointer :: fmap
   integer, dimension(:), allocatable :: fmap_alloc

   ! Deal with case where only one of val_in or val_out is present
   if(C_ASSOCIATED(cval_in) .neqv. C_ASSOCIATED(cval_out)) then
      flag = -15
      return
   endif
   ! Deal with case where only one of map or lmap is present
   if(C_ASSOCIATED(cmap) .neqv. C_ASSOCIATED(clmap)) then
      flag = -16
      return
   endif

   ! Copy data in and associate pointers correctly
   call C_F_POINTER(crow_in, frow_in, shape = (/ ne /))
   if(findex.eq.0) then
      allocate(frow_in_alloc(ne))
      frow_in_alloc(:) = frow_in(:) + 1
      frow_in => frow_in_alloc
   endif
   call C_F_POINTER(ccol_in, fcol_in, shape = (/ ne /))
   if(findex.eq.0) then
      allocate(fcol_in_alloc(ne))
      fcol_in_alloc(:) = fcol_in(:) + 1
      fcol_in => fcol_in_alloc
   endif
   if(C_ASSOCIATED(cval_in)) then
      call C_F_POINTER(cval_in, fval_in, shape = (/ ne /))
      call C_F_POINTER(cval_out, fval_out, shape = (/ lrow /))
   else
      nullify(fval_in)
      nullify(fval_out)
   endif
   if(C_ASSOCIATED(clmap)) then
      ! Note: If clmap is present, then so is cmap by check above
      call C_F_POINTER(clmap, flmap)
      call C_F_POINTER(cmap, fmap, shape = (/ flmap /))
   else
      nullify(flmap)
      nullify(fmap)
   endif
   call C_F_POINTER(cptr_out, fptr_out, shape = (/ n+1 /))
   call C_F_POINTER(crow_out, frow_out, shape = (/ lrow /))

   ! Call the Fortran routine
   if(associated(fval_in)) then
      if(associated(flmap)) then
         call f_mc69_coord_convert(matrix_type, m, n, ne, frow_in, fcol_in, &
            fptr_out, frow_out_alloc, flag, val_in=fval_in, &
            val_out=fval_out_alloc, lmap=flmap, map=fmap_alloc, lp=unit, &
            noor=noor, ndup=ndup)
      else
         call f_mc69_coord_convert(matrix_type, m, n, ne, frow_in, fcol_in, &
            fptr_out, frow_out_alloc, flag, val_in=fval_in, &
            val_out=fval_out_alloc, lp=unit, noor=noor, ndup=ndup)
      endif
   else
      if(associated(flmap)) then
         call f_mc69_coord_convert(matrix_type, m, n, ne, frow_in, fcol_in, &
            fptr_out, frow_out_alloc, flag, lmap=flmap, map=fmap_alloc, &
            lp=unit, noor=noor, ndup=ndup)
      else
         call f_mc69_coord_convert(matrix_type, m, n, ne, frow_in, fcol_in, &
            fptr_out, frow_out_alloc, flag, lp=unit, noor=noor, &
            ndup=ndup)
      endif
   endif

   ! Return without copying out data if there was an error
   if(flag .lt. 0) return

   ! Unpack the data back to C, adjusting index findex if necessary
   ne_out = fptr_out(n+1)-1
   if(findex.eq.0) fptr_out(:) = fptr_out(:) - 1
   if(ne_out .gt. lrow) then
      ! row_out and val_out are too small
      flag = MC69_ERROR_C_ROW_TOO_SHORT
      return
   endif
   if(findex.eq.0) then
      frow_out(1:ne) = frow_out_alloc(1:ne) - 1
   else
      frow_out(1:ne) = frow_out_alloc(1:ne)
   endif
   if(associated(fval_out)) then
      fval_out(1:ne_out) = fval_out_alloc(1:ne_out)
   endif
   ! Note: map is always 1-based
   if(associated(flmap)) then
      if(flmap .gt. size(fmap)) then
         ! map array allocated by C is too small
         flag = MC69_ERROR_C_MAP_TOO_SHORT
         return
      endif
      fmap(1:flmap) = fmap_alloc(1:flmap)
   endif
end function mc69_coord_convert_s
