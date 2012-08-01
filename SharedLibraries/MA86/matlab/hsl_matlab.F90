!-*-*-*-*-*-*-*-*-        HSL_MATLAB MODULE        -*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  Copyright (c) 2010,2011
!  Science and Technology Facilities Council (STFC)
!  Additional authors: Sue Dollar, Mario Arioli, Jonathan Hogg

!  Version 2.1.0     27 September 2011
!  History - See ChangeLog

! Warning: mathworks only support gfortran under linux (g95 prior to 2011a).
! Ensure that you use the support version as well.
! Use with unsupported compilers at your own risk!

#include <fintrf.h>
    MODULE HSL_MATLAB

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: MATLAB_copy_from_ptr,                                          &
                MATLAB_copy_to_ptr,                                            &
                MATLAB_create_integer,                                         &
                MATLAB_create_real,                                            &
                MATLAB_create_complex,                                         &
                MATLAB_create_sparse,                                          &
                MATLAB_create_structure_expert,                                &
                MATLAB_create_substructure,                                    &
                MATLAB_create_integer_component,                               &
                MATLAB_create_real_component,                                  &
                MATLAB_create_character_component,                             &
                MATLAB_create_logical_component,                               &
                MATLAB_get_field,                                              &
                MATLAB_get_field_name_by_no,                                   &
                MATLAB_get_ir,                                                 &
                MATLAB_get_jc,                                                 &
                MATLAB_get_m,                                                  &
                MATLAB_get_n,                                                  &
                MATLAB_get_no_fields,                                          &
                MATLAB_get_nzmax,                                              &
                MATLAB_get_ptr,                                                &
                MATLAB_get_imag_ptr,                                           &
                MATLAB_get_string,                                             &
                MATLAB_get_value,                                              &
                MATLAB_is_character,                                           &
                MATLAB_is_sparse,                                              &
                MATLAB_is_numeric,                                             &
                MATLAB_is_structure,                                           &
                MATLAB_is_complex,                                             &
                MATLAB_error,                                                  &
                MATLAB_context_error,                                          &
                MATLAB_warning

      ! Simple interface
      public :: matlab_to_fortran,           &
                fortran_to_matlab,           &
                matlab_create_structure,     &
                matlab_set_field

      ! Expose a few mx call interfaces directly for complicated things
      public :: mxSetField

      ! Define some kind parameters
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      mwPointer :: dummy_mwPointer__
      mwSize :: dummy_mwSize__
      mwIndex :: dummy_mwIndex__
      integer*1 :: dummy_int1__
      integer*2 :: dummy_int2__
      integer*4 :: dummy_int4__
      integer*8 :: dummy_int8__
      real*4 :: dummy_real4__
      real*8 :: dummy_real8__
      integer :: dummy_default_integer__
      real :: dummy_default_real__
      integer, parameter, public :: mwp_ = kind(dummy_mwPointer__)
      integer, parameter, public :: mws_ = kind(dummy_mwSize__)
      integer, parameter, public :: mwi_ = kind(dummy_mwIndex__)
      integer, parameter, public :: int1_ = kind(dummy_int1__)
      integer, parameter, public :: int2_ = kind(dummy_int2__)
      integer, parameter, public :: int4_ = kind(dummy_int4__)
      integer, parameter, public :: int8_ = kind(dummy_int8__)
      integer, parameter, public :: real4_ = kind(dummy_real4__)
      integer, parameter, public :: real8_ = kind(dummy_real8__)
      integer, parameter, public :: di_ = kind(dummy_default_integer__)
      integer, parameter, public :: dr_ = kind(dummy_default_real__)

!---------------------------------
!   I n t e r f a c e  B l o c k s
!---------------------------------

      INTERFACE MATLAB_error
        SUBROUTINE mexErrMsgTxt(errormsg)
          character*(*):: errormsg
        END SUBROUTINE mexErrMsgTxt
      END INTERFACE MATLAB_error

      INTERFACE MATLAB_context_error
        SUBROUTINE mexErrMsgIdAndTxt(errorid, errormsg)
          character*(*):: errorid, errormsg
        END SUBROUTINE mexErrMsgIdAndTxt
      END INTERFACE MATLAB_context_error

      INTERFACE MATLAB_create_integer_component
        MODULE PROCEDURE hslmxCreateIntegerComponent,                          &
                         hslmxCreateIntegerArrayComponent,                     &
                         hslmxCreateIntegerMatrixComponent
      END INTERFACE MATLAB_create_integer_component

      INTERFACE MATLAB_create_real_component
        MODULE PROCEDURE hslmxCreateRealComponent,                             &
                         hslmxCreateRealArrayComponent,                        &
                         hslmxCreateRealMatrixComponent
      END INTERFACE MATLAB_create_real_component

      INTERFACE MATLAB_get_value
        MODULE PROCEDURE hslmxGetInteger,                                      &
                         hslmxGetReal,                                         &
                         hslmxGetLogical,                                      &
                         hslmxGetCharacter
      END INTERFACE MATLAB_get_value

      INTERFACE matlab_to_fortran
         MODULE PROCEDURE matlab_to_fortran_scalar_int4
         MODULE PROCEDURE matlab_to_fortran_scalar_int8
         MODULE PROCEDURE matlab_to_fortran_scalar_real4
         MODULE PROCEDURE matlab_to_fortran_scalar_real8
         MODULE PROCEDURE matlab_to_fortran_scalar_complex8
         MODULE PROCEDURE matlab_to_fortran_scalar_complex16
         MODULE PROCEDURE matlab_to_fortran_array1d_int4
         MODULE PROCEDURE matlab_to_fortran_array1d_int8
         MODULE PROCEDURE matlab_to_fortran_array1d_real4
         MODULE PROCEDURE matlab_to_fortran_array1d_real8
         MODULE PROCEDURE matlab_to_fortran_array1d_complex8
         MODULE PROCEDURE matlab_to_fortran_array1d_complex16
         MODULE PROCEDURE matlab_to_fortran_array2d_int4
         MODULE PROCEDURE matlab_to_fortran_array2d_int8
         MODULE PROCEDURE matlab_to_fortran_array2d_real4
         MODULE PROCEDURE matlab_to_fortran_array2d_real8
         MODULE PROCEDURE matlab_to_fortran_array2d_complex8
         MODULE PROCEDURE matlab_to_fortran_array2d_complex16
         MODULE PROCEDURE matlab_to_fortran_sparse_int4_int4
         MODULE PROCEDURE matlab_to_fortran_sparse_int4_int8
         MODULE PROCEDURE matlab_to_fortran_sparse_int4_real4
         MODULE PROCEDURE matlab_to_fortran_sparse_int4_real8
         MODULE PROCEDURE matlab_to_fortran_sparse_int4_complex8
         MODULE PROCEDURE matlab_to_fortran_sparse_int4_complex16
         MODULE PROCEDURE matlab_to_fortran_sparse_int8_int4
         MODULE PROCEDURE matlab_to_fortran_sparse_int8_int8
         MODULE PROCEDURE matlab_to_fortran_sparse_int8_real4
         MODULE PROCEDURE matlab_to_fortran_sparse_int8_real8
         MODULE PROCEDURE matlab_to_fortran_sparse_int8_complex8
         MODULE PROCEDURE matlab_to_fortran_sparse_int8_complex16
      END INTERFACE

      INTERFACE fortran_to_matlab
         MODULE PROCEDURE fortran_to_matlab_scalar_int4
         MODULE PROCEDURE fortran_to_matlab_array1d_int4
         MODULE PROCEDURE fortran_to_matlab_array2d_int4
         MODULE PROCEDURE fortran_to_matlab_sparse_int4_int4
         MODULE PROCEDURE fortran_to_matlab_sparse_int8_int4
         MODULE PROCEDURE fortran_to_matlab_scalar_int8
         MODULE PROCEDURE fortran_to_matlab_array1d_int8
         MODULE PROCEDURE fortran_to_matlab_array2d_int8
         MODULE PROCEDURE fortran_to_matlab_sparse_int4_int8
         MODULE PROCEDURE fortran_to_matlab_sparse_int8_int8
         MODULE PROCEDURE fortran_to_matlab_scalar_real4
         MODULE PROCEDURE fortran_to_matlab_array1d_real4
         MODULE PROCEDURE fortran_to_matlab_array2d_real4
         MODULE PROCEDURE fortran_to_matlab_sparse_int4_real4
         MODULE PROCEDURE fortran_to_matlab_sparse_int8_real4
         MODULE PROCEDURE fortran_to_matlab_scalar_real8
         MODULE PROCEDURE fortran_to_matlab_array1d_real8
         MODULE PROCEDURE fortran_to_matlab_array2d_real8
         MODULE PROCEDURE fortran_to_matlab_sparse_int4_real8
         MODULE PROCEDURE fortran_to_matlab_sparse_int8_real8
         MODULE PROCEDURE fortran_to_matlab_scalar_complex8
         MODULE PROCEDURE fortran_to_matlab_array1d_complex8
         MODULE PROCEDURE fortran_to_matlab_array2d_complex8
         MODULE PROCEDURE fortran_to_matlab_sparse_int4_complex8
         MODULE PROCEDURE fortran_to_matlab_sparse_int8_complex8
         MODULE PROCEDURE fortran_to_matlab_scalar_complex16
         MODULE PROCEDURE fortran_to_matlab_array1d_complex16
         MODULE PROCEDURE fortran_to_matlab_array2d_complex16
         MODULE PROCEDURE fortran_to_matlab_sparse_int4_complex16
         MODULE PROCEDURE fortran_to_matlab_sparse_int8_complex16
      END INTERFACE fortran_to_matlab

      interface matlab_set_field
         module procedure matlab_set_field_scalar_int4
         module procedure matlab_set_field_scalar_int8
         module procedure matlab_set_field_scalar_real4
         module procedure matlab_set_field_scalar_real8
         module procedure matlab_set_field_scalar_complex8
         module procedure matlab_set_field_scalar_complex16
         module procedure matlab_set_field_character
         module procedure matlab_set_field_logical
      end interface matlab_set_field

      interface matlab_get_field
         module procedure matlab_get_field_scalar_int4
         module procedure matlab_get_field_scalar_int8
         module procedure matlab_get_field_scalar_real4
         module procedure matlab_get_field_scalar_real8
         module procedure matlab_get_field_scalar_complex8
         module procedure matlab_get_field_scalar_complex16
         module procedure matlab_get_field_character
      end interface matlab_get_field

      integer*4, parameter :: mx_class_int4      = 12,   &
                              mx_class_int8      = 14,   &
                              mx_class_real4     = 7,    &
                              mx_class_real8     = 6

      interface matlab_copy_from_ptr
         subroutine mxCopyPtrToInteger4(px, y, n)
            mwPointer :: px
            mwSize :: n
            integer*4 :: y(n)
         end subroutine mxCopyPtrToInteger4
         subroutine mxCopyPtrToInteger8(px, y, n)
            mwPointer :: px
            mwSize :: n
            integer*8 :: y(n)
         end subroutine mxCopyPtrToInteger8

         subroutine mxCopyPtrToReal4(px, y, n)
            mwPointer :: px
            mwSize :: n
            real*4 :: y(n)
         end subroutine mxCopyPtrToReal4
         subroutine mxCopyPtrToReal8(px, y, n)
            mwPointer :: px
            mwSize :: n
            real*8 :: y(n)
         end subroutine mxCopyPtrToReal8

         subroutine mxCopyPtrToComplex8(pr, pi, y, n)
            mwPointer :: pr, pi
            mwSize :: n
            complex*8 :: y(n)
         end subroutine mxCopyPtrToComplex8
         subroutine mxCopyPtrToComplex16(pr, pi, y, n)
            mwPointer :: pr, pi
            mwSize :: n
            complex*16 :: y(n)
         end subroutine mxCopyPtrToComplex16
      end interface matlab_copy_from_ptr

      interface MATLAB_copy_to_ptr
         subroutine mxCopyReal4ToPtr(y, px, n)
            mwSize :: n
            real*4 :: y(n)
            mwPointer :: px
         end subroutine mxCopyReal4ToPtr
         subroutine mxCopyReal8ToPtr(y, px, n)
            mwSize :: n
            real*8 :: y(n)
            mwPointer :: px
         end subroutine mxCopyReal8ToPtr

         subroutine mxCopyInteger1ToPtr(y, px, n)
            mwSize :: n
            integer*1 :: y(n)
            mwPointer :: px
         end subroutine mxCopyInteger1ToPtr
         subroutine mxCopyInteger2ToPtr(y, px, n)
            mwSize :: n
            integer*2 :: y(n)
            mwPointer :: px
         end subroutine mxCopyInteger2ToPtr
         subroutine mxCopyInteger4ToPtr(y, px, n)
            mwSize :: n
            integer*4 :: y(n)
            mwPointer :: px
         end subroutine mxCopyInteger4ToPtr
         subroutine mxCopyInteger8ToPtr(y, px, n)
            mwSize :: n
            integer*8 :: y(n)
            mwPointer :: px
         end subroutine mxCopyInteger8ToPtr

         subroutine mxCopyComplex8ToPtr( Y, pr, pi, n )
            mwSize :: n
            complex*8 :: y(n)
            mwPointer :: pr, pi
         end subroutine mxCopyComplex8ToPtr
         subroutine mxCopyComplex16ToPtr( Y, pr, pi, n )
            mwSize :: n
            complex*16 :: y(n)
            mwPointer :: pr, pi
         end subroutine mxCopyComplex16ToPtr
      end interface MATLAB_copy_to_ptr

      interface MATLAB_get_field_expert
         mwPointer function mxGetField(pm, index, fieldname)
            mwPointer:: pm
            mwIndex :: index
            character*(*) :: fieldname
         end function mxGetField
      end interface MATLAB_get_field_expert


      interface MATLAB_get_field_name_by_no
         character(80) function mxGetFieldNameByNumber(pm, fieldnumber)
            mwPointer :: pm
            integer*4 :: fieldnumber
         end function mxGetFieldNameByNumber
      end interface MATLAB_get_field_name_by_no

      interface MATLAB_warning
         subroutine mexWarnMsgTxt(warningmsg)
            character*(*) :: warningmsg
         end subroutine mexWarnMsgTxt
      end interface MATLAB_warning

      interface MATLAB_get_jc
         mwPointer function mxGetJc(pm)
            mwPointer :: pm
         end function mxGetJc
      end interface MATLAB_get_jc

      interface MATLAB_get_ir
         mwPointer function mxGetIr(pm)
            mwPointer :: pm
         end function mxGetIr
      end interface MATLAB_get_ir

      interface MATLAB_get_n
         mwPointer function mxGetN(pm)
            mwPointer :: pm
         end function mxGetN
      end interface MATLAB_get_n

      interface MATLAB_get_m
         mwPointer function mxGetM(pm)
            mwPointer :: pm
         end function mxGetM
      end interface MATLAB_get_m

      interface MATLAB_get_no_fields
         integer*4 function mxGetNumberOfFields(pm)
            mwPointer :: pm
         end function mxGetNumberOfFields
      end interface MATLAB_get_no_fields

      interface MATLAB_get_nzmax
         mwSize function mxGetNzmax(pm)
            mwPointer :: pm
         end function mxGetNzmax
      end interface MATLAB_get_nzmax

      interface MATLAB_get_ptr
         mwPointer function mxGetPr(pm)
            mwPointer :: pm
         end function mxGetPr
      end interface MATLAB_get_ptr

      interface MATLAB_get_imag_ptr
         mwPointer function mxGetPi(pm)
            mwPointer :: pm
         end function mxGetPi
      end interface MATLAB_get_imag_ptr

      interface MATLAB_create_structure_expert
         mwPointer function mxCreateStructMatrix(m, n, nfields, fieldnames)
            mwSize :: m, n
            integer*4 :: nfields
            character*(*) :: fieldnames(nfields)
         end function mxCreateStructMatrix
      end interface MATLAB_create_structure_expert

      interface MATLAB_create_sparse
         mwPointer function mxCreateSparse(m, n, nzmax, complexFlag)
            mwSize :: m, n, nzmax
            integer*4 :: complexFlag
         end function mxCreateSparse
      end interface MATLAB_create_sparse

      interface
         integer*4 function mxGetString(pm, str, strlen)
            mwPointer :: pm
            character*(*) :: str
            mwSize :: strlen
         end function mxGetString

         integer*4 function mxIsStruct(pm)
            mwPointer :: pm
         end function mxIsStruct

         integer*4 function mxIsNumeric(pm)
            mwPointer :: pm
         end function mxIsNumeric

         integer*4 function mxIsSparse(pm)
            mwPointer :: pm
         end function mxIsSparse

         integer*4 function mxIsChar(pm)
            mwPointer :: pm
         end function mxIsChar

         integer*4 function mxIsComplex(pm)
            mwPointer :: pm
         end function mxIsComplex

         subroutine mxSetField(pm, index, fieldname, pvalue)
            mwPointer :: pm, pvalue
            mwIndex :: index
            character*(*) :: fieldname
         end subroutine mxSetField

         mwPointer function mxCreateCharMatrixFromStrings(m, str)
            mwSize :: m
            character*(*) :: str
         end function mxCreateCharMatrixFromStrings

         real*8 function mxGetScalar(pm)
            mwPointer :: pm
         end function mxGetScalar

         mwPointer function mxCreateString(str)
            character*(*) :: str
         end function mxCreateString

         mwPointer function mxCreateNumericMatrix(m, n, classid, complexFlag)
            mwSize :: m, n
            integer*4 :: classid, complexFlag
         end function mxCreateNumericMatrix

         mwPointer function mxCreateNumericArray(ndim, dims, classid, &
               complexFlag)
            mwSize :: ndim
            mwSize :: dims(ndim)
            integer*4 :: classid, complexFlag
         end function mxCreateNumericArray

         integer*4 function mxClassIDFromClassName(classname)
            character*(*) :: classname
         end function mxClassIDFromClassName

         integer*4 function mxGetClassID(pm)
            mwPointer :: pm
         end function mxGetClassID
            
      end interface

CONTAINS

!  -*-*-*-  M A T L A B _ g e t _ s t r i n g  -*-*-*-*-
!  -----------------------------------------------------
!  Copy string from Matlab pointer to Fortran string 
!
!  Arguments
!
!  pr - pointer to a string 
!  str - starting location into which the string should be written
!  strlen - maximum number of characters to read into str
!  ------------------------------------------------------
      SUBROUTINE MATLAB_get_string( pr, str, strlen)
         mwPointer :: pr
         CHARACTER ( len = * ) :: str
         mwSize :: strlen
         mwSize :: i

         i = mxGetString( pr, str, strlen)
      END SUBROUTINE MATLAB_get_string

!  -*-*-*-  M A T L A B _ i s _ s t r u c t u r e  -*-*-*-*-
!  -----------------------------------------------------
!  Check whether Matlab pointer is associated with a NUMERIC object
!
!  Arguments
!
!  pr - pointer to object
!  ------------------------------------------------------
      logical FUNCTION MATLAB_is_structure( pr )
         mwPointer :: pr

         integer(int4_) :: rval

         rval = mxIsStruct(pr)
         MATLAB_is_structure = (rval.eq.1)
      END FUNCTION MATLAB_is_structure

!  -*-*-*-  M A T L A B _ i s _ n u m e r i c  -*-*-*-*-
!  -----------------------------------------------------
!  Check whether Matlab pointer is associated with a NUMERIC object
!
!  Arguments
!
!  pr - pointer to object
!  ------------------------------------------------------
      logical FUNCTION MATLAB_is_numeric( pr )
         mwPointer :: pr

         integer(int4_) :: rval

         rval = mxIsNumeric(pr)
         MATLAB_is_numeric = (rval.eq.1)
      END FUNCTION MATLAB_is_numeric

!  -*-*-*-  M A T L A B _ i s _ s p a r s e  -*-*-*-*-
!  -----------------------------------------------------
!  Check whether Matlab pointer is associated with a SPARSE object
!
!  Arguments
!
!  pr - pointer to object
!  ------------------------------------------------------
      logical FUNCTION MATLAB_is_sparse( pr )
         mwPointer :: pr

         integer(int4_) :: rval

         rval = mxIsSparse(pr)
         MATLAB_is_sparse = (rval.eq.1)
      END FUNCTION MATLAB_is_sparse

!  -*-*-*-  M A T L A B _ i s _ c h a r  -*-*-*-*-
!  -----------------------------------------------------
!  Check whether Matlab pointer is associated with a CHARACTER object
!
!  Arguments
!
!  pr - pointer to object
!  ------------------------------------------------------
      logical FUNCTION MATLAB_is_character( pr )
         mwPointer :: pr

         integer(int4_) :: rval

         rval = mxIsChar(pr)
         MATLAB_is_character = (rval .eq. 1)
      END FUNCTION MATLAB_is_character

!  -*-*-*-  M A T L A B _ c r e a t e _ s u b s t r u c t u r e  -*-*-*-*-
!  -----------------------------------------------------
!  Create a named sub-structure of a structure
!
!  Arguments
!
!  struct - existing pointer to the structure
!  name - name of component of the structure
!  pr - pointer to the sub-tructure
!  ninform - number of components of the substructure
!  finform - names of components of the substructure
!  ------------------------------------------------------
      SUBROUTINE MATLAB_create_substructure(struct, name, pr, ninform, finform)
         mwPointer ::struct, pr
         CHARACTER ( len = * ) :: name
         integer(int4_) :: ninform
         CHARACTER ( LEN = * ), DIMENSION( ninform ) :: finform

         pr = mxCreateStructMatrix( 1_mws_, 1_mws_, ninform, finform )
         CALL mxSetField( struct, 1_mwi_, name, pr )

      END SUBROUTINE MATLAB_create_substructure

!  -*-*-*- h s l m x  C r e a t e  I n t e g e r  C o m p o n e n t -*-*-*-*-
!  -----------------------------------------------
!  Create a named INTEGER component of a structure
!
!  Arguments
!
!  struct - structure
!  name - name of component
!  pr - pointer to the structure
!  -----------------------------------------------
      SUBROUTINE hslmxCreateIntegerComponent( struct, name, pr )
         mwPointer ::struct
         CHARACTER ( len = * ) :: name
         mwPointer :: pr

         pr = MATLAB_create_integer( ) 
         CALL mxSetField( struct, 1_mwi_, name, pr )
      END SUBROUTINE hslmxCreateIntegerComponent

!  -*- h s l m x  C r e a t e  I n t e g e r  A r r a y  C o m p o n e n t -*-
!  -----------------------------------------------------
!  Create a named INTEGER array component of a structure
!
!  Arguments
!
!  struct - structure
!  name - name of component
!  n - dimension of array
!  pr - pointer to the structure
!  -----------------------------------------------------
      SUBROUTINE hslmxCreateIntegerArrayComponent( struct, name, n, pr )
         mwPointer ::struct
         CHARACTER ( len = * ) :: name
         mwSize :: n
         mwPointer :: pr

         pr = MATLAB_create_integer( m=n ) 
         CALL mxSetField( struct, 1_mwi_, name, pr )
      END SUBROUTINE hslmxCreateIntegerArrayComponent

!  -*- h s l m x  C r e a t e  I n t e g e r  M a t r i x  C o m p o n e n t -*-
!  -----------------------------------------------------
!  Create a named INTEGER matrix component of a structure
!
!  Arguments
!
!  struct - structure
!  name - name of component
!  m - row dimension of matrix
!  n - column dimension of matrix
!  pr - pointer to the structure
!  -----------------------------------------------------
      SUBROUTINE hslmxCreateIntegerMatrixComponent( struct, name, m, n, pr )
         mwPointer ::struct
         CHARACTER ( len = * ) :: name
         mwSize :: m, n
         mwPointer :: pr

         pr = MATLAB_create_integer( m=m, n=n ) 
         CALL mxSetField( struct, 1_mwi_, name, pr )
      END SUBROUTINE hslmxCreateIntegerMatrixComponent

!  -*-*-*-*-*- h s l m x  C r e a t e  R e a l  C o m p o n e n t -*-*-*-*-*-*-
!  --------------------------------------------
!  Create a named REAL component of a structure
!
!  Arguments
!
!  struct - structure
!  name - name of component
!  pr - pointer to the structure
!  --------------------------------------------
      SUBROUTINE hslmxCreateRealComponent( struct, name, pr )
         mwPointer ::struct
         CHARACTER ( len = * ) :: name
         mwPointer :: pr

         pr = MATLAB_create_real( ) 
         CALL mxSetField( struct, 1_mwi_, name, pr )
      END SUBROUTINE hslmxCreateRealComponent


!  -*-*-*- h s l m x  C r e a t e  R e a l  A r r a y  C o m p o n e n t -*-*-
!  --------------------------------------------------
!  Create a named REAL array component of a structure
!
!  Arguments
!
!  struct - structure
!  name - name of component
!  n - dimension of array
!  pr - pointer to the structure
!  --------------------------------------------------
      SUBROUTINE hslmxCreateRealArrayComponent( struct, name, n, pr )
         mwPointer ::struct
         CHARACTER ( len = * ) :: name
         mwSize :: n
         mwPointer :: pr

         pr = MATLAB_create_real( m=n ) 
         CALL mxSetField( struct, 1_mwi_, name, pr )
      END SUBROUTINE hslmxCreateRealArrayComponent


!  -*-*-*- h s l m x  C r e a t e  R e a l  M a t r i x  C o m p o n e n t -*-*-
!  --------------------------------------------------
!  Create a named REAL matrix component of a structure
!
!  Arguments
!
!  struct - structure
!  name - name of component
!  m - row dimension of matrix
!  n - column dimension of matrix
!  pr - pointer to the structure
!  --------------------------------------------------
      SUBROUTINE hslmxCreateRealMatrixComponent( struct, name, m, n, pr )
         mwPointer ::struct
         CHARACTER ( len = * ) :: name
         mwSize :: m, n
         mwPointer :: pr

         pr = MATLAB_create_real( m=m, n=n ) 
         CALL mxSetField( struct, 1_mwi_, name, pr )
      END SUBROUTINE hslmxCreateRealMatrixComponent

!  -*-*-*-  M A T L A B _ c r e a t e _ c h a r _ c o m p o n e n t  -*-*-*-
!  -------------------------------------------------
!  Create a named CHARACTER component of a structure
!
!  Arguments
!
!  struct - structure
!  name - name of component
!  pr - pointer to the structure
!  -------------------------------------------------
      SUBROUTINE MATLAB_create_character_component( struct, name, pr )
         mwPointer ::struct
         CHARACTER ( len = * ) :: name
         mwPointer :: pr

         INTEGER, PARAMETER :: len_blank = 80
         CHARACTER ( len = len_blank ) :: blank

         blank = REPEAT( ' ', len_blank )

         pr = mxCreateCharMatrixFromStrings( 1_mws_, blank )
         CALL mxSetField( struct, 1_mwi_, name, pr )
      END SUBROUTINE MATLAB_create_character_component

!  -*-*-  M A T L A B _ c r e a t e _ l o g i c a l _ c o m p o n e n t  -*-*-
!  -----------------------------------------------------
!  Create a named LOGICAL component of a structure
!
!  Arguments
!
!  struct - structure
!  name - name of component
!  pr - pointer to the structure
!
! ** - NB - ** This is a bodge since Mex doesn't appear 
!              to handle Fortran logicals ** - NB - **
!  ----------------------------------------------------
      SUBROUTINE MATLAB_create_logical_component( struct, name, pr )
         mwPointer ::struct
         CHARACTER ( len = * ) :: name
         mwPointer :: pr

         pr = MATLAB_create_integer( )
         CALL mxSetField( struct, 1_mwi_, name, pr )
      END SUBROUTINE MATLAB_create_logical_component


!  -*-*-*-*-*-*-*-*-*-*- h s l m x  G e t  I n t e g e r  -*-*-*-*-*-*-*-*-*-*-
!  ---------------------------------------------------------
!  Obtain an integer value from the component of a structure
!
!  Arguments
!
!  ps - given pointer to the structure
!  name - name of component
!  pc - pointer to the component
!  value - value of the component
!  ---------------------------------------------------------
      SUBROUTINE hslmxGetInteger( ps, name, pc, value )
         mwPointer :: ps, pc
         CHARACTER ( LEN = * ) :: name
         INTEGER :: value

         pc = mxGetField( ps, 1_mwi_, name )
         value = INT( mxGetScalar( pc ) )
      END SUBROUTINE hslmxGetInteger


!  -*-*-*-*-*-*-*-*-*-*- h s l m x  G e t  R e a l  -*-*-*-*-*-*-*-*-*-*-*-
!  -----------------------------------------------------
!  Obtain a real value from the component of a structure
!
!  Arguments
!
!  ps - given pointer to the structure
!  name - name of component
!  pc - pointer to the component
!  value - value of the component
!  -----------------------------------------------------
      SUBROUTINE hslmxGetReal( ps, name, pc, value )
         mwPointer :: ps, pc
         CHARACTER ( LEN = * ) :: name
         real(real8_) :: value

         pc = mxGetField( ps, 1_mwi_, name )
         value = mxGetScalar( pc )
      END SUBROUTINE hslmxGetReal


!  -*-*-*-*-*-*-*-*-*-*- h s l m x  G e t  L o g i c a l  -*-*-*-*-*-*-*-*-*-*-
!  ---------------------------------------------------------
!  Obtain a logical value from the component of a structure
!
!  Arguments
!
!  ps - given pointer to the structure
!  name - name of component
!  pc - pointer to the component
!  value - value of the component
!
! ** - NB - ** This is a bodge since Mex doesn't appear 
!              to handle Fortran logicals ** - NB - **
!  ---------------------------------------------------------
      SUBROUTINE hslmxGetLogical( ps, name, pc, value )
         mwPointer :: ps, pc
         CHARACTER ( LEN = * ) :: name
         LOGICAL :: value

         pc = mxGetField( ps, 1_mwi_, name )
         value = ( 1 .eq. INT( mxGetScalar( pc ) ) )
      END SUBROUTINE hslmxGetLogical


!  -*-*-*-*-*-*-*-*-*- h s l m x  G e t  C h a r a c t e r  -*-*-*-*-*-*-*-*-*-
!  -----------------------------------------------------------------
!  Obtain a character string value from the component of a structure
!
!  Arguments
!
!  ps - given pointer to the structure
!  name - name of component
!  pc - pointer to the component
!  value - value of the character string
!  len - maximum length of the string
!  ------------------------------------------------------------------
      SUBROUTINE hslmxGetCharacter( ps, name, pc, value, len )
         mwPointer :: ps, pc
         CHARACTER ( LEN = * ) :: name
         CHARACTER ( LEN = * ) :: value
         mwSize :: len

         integer(int4_) :: flag

         pc = mxGetField( ps, 1_mwi_, name )
         flag = mxGetString( pc, value, len )
      END SUBROUTINE hslmxGetCharacter

!  -*-*-*-  h s l m x  C o p y  L o g i c a l  T o  P t r  -*-*-*-*-
!  --------------------------------------------------------
!  Copy LOGICAL value from Fortran scalar to Matlab pointer
!
!  Arguments
!
!  Y   - Logical Fortran scalar
!  px  - Pointer to ir or jc array
!  --------------------------------------------------------
      SUBROUTINE hslmxCopyLogicalToPtr( Y, px )
         mwPointer :: px
         LOGICAL :: Y

         CALL hslmxCopyLogicalArrayToPtr( (/y/), px, 1_mws_ )
      END SUBROUTINE hslmxCopyLogicalToPtr


!  -*-*-*-  h s l m x  C o p y  L o g i c a l  A r r a y  T o  P t r  -*-*-*-*-
!  --------------------------------------------------------------
!  Copy LOGICAL values from Fortran array to Matlab pointer array
!
!  Arguments
!
!  Y   - Logical Fortran array
!  px  - Pointer to ir or jc array
!  n   - number of elements to copy
!  --------------------------------------------------------------
      SUBROUTINE hslmxCopyLogicalArrayToPtr( Y, px, n )
         mwPointer :: px
         mwSize :: n
         LOGICAL, DIMENSION( n ) :: Y

         INTEGER :: i
         integer, dimension(n) :: ly ! Deliberatly left as default integer

         DO i = 1, n
           IF ( Y( i ) ) THEN
             LY( i ) = 1
           ELSE
             LY( i ) = 0
           END IF
         END DO
         CALL MATLAB_copy_to_ptr( LY, px, n )
      END SUBROUTINE hslmxCopyLogicalArrayToPtr


!  -*-*-*-  h s l m x  S e t  C h a r a c t e r  C o m p o n e n t  -*-*-*-*-
!  ---------------------------------------------------------
!  Set a named CHARACTER component of a structure to a value
!
!  Arguments
!
!  struct - structure
!  name - name of component
!  value - character value to be assigned
!  ---------------------------------------------------------
      SUBROUTINE hslmxSetCharacterComponent( struct, name, value )
         mwPointer ::struct
         CHARACTER ( len = * ) :: name
         CHARACTER ( len = * ) :: value

         CALL mxSetField( struct, 1_mwi_, name, mxCreateString( value ) )
      END SUBROUTINE hslmxSetCharacterComponent

!  -*-*-*-  h s l m x  C r e a t e  I n t e g e r  M a t r i x -*-*-*-*-
!  --------------------------------------------
!  Create an unpopulated INTEGER pointer array
!
!  Arguments
!
!  m - number of rows
!  n - number of columns
!  --------------------------------------------

      mwPointer FUNCTION MATLAB_create_integer( m, n, kind )
         mwSize, optional :: m, n
         integer, optional :: kind

         INTEGER :: mykind
         mwSize :: mym, myn
         integer(int4_) :: classID

         mykind = di_
         if(present(kind)) mykind = kind
         mym = 1
         if(present(m)) mym = m
         myn = 1
         if(present(n)) myn = n

         SELECT CASE ( mykind )
         CASE ( int8_ )
            classID = mxClassIDFromClassName('int64')
         CASE ( int2_ )
            classID = mxClassIDFromClassName('int16')
         CASE ( int1_ )
            classID = mxClassIDFromClassName('int8')
         CASE default
            classID = mxClassIDFromClassName('int32')
         END SELECT

         if(mym.eq.1 .and. myn.eq.1) then
            MATLAB_create_integer = mxCreateNumericMatrix( 1_mws_, 1_mws_, &
               classID, 0_int4_ ) 
         else
            MATLAB_create_integer = mxCreateNumericArray( 2_mws_, &
               (/ mym, myn /), classID, 0_int4_ ) 
         endif
      END FUNCTION MATLAB_create_integer

!  -*-*-*-  h s l m x  C r e a t e  R e a l  M a t r i x  -*-*-*-*-
!  -----------------------------------------
!  Create an unpopulated REAL pointer matrix
!
!  Arguments
!
!  m - number of rows
!  n - number of columns
!  -----------------------------------------
      mwPointer FUNCTION MATLAB_create_real( m, n, kind )
         mwSize, optional :: m, n
         integer, optional :: kind

         integer :: mykind
         mwSize :: mym, myn
         integer(int4_) :: classID

         mykind = real8_
         if(present(kind)) mykind = kind
         mym = 1
         if(present(m)) mym = m
         myn = 1
         if(present(n)) myn = n

         classID = mxClassIDFromClassName('single')
         if(mykind .eq. real8_) classID = mxClassIDFromClassName('double')

         if(mym.eq.1 .and. myn.eq.1) then
            MATLAB_create_real = mxCreateNumericMatrix( 1_mws_, 1_mws_, &
               classID, 0_int4_ ) 
         else
            MATLAB_create_real = mxCreateNumericArray( 2_mws_, &
               (/ mym, myn /), classID, 0_int4_ ) 
         endif
      END FUNCTION MATLAB_create_real

!  -*-*-*-  h s l m x  C r e a t e  C o m p l e x   -*-*-*-*-
!  ----------------------------------
!  Create an unpopulated COMPLEX pointer
!
!  Arguments:
!
!  kind  optional argument specifying kind of complex to create. If not
!        specified a default real is created
!  ----------------------------------
      mwPointer FUNCTION MATLAB_create_complex( m, n, kind )
         mwSize, optional, intent(in) :: m
         mwSize, optional, intent(in) :: n
         integer, optional, intent(in) :: kind

         integer :: mykind
         mwSize :: mym, myn
         integer(int4_) :: classID

         mykind = real8_
         if(present(kind)) mykind = kind
         mym = 1
         if(present(m)) mym = m
         myn = 1
         if(present(n)) myn = n

         classID = mxClassIDFromClassName('single')
         if(mykind .eq. real8_) classID = mxClassIDFromClassName('double')

         if(mym.eq.1 .and. myn.eq.1) then
            MATLAB_create_complex = mxCreateNumericMatrix( 1_mws_, 1_mws_, &
               classID, 1_int4_ ) 
         else
            MATLAB_create_complex = mxCreateNumericArray( 2_mws_, &
               (/ mym, myn /), classID, 1_int4_ ) 
         endif
      END FUNCTION MATLAB_create_complex

!  -*-*-*-  M A T L A B _ c r e a t e _ sparse -*-*-*-*-
!  -----------------------------------------------------
!  Create a Matlab sparse matrix
!
!  Arguments
!
!  m - the desired number of rows. This must be a positive integer.
!  n - the desired number of columns. This must be a positive integer.
!  nzmax - The number of elements that MATLAB_create_sparse should allocate
!          to hold the pr, ir, and, if ComplexFlag is 1 in , pi arrays.
!          Set the  value of nzmax to be greater than or equal to the number of
!          nonzero elements you plan to put into the mxArray, but make sure
!          that nzmax is less than or equal to m*n.the desired number of fields
!          in each element.
! ComplexFlag - If the mxArray you are creating is to contain imaginary data,
!               set ComplexFlag to 1 . Otherwise, set ComplexFlag to 0.
!
!  ------------------------------------------------------

      logical function MATLAB_is_complex(pm)
         mwPointer :: pm

         integer(int4_) :: temp

         temp = mxIsComplex(pm)
         MATLAB_is_complex = (temp .eq. 1)
      end function MATLAB_is_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine matlab_to_fortran_scalar_int4(mp, scalar, argname)
         mwPointer :: mp
         integer(int4_), intent(out) :: scalar
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwSize :: m, n
         mwPointer :: dataptr
         integer(int4_) :: tempint4(1)
         integer(int8_) :: tempint8(1)
         real(real4_) :: tempreal4(1)
         real(real8_) :: tempreal8(1)

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)
         if(m.ne.1 .or. n.ne.1) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected 1x1 matrix."
            call MATLAB_context_error("hsl_matlab:Expects1x1", errmsg)
         endif

         if(MATLAB_is_complex(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected integer matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         endif

         dataptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_int4)
            call MATLAB_copy_from_ptr(dataptr, tempint4, 1_mws_)
            scalar = tempint4(1)
         case(mx_class_int8)
            call MATLAB_copy_from_ptr(dataptr, tempint8, 1_mws_)
            scalar = tempint8(1)
         case(mx_class_real4)
            call MATLAB_copy_from_ptr(dataptr, tempreal4, 1_mws_)
            scalar = int(tempreal4(1))
            if(real(scalar) .ne. tempreal4(1)) then
               write(errmsg, "(3a)") &
                  "Error in argument ", argname, ". Expected integer data."
               call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
            endif
         case(mx_class_real8)
            call MATLAB_copy_from_ptr(dataptr, tempreal8, 1_mws_)
            scalar = int(tempreal8(1))
            if(real(scalar) .ne. tempreal8(1)) then
               write(errmsg, "(3a)") &
                  "Error in argument ", argname, ". Expected integer data."
               call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
            endif
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected integer data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select
      end subroutine matlab_to_fortran_scalar_int4
      subroutine matlab_to_fortran_scalar_int8(mp, scalar, argname)
         mwPointer :: mp
         integer(int8_), intent(out) :: scalar
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwSize :: m, n
         mwPointer :: dataptr
         integer(int4_) :: tempint4(1)
         integer(int8_) :: tempint8(1)
         real(real4_) :: tempreal4(1)
         real(real8_) :: tempreal8(1)

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)
         if(m.ne.1 .or. n.ne.1) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected 1x1 matrix."
            call MATLAB_context_error("hsl_matlab:Expects1x1", errmsg)
         endif

         if(MATLAB_is_complex(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected integer matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         endif

         dataptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_int4)
            call MATLAB_copy_from_ptr(dataptr, tempint4, 1_mws_)
            scalar = tempint4(1)
         case(mx_class_int8)
            call MATLAB_copy_from_ptr(dataptr, tempint8, 1_mws_)
            scalar = tempint8(1)
         case(mx_class_real4)
            call MATLAB_copy_from_ptr(dataptr, tempreal4, 1_mws_)
            scalar = int(tempreal4(1))
            if(real(scalar) .ne. tempreal4(1)) then
               write(errmsg, "(3a)") &
                  "Error in argument ", argname, ". Expected integer data."
               call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
            endif
         case(mx_class_real8)
            call MATLAB_copy_from_ptr(dataptr, tempreal8, 1_mws_)
            scalar = int(tempreal8(1))
            if(real(scalar) .ne. tempreal8(1)) then
               write(errmsg, "(3a)") &
                  "Error in argument ", argname, ". Expected integer data."
               call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
            endif
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected integer data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select
      end subroutine matlab_to_fortran_scalar_int8
      subroutine matlab_to_fortran_scalar_real4(mp, scalar, argname)
         mwPointer :: mp
         real(real4_), intent(out) :: scalar
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwSize :: m, n
         mwPointer :: dataptr
         integer(int4_) :: tempint4(1)
         integer(int8_) :: tempint8(1)
         real(real4_) :: tempreal4(1)
         real(real8_) :: tempreal8(1)

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)
         if(m.ne.1 .or. n.ne.1) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected 1x1 matrix."
            call MATLAB_context_error("hsl_matlab:Expects1x1", errmsg)
         endif

         if(MATLAB_is_complex(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected real matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsReal", errmsg)
         endif

         dataptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_int4)
            call MATLAB_copy_from_ptr(dataptr, tempint4, 1_mws_)
            scalar = tempint4(1)
         case(mx_class_int8)
            call MATLAB_copy_from_ptr(dataptr, tempint8, 1_mws_)
            scalar = tempint8(1)
         case(mx_class_real4)
            call MATLAB_copy_from_ptr(dataptr, tempreal4, 1_mws_)
            scalar = tempreal4(1)
         case(mx_class_real8)
            call MATLAB_copy_from_ptr(dataptr, tempreal8, 1_mws_)
            scalar = tempreal8(1)
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected real data."
            call MATLAB_context_error("hsl_matlab:ExpectsReal", errmsg)
         end select
      end subroutine matlab_to_fortran_scalar_real4
      subroutine matlab_to_fortran_scalar_real8(mp, scalar, argname)
         mwPointer :: mp
         real(real8_), intent(out) :: scalar
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwSize :: m, n
         mwPointer :: dataptr
         integer(int4_) :: tempint4(1)
         integer(int8_) :: tempint8(1)
         real(real4_) :: tempreal4(1)
         real(real8_) :: tempreal8(1)

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)
         if(m.ne.1 .or. n.ne.1) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected 1x1 matrix."
            call MATLAB_context_error("hsl_matlab:Expects1x1", errmsg)
         endif

         if(MATLAB_is_complex(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected real matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsReal", errmsg)
         endif

         dataptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_int4)
            call MATLAB_copy_from_ptr(dataptr, tempint4, 1_mws_)
            scalar = tempint4(1)
         case(mx_class_int8)
            call MATLAB_copy_from_ptr(dataptr, tempint8, 1_mws_)
            scalar = tempint8(1)
         case(mx_class_real4)
            call MATLAB_copy_from_ptr(dataptr, tempreal4, 1_mws_)
            scalar = tempreal4(1)
         case(mx_class_real8)
            call MATLAB_copy_from_ptr(dataptr, tempreal8, 1_mws_)
            scalar = tempreal8(1)
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected real data."
            call MATLAB_context_error("hsl_matlab:ExpectsReal", errmsg)
         end select
      end subroutine matlab_to_fortran_scalar_real8
      subroutine matlab_to_fortran_scalar_complex8(mp, scalar, argname)
         mwPointer :: mp
         complex(real4_), intent(out) :: scalar
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwSize :: m, n
         mwPointer :: realptr, imagptr
         real(real4_) :: tempreal4(1)
         real(real8_) :: tempreal8(1)
         complex(real4_) :: tempcomplex8(1)
         complex(real8_) :: tempcomplex16(1)

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)
         if(m.ne.1 .or. n.ne.1) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected 1x1 matrix."
            call MATLAB_context_error("hsl_matlab:Expects1x1", errmsg)
         endif

         realptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_real4)
            if(MATLAB_is_complex(mp)) then
               imagptr = MATLAB_get_imag_ptr(mp)
               call MATLAB_copy_from_ptr(realptr, imagptr, tempcomplex8, 1_mws_)
               scalar = tempcomplex8(1)
            else
               call MATLAB_copy_from_ptr(realptr, tempreal4, 1_mws_)
               scalar = tempreal4(1)
            endif
         case(mx_class_real8)
            if(MATLAB_is_complex(mp)) then
               imagptr = MATLAB_get_imag_ptr(mp)
               call MATLAB_copy_from_ptr(realptr, imagptr, tempcomplex16, &
                  1_mws_)
               scalar = tempcomplex16(1)
            else
               call MATLAB_copy_from_ptr(realptr, tempreal8, 1_mws_)
               scalar = tempreal8(1)
            endif
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected complex matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsComplex", errmsg)
         end select
      end subroutine matlab_to_fortran_scalar_complex8
      subroutine matlab_to_fortran_scalar_complex16(mp, scalar, argname)
         mwPointer :: mp
         complex(real8_), intent(out) :: scalar
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwSize :: m, n
         mwPointer :: realptr, imagptr
         real(real4_) :: tempreal4(1)
         real(real8_) :: tempreal8(1)
         complex(real4_) :: tempcomplex8(1)
         complex(real8_) :: tempcomplex16(1)

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)
         if(m.ne.1 .or. n.ne.1) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected 1x1 matrix."
            call MATLAB_context_error("hsl_matlab:Expects1x1", errmsg)
         endif

         realptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_real4)
            if(MATLAB_is_complex(mp)) then
               imagptr = MATLAB_get_imag_ptr(mp)
               call MATLAB_copy_from_ptr(realptr, imagptr, tempcomplex8, 1_mws_)
               scalar = tempcomplex8(1)
            else
               call MATLAB_copy_from_ptr(realptr, tempreal4, 1_mws_)
               scalar = tempreal4(1)
            endif
         case(mx_class_real8)
            if(MATLAB_is_complex(mp)) then
               imagptr = MATLAB_get_imag_ptr(mp)
               call MATLAB_copy_from_ptr(realptr, imagptr, tempcomplex16, &
                  1_mws_)
               scalar = tempcomplex16(1)
            else
               call MATLAB_copy_from_ptr(realptr, tempreal8, 1_mws_)
               scalar = tempreal8(1)
            endif
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected complex matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsComplex", errmsg)
         end select
      end subroutine matlab_to_fortran_scalar_complex16

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine matlab_to_fortran_array1d_int4(mp, array1d, n, argname)
         mwPointer :: mp
         integer(int4_), dimension(:), allocatable, intent(out) :: array1d
         mwSize :: n
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwSize :: m
         mwPointer :: dataptr
         integer(int8_), dimension(:), allocatable :: tempint8
         real(real4_), dimension(:), allocatable :: tempreal4
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         integer :: i

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)
         if(m.ne.1 .and. n.ne.1) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected 1-dimensional matrix."
            call MATLAB_context_error("hsl_matlab:Expects1D", errmsg)
         endif
         if(n.eq.1) n = m

         if(MATLAB_is_complex(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected integer matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         endif

         if(MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected dense matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsDense", errmsg)
         endif

         deallocate(array1d, stat=st)
         allocate(array1d(n), stat=st)
         if(st.ne.0) goto 100

         dataptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_int4)
            call MATLAB_copy_from_ptr(dataptr, array1d, n)
         case(mx_class_int8)
            allocate(tempint8(n), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(dataptr, tempint8, n)
            array1d(:) = tempint8(:)
         case(mx_class_real4)
            allocate(tempreal4(n), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(dataptr, tempreal4, n)
            do i = 1, n
               array1d(i) = int(tempreal4(i))
               if(real(array1d(i)) .ne. tempreal4(i)) then
                  write(errmsg, "(3a)") &
                     "Error in argument ", argname, &
                     ". Expected integer data."
                  call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
               endif
            end do
         case(mx_class_real8)
            allocate(tempreal8(n), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(dataptr, tempreal8, n)
            do i = 1, n
               array1d(i) = int(tempreal8(i))
               if(real(array1d(i)) .ne. tempreal8(i)) then
                  write(errmsg, "(3a)") &
                     "Error in argument ", argname, &
                     ". Expected integer data."
                  call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
               endif
            end do
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected integer data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select

         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_array1d_int4
      subroutine matlab_to_fortran_array1d_int8(mp, array1d, n, argname)
         mwPointer :: mp
         integer(int8_), dimension(:), allocatable, intent(out) :: array1d
         mwSize :: n
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwSize :: m
         mwPointer :: dataptr
         integer(int4_), dimension(:), allocatable :: tempint4
         real(real4_), dimension(:), allocatable :: tempreal4
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         integer :: i

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)
         if(m.ne.1 .and. n.ne.1) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected 1-dimensional matrix."
            call MATLAB_context_error("hsl_matlab:Expects1D", errmsg)
         endif
         if(n.eq.1) n = m

         if(MATLAB_is_complex(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected integer matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         endif

         if(MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected dense matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsDense", errmsg)
         endif

         deallocate(array1d, stat=st)
         allocate(array1d(n), stat=st)
         if(st.ne.0) goto 100

         dataptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_int8)
            call MATLAB_copy_from_ptr(dataptr, array1d, n)
         case(mx_class_int4)
            allocate(tempint4(n), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(dataptr, tempint4, n)
            array1d(:) = tempint4(:)
         case(mx_class_real4)
            allocate(tempreal4(n), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(dataptr, tempreal4, n)
            do i = 1, n
               array1d(i) = int(tempreal4(i))
               if(real(array1d(i)) .ne. tempreal4(i)) then
                  write(errmsg, "(3a)") &
                     "Error in argument ", argname, &
                     ". Expected integer data."
                  call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
               endif
            end do
         case(mx_class_real8)
            allocate(tempreal8(n), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(dataptr, tempreal8, n)
            do i = 1, n
               array1d(i) = int(tempreal8(i))
               if(real(array1d(i)) .ne. tempreal8(i)) then
                  write(errmsg, "(3a)") &
                     "Error in argument ", argname, &
                     ". Expected integer data."
                  call MATLAB_context_error("hsl_matlabLExpectsInteger", errmsg)
               endif
            end do
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected integer data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select

         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)

      end subroutine matlab_to_fortran_array1d_int8
      subroutine matlab_to_fortran_array1d_real4(mp, array1d, n, argname)
         mwPointer :: mp
         real(real4_), dimension(:), allocatable, intent(out) :: array1d
         mwSize :: n
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwSize :: m
         mwPointer :: dataptr
         integer(int4_), dimension(:), allocatable :: tempint4
         integer(int8_), dimension(:), allocatable :: tempint8
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)
         if(m.ne.1 .and. n.ne.1) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected 1-dimensional matrix."
            call MATLAB_context_error("hsl_matlab:Expects1D", errmsg)
         endif
         if(n.eq.1) n = m

         if(MATLAB_is_complex(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected real matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsReal", errmsg)
         endif

         if(MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected dense matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsDense", errmsg)
         endif

         deallocate(array1d, stat=st)
         allocate(array1d(n), stat=st)
         if(st.ne.0) goto 100

         dataptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_int4)
            allocate(tempint4(n), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(dataptr, tempint4, n)
            array1d(:) = tempint4(:)
         case(mx_class_int8)
            allocate(tempint8(n), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(dataptr, tempint8, n)
            array1d(:) = tempint8(:)
         case(mx_class_real4)
            call MATLAB_copy_from_ptr(dataptr, array1d, n)
         case(mx_class_real8)
            allocate(tempreal8(n), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(dataptr, tempreal8, n)
            array1d(:) = tempreal8(:)
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected real data."
            call MATLAB_context_error("hsl_matlab:ExpectsReal", errmsg)
         end select

         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_array1d_real4
      subroutine matlab_to_fortran_array1d_real8(mp, array1d, n, argname)
         mwPointer :: mp
         real(real8_), dimension(:), allocatable, intent(out) :: array1d
         mwSize :: n
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwSize :: m
         mwPointer :: dataptr
         integer(int4_), dimension(:), allocatable :: tempint4
         integer(int8_), dimension(:), allocatable :: tempint8
         real(real4_), dimension(:), allocatable :: tempreal4
         integer :: st

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)
         if(m.ne.1 .and. n.ne.1) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected 1-dimensional matrix."
            call MATLAB_context_error("hsl_matlab:Expects1D", errmsg)
         endif
         if(n.eq.1) n = m

         if(MATLAB_is_complex(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected real matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsReal", errmsg)
         endif

         if(MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected dense matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsDense", errmsg)
         endif

         deallocate(array1d, stat=st)
         allocate(array1d(n), stat=st)
         if(st.ne.0) goto 100

         dataptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_int4)
            allocate(tempint4(n), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(dataptr, tempint4, n)
            array1d(:) = tempint4(:)
         case(mx_class_int8)
            allocate(tempint8(n), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(dataptr, tempint8, n)
            array1d(:) = tempint8(:)
         case(mx_class_real4)
            allocate(tempreal4(n), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(dataptr, tempreal4, n)
            array1d(:) = tempreal4(:)
         case(mx_class_real8)
            call MATLAB_copy_from_ptr(dataptr, array1d, n)
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected real data."
            call MATLAB_context_error("hsl_matlab:ExpectsReal", errmsg)
         end select

         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_array1d_real8
      subroutine matlab_to_fortran_array1d_complex8(mp, array1d, n, argname)
         mwPointer :: mp
         complex(real4_), dimension(:), allocatable, intent(out) :: array1d
         mwSize :: n
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwSize :: m
         mwPointer :: realptr, imagptr
         real(real4_), dimension(:), allocatable :: tempreal4
         real(real8_), dimension(:), allocatable :: tempreal8
         complex(real8_), dimension(:), allocatable :: tempcomplex16
         integer :: st

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)
         if(m.ne.1 .and. n.ne.1) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected 1-dimensional matrix."
            call MATLAB_context_error("hsl_matlab:Expects1D", errmsg)
         endif
         if(n.eq.1) n = m

         if(MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected dense matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsDense", errmsg)
         endif


         deallocate(array1d, stat=st)
         allocate(array1d(n), stat=st)
         if(st.ne.0) goto 100

         realptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_real4)
            if(MATLAB_is_complex(mp)) then
               imagptr = MATLAB_get_imag_ptr(mp)
               call MATLAB_copy_from_ptr(realptr, imagptr, array1d, n)
            else
               allocate(tempreal4(n), stat=st)
               if(st.ne.0) goto 100
               call MATLAB_copy_from_ptr(realptr, tempreal4, n)
               array1d(:) = tempreal4(:)
            endif
         case(mx_class_real8)
            if(MATLAB_is_complex(mp)) then
               imagptr = MATLAB_get_imag_ptr(mp)
               allocate(tempcomplex16(n), stat=st)
               if(st.ne.0) goto 100
               call MATLAB_copy_from_ptr(realptr, imagptr, tempcomplex16, n)
               array1d(:) = tempcomplex16(:)
            else
               allocate(tempreal8(n), stat=st)
               if(st.ne.0) goto 100
               call MATLAB_copy_from_ptr(realptr, tempreal8, n)
               array1d(:) = tempreal8(:)
            endif
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected complex matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsComplex", errmsg)
         end select

         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_array1d_complex8
      subroutine matlab_to_fortran_array1d_complex16(mp, array1d, n, argname)
         mwPointer :: mp
         complex(real8_), dimension(:), allocatable, intent(out) :: array1d
         mwSize :: n
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwSize :: m
         mwPointer :: realptr, imagptr
         real(real4_), dimension(:), allocatable :: tempreal4
         real(real8_), dimension(:), allocatable :: tempreal8
         complex(real4_), dimension(:), allocatable :: tempcomplex8
         integer :: st

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)
         if(m.ne.1 .and. n.ne.1) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected 1-dimensional matrix."
            call MATLAB_context_error("hsl_matlab:Expects1D", errmsg)
         endif
         if(n.eq.1) n = m

         if(MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected dense matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsDense", errmsg)
         endif


         deallocate(array1d, stat=st)
         allocate(array1d(n), stat=st)
         if(st.ne.0) goto 100

         realptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_real4)
            if(MATLAB_is_complex(mp)) then
               imagptr = MATLAB_get_imag_ptr(mp)
               allocate(tempcomplex8(n), stat=st)
               if(st.ne.0) goto 100
               call MATLAB_copy_from_ptr(realptr, imagptr, tempcomplex8, n)
               array1d(:) = tempcomplex8(:)
            else
               allocate(tempreal4(n), stat=st)
               if(st.ne.0) goto 100
               call MATLAB_copy_from_ptr(realptr, tempreal4, n)
               array1d(:) = tempreal4(:)
            endif
         case(mx_class_real8)
            if(MATLAB_is_complex(mp)) then
               imagptr = MATLAB_get_imag_ptr(mp)
               call MATLAB_copy_from_ptr(realptr, imagptr, array1d, n)
            else
               allocate(tempreal8(n), stat=st)
               if(st.ne.0) goto 100
               call MATLAB_copy_from_ptr(realptr, tempreal8, n)
               array1d(:) = tempreal8(:)
            endif
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected complex matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsComplex", errmsg)
         end select

         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_array1d_complex16

      subroutine matlab_to_fortran_array2d_int4(mp, array2d, m, n, argname)
         mwPointer :: mp
         integer(int4_), dimension(:,:), allocatable, intent(out) :: array2d
         mwSize, intent(out) :: m, n
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwPointer :: dataptr
         integer(int8_), dimension(:,:), allocatable :: tempint8
         real(real4_), dimension(:,:), allocatable :: tempreal4
         real(real8_), dimension(:,:), allocatable :: tempreal8
         integer :: st
         integer :: i, j

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         if(MATLAB_is_complex(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected integer matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         endif

         if(MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected dense matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsDense", errmsg)
         endif

         deallocate(array2d, stat=st)
         allocate(array2d(m,n), stat=st)
         if(st.ne.0) goto 100

         dataptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_int4)
            call mxCopyPtrToInteger4(dataptr, array2d, m*n)
         case(mx_class_int8)
            allocate(tempint8(m,n), stat=st)
            if(st.ne.0) goto 100
            call mxCopyPtrToInteger8(dataptr, tempint8, m*n)
            array2d(:,:) = tempint8(:,:)
         case(mx_class_real4)
            allocate(tempreal4(m,n), stat=st)
            if(st.ne.0) goto 100
            call mxCopyPtrToReal4(dataptr, tempreal4, m*n)
            do i = 1, m
               do j = 1, n
                  array2d(i,j) = int(tempreal4(i,j))
                  if(real(array2d(i,j)) .ne. tempreal4(i,j)) then
                     write(errmsg, "(3a)") &
                        "Error in argument ", argname, &
                        ". Expected integer data."
                     call MATLAB_context_error("hsl_matlab:ExpectsInteger", &
                        errmsg)
                  endif
               end do
            end do
         case(mx_class_real8)
            allocate(tempreal8(m,n), stat=st)
            if(st.ne.0) goto 100
            call mxCopyPtrToReal8(dataptr, tempreal8, m*n)
            do i = 1, m
               do j = 1, n
                  array2d(i,j) = int(tempreal8(i,j))
                  if(real(array2d(i,j)) .ne. tempreal8(i,j)) then
                     write(errmsg, "(3a)") &
                        "Error in argument ", argname, &
                        ". Expected integer data."
                     call MATLAB_context_error("hsl_matlabLExpectsInteger", &
                        errmsg)
                  endif
               end do
            end do
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected integer data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select

         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)

      end subroutine matlab_to_fortran_array2d_int4
      subroutine matlab_to_fortran_array2d_int8(mp, array2d, m, n, argname)
         mwPointer :: mp
         integer(int8_), dimension(:,:), allocatable, intent(out) :: array2d
         mwSize, intent(out) :: m, n
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwPointer :: dataptr
         integer(int4_), dimension(:,:), allocatable :: tempint4
         real(real4_), dimension(:,:), allocatable :: tempreal4
         real(real8_), dimension(:,:), allocatable :: tempreal8
         integer :: st
         integer :: i, j

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         if(MATLAB_is_complex(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected integer matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         endif

         if(MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected dense matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsDense", errmsg)
         endif

         deallocate(array2d, stat=st)
         allocate(array2d(m,n), stat=st)
         if(st.ne.0) goto 100

         dataptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_int8)
            call mxCopyPtrToInteger8(dataptr, array2d, m*n)
         case(mx_class_int4)
            allocate(tempint4(m,n), stat=st)
            if(st.ne.0) goto 100
            call mxCopyPtrToInteger4(dataptr, tempint4, m*n)
            array2d(:,:) = tempint4(:,:)
         case(mx_class_real4)
            allocate(tempreal4(m,n), stat=st)
            if(st.ne.0) goto 100
            call mxCopyPtrToReal4(dataptr, tempreal4, m*n)
            do i = 1, m
               do j = 1, n
                  array2d(i,j) = int(tempreal4(i,j))
                  if(real(array2d(i,j)) .ne. tempreal4(i,j)) then
                     write(errmsg, "(3a)") &
                        "Error in argument ", argname, &
                        ". Expected integer data."
                     call MATLAB_context_error("hsl_matlab:ExpectsInteger", &
                        errmsg)
                  endif
               end do
            end do
         case(mx_class_real8)
            allocate(tempreal8(m,n), stat=st)
            if(st.ne.0) goto 100
            call mxCopyPtrToReal8(dataptr, tempreal8, m*n)
            do i = 1, m
               do j = 1, n
                  array2d(i,j) = int(tempreal8(i,j))
                  if(real(array2d(i,j)) .ne. tempreal8(i,j)) then
                     write(errmsg, "(3a)") &
                        "Error in argument ", argname, &
                        ". Expected integer data."
                     call MATLAB_context_error("hsl_matlabLExpectsInteger", &
                        errmsg)
                  endif
               end do
            end do
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected integer data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select

         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_array2d_int8
      subroutine matlab_to_fortran_array2d_real4(mp, array2d, m, n, argname)
         mwPointer :: mp
         real(real4_), dimension(:,:), allocatable, intent(out) :: array2d
         mwSize, intent(out) :: m, n
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwPointer :: dataptr
         integer(int4_), dimension(:,:), allocatable :: tempint4
         integer(int8_), dimension(:,:), allocatable :: tempint8
         real(real8_), dimension(:,:), allocatable :: tempreal8
         integer :: st

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         if(MATLAB_is_complex(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected real matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsReal", errmsg)
         endif

         if(MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected dense matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsDense", errmsg)
         endif

         deallocate(array2d, stat=st)
         allocate(array2d(m,n), stat=st)
         if(st.ne.0) goto 100

         dataptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_int4)
            allocate(tempint4(m,n), stat=st)
            if(st.ne.0) goto 100
            call mxCopyPtrToInteger4(dataptr, tempint4, m*n)
            array2d(:,:) = tempint4(:,:)
         case(mx_class_int8)
            allocate(tempint8(m,n), stat=st)
            if(st.ne.0) goto 100
            call mxCopyPtrToInteger8(dataptr, tempint8, m*n)
            array2d(:,:) = tempint8(:,:)
         case(mx_class_real4)
            call mxCopyPtrToReal4(dataptr, array2d, m*n)
         case(mx_class_real8)
            allocate(tempreal8(m,n), stat=st)
            if(st.ne.0) goto 100
            call mxCopyPtrToReal8(dataptr, tempreal8, m*n)
            array2d(:,:) = tempreal8(:,:)
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected real data."
            call MATLAB_context_error("hsl_matlab:ExpectsReal", errmsg)
         end select

         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_array2d_real4
      subroutine matlab_to_fortran_array2d_real8(mp, array2d, m, n, argname)
         mwPointer :: mp
         real(real8_), dimension(:,:), allocatable, intent(out) :: array2d
         mwSize, intent(out) :: m, n
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwPointer :: dataptr
         integer(int4_), dimension(:,:), allocatable :: tempint4
         integer(int8_), dimension(:,:), allocatable :: tempint8
         real(real4_), dimension(:,:), allocatable :: tempreal4
         integer :: st

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         if(MATLAB_is_complex(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected real matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsReal", errmsg)
         endif

         if(MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected dense matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsDense", errmsg)
         endif

         deallocate(array2d, stat=st)
         allocate(array2d(m,n), stat=st)
         if(st.ne.0) goto 100

         dataptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_int4)
            allocate(tempint4(m,n), stat=st)
            if(st.ne.0) goto 100
            call mxCopyPtrToInteger4(dataptr, tempint4, m*n)
            array2d(:,:) = tempint4(:,:)
         case(mx_class_int8)
            allocate(tempint8(m,n), stat=st)
            if(st.ne.0) goto 100
            call mxCopyPtrToInteger8(dataptr, tempint8, m*n)
            array2d(:,:) = tempint8(:,:)
         case(mx_class_real4)
            allocate(tempreal4(m,n), stat=st)
            if(st.ne.0) goto 100
            call mxCopyPtrToReal4(dataptr, tempreal4, m*n)
            array2d(:,:) = tempreal4(:,:)
         case(mx_class_real8)
            call mxCopyPtrToReal8(dataptr, array2d, m*n)
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected real data."
            call MATLAB_context_error("hsl_matlab:ExpectsReal", errmsg)
         end select

         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_array2d_real8
      subroutine matlab_to_fortran_array2d_complex8(mp, array2d, m, n, argname)
         mwPointer :: mp
         complex(real4_), dimension(:,:), allocatable, intent(out) :: array2d
         mwSize, intent(out) :: m, n
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwPointer :: realptr, imagptr
         real(real4_), dimension(:,:), allocatable :: tempreal4
         real(real8_), dimension(:,:), allocatable :: tempreal8
         complex(real8_), dimension(:,:), allocatable :: tempcomplex16
         integer :: st

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         if(MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected dense matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsDense", errmsg)
         endif

         deallocate(array2d, stat=st)
         allocate(array2d(m,n), stat=st)
         if(st.ne.0) goto 100

         realptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_real4)
            if(MATLAB_is_complex(mp)) then
               imagptr = MATLAB_get_imag_ptr(mp)
               call mxCopyPtrToComplex8(realptr, imagptr, array2d, m*n)
            else
               allocate(tempreal4(m,n), stat=st)
               if(st.ne.0) goto 100
               call mxCopyPtrToReal4(realptr, tempreal4, m*n)
               array2d(:,:) = tempreal4(:,:)
            endif
         case(mx_class_real8)
            if(MATLAB_is_complex(mp)) then
               imagptr = MATLAB_get_imag_ptr(mp)
               allocate(tempcomplex16(m,n), stat=st)
               if(st.ne.0) goto 100
               call mxCopyPtrToComplex16(realptr, imagptr, tempcomplex16, m*n)
               array2d(:,:) = tempcomplex16(:,:)
            else
               allocate(tempreal8(m,n), stat=st)
               if(st.ne.0) goto 100
               call mxCopyPtrToReal8(realptr, tempreal8, m*n)
               array2d(:,:) = tempreal8(:,:)
            endif
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected complex matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsComplex", errmsg)
         end select

         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_array2d_complex8
      subroutine matlab_to_fortran_array2d_complex16(mp, array2d, m, n, &
            argname)
         mwPointer :: mp
         complex(real8_), dimension(:,:), allocatable, intent(out) :: array2d
         mwSize, intent(out) :: m, n
         character(len=*), intent(in) :: argname

         character(len=200) :: errmsg
         mwPointer :: realptr, imagptr
         real(real4_), dimension(:,:), allocatable :: tempreal4
         real(real8_), dimension(:,:), allocatable :: tempreal8
         complex(real4_), dimension(:,:), allocatable :: tempcomplex8
         integer :: st

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif
         
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         if(MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a,2i4)") &
               "Error in argument ", argname, ". Expected dense matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsDense", errmsg)
         endif

         deallocate(array2d, stat=st)
         allocate(array2d(m,n), stat=st)
         if(st.ne.0) goto 100

         realptr = MATLAB_get_ptr(mp)
         select case(mxGetClassID(mp))
         case(mx_class_real4)
            if(MATLAB_is_complex(mp)) then
               imagptr = MATLAB_get_imag_ptr(mp)
               allocate(tempcomplex8(m,n), stat=st)
               if(st.ne.0) goto 100
               call mxCopyPtrToComplex8(realptr, imagptr, tempcomplex8, m*n)
               array2d(:,:) = tempcomplex8(:,:)
            else
               allocate(tempreal4(m,n), stat=st)
               if(st.ne.0) goto 100
               call mxCopyPtrToReal4(realptr, tempreal4, m*n)
               array2d(:,:) = tempreal4(:,:)
            endif
         case(mx_class_real8)
            if(MATLAB_is_complex(mp)) then
               imagptr = MATLAB_get_imag_ptr(mp)
               call mxCopyPtrToComplex16(realptr, imagptr, array2d, m*n)
            else
               allocate(tempreal8(m,n), stat=st)
               if(st.ne.0) goto 100
               call mxCopyPtrToReal8(realptr, tempreal8, m*n)
               array2d(:,:) = tempreal8(:,:)
            endif
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected complex matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsComplex", errmsg)
         end select

         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_array2d_complex16

      subroutine matlab_to_fortran_sparse_int4_int4(mp, m, n, ptr, row, val, &
            argname)
         mwPointer :: mp
         mwSize, intent(out) :: m, n
         integer(int4_), dimension(:), allocatable :: ptr
         integer, dimension(:), allocatable :: row ! default integer
         integer(int4_), dimension(:), allocatable :: val
         character(len=*), intent(in) :: argname

         mwPointer :: ptrptr, rowptr, valptr
         mwIndex :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         integer :: i
         character(len=200) :: errmsg

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif

         if(.not. MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, &
               ". Expected sparse matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsSparse", errmsg)
         endif

         ! Ensure arrays are deallocated before we start
         deallocate(ptr, stat=st)
         deallocate(row, stat=st)
         deallocate(val, stat=st)

         ! Get matrix dimensions
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         ! Get matrix pointers
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         valptr = MATLAB_get_ptr(mp)

         ! Copy out column pointers
         allocate(ptr(n+1), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(ptr))
            call MATLAB_copy_from_ptr(ptrptr, ptr, n+1_mws_)
         case default
            allocate(tempmwi(n+1), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(ptrptr, tempmwi, n+1_mws_)
            ptr(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         ptr(:) = ptr(:) + 1
         nz = ptr(n+1)-1

         ! Copy row indices
         allocate(row(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(row))
            call MATLAB_copy_from_ptr(rowptr, row, nz)
         case default
            allocate(tempmwi(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(rowptr, tempmwi, nz)
            row(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         row(:) = row(:) + 1

         ! Copy values
         ! Note, at present MATLAB sparse only offered in double and logical
         allocate(val(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mxGetClassID(mp))
         case(mx_class_real8)
            allocate(tempreal8(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(valptr, tempreal8, nz)
            do i = 1, nz
               val(i) = int(tempreal8(i))
               if(real(val(i)) .ne. tempreal8(i)) then
                  write(errmsg, "(3a)") &
                     "Error in argument ", argname, &
                     ". Expected integer data."
                  call MATLAB_context_error("hsl_matlabLExpectsInteger", errmsg)
               endif
            end do
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected integer data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select
         
         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_sparse_int4_int4
      subroutine matlab_to_fortran_sparse_int4_int8(mp, m, n, ptr, row, val, &
            argname)
         mwPointer :: mp
         mwSize, intent(out) :: m, n
         integer(int4_), dimension(:), allocatable :: ptr
         integer, dimension(:), allocatable :: row ! default integer
         integer(int8_), dimension(:), allocatable :: val
         character(len=*), intent(in) :: argname

         mwPointer :: ptrptr, rowptr, valptr
         mwIndex :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         integer :: i
         character(len=200) :: errmsg

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif

         if(.not. MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, &
               ". Expected sparse matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsSparse", errmsg)
         endif

         ! Ensure arrays are deallocated before we start
         deallocate(ptr, stat=st)
         deallocate(row, stat=st)
         deallocate(val, stat=st)

         ! Get matrix dimensions
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         ! Get matrix pointers
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         valptr = MATLAB_get_ptr(mp)

         ! Copy out column pointers
         allocate(ptr(n+1), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(ptr))
            call MATLAB_copy_from_ptr(ptrptr, ptr, n+1_mws_)
         case default
            allocate(tempmwi(n+1), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(ptrptr, tempmwi, n+1_mws_)
            ptr(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         ptr(:) = ptr(:) + 1
         nz = ptr(n+1)-1

         ! Copy row indices
         allocate(row(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(row))
            call MATLAB_copy_from_ptr(rowptr, row, nz)
         case default
            allocate(tempmwi(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(rowptr, tempmwi, nz)
            row(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         row(:) = row(:) + 1

         ! Copy values
         ! Note, at present MATLAB sparse only offered in double and logical
         allocate(val(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mxGetClassID(mp))
         case(mx_class_real8)
            allocate(tempreal8(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(valptr, tempreal8, nz)
            do i = 1, nz
               val(i) = int(tempreal8(i))
               if(real(val(i)) .ne. tempreal8(i)) then
                  write(errmsg, "(3a)") &
                     "Error in argument ", argname, &
                     ". Expected integer data."
                  call MATLAB_context_error("hsl_matlabLExpectsInteger", errmsg)
               endif
            end do
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected integer data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select
         
         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_sparse_int4_int8
      subroutine matlab_to_fortran_sparse_int4_real4(mp, m, n, ptr, row, val, &
            argname)
         mwPointer :: mp
         mwSize, intent(out) :: m, n
         integer(int4_), dimension(:), allocatable :: ptr
         integer, dimension(:), allocatable :: row ! default integer
         real(real4_), dimension(:), allocatable :: val
         character(len=*), intent(in) :: argname

         mwPointer :: ptrptr, rowptr, valptr
         mwIndex :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         character(len=200) :: errmsg

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif

         if(.not. MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, &
               ". Expected sparse matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsSparse", errmsg)
         endif

         ! Ensure arrays are deallocated before we start
         deallocate(ptr, stat=st)
         deallocate(row, stat=st)
         deallocate(val, stat=st)

         ! Get matrix dimensions
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         ! Get matrix pointers
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         valptr = MATLAB_get_ptr(mp)

         ! Copy out column pointers
         allocate(ptr(n+1), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(ptr))
            call MATLAB_copy_from_ptr(ptrptr, ptr, n+1_mws_)
         case default
            allocate(tempmwi(n+1), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(ptrptr, tempmwi, n+1_mws_)
            ptr(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         ptr(:) = ptr(:) + 1
         nz = ptr(n+1)-1

         ! Copy row indices
         allocate(row(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(row))
            call MATLAB_copy_from_ptr(rowptr, row, nz)
         case default
            allocate(tempmwi(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(rowptr, tempmwi, nz)
            row(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         row(:) = row(:) + 1

         ! Copy values
         ! Note, at present MATLAB sparse only offered in double and logical
         allocate(val(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mxGetClassID(mp))
         case(mx_class_real8)
            allocate(tempreal8(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(valptr, tempreal8, nz)
            val(:) = tempreal8(:)
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected real data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select
         
         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_sparse_int4_real4
      subroutine matlab_to_fortran_sparse_int4_real8(mp, m, n, ptr, row, val, &
            argname)
         mwPointer :: mp
         mwSize, intent(out) :: m, n
         integer(int4_), dimension(:), allocatable :: ptr
         integer, dimension(:), allocatable :: row ! default integer
         real(real8_), dimension(:), allocatable :: val
         character(len=*), intent(in) :: argname

         mwPointer :: ptrptr, rowptr, valptr
         mwIndex :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         integer :: st
         character(len=200) :: errmsg

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif

         if(.not. MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, &
               ". Expected sparse matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsSparse", errmsg)
         endif

         ! Ensure arrays are deallocated before we start
         deallocate(ptr, stat=st)
         deallocate(row, stat=st)
         deallocate(val, stat=st)

         ! Get matrix dimensions
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         ! Get matrix pointers
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         valptr = MATLAB_get_ptr(mp)

         ! Copy out column pointers
         allocate(ptr(n+1), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(ptr))
            call MATLAB_copy_from_ptr(ptrptr, ptr, n+1_mws_)
         case default
            allocate(tempmwi(n+1), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(ptrptr, tempmwi, n+1_mws_)
            ptr(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         ptr(:) = ptr(:) + 1
         nz = ptr(n+1)-1

         ! Copy row indices
         allocate(row(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(row))
            call MATLAB_copy_from_ptr(rowptr, row, nz)
         case default
            allocate(tempmwi(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(rowptr, tempmwi, nz)
            row(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         row(:) = row(:) + 1

         ! Copy values
         ! Note, at present MATLAB sparse only offered in double and logical
         allocate(val(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mxGetClassID(mp))
         case(mx_class_real8)
            call MATLAB_copy_from_ptr(valptr, val, nz)
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected real data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select
         
         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_sparse_int4_real8
      subroutine matlab_to_fortran_sparse_int4_complex8(mp, m, n, ptr, row, &
            val, argname)
         mwPointer :: mp
         mwSize, intent(out) :: m, n
         integer(int4_), dimension(:), allocatable :: ptr
         integer, dimension(:), allocatable :: row ! default integer
         complex(real4_), dimension(:), allocatable :: val
         character(len=*), intent(in) :: argname

         mwPointer :: ptrptr, rowptr, realptr, imagptr
         mwIndex :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         complex(real8_), dimension(:), allocatable :: tempcomplex16
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         character(len=200) :: errmsg

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif

         if(.not. MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, &
               ". Expected sparse matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsSparse", errmsg)
         endif

         ! Ensure arrays are deallocated before we start
         deallocate(ptr, stat=st)
         deallocate(row, stat=st)
         deallocate(val, stat=st)

         ! Get matrix dimensions
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         ! Get matrix pointers
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         realptr = MATLAB_get_ptr(mp)

         ! Copy out column pointers
         allocate(ptr(n+1), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(ptr))
            call MATLAB_copy_from_ptr(ptrptr, ptr, n+1_mws_)
         case default
            allocate(tempmwi(n+1), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(ptrptr, tempmwi, n+1_mws_)
            ptr(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         ptr(:) = ptr(:) + 1
         nz = ptr(n+1)-1

         ! Copy row indices
         allocate(row(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(row))
            call MATLAB_copy_from_ptr(rowptr, row, nz)
         case default
            allocate(tempmwi(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(rowptr, tempmwi, nz)
            row(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         row(:) = row(:) + 1

         ! Copy values
         ! Note, at present MATLAB sparse only offered in double and logical
         allocate(val(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mxGetClassID(mp))
         case(mx_class_real8)
            if(MATLAB_is_complex(mp)) then
               imagptr = MATLAB_get_imag_ptr(mp)
               allocate(tempcomplex16(nz), stat=st)
               if(st.ne.0) goto 100
               call MATLAB_copy_from_ptr(realptr, imagptr, tempcomplex16, nz)
               val(:) = tempcomplex16(:)
            else
               allocate(tempreal8(nz), stat=st)
               if(st.ne.0) goto 100
               call MATLAB_copy_from_ptr(realptr, tempreal8, nz)
               val(:) = tempreal8(:)
            endif
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected real data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select
         
         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_sparse_int4_complex8
      subroutine matlab_to_fortran_sparse_int4_complex16(mp, m, n, ptr, row, &
            val, argname)
         mwPointer :: mp
         mwSize, intent(out) :: m, n
         integer(int4_), dimension(:), allocatable :: ptr
         integer, dimension(:), allocatable :: row ! default integer
         complex(real8_), dimension(:), allocatable :: val
         character(len=*), intent(in) :: argname

         mwPointer :: ptrptr, rowptr, realptr, imagptr
         mwIndex :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         character(len=200) :: errmsg

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif

         if(.not. MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, &
               ". Expected sparse matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsSparse", errmsg)
         endif

         ! Ensure arrays are deallocated before we start
         deallocate(ptr, stat=st)
         deallocate(row, stat=st)
         deallocate(val, stat=st)

         ! Get matrix dimensions
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         ! Get matrix pointers
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         realptr = MATLAB_get_ptr(mp)

         ! Copy out column pointers
         allocate(ptr(n+1), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(ptr))
            call MATLAB_copy_from_ptr(ptrptr, ptr, n+1_mws_)
         case default
            allocate(tempmwi(n+1), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(ptrptr, tempmwi, n+1_mws_)
            ptr(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         ptr(:) = ptr(:) + 1
         nz = ptr(n+1)-1

         ! Copy row indices
         allocate(row(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(row))
            call MATLAB_copy_from_ptr(rowptr, row, nz)
         case default
            allocate(tempmwi(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(rowptr, tempmwi, nz)
            row(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         row(:) = row(:) + 1

         ! Copy values
         ! Note, at present MATLAB sparse only offered in double and logical
         allocate(val(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mxGetClassID(mp))
         case(mx_class_real8)
            if(MATLAB_is_complex(mp)) then
               imagptr = MATLAB_get_imag_ptr(mp)
               call MATLAB_copy_from_ptr(realptr, imagptr, val, nz)
            else
               allocate(tempreal8(nz), stat=st)
               if(st.ne.0) goto 100
               call MATLAB_copy_from_ptr(realptr, tempreal8, nz)
               val(:) = tempreal8(:)
            endif
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected real data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select
         
         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_sparse_int4_complex16

      subroutine matlab_to_fortran_sparse_int8_int4(mp, m, n, ptr, row, val, &
            argname)
         mwPointer :: mp
         mwSize, intent(out) :: m, n
         integer(int8_), dimension(:), allocatable :: ptr
         integer, dimension(:), allocatable :: row ! default integer
         integer(int4_), dimension(:), allocatable :: val
         character(len=*), intent(in) :: argname

         mwPointer :: ptrptr, rowptr, valptr
         mwIndex :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         integer :: i
         character(len=200) :: errmsg

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif

         if(.not. MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, &
               ". Expected sparse matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsSparse", errmsg)
         endif

         ! Ensure arrays are deallocated before we start
         deallocate(ptr, stat=st)
         deallocate(row, stat=st)
         deallocate(val, stat=st)

         ! Get matrix dimensions
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         ! Get matrix pointers
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         valptr = MATLAB_get_ptr(mp)

         ! Copy out column pointers
         allocate(ptr(n+1), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(ptr))
            call MATLAB_copy_from_ptr(ptrptr, ptr, n+1_mws_)
         case default
            allocate(tempmwi(n+1), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(ptrptr, tempmwi, n+1_mws_)
            ptr(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         ptr(:) = ptr(:) + 1
         nz = ptr(n+1)-1

         ! Copy row indices
         allocate(row(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(row))
            call MATLAB_copy_from_ptr(rowptr, row, nz)
         case default
            allocate(tempmwi(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(rowptr, tempmwi, nz)
            row(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         row(:) = row(:) + 1

         ! Copy values
         ! Note, at present MATLAB sparse only offered in double and logical
         allocate(val(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mxGetClassID(mp))
         case(mx_class_real8)
            allocate(tempreal8(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(valptr, tempreal8, nz)
            do i = 1, nz
               val(i) = int(tempreal8(i))
               if(real(val(i)) .ne. tempreal8(i)) then
                  write(errmsg, "(3a)") &
                     "Error in argument ", argname, &
                     ". Expected integer data."
                  call MATLAB_context_error("hsl_matlabLExpectsInteger", errmsg)
               endif
            end do
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected integer data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select
         
         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_sparse_int8_int4
      subroutine matlab_to_fortran_sparse_int8_int8(mp, m, n, ptr, row, val, &
            argname)
         mwPointer :: mp
         mwSize, intent(out) :: m, n
         integer(int8_), dimension(:), allocatable :: ptr
         integer, dimension(:), allocatable :: row ! default integer
         integer(int8_), dimension(:), allocatable :: val
         character(len=*), intent(in) :: argname

         mwPointer :: ptrptr, rowptr, valptr
         mwIndex :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         integer :: i
         character(len=200) :: errmsg

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif

         if(.not. MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, &
               ". Expected sparse matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsSparse", errmsg)
         endif

         ! Ensure arrays are deallocated before we start
         deallocate(ptr, stat=st)
         deallocate(row, stat=st)
         deallocate(val, stat=st)

         ! Get matrix dimensions
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         ! Get matrix pointers
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         valptr = MATLAB_get_ptr(mp)

         ! Copy out column pointers
         allocate(ptr(n+1), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(ptr))
            call MATLAB_copy_from_ptr(ptrptr, ptr, n+1_mws_)
         case default
            allocate(tempmwi(n+1), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(ptrptr, tempmwi, n+1_mws_)
            ptr(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         ptr(:) = ptr(:) + 1
         nz = ptr(n+1)-1

         ! Copy row indices
         allocate(row(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(row))
            call MATLAB_copy_from_ptr(rowptr, row, nz)
         case default
            allocate(tempmwi(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(rowptr, tempmwi, nz)
            row(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         row(:) = row(:) + 1

         ! Copy values
         ! Note, at present MATLAB sparse only offered in double and logical
         allocate(val(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mxGetClassID(mp))
         case(mx_class_real8)
            allocate(tempreal8(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(valptr, tempreal8, nz)
            do i = 1, nz
               val(i) = int(tempreal8(i))
               if(real(val(i)) .ne. tempreal8(i)) then
                  write(errmsg, "(3a)") &
                     "Error in argument ", argname, &
                     ". Expected integer data."
                  call MATLAB_context_error("hsl_matlabLExpectsInteger", errmsg)
               endif
            end do
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected integer data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select
         
         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_sparse_int8_int8
      subroutine matlab_to_fortran_sparse_int8_real4(mp, m, n, ptr, row, val, &
            argname)
         mwPointer :: mp
         mwSize, intent(out) :: m, n
         integer(int8_), dimension(:), allocatable :: ptr
         integer, dimension(:), allocatable :: row ! default integer
         real(real4_), dimension(:), allocatable :: val
         character(len=*), intent(in) :: argname

         mwPointer :: ptrptr, rowptr, valptr
         mwIndex :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         character(len=200) :: errmsg

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif

         if(.not. MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, &
               ". Expected sparse matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsSparse", errmsg)
         endif

         ! Ensure arrays are deallocated before we start
         deallocate(ptr, stat=st)
         deallocate(row, stat=st)
         deallocate(val, stat=st)

         ! Get matrix dimensions
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         ! Get matrix pointers
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         valptr = MATLAB_get_ptr(mp)

         ! Copy out column pointers
         allocate(ptr(n+1), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(ptr))
            call MATLAB_copy_from_ptr(ptrptr, ptr, n+1_mws_)
         case default
            allocate(tempmwi(n+1), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(ptrptr, tempmwi, n+1_mws_)
            ptr(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         ptr(:) = ptr(:) + 1
         nz = ptr(n+1)-1

         ! Copy row indices
         allocate(row(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(row))
            call MATLAB_copy_from_ptr(rowptr, row, nz)
         case default
            allocate(tempmwi(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(rowptr, tempmwi, nz)
            row(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         row(:) = row(:) + 1

         ! Copy values
         ! Note, at present MATLAB sparse only offered in double and logical
         allocate(val(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mxGetClassID(mp))
         case(mx_class_real8)
            allocate(tempreal8(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(valptr, tempreal8, nz)
            val(:) = tempreal8(:)
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected integer data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select
         
         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_sparse_int8_real4
      subroutine matlab_to_fortran_sparse_int8_real8(mp, m, n, ptr, row, val, &
            argname)
         mwPointer :: mp
         mwSize, intent(out) :: m, n
         integer(int8_), dimension(:), allocatable :: ptr
         integer, dimension(:), allocatable :: row ! default integer
         real(real8_), dimension(:), allocatable :: val
         character(len=*), intent(in) :: argname

         mwPointer :: ptrptr, rowptr, valptr
         mwIndex :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         integer :: st
         character(len=200) :: errmsg

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif

         if(.not. MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, &
               ". Expected sparse matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsSparse", errmsg)
         endif

         ! Ensure arrays are deallocated before we start
         deallocate(ptr, stat=st)
         deallocate(row, stat=st)
         deallocate(val, stat=st)

         ! Get matrix dimensions
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         ! Get matrix pointers
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         valptr = MATLAB_get_ptr(mp)

         ! Copy out column pointers
         allocate(ptr(n+1), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(ptr))
            call MATLAB_copy_from_ptr(ptrptr, ptr, n+1_mws_)
         case default
            allocate(tempmwi(n+1), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(ptrptr, tempmwi, n+1_mws_)
            ptr(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         ptr(:) = ptr(:) + 1
         nz = ptr(n+1)-1

         ! Copy row indices
         allocate(row(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(row))
            call MATLAB_copy_from_ptr(rowptr, row, nz)
         case default
            allocate(tempmwi(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(rowptr, tempmwi, nz)
            row(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         row(:) = row(:) + 1

         ! Copy values
         ! Note, at present MATLAB sparse only offered in double and logical
         allocate(val(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mxGetClassID(mp))
         case(mx_class_real8)
            call MATLAB_copy_from_ptr(valptr, val, nz)
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected integer data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select
         
         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_sparse_int8_real8
      subroutine matlab_to_fortran_sparse_int8_complex8(mp, m, n, ptr, row, &
            val, argname)
         mwPointer :: mp
         mwSize, intent(out) :: m, n
         integer(int8_), dimension(:), allocatable :: ptr
         integer, dimension(:), allocatable :: row ! default integer
         complex(real4_), dimension(:), allocatable :: val
         character(len=*), intent(in) :: argname

         mwPointer :: ptrptr, rowptr, realptr, imagptr
         mwIndex :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         complex(real8_), dimension(:), allocatable :: tempcomplex16
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         character(len=200) :: errmsg

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif

         if(.not. MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, &
               ". Expected sparse matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsSparse", errmsg)
         endif

         ! Ensure arrays are deallocated before we start
         deallocate(ptr, stat=st)
         deallocate(row, stat=st)
         deallocate(val, stat=st)

         ! Get matrix dimensions
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         ! Get matrix pointers
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         realptr = MATLAB_get_ptr(mp)

         ! Copy out column pointers
         allocate(ptr(n+1), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(ptr))
            call MATLAB_copy_from_ptr(ptrptr, ptr, n+1_mws_)
         case default
            allocate(tempmwi(n+1), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(ptrptr, tempmwi, n+1_mws_)
            ptr(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         ptr(:) = ptr(:) + 1
         nz = ptr(n+1)-1

         ! Copy row indices
         allocate(row(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(row))
            call MATLAB_copy_from_ptr(rowptr, row, nz)
         case default
            allocate(tempmwi(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(rowptr, tempmwi, nz)
            row(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         row(:) = row(:) + 1

         ! Copy values
         ! Note, at present MATLAB sparse only offered in double and logical
         allocate(val(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mxGetClassID(mp))
         case(mx_class_real8)
            if(MATLAB_is_complex(mp)) then
               imagptr = MATLAB_get_imag_ptr(mp)
               allocate(tempcomplex16(nz), stat=st)
               if(st.ne.0) goto 100
               call MATLAB_copy_from_ptr(realptr, imagptr, tempcomplex16, nz)
               val(:) = tempcomplex16(:)
            else
               allocate(tempreal8(nz), stat=st)
               if(st.ne.0) goto 100
               call MATLAB_copy_from_ptr(realptr, tempreal8, nz)
               val(:) = tempreal8(:)
            endif
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected real data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select
         
         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_sparse_int8_complex8
      subroutine matlab_to_fortran_sparse_int8_complex16(mp, m, n, ptr, row, &
            val, argname)
         mwPointer :: mp
         mwSize, intent(out) :: m, n
         integer(int8_), dimension(:), allocatable :: ptr
         integer, dimension(:), allocatable :: row ! default integer
         complex(real8_), dimension(:), allocatable :: val
         character(len=*), intent(in) :: argname

         mwPointer :: ptrptr, rowptr, realptr, imagptr
         mwIndex :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         character(len=200) :: errmsg

         if(mp.eq.0) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Passed null pointer."
            call MATLAB_context_error("hsl_matlab:NullPtr", errmsg)
         endif

         if(.not. MATLAB_is_sparse(mp)) then
            write(errmsg, "(3a)") &
               "Error in argument ", argname, &
               ". Expected sparse matrix."
            call MATLAB_context_error("hsl_matlab:ExpectsSparse", errmsg)
         endif

         ! Ensure arrays are deallocated before we start
         deallocate(ptr, stat=st)
         deallocate(row, stat=st)
         deallocate(val, stat=st)

         ! Get matrix dimensions
         m = MATLAB_get_m(mp)
         n = MATLAB_get_n(mp)

         ! Get matrix pointers
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         realptr = MATLAB_get_ptr(mp)

         ! Copy out column pointers
         allocate(ptr(n+1), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(ptr))
            call MATLAB_copy_from_ptr(ptrptr, ptr, n+1_mws_)
         case default
            allocate(tempmwi(n+1), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(ptrptr, tempmwi, n+1_mws_)
            ptr(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         ptr(:) = ptr(:) + 1
         nz = ptr(n+1)-1

         ! Copy row indices
         allocate(row(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mwi_)
         case(kind(row))
            call MATLAB_copy_from_ptr(rowptr, row, nz)
         case default
            allocate(tempmwi(nz), stat=st)
            if(st.ne.0) goto 100
            call MATLAB_copy_from_ptr(rowptr, tempmwi, nz)
            row(:) = tempmwi(:)
            deallocate(tempmwi)
         end select
         ! Convert from 1-based to 0-based
         row(:) = row(:) + 1

         ! Copy values
         ! Note, at present MATLAB sparse only offered in double and logical
         allocate(val(nz), stat=st)
         if(st.ne.0) goto 100
         select case(mxGetClassID(mp))
         case(mx_class_real8)
            if(MATLAB_is_complex(mp)) then
               imagptr = MATLAB_get_imag_ptr(mp)
               call MATLAB_copy_from_ptr(realptr, imagptr, val, nz)
            else
               allocate(tempreal8(nz), stat=st)
               if(st.ne.0) goto 100
               call MATLAB_copy_from_ptr(realptr, tempreal8, nz)
               val(:) = tempreal8(:)
            endif
         case default
            write(errmsg, "(3a)") &
               "Error in argument ", argname, ". Expected real data."
            call MATLAB_context_error("hsl_matlab:ExpectsInteger", errmsg)
         end select
         
         return
         100 continue
         write(errmsg, "(3a)") &
            "Error processing argument ", argname, ". Allocation failed."
         call MATLAB_context_error("hsl_matlab:AllocFail", errmsg)
      end subroutine matlab_to_fortran_sparse_int8_complex16

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      mwPointer function fortran_to_matlab_scalar_int4(scalar)
         integer(int4_), intent(in) :: scalar

         fortran_to_matlab_scalar_int4 = &
            fortran_to_matlab_array1d_int4( (/scalar/) )
         end function fortran_to_matlab_scalar_int4
         mwPointer function fortran_to_matlab_array1d_int4(array1d, n)
            integer(int4_), dimension(:), intent(in) :: array1d
            mwSize, optional :: n

            mwPointer :: mp, dataptr
            mwSize :: len

            len = size(array1d)
            if(present(n)) len = n
            
            mp = mxCreateNumericArray( 2_mws_, (/ len, 1_mws_ /), &
               mxClassIDFromClassName('int32'), 0_int4_ ) 
            dataptr = MATLAB_get_ptr(mp)
            call MATLAB_copy_to_ptr(array1d, dataptr, len)

            fortran_to_matlab_array1d_int4 = mp
         end function fortran_to_matlab_array1d_int4
         mwPointer function fortran_to_matlab_array2d_int4(array2d, m, n)
            integer(int4_), dimension(:,:), intent(in) :: array2d
            mwSize, optional :: m
         mwSize, optional :: n

         mwPointer :: mp, dataptr
         mwSize :: len1, len2

         len1 = size(array2d,1)
         if(present(m)) len1 = m
         len2 = size(array2d,2)
         if(present(n)) len2 = n
         
         mp = mxCreateNumericArray( 2_mws_, (/ len1, len2 /), &
            mxClassIDFromClassName('int32'), 0_int4_ ) 
         dataptr = MATLAB_get_ptr(mp)
         call mxCopyInteger4ToPtr(array2d, dataptr, len1*len2)

         fortran_to_matlab_array2d_int4 = mp
      end function fortran_to_matlab_array2d_int4
      mwPointer function fortran_to_matlab_sparse_int4_int4(m, n, ptr, row, val)
         mwSize, intent(in) :: m, n
         integer(int4_), dimension(n+1), intent(in) :: ptr
         integer, dimension(ptr(n+1)-1), intent(in) :: row ! default integer
         integer(int4_), dimension(ptr(n+1)-1), intent(in) :: val

         mwPointer :: mp, ptrptr, rowptr, valptr
         mwSize :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         
         ! Setup a sparse matrix. This allocates ptr, row, val on matlab side
         nz = ptr(n+1)-1
         mp = mxCreateSparse(m, n, nz, 0_int4_)
         
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         valptr = MATLAB_get_ptr(mp)

         ! Convert ptr from 1-based to 0-based then copy out
         allocate(tempmwi(n+1), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = ptr(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, ptrptr, n+1_mws_)
         deallocate(tempmwi)

         ! Convert row from 1-based to 0-based then copy out
         allocate(tempmwi(nz), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = row(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, rowptr, nz)
         deallocate(tempmwi)

         ! Copy out val. mxCreateSparse only creates double prec matrices
         allocate(tempreal8(nz), stat=st)
         if(st.ne.0) goto 100
         tempreal8(:) = val(:)
         call MATLAB_copy_to_ptr(tempreal8, valptr, nz)

         fortran_to_matlab_sparse_int4_int4 = mp

         return
         100 continue
         call MATLAB_context_error("hsl_matlab:AllocFail", "Allocation failed")
      end function fortran_to_matlab_sparse_int4_int4
      mwPointer function fortran_to_matlab_sparse_int8_int4(m, n, ptr, row, val)
         mwSize, intent(in) :: m, n
         integer(int8_), dimension(n+1), intent(in) :: ptr
         integer, dimension(ptr(n+1)-1), intent(in) :: row ! default integer
         integer(int4_), dimension(ptr(n+1)-1), intent(in) :: val

         mwPointer :: mp, ptrptr, rowptr, valptr
         mwSize :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         
         ! Setup a sparse matrix. This allocates ptr, row, val on matlab side
         nz = ptr(n+1)-1
         mp = mxCreateSparse(m, n, nz, 0_int4_)
         
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         valptr = MATLAB_get_ptr(mp)

         ! Convert ptr from 1-based to 0-based then copy out
         allocate(tempmwi(n+1), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = ptr(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, ptrptr, n+1_mws_)
         deallocate(tempmwi)

         ! Convert row from 1-based to 0-based then copy out
         allocate(tempmwi(nz), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = row(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, rowptr, nz)
         deallocate(tempmwi)

         ! Copy out val. mxCreateSparse only creates double prec matrices
         allocate(tempreal8(nz), stat=st)
         if(st.ne.0) goto 100
         tempreal8(:) = val(:)
         call MATLAB_copy_to_ptr(tempreal8, valptr, nz)

         fortran_to_matlab_sparse_int8_int4 = mp

         return
         100 continue
         call MATLAB_context_error("hsl_matlab:AllocFail", "Allocation failed")
      end function fortran_to_matlab_sparse_int8_int4
         
      mwPointer function fortran_to_matlab_scalar_int8(scalar)
         integer(int8_), intent(in) :: scalar

         fortran_to_matlab_scalar_int8 = &
            fortran_to_matlab_array1d_int8( (/scalar/) )
      end function fortran_to_matlab_scalar_int8
      mwPointer function fortran_to_matlab_array1d_int8(array1d, n)
         integer(int8_), dimension(:), intent(in) :: array1d
         mwSize, optional :: n

         mwPointer :: mp, dataptr
         mwSize :: len

         len = size(array1d)
         if(present(n)) len = n
         
         mp = mxCreateNumericArray( 2_mws_, (/ len, 1_mws_ /), &
            mxClassIDFromClassName('int64'), 0_int4_ ) 
         dataptr = MATLAB_get_ptr(mp)
         call MATLAB_copy_to_ptr(array1d, dataptr, len)

         fortran_to_matlab_array1d_int8 = mp
      end function fortran_to_matlab_array1d_int8
      mwPointer function fortran_to_matlab_array2d_int8(array2d, m, n)
         integer(int8_), dimension(:,:), intent(in) :: array2d
         mwSize, optional :: m
         mwSize, optional :: n

         mwPointer :: mp, dataptr
         mwSize :: len1, len2

         len1 = size(array2d,1)
         if(present(m)) len1 = m
         len2 = size(array2d,2)
         if(present(n)) len2 = n
         
         mp = mxCreateNumericArray( 2_mws_, (/ len1, len2 /), &
            mxClassIDFromClassName('int64'), 0_int4_ ) 
         dataptr = MATLAB_get_ptr(mp)
         call mxCopyInteger8ToPtr(array2d, dataptr, len1*len2)

         fortran_to_matlab_array2d_int8 = mp
      end function fortran_to_matlab_array2d_int8
      mwPointer function fortran_to_matlab_sparse_int4_int8(m, n, ptr, row, val)
         mwSize, intent(in) :: m, n
         integer(int4_), dimension(n+1), intent(in) :: ptr
         integer, dimension(ptr(n+1)-1), intent(in) :: row ! default integer
         integer(int8_), dimension(ptr(n+1)-1), intent(in) :: val

         mwPointer :: mp, ptrptr, rowptr, valptr
         mwSize :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         
         ! Setup a sparse matrix. This allocates ptr, row, val on matlab side
         nz = ptr(n+1)-1
         mp = mxCreateSparse(m, n, nz, 0_int4_)
         
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         valptr = MATLAB_get_ptr(mp)

         ! Convert ptr from 1-based to 0-based then copy out
         allocate(tempmwi(n+1), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = ptr(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, ptrptr, n+1_mws_)
         deallocate(tempmwi)

         ! Convert row from 1-based to 0-based then copy out
         allocate(tempmwi(nz), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = row(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, rowptr, nz)
         deallocate(tempmwi)

         ! Copy out val. mxCreateSparse only creates double prec matrices
         allocate(tempreal8(nz), stat=st)
         if(st.ne.0) goto 100
         tempreal8(:) = val(:)
         call MATLAB_copy_to_ptr(tempreal8, valptr, nz)

         fortran_to_matlab_sparse_int4_int8 = mp

         return
         100 continue
         call MATLAB_context_error("hsl_matlab:AllocFail", "Allocation failed")
      end function fortran_to_matlab_sparse_int4_int8
      mwPointer function fortran_to_matlab_sparse_int8_int8(m, n, ptr, row, val)
         mwSize, intent(in) :: m, n
         integer(int8_), dimension(n+1), intent(in) :: ptr
         integer, dimension(ptr(n+1)-1), intent(in) :: row ! default integer
         integer(int8_), dimension(ptr(n+1)-1), intent(in) :: val

         mwPointer :: mp, ptrptr, rowptr, valptr
         mwSize :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         
         ! Setup a sparse matrix. This allocates ptr, row, val on matlab side
         nz = ptr(n+1)-1
         mp = mxCreateSparse(m, n, nz, 0_int4_)
         
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         valptr = MATLAB_get_ptr(mp)

         ! Convert ptr from 1-based to 0-based then copy out
         allocate(tempmwi(n+1), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = ptr(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, ptrptr, n+1_mws_)
         deallocate(tempmwi)

         ! Convert row from 1-based to 0-based then copy out
         allocate(tempmwi(nz), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = row(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, rowptr, nz)
         deallocate(tempmwi)

         ! Copy out val. mxCreateSparse only creates double prec matrices
         allocate(tempreal8(nz), stat=st)
         if(st.ne.0) goto 100
         tempreal8(:) = val(:)
         call MATLAB_copy_to_ptr(tempreal8, valptr, nz)

         fortran_to_matlab_sparse_int8_int8 = mp

         return
         100 continue
         call MATLAB_context_error("hsl_matlab:AllocFail", "Allocation failed")
      end function fortran_to_matlab_sparse_int8_int8

      mwPointer function fortran_to_matlab_scalar_real4(scalar)
         real(real4_), intent(in) :: scalar

         fortran_to_matlab_scalar_real4 = &
            fortran_to_matlab_array1d_real4( (/scalar/) )
      end function fortran_to_matlab_scalar_real4
      mwPointer function fortran_to_matlab_array1d_real4(array1d, n)
         real(real4_), dimension(:), intent(in) :: array1d
         mwSize, optional :: n

         mwPointer :: mp, dataptr
         mwSize :: len

         len = size(array1d)
         if(present(n)) len = n
         
         mp = mxCreateNumericArray( 2_mws_, (/ len, 1_mws_ /), &
            mxClassIDFromClassName('single'), 0_int4_ ) 
         dataptr = MATLAB_get_ptr(mp)
         call MATLAB_copy_to_ptr(array1d, dataptr, len)

         fortran_to_matlab_array1d_real4 = mp
      end function fortran_to_matlab_array1d_real4
      mwPointer function fortran_to_matlab_array2d_real4(array2d, m, n)
         real(real4_), dimension(:,:), intent(in) :: array2d
         mwSize, optional :: m
         mwSize, optional :: n

         mwPointer :: mp, dataptr
         mwSize :: len1, len2

         len1 = size(array2d,1)
         if(present(m)) len1 = m
         len2 = size(array2d,2)
         if(present(n)) len2 = n
         
         mp = mxCreateNumericArray( 2_mws_, (/ len1, len2 /), &
            mxClassIDFromClassName('single'), 0_int4_ ) 
         dataptr = MATLAB_get_ptr(mp)
         call mxCopyReal4ToPtr(array2d, dataptr, len1*len2)

         fortran_to_matlab_array2d_real4 = mp
      end function fortran_to_matlab_array2d_real4
      mwPointer function fortran_to_matlab_sparse_int4_real4(m, n, ptr, row, val)
         mwSize, intent(in) :: m, n
         integer(int4_), dimension(n+1), intent(in) :: ptr
         integer, dimension(ptr(n+1)-1), intent(in) :: row ! default integer
         real(real4_), dimension(ptr(n+1)-1), intent(in) :: val

         mwPointer :: mp, ptrptr, rowptr, valptr
         mwSize :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         
         ! Setup a sparse matrix. This allocates ptr, row, val on matlab side
         nz = ptr(n+1)-1
         mp = mxCreateSparse(m, n, nz, 0_int4_)
         
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         valptr = MATLAB_get_ptr(mp)

         ! Convert ptr from 1-based to 0-based then copy out
         allocate(tempmwi(n+1), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = ptr(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, ptrptr, n+1_mws_)
         deallocate(tempmwi)

         ! Convert row from 1-based to 0-based then copy out
         allocate(tempmwi(nz), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = row(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, rowptr, nz)
         deallocate(tempmwi)

         ! Copy out val. mxCreateSparse only creates double prec matrices
         allocate(tempreal8(nz), stat=st)
         if(st.ne.0) goto 100
         tempreal8(:) = val(:)
         call MATLAB_copy_to_ptr(tempreal8, valptr, nz)

         fortran_to_matlab_sparse_int4_real4 = mp

         return
         100 continue
         call MATLAB_context_error("hsl_matlab:AllocFail", "Allocation failed")
      end function fortran_to_matlab_sparse_int4_real4
      mwPointer function fortran_to_matlab_sparse_int8_real4(m, n, ptr, row, val)
         mwSize, intent(in) :: m, n
         integer(int8_), dimension(n+1), intent(in) :: ptr
         integer, dimension(ptr(n+1)-1), intent(in) :: row ! default integer
         real(real4_), dimension(ptr(n+1)-1), intent(in) :: val

         mwPointer :: mp, ptrptr, rowptr, valptr
         mwSize :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         real(real8_), dimension(:), allocatable :: tempreal8
         integer :: st
         
         ! Setup a sparse matrix. This allocates ptr, row, val on matlab side
         nz = ptr(n+1)-1
         mp = mxCreateSparse(m, n, nz, 0_int4_)
         
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         valptr = MATLAB_get_ptr(mp)

         ! Convert ptr from 1-based to 0-based then copy out
         allocate(tempmwi(n+1), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = ptr(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, ptrptr, n+1_mws_)
         deallocate(tempmwi)

         ! Convert row from 1-based to 0-based then copy out
         allocate(tempmwi(nz), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = row(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, rowptr, nz)
         deallocate(tempmwi)

         ! Copy out val. mxCreateSparse only creates double prec matrices
         allocate(tempreal8(nz), stat=st)
         if(st.ne.0) goto 100
         tempreal8(:) = val(:)
         call MATLAB_copy_to_ptr(tempreal8, valptr, nz)

         fortran_to_matlab_sparse_int8_real4 = mp

         return
         100 continue
         call MATLAB_context_error("hsl_matlab:AllocFail", "Allocation failed")
      end function fortran_to_matlab_sparse_int8_real4
      mwPointer function fortran_to_matlab_sparse_int4_real8(m, n, ptr, row, val)
         mwSize, intent(in) :: m, n
         integer(int4_), dimension(n+1), intent(in) :: ptr
         integer, dimension(ptr(n+1)-1), intent(in) :: row ! default integer
         real(real8_), dimension(ptr(n+1)-1), intent(in) :: val

         mwPointer :: mp, ptrptr, rowptr, valptr
         mwSize :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         integer :: st
         
         ! Setup a sparse matrix. This allocates ptr, row, val on matlab side
         nz = ptr(n+1)-1
         mp = mxCreateSparse(m, n, nz, 0_int4_)
         
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         valptr = MATLAB_get_ptr(mp)

         ! Convert ptr from 1-based to 0-based then copy out
         allocate(tempmwi(n+1), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = ptr(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, ptrptr, n+1_mws_)
         deallocate(tempmwi)

         ! Convert row from 1-based to 0-based then copy out
         allocate(tempmwi(nz), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = row(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, rowptr, nz)
         deallocate(tempmwi)

         ! Copy out val. mxCreateSparse only creates double prec matrices
         call MATLAB_copy_to_ptr(val, valptr, nz)

         fortran_to_matlab_sparse_int4_real8 = mp

         return
         100 continue
         call MATLAB_context_error("hsl_matlab:AllocFail", "Allocation failed")
      end function fortran_to_matlab_sparse_int4_real8
      mwPointer function fortran_to_matlab_sparse_int8_real8(m, n, ptr, row, &
            val)
         mwSize, intent(in) :: m, n
         integer(int8_), dimension(n+1), intent(in) :: ptr
         integer, dimension(ptr(n+1)-1), intent(in) :: row ! default integer
         real(real8_), dimension(ptr(n+1)-1), intent(in) :: val

         mwPointer :: mp, ptrptr, rowptr, valptr
         mwSize :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         integer :: st
         
         ! Setup a sparse matrix. This allocates ptr, row, val on matlab side
         nz = ptr(n+1)-1
         mp = mxCreateSparse(m, n, nz, 0_int4_)
         
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         valptr = MATLAB_get_ptr(mp)

         ! Convert ptr from 1-based to 0-based then copy out
         allocate(tempmwi(n+1), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = ptr(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, ptrptr, n+1_mws_)
         deallocate(tempmwi)

         ! Convert row from 1-based to 0-based then copy out
         allocate(tempmwi(nz), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = row(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, rowptr, nz)
         deallocate(tempmwi)

         ! Copy out val. mxCreateSparse only creates double prec matrices
         call MATLAB_copy_to_ptr(val, valptr, nz)

         fortran_to_matlab_sparse_int8_real8 = mp

         return
         100 continue
         call MATLAB_context_error("hsl_matlab:AllocFail", "Allocation failed")
      end function fortran_to_matlab_sparse_int8_real8

      mwPointer function fortran_to_matlab_scalar_real8(scalar)
         real(real8_), intent(in) :: scalar

         fortran_to_matlab_scalar_real8 = &
            fortran_to_matlab_array1d_real8( (/scalar/) )
      end function fortran_to_matlab_scalar_real8
      mwPointer function fortran_to_matlab_array1d_real8(array1d, n)
         real(real8_), dimension(:), intent(in) :: array1d
         mwSize, optional :: n

         mwPointer :: mp, dataptr
         mwSize :: len

         len = size(array1d)
         if(present(n)) len = n
         
         mp = mxCreateNumericArray( 2_mws_, (/ len, 1_mws_ /), &
            mxClassIDFromClassName('double'), 0_int4_ ) 
         dataptr = MATLAB_get_ptr(mp)
         call MATLAB_copy_to_ptr(array1d, dataptr, len)

         fortran_to_matlab_array1d_real8 = mp
      end function fortran_to_matlab_array1d_real8
      mwPointer function fortran_to_matlab_array2d_real8(array2d, m, n)
         real(real8_), dimension(:,:), intent(in) :: array2d
         mwSize, optional :: m
         mwSize, optional :: n

         mwPointer :: mp, dataptr
         mwSize :: len1, len2

         len1 = size(array2d,1)
         if(present(m)) len1 = m
         len2 = size(array2d,2)
         if(present(n)) len2 = n
         
         mp = mxCreateNumericArray( 2_mws_, (/ len1, len2 /), &
            mxClassIDFromClassName('double'), 0_int4_ ) 
         dataptr = MATLAB_get_ptr(mp)
         call mxCopyReal8ToPtr(array2d, dataptr, len1*len2)

         fortran_to_matlab_array2d_real8 = mp
      end function fortran_to_matlab_array2d_real8

      mwPointer function fortran_to_matlab_scalar_complex8(scalar)
         complex(real4_), intent(in) :: scalar

         fortran_to_matlab_scalar_complex8 = &
            fortran_to_matlab_array1d_complex8( (/scalar/) )
      end function fortran_to_matlab_scalar_complex8
      mwPointer function fortran_to_matlab_array1d_complex8(array1d, n)
         complex(real4_), dimension(:), intent(in) :: array1d
         mwSize, optional :: n

         mwPointer :: mp, realptr, imagptr
         mwSize :: len

         len = size(array1d)
         if(present(n)) len = n
         
         mp = mxCreateNumericArray( 2_mws_, (/ len, 1_mws_ /), &
            mxClassIDFromClassName('single'), 1_int4_ ) 
         realptr = MATLAB_get_ptr(mp)
         imagptr = MATLAB_get_imag_ptr(mp)
         call MATLAB_copy_to_ptr(array1d, realptr, imagptr, len)

         fortran_to_matlab_array1d_complex8 = mp
      end function fortran_to_matlab_array1d_complex8
      mwPointer function fortran_to_matlab_array2d_complex8(array2d, m, n)
         complex(real4_), dimension(:,:), intent(in) :: array2d
         mwSize, optional :: m
         mwSize, optional :: n

         mwPointer :: mp, realptr, imagptr
         mwSize :: len1, len2

         len1 = size(array2d,1)
         if(present(m)) len1 = m
         len2 = size(array2d,2)
         if(present(n)) len2 = n
         
         mp = mxCreateNumericArray( 2_mws_, (/ len1, len2 /), &
            mxClassIDFromClassName('single'), 1_int4_ ) 
         realptr = MATLAB_get_ptr(mp)
         imagptr = MATLAB_get_imag_ptr(mp)
         call mxCopyComplex8ToPtr(array2d, realptr, imagptr, len1*len2)

         fortran_to_matlab_array2d_complex8 = mp
      end function fortran_to_matlab_array2d_complex8
      mwPointer function fortran_to_matlab_sparse_int4_complex8(m, n, ptr, &
            row, val)
         mwSize, intent(in) :: m, n
         integer(int4_), dimension(n+1), intent(in) :: ptr
         integer, dimension(ptr(n+1)-1), intent(in) :: row ! default integer
         complex(real4_), dimension(ptr(n+1)-1), intent(in) :: val

         mwPointer :: mp, ptrptr, rowptr, realptr, imagptr
         mwSize :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         complex(real8_), dimension(:), allocatable :: tempcomplex16
         integer :: st
         
         ! Setup a sparse matrix. This allocates ptr, row, val on matlab side
         nz = ptr(n+1)-1
         mp = mxCreateSparse(m, n, nz, 1_int4_)
         
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         realptr = MATLAB_get_ptr(mp)
         imagptr = MATLAB_get_imag_ptr(mp)

         ! Convert ptr from 1-based to 0-based then copy out
         allocate(tempmwi(n+1), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = ptr(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, ptrptr, n+1_mws_)
         deallocate(tempmwi)

         ! Convert row from 1-based to 0-based then copy out
         allocate(tempmwi(nz), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = row(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, rowptr, nz)
         deallocate(tempmwi)

         ! Copy out val. mxCreateSparse only creates double prec matrices
         allocate(tempcomplex16(nz), stat=st)
         if(st.ne.0) goto 100
         tempcomplex16(:) = val(:)
         call MATLAB_copy_to_ptr(tempcomplex16, realptr, imagptr, nz)

         fortran_to_matlab_sparse_int4_complex8 = mp

         return
         100 continue
         call MATLAB_context_error("hsl_matlab:AllocFail", "Allocation failed")
      end function fortran_to_matlab_sparse_int4_complex8
      mwPointer function fortran_to_matlab_sparse_int8_complex8(m, n, ptr, &
            row, val)
         mwSize, intent(in) :: m, n
         integer(int8_), dimension(n+1), intent(in) :: ptr
         integer, dimension(ptr(n+1)-1), intent(in) :: row ! default integer
         complex(real4_), dimension(ptr(n+1)-1), intent(in) :: val

         mwPointer :: mp, ptrptr, rowptr, realptr, imagptr
         mwSize :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         complex(real8_), dimension(:), allocatable :: tempcomplex16
         integer :: st
         
         ! Setup a sparse matrix. This allocates ptr, row, val on matlab side
         nz = ptr(n+1)-1
         mp = mxCreateSparse(m, n, nz, 1_int4_)
         
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         realptr = MATLAB_get_ptr(mp)
         imagptr = MATLAB_get_imag_ptr(mp)

         ! Convert ptr from 1-based to 0-based then copy out
         allocate(tempmwi(n+1), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = ptr(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, ptrptr, n+1_mws_)
         deallocate(tempmwi)

         ! Convert row from 1-based to 0-based then copy out
         allocate(tempmwi(nz), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = row(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, rowptr, nz)
         deallocate(tempmwi)

         ! Copy out val. mxCreateSparse only creates double prec matrices
         allocate(tempcomplex16(nz), stat=st)
         if(st.ne.0) goto 100
         tempcomplex16(:) = val(:)
         call MATLAB_copy_to_ptr(tempcomplex16, realptr, imagptr, nz)

         fortran_to_matlab_sparse_int8_complex8 = mp

         return
         100 continue
         call MATLAB_context_error("hsl_matlab:AllocFail", "Allocation failed")
      end function fortran_to_matlab_sparse_int8_complex8

      mwPointer function fortran_to_matlab_scalar_complex16(scalar)
         complex(real8_), intent(in) :: scalar

         fortran_to_matlab_scalar_complex16 = &
            fortran_to_matlab_array1d_complex16( (/scalar/) )
      end function fortran_to_matlab_scalar_complex16
      mwPointer function fortran_to_matlab_array1d_complex16(array1d, n)
         complex(real8_), dimension(:), intent(in) :: array1d
         mwSize, optional :: n

         mwPointer :: mp, realptr, imagptr
         mwSize :: len

         len = size(array1d)
         if(present(n)) len = n
         
         mp = mxCreateNumericArray( 2_mws_, (/ len, 1_mws_ /), &
            mxClassIDFromClassName('double'), 1_int4_ ) 
         realptr = MATLAB_get_ptr(mp)
         imagptr = MATLAB_get_imag_ptr(mp)
         call MATLAB_copy_to_ptr(array1d, realptr, imagptr, len)

         fortran_to_matlab_array1d_complex16 = mp
      end function fortran_to_matlab_array1d_complex16
      mwPointer function fortran_to_matlab_array2d_complex16(array2d, m, n)
         complex(real8_), dimension(:,:), intent(in) :: array2d
         mwSize, optional :: m
         mwSize, optional :: n

         mwPointer :: mp, realptr, imagptr
         mwSize :: len1, len2

         len1 = size(array2d,1)
         if(present(m)) len1 = m
         len2 = size(array2d,2)
         if(present(n)) len2 = n
         
         mp = mxCreateNumericArray( 2_mws_, (/ len1, len2 /), &
            mxClassIDFromClassName('double'), 1_int4_ ) 
         realptr = MATLAB_get_ptr(mp)
         imagptr = MATLAB_get_imag_ptr(mp)
         call mxCopyComplex16ToPtr(array2d, realptr, imagptr, len1*len2)

         fortran_to_matlab_array2d_complex16 = mp
      end function fortran_to_matlab_array2d_complex16
      mwPointer function fortran_to_matlab_sparse_int4_complex16(m, n, ptr, &
            row, val)
         mwSize, intent(in) :: m, n
         integer(int4_), dimension(n+1), intent(in) :: ptr
         integer, dimension(ptr(n+1)-1), intent(in) :: row ! default integer
         complex(real8_), dimension(ptr(n+1)-1), intent(in) :: val

         mwPointer :: mp, ptrptr, rowptr, realptr, imagptr
         mwSize :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         integer :: st
         
         ! Setup a sparse matrix. This allocates ptr, row, val on matlab side
         nz = ptr(n+1)-1
         mp = mxCreateSparse(m, n, nz, 1_int4_)
         
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         realptr = MATLAB_get_ptr(mp)
         imagptr = MATLAB_get_imag_ptr(mp)

         ! Convert ptr from 1-based to 0-based then copy out
         allocate(tempmwi(n+1), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = ptr(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, ptrptr, n+1_mws_)
         deallocate(tempmwi)

         ! Convert row from 1-based to 0-based then copy out
         allocate(tempmwi(nz), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = row(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, rowptr, nz)
         deallocate(tempmwi)

         ! Copy out val. mxCreateSparse only creates double prec matrices
         call MATLAB_copy_to_ptr(val, realptr, imagptr, nz)

         fortran_to_matlab_sparse_int4_complex16 = mp

         return
         100 continue
         call MATLAB_context_error("hsl_matlab:AllocFail", "Allocation failed")
      end function fortran_to_matlab_sparse_int4_complex16
      mwPointer function fortran_to_matlab_sparse_int8_complex16(m, n, ptr, &
            row, val)
         mwSize, intent(in) :: m, n
         integer(int8_), dimension(n+1), intent(in) :: ptr
         integer, dimension(ptr(n+1)-1), intent(in) :: row ! default integer
         complex(real8_), dimension(ptr(n+1)-1), intent(in) :: val

         mwPointer :: mp, ptrptr, rowptr, realptr, imagptr
         mwSize :: nz
         mwIndex, dimension(:), allocatable :: tempmwi
         integer :: st
         
         ! Setup a sparse matrix. This allocates ptr, row, val on matlab side
         nz = ptr(n+1)-1
         mp = mxCreateSparse(m, n, nz, 1_int4_)
         
         ptrptr = MATLAB_get_jc(mp)
         rowptr = MATLAB_get_ir(mp)
         realptr = MATLAB_get_ptr(mp)
         imagptr = MATLAB_get_imag_ptr(mp)

         ! Convert ptr from 1-based to 0-based then copy out
         allocate(tempmwi(n+1), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = ptr(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, ptrptr, n+1_mws_)
         deallocate(tempmwi)

         ! Convert row from 1-based to 0-based then copy out
         allocate(tempmwi(nz), stat=st)
         if(st.ne.0) goto 100
         tempmwi(:) = row(:) - 1
         call MATLAB_copy_to_ptr(tempmwi, rowptr, nz)
         deallocate(tempmwi)

         ! Copy out val. mxCreateSparse only creates double prec matrices
         call MATLAB_copy_to_ptr(val, realptr, imagptr, nz)

         fortran_to_matlab_sparse_int8_complex16 = mp

         return
         100 continue
         call MATLAB_context_error("hsl_matlab:AllocFail", "Allocation failed")
      end function fortran_to_matlab_sparse_int8_complex16

      mwPointer function MATLAB_create_structure(fieldnames)
         character(len=*), dimension(:), intent(in) :: fieldnames

         integer(int4_) :: nfields

         nfields = size(fieldnames)
         
         MATLAB_create_structure = MATLAB_create_structure_expert( &
            1_mws_, 1_mws_, nfields, fieldnames)
      end function MATLAB_create_structure

      subroutine matlab_set_field_scalar_int4(mp, field, scalar)
         mwPointer :: mp
         character(len=*), intent(in) :: field
         integer(int4_), intent(in) :: scalar

         mwPointer :: fptr

         fptr = fortran_to_matlab(scalar)
         call mxSetField( mp, 1_mwi_, field, fptr )
      end subroutine matlab_set_field_scalar_int4
      subroutine matlab_set_field_scalar_int8(mp, field, scalar)
         mwPointer :: mp
         character(len=*), intent(in) :: field
         integer(int8_), intent(in) :: scalar

         mwPointer :: fptr

         fptr = fortran_to_matlab(scalar)
         call mxSetField( mp, 1_mwi_, field, fptr )
      end subroutine matlab_set_field_scalar_int8
      subroutine matlab_set_field_scalar_real4(mp, field, scalar)
         mwPointer :: mp
         character(len=*), intent(in) :: field
         real(real4_), intent(in) :: scalar

         mwPointer :: fptr

         fptr = fortran_to_matlab(scalar)
         call mxSetField( mp, 1_mwi_, field, fptr )
      end subroutine matlab_set_field_scalar_real4
      subroutine matlab_set_field_scalar_real8(mp, field, scalar)
         mwPointer :: mp
         character(len=*), intent(in) :: field
         real(real8_), intent(in) :: scalar

         mwPointer :: fptr

         fptr = fortran_to_matlab(scalar)
         call mxSetField( mp, 1_mwi_, field, fptr )
      end subroutine matlab_set_field_scalar_real8
      subroutine matlab_set_field_scalar_complex8(mp, field, scalar)
         mwPointer :: mp
         character(len=*), intent(in) :: field
         complex(real4_), intent(in) :: scalar

         mwPointer :: fptr

         fptr = fortran_to_matlab(scalar)
         call mxSetField( mp, 1_mwi_, field, fptr )
      end subroutine matlab_set_field_scalar_complex8
      subroutine matlab_set_field_scalar_complex16(mp, field, scalar)
         mwPointer :: mp
         character(len=*), intent(in) :: field
         complex(real8_), intent(in) :: scalar

         mwPointer :: fptr

         fptr = fortran_to_matlab(scalar)
         call mxSetField( mp, 1_mwi_, field, fptr )
      end subroutine matlab_set_field_scalar_complex16
      subroutine matlab_set_field_character(mp, field, scalar)
         mwPointer :: mp
         character(len=*), intent(in) :: field
         character(len=*), intent(in) :: scalar

         mwPointer :: fptr

         fptr = mxCreateCharMatrixFromStrings( 1_mws_, trim(scalar) )
         call mxSetField(mp, 1_mwi_, field, fptr)
      end subroutine matlab_set_field_character
      subroutine matlab_set_field_logical(mp, field, scalar)
         mwPointer :: mp
         character(len=*), intent(in) :: field
         logical, intent(in) :: scalar

         mwPointer :: fptr
         integer :: sint

         sint = 0
         if(scalar) sint = 1

         fptr = fortran_to_matlab(sint)
         call mxSetField(mp, 1_mwi_, field, fptr)
      end subroutine matlab_set_field_logical

      subroutine matlab_get_field_scalar_int4(mp, field, scalar, argname, &
            isset, default)
         mwPointer :: mp
         character(len=*), intent(in) :: field
         integer(int4_), intent(out) :: scalar
         character(len=*), intent(in) :: argname
         logical, optional, intent(out) :: isset
         integer(int4_), optional, intent(in) :: default

         mwPointer :: fptr
         character(len=200) :: errmsg

         if(.not. MATLAB_is_structure(mp)) then
            write(errmsg, "(3a)") "Error in argument ", argname, &
               ". Expected structure."
            call MATLAB_context_error('hsl_matlab:ExpectsStructure', errmsg)
         endif

         fptr = mxGetField(mp, 1_mwi_, field)
         if(fptr.eq.0) then
            if(present(isset)) isset = .false.
            if(present(default)) scalar = default
            return
         elseif(present(isset)) then
            isset = .true.
         endif

         write(errmsg, "(3a)") trim(argname), ".", trim(field)
         call matlab_to_fortran(fptr, scalar, trim(errmsg))
      end subroutine matlab_get_field_scalar_int4
      subroutine matlab_get_field_scalar_int8(mp, field, scalar, argname, &
            isset, default)
         mwPointer :: mp
         character(len=*), intent(in) :: field
         integer(int8_), intent(out) :: scalar
         character(len=*), intent(in) :: argname
         logical, optional, intent(out) :: isset
         integer(int8_), optional, intent(in) :: default

         mwPointer :: fptr
         character(len=200) :: errmsg

         if(.not. MATLAB_is_structure(mp)) then
            write(errmsg, "(3a)") "Error in argument ", argname, &
               ". Expected structure."
            call MATLAB_context_error('hsl_matlab:ExpectsStructure', errmsg)
         endif

         fptr = mxGetField(mp, 1_mwi_, field)
         if(fptr.eq.0) then
            if(present(isset)) isset = .false.
            if(present(default)) scalar = default
            return
         elseif(present(isset)) then
            isset = .true.
         endif

         write(errmsg, "(3a)") trim(argname), ".", trim(field)
         call matlab_to_fortran(fptr, scalar, trim(errmsg))
      end subroutine matlab_get_field_scalar_int8
      subroutine matlab_get_field_scalar_real4(mp, field, scalar, argname, &
            isset, default)
         mwPointer :: mp
         character(len=*), intent(in) :: field
         real(real4_), intent(out) :: scalar
         character(len=*), intent(in) :: argname
         logical, optional, intent(out) :: isset
         real(real4_), optional, intent(in) :: default

         mwPointer :: fptr
         character(len=200) :: errmsg

         if(.not. MATLAB_is_structure(mp)) then
            write(errmsg, "(3a)") "Error in argument ", argname, &
               ". Expected structure."
            call MATLAB_context_error('hsl_matlab:ExpectsStructure', errmsg)
         endif

         fptr = mxGetField(mp, 1_mwi_, field)
         if(fptr.eq.0) then
            if(present(isset)) isset = .false.
            if(present(default)) scalar = default
            return
         elseif(present(isset)) then
            isset = .true.
         endif

         write(errmsg, "(3a)") trim(argname), ".", trim(field)
         call matlab_to_fortran(fptr, scalar, trim(errmsg))
      end subroutine matlab_get_field_scalar_real4
      subroutine matlab_get_field_scalar_real8(mp, field, scalar, argname, &
            isset, default)
         mwPointer :: mp
         character(len=*), intent(in) :: field
         real(real8_), intent(out) :: scalar
         character(len=*), intent(in) :: argname
         logical, optional, intent(out) :: isset
         real(real8_), optional, intent(in) :: default

         mwPointer :: fptr
         character(len=200) :: errmsg

         if(.not. MATLAB_is_structure(mp)) then
            write(errmsg, "(3a)") "Error in argument ", argname, &
               ". Expected structure."
            call MATLAB_context_error('hsl_matlab:ExpectsStructure', errmsg)
         endif

         fptr = mxGetField(mp, 1_mwi_, field)
         if(fptr.eq.0) then
            if(present(isset)) isset = .false.
            if(present(default)) scalar = default
            return
         elseif(present(isset)) then
            isset = .true.
         endif

         write(errmsg, "(3a)") trim(argname), ".", trim(field)
         call matlab_to_fortran(fptr, scalar, trim(errmsg))
      end subroutine matlab_get_field_scalar_real8
      subroutine matlab_get_field_scalar_complex8(mp, field, scalar, argname, &
            isset, default)
         mwPointer :: mp
         character(len=*), intent(in) :: field
         complex(real4_), intent(out) :: scalar
         character(len=*), intent(in) :: argname
         logical, optional, intent(out) :: isset
         complex(real4_), optional, intent(in) :: default

         mwPointer :: fptr
         character(len=200) :: errmsg

         if(.not. MATLAB_is_structure(mp)) then
            write(errmsg, "(3a)") "Error in argument ", argname, &
               ". Expected structure."
            call MATLAB_context_error('hsl_matlab:ExpectsStructure', errmsg)
         endif

         fptr = mxGetField(mp, 1_mwi_, field)
         if(fptr.eq.0) then
            if(present(isset)) isset = .false.
            if(present(default)) scalar = default
            return
         elseif(present(isset)) then
            isset = .true.
         endif

         write(errmsg, "(3a)") trim(argname), ".", trim(field)
         call matlab_to_fortran(fptr, scalar, trim(errmsg))
      end subroutine matlab_get_field_scalar_complex8
      subroutine matlab_get_field_scalar_complex16(mp, field, scalar, argname, &
            isset, default)
         mwPointer :: mp
         character(len=*), intent(in) :: field
         complex(real8_), intent(out) :: scalar
         character(len=*), intent(in) :: argname
         logical, optional, intent(out) :: isset
         complex(real8_), optional, intent(in) :: default

         mwPointer :: fptr
         character(len=200) :: errmsg

         if(.not. MATLAB_is_structure(mp)) then
            write(errmsg, "(3a)") "Error in argument ", argname, &
               ". Expected structure."
            call MATLAB_context_error('hsl_matlab:ExpectsStructure', errmsg)
         endif

         fptr = mxGetField(mp, 1_mwi_, field)
         if(fptr.eq.0) then
            if(present(isset)) isset = .false.
            if(present(default)) scalar = default
            return
         elseif(present(isset)) then
            isset = .true.
         endif

         write(errmsg, "(3a)") trim(argname), ".", trim(field)
         call matlab_to_fortran(fptr, scalar, errmsg)
      end subroutine matlab_get_field_scalar_complex16
      subroutine matlab_get_field_character(mp, field, string, argname, &
            isset, default)
         mwPointer :: mp
         character(len=*), intent(in) :: field
         character(len=*), intent(out) :: string
         character(len=*), intent(in) :: argname
         logical, optional, intent(out) :: isset
         character(len=*), optional, intent(in) :: default

         mwPointer :: fptr
         character(len=200) :: errmsg
         mwSize :: strlen

         if(.not. MATLAB_is_structure(mp)) then
            write(errmsg, "(3a)") "Error in argument ", argname, &
               ". Expected structure."
            call MATLAB_context_error('hsl_matlab:ExpectsStructure', errmsg)
         endif

         fptr = mxGetField(mp, 1_mwi_, field)
         if(fptr.eq.0) then
            if(present(isset)) isset = .false.
            if(present(default)) string = default
            return
         elseif(present(isset)) then
            isset = .true.
         endif

         strlen = len(string)
         call MATLAB_get_string(fptr, string, strlen)
      end subroutine matlab_get_field_character


    END MODULE HSL_MATLAB
!-*-*-*-*-*-*- E N D  o f  H S L _ M A T L A B   M O D U L E -*-*-*-*-*-
