! COPYRIGHT (c) 2010 Science and Technology Facilities Council
! Original date 13 December 2010, Version 1.0.0
!
! Written by: Jonathan Hogg and Jennifer Scott
!
! Version 1.3.0
! See ChangeLog for history
!

!
! To convert from single:
! * Change wp
! * Change _complex [mc34_single]
! * Change BLAS calls: cgemv, cgemm, ctrsm, ctrsv, ccopy, cswap, caxpy
! * Change HSL calls: mc77i, mc77a, mc64w
! * Change control%u default to 0.01
!

! Note: in complex case, calls to factorize
! require an additional argument matrix_type.

module hsl_MA86_complex
!$ use omp_lib
   use hsl_mc78_integer
   use hsl_mc34_single
   implicit none

   private
   public :: ma86_keep, ma86_control, ma86_info
   public :: ma86_analyse, ma86_factor, ma86_factor_solve, ma86_solve, &
             ma86_finalise
   ! Following line contains symbols that should be used by HSL interfaces only
   public :: ma86_get_n__

!*************************************************
   ! Parameters (all private)
   ! Data kinds
   integer, parameter :: wp   = kind(0e0)
   integer, parameter :: long = selected_int_kind(18)

   ! Numerical constants
   complex(wp), parameter :: cone  = (1.0_wp,0.0_wp)
   complex(wp), parameter :: czero = (0.0_wp,0.0_wp)
   real(wp), parameter :: rone  = 1.0_wp
   real(wp), parameter :: rzero = 0.0_wp

   ! Default values
   integer, parameter :: nemin_default = 32
      ! node amalgamation parameter
   integer, parameter :: nb_default = 256
      ! Block size with dense kernel
   integer, parameter :: nbi_default = 16 ! Default inner block size.
   integer, parameter :: pool_default = 25000
      ! size of task pool

   ! Symbolic constants
   ! These flag the different tasks within factor and solve
   integer, parameter :: TASK_DONE             = -1
   integer, parameter :: TASK_NONE             = 0
   integer, parameter :: TASK_FACTORIZE_COLUMN = 2
   integer, parameter :: TASK_UPDATE_INTERNAL  = 3
   integer, parameter :: TASK_UPDATE_BETWEEN   = 4
   integer, parameter :: TASK_SLV_FSLV         = 6
      ! Fwds solve on diag block
   integer, parameter :: TASK_SLV_BSLV         = 8
      ! Bwds solve on diag block

   ! Types of solve job
   integer, parameter :: SOLVE_JOB_ALL         = 0
   integer, parameter :: SOLVE_JOB_FWD         = 1
   integer, parameter :: SOLVE_JOB_D           = 2
   integer, parameter :: SOLVE_JOB_BWD         = 3
   integer, parameter :: SOLVE_JOB_D_AND_BWD   = 4

   ! How processors share cache                    Example
   integer, parameter :: CACHE_COMPACT       = 1
      ! [0,1], [2,3], [4,5], [6,7]
   integer, parameter :: CACHE_SCATTER       = 2
      ! [0,4]. [1,5], [2,6], [3,7]
   integer, parameter :: CACHE_IDENTITY      = 3
      ! 0, 1, 2, 3, 4, 5, 6, 7

   ! Error flags
   integer, parameter :: MA86_SUCCESS               = 0
   integer, parameter :: MA86_ERROR_ALLOCATION      = -1
   integer, parameter :: MA86_ERROR_ORDER           = -2
   integer, parameter :: MA86_ERROR_SINGULAR        = -3
   integer, parameter :: MA86_ERROR_X_SIZE          = -4
   integer, parameter :: MA86_ERROR_INFINITY        = -5
   integer, parameter :: MA86_ERROR_JOB_OOR         = -6
   integer, parameter :: MA86_ERROR_STATIC_SMALL    = -7
   integer, parameter :: MA86_ERROR_MATRIX_TYPE     = -8
   !integer, parameter :: MA86_ERROR_NO_SCALE        = -9
   integer, parameter :: MA86_ERROR_UNKNOWN         = -99

   ! warning flags
   integer, parameter :: MA86_WARNING_POOL_SMALL       = 1
   integer, parameter :: MA86_WARNING_SINGULAR         = 2
   integer, parameter :: MA86_WARNING_POOL_SING        = 3

   ! matrix types (Hermitian and complex symmetric)
   integer, parameter :: HSL_MATRIX_CPLX_HERM_INDEF = -4
   integer, parameter :: HSL_MATRIX_CPLX_SYM = -5

   !*************************************************

   interface MA86_analyse
      module procedure MA86_analyse_complex
   end interface

   interface MA86_factor
      module procedure MA86_factor_complex
   end interface

   interface MA86_factor_solve
      module procedure MA86_factor_solve_one_complex, &
                       MA86_factor_solve_mult_complex
   end interface

   interface MA86_solve
      module procedure MA86_solve_one_complex, &
         MA86_solve_mult_complex
   end interface

   interface MA86_finalise 
      module procedure MA86_finalise_complex
   end interface

   interface MA86_get_n__
      module procedure ma86_get_n_complex
   end interface MA86_get_n__

   !*************************************************

   ! Data type for storing information for each block (BLK)
   ! The blocks are numbered 1,2,..., keep%final_blk
   type block_type
      ! Static info, which is set in ma86_analayse
      integer :: bcol            ! block column that blk belongs to
      integer :: blkm            ! height of block (number of rows in blk)
      integer :: blkn            ! width of block (number of columns in blk)
      integer(long) :: dblk      ! id of the block on the diagonal within the 
         ! block column to which blk belongs
      integer :: dep_initial     ! initial dependency count for block,
         ! In indef case, countdown occurs in bcol on a block column basis
      integer(long) :: id        ! The block identitifier (ie, its number blk)
      integer(long) :: last_blk  ! id of the last block within the
         ! block column to which blk belongs
      integer :: node            ! node to which blk belongs
      integer :: sa              ! posn of the first entry of the
         ! block blk within the array that holds the block column of L
         ! that blk belongs to.

      ! Non-static info
      logical :: touched ! is this the first time block is touched
      integer :: sa_new          ! posn of the first entry of the
         ! block blk within the array that holds the block column of L
         ! that blk belongs to, after delays have been allowed for. 
         ! This is computed during factorize. 
!$    integer(omp_lock_kind) :: lock   ! Lock for altering dep
!$    integer(omp_lock_kind) :: alock  ! Lock for altering values in keep%lfact 
         ! for this block.
         ! Note: locks initialised in ma86_analyse and destroyed
         !       in ma86_finalise
   end type block_type

   !*************************************************

   ! Derived type for holding data for each node.
   ! This information is set up by ma86_analyse once the assembly tree
   ! has been constructed.
   type node_type
      integer(long) :: blk_sa ! identifier of the first block in node
      integer(long) :: blk_en ! identifier of the last block in node

      integer :: nb ! Block size for nodal matrix
         ! If number of cols nc in nodal matrix is less than control%nb but 
         ! number of rows is large, the block size for the node is taken as 
         ! control%nb**2/nc, rounded up to a multiple of 8. The aim is for
         ! the numbers of entries in the blocks to be similar to those in the 
         ! normal case. 

      integer :: sa ! index (in pivotal order) of the first column of the node
      integer :: en ! index (in pivotal order) of the last column of the node

      integer, allocatable :: index(:) ! holds the permuted variable
         ! list for node. They are sorted into increasing order.
         ! index is set up by ma86_analyse

      integer :: nchild ! number of children node has in assembly tree
      integer, allocatable :: child(:) ! holds children of node
      integer :: num_delay = 0 ! number of delayed eliminations to
         ! pass to parent
      integer :: parent ! Parent of node in assembly tree
      integer :: least_desc ! Least descendant in assembly tree

      ! next pointer in linked list of delay sources
      integer :: delay_next
   end type node_type

   !*************************************************

   ! This type stores overall values that this thread has encountered. It
   ! must be combined across all threads to be meaningful.
   type thread_info
      integer :: num_delay     = 0 ! number of delayed pivots
      integer(long) :: num_flops = 0 ! number of floating point operations
      integer(long) :: num_factor = 0 ! number of entries in factors
      integer :: num_neg       = 0 ! number of negative pivots in the real
                                   ! or complex Hermitian case.
      integer :: num_nothresh  = 0 ! number of pivots not satisfying
         ! pivot threshold test with control%u
      integer :: num_perturbed = 0 ! number of perturbed pivots
      integer :: num_two       = 0 ! number of 2x2 pivots
      integer :: num_zero_pivots  = 0 ! number of zero pivots
      real(wp) :: usmall       = -rone ! Set to zero if num_perturbed > 0.
         ! Otherwise, if q < p, it holds the value of cntl%umin that 
         ! would have led to a greater value of q and if q = p, it holds 
         ! the smallest relative pivot value of the chosen pivots.

      real(wp) :: detlog       = rzero ! logarithm of abs value of det A
      integer :: detsign       = 1 ! in the real or complex Hermitian case,
         ! holds sign of determinant or 0 if A is singular
      complex(wp) :: detarg = cone ! complex symmetric case, holds determinant 
         ! of A divided by its absolute value or one if A singular. 
   end type thread_info


   !*************************************************

   ! Data type that represents a single block column in L and corresponding
   ! entries in the indefinite case (allocated by ma86_analyse)
   type lfactor
      integer :: blkn_new ! number of columns in lcol_new
      integer(long) :: dblk      ! identifier of diagonal block in block col.
      integer :: dep             ! dependency count
      integer :: local ! local index of block column within node.
      integer :: nelim ! number of eliminations performed in block column.
      integer :: nrow ! number of rows in block column at start of factorize
         ! ie number of rows in lcol. Number of rows in lcol_new is
         ! nrow + num_delay   (num_delay = blkn_new - blkn, where blkn is no.
         ! cols in lcol and can be found from blocks using dblk)
      integer :: col ! start position of bcol in rhs vector for solve

      ! Linked list of columns to draw delays from
      integer :: delay_head

      integer, allocatable :: index_new(:) ! holds index list within
         ! factorize (will include delays)
!$    integer(omp_lock_kind) :: lock   ! lock so only one thread at a time
         ! can alter the block column
      complex(wp), dimension(:), allocatable :: lcol ! holds block column
      complex(wp) , dimension(:), allocatable :: d ! holds block of D.
   end type lfactor

   !*************************************************

   ! Data type for storing mapping from user's matrix values into
   ! block columns of L
   type lmap_type
      integer(long) :: len_map ! length of map
      integer(long), allocatable :: map(:,:) ! holds map from user's val
         ! array into lfact(:)%lcol values as follows:
         ! lcol( map(1,i) ) += val( map(2,i) )     i = 1:lmap
         ! map is set at end of analyse phase using subroutines make_map
         ! and lcol_map, and is used in factorization phase by blk_col_add_a
   end type lmap_type

   !*************************************************

   ! Data type that contains counts and locks for the solve
   ! (one per a block column)
   type slv_count_type
      integer :: dep
      integer :: dblk
!$    integer(kind=omp_lock_kind) :: lock
   end type slv_count_type

   !*************************************************

   ! Data type for user controls
   type MA86_control

      logical :: action          = .true. ! Do we keep
         ! going even if matrix is singular (abort if .false.)
      integer :: diagnostics_level = 0      ! Controls diagnostic printing.
         ! Possible values are:
         !  < 0: no printing.
         !    0: error and warning messages only.
         !    1: as 0 plus basic diagnostic printing.
         !    2: as 1 plus some more detailed diagnostic messages.
         !    3: as 2 plus all entries of user-supplied arrays.
      integer :: nb    = nb_default ! Controls the size of the
         ! blocks used within each node (used to set nb within node_type)
      integer :: nbi   = nbi_default ! Inner block size for use with ma64
      integer :: nemin = nemin_default    
         ! Node amalgamation parameter. A child node is merged with its parent 
         ! if they both involve fewer than nemin eliminations.
      integer :: pool_size       = pool_default ! Size of task pool arrays
      real(wp) :: small          = 1e-20 ! Pivots less than small are
         ! treated as zero
      real(wp) :: static         = rzero ! Control static pivoting
      real(wp) :: u              = 0.1 ! Pivot tolerance
      real(wp) :: umin           = rone ! Minimum pivot tolerance
      integer :: unit_diagnostics = 6    ! unit for diagnostic messages
         ! Printing is suppressed if unit_diagnostics  <  0.
      integer :: unit_error       = 6    ! unit for error messages
         ! Printing is suppressed if unit_error  <  0.
      integer :: unit_warning     = 6    ! unit for warning messages
         ! Printing is suppressed if unit_warning  <  0.
      integer :: scaling          = 0    ! scaling routine to use
         ! 0 = none or user defined (user iff scale is present)
         ! 1 = mc64
         ! 2 = mc77

      !!!! Undocumented
   !**   integer :: time_out        = -1     ! If >= 0 some internal timings
   !**      are printed on unit time_out. For HSL 2011 these are commented
   !**      using comments !** so easy to search and uncomment
   !%%%  integer :: unit_log        = -1     ! For profiling log output
   !%%%     commented out for HSL 2011 using !%%%
   !%%%  integer :: log_level       = 1      ! Level of profiling information
   !!! Note: commenting out use of time_out and unit_log means
   !%%%     commented out for HSL 2011 using !%%%

   !!! As a result, some subroutines have unused dummy arguments that
   !!! give warnings at compile time. We have not removed them
   !!! since uncommenting the above controls would then be more tedious.

      integer :: cache_tq_sz     = 100    ! Size of local task stack
      integer :: cache_layout    = CACHE_COMPACT ! Proc <=> cache mapping
      integer :: cache_cores     = 2      ! Number of cores per cache
      integer :: min_width_blas  = 8      ! Minimum width of source block
         ! before we use an indirect update_between
    
   end type MA86_control

   !*************************************************

   ! data type for returning information to user.
   type MA86_info 
      real(wp) :: detlog = rzero        ! Holds logarithm of abs det A (or 0)
      integer :: detsign = 0            ! Holds sign of determinant (+/-1 or 0)
      integer :: flag = 0               ! Error return flag (0 on success)
      integer :: matrix_rank = 0        ! Rank of matrix
      integer :: maxdepth = 0           ! Maximum depth of the tree.
      integer :: num_delay = 0          ! Number of delayed pivots
      integer(long) :: num_factor = 0_long ! Number of entries in the factor.
      integer(long) :: num_flops = 0_long  ! Number of flops for factor.
      integer :: num_neg = 0            ! Number of negative pivots
      integer :: num_nodes = 0          ! Number of nodes
      integer :: num_nothresh = 0       ! Number of pivots not satisfying u
      integer :: num_perturbed = 0      ! Number of perturbed pivots
      integer :: num_two = 0            ! Number of 2x2 pivots
      integer :: pool_size = pool_default  ! Maximum size of task pool used
      integer :: stat = 0               ! STAT value on error return -1.
      real(wp) :: usmall = rzero        ! smallest threshold parameter used
   end type MA86_info

   !*************************************************

   ! Data type for communication between threads and routines
   type ma86_keep
      private
      type(block_type), dimension(:), allocatable :: blocks ! block info
      integer, dimension(:), allocatable :: flag_array ! allocated to
         ! have size equal to the number of threads. For each thread, holds
         ! error flag
      integer(long) :: final_blk = 0 ! Number of blocks. Used for destroying
         ! locks in finalise
      type(ma86_info) :: info ! Holds copy of info
      integer :: maxm ! holds largest block row dimension
      integer :: maxn ! holds largest block col dimension
      integer :: matrix_type ! copy of matrix_type input to factorize
      integer :: n  ! Order of the system.
      type(node_type), dimension(:), allocatable :: nodes ! nodal info
      integer :: nbcol = 0 ! number of block columns in L
      type(lfactor), dimension(:), allocatable :: lfact
         ! holds block cols of L
      type(lmap_type), dimension(:), allocatable :: lmap
         ! holds mapping from matrix values into lfact
      real(wp), dimension(:), allocatable :: scaling
   end type ma86_keep

   !*************************************************

   ! Data type for a task
   type dagtask
      integer :: task_type    ! One of TASK_FACTORIZE_COLUMN, ...
      integer(long) :: dest   ! id of the target (destination) block
      integer(long) :: src1   ! 
         ! if task_type = TASK_UPDATE_INTERNAL, src1 holds the id of the first 
         ! source block
         ! if task_type = TASK_UPDATE_BETWEEN, src1 holds the id of a block 
         ! in the block column of the source node that is used
         ! in updating dest.
      integer(long) :: src2   
         ! if task_type = TASK_UPDATE_INTERNAL, src2 holds the id of the second 
         ! source block
         ! (src1 and src2 are blocks belonging to the same block column
         ! of the source node with src1 .le. src2)
         ! src2 is not used by the other tasks
      integer :: csrc(2)
      integer :: rsrc(2)
         ! for an UPDATE_BETWEEN task, we need to hold some additional
         ! information, which locates the source blocks rsrc and csrc
         ! within the source block col.
         ! This info. is set up the subroutine add_between_updates
   end type dagtask

   !*************************************************

   ! Data type for storing tasks we need to do.
   ! Task pool is held as a collection of 4 stacks for different priorities 
   ! using the same space. 
   ! Each stack is a linked list with its head given by an element of prihead.
   ! There are also local stacks for each cache. 

   type taskstack
      integer :: max_pool_size = 0 ! max. number of tasks that are in
         ! the task pool at any one time during the factorization.
      logical :: abort = .false.   ! true if we have aborted
      integer :: active ! Number of active threads 
         ! (number of tasks in execution)
      type(dagtask), dimension(:,:), allocatable :: ctasks ! local task stacks.
         ! allocated to have size (control%cache_tq_sz, ncache), where
         ! ncache is number of caches
      integer, dimension(:), allocatable :: cheads   ! Heads for local stacks.
         ! allocated to have size equal to number of caches
!$    integer(omp_lock_kind), dimension(:), allocatable :: clocks 
         ! Locks for local stacks.
      integer :: freehead  ! Holds the head of linked list of
         ! entries in the task pool that are free
!$    integer(omp_lock_kind) :: lock   ! lock so only one thread at a time
         ! can read/alter the task pool 
      integer :: lowest_priority_value = huge(0) ! 
         ! lowest priority value of the tasks in the pool.
         ! The priority value for each of the different types of task is
         !  1. factor             Highest priority 
         !  2. solve 
         !  3. update_internal 
         !  4. update_between     Lowest priority
      integer, dimension(:), allocatable :: next  ! next task in linked list.
         ! allocated to have size pool_size. Reallocated if initial setting
         ! for pool_size found to be too small.
      integer :: pool_size   ! sizes of task pool arrays next and tasks. 
         ! Initialised to control%pool_size
      integer :: prihead(4)  ! Holds the heads of the linked lists for tasks
         ! with priority values 1,2,3,4. 
      type(dagtask), dimension(:), allocatable :: tasks ! Holds tasks.
         ! allocated to have size pool_size. Reallocated if initial setting
         ! for pool_size found to be too small.
      integer :: total       ! Total number of tasks in pool
!**   real, dimension(:), allocatable :: waiting  ! Allocated to have size
!**      ! equal to the number of threads. Used to hold times the threads spent 
!**      ! waiting if control%time_out >= 0

   end type taskstack

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Analyse phase.
! The user inputs the pivot order and lower
! triangular parts of A. Structure is expanded.
! Supervariables are computed
! and then the assembly tree is constructed and the data structures
! required by the factorization are set up.
! There is no checking of the user's data.
!
subroutine MA86_analyse_complex(n, ptr, row, order, keep, control, info)
   integer, intent(in) :: n ! order of A
   integer, intent(in) :: row(:) ! row indices of lower triangular part
   integer, intent(in) :: ptr(:) ! col pointers for lower triangular part
   integer, intent(inout), dimension(:) :: order
      ! order(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by MA86_factor.
   ! For details of keep, control, info : see derived type descriptions
   type(MA86_keep), intent(out) :: keep
   type(MA86_control), intent(in) :: control
   type(MA86_info), intent(inout) :: info

   ! Local arrays
   integer, allocatable :: amap(:) ! map from user a to reordered a
   integer, allocatable :: aptr(:) ! column pointers of expanded matrix
   integer, allocatable :: arow(:) ! row pointers of expanded matrix
   integer, allocatable :: iw(:) ! work array
   integer, allocatable :: map(:) ! Allocated to have size n.
      ! used in computing dependency counts. For each row k in 
      ! j-th block column of a node, map(k1) is set to j
   integer, allocatable :: perm(:) ! inverse permutation.
      ! perm(i) holds the variable that is i-th in the pivot sequence.
      ! Also used for checking user-supplied permutation.

   ! Local scalars.
   integer :: a_nb ! block size of anode
   integer :: anode ! ancestoral node of snode in tree
   integer(long) :: blk ! temporary variable for holding block number
   integer :: blkn ! number of columns within block column
   integer :: cptr ! pointer into rows of snode that correspond
      ! to columns of an ancestor anode
   integer :: cb ! index of block column within node
   integer :: col_used ! used in computing number of cols in block col.
   integer :: ci ! do loop variable. current block col.
   integer(long) :: dblk ! diagonal block within block column
   integer :: en ! holds keep%nodes(snode)%en
   integer :: i ! temporary variable
   integer :: j ! temporary variable
   integer :: jj ! temporary variable
   integer :: jb ! block index in anode
   integer :: k
   integer :: k1 ! temporary variable
   integer :: l_nb ! set to block size of snode (keep%nodes(snode)%nb)
   integer :: mp  ! copy of control%unit_diagnostics
   integer :: ne ! set to ptr(n+1) - 1
   integer :: nemin ! min. number of eliminations (see control%nemin)
   integer :: node ! a node in tree
   integer :: num_nodes ! number of nodes in tree
   integer :: numcol ! number of cols in node (en-sa+1)
   integer :: numrow ! number of rows in a block column
   integer :: row_used ! used in computing number of rows in a block.
   integer :: sa ! holds keep%nodes(snode)%sa
   integer :: size_anode ! size(keep%nodes(anode)%index)
   integer :: st ! stat parameter
   integer :: sz ! number of blocks in a block column of node
   integer :: swidth ! number of block columns in node
!**integer :: t_start, t_end, t_rate

   type(mc78_control) :: control78
   integer :: par
   integer :: info78
   integer, dimension(:), allocatable :: sptr, sparent, rlist
   integer(long), dimension(:), allocatable :: rptr

   ! Possible error returns:
   !  MA86_ERROR_ALLOCATION   Allocation error
   !  MA86_ERROR_ORDER        Error in order

   ! initialise
   info%flag = 0
   info%num_factor = 0_long
   info%num_flops = 0_long
   info%num_nodes = 0
   info%maxdepth = 0
   info%stat = 0

   keep%n = n
   ne = ptr(n+1) - 1

   ! Perform appropriate printing
   mp = control%unit_diagnostics
   if (control%diagnostics_level>=1 .and. mp>=0) then
      write (mp,'(/a)') ' On entry to MA86_analyse:'
      write (mp,'(a,i15)') ' control%diagnostics_level =  ', &
         control%diagnostics_level
      write (mp,'(a,i15)') ' control%unit_diagnostics  =  ',mp
      write (mp,'(a,i15)') ' control%unit_error        =  ',control%unit_error
      write (mp,'(a,i15)') ' control%unit_warning      =  ',control%unit_warning
      write (mp,'(a,i15)') ' control%nemin             =  ',control%nemin
      write (mp,'(a,i15)') ' control%nb                =  ',control%nb
      write (mp,'(a,i15)') ' n                         =  ',n
   end if

   if (control%diagnostics_level > 2 .and. mp>=0) then

      write (mp,'(a)') ' ptr = '
      write (mp,'(5i15)') ptr(1:n+1)

      write (mp,'(a)') ' row = '
      write (mp,'(5i15)') row(1:ne)

      ! Print out pivot order.
      write (mp,'(a)') ' User-supplied elimination order :'
      i = min(size(order),n)
      write (mp,'(5i15)') order(1:i)

   else if (control%diagnostics_level==2 .and. mp>=0) then

      write (mp,'(a)') ' ptr(1:min(5,n+1)) = '
      write (mp,'(5i15)') ptr(1:min(5,n+1))

      write (mp,'(a)') ' row(1:min(5,ne)) =  '
      write (mp,'(5i15)') row(1:min(5,ne))

      i = min(10,n)
      i = min(size(order),i)
      write (mp,'(a,2(/5i12))')  &
         ' User-supplied elimination order :', order(1:i)
      if (i < n) write (mp,'(a)') '  . . . . . .'

   end if

   ! immediate return if n = 0
   if (n == 0) return

   ! expand the matrix

   ! allocate space for expanded matrix (held in aptr,arow)
   allocate (arow(2*ne),aptr(n+3),iw(n+1),stat=st)
   if (st /= 0) go to 490

   arow(1:ne) = row(1:ne)
   aptr(1:n+1) = ptr(1:n+1)
   call mc34_expand(n, arow, aptr, iw)

   deallocate(iw,stat=st)

   nemin = control%nemin
   ! Check nemin (a node is merged with its parent if both involve
   ! fewer than nemin eliminations). If out of range, use the default
   if (nemin < 1) nemin = nemin_default

   ! Check the user-supplied array order and set the inverse in perm.
   if (size(order).lt.n) then
      info%flag = MA86_ERROR_ORDER
      call MA86_print_flag(info%flag, control, context='MA86_analyse')
      go to 500
   end if

   deallocate (perm,stat=st)
   allocate (perm(n),stat=st)
   if (st /= 0) go to 490
   perm(:) = 0
   k1 = 0
   do i = 1, n
      jj = order(i)
      if (jj < 1 .or. jj > n) exit
      if (perm(jj) /= 0) exit ! Duplicate found
      perm(jj) = i
   end do
   if (i-1 /= n) then
      info%flag = MA86_ERROR_ORDER
      call MA86_print_flag(info%flag, control, context='MA86_analyse')
      go to 500
   end if

   control78%nemin = nemin
   control78%sort = .true.
   control78%lopt = .true.
   call mc78_analyse(n, aptr, arow, order, num_nodes, &
      sptr, sparent, rptr, rlist, control78, info78, nfact=info%num_factor, &
      nflops=info%num_flops, stat=st)
   select case (info78)
   case(0)
      ! success, do nothing
   case(-1)
      ! allocation error
      goto 490
   case(1)
      ! symbolically singular
      if(control%action) then
         ! issue warning and continue
         info%flag = MA86_WARNING_SINGULAR
      else
         ! abort with error
         info%flag = MA86_ERROR_SINGULAR
         goto 500
      endif
   case default
      ! panic
      info%flag = MA86_ERROR_UNKNOWN
      goto 500
   end select

   info%num_nodes = num_nodes
   !**************************************
   ! Set up nodal data structures
   ! For each node, hold data in keep%nodes(node) 
   deallocate(keep%nodes, stat=st)

   allocate(keep%nodes(-1:num_nodes),stat=st)
   if (st /= 0) go to 490

   keep%nodes(0)%blk_en = 0
   keep%nodes(1)%blk_sa = 1
   keep%nodes(1)%sa = 1

   ! loop over nodes
   keep%nodes(:)%nchild = 0
   do node = 1, num_nodes
      keep%nodes(node)%sa = sptr(node)
      keep%nodes(node)%en = sptr(node+1)-1

      par = sparent(node)
      keep%nodes(node)%parent = par
      if(par .le. num_nodes) then
         keep%nodes(par)%nchild = keep%nodes(par)%nchild + 1
      else
         keep%nodes(node)%parent = -1
      endif

      ! determine and record the block size for node
      ! note we are careful in case l_nb**2 overflows (in fact 1+l_nb must
      ! not overflow at the end), and limit the answer to huge(l_nb)/2
      l_nb = control%nb
      if (l_nb < 1) l_nb = nb_default
      l_nb = min(huge(l_nb)/2_long, &
         (l_nb**2_long) / min(sptr(node+1)-sptr(node), l_nb) )
      l_nb = (l_nb-1) / 8 + 1
      l_nb = 8 * l_nb
      keep%nodes(node)%nb = l_nb

      ! Copy row list into keep%nodes
      allocate(keep%nodes(node)%index(rptr(node+1)-rptr(node)),stat=st)
      if (st /= 0) go to 490
      j = 1
      do i = rptr(node), rptr(node+1)-1
         keep%nodes(node)%index(j) = rlist(i)
         j = j + 1
      end do

      ! Allocate space to store child nodes
      allocate(keep%nodes(node)%child(keep%nodes(node)%nchild), stat=st)
      if(st.ne.0) goto 490

      ! Calculate number j of blocks in node and set
      ! keep%nodes(node)%blk_en
      sz = (rptr(node+1)-rptr(node) - 1) / l_nb + 1
      j = 0
      do i = keep%nodes(node)%sa, keep%nodes(node)%en, l_nb
         j = j + sz
         sz = sz - 1
      end do
      keep%nodes(node)%blk_en = keep%nodes(node-1)%blk_en + j

      ! if node is not the final node, record first block
      ! for the next node (node+1)
      if (node < num_nodes)  &
         keep%nodes(node+1)%blk_sa = keep%nodes(node)%blk_en + 1
   end do

   ! set keep%final_blk to hold total number of blocks.
   keep%final_blk = keep%nodes(num_nodes)%blk_en

   ! Add children to nodes, use sptr as a counter as it has fufilled its purpose
   sptr(:) = 0
   do node = 1, num_nodes
      par = sparent(node)
      if(par.gt.num_nodes) cycle
      sptr(par) = sptr(par) + 1
      keep%nodes(par)%child(sptr(par)) = node
   end do

   ! Setup least descendants, to allow easy walk of subtrees
   do node = -1, num_nodes
      ! initialise least descendat to self
      keep%nodes(node)%least_desc = node
   end do
   do node = 1, num_nodes
      ! walk up tree from leaves. A parent's least descendant is either this
      ! nodes least descendant (itself in case of a leaf), or its existing
      ! one if that is smaller.
      anode = keep%nodes(node)%parent
      keep%nodes(anode)%least_desc = &
         min(keep%nodes(node)%least_desc, keep%nodes(anode)%least_desc)
   end do

   !**************************************   
   ! Fill out block information. 

!**call system_clock(t_start)

   deallocate(keep%blocks,stat=st)
   allocate(keep%blocks(keep%final_blk),stat=st)
   if(st.ne.0) go to 490

   ! Loop over the nodes. Number the blocks in the first node
   ! contiguously, then those in the second node, and so on.
   ! Each node has a number of block columns; the blocks within
   ! each block column are numbered contiguously.
   ! Also add up the number of block columns and store largest block dimension.
   blk = 1
   keep%nbcol = 0
   keep%maxm = 0
   keep%maxn = 0
   do node = 1, num_nodes

      sa = keep%nodes(node)%sa
      en = keep%nodes(node)%en
      numcol = en - sa + 1
      numrow = size(keep%nodes(node)%index)

      ! l_nb is the size of the blocks
      l_nb = keep%nodes(node)%nb

      ! sz is number of blocks in the current block column
      sz = (numrow - 1) / l_nb + 1

      ! cb is the index of the block col. within node
      cb = 0
      col_used = 0

      ! Loop over the block columns in node. 
      do ci = sa, en, l_nb
         k = 1 ! use k to hold position of block within block column
         ! increment count of block columns
         keep%nbcol = keep%nbcol + 1

         cb = cb + 1

         ! blkn is the number of columns in the block column.
         ! For all but the last block col. blkn = l_nb.
         blkn = min(l_nb, numcol-col_used)
         col_used = col_used + blkn

         dblk = blk

         ! loop over the row blocks (that is, loop over blocks in block col)
         row_used = 0 
         do blk = dblk, dblk+sz-1
            ! store identity of block
            keep%blocks(blk)%id       = blk

            ! store number of rows in the block.
            ! For all but the last block, the number of rows is l_nb
            keep%blocks(blk)%blkm     = min(l_nb, numrow-row_used)
            row_used = row_used + keep%blocks(blk)%blkm

            ! store number of columns in the block.
            keep%blocks(blk)%blkn     = blkn

            keep%maxm = max(keep%maxm, keep%blocks(blk)%blkm)
            keep%maxn = max(keep%maxn, keep%blocks(blk)%blkn)

            ! store position of the first entry of the block within the
            ! block column of L
            keep%blocks(blk)%sa       = k

            ! store identity of diagonal block within current block column
            keep%blocks(blk)%dblk     = dblk

            ! store identity of last block within current block column
            keep%blocks(blk)%last_blk = dblk + sz - 1

            ! store node the blk belongs to
            keep%blocks(blk)%node     = node

            ! initialise dependency count
            keep%blocks(blk)%dep_initial = cb

            ! store identity of block column that blk belongs to
            keep%blocks(blk)%bcol     = keep%nbcol

!$          call omp_init_lock(keep%blocks(blk)%alock)

            ! increment k by number of entries in block
            k = k + keep%blocks(blk)%blkm * keep%blocks(blk)%blkn

         end do

         ! Diagonal block has no dependency for factor(dblk)
         keep%blocks(dblk)%dep_initial = cb - 1 

         ! decrement number of row blocks and rows in next block column
         sz = sz - 1
         numrow = numrow - l_nb
      end do
   end do

!**if(control%time_out.ge.0) then
!**   call system_clock(t_end, t_rate)
!**   write(control%time_out,"(a,es12.4)") "fill block took ", &
!**   (t_end - t_start) / real(t_rate)
!**end if

   !
   ! Compute dependency counts
   ! Note: This might be more efficient if implemented in left-looking
   ! (all descendants) rather than right-looking (all ancestors) fashion.
   !
!**call system_clock(t_start)
   allocate (map(n),stat=st)
   if(st.ne.0) go to 490

   ! loop over nodes
   ! FIXME: this loop can be particuarly expensive, and should perhaps be
   ! reordered so the map biulding loop only occurs once for each node.
   ! (eg PARSEC/Ga41As41H72 is particularly bad)
   do node = 1, num_nodes
      ! anode is initially the parent of node. Later it will be the
      ! grandparent, then great grandparent as we move up the tree to the root
      anode = keep%nodes(node)%parent
      ! numcol is number of columns in node
      numcol = keep%nodes(node)%en - keep%nodes(node)%sa + 1

      ! initialise cptr 
      cptr = 1 + numcol

      ! set swidth to number of block columns in node
      l_nb = keep%nodes(node)%nb
      swidth = (numcol-1)/l_nb + 1

      ! loop over ancestors of node
      do while(anode.gt.0)
         ! if we have finished with node, move to next node
         if(cptr.gt.size(keep%nodes(node)%index)) exit

         ! If we have skipped an anode (eg if its only a parent because of 
         ! other nodes in the subtree) we skip the current anode
         if(keep%nodes(node)%index(cptr).gt.keep%nodes(anode)%en) then
            anode = keep%nodes(anode)%parent
            cycle
         endif

         ! Build a map of anode's blocks. 
         ! Within the matrix for anode, the block columns are numbered
         ! 1,2,3... For each row k1 in jb-th block column,
         ! map(k1) is set to jb.

         a_nb = keep%nodes(anode)%nb
         jb = 1 ! Block
         ! loop over the block columns in anode
         size_anode = size(keep%nodes(anode)%index)
         do i = 1, size_anode, a_nb
            ! loop over the rows in the block column
            do k = i, min(i+a_nb-1, size_anode)
               k1 = keep%nodes(anode)%index(k)
               map(k1) = jb
            end do
            jb = jb + 1
         end do
         !print *, "   map = ", map

         ! Loop over affected block columns
         call calc_dep(cptr, node, anode, keep%nodes, keep%blocks, &
            swidth, map)

         ! Move up the tree to the parent of anode
         anode = keep%nodes(anode)%parent 
      end do
   end do

   allocate(amap(ptr(n+1)-1), keep%lmap(keep%nbcol), stat=st)
   if(st.ne.0) goto 490
   ! The output is keep%lmap, which holds mapping from user's matrix values
   ! into the block columns of the matrix factor
   call make_map(n, order, ptr, row, aptr, arow, amap)
   call lcol_map(aptr, arow, num_nodes, keep%nodes, keep%blocks, &
      keep%lmap, map, amap, st)
   if(st.ne.0) goto 490

!**if(control%time_out.ge.0) then
!**   call system_clock(t_end)
!**   write(control%time_out,"(a,es12.4)") &
!**      "calculating initial dependencies took ", &
!**       (t_end - t_start) / real(t_rate)
!**end if

   if (mp < 0) go to 500

   if (control%diagnostics_level >= 1) then
      write (mp,'(/a)') ' Leaving MA86_analyse with:'
      write (mp,'(a,i15)')    ' flag              = ',info%flag
      write (mp,'(a,i15)')    ' maxdepth          = ',info%maxdepth
      write (mp,'(a,es15.5)') ' num_factor        = ',real(info%num_factor)
      write (mp,'(a,i15)')    ' num_nodes         = ',info%num_nodes
      write (mp,'(a,i15)')    ' stat              = ',info%stat
   end if
   if (control%diagnostics_level>2) then
      ! Print out pivot order.
      write (mp,'(a)') ' On exit, elimination order :'
      write (mp,'(5i15)') order(1:n)
   else if (control%diagnostics_level==2) then
      i = min(10,n)
      write (mp,'(a,2(/5i12))')  &
         ' On exit, elimination order :', order(1:i)
      if (i < n) write (mp,'(a)') '  . . . . . .'
   end if

   go to 500

   490 info%flag = MA86_ERROR_ALLOCATION
       info%stat = st
       call MA86_print_flag(info%flag, control, context='MA86_analyse',st=st)

   500 continue
   ! before returning take copy of components of info set by MA86_analyse
   keep%info%flag         = info%flag
   keep%info%num_factor   = info%num_factor
   keep%info%num_flops    = info%num_flops
   keep%info%num_nodes    = info%num_nodes
   keep%info%maxdepth     = info%maxdepth
   keep%info%stat         = info%stat

   deallocate (arow,stat=st)
   deallocate (aptr,stat=st)
   deallocate (amap,stat=st)
   deallocate (iw,stat=st)
   deallocate (map,stat=st)
   deallocate (perm,stat=st)

end subroutine MA86_analyse_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Make a map from original A to reordered half matrix A
! The reordered half matrix's pattern is returned in nptr and nrow
subroutine make_map(n, perm, optr, orow, nptr, nrow, map)
   integer, intent(in) :: n
   integer, dimension(n), intent(in) :: perm
   integer, dimension(n+1), intent(in) :: optr
   integer, dimension(optr(n+1)-1), intent(in) :: orow
   integer, dimension(n+3), intent(out) :: nptr ! extra space used for tricks
   integer, dimension(optr(n+1)-1), intent(out) :: nrow
   integer, dimension(optr(n+1)-1), intent(out) :: map

   integer :: i, k, l
   integer(long) :: j

   nptr(:) = 0

   ! Count number of entries in each column of new matrix (at offset 2)
   do i = 1, n
      l = perm(i)
      do j = optr(i), optr(i+1)-1
         k = perm(orow(j))
         if(k<l) then
            nptr(k+2) = nptr(k+2) + 1
         else
            nptr(l+2) = nptr(l+2) + 1
         endif
      end do
   end do

   ! Calculate column starts (at offset 1)
   nptr(1:2) = 1
   do i = 2, n
      nptr(i+1) = nptr(i) + nptr(i+1)
   end do

   ! Now build map
   do i = 1, n
      l = perm(i)
      do j = optr(i), optr(i+1)-1
         k = perm(orow(j))
         if(k<l) then
            ! in upper triangle of reordered matrix - need to take conjugate
            map(nptr(k+1)) = -j ! use minus sign to indicate conjugate
            nrow(nptr(k+1)) = l
            nptr(k+1) = nptr(k+1) + 1
         else
            ! in lower triangle of reordered matrix
            map(nptr(l+1)) = j
            nrow(nptr(l+1)) = k
            nptr(l+1) = nptr(l+1) + 1
         endif
      end do
   end do
end subroutine make_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Build mapping on per block column basis from user's val to block col's lcol
! This routine uses the reordered half matrix and map from the make_map routine
subroutine lcol_map(aptr, arow, num_nodes, nodes, blocks, lmap, map, amap, st)
   ! Reordered lower triangle of reordered matrix held using aptr and arow
   integer, dimension(:), intent(in) :: aptr
   integer, dimension(:), intent(in) :: arow
   integer, intent(in) :: num_nodes
   type(node_type), dimension(-1:), intent(in) :: nodes ! Node info
   type(block_type), dimension(:), intent(in) :: blocks ! block info
   type(lmap_type), dimension(:), intent(out) :: lmap ! output lcol map
   integer, dimension(:), intent(out) :: map ! work array
   integer, dimension(:), intent(in) :: amap ! map set up by make_map
   integer, intent(out) :: st

   ! Local scalars
   integer :: bcol ! block column
   integer :: cb ! Temporary variable
   integer :: col ! do loop index
   integer(long) :: dblk ! set to keep%nodes(snode)%blk_sa (first block
     ! in snode which is, of course, a diagonal block)
   integer :: en ! set keep%nodes(snode)%en
   integer(long) :: i ! Temporary variable. global row index
   integer(long) :: j ! Temporary variable
   integer :: l_nb ! set to keep%nodes(snode)%nb
   integer :: nrow
   integer(long) :: offset
   integer :: sa ! set keep%nodes(snode)%sa
   integer :: snode ! node
   integer :: swidth ! set to keep%blocks(dblk)%blkn (number of columns
     ! in block column to which dblk belongs)

   integer :: k

   st = 0

   do snode = 1, num_nodes ! loop over nodes
      ! Build a map from global to local indices
      do j = 1, size(nodes(snode)%index)
        i = nodes(snode)%index(j)
        map(i) = j - 1
      end do

      ! Fill in lfact by block columns
      dblk = nodes(snode)%blk_sa

      l_nb = nodes(snode)%nb
      sa = nodes(snode)%sa
      en = nodes(snode)%en

      do cb = sa, en, l_nb
         bcol = blocks(dblk)%bcol

         offset = blocks(dblk)%sa - (cb-sa)*blocks(dblk)%blkn
         nrow = size(nodes(snode)%index) - (cb-sa)
         swidth = blocks(dblk)%blkn

         k = 1
         lmap(bcol)%len_map = aptr(min(cb+l_nb-1,en)+1) - aptr(cb)
         allocate(lmap(bcol)%map(2,lmap(bcol)%len_map), stat=st)
         if(st.ne.0) return

         ! Loop over columns in the block column bcol
         if(nodes(snode)%nchild.eq.0 .and. cb.eq.sa) then
            ! Leaf node, go colwise
            do col = cb, min(cb+l_nb-1, en)
               ! loop over rows in column col
               do j = aptr(col), aptr(col+1)-1
                  i = arow(j)
                  i = map(i)
                  i = offset + i
                  lmap(bcol)%map(1,k) = i ! destination in lfact(:)%lcol
                  lmap(bcol)%map(2,k) = amap(j) ! source
                  k = k + 1
               end do
               offset = offset + nrow
            end do
         else
            ! Non-leaf node, go rowwise
            do col = cb, min(cb+l_nb-1, en)
               ! loop over rows in column col
               do j = aptr(col), aptr(col+1)-1
                  i = arow(j)
                  i = map(i)
                  i = offset + i*swidth
                  lmap(bcol)%map(1,k) = i ! destination in lfact(:)%lcol
                  lmap(bcol)%map(2,k) = amap(j) ! source
                  k = k + 1
               end do
               offset = offset + 1
            end do
         endif
         ! move to next block column in snode
         dblk = blocks(dblk)%last_blk + 1
      end do
   end do
end subroutine lcol_map

!****************************************************************************

!
! Factorisation phase.
!
subroutine MA86_factor_complex(matrix_type, n, ptr, row, val, order, &
      keep, control, info, scale)
   integer, intent(in) :: matrix_type ! must be set to -4 for
      ! Hermitian and -5 for complex symmetric
   integer, intent(in) :: n ! order of A
   integer, intent(in) :: row(:) ! row indices of lower triangular part
   integer, intent(in) :: ptr(:) ! col pointers for lower triangular part
   complex(wp), intent(in) :: val(:) ! matrix values for lower triangular part
   integer, intent(in) :: order(:) ! holds pivot order (must be unchanged
      ! since the analyse phase)
   type(MA86_keep), intent(inout) :: keep ! see description of derived type
   type(MA86_control), intent(in) :: control ! see description of derived type
   type(MA86_info), intent(out) :: info ! see description of derived type
   real(wp), optional, intent(inout) :: scale(*) ! optional scaling factors

   integer :: i ! temporary variable
   integer :: mp ! copy of control%unit_diagnostics
   integer :: ne ! set to ptr(n+1) - 1
!**integer :: t_start, t_end, t_rate
   complex(wp) :: soln(0)
   integer :: st


   ! Reset components of info and return if an error was encountered earlier
   select case (keep%info%flag)
   case(MA86_ERROR_ALLOCATION, MA86_ERROR_ORDER)
      ! Previous call was an analyse or unrecoverable error
      return ! Cannot do factorize
   case default
      ! In all other cases we can do a factorize, reset flag to 0 and proceed
      info%flag = 0
   end select
   info%num_factor   = keep%info%num_factor
   info%num_flops    = keep%info%num_flops
   info%num_nodes    = keep%info%num_nodes
   info%maxdepth     = keep%info%maxdepth
   info%stat         = keep%info%stat

   ne = ptr(n+1) - 1
   ! Perform appropriate printing
   mp = control%unit_diagnostics
   if (control%diagnostics_level>=1 .and. mp>=0) then
      write (mp,'(/a)') ' On entry to MA86_factor:'
      write (mp,'(a,i15)') ' control%diagnostics_level =  ', &
         control%diagnostics_level
      write (mp,'(a,i15)') ' control%unit_diagnostics  =  ',mp
      write (mp,'(a,i15)') ' control%unit_error        =  ',control%unit_error
      write (mp,'(a,i15)') ' control%unit_warning      =  ',control%unit_warning
      write (mp,'(a,i15)') ' control%pool_size         =  ', &
        control%pool_size
      write (mp,'(a,i15)')    ' control%nbi               =  ', control%nbi
      write (mp,'(a,i15)') ' control%scaling           =  ', control%scaling
      write (mp,'(a,es15.5)') ' control%small             =  ', control%small
      write (mp,'(a,es15.5)') ' control%static            =  ', control%static
      write (mp,'(a,i15)') ' matrix_type               =  ',matrix_type
      write (mp,'(a,i15)') ' n                         =  ',n
   end if

   if (control%diagnostics_level > 2 .and. mp>=0) then

      write (mp,'(a)') ' ptr = '
      write (mp,'(5i15)') ptr(1:n+1)

      write (mp,'(a)') ' row = '
      write (mp,'(5i15)') row(1:ne)

      write (mp,'(a)') ' val = '
      write (mp,'(4es14.6)') val(1:ne)

      write (mp,'(a)') ' Elimination order :'
      write (mp,'(5i15)') order(1:n)

   else if (control%diagnostics_level==2 .and. mp>=0) then

      write (mp,'(a)') ' ptr(1:min(5,n+1)) = '
      write (mp,'(5i15)') ptr(1:min(5,n+1))

      write (mp,'(a)') ' row(1:min(5,ne)) =  '
      write (mp,'(5i15)') row(1:min(5,ne))

      write (mp,'(a)') ' val(1:min(5,ne)) =  '
      write (mp,'(4es14.6)') val(1:min(5,ne))

      i = min(10,n)
      write (mp,'(a,2(/5i12))')  &
         ' Elimination order :', order(1:i)
      if (i < n) write (mp,'(a)') '  . . . . . .'

   end if

   if(matrix_type .ne. HSL_MATRIX_CPLX_HERM_INDEF .and. &
         matrix_type .ne. HSL_MATRIX_CPLX_SYM) then
      info%flag = MA86_ERROR_MATRIX_TYPE
      call MA86_print_flag(info%flag, control, context='MA86_factor')
      return
   endif
   keep%matrix_type = matrix_type

   if(control%static.lt.abs(control%small) .and. control%static.ne.rzero) then
      info%flag = MA86_ERROR_STATIC_SMALL
      call MA86_print_flag(info%flag, control, context='MA86_factor')
      return
   endif

   ! immediate return if n = 0
   if (n == 0) return

   ! Handle scaling
   deallocate(keep%scaling, stat=st)
   if(control%scaling.gt.0 .or. present(scale)) then
      allocate(keep%scaling(n), stat=st)
      if(st.ne.0) then
         info%flag = MA86_ERROR_ALLOCATION
         info%stat = st
         call MA86_print_flag(info%flag, control, &
            context='MA86_factor', st=st)
         return
      endif
   endif
   select case (control%scaling)
   case(:0) ! user supplied (scale present) or none (scale not present)
      st = 0
      if(present(scale)) keep%scaling(1:n) = scale(1:n)
   case(1) ! mc64 scaling
      call mc64_scale(n, ptr, row, val, keep%scaling, control, info%flag, st)
      if(info%flag.lt.0) return
      if(present(scale)) scale(1:n) = keep%scaling(1:n)
   case(2:) ! mc77 scaling
      call mc77_scale(n, ptr, row, val, keep%scaling, st)
      if(present(scale)) scale(1:n) = keep%scaling(1:n)
   end select
   if(st.ne.0) then
      info%flag = MA86_ERROR_ALLOCATION
      info%stat = st
      call MA86_print_flag(info%flag, control, &
         context='MA86_factor', st=st)
      return
   endif

   ! Ready to perform the sparse factorization
!**call system_clock(t_start)

   if(allocated(keep%scaling)) then
      call factorize_indef(matrix_type, val, order, keep, control, info, &
         0, 0, soln, scale=keep%scaling)
   else
      call factorize_indef(matrix_type, val, order, keep, control, info, &
         0, 0, soln)
   endif

   if (info%flag < 0) then
      keep%info%flag  = info%flag
      return
   end if

!**if(control%time_out.ge.0) then
!**   call system_clock(t_end, t_rate)
!**   write(control%time_out,"(a,es12.4)") "factorization took ", &
!**   (t_end - t_start) / real(t_rate)
!**end if

   if (control%diagnostics_level >= 1 .and. mp >= 0) then
      write (mp,'(/a)')       ' Leaving MA86_factor with:'
      write (mp,'(a,i15)')    ' flag              = ',info%flag
      write (mp,'(a,i15)')    ' matrix_rank       = ',info%matrix_rank
      write (mp,'(a,i15)')    ' num_delay         = ',info%num_delay
      write (mp,'(a,i15)')    ' num_nodes         = ',info%num_nodes
      write (mp,'(a,i15)')    ' num_factor        = ',info%num_factor
      write (mp,'(a,i15)')    ' num_flops         = ',info%num_flops
      write (mp,'(a,i15)')    ' num_two           = ',info%num_two
      write (mp,'(a,i15)')    ' num_neg           = ',info%num_neg
      write (mp,'(a,i15)')    ' num_perturbed     = ',info%num_perturbed
      write (mp,'(a,i15)')    ' pool_size         = ',info%pool_size
      write (mp,'(a,i15)')    ' stat              = ',info%stat
      write (mp,'(a,es12.5)') ' usmall            = ',info%usmall
   end if

   ! Take a copy of any components of info that may have changed
   keep%info%flag          = info%flag
   keep%info%matrix_rank   = info%matrix_rank
   keep%info%num_delay     = info%num_delay
   keep%info%num_nodes     = info%num_nodes
   keep%info%num_factor    = info%num_factor
   keep%info%num_flops     = info%num_flops
   keep%info%num_two       = info%num_two
   keep%info%num_neg       = info%num_neg 
   keep%info%num_perturbed = info%num_perturbed
   keep%info%pool_size     = info%pool_size
   keep%info%stat          = info%stat
   keep%info%usmall        = info%usmall

end subroutine MA86_factor_complex

!****************************************************************************

!
! Combined factor+solve (simplified interface for single rhs).
!
subroutine MA86_factor_solve_one_complex(matrix_type, n, ptr, row, &
      val, order, keep, control, info, x, scale)
   integer, intent(in) :: matrix_type ! must be set to -4 or -5
   integer, intent(in) :: n ! order of A
   integer, intent(in) :: row(:) ! row indices of lower triangular part
   integer, intent(in) :: ptr(:) ! col pointers for lower triangular part
   complex(wp), intent(in) :: val(:) ! matrix values
   integer, intent(in) :: order(:) ! holds pivot order (must be unchanged
      ! since the analyse phase)
   type(MA86_keep), intent(inout) :: keep ! see description of derived type
   type(MA86_control), intent(in) :: control ! see description of derived type
   type(MA86_info), intent(out) :: info ! see description of derived type
   complex(wp), intent(inout) :: x(keep%n) ! On entry, x must
      ! be set so that if i has been used to index a variable,
      ! x(i) is the corresponding component of the right-hand side.
      ! On exit, if i has been used to index a variable,
      ! x(i) holds solution for variable i.
   real(wp), optional, intent(inout) :: scale(*) ! optional scaling factors

   call MA86_factor_solve_mult_complex(matrix_type, n, ptr, row, &
      val, order, keep, control, info, 1, keep%n, x, scale=scale)

end subroutine MA86_factor_solve_one_complex

!****************************************************************************

!
! Combined factor+solve.
!
subroutine MA86_factor_solve_mult_complex(matrix_type, n, ptr, row, &
      val, order, keep, control, info, nrhs, lx, x, scale)
   integer, intent(in) :: matrix_type ! must be set to -4 or -5
   integer, intent(in) :: n ! order of A
   integer, intent(in) :: row(:) ! row indices of lower triangular part
   integer, intent(in) :: ptr(:) ! col pointers for lower triangular part
   complex(wp), intent(in) :: val(:) ! matrix values
   integer, intent(in) :: order(:) ! holds pivot order (must be unchanged
      ! since the analyse phase)
   type(MA86_keep), intent(inout) :: keep ! see description of derived type
   type(MA86_control), intent(in) :: control ! see description of derived type
   type(MA86_info), intent(out) :: info ! see description of derived type
   integer, intent(in) :: nrhs ! number of right-hand sides to solver for.
   integer, intent(in) :: lx ! first dimension of x
   complex(wp), intent(inout) :: x(lx,nrhs) ! On entry, x must
      ! be set so that if i has been used to index a variable,
      ! x(i,j) is the corresponding component of the
      ! right-hand side for the jth system (j = 1,2,..., nrhs).
      ! On exit, if i has been used to index a variable,
      ! x(i,j) holds solution for variable i to system j
   real(wp), optional, intent(inout) :: scale(*) ! optional scaling factors

   integer :: i ! temporary variable
   integer :: j ! temporary variable
   integer :: mp ! copy of control%unit_diagnostics
   integer :: ne ! set to ptr(n+1) - 1
!**integer :: t_start, t_end, t_rate
   integer :: st ! stat parameter

   complex(wp), dimension(:), allocatable :: soln ! allocated to have 
     ! size n*nrhs.
     ! used to hold reordered rhs and then overwritten by reordered solution.

   ! Reset components of info and return if an error was encountered earlier
   select case (keep%info%flag)
   case(MA86_ERROR_ALLOCATION, MA86_ERROR_ORDER)
      ! Previous call was an analyse or unrecoverable error
      return ! Cannot do factorize
   case default
      ! In all other cases we can do a factorize, reset flag to 0 and proceed
      info%flag = 0
   end select
   info%num_factor   = keep%info%num_factor
   info%num_flops    = keep%info%num_flops
   info%num_nodes    = keep%info%num_nodes
   info%maxdepth     = keep%info%maxdepth
   info%stat         = keep%info%stat

   ne = ptr(n+1) - 1
   ! Perform appropriate printing
   mp = control%unit_diagnostics
   if (control%diagnostics_level>=1 .and. mp>=0) then
      write (mp,'(/a)') ' On entry to MA86_factor:'
      write (mp,'(a,i15)') ' control%diagnostics_level =  ', &
         control%diagnostics_level
      write (mp,'(a,i15)') ' control%unit_diagnostics  =  ',mp
      write (mp,'(a,i15)') ' control%unit_error        =  ',control%unit_error
      write (mp,'(a,i15)') ' control%unit_warning      =  ',control%unit_warning
      write (mp,'(a,i15)') ' control%pool_size         =  ', &
        control%pool_size
      write (mp,'(a,i15)')    ' control%nbi               =  ', control%nbi
      write (mp,'(a,i15)') ' control%scaling           =  ', control%scaling
      write (mp,'(a,es15.5)') ' control%small             =  ', control%small
      write (mp,'(a,es15.5)') ' control%static            =  ', control%static
      write (mp,'(a,i15)') ' matrix_type               =  ',matrix_type
      write (mp,'(a,i15)') ' n                         =  ',n
   end if

   if (control%diagnostics_level > 2 .and. mp>=0) then

      write (mp,'(a)') ' ptr = '
      write (mp,'(5i15)') ptr(1:n+1)

      write (mp,'(a)') ' row = '
      write (mp,'(5i15)') row(1:ne)

      write (mp,'(a)') ' val = '
      write (mp,'(4es14.6)') val(1:ne)

      write (mp,'(a)') ' Elimination order :'
      write (mp,'(5i15)') order(1:n)

   else if (control%diagnostics_level==2 .and. mp>=0) then

      write (mp,'(a)') ' ptr(1:min(5,n+1)) = '
      write (mp,'(5i15)') ptr(1:min(5,n+1))

      write (mp,'(a)') ' row(1:min(5,ne)) =  '
      write (mp,'(5i15)') row(1:min(5,ne))

      write (mp,'(a)') ' val(1:min(5,ne)) =  '
      write (mp,'(4es14.6)') val(1:min(5,ne))

      i = min(10,n)
      write (mp,'(a,2(/5i12))')  &
         ' Elimination order :', order(1:i)
      if (i < n) write (mp,'(a)') '  . . . . . .'

   end if

   if(matrix_type .ne. HSL_MATRIX_CPLX_HERM_INDEF .and. &
         matrix_type .ne. HSL_MATRIX_CPLX_SYM) then
      info%flag = MA86_ERROR_MATRIX_TYPE
      call MA86_print_flag(info%flag, control, context='MA86_factor')
      return
   endif
   keep%matrix_type = matrix_type

   if(control%static.lt.abs(control%small) .and. control%static.ne.rzero) then
      info%flag = MA86_ERROR_STATIC_SMALL
      call MA86_print_flag(info%flag, control, context='MA86_factor')
      return
   endif

   ! immediate return if n = 0
   if (n == 0) return

   if (lx < n .or. nrhs < 1) then
      info%flag = MA86_ERROR_X_SIZE
      call MA86_print_flag(info%flag, control, context='MA86_factor_solve')
      return
   end if

   deallocate (soln,stat=st)

   allocate (soln(n*nrhs),stat=st)
   if (st.ne.0) then
      info%flag = MA86_ERROR_ALLOCATION
      info%stat = st
      call MA86_print_flag(info%flag, control, context='MA86_factor_solve', &
        st=st)
      return
   end if

   ! Handle scaling
   deallocate(keep%scaling, stat=st)
   if(control%scaling.gt.0 .or. present(scale)) then
      allocate(keep%scaling(n), stat=st)
      if(st.ne.0) then
         info%flag = MA86_ERROR_ALLOCATION
         info%stat = st
         call MA86_print_flag(info%flag, control, &
            context='MA86_factor', st=st)
         return
      endif
   endif
   select case (control%scaling)
   case(:0) ! user supplied (scale present) or none (scale not present)
      st = 0
      if(present(scale)) keep%scaling(1:n) = scale(1:n)
   case(1) ! mc64 scaling
      call mc64_scale(n, ptr, row, val, keep%scaling, control, info%flag, st)
      if(info%flag.lt.0) return
      if(present(scale)) scale(1:n) = keep%scaling(1:n)
   case(2:) ! mc77 scaling
      call mc77_scale(n, ptr, row, val, keep%scaling, st)
      if(present(scale)) scale(1:n) = keep%scaling(1:n)
   end select
   if(st.ne.0) then
      info%flag = MA86_ERROR_ALLOCATION
      info%stat = st
      call MA86_print_flag(info%flag, control, &
         context='MA86_factor', st=st)
      return
   endif

   !
   ! Scale right hand side if required
   !
   if(allocated(keep%scaling)) then
      do i = 1, nrhs
         x(1:keep%n,i) = x(1:keep%n,i)*keep%scaling(1:keep%n)
      end do
   endif

   ! Reorder rhs

   do i = 1, nrhs
      do j = 1, n
         soln((i-1)*n + order(j)) = x(j, i)
      end do
   end do

   ! Ready to perform the sparse factorization
!**call system_clock(t_start)

   if(allocated(keep%scaling)) then
      call factorize_indef(matrix_type, val, order, keep, control, info, &
         nrhs, n, soln, scale=keep%scaling)
   else
      call factorize_indef(matrix_type, val, order, keep, control, info, &
         nrhs, n, soln)
   endif
   if (info%flag < 0) go to 10

!**if(control%time_out.ge.0) then
!**   call system_clock(t_end, t_rate)
!**   write(control%time_out,"(a,es12.4)") "factorization took ", &
!**   (t_end - t_start) / real(t_rate)
!**end if

   ! Perform back substitution
   call solve_indef(matrix_type, SOLVE_JOB_D_AND_BWD, nrhs, soln, n, keep, &
      control, info)

   if (info%flag.lt.0) go to 10

   ! Reorder soln
   
   do i = 1, nrhs
      do j = 1, n
         x(j, i) = soln((i-1)*n + order(j))
      end do
   end do

   !
   ! Scale right hand side if required
   !
   if(allocated(keep%scaling)) then
      do i = 1, nrhs
         x(1:keep%n,i) = x(1:keep%n,i)*keep%scaling(1:keep%n)
      end do
   endif

   if (control%diagnostics_level >= 1 .and. mp >= 0) then
      write (mp,'(/a)')       ' Leaving MA86_factor with:'
      write (mp,'(a,i15)')    ' flag              = ',info%flag
      write (mp,'(a,i15)')    ' matrix_rank       = ',info%matrix_rank
      write (mp,'(a,i15)')    ' num_delay         = ',info%num_delay
      write (mp,'(a,i15)')    ' num_nodes         = ',info%num_nodes
      write (mp,'(a,i15)')    ' num_factor        = ',info%num_factor
      write (mp,'(a,i15)')    ' num_flops         = ',info%num_flops
      write (mp,'(a,i15)')    ' num_two           = ',info%num_two
      write (mp,'(a,i15)')    ' num_neg           = ',info%num_neg
      write (mp,'(a,i15)')    ' num_perturbed     = ',info%num_perturbed
      write (mp,'(a,i15)')    ' pool_size         = ',info%pool_size
      write (mp,'(a,i15)')    ' stat              = ',info%stat
      write (mp,'(a,es15.5)') ' usmall            = ',info%usmall
   end if

   if (control%diagnostics_level > 2 .and. mp>=0) then
      write (mp,'(a)') ' Computed solution for first right-hand side :'
      write (mp,'(4es14.6)') x(1:n,1)
   else if (control%diagnostics_level==2 .and. mp>=0) then
      i = min(10,n)
      write (mp,'(a)') ' Computed solution for first right-hand side :'
      write (mp,'(4es14.6)') x(1:i,1)
      if (i < n) write (mp,'(a)') '  . . . . . .'
   end if

 10 deallocate(soln,stat=st)

   ! Take a copy of any components of info that may have changed
   keep%info%flag          = info%flag
   keep%info%matrix_rank   = info%matrix_rank
   keep%info%num_delay     = info%num_delay
   keep%info%num_nodes     = info%num_nodes
   keep%info%num_factor    = info%num_factor
   keep%info%num_flops     = info%num_flops
   keep%info%num_two       = info%num_two
   keep%info%num_neg       = info%num_neg 
   keep%info%num_perturbed = info%num_perturbed
   keep%info%pool_size     = info%pool_size
   keep%info%stat          = info%stat
   keep%info%usmall        = info%usmall

end subroutine MA86_factor_solve_mult_complex

!*************************************************

!
! Solve phase. simplified interface for a single rhs
!
subroutine MA86_solve_one_complex(x,order,keep,control,info,job,scale)
   type(MA86_keep), intent(inout) :: keep
   complex(wp), intent(inout) :: x(keep%n) ! On entry, x must
      ! be set so that if i has been used to index a variable,
      ! x(i) is the corresponding component of the right-hand side.
      ! On exit, if i has been used to index a variable,
      ! x(i) holds solution for variable i.
   integer, intent(in) :: order(:) ! pivot order. must be unchanged
   ! For details of keep, control, info : see derived type description
   type(MA86_control), intent(in) :: control
   type(MA86_info), intent(out) :: info
   integer, optional, intent(in) :: job  ! used to indicate whether
      ! partial solution required
      ! job = 0 or absent: complete solve performed
      ! job = 1 : forward eliminations only (PLx = b)
      ! job = 2 : backsubs only ((PL)^Tx = b)
   real(wp), optional, intent(in) :: scale(*) ! deprecated, ignored

   call MA86_solve_mult_complex(1, keep%n, x, order, keep, &
      control, info, job)

end subroutine MA86_solve_one_complex

!*************************************************

!
! Solve phase. Optionally performs only the forward, diagonal or backward sub.
!
subroutine MA86_solve_mult_complex(nrhs,lx,x,order,keep, &
      control,info,job,scale)
   integer, intent(in) :: nrhs ! number of right-hand sides to solver for.
   integer, intent(in) :: lx ! first dimension of x
   complex(wp), intent(inout) :: x(lx,nrhs) ! On entry, x must
      ! be set so that if i has been used to index a variable,
      ! x(i,j) is the corresponding component of the
      ! right-hand side for the jth system (j = 1,2,..., nrhs).
      ! On exit, if i has been used to index a variable,
      ! x(i,j) holds solution for variable i to system j
   integer, intent(in) :: order(:) ! pivot order. must be unchanged
   ! For details of keep, control, info : see derived type description
   type(MA86_keep), intent(inout) :: keep
   type(MA86_control), intent(in) :: control
   type(MA86_info), intent(out) :: info
   integer, optional, intent(in) :: job  ! used to indicate whether
      ! partial solution required
      ! job = 0 or absent: complete solve performed
      ! job = 1 : forward eliminations only (PLX = B). 
      ! job = 2 : diagonal solve (DX = B)
      ! job = 3 : backsubs only ((PL)^TX = B)
      ! job = 4 : backsubs and diag solve only (D(PL)^TX = B).
   real(wp), optional, intent(in) :: scale(*) ! deprecated, ignored

   integer :: i
   integer :: j
   integer :: local_job ! set to job or 0 if job not present
   integer :: mp ! set to control%unit_diagnostics
   integer :: n ! order of system
   integer :: st ! stat parameter
   complex(wp), dimension(:), allocatable :: soln ! allocated to have size
     ! n*nrhs. used to hold reordered rhs and then overwritten by reorder
     ! solution. 

   ! Reset components of info and return if an error was encountered earlier
   info%flag            = keep%info%flag
   if (info%flag < 0) return
   info%matrix_rank     = keep%info%matrix_rank
   info%maxdepth        = keep%info%maxdepth
   info%num_delay       = keep%info%num_delay
   info%num_factor      = keep%info%num_factor
   info%num_flops       = keep%info%num_flops
   info%num_nodes       = keep%info%num_nodes
   info%num_two         = keep%info%num_two
   info%num_neg         = keep%info%num_neg
   info%num_perturbed   = keep%info%num_perturbed
   info%pool_size       = keep%info%pool_size
   info%stat            = keep%info%stat
   info%usmall          = keep%info%usmall

   ! Perform appropriate printing
   mp = control%unit_diagnostics
   if (control%diagnostics_level>=1 .and. mp>=0) then
      write (mp,'(/a)') ' On entry to MA86_solve:'
      write (mp,'(a,i15)') ' control%diagnostics_level =  ', &
         control%diagnostics_level
      write (mp,'(a,i15)') ' control%unit_diagnostics  =  ',mp
      write (mp,'(a,i15)') ' control%unit_error        =  ',control%unit_error
      write (mp,'(a,i15)') ' control%unit_warning      =  ',control%unit_warning
      write (mp,'(a,i15)') ' control%pool_size         =  ', &
         control%pool_size
      write (mp,'(a,i15)') ' nrhs                      =  ',nrhs
      write (mp,'(a,i15)') ' lx                        =  ',lx
      if (present(job)) &
      write (mp,'(a,i15)') ' job                       =  ',job
   end if

   local_job = 0
   if (present(job)) then
      select case (job)
      case (SOLVE_JOB_ALL,SOLVE_JOB_FWD,SOLVE_JOB_BWD)
         ! everything OK
      case (SOLVE_JOB_D,SOLVE_JOB_D_AND_BWD)
         ! also OK
      case default
         info%flag = MA86_ERROR_JOB_OOR
         call MA86_print_flag(info%flag, control, context='MA86_solve')
         return
      end select
      local_job = job
   end if

   n = keep%n

   ! immediate return if n = 0
   if (n == 0) return

   if (lx < n .or. nrhs < 1) then
      info%flag = MA86_ERROR_X_SIZE
      call MA86_print_flag(info%flag, control, context='MA86_solve')
      return
   end if

   !
   ! Scale right hand side if required
   !
   if(allocated(keep%scaling) .and. &
         (local_job.eq.SOLVE_JOB_ALL.or.local_job.eq.SOLVE_JOB_FWD)) then
      do i = 1, nrhs
         x(1:keep%n,i) = x(1:keep%n,i)*keep%scaling(1:keep%n)
      end do
   endif

   !
   ! Reorder rhs
   !
   deallocate(soln,stat=st)
   allocate(soln(n*nrhs),stat=st)
   if (st.ne.0) then
      info%flag = MA86_ERROR_ALLOCATION
      info%stat = st
      call MA86_print_flag(info%flag, control, context='MA86_solve',st=st)
      return
   end if

   do i = 1, nrhs
      do j = 1, n
         soln((i-1)*n + order(j)) = x(j, i)
      end do
   end do

   call solve_indef(keep%matrix_type, local_job, nrhs, soln, n, keep, control, &
      info)

   if (info%flag.lt.0) go to 10

   !
   ! Reorder soln
   !
   do i = 1, nrhs
      do j = 1, n
         x(j, i) = soln((i-1)*n + order(j))
      end do
   end do

   !
   ! Scale right hand side if required
   !
   if(allocated(keep%scaling) .and. &
         (local_job.eq.SOLVE_JOB_ALL .or. local_job.eq.SOLVE_JOB_BWD &
          .or. local_job.eq.SOLVE_JOB_D_AND_BWD) ) then
      do i = 1, nrhs
         x(1:keep%n,i) = x(1:keep%n,i)*keep%scaling(1:keep%n)
      end do
   endif

   if (control%diagnostics_level >= 1 .and. mp >= 0) then
      write (mp,'(/a)') ' Leaving MA86_solve with:'
      write (mp,'(a,i15)') ' flag              = ',info%flag
      write (mp,'(a,i15)') ' stat              = ',info%stat
   end if

   if (control%diagnostics_level > 2 .and. mp>=0) then
      write (mp,'(a)') ' Computed solution for first right-hand side :'
      write (mp,'(4es14.6)') x(1:n,1)
   else if (control%diagnostics_level==2 .and. mp>=0) then
      i = min(10,n)
      write (mp,'(a)') ' Computed solution for first right-hand side :'
      write (mp,'(4es14.6)') x(1:i,1)
      if (i < n) write (mp,'(a)') '  . . . . . .'
   end if

 10 deallocate(soln,stat=st)

   ! Take a copy of any components of info that may have changed
   keep%info%flag       = info%flag
   keep%info%stat       = info%stat

end subroutine MA86_solve_mult_complex

!*************************************************

!
! This routine must be called after all other calls to routines
! in the package.
!
subroutine MA86_finalise_complex(keep, control)
   type(MA86_keep), intent(inout) :: keep    ! See derived-type declaration
   type(MA86_control), intent(in) :: control ! See derived-type declaration

   integer :: i
   integer :: st     ! stat parameter

   if (control%diagnostics_level >= 1 .and. control%unit_diagnostics >= 0) then
      write (control%unit_diagnostics,'(/a)') &
         ' Entering MA86_finalise'
   end if

   if(allocated(keep%lfact)) then
     do i = 1, keep%nbcol
        ! Note: we need to be careful if we've failed to allocate sufficient
        ! memory then we may not have initialized all the lfact locks.
!$      if(allocated(keep%lfact(i)%lcol)) &
!$        call omp_destroy_lock(keep%lfact(i)%lock)
        deallocate(keep%lfact(i)%lcol,stat=st)
        deallocate(keep%lfact(i)%d,stat=st)
     end do
     deallocate(keep%lfact,stat=st)
     keep%nbcol=0    
   endif

   if(allocated(keep%blocks)) then
!$    do i = 1, keep%final_blk
!$       call omp_destroy_lock(keep%blocks(i)%alock)
!$    end do
      keep%final_blk = 0
      deallocate(keep%blocks)
   endif

   deallocate(keep%nodes,stat=st)
   deallocate(keep%flag_array,stat=st)

end subroutine MA86_finalise_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Call mc64 to get a scaling, then symmetrize it
!
subroutine mc64_scale(n, ptr, row, val, scaling, control, flag, st)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   complex(wp), dimension(ptr(n+1)-1), intent(in) :: val
   real(wp), dimension(n), intent(out) :: scaling
   type(ma86_control), intent(in) :: control
   integer, intent(out) :: flag
   integer, intent(out) :: st

   integer :: i, j, k, ndiag
   integer, dimension(:), allocatable :: ptr2, row2, iw, cperm
   real(wp), dimension(:), allocatable :: val2, dw

   real(wp) :: colmax

   ! Expand out, removing any explicit zeroes
   allocate(ptr2(n+1), row2(2*ptr(n+1)), val2(2*ptr(n+1)), stat=st)
   if(st.ne.0) return
   allocate(iw(5*n), cperm(n), dw(3*n+2*ptr(n+1)), stat=st)
   if(st.ne.0) return

   k = 1
   do i = 1, n
      ptr2(i) = k
      do j = ptr(i), ptr(i+1)-1
         if(abs(val(j)).eq.rzero) cycle
         row2(k) = row(j)
         val2(k) = abs(val(j))
         k = k + 1
      end do
   end do
   ptr2(n+1) = k
   call mc34_expand(n, row2, ptr2, iw, a=val2)

   ! Call mc64
   do i = 1,n
      colmax = maxval(val2(ptr2(i):ptr2(i+1)-1))
      if (colmax.ne.0) colmax = log(colmax)
      dw(2*n+i) = colmax
      val2(ptr2(i):ptr2(i+1)-1) = colmax - log(val2(ptr2(i):ptr2(i+1)-1))
   end do
   call mc64w(n,ptr2(n+1)-1,ptr2,row2,val2,cperm,ndiag, &
      iw(1),iw(n+1),iw(2*n+1),iw(3*n+1),iw(4*n+1), &
      dw(1),dw(n+1))
   if (ndiag.eq.n) then
      do i = 1,n
         if (dw(2*n+i).ne.rzero) then
            dw(n+i) = dw(n+i) - dw(2*n+i)
         else
            dw(n+i) = rzero
         endif
      end do
   elseif(.not. control%action) then
      ! Matrix is singular, abort
      flag = MA86_ERROR_SINGULAR
      call MA86_print_flag(flag, control, context='MA86_factor')
      return
   endif

   do i = 1, n
      scaling(i) = exp( ( dw(i) + dw(n+i) ) / 2 )
      ! Just ignore badly scaled rows/columns. (Probably due to singular matrix)
      !if(scaling(i).ge.1e10) scaling(i) = 1.0
   end do
end subroutine mc64_scale

! Following Ruiz and Ucar
! "A symmetry preserving algorithm for matrix scaling"
! we do one iteration of the inf norm and then 3 of the one norm.
subroutine mc77_scale(n, ptr, row, val, scaling, st)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   complex(wp), dimension(ptr(n+1)-1), intent(in) :: val
   real(wp), dimension(n), intent(out) :: scaling
   integer, intent(out) :: st

   integer :: i, j, k

   integer, allocatable, dimension(:) :: iw
   real(wp), allocatable, dimension(:) :: dw
   integer :: icntl(10), info(10)
   real(wp) :: cntl(10), rinfo(10)

   real(wp), dimension(:), allocatable :: val2

   ! Take absolute value of matrix to avoid overheads
   allocate(val2(ptr(n+1)-1), stat=st)
   if(st.ne.0) return
   val2(1:ptr(n+1)-1) = abs(val(1:ptr(n+1)-1))

   ! Set controls
   call mc77i(icntl, cntl)
   !icntl(1) = -1 ! error messages
   !icntl(2) = -1 ! warning messages
   !icntl(3) = -1 ! diagnostic messages
   icntl(4) = -1 ! disable checking
   icntl(5) = -1 ! absolute value precomputed
   icntl(6) = -1 ! symmetric matrix

   allocate(iw(2*n), dw(2*n), stat=st)
   if(st.ne.0) return

   ! Single iteration of inf norm
   icntl(7) = 1 ! max number of iterations
   call mc77a(0, n, n, ptr(n+1)-1, ptr, row, val2, iw, size(iw), &
      dw, size(dw), icntl, cntl, info, rinfo)

   ! Apply scaling
   ! Note: we could modify mc77 to take an input vector of scaling to make
   ! this step unnesscary
   do i = 1, n
      do j = ptr(i), ptr(i+1)-1
         k = row(j)
         val2(j) = val2(j) / (dw(i) * dw(k))
      end do
   end do
   scaling(1:n) = 1/dw(1:n)

   ! Up to 3 iterations of one norm
   icntl(7) = 3 ! max number of iterations
   call mc77a(1, n, n, ptr(n+1)-1, ptr, row, val2, iw, size(iw), &
      dw, size(dw), icntl, cntl, info, rinfo)

   scaling(1:n) = scaling(1:n) / dw(1:n)
   
end subroutine mc77_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Perform a sparse factorization. Called by the driver MA86_factor.
!
subroutine factorize_indef(matrix_type, val, order, keep, control, &
      info, nrhs, ldr, rhs, scale)
   integer, intent(in) :: matrix_type  ! = -4 (Hermitian) or
                                       ! = -5 (symmetric)
   complex(wp), dimension(*), intent(in) :: val ! matrix values (lower triangle)
   integer, intent(in) :: order(:)  ! holds pivot order
   type(MA86_keep), intent(inout) :: keep ! see description of derived type
   type(MA86_control), intent(in) :: control ! see description of derived type
   type(MA86_info), intent(inout) :: info ! see description of derived type
   integer, intent(in) :: nrhs  ! number of right-hand sides (maybe = 0)
   integer, intent(in) :: ldr  ! leading extent of rhs
   complex(wp), intent(inout) :: rhs(ldr*nrhs)  ! On entry holds rhs data. 
      ! Overwritten by partial solution (forward substitution performed).
   real(wp), optional, intent(in) :: scale(*)
   
   ! local derived types
   type(dagtask) :: task ! see description of derived type
   type(taskstack) :: stack ! see description of derived type
   
   ! local arrays
   integer, dimension(:), allocatable ::  invp ! used to hold inverse ordering
   integer, dimension(:), allocatable ::  map ! allocated to have size n.
     ! used in copying entries of user's matrix a into factor storage 
     ! (keep%fact).
   integer, dimension(:), allocatable ::  pos ! allocated to have size n.
     ! used in copying from a child into its parent.
   complex(wp), dimension(:,:), allocatable ::  rhs_local ! Local right-hand 
     ! side arrays. allocated to have size (nrhs*ldr,0:total_threads)

   ! local scalars
   integer :: bcol ! block col.
   integer(long) :: blk ! block identity
   integer(long) :: dblk ! diagonal block within block column
   integer :: dest ! destination column for calculating delay destination
   integer :: en ! holds keep%nodes(snode)%en
   integer :: flag ! Error flag
   integer(long) :: i
   integer :: l_nb ! set to block size of snode (keep%nodes(snode)%nb)
   integer :: local ! used to set local index of block col. within snode
   integer :: nbcol ! number of block column
   integer :: nnodes ! number of nodes
   integer :: nrow ! used to count no. of rows in block col.
   integer :: parent ! parent of snode
   integer :: pool_size ! Normally set to control%pool_size
   integer :: sa ! holds keep%nodes(snode)%sa
   integer :: size_bcol ! number of entries in the block column (sum over the
     ! row blocks in the column)
   integer :: snode 
   integer :: st ! stat parameter
   integer :: sz ! number of blocks in a block column of snode
!**integer :: t_start, t_end, t_rate
   integer :: this_thread
   integer :: total_threads ! number of threads being used
   real(wp) :: u

   ! array for accumulating ma64 statistics
   type(thread_info), dimension(:), allocatable :: thread_stats

   ! Initialise
   flag = 0

   total_threads = 1
!$ total_threads = omp_get_max_threads()

   deallocate(thread_stats,stat=st)
   allocate(thread_stats(0:total_threads-1),stat=st)
   if(st.ne.0) go to 10

   u = control%u
   u = min(u, 0.5_wp)
   u = max(u, rzero)
   thread_stats(:)%usmall = u

   call zero_task(task)

   nnodes = keep%info%num_nodes

!**call system_clock(t_start)

   ! Initialize task pool
   pool_size = control%pool_size
   if (pool_size < 1) pool_size = pool_default

   call init_stack(stack, pool_size, control, flag, st)

   ! check for allocation error
   if (flag == MA86_ERROR_ALLOCATION) go to 10

   ! Set up inverse permutation
   deallocate (invp,stat=st)
   allocate (invp(keep%n),stat=st)
   if(st.ne.0) go to 10

   do i = 1, keep%n
      invp(order(i)) = i
   end do

   ! Allocate factor storage (in keep%lfact)
   ! Be careful to unitialise locks in any preexisting factorization!
   if(allocated(keep%lfact)) then
!$   do i = 1, keep%nbcol
!$      if(allocated(keep%lfact(i)%lcol)) &
!$        call omp_destroy_lock(keep%lfact(i)%lock)
!$   end do
     deallocate (keep%lfact,stat=st)
   endif
   allocate (keep%lfact(keep%nbcol),stat=st)
   if(st.ne.0) go to 10

   blk = 1
   nbcol = 0
   keep%lfact(:)%delay_head = -1
   ! loop over the nodes
   do snode = 1, nnodes
      ! Loop over the block columns in snode, allocating space 
      ! l_nb is the size of the blocks and sz is number of
      ! blocks in the current block column
      l_nb = keep%nodes(snode)%nb
      sz = (size(keep%nodes(snode)%index) - 1) / l_nb + 1
      sa = keep%nodes(snode)%sa
      en = keep%nodes(snode)%en

      local = 0
      do i = sa, en, l_nb
         nbcol = nbcol + 1
         local = local + 1
         ! store the local index of the block col within the node
         keep%lfact(nbcol)%local = local
         
         ! if this is last blk col of node, determine where delays end up
         parent = keep%nodes(snode)%parent
         if(i+l_nb.gt.en .and. parent.ne.-1) then
            ! go to some column of parent node
            dest = keep%nodes(snode)%index(en-sa+2) ! first row below diagonal
            ! now calculate which block column of parent dest corresponds to
            dest = (dest - keep%nodes(parent)%sa) / keep%nodes(parent)%nb + 1
            bcol = keep%blocks(keep%nodes(parent)%blk_sa)%bcol + dest - 1
            keep%nodes(snode)%delay_next = keep%lfact(bcol)%delay_head
            keep%lfact(bcol)%delay_head = snode
         endif

         ! store diagonal block
         dblk = blk
         keep%lfact(nbcol)%dblk = dblk

         ! loop over the row blocks adding up no. of rows
         nrow = 0
         do blk = dblk, dblk+sz-1
            nrow = nrow + keep%blocks(blk)%blkm
         end do
         keep%lfact(nbcol)%nrow = nrow

         ! allocate storage for L
         deallocate (keep%lfact(nbcol)%lcol, stat=st)
         size_bcol = nrow*keep%blocks(dblk)%blkn
         allocate (keep%lfact(nbcol)%lcol(size_bcol),stat=st)
         if(st.ne.0) go to 10
!$       call omp_init_lock(keep%lfact(nbcol)%lock)

         sz = sz - 1
      end do

      ! initialise number of delayed pivots to 0
      keep%nodes(snode)%num_delay = 0
   end do

   ! compute block column dependency counts by
   ! adding up all the dependencies of the blocks in each block col.

   keep%lfact(1:keep%nbcol)%dep = 0
   do i = 1,keep%final_blk
     bcol = keep%blocks(i)%bcol
     if (keep%blocks(i)%id.eq.keep%blocks(i)%dblk) then
       keep%lfact(bcol)%dep = keep%lfact(bcol)%dep + &
          keep%blocks(i)%dep_initial
     else
       ! for off diag. blocks, reduce count by one since no solve
       ! for combined factor-solve task
       keep%lfact(bcol)%dep = keep%lfact(bcol)%dep + &
         keep%blocks(i)%dep_initial - 1
     end if
     ! flag block as not yet touched
     keep%blocks(i)%touched = .false.
   end do

   ! Add initial tasks into the global task pool
   ! (i.e. those tasks with block col. dependency count equal to 0)
   task%task_type = TASK_FACTORIZE_COLUMN
   do i = 1, keep%nbcol
      if(keep%lfact(i)%dep.ne.0) cycle
      task%dest = keep%lfact(i)%dblk
      flag = 0
      call add_task_g(stack, task, control, flag, st)
      ! check for allocation error
      if(flag < 0) go to 10
   end do

!**if(control%time_out.ge.0) then
!**   call system_clock(t_end)
!**   write(control%time_out,"(a,es12.4)") "init took ", &
!**      (t_end - t_start) / real(t_rate)
!**end if

!**call system_clock(t_start, t_rate)

   ! Allocate parallel error array
   deallocate(keep%flag_array,stat=st)
   allocate(keep%flag_array(0:total_threads-1),stat=st)
   if(st.ne.0) go to 10

   ! Allocate local right-hand side arrays
   deallocate(rhs_local,stat=st)
   allocate(rhs_local(nrhs*ldr,0:total_threads-1),stat=st)

   10 if(st.ne.0) then
      info%flag = MA86_ERROR_ALLOCATION
      info%stat = st
      call cleanup_stack(stack)
      call MA86_print_flag(info%flag, control, context='MA86_factor',st=st)
      return
   endif

   ! initialise local error flags
   keep%flag_array = 0

   ! initialise rhs_local
   rhs_local(1:nrhs*ldr,0:total_threads-1) = czero

   !
   ! Copy matrix values across from a into keep%lfact
   !
!**call system_clock(t_start)
   st = 0
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(map, pos, flag, this_thread, st) &
!$OMP SHARED(control,info,invp,keep,ldr,matrix_type,nnodes,nrhs,order,rhs, &
!$OMP    rhs_local,scale,stack,thread_stats,total_threads,val &
!**!$OMP        ,t_start,t_end,t_rate &
!$OMP        )
      this_thread = 0
!$    this_thread = omp_get_thread_num()

      deallocate(map,stat=st)
      deallocate(pos,stat=st)
      allocate(map(keep%n),pos(keep%n),stat=st)
      if (st.ne.0) then
         keep%flag_array(this_thread) = MA86_ERROR_ALLOCATION
         info%stat = st
         call ma86_print_flag(flag, control, context='MA86_factor',st=st)
         go to 15
      end if
!$OMP BARRIER
      if(any(keep%flag_array(:).lt.0)) go to 15

      pos(1:keep%n) = 0

!$OMP BARRIER
!**!$OMP SINGLE
!**   if(control%time_out.ge.0) then
!**      call system_clock(t_end,t_rate)
!**      write(control%time_out,"(a,es12.4)") "copy matrix took ", &
!**      (t_end - t_start) / real(t_rate)
!**   end if

!**   call system_clock(t_start)
!**!$OMP END SINGLE NOWAIT

      !
      ! Perform actual factorization
      !
      keep%flag_array(this_thread) = info%flag
      call task_dispatch(matrix_type, val, invp, keep%nbcol, keep%maxm, &
         keep%maxn, keep%lfact, keep%lmap, map, pos, stack, keep%blocks, &
         keep%nodes, control, &
         thread_stats(this_thread), keep%flag_array(this_thread), st, nrhs, &
         rhs, ldr, total_threads, rhs_local, scale=scale)
      if(keep%flag_array(this_thread).lt.0) call set_abort(stack)

      if(keep%flag_array(this_thread).eq.MA86_ERROR_ALLOCATION) &
         info%stat = st

      ! Reductions
   15 continue
!$OMP BARRIER
!$OMP SINGLE
         deallocate(invp,stat=st)
         ! Carefully merge warnings or select most negative error
         do i = 0, total_threads-1
            if(keep%flag_array(i) .lt. 0) &
               info%flag = min(info%flag, keep%flag_array(i))
            if(info%flag.lt.0) cycle ! just trying to find minimum error
            flag = keep%flag_array(i)
            select case(flag)
            case(MA86_WARNING_SINGULAR)
               ! already found singular, was pool small also?
               if(info%flag.eq.MA86_WARNING_POOL_SMALL .or. &
                     info%flag.eq.MA86_WARNING_POOL_SING) &
                  flag = MA86_WARNING_POOL_SING
            case(MA86_WARNING_POOL_SMALL)
               ! already has a small pool, was it singular also?
               if(info%flag.eq.MA86_WARNING_SINGULAR .or. &
                     info%flag.eq.MA86_WARNING_POOL_SING) &
                  flag = MA86_WARNING_POOL_SING
            end select
            info%flag = max(info%flag, flag)
         end do
         if(info%flag.eq.MA86_WARNING_POOL_SING .or. &
               info%flag.eq.MA86_WARNING_SINGULAR ) &
            call MA86_print_flag(info%flag, control, context='MA86_factor')

         info%pool_size = stack%max_pool_size
         info%matrix_rank = keep%n - sum(thread_stats(:)%num_zero_pivots)
         info%num_factor = sum(thread_stats(:)%num_factor)
         info%num_flops = sum(thread_stats(:)%num_flops)
         info%num_delay = sum(thread_stats(:)%num_delay)
         info%num_two = sum(thread_stats(:)%num_two)
         info%num_neg = sum(thread_stats(:)%num_neg)
         info%num_perturbed = sum(thread_stats(:)%num_perturbed)
         info%num_nothresh = sum(thread_stats(:)%num_nothresh)
         info%usmall = minval(thread_stats(:)%usmall)
         info%detlog = sum(thread_stats(:)%detlog)
         info%detsign = product(thread_stats(:)%detsign)
         if(info%detsign.eq.0) info%detlog = 0.0_wp
!$OMP END SINGLE
      deallocate(map, stat=st)
      deallocate(pos, stat=st)
!$OMP END PARALLEL

!**if(control%time_out.ge.0) then
!**   call system_clock(t_end)
!**   write(control%time_out,"(a,es12.4)") "task_dispatch took ", &
!**      (t_end - t_start) / real(t_rate)
!**   write(control%time_out,"(a,50es12.4)")  &
!**      "Waiting time = ", stack%waiting(:)
!**end if

   ! write (6,*) 'max task pool size ', stack%max_pool_size

   call cleanup_stack(stack)
   deallocate(invp, stat=st)
   deallocate(rhs_local, stat=st)
   deallocate(thread_stats, stat=st)
end subroutine factorize_indef

!*************************************************

!
! Add entries from A into L.
! At a non-leaf node there are contributions from children to worry about,
! so these entries are added in to a rowwise data structure
!
subroutine blk_col_add_a(matrix_type, local, val, lmap, invp, snode, dblk, &
      lcol, scale)
   integer, intent(in) :: matrix_type
   integer, intent(in) :: local ! Block column position
   complex(wp), dimension(*), intent(in) :: val
   type(lmap_type), intent(in) :: lmap ! mapping set up by analyse phase
   integer, dimension(*), intent(in) :: invp
   type(node_type), intent(in) :: snode
   type(block_type), intent(in) :: dblk
   complex(wp), dimension(*), intent(inout) :: lcol
   real(wp), dimension(*), optional, intent(in) :: scale

   integer :: j, k, p, q, jcol, swidth, offset, offset2
   integer(long) :: i
   complex(wp) :: v

   swidth = dblk%blkn
   offset = snode%sa + (local-1)*snode%nb
   offset2 = (local-1)*snode%nb

   ! Loop over columns in the block column, adding in original matrix values
   select case(matrix_type)
   case(HSL_MATRIX_CPLX_HERM_INDEF)
      if(present(scale)) then
         do i = 1, lmap%len_map
            j = lmap%map(1,i)
            ! Calculate row and column of original matrix
            jcol = invp( offset + mod(j-1,swidth) ) ! col of original
            k = (j-1)/swidth + 1 ! row of local
            k = invp( snode%index(offset2+k) ) ! row of original
            p = lmap%map(2,i)
            q = abs(p)
            v = cmplx( real(val(q)), sign(1,-p)*aimag(val(q)), &
               kind=wp ) ! Take complex conjugate if required
            lcol(j) = lcol(j) + scale(jcol) * v * scale(k)
         end do
      else
         do i = 1, lmap%len_map
            p = lmap%map(2,i)
            q = abs(p)
            v = cmplx( real(val(q)), sign(1,-p)*aimag(val(q)), &
               kind=wp ) ! Take complex conjugate if required
            lcol(lmap%map(1,i)) = lcol(lmap%map(1,i)) + v
         end do
      endif
   case(HSL_MATRIX_CPLX_SYM)
      if(present(scale)) then
         do i = 1, lmap%len_map
            j = lmap%map(1,i)
            ! Calculate row and column of original matrix
            jcol = invp( offset + mod(j-1,swidth) ) ! col of original
            k = (j-1)/swidth + 1 ! row of local
            k = invp( snode%index(offset2+k) ) ! row of original
            q = abs(lmap%map(2,i))
            lcol(j) = lcol(j) + scale(jcol) * val(q) * scale(k)
         end do
      else
         do i = 1, lmap%len_map
            lcol(lmap%map(1,i)) = lcol(lmap%map(1,i)) + val(abs(lmap%map(2,i)))
         end do
      endif
   end select
end subroutine blk_col_add_a

!
! Add entries from A into L.
! At a leaf node there are no contributions from children to worry about,
! so this is built colwise
!
subroutine blk_col_add_a_leaf(matrix_type, val, lmap, invp, snode, dblk, &
      colwork, blkm, scale)
   integer, intent(in) :: matrix_type
   complex(wp), dimension(*), intent(in) :: val
   type(lmap_type), intent(in) :: lmap ! mapping set up by analyse phase
   integer, dimension(*), intent(in) :: invp
   type(node_type), intent(in) :: snode
   type(block_type), intent(in) :: dblk
   complex(wp), dimension(*), intent(out) :: colwork
   integer, intent(in) :: blkm
   real(wp), dimension(*), optional, intent(in) :: scale

   integer :: j, k, p, q, jcol, blkn, offset
   integer(long) :: i
   complex(wp) :: v

   blkn = dblk%blkn
   offset = snode%sa

   colwork(1:blkn*blkm) = czero

   ! Loop over columns in the block column, adding in original matrix values
   select case(matrix_type)
   case(HSL_MATRIX_CPLX_HERM_INDEF)
      if(present(scale)) then
         do i = 1, lmap%len_map
            j = lmap%map(1,i)
            ! Calculate row and column of original matrix
            jcol = invp( offset + (j-1)/blkm ) ! col of original
            k = mod(j-1,blkm) + 1 ! row of local
            k = invp( snode%index(k) ) ! row of original
            p = lmap%map(2,i)
            q = abs(p)
            v = cmplx( real(val(q)), sign(1,p)*aimag(val(q)), &
               kind=wp ) ! Take complex conjugate if required
            colwork(j) = colwork(j) + scale(jcol) * v * scale(k)
         end do
      else
         do i = 1, lmap%len_map
            j = lmap%map(1,i)
            p = lmap%map(2,i)
            q = abs(p)
            v = cmplx( real(val(q)), sign(1,p)*aimag(val(q)), &
               kind=wp ) ! Take complex conjugate if required
            colwork(j) = colwork(j) + v
         end do
      endif
   case(HSL_MATRIX_CPLX_SYM)
      if(present(scale)) then
         do i = 1, lmap%len_map
            j = lmap%map(1,i)
            ! Calculate row and column of original matrix
            jcol = invp( offset + (j-1)/blkm ) ! col of original
            k = mod(j-1,blkm) + 1 ! row of local
            k = invp( snode%index(k) ) ! row of original
            q = abs(lmap%map(2,i))
            colwork(j) = colwork(j) + scale(jcol) * val(q) * scale(k)
         end do
      else
         do i = 1, lmap%len_map
            j = lmap%map(1,i)
            colwork(j) = colwork(j) + val(abs(lmap%map(2,i)))
         end do
      endif
   end select

end subroutine blk_col_add_a_leaf

!*************************************************

!
! Performs solve. Called by driver MA86_solve.
!
subroutine solve_indef(matrix_type, job, nrhs, rhs, ldr, &
      keep, control, info)
   integer, intent(in) :: matrix_type        ! = -4 (Hermitian) or
                                             ! = -5 (symmetric)
   integer, intent(in) :: job                ! controls full or partial solve
   integer, intent(in) :: nrhs               ! number of rhs
   integer, intent(in) :: ldr                ! leading extent of rhs
   complex(wp), intent(inout) :: rhs(ldr*nrhs)  ! On entry holds rhs data. 
      ! Overwritten by solution.
   type(MA86_keep), intent(inout) :: keep
   type(MA86_control), intent(in) :: control
   type(MA86_info), intent(inout) :: info

   integer :: bcol ! block column (blocks(blk)%bcol)
   integer(long) :: blk ! block identifier
   type(slv_count_type), dimension(:), allocatable :: counts
   integer :: flag ! temporary store for info%flag
   integer :: i ! loop variable
   integer :: j ! temporary variable
   integer :: k ! index of rhs
   integer :: m   ! number of rows in block blk
   integer :: maxmn ! holds largest block dimension
   integer :: n   ! number of cols in block col.
   integer :: nelim  ! number of eliminations performed in block col.
   integer :: node ! temporary variable (for looping over nodes)
   integer :: num_nodes ! number of nodes
   integer :: offset_rhs
   integer :: pool_size ! inital size of task pool
   complex(wp), dimension(:), allocatable :: rhslocal ! allocated to have
      ! size (maxmn*nrhs).
   complex(wp), dimension(:,:), allocatable :: rhs_local
   integer :: s ! number of extra rows/cols in bcol because of delays
   integer :: st
   type(taskstack) :: stack
   type(dagtask) :: task
   integer :: this_thread
   integer :: total_threads

   integer :: col

   total_threads = 1
!$ total_threads = omp_get_max_threads()

   ! Find largest block dimension (m or n), allowing for delays
   num_nodes = keep%info%num_nodes
   maxmn = 0
   do node = 1, num_nodes
      do blk = keep%nodes(node)%blk_sa, keep%nodes(node)%blk_en
         bcol  = keep%blocks(blk)%bcol
         n     = keep%lfact(bcol)%blkn_new
         ! s is the number of extra rows/cols in bcol because of delays
         s     = n - keep%blocks(blk)%blkn
         m     = keep%blocks(blk)%blkm
         if (blk.eq.keep%blocks(blk)%dblk) m = m + s
         maxmn = max(maxmn, m,  n)
      end do
   end do

   ! Allocate workspace
   allocate(rhslocal(maxmn*nrhs),stat=st)
   if(st.ne.0) then
      info%flag = MA86_ERROR_ALLOCATION
      info%stat = st
      call MA86_print_flag(info%flag, control, context='MA86_solve',st=st)
      return
   endif

   ! Initialize task pool
   pool_size = control%pool_size
   if (pool_size < 1) pool_size = pool_default
   call init_stack(stack, pool_size, control, flag, st)
   if(st.ne.0) go to 10

   ! Allocate parallel error array
   deallocate(keep%flag_array,stat=st)
   allocate(keep%flag_array(0:total_threads-1),stat=st)
   if(st.ne.0) go to 10
   keep%flag_array(:) = 0

   col = 1
   do node = 1, num_nodes
      blk = keep%nodes(node)%blk_sa
      do while(blk.le.keep%nodes(node)%blk_en)
         bcol = keep%blocks(blk)%bcol

         keep%lfact(bcol)%col = col
         col = col + keep%lfact(bcol)%nelim

         blk = keep%blocks(blk)%last_blk+1
      end do
   end do

   if (job == SOLVE_JOB_ALL .or. job == SOLVE_JOB_FWD) then
      
      !
      ! Setup for forwards solve
      !

      allocate(counts(num_nodes),stat=st)
      if(st.ne.0) go to 10
      ! Next loop not integrated with add_task_g() loop to ensure locks
      ! can be destroyed safely if we get an error.
!$    do node = 1, num_nodes
!$       call omp_init_lock(counts(node)%lock)
!$    end do
      call zero_task(task)
      task%task_type = TASK_SLV_FSLV
      do node = 1, num_nodes
         counts(node)%dep = keep%nodes(node)%nchild
         if(counts(node)%dep.eq.0) then
            task%dest = node
            flag = 0
            call add_task_g(stack, task, control, flag, st)
            ! check for allocation error
            if(flag < 0) go to 10
         endif
      end do
   else
      allocate(counts(1),stat=st)
      if(st.ne.0) go to 10
   end if ! job == all or fwd

   if (job.eq.SOLVE_JOB_D) then
      ! Diagonal Solve only (otherwise combined with bwd solve)
      ! Loop over nodes.
      do node = 1, num_nodes
         ! set blk to be first block in node
         blk = keep%nodes(node)%blk_sa
         ! each pass of this loop deals with a block column in node
         do while(blk.le.keep%nodes(node)%blk_en)
            bcol  = keep%blocks(blk)%bcol
            nelim = keep%lfact(bcol)%nelim

            do k = 1, nrhs
               offset_rhs = (k-1)*ldr
               do i = 1, nelim
                  j = keep%lfact(bcol)%index_new(i)
                  rhslocal(i) = rhs(j+offset_rhs)
               end do

               call solveD(nelim, nelim, rhslocal, keep%lfact(bcol)%d, &
                  matrix_type)

               do i = 1, nelim
                  j = keep%lfact(bcol)%index_new(i)
                  rhs(j+offset_rhs) = rhslocal(i)
               end do
            end do

            blk = keep%blocks(blk)%last_blk + 1
         end do ! block column
      end do ! nodes
   endif ! diagonal solve

   if (job.eq.SOLVE_JOB_BWD .or. job.eq.SOLVE_JOB_D_AND_BWD) then
      
      !
      ! Setup for Backwards Solve (only if no fwds)
      !

      ! Add root nodes to task stack
      call zero_task(task)
      task%task_type = TASK_SLV_BSLV
      do node = 1, num_nodes
         if(keep%nodes(node)%parent.gt.0) cycle ! not a root

         ! If we reach here, then node is a root.
         task%dest = node
         flag = 0
         call add_task_g(stack, task, control, flag, st)
         ! check for allocation error
         if(flag < 0) go to 10
      end do
   end if ! bwd

   !
   ! Now execute the establish task graph in parallel
   !
   deallocate(rhs_local,stat=st)
   allocate(rhs_local(ldr*nrhs, 0:total_threads-1),stat=st)

   10 if(st.ne.0) then
      info%flag = MA86_ERROR_ALLOCATION
      info%stat = st
      call cleanup_stack(stack)
!$    if(allocated(counts) .and. &
!$          (job.eq.SOLVE_JOB_ALL .or. job.eq.SOLVE_JOB_FWD)) then
!$       do i = 1, size(counts)
!$          call omp_destroy_lock(counts(i)%lock)
!$       end do
!$    endif
      call MA86_print_flag(info%flag, control, context='MA86_solve',st=st)
      return
   endif
   rhs_local(:,:) = czero

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(this_thread) &
!$OMP SHARED(matrix_type, control, counts, info, job, keep, ldr, maxmn, nrhs, &
!$OMP    rhs, rhs_local, st, stack, total_threads)
   this_thread = 0
!$ this_thread = omp_get_thread_num()

   call solve_task_dispatch(matrix_type, keep%nbcol, keep%lfact, stack,       &
      keep%blocks, keep%nodes, counts, control, keep%flag_array(this_thread), &
      st, nrhs, ldr, rhs, total_threads, rhs_local, maxmn, job)

   if(keep%flag_array(this_thread).lt.0) call set_abort(stack)

   if(keep%flag_array(this_thread).eq.MA86_ERROR_ALLOCATION) &
      info%stat = st
!$OMP END PARALLEL

   ! Reduce flag_array nicely
   flag = info%flag
   info%flag = minval(keep%flag_array(:))
   if(info%flag.ge.0) &
      info%flag = max(flag, maxval(keep%flag_array(:)))

   !
   ! Dependency error checking - reneable by commenting if statement if
   ! needed for debugging
   !
   !if(.true.) then
   !   ! Check for errors in dep counting, fwd solve
   !   if(job.eq.SOLVE_JOB_ALL .or. job.eq.SOLVE_JOB_FWD) then
   !      do i = 1, size(counts)
   !         if(counts(i)%dep.ne.0) then
   !            print *, "fdep(", i, ") = ", counts(i)%dep
   !         endif
   !      end do
   !   endif
   !endif

   !
   ! Cleanup things that require explicit stuff
   !
   call cleanup_stack(stack)

   ! Destroy locks
!$ if(job.eq.SOLVE_JOB_ALL .or. job.eq.SOLVE_JOB_FWD) then
!$    do i = 1, size(counts)
!$       call omp_destroy_lock(counts(i)%lock)
!$    end do
!$ endif

end subroutine solve_indef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main task dispatch routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Main worker routine
!
subroutine task_dispatch(matrix_type, val, invp, nbcol, maxm, maxn, lfact, &
      lmap, map, pos, stack, blocks, nodes, control, tstats, info, st, &
      nrhs, rhs, ldr, total_threads, rhs_local, scale)
   integer, intent(in) :: matrix_type
   complex(wp), dimension(*), intent(in) :: val
   integer, dimension(*), intent(in) :: invp
   integer, intent(in) :: nbcol  ! Size of lfact (no. block cols)
   integer, intent(in) :: maxm ! Maximum value of blkm per analyse
   integer, intent(in) :: maxn ! Maximum value of blkn per analyse
   type(lfactor), dimension(nbcol), intent(inout) :: lfact ! Entries in 
      ! block columns of L
   type(lmap_type), dimension(nbcol), intent(in) :: lmap
   integer, dimension(:), intent(inout) :: map     ! private work array
      ! and used for mapping from child nodes into node
   integer, dimension(:), intent(inout) :: pos     ! 
      ! used for mapping from child nodes into node
   type(taskstack), intent(inout) :: stack         ! task pool
   type(block_type), dimension(:), intent(inout) :: blocks ! block info
   type(node_type), dimension(-1:), intent(inout) :: nodes ! Node info
      ! Components associated with indefinite case (delayed pivots) will
      ! be altered.
   type(MA86_control), intent(in) :: control
   type(thread_info), intent(inout) :: tstats
   integer, intent(out) :: info ! error flag
   integer, intent(out) :: st ! stat parameter
   integer, intent(in) :: nrhs  ! number of right-hand sides (maybe = 0)
   integer, intent(in) :: ldr  ! leading extent of rhs
   complex(wp), intent(inout) :: rhs(ldr*nrhs)  ! On entry holds rhs data. 
      ! Overwritten by partial solution (forward substitution performed).
   integer, intent(in) :: total_threads  ! number of threads
   complex(wp), intent(inout) :: rhs_local(ldr*nrhs,0:total_threads-1)  ! Local
      ! right-hand side arrays (one for each thread).
   real(wp), dimension(*), optional, intent(in) :: scale

   character(len=1) :: transa ! set to 'C' for Hermitian and 'T' for symmetric
   integer :: bcol ! block col that task%dest belongs to
   integer :: bcol_src ! block col that task%src1 belongs to
   integer(long) :: blk ! task%dest
   integer :: blkm ! set to number of rows in block
   integer :: blkn ! set to number of cols in block
   complex(wp), dimension(:), allocatable :: buffer ! work array
   complex(wp), dimension(:), allocatable :: colwork! buffer for storing colwise
   integer(long) :: dblk ! block on diagonal
   integer :: delay_col ! number of delayed cols added to the block col. that
     ! is passed to factor_solve
   integer :: flag
   type(dagtask) :: task
   integer :: i
   integer :: m ! set to blocks(blk)%blkm (number of rows in blk)
   integer :: m1 ! set to blocks(dblk)%blkm (number of rows in dblk)
   integer :: m2 ! set to blocks(task%src2)%blkm 
   integer :: maxm_local
   integer :: n ! set to blocks(blk)%blkn (number of columns in blk and dblk)
   integer :: n1 ! number of columns in source block col.
   integer :: nelim ! number of eliminations
   integer :: node 
   integer :: sa ! set to blocks(blk)%sa 
   integer :: sa1 ! set to blocks(dblk)%sa + s
   integer :: sa2 ! set to blocks(task%src2)%sa + s
   integer :: slocal ! local index of block column that source block belongs to 
!%%%   integer :: t_start, t_end
   integer :: this_thread
   complex(wp), dimension(:), allocatable :: work
   complex(wp), dimension(:), allocatable :: xlocal
   integer, dimension(:), allocatable :: col_list 
      ! work array -> update_between
   integer, dimension(:), allocatable :: row_list 
      ! work array -> update_between

   this_thread = 0
!$ this_thread = omp_get_thread_num()

   !if(control%diagnostics_level.gt.2 .and. control%unit_diagnostics.ge.0) &
   !   write(control%unit_diagnostics, "(a,i4,a)") &
   !      "Thread ", this_thread, " joined worker pool."

   ! Initialize things
   info = 0; st = 0             ! By default everything went OK
   call zero_task(task)         ! Set everything to zero
   task%task_type = TASK_NONE   ! Needs to be set to prevent 
                                ! dispatched incrementing
   maxm_local = 1               ! maximum block height on this thread

   transa = 'C'
   if (matrix_type == HSL_MATRIX_CPLX_SYM) transa = 'T'

   allocate(col_list(maxn), row_list(maxm), buffer(maxm*maxn), &
      colwork(1), work(1), xlocal(1), stat=st)
   if (st.ne.0) then
      info = MA86_ERROR_ALLOCATION
      call MA86_print_flag(info, control, context='MA86_factor',st=st)
      return
   end if

   ! Mark thread as active
!$ call omp_set_lock(stack%lock)
   stack%active = stack%active + 1
!$ call omp_unset_lock(stack%lock)

   !
   ! Main loop
   !
   do
      !
      ! Retrieve next task to perform and wait till all is ready
      !
!%%%  if(control%unit_log.gt.0) call system_clock(t_start)

      call get_task(stack, task, control, info, st)
      if(info.lt.0) return

!%%%  if(control%unit_log.gt.0) then
!%%%     call system_clock(t_end)
!%%%     call log_task(control, this_thread, t_start, t_end, "GT", &
!%%%        int(task%task_type,long))
!%%%  endif

!$OMP FLUSH(lfact)

      !if(control%diagnostics_level.gt.2 .and. control%unit_diagnostics.ge.0 &
      !      .and. task%task_type.ne.TASK_NONE .and. task%type.ne.TASK_DONE)  &
      !   then
      !   write(control%unit_diagnostics, "(i3,2a)") &
      !      this_thread, " got task ",print_job(task) 
      !endif

      !
      ! Perform task and update dependencies
      !
      select case(task%task_type)
      case(TASK_DONE)
         exit

      case(TASK_NONE) ! Job not finished but no tasks available, spin.
         cycle

      case(TASK_FACTORIZE_COLUMN) ! Factorize and solve
         dblk = task%dest
         bcol = blocks(dblk)%bcol
         node = blocks(dblk)%node

         !
         ! We handle leaf nodes seperatly to non-leaf nodes. The lack of
         ! children from which delays and updates can be inherited allows
         ! considerable simplifications compared to the case of a non-leaf
         ! node.
         ! In particular we construct the entries in a columnwise fashion
         ! for a leaf node, whereas a non-leaf node has the entries added
         ! in a rowwise fashion, and is copied to colwork in a columnwise
         ! fashion at the same time delays are incoporated.
         ! (Delays are held columnwise after the rowwise portion of lcol
         ! that holds the eliminated columns)
         !
         if(nodes(node)%nchild.eq.0 .and. lfact(bcol)%local.eq.1) then
            ! Leaf node
            blkm = lfact(bcol)%nrow
            blkn = blocks(dblk)%blkn

            if(size(colwork).lt.blkm*blkn) then
               deallocate(colwork,stat=st)
               allocate(colwork(blkm*blkn),stat=st)
               if (st.ne.0) goto 100
            endif

            ! Leaf node - no blocks touched yet, mark them as touched and zero
            ! them in blk_col_add_a_leaf.
            do blk = dblk, blocks(dblk)%last_blk
               m = blocks(blk)%blkm
               maxm_local = max(m, maxm_local)
               blocks(blk)%touched = .true.
            end do

            lfact(bcol)%blkn_new = blocks(dblk)%blkn

            deallocate (lfact(bcol)%d,stat=st)
            deallocate (lfact(bcol)%index_new,stat=st)
            allocate (lfact(bcol)%d(2*blkm), &
                      lfact(bcol)%index_new(blkm),stat=st)
            if (st.ne.0) goto 100
            lfact(bcol)%index_new(:) = nodes(node)%index(:)
             
            i = blocks(dblk)%last_blk
            blocks(dblk:i)%sa_new = blocks(dblk:i)%sa

            ! delay_col is number of delayed columns that have been added to
            ! the block col.
            delay_col = 0

            ! Expand original matrix entries into colwork
            call blk_col_add_a_leaf(matrix_type, val, lmap(bcol), invp, &
               nodes(node), blocks(dblk), colwork, blkm, scale=scale)
         else ! node has children
            ! If any blocks haven't been touched, zero them (eg first column)
            do blk = dblk, blocks(dblk)%last_blk
               m = blocks(blk)%blkm
               maxm_local = max(m, maxm_local)
               if(blocks(blk)%touched) cycle
               n = blocks(blk)%blkn
               lfact(bcol)%lcol(blocks(blk)%sa:blocks(blk)%sa+m*n-1) = czero
               blocks(blk)%touched = .true.
            end do

            ! Add original matrix entries to updates accummulated in lcol
            call blk_col_add_a(matrix_type, lfact(bcol)%local, val, &
               lmap(bcol), invp, nodes(node), blocks(dblk), lfact(bcol)%lcol, &
               scale=scale)

            ! Merge any delays into column
            ! On exit, blkm and blkn hold the number of rows and cols
            ! in the enlarged block col.
            call merge_delays(matrix_type, nodes, blocks, lfact, pos, dblk, &
               blkm, blkn, tstats, colwork, st)
            if (st.ne.0) then
               info = MA86_ERROR_ALLOCATION
               call MA86_print_flag(info, control, context='MA86_factor',st=st)
               return
            end if

            ! delay_col is number of delayed columns that have been added to
            ! the block col.
            delay_col = blkn - blocks(dblk)%blkn

            ! Handle diagonal block maxm carefully to include delays
            m = blocks(dblk)%blkm + delay_col
            maxm_local = max(maxm_local, m)
         endif

         if(size(work).lt.blkm*blkn) then
            deallocate(work,stat=st)
            allocate(work(blkm*blkn),stat=st)
            if(st.ne.0) goto 100
         endif
         ! Perform actual factorization
         call factor_solve_block(matrix_type, blkm, blkn, blkn, delay_col, &
            colwork, blkm*blkn, work, lfact(bcol)%d, map, lfact(bcol)%nelim, &
            tstats, control, flag)
         if(flag.lt.0) then
            if (flag.eq.-13) then
               ! IEEE infinity detected
               info = MA86_ERROR_INFINITY
            else
               info = MA86_ERROR_UNKNOWN
            end if
            call MA86_print_flag(info, control, context='MA86_factor')
            return
         endif

         if (matrix_type == HSL_MATRIX_CPLX_HERM_INDEF) then
            call col_to_row_herm(blkm,lfact(bcol)%nelim,lfact(bcol),colwork,st)
         else
            call col_to_row(blkm,lfact(bcol)%nelim,lfact(bcol),colwork,st)
         end if
         if(st.ne.0) goto 100

         ! Work out number of flops and entries in factor
         do i = 0, lfact(bcol)%nelim-1
            tstats%num_flops = tstats%num_flops + (blkm-i)**2
            tstats%num_factor = tstats%num_factor + (blkm-i)
         end do

! Use pos to take temporary copy (because f95 compiler with -O is giving me
! wrong answer if I use following) FIXME check this does cause problems
! (it did in ma77)
!   lfact(bcol)%index_new(1:blkn) = lfact(bcol)%index_new(map(1:blkn))
         pos(1:blkn) = lfact(bcol)%index_new(1:blkn)
         do i = 1,blkn
            lfact(bcol)%index_new(i) = pos(map(i))
         end do
         pos(1:blkn) = 0

         ! store the number of delayed pivots (note: nodes(node)%num_delay will
         ! get overwritten by value for the last block col. in node
         ! and this is what we have to pass to parent node)
         ! This is accumulated into tstats%num_delay in merge_delays
         ! if (blkn - lfact(bcol)%nelim.ne.0) &
         ! write(6,*) 'num_delay',node,blkn-lfact(bcol)%nelim
         nodes(node)%num_delay = blkn - lfact(bcol)%nelim

         ! Do the right thing if we're singular
         if(tstats%num_zero_pivots.ne.0) then
            if(control%action) then
               ! Continue but raise a warning
               if(info.ne.MA86_WARNING_SINGULAR .and. &
                     info.ne.MA86_WARNING_POOL_SING) then
                  if(info.eq.MA86_WARNING_POOL_SMALL) then
                     info = MA86_WARNING_POOL_SING
                  else
                     info = MA86_WARNING_SINGULAR
                  endif
                  ! Note: we don't print out a warning here as we only want to
                  ! print it once (not once per thread). We do so instead
                  ! at the end of the parallel section
               endif
            else
               ! Abort with an error
               info = MA86_ERROR_SINGULAR
               call MA86_print_flag(info, control, context='MA86_factor')
               return
            endif
         endif

         ! if we are doing a combined factor_solve, use the factorization
         ! to do a solve
         if(nrhs.gt.0) then
            ! Ensure work is large enough
            if(size(work).lt.lfact(bcol)%nelim*nrhs) then
               deallocate(work, stat=st)
               allocate(work(lfact(bcol)%nelim*nrhs), stat=st)
               if(st.ne.0) goto 100
            endif

            ! Ensure xlocal is large enough
            if(size(xlocal).lt.nrhs*maxm_local) then
               deallocate(xlocal, stat=st)
               allocate(xlocal(nrhs*maxm_local), stat=st)
               if(st.ne.0) goto 100
            endif

            blk = dblk
            call fwd_solve_bcol(matrix_type, bcol, blk, 1, lfact, blocks, &
               nrhs, rhs, ldr, work, lfact(bcol)%nelim, xlocal, rhs_local, &
               control,  this_thread, total_threads)
         endif

         !
         ! add update tasks
         !

!$OMP FLUSH 
  
         call add_updates_new(stack, map, nodes, blocks, &
            dblk, blocks(dblk)%node, control, info, st)
         if(info.lt.0) return

      case(TASK_UPDATE_INTERNAL) ! update one block by two others
         ! in the same node. If pivots have been delayed, we can 
         ! only update by the number of pivots selected.

         blk = task%dest
         n   = blocks(blk)%blkn
         m   = blocks(blk)%blkm
         sa  = blocks(blk)%sa
         bcol = blocks(blk)%bcol

         ! src1 and src2 belong to same block column. 
         ! They have m1 and m2 rows.

         m1  = blocks(task%src1)%blkm
         m2  = blocks(task%src2)%blkm

         bcol_src = blocks(task%src1)%bcol

         ! n1 is the number of columns in source block column bcol_src.
         n1  = lfact(bcol_src)%blkn_new

         ! nelim is the number of update operations (ie number of
         ! eliminations performed in the source block). Note that it can be 0.
         nelim  = lfact(bcol_src)%nelim

         ! determine position of the source blocks within bcol_src
         sa1 = swizzle_sa(blocks(task%src1)%sa_new, n1, nelim)
         sa2 = swizzle_sa(blocks(task%src2)%sa_new, n1, nelim)

         if (nelim.gt.0) then
            if(size(work).lt.n1*m1) then ! Resize work array for calc_ld
               deallocate(work)
               allocate(work(n1*m1),stat=st)
               if(st.ne.0) goto 100
            endif
            call update_block_block(transa, m, n, nelim, &
               lfact(bcol)%lcol(sa:sa+n*m-1), blocks(blk), nelim, &
               lfact(bcol_src)%lcol(sa1:sa1+nelim*m1-1), &
               lfact(bcol_src)%lcol(sa2:sa2+nelim*m2-1), &
               lfact(bcol_src)%d, work, &
               control)
         endif

         ! reduce dependency count for bcol
         call reduce_block_dep(stack, lfact(bcol), control, info, st)
         if(info.lt.0) return

      case(TASK_UPDATE_BETWEEN) ! update a block with information from 
                                ! another node

         blk = task%dest
         bcol = blocks(task%dest)%bcol

         ! task%src1 is a block in the source node
         ! bcol_src is the block col. it belongs to
         ! slocal is the the local index of the block col. within source node

         bcol_src = blocks(task%src1)%bcol
         slocal = lfact(bcol_src)%local
         n1  = lfact(bcol_src)%blkn_new

         nelim  = lfact(bcol_src)%nelim
         if (nelim.gt.0) then
            call update_between(transa, blk, nodes(blocks(blk)%node),   &
               n1, nelim, lfact(bcol_src)%index_new, lfact(bcol)%lcol,  &
               lfact(bcol_src)%lcol, lfact(bcol_src)%d,                 &
               blocks, col_list, row_list, buffer, control, info, st,   &
               slocal, nodes(blocks(task%src1)%node), work)
         endif
         if(info.lt.0) return

         ! reduce dependency count for bcol
         call reduce_block_dep(stack, lfact(bcol), control, info, st)
         if(info.lt.0) return

      !case default
      !   info = MA86_ERROR_UNKNOWN
      !   call MA86_print_flag(info, control, context='MA86_factor')
      !   exit
      end select

   end do

   return
   100 continue
   info = MA86_ERROR_ALLOCATION
   call ma86_print_flag(info, control, context='MA86_factor',st=st)
   return

end subroutine task_dispatch

!************************************************

!
! Before we can do the factor_solve, we have to add in columns
! from delayed pivots. Treat the case of first block column in node
! separately (since must pick up delays from children)
!
subroutine merge_delays(matrix_type, nodes, blocks, lfact, pos, dblk, blkm_new,&
      blkn_new, tstats, work, st)
   integer, intent(in) :: matrix_type
   type(node_type), dimension(-1:), intent(in) :: nodes ! Node info
   type(block_type), dimension(:), intent(inout) :: blocks ! block info
      ! For each block in the block column, the component sa_new is set
      ! and other components are unchanged.
   type(lfactor), dimension(:), intent(inout) :: lfact
   integer, dimension(:), intent(inout) :: pos     ! 
      ! used for mapping from child nodes into node
   integer(long), intent(in) :: dblk ! block on diagonal
   integer, intent(out) :: blkm_new ! new number of rows in the block col.
   integer, intent(out) :: blkn_new ! new number of columns in the block col.
   type(thread_info), intent(inout) :: tstats ! thread stats
   complex(wp), dimension(:), allocatable, target, intent(inout) :: work
   integer, intent(out) :: st ! stat parameter

   integer :: bcol
   integer :: bcol_cnode ! last block column in child node
   integer :: blk
   integer :: blkn ! set to blocks(dblk)%blkn
   integer :: cblkn ! 
   integer :: cblkm ! 
   integer :: cnum_delay ! number of delays at child node
   integer :: cnelim ! no. of eliminations at child node
   integer :: cnode ! a child of node
   integer :: cnrow
   integer(long) :: en_child ! identifier of last block in child node
   integer :: i
   integer :: ip
   integer :: ip1
   integer :: ipd
   integer :: j
   integer :: jj
   integer :: k
   integer :: lindex ! size of array index_new (holds indices in block col).
   integer :: l_nb ! block size
   integer :: local ! local index of block col.
   integer :: nd ! number of delayed pivots (set once only and not changed)
   integer :: num_delay ! number of delayed pivots (used as temporary)
   integer :: nelim ! number of eliminations
   integer :: nrow ! number of rows in block col.
   integer :: node
   integer :: nvar
   integer :: nzl ! size of nodes(node)%index(:)
   integer :: sa
   integer(long) :: size_bcol ! size of block col. when reallocated to
      ! allow for delays
   integer :: size_d ! number of entries in D
   logical :: do_copy ! flag indicating if delays are present

   ! bcol and node are block column and node that task%dest belongs to
   bcol = blocks(dblk)%bcol
   node = blocks(dblk)%node
   blkn = blocks(dblk)%blkn
   nrow = lfact(bcol)%nrow

   num_delay = 0
   local = lfact(bcol)%local

   ! Count any delays from previous blk col of this node
   if(local .ne. 1) then
      ! blkn_new is no. of columns in previous block col.
      ! Set num_delay to number of delayed pivots to pass from bcol-1 to bcol
      blkn_new = lfact(bcol-1)%blkn_new
      nelim = lfact(bcol-1)%nelim
      num_delay = blkn_new - nelim
   endif

   ! Count any delays from child nodes (don't always hit first col of a node!)
   cnode = lfact(bcol)%delay_head
   do while(cnode.ne.-1)
      num_delay = num_delay + nodes(cnode)%num_delay
      tstats%num_delay = tstats%num_delay + nodes(cnode)%num_delay
      cnode = nodes(cnode)%delay_next
   end do
   nd = num_delay

   do_copy = (nd.ne.0) ! If there are no delays things are easier
   lfact(bcol)%blkn_new = blkn + nd
   blkn_new = lfact(bcol)%blkn_new
   blkm_new = nrow + nd
   size_bcol = blkm_new * blkn_new

   !
   ! allocate storage for L and D
   ! allocate new array for permuted variable list for block col.
   ! (including delays from children)
   !
   ! Array index_new holds index list for block col.

   size_d     = 2*(blocks(dblk)%blkm + nd) ! diagonal pivots

   lindex = blkm_new

   deallocate (lfact(bcol)%d,stat=st)
   deallocate (lfact(bcol)%index_new,stat=st)
   allocate (lfact(bcol)%d(size_d),             &
             lfact(bcol)%index_new(lindex),stat=st)
   if (st.ne.0) return
             
   if(size(work).lt.size_bcol) then
      deallocate(work,stat=st)
      allocate(work(size_bcol),stat=st)
      if (st.ne.0) return
   endif

   ! Copy any delays from previous column of this node into index list
   num_delay = 0
   if(local .ne. 1) then
      cnelim = lfact(bcol-1)%nelim
      cnum_delay = lfact(bcol-1)%blkn_new - cnelim
      lfact(bcol)%index_new(num_delay+1:num_delay+cnum_delay) =  &
         lfact(bcol-1)%index_new(cnelim+1:cnelim+cnum_delay)
      num_delay = num_delay + cnum_delay
   endif

   ! Copy any delays from child nodes into index list
   cnode = lfact(bcol)%delay_head
   do while(cnode .ne. -1)
      cnum_delay = nodes(cnode)%num_delay

      ! determine which block column was the final block column of
      ! the child cnode. en_child is last block in child node.
      en_child = nodes(cnode)%blk_en
      bcol_cnode = blocks(en_child)%bcol
      cnelim = lfact(bcol_cnode)%nelim

      lfact(bcol)%index_new(num_delay+1:num_delay+cnum_delay) = &
         lfact(bcol_cnode)%index_new(cnelim+1:cnelim+cnum_delay)

      num_delay = num_delay + cnum_delay

      ! Move to next child with potential delays
      cnode = nodes(cnode)%delay_next
   end do

   nzl = size(nodes(node)%index)
   lfact(bcol)%index_new(nd+1:lindex) = &
      nodes(node)%index(nzl-lindex+nd+1:nzl)

   ! Compute new block starts so that block can be located within expanded
   ! block column.
   sa = 1
   blocks(dblk)%sa_new = sa
   sa = sa + (blocks(dblk)%blkm+nd) * blkn_new
   if (blocks(dblk)%last_blk.gt.dblk) then
     blocks(dblk+1)%sa_new = sa
     sa = sa + blocks(dblk+1)%blkm * blkn_new
     do blk = dblk+2, blocks(dblk)%last_blk
       blocks(blk)%sa_new = sa
       sa = sa + blocks(blk)%blkm * blkn_new
     end do
   end if

   !
   ! The remainder only applies if a copy is required (i.e. there are delays).
   !
   if(.not. do_copy) then
      if (matrix_type == HSL_MATRIX_CPLX_HERM_INDEF) then
         call row_to_col_herm(blkm_new,blkn_new,lfact(bcol)%lcol,work)
      else
         call row_to_col(blkm_new,blkn_new,lfact(bcol)%lcol,work)
      endif
      return
   endif

   ! copy contents of lcol into work (held by cols)
   ip = nd*blkm_new
   ip1 = 0
   ! Note: If we ONLY recieve entries from previous column of same node
   ! then all entries are overwritten and don't need to be set to zero.
   if(lfact(bcol)%delay_head.ne.-1) then
      ! Entries coming from child nodes in addition to any from this node's
      ! previous columns. Not all locations are overwritten, so set all
      ! entries in delayed columns to zero.
      work(1:ip) = czero
   endif

   ! copy columns from original version of block column
   if (matrix_type == HSL_MATRIX_CPLX_HERM_INDEF) then
      do i = 1,blkn
         ! Following line zeroes part of column above diagonal, it is
         ! not strictly required, but avoids various innocous undefined
         ! variable errors.
         work(ip+1:ip+nd) = czero
         ip  = ip + nd
         work(ip+1:ip+nrow) = &
            conjg(lfact(bcol)%lcol(ip1+1:blkn*nrow:blkn))
         ip  = ip + nrow
         ip1 = ip1 + 1
      end do
   else
      do i = 1,blkn
         ! Following line zeroes part of column above diagonal, it is
         ! not strictly required, but avoids various innocous undefined
         ! variable errors.
         work(ip+1:ip+nd) = czero
         ip  = ip + nd
            work(ip+1:ip+nrow) = &
               lfact(bcol)%lcol(ip1+1:blkn*nrow:blkn)
         ip  = ip + nrow
         ip1 = ip1 + 1
      end do
   endif

   num_delay = 0
   if(local.ne.1) then
      ! copy delayed columns from previous block column (bcol-1)
      cnelim = lfact(bcol-1)%nelim
      cnum_delay = lfact(bcol-1)%blkn_new - cnelim
      cblkm = lfact(bcol-1)%nrow + &
         ( lfact(bcol-1)%blkn_new - blocks(lfact(bcol-1)%dblk)%blkn )
      ip = 0
      ip1 = (cblkm+1) * cnelim
      do i = 1, cnum_delay
         ! Add rows matching expected elimination in bcol
         work(ip+1:ip+cnum_delay) = &
            lfact(bcol-1)%lcol(ip1+1:ip1+cnum_delay)
         ! Note: If we inherit from children there may be some zero rows
         ! to skip before adding remaining rows
         work(ip+nd+1:ip+nd+nrow) = &
            lfact(bcol-1)%lcol(ip1+cnum_delay+1:ip1+cnum_delay+nrow)
         ip = ip + blkm_new
         ip1 = ip1 + cblkm
      end do
      num_delay = num_delay + cnum_delay
   endif

   if(lfact(bcol)%delay_head .ne. -1) then
      ! copy delayed columns from each of the children of node

      ! set up mapping array pos
      nvar = size(nodes(node)%index)
      l_nb = nodes(node)%nb
      do i = 1+l_nb*(local-1),nvar
         j = nodes(node)%index(i)
         pos(j) = i + nd - l_nb*(local-1)
      end do

      ! Iterate over children
      cnode = lfact(bcol)%delay_head
      do while(cnode.ne.-1)
         cnum_delay = nodes(cnode)%num_delay
         if (cnum_delay == 0) then
            cnode = nodes(cnode)%delay_next
            cycle
         endif
         ! determine which block column was the final block column of
         ! the child cnode. en_child is last block in child node.
         en_child = nodes(cnode)%blk_en
         bcol_cnode = blocks(en_child)%bcol
         cnelim = lfact(bcol_cnode)%nelim
         cblkn = lfact(bcol_cnode)%blkn_new
         cblkm = lfact(bcol_cnode)%nrow + &
            (lfact(bcol_cnode)%blkn_new-blocks(lfact(bcol_cnode)%dblk)%blkn)

         ! set cnrow to hold total number of rows in block column
         cnrow = size(lfact(bcol_cnode)%index_new)

         ip = num_delay*blkm_new
         ip1 = (cblkm+1) * cnelim
         do j = 1, cnum_delay
            ! Rows corresponding to delays
            ipd = ip + num_delay
            work(ipd+j:ipd+cnum_delay) = &
               lfact(bcol_cnode)%lcol(ip1+1:ip1+(cnum_delay-j+1))
            ! Rows in child corresponding to expected rows of node
            ip1 = ip1 + (cnum_delay-j+1) + 1
            do i = cnelim+cnum_delay+1, cnrow
               k = lfact(bcol_cnode)%index_new(i)
               jj = pos(k) ! ith row of child maps to row jj in node
               work(ip+jj) = &
                  lfact(bcol_cnode)%lcol(ip1)
               ip1 = ip1 + 1
            end do
            ip = ip + blkm_new
            ip1 = ip1 + cnelim + j-1
         end do

         num_delay = num_delay + cnum_delay
         cnode = nodes(cnode)%delay_next
      end do

      ! reset pos
      do i = 1+l_nb*(local-1),nvar
         j = nodes(node)%index(i)
         pos(j) = 0
      end do

   end if
   ! free up space that was used by lcol
   deallocate(lfact(bcol)%lcol,stat=st)

end subroutine merge_delays

integer function swizzle_sa(sa, blkn, nelim)
   integer, intent(in) :: sa
   integer, intent(in) :: blkn
   integer, intent(in) :: nelim

   swizzle_sa = (sa-1) / blkn
   swizzle_sa = nelim *swizzle_sa + 1
end function swizzle_sa

!*************************************************

!
! This routine adds internal update tasks before 
! adding inter-nodal updates.
! It is called after the factor_solve has been performed
! for the block column lcol that dblk belongs to (lcol
! block col. within node).
!
subroutine add_updates_new(stack, map, nodes, blocks, dblk, snode, &
      control, info, st)
   type(taskstack), intent(inout) :: stack ! holds tasks that have to be done.
   integer, dimension(*), intent(inout) :: map ! used for inter-nodal updates.
   type(node_type), dimension(-1:), intent(in) :: nodes
   type(block_type), dimension(*), intent(inout) :: blocks
   integer(long), intent(in) :: dblk ! diagonal block in block col.
   integer, intent(in) :: snode ! node to which block col. belongs
   type(MA86_control), intent(in) :: control
   integer, intent(inout) :: info ! error flag
   integer, intent(inout) :: st ! stat parameter

   integer(long) :: i
   integer(long) :: last ! id of final block that generates an update 
      ! within node
   integer :: j
   integer :: nb ! block size nodes(node)%nb
   integer :: nbcol ! number of block cols in node
   integer :: numcol ! number of cols in node
   integer :: numrow ! number of rows in node
   integer :: sz ! number of blocks in first block col of node
   type(dagtask) :: task

   ! Work out useful information on blocks in snode
   nb = nodes(snode)%nb
   numrow = size(nodes(snode)%index)
   numcol = nodes(snode)%en - nodes(snode)%sa + 1

   !
   ! Determine final block that generates an update within node
   ! (it is the last block in a block row that has diagonal elements)
   !
   ! sz is number of (row) blocks in the first block column of node
   ! nbcol is the number of block cols in node

   sz = (numrow - 1) / nb + 1
   nbcol = (numcol - 1) / nb + 1 
   last = blocks(dblk)%last_blk - (sz - nbcol)

   ! add the UPDATE_INTERNAL tasks (if any) from this block column.

   call zero_task(task)
   task%task_type = TASK_UPDATE_INTERNAL

   do j = dblk+1, last
     task%src1 = j
     task%dest = get_dest_block(blocks(j), blocks(j))
     do i = j, blocks(dblk)%last_blk
       task%src2 = i
       call add_task(stack, task, control, info, st)
       if(info.lt.0) return
       task%dest = task%dest + 1
     end do
   end do

   ! add UPDATE_BETWEEN tasks.

   call add_between_updates_simple(nodes, blocks, dblk, snode, &
      stack, map, control, info, st)

end subroutine add_updates_new

!*************************************************

!
! Add internodal updates
! Do this by moving up the elimination tree comparing row indices each time.
!
subroutine add_between_updates_simple(nodes, blocks, src, snode,  &
      stack, map, control, info, st)
   type(node_type), dimension(-1:), intent(in) :: nodes
   type(block_type), dimension(*), intent(inout) :: blocks
   integer(long), intent(in) :: src ! Id of a block in the
     ! block col where factor_solve just done (can be any block
     ! in the block col ... we in fact pass the diagonal block).
   integer, intent(in) :: snode ! node that src belongs to.
   type(taskstack), intent(inout) :: stack ! holds tasks that have to be done.
   integer, dimension(*), intent(inout) :: map ! Workarray to hold map from row 
     ! indices to block indices in ancestor node. 
   type(MA86_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(inout) :: st

   type(dagtask) :: task ! Used to send a task to the pool
   integer :: a_nb  ! Block size of anode
   integer :: anode ! Ancestor of snode
   integer :: cb    ! Local index of column block in anode
   integer :: cptr  ! Position in snode of the first row 
     ! matching a column of the current block column of anode.
   integer(long) :: d_anode ! id of diagonal block of anode
   integer :: i
   integer :: jb ! Block index in anode
   integer :: jlast ! Last column in the cb-th block column of anode
   integer :: k
   integer :: k1
   integer :: numcol ! number of cols in snode
   integer :: size_anode ! size(nodes(anode)%index)
   integer :: size_snode ! size(nodes(snode)%index)

   logical :: map_done ! True if map has been build for anode

   call zero_task(task)
   task%task_type = TASK_UPDATE_BETWEEN
   task%src1 = src

   ! cache some values in variables
   size_snode = size(nodes(snode)%index)

   anode = nodes(snode)%parent
   numcol = nodes(snode)%en - nodes(snode)%sa + 1
   ! initialise cptr to point to first row below triangular part of snode.
   cptr = 1 + numcol

   do while(anode.gt.0)

      ! Skip columns that come from other children
      do cptr = cptr, size_snode
         if(nodes(snode)%index(cptr).ge.nodes(anode)%sa) exit
      end do
      if(cptr.gt.size_snode) exit ! finished with snode

      map_done = .false. ! We will only build a map when we need it
      a_nb = nodes(anode)%nb

      ! Loop over affected block columns of anode
      bcols: do

         if(nodes(snode)%index(cptr).gt.nodes(anode)%en) exit

         ! compute local index of block column in anode and find the id of 
         ! its diagonal block
         cb = (nodes(snode)%index(cptr) - nodes(anode)%sa)/a_nb + 1
         d_anode = nodes(anode)%blk_sa
         do jb = 2, cb
            d_anode = blocks(d_anode)%last_blk + 1
         end do

         ! Build a map of ancestor's blocks
         if(.not.map_done) then
            ! The indices for each row block in anode are mapped to a local row
            ! block index.
            size_anode = size(nodes(anode)%index)
            jb = 1
            do i = 1, size_anode, a_nb
               do k = i, min(i+a_nb-1, size_anode)
                  k1 = nodes(anode)%index(k)
                  map(k1) = jb
               end do
               jb = jb + 1 
            end do
         endif

         ! Loop over remaining rows of snode, adding update_between
         ! task for each row block in current block col. of anode
         ! corresponding to rows of snode
         jb = -1 
         do i = cptr, size_snode
            k1 = nodes(snode)%index(i)
            k = map(k1)
            if(k.ne.jb) then
               task%dest = d_anode + k - cb
               call add_task(stack, task, control, info, st)
               if(info.lt.0) return
               ! block k in anode dealt with so set jb to avoid reconsideration
               jb = k
            endif
         end do

         ! Move cptr to first row in another block of anode
         jlast = min(nodes(anode)%sa + cb*a_nb - 1, nodes(anode)%en)
         do cptr = cptr, size_snode
            if(nodes(snode)%index(cptr) > jlast) exit
         end do
         if(cptr.gt.size_snode) exit 

      end do bcols

      ! Move up the tree
      anode = nodes(anode)%parent
   end do

end subroutine add_between_updates_simple

!*************************************************
!
! Reduce the given dependency of a block column lcol by 1; 
! if it is then zero add the resultant factor_solve task
! to the local task stack (or, if full, to task pool)
!
subroutine reduce_block_dep(stack, lcol, control, info, st)
   type(taskstack), intent(inout) :: stack ! holds tasks that have to be done.
      ! it is altered if the dependency for block column is reduced to 0.
   type(lfactor), intent(inout) :: lcol ! block column for which
      ! dependency count is to be reduced
   type(MA86_control), intent(in) :: control
   integer, intent(inout) :: info ! error flag
   integer, intent(inout) :: st ! stat parameter

   ! Local
   type(dagtask) :: task ! used to hold the task to be added

   ! acquire the lock for the block column
!$ call omp_set_lock(lcol%lock)

   ! Decrement dep count
   lcol%dep = lcol%dep - 1
   !
   ! If dep count is zero then add the factor_solve task
   !
   if(lcol%dep.ne.0) then
!$    call omp_unset_lock(lcol%lock)
      return
   endif

   lcol%dep = - 1

   call zero_task(task)

   task%task_type = TASK_FACTORIZE_COLUMN

   task%dest = lcol%dblk

   call add_task(stack, task, control, info, st)
   ! If info is non-zero we're about to return anyway

!$ call omp_unset_lock(lcol%lock)

end subroutine reduce_block_dep

!*************************************************

subroutine solve_task_dispatch(matrix_type, nbcol, lfact, stack, blocks, &
      nodes, counts, control, info, st, nrhs, ldr, rhs, total_threads, &
      rhs_local, maxmn, job)
   integer, intent(in) :: matrix_type ! (-4 Hermitian, -5 symmetric)
   integer, intent(in) :: nbcol  ! Size of lfact (no. block cols)
   type(lfactor), dimension(nbcol), intent(inout) :: lfact ! Entries in 
      ! block columns of L
   type(taskstack), intent(inout) :: stack         ! task pool
   type(block_type), dimension(:), intent(inout) :: blocks ! block info
   type(node_type), dimension(-1:), intent(in) :: nodes ! Node info
   type(slv_count_type), dimension(:), intent(inout) :: counts
   type(MA86_control), intent(in) :: control
   integer, intent(out) :: info ! error flag
   integer, intent(out) :: st ! stat parameter
   integer, intent(in) :: nrhs  ! number of right-hand sides (maybe = 0)
   integer, intent(in) :: ldr  ! leading extent of rhs
   complex(wp), intent(inout) :: rhs(ldr*nrhs)  ! On entry holds rhs data. 
      ! Overwritten by partial solution (forward substitution performed).
   integer, intent(in) :: total_threads  ! number of threads
   complex(wp), intent(inout) :: rhs_local(ldr*nrhs,0:total_threads-1)  ! Local
      ! right-hand side arrays (one for each thread).
   integer, intent(in) :: maxmn ! max block dimension (=maxmn)
   integer, intent(in) :: job

   integer :: bcol ! block col that task%dest belongs to
   integer(long) :: blk ! task%dest
   integer :: col ! global index of fist col. in dcol
   integer :: i ! loop variable
   integer :: j ! loop variable
   integer :: m ! set to blocks(blk)%blkm (number of rows in blk)
   integer :: n ! set to blocks(blk)%blkn (number of columns in blk and dblk)
   integer :: node
   integer :: r ! current right hand side, 0 indexed
   integer :: sa ! set to blocks(blk)%sa 
   type(dagtask) :: task
 !%%%  integer :: t_start, t_end
   integer :: this_thread

   complex(wp), dimension(:), allocatable :: evars ! elimination variables
   complex(wp), dimension(:), allocatable :: xlocal ! update_buffer workspace

   integer :: idx_offset, nelim, offset_local, offset_rhs, s

   this_thread = 0
!$ this_thread = omp_get_thread_num()

   ! Initialize things
   info = 0; st = 0        ! By default everything went OK
   call zero_task(task)    ! Set everything to zero
   task%task_type = TASK_NONE   ! Needs to be set to prevent 
                                ! dispatched incrementing

   allocate(xlocal(maxmn*nrhs), evars(ldr*nrhs), stat=st)
   if (st.ne.0) then
      info = MA86_ERROR_ALLOCATION
      call MA86_print_flag(info, control, context='MA86_solve',st=st)
      return
   end if

   ! Mark thread as active
!$ call omp_set_lock(stack%lock)
   stack%active = stack%active + 1
!$ call omp_unset_lock(stack%lock)

   !
   ! Main loop
   !
   do
      !
      ! Retrieve next task to perform and wait till all is ready
      !
!%%%  if(control%unit_log.gt.0) call system_clock(t_start)

      call get_task(stack, task, control, info, st)
      if(info.lt.0) return

      !write(6, "(i3,2a)") &
      !   this_thread, " got task ",print_job(task) 

!%%%  if(control%unit_log.gt.0) then
!%%%     call system_clock(t_end)
!%%%     call log_task(control, this_thread, t_start, t_end, "GT", &
!%%%        int(task%task_type,long))
!%%%  endif

!$OMP FLUSH(stack,rhs)

      !
      ! Perform task and update dependencies
      !
      select case(task%task_type)
      case(TASK_DONE)
         exit

      case(TASK_NONE) ! Job not finished but no tasks available, spin.
         cycle

      case(TASK_SLV_FSLV) ! Forwards solve with node
         node = task%dest

         ! set blk to be first block in node
         blk = nodes(node)%blk_sa

         ! each pass of this loop deals with a block column in node
         do while(blk.le.nodes(node)%blk_en)
            ! Establish variables describing block column
            bcol  = blocks(blk)%bcol
            col   = lfact(bcol)%col

            call fwd_solve_bcol(matrix_type, bcol, blk, col, lfact,        &
               blocks, nrhs, rhs, ldr, evars, ldr, xlocal, rhs_local,      &
               control, this_thread, total_threads)
         end do ! block column

!$OMP FLUSH

         if(nodes(node)%parent .gt. 0) then
!$          call omp_set_lock(counts(nodes(node)%parent)%lock)
            counts(nodes(node)%parent)%dep = &
               counts(nodes(node)%parent)%dep - 1
            if(counts(nodes(node)%parent)%dep .eq. 0) then
               task%task_type = TASK_SLV_FSLV
               task%dest = nodes(node)%parent
               call add_task(stack, task, control, info, st)
               if(info.lt.0) return
            endif
!$          call omp_unset_lock(counts(nodes(node)%parent)%lock)
         endif

         ! If we are a root node and need to do a backwards solve, add the
         ! relevant task
         if(job.ne.SOLVE_JOB_FWD .and. nodes(node)%parent.eq.-1) then
            task%task_type = TASK_SLV_BSLV
            task%dest = node
            call add_task(stack, task, control, info, st)
            if(info.lt.0) return
         endif

      case(TASK_SLV_BSLV) ! Backward solve with node
         node = task%dest

         ! Loop over block columns
         blk = nodes(node)%blk_en
         do while(blk.ge.nodes(node)%blk_sa)
            ! Establish variables describing block column
            bcol     = blocks(blk)%bcol
            n        = lfact(bcol)%blkn_new
            nelim    = lfact(bcol)%nelim
            s        = n - blocks(blk)%blkn
            col      = lfact(bcol)%col
         
            ! Zero relevant part of rhs_local
            do r = 1, nrhs
               offset_local = (r-1)*ldr + col - 1
               do i = 1, nelim
                  rhs_local(i+offset_local,0) = czero
               end do
            end do

            ! Loop over blocks in column
            do blk = blk, blocks(blk)%dblk+1, -1
               m        = blocks(blk)%blkm
               sa       = swizzle_sa(blocks(blk)%sa_new, n, nelim)

               idx_offset = blk - blocks(blk)%dblk ! block row, 0 indexed
               idx_offset = 1 + s + idx_offset * nodes(node)%nb

               call slv_bwd_update(m, nelim, col, idx_offset, &
                  lfact(bcol)%index_new, &
                  lfact(bcol)%lcol(sa:sa+m*nelim-1), &
                  nelim, nrhs, rhs, rhs_local(:,0), ldr, xlocal, control, &
                  blocks(blk)%id)
            end do

            ! Determine properties of diagonal block
            m        = blocks(blk)%blkm + s 
            col      = lfact(bcol)%col

            ! Perform any retangular update from diagonal block
            if(m.gt.nelim) then
               sa = 1 + nelim*nelim
               call slv_bwd_update(m-nelim, nelim, col, nelim+1,        &
                  lfact(bcol)%index_new,                                &
                  lfact(bcol)%lcol(sa:sa+(m-nelim)*nelim-1), nelim,     &
                  nrhs, rhs, rhs_local(:,0), ldr, xlocal, control,      &
                  blocks(blk)%id)
            endif

            ! Add in b and perform diagonal solve (if required)
            if(job.eq.SOLVE_JOB_BWD) then
               ! No diagonal solve
               do r = 1, nrhs
                  offset_rhs = (r-1)*ldr
                  offset_local = offset_rhs + col - 1
                  do i = 1, nelim
                     j = lfact(bcol)%index_new(i)
                     rhs_local(i+offset_local,0) = &
                        rhs_local(i+offset_local,0) + rhs(j+offset_rhs)
                  end do
               end do
            else
               ! Add and do diagonal solve simultaneously
               do r = 1, nrhs
                  offset_rhs = (r-1)*ldr
                  offset_local = offset_rhs + col - 1
                  if(nelim.gt.0) &
                     call solve_add_D(nelim, lfact(bcol)%d, rhs(offset_rhs+1), &
                        rhs_local(offset_local+1,0), lfact(bcol)%index_new, &
                        matrix_type)
               end do
            endif

            ! Perform triangular solve
            call slv_solve_bwd(nelim, nelim, col, &
               lfact(bcol)%lcol(1:nelim*nelim),   &
               nrhs, rhs_local(:,0), ldr, control, blocks(blk)%id)

            ! Permute rhs_local into rhs
            do r = 1, nrhs
               offset_rhs = (r-1)*ldr
               offset_local = offset_rhs + col - 1
               do i = 1, nelim
                  j = lfact(bcol)%index_new(i)
                  rhs(j+offset_rhs) = rhs_local(i+offset_local,0)
               end do
            end do

            blk = blk - 1
         end do ! loop over block columns

!$OMP FLUSH

         ! Add bwd solves for child nodes
         task%task_type = TASK_SLV_BSLV
         do j = 1, nodes(node)%nchild
            task%dest = nodes(node)%child(j)
            call add_task(stack, task, control, info, st)
            if(info.lt.0) return
         end do
      case default
         info = MA86_ERROR_UNKNOWN
         call MA86_print_flag(info, control, context='MA86_factor')
         !if(control%diagnostics_level.ge.0 .and. control%unit_error.ge.0) &
         !   write(control%unit_error, "(/a,i3,a,i8)") &
         !   " MA86_factor: Internal Error ", info, &
         !   " Unknown task type encountered type = ", task%task_type
         exit
      end select

   end do
end subroutine solve_task_dispatch

!*************************************************

subroutine fwd_solve_bcol(matrix_type,  &
      bcol, blk, col, lfact, blocks, nrhs, rhs, ldr, &
      evars, lde, xlocal, rhs_local, control, this_thread, total_threads)
   integer, intent(in) :: matrix_type ! (-4 Hermitian, -5 symmetric)
   integer, intent(in) :: bcol ! block column to perform fwd solve on
   integer(long), intent(inout) :: blk ! on entry: first block of bcol,
      ! on exit: (final block of bcol) + 1
   integer, intent(in) :: col ! offset into rhs workspace
   type(lfactor), dimension(*), intent(inout) :: lfact ! Entries in 
      ! block columns of L
   type(block_type), dimension(:), intent(inout) :: blocks ! block info
   integer, intent(in) :: nrhs ! number of right hand sides
   integer, intent(in) :: ldr ! leading dimension of rhs
   complex(wp), dimension(nrhs*ldr), intent(inout) :: rhs ! right hand side
      ! vector
   integer, intent(in) :: lde ! leading dimension of evars
   complex(wp), dimension(nrhs*lde), intent(inout) :: evars ! workspace to
      ! store contigous version of rhs in (potentially shared)
   complex(wp), dimension(*), intent(out) :: xlocal ! local workspace. Must be
      ! of size at least nrhs*nelim
   integer, intent(in) :: total_threads ! total number of threads running in
      ! parallel
   complex(wp), dimension(nrhs*ldr, 0:total_threads-1), intent(inout) :: &
      rhs_local ! per thread accumulation vectors
   type(MA86_control), intent(in) :: control
   integer, intent(in) :: this_thread

   character(len=1) :: transa ! 'C' Hermitian, 'T' symmetric
   integer :: i ! temporary variable
   integer :: j ! temporary variable
   integer :: idx_offset ! offset into index
   integer :: m ! number of rows in block
   integer :: n ! leading dimension of column
   integer :: nelim ! number of eliminated variables
   integer :: offset_rhs ! offset into rhs
   integer :: offset_local ! offset into local workspace
   integer :: r ! loop index (over right hand sides)
   integer :: s ! number of extra variables in block due to delays
   integer :: sa ! start of block within bcol
   integer :: t ! loop index (over threads)

   transa = 'C'
   if (matrix_type == HSL_MATRIX_CPLX_SYM) transa = 'T'

   ! Establish variables describing block column
   n     = lfact(bcol)%blkn_new
   nelim = lfact(bcol)%nelim
   idx_offset = 1

   do r = 1, nrhs
      offset_rhs = (r-1)*ldr
      offset_local = (r-1)*lde + col - 1
      do i = 1, nelim
         j = lfact(bcol)%index_new(i)
         evars(i+offset_local) = rhs(j+offset_rhs)
      end do
      do t = 0, total_threads-1
         do i = 1, nelim
            j = lfact(bcol)%index_new(i)
            evars(i+offset_local) = evars(i+offset_local) + &
               rhs_local(j+offset_rhs, t)
         end do
      end do
   end do

   ! s is the number of extra rows/cols in bcol because of delays
   s     = n - blocks(blk)%blkn
   m     = blocks(blk)%blkm + s
   call slv_solve_fwd(transa, nelim, nelim, col, &
      lfact(bcol)%lcol(1:nelim*nelim), &
      nrhs, evars, lde, control, blocks(blk)%id)
   idx_offset = idx_offset + nelim

   ! Permute rhslocal back in to rhs (rhslocal is still used later)
   do r = 1, nrhs
      offset_rhs = (r-1)*ldr
      offset_local = (r-1)*lde + col - 1
      do i = 1, nelim
         j = lfact(bcol)%index_new(i)
         rhs(j+offset_rhs) = evars(i+offset_local)
      end do
   end do

   m = m - nelim
   if(m.gt.0) then
      sa = 1 + nelim*nelim
      call slv_fwd_update(transa, m, nelim, col, idx_offset,      &
         lfact(bcol)%index_new,                                   &
         lfact(bcol)%lcol(sa:sa+m*nelim-1), nelim, nrhs,          &
         rhs_local(:,this_thread), ldr, evars, lde, xlocal,       &
         control, blocks(blk)%id)
      idx_offset = idx_offset + m
   endif

   ! loop over remaining blocks in the block column bcol
   do blk = blocks(blk)%dblk+1, blocks(blk)%last_blk
      m    = blocks(blk)%blkm
      sa   = swizzle_sa(blocks(blk)%sa_new, n, nelim)

      call slv_fwd_update(transa, m, nelim, col, idx_offset,      &
         lfact(bcol)%index_new,                                   &
         lfact(bcol)%lcol(sa:sa+m*nelim-1), nelim, nrhs,          &
         rhs_local(:,this_thread), ldr, evars, lde, xlocal,       &
         control, blocks(blk)%id)
      idx_offset = idx_offset + m
   end do
end subroutine fwd_solve_bcol

!*************************************************

!
! Returns the destination block of an internal update task.
! Called by add_updates.
!
integer(long) function get_dest_block(src1, src2)
   type(block_type), intent(in) :: src1
   type(block_type), intent(in) :: src2

   integer(long) :: i
   integer :: sz

   ! Move to diagonal block of target column
   ! sz is the number of (row) blocks in src1
   sz = src1%last_blk - src1%dblk + 1 
   get_dest_block = src1%dblk
   do i = src1%dblk+1, src1%id
      get_dest_block = get_dest_block + sz
      sz = sz - 1
   end do

   ! Move to relevant row block in target col.
   get_dest_block = get_dest_block + src2%id - src1%id

end function get_dest_block

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Numerical block operation routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! performs factorization of trapezoidal matrix.
! based on hsl_ma64 version 6.0.0 (ported 18th January 2011)
!
subroutine factor_solve_block(matrix_type, n, p, nb, s, a, la, buf, d, perm, &
     q, tstats, control, flag)
   integer, intent(in) :: matrix_type
   integer, intent(in) :: n ! number of rows in trapeziodal matrix
   integer, intent(in) :: p ! number of cols in trapeziodal matrix
   integer, intent (in) :: nb ! Block size
   integer, intent(in) :: s ! number of cols in trapeziodal matrix that
     ! are from delayed pivots. These cols are searched last for pivots.
   integer :: la  ! size of array lcol
   complex(wp), intent(inout) :: a(la) ! holds trapezoidal matrix 
     ! to be factorized. Holds the factorization on exit.
      ! Delayed columns are permuted to the final columns. Held rowwise
      ! at all times.
   complex (wp) :: buf(nb*n) ! Work array
   complex (wp), intent (out) :: d(2*p) ! d(1:2*q) is set to hold the inverse
      ! of D, except that zero diagonal blocks (pivots) are not inverted.
      ! Diagonal entries are in d(1:2*q-1:2) and entries to the right of
      ! the diagonal are in d(2:2*q-2:2). d(2*q) is set to zero.
   integer, intent(inout) :: perm(:) ! perm array. Input ignored.
      ! On return, for i = 1, 2, ..., p, perm(i) is set to the index
      ! of the row of A that is permuted to row i.
   integer, intent(out) :: q ! number of pivots chosen
   type(thread_info), intent(inout) :: tstats
   type(MA86_control), intent(in) :: control
   integer, intent(out) :: flag

   !     .. Local Scalars ..
   real(wp) abs_amax ! Largest absolute value to left of diagonal in row m.
   complex(wp) amax ! Entry of largest absolute value to left of diagonal in
      ! row m.
   real(wp) amax2 ! Second largest absolute value to left of diagonal in row m
   real(wp) amaxm1 ! Largest absolute value below diagonal in column m, not
      ! including row m+1
   real(wp) amaxb ! Largest absolute value below diagonal in column m
   real(wp) amaxt ! Largest absolute value in column t.
   real(wp) amaxt_cache ! Largest absolute value in column m-1 (used as cache).
   integer  deti   ! Determinant held as detr*radix**deti
   complex(wp) detpiv ! Determinant of candidate pivot is detpiv/detscale.
   complex(wp) detpiv0 ! First term in calculation of detpiv.
   complex(wp) detpiv1 ! Second term in calculation of detpiv.
   complex(wp) detpiv2 ! Value of detpiv for best candidate 2x2 pivot.
   real(wp) detr   ! Scaled determinant
   real(wp) detscale ! Inverse of the largest entry in the candidate pivot.
   real(wp) detscale2 ! Value of detscale for best candidate 2x2 pivot.
   integer i ! Row index
   integer j ! Column index
   integer(long) k  ! Position in a
   integer(long) kkj ! Position of diagonal of column j
   integer(long) kkt ! Position of diagonal of column t
   integer(long) kkr1 ! Position of diagonal of column r+1
   integer(long) kkq !  Start of the pivotal block column
   integer(long) kkq1 ! Position of diagonal of column q+1
   integer(long) kkm ! Position of diagonal of column m
   integer(long) kq1 ! Position in buf of diagonal of column q+1
   integer lq ! Height of the pivot block
   integer m  ! Column searched most recently
   integer mbest ! Column with best relative pivot value for a 1x1 pivot
   integer mbest1, mbest2 ! Pair with best relative pivot value for a 2x2 pivot
   integer mdummy  ! Loop execution count
   integer mlast ! Last column of the block in which m appears
   integer nbi ! Inner block size
   integer pivsiz ! Size of the chosen pivot, 0 if none chosen, or -1
      ! if column is essentially zero
   complex(wp) pivval ! temporary variable for storing pivotal value
   integer qlast ! Last column of the inner block containing column q+1
   integer r ! Number of pivot operations applied to columns q+1:p
   integer rm ! Number of pivot operations applied to the (outer) block
              ! containing column m
   real(wp) rmax ! Largest entry in column m
   real(wp) rmax2 ! Largest entry in column m outwith rows m,m-1.
   integer t ! Candidate 2x2 pivot is in columns t and m
   real(wp) u ! Relative pivot threshold
   real(wp) ubest1 ! Relative pivot value of best candidate 1x1 pivot
   real(wp) ubest2 ! Relative pivot value of best candidate 2x2 pivot
   real(wp) urel ! Relative pivot value
   real(wp) umin ! Minimum relative pivot threshold
   integer :: nzero ! Number of zero pivots since last non-zero pivot

   flag = 0
   !if (n < 0) then
   !   flag = -1
   !else if (p < 0) then
   !   flag = -2
   !else if (p > n) then
   !   flag = -3
   !else if (nbi <= 1) then
   !   flag = -4
   !else if ( la < min(n*int(n,long),(n*(n+nb+1_long))/2) ) then
   !   flag = -7
   !else if (control%static < control%small .and. control%static/=rzero) then
   !   flag = -10
   !else if (mod(nb,nbi) /= 0) then
   !   flag = -12
   !end if

   !info%detlog = rzero
   !info%num_neg = 0
   !info%num_nothresh = 0
   !info%num_perturbed = 0
   !info%num_zero = 0
   !info%num_2x2 = 0

   ! u is reset to control%u for each block column, but may then be relaxed to
   ! umin on a block column by block column basis
   u = min(max(control%u,rzero),rone)

   nbi = control%nbi
   if(nbi.le.1) nbi = nbi_default ! Ensure positive

   umin = min(max(control%umin,rzero),u)
   if (p==n) umin = min(umin,0.5_wp)
   q = 0
   if (flag/=0 .or. p==0) return
   deti = 0
   detr = 1.0_wp
   lq = n
   kkq1 = 1
   kkq = 1
   m = p
   qlast = min(nbi,p)
   ! m is updated at the start of the main loop so initializing it to p causes
   ! it to have the value 1 during the first execution of that loop.

   do j = 1, p
      perm(j) = j
   end do

   !print *
   !print *, "======================="
   !print *, "factor n, p, s = ", n, p, s
   !print *, "a input = "
   !print "(11es12.4)", a

   if (s>0 .and. s<p) then
      i = min(s,p-s)
      ! Make first i columns be last
      do j = 1,i
         call swap_cols(matrix_type,n,nb,q,1,0,a,buf,perm,j,p+1-j)
      end do
      !! Ensure that the first index of a 2x2 pair is labelled
      !do j = 2,i
      !   if (perm(j)<0) then
      !      perm(j) = -perm(j)
      !      perm(j-1) = -perm(j-1)
      !   end if
      !end do
      !do j = p+2-i,p
      !   if (perm(j)<0) then
      !      perm(j) = -perm(j)
      !      perm(j-1) = -perm(j-1)
      !   end if
      !end do
   end if

   nzero = 0 ! number of zero columns encountered since last proper pivot
   pivot: do
      ! Perform a pivotal operation
      pivsiz = 0
      ubest1 = rzero
      ubest2 = rzero
      sweep: do mdummy = 1, p-q ! Look for a pivotal column or pair of columns
         ! Update m and the scalars associated with column m
         m = m+1
         if (m>p) then
            ! Go back to column q+1
            m = q+1
            kkm = kkq1
            r = q
            rm = q
         else if (m<1+nb) then
            ! Within the current block column
            kkm = kkm + lq + 1
         end if

         ! Update column m
         kkr1 = kkq1 + (rm-q)*(lq+1_long)
         k = kkr1+m-rm-1
         if(q>rm) then
            call cgemv('NoTrans',n-m+1,q-rm,cone,a(k),lq, &
               buf(n*(rm+0_long)+m),n,cone,a(kkm),1)
         end if
         if(matrix_type .eq. HSL_MATRIX_CPLX_HERM_INDEF) &
            a(kkm) = real(a(kkm),wp)

         ! Find largest and second largest entry to left of diagonal in row m.
         j = q + 1
         k = kkq1 + m - j - n
         amax = rzero
         abs_amax = rzero
         amax2 = rzero
         t = 0
         if (j<m) then
            t = j
            kkt = kkq1
         end if
         do j = j, m-1
            k = k + lq
            if (abs(a(k))>abs_amax) then
               t = j
               amax2 = abs(amax)
               amax = a(k)
               abs_amax = abs(a(k))
               kkt = k - (m-j)
            else
               !amax2 = max(abs(a(k)),amax2)
               if(abs(a(k)).gt.amax2) amax2 = abs(a(k))
            end if
         end do

         ! Now calculate largest entry below the diagonal of column m.
         amaxm1 = rzero
         amaxb = rzero
         !do i = m+1,n
         !   amaxb = max(abs(a(kkm+i-m)),amaxb)
         !end do
         if(n-m.ge.1) amaxb = abs(a(kkm+1))
         do k = kkm+2,kkm+n-m
            if(abs(a(k)).gt.amaxm1) amaxm1 = abs(a(k))
         end do
         if(amaxm1.gt.amaxb) amaxb = amaxm1

         ! Now calculate largest entry in the whole of column m and make sure
         ! that it is neither small nor infinity.
         !rmax = max(abs_amax,abs(a(kkm)),amaxb)
         rmax = abs_amax
         if(abs(a(kkm)).gt.rmax) rmax = abs(a(kkm))
         if(amaxb.gt.rmax) rmax = amaxb
         if (rmax<=control%small) then
            ! All entries of the column are small
            a(kkm) = czero
            tstats%num_zero_pivots = tstats%num_zero_pivots + 1
            pivsiz = -1
            perm(m) = abs(perm(m))
            exit sweep
         else if (rmax > huge(rzero)) then
            ! There is an infinity in the column
            flag = -13
            return
         end if

         ! Calculate the relative pivot value and see if it is the best so far
         if (abs(a(kkm))>control%small) then
            urel = abs(a(kkm))/rmax
         else
            urel = rzero
         end if
         if (urel >= ubest1) then
            ubest1 = urel
            mbest = m
         end if

         ! If there is a candidate 2x2 pivot, try it.
         ! Look for a 1x1 pivot only if the 2x2 pivot is unacceptable.
         tgt0: if (t>0) then
            if( min(abs(a(kkm)),abs(a(kkt))).gt.control%small .or. &
                  abs_amax.ge.control%small) then

               ! Store value of the largest entry in whole of column m outwith
               ! rows m and t
               rmax2 = max(amax2,amaxb)

               if (matrix_type == HSL_MATRIX_CPLX_HERM_INDEF) &
                  a(kkt) = real(a(kkt),wp)
               detscale = rone/max( abs(a(kkm)), abs(a(kkt)), abs_amax )
               if (matrix_type == HSL_MATRIX_CPLX_HERM_INDEF) then
                  detpiv1 =  (abs_amax*detscale)*abs_amax
               else if (matrix_type == HSL_MATRIX_CPLX_SYM) then
                  detpiv1 =  (amax*detscale)*amax
               end if
               detpiv0 = a(kkm)*detscale*a(kkt)
               detpiv = detpiv0 - detpiv1
               ! Make sure that the 2x2 pivot is not singular and that there is
               ! little cancellation in calculating detpiv. Bearing in mind the
               ! way detscale is calculated, if the largest entry of the matrix
               ! is chosen as pivot, the one entry of the reduced matrix has
               ! absolute value abs(detpiv).
               !print *, "test 2x2 = ", detpiv, detpiv0, detpiv1
               left2x2:if (abs(detpiv)> &
                     max(control%small,abs(detpiv0)/2,abs(detpiv1)/2)) then

                  ! Find largest entry in column t outwith rows m and t
                  if(t.eq.m-1 .and. amaxt_cache.ne.-1) then
                     ! Use cached answer from scan of previous column
                     amaxt = amaxt_cache
                  else
                     amaxt = rzero
                     j = q + 1
                     k = kkq1 + t - j - lq
                     do j = q+1, t-1
                        k = k + lq
                        amaxt = max(abs(a(k)),amaxt)
                     end do
                     k = k + lq
                     do i = t+1,m-1
                        amaxt = max(abs(a(k+i-t)),amaxt)
                     end do
                     do i = m+1,n
                        amaxt = max(abs(a(k+i-t)),amaxt)
                     end do
                  endif

                  ! OK as 2x2 pivot if all entries in the rest of the columns
                  ! are small
                  if(max(rmax2,amaxt)<=control%small) then
                     pivsiz = 2
                     exit sweep
                  end if

                  ! Calculate the relative pivot value (as 2x2 pivot)
                  urel = abs(detpiv)/max( &
                     abs(a(kkm)*detscale)*amaxt+(abs_amax*detscale)*rmax2, &
                     abs(a(kkt)*detscale)*rmax2+(abs_amax*detscale)*amaxt )
                  !print *, "urel2x2 = ", urel

                  ! OK as 2x2 pivot if relative pivot value is big enough
                  if (urel>u) then
                     pivsiz = 2
                     tstats%usmall = min(urel,tstats%usmall)
                     exit sweep
                  end if

                  ! If this has the best relative pivot value so far, record
                  ! this
                  if(urel>ubest2)then
                     ubest2 = urel
                     detpiv2 = detpiv
                     detscale2 = detscale
                     mbest1 = m
                     mbest2 = t
                     !print *, "store 2x2", m, t, urel
                  end if
               end if left2x2
            end if
         end if tgt0
         amaxt_cache = max(abs_amax, amaxm1)

         ! If 2x2 pivot rejected or only one column left, take best 1x1 pivot
         ! if it is OK.
         if(t>0 .or. m==p) then
            !print *, "   test 1x1 ubest1 = ", ubest1
            if (ubest1>u) then
               pivsiz = 1
               tstats%usmall = min(ubest1,tstats%usmall)
               if (mbest/=m) &
                  call swap_cols(matrix_type,n,nb,q,1,rm,a,buf,perm,mbest,m)
               exit sweep
            end if
         end if
      end do sweep

      !
      ! At this stage following variables should have good values:
      ! m - pivot column (wil be swapped to posn q+1 for 1x1 or q+2 for 2x2)
      ! q - number of pivots already performed
      ! kkq1 - posn in a of diagonal of col q+1
      ! t - other pivot column for 2x2 (t<m) (swapped to posn q+q
      ! detscale & detpiv - used for stats
      !

      pivsiz0: if (pivsiz==0) then
         !print *, "pivsize0 dropout", ubest1, ubest2, umin
         ! No pivot found in search of all available columns
         ! Since all the columns have been updated, revise m to q+1
         m = q+1
         kkm = kkq1
         rm = q
         r = q
         ! Perform relaxed pivoting if the best pivot is good enough
         if (max(ubest1,ubest2)>=umin) then
            if (ubest1>=ubest2) then
               ! Accept 1x1 pivot
               u = min(u,ubest1)
               tstats%usmall = min(ubest1,tstats%usmall)
               pivsiz = 1
               if (mbest/=m) &
                  call swap_cols(matrix_type,n,nb,q,1,rm,a,buf,perm,m,mbest)
            else
               ! Accept 2x2 pivot
               u = min(u,ubest2)
               tstats%usmall = min(ubest2,tstats%usmall)
               pivsiz = 2
               ! Revise m to q+1
               m = q+2
               kkm = kkq1+lq+1
               detpiv = detpiv2
               detscale = detscale2
               t = min(mbest1,mbest2)
               if (t/=q+1) &
                  call swap_cols(matrix_type,n,nb,q,1,rm,a,buf,perm,q+1,t)
               t = max(mbest1,mbest2)
               if (t/=q+2) &
                  call swap_cols(matrix_type,n,nb,q,1,rm,a,buf,perm,q+2,t)
               t = q+1
            end if

         else if (control%static>0) then
            ! Perform static pivoting if this has been requested
            if (mbest/=m) &
               call swap_cols(matrix_type,n,nb,q,1,rm,a,buf,perm,m,mbest)
            pivsiz = 1
            tstats%num_nothresh = tstats%num_nothresh + 1
            if (abs(a(kkm)) < control%static) then
               if (matrix_type == HSL_MATRIX_CPLX_HERM_INDEF) then
                  a(kkm) = sign(control%static,real(a(kkm),wp))
               else
                  if(abs(a(kkm))>tiny(rone)) then
                     a(kkm) = control%static*(a(kkm)/abs(a(kkm)))
                  else
                     a(kkm) = control%static
                  end if
               end if
               tstats%num_perturbed = tstats%num_perturbed + 1
               tstats%usmall = -rone
            else
               tstats%usmall = min(ubest1,tstats%usmall)
            end if
         end if
      end if pivsiz0

      pivsizes: if (pivsiz==1) then
         nzero = 0 ! pivot found, reset recent zero count
         amaxt_cache = -1 ! mark cache as invalid
         ! Swap columns if m not q+1
         if (q+1/=m) call swap_cols(matrix_type,n,nb,q,1,rm,a,buf,perm,q+1,m)
         ! We store D**-1.
         !print *, "1x1 pivot ", a(kkq1), "(col ", m, ")"
         d(2*q+1) = cone/a(kkq1)
         d(2*q+2) = czero
         ! Update tstats
         !tstats%detlog = tstats%detlog + log(abs(a(kkq1)))
         deti = deti + exponent(detr) + exponent(abs(a(kkq1)))
         detr = fraction(detr)*fraction(abs(a(kkq1)))
         if (matrix_type == HSL_MATRIX_CPLX_HERM_INDEF) then
            if (real(a(kkq1))<rzero) tstats%num_neg = tstats%num_neg + 1
         else
            tstats%detarg = tstats%detarg*(a(kkq1)/abs(a(kkq1)))
         end if
         ! Store L in A and -LD or -conjg(LD) in buf
         kq1 = q+1+(q+0_long)*n
         buf(kq1) = -a(kkq1)
         a(kkq1) = cone
         if( matrix_type==HSL_MATRIX_CPLX_SYM) then
            !do i = 1, n-q-1
            !   buf(kq1+i) = -a(kkq1+i)
            !   a(kkq1+i) = d(2*q+1)*a(kkq1+i)
            !end do
            pivval = d(2*q+1)
            do i = kq1+1, kq1+n-q-1
               buf(i) = -a(i)
               a(i) = pivval*a(i)
            end do
         else
            !do i = 1, n-q-1
            !   buf(kq1+i) = -conjg(a(kkq1+i))
            !   a(kkq1+i) = d(2*q+1)*a(kkq1+i)
            !end do
            pivval = d(2*q+1)
            do i = kq1+1, kq1+n-q-1
               buf(i) = -conjg(a(i))
               a(i) = pivval*a(i)
            end do
         end if
         ! Update columns q+2 to m
         kkj = kkq1
         do j = q+2, min(nb,m)
            kkj = kkj + lq + 1
            call caxpy(n-j+1,buf(kq1+j-q-1),a(kkq1+j-q-1),1,a(kkj),1)
         end do
         ! Update q and kkq1
         kkq1 = kkq1 + lq + 1
         q = q+1

      else if (pivsiz==2) then pivsizes
         nzero = 0 ! pivot found, reset recent zero count
         amaxt_cache = -1 ! mark cache as invalid
         ! Swap columns unless t==q+1 and m==q+2
         if (q+2/=m) then
            if (q+1/=t) call swap_cols(matrix_type,n,nb,q,1,rm,a,buf,perm,q+1,t)
            call swap_cols(matrix_type,n,nb,q,1,rm,a,buf,perm,q+2,m)
         end if
         ! We store D**-1
         k = kkq1 + lq + 1
         d(2*q+1) = (a(k)*detscale)/detpiv
         d(2*q+3) = (a(kkq1)*detscale)/detpiv
         d(2*q+2) = -(a(kkq1+1)*detscale)/detpiv
         d(2*q+4) = czero
         !print *, "2x2 pivot ", d(2*q+1:2*q+3)
         !print *, "   (cols ", t, m, ")"
         !print *, "   [vals ", a(k), a(kkq1), a(kkq1+1), &
         !   "dp", detscale, detpiv, "]"
         ! Update info
         tstats%num_two = tstats%num_two + 1
         !tstats%detlog = tstats%detlog + log(abs(detpiv)) - log(detscale)
         deti = deti + exponent(detr) - exponent(detscale)
         detr = fraction(detr)*abs(detpiv)/fraction(detscale)
         if (matrix_type == HSL_MATRIX_CPLX_HERM_INDEF) then
            if (real(detpiv)<rzero) then
               tstats%num_neg = tstats%num_neg + 1
            else if (real(a(kkq1)+a(k))<rzero) then
               tstats%num_neg = tstats%num_neg + 2
            end if
         else
            tstats%detarg = tstats%detarg*(detpiv/abs(detpiv))
         end if
         ! Store L in A and -LD or -conjg(LD) in buf
         kq1 = q+1+(q+0_long)*n
         buf(kq1) = -a(kkq1)
         if(matrix_type==HSL_MATRIX_CPLX_HERM_INDEF) then
            buf(kq1+1) = -conjg(a(kkq1+1))
         else
            buf(kq1+1) = -a(kkq1+1)
         end if
         buf(kq1+n) = -a(kkq1+1)
         buf(kq1+n+1) = -a(k)
         a(kkq1) = cone
         a(k) = cone
         a(kkq1+1) = czero
         if(matrix_type==HSL_MATRIX_CPLX_HERM_INDEF) then
            do i = 2, n-q-1
               buf(kq1+i  ) = -conjg(a(kkq1+i))
               buf(kq1+i+n) = -conjg(a(k+i-1))
               a(kkq1+i) = d(2*q+1)*a(kkq1+i) + d(2*q+2)*a(k+i-1)
               a(k+i-1) = d(2*q+3)*a(k+i-1) - conjg(d(2*q+2)*buf(kq1+i))
            end do
         else
            do i = 2, n-q-1
               buf(kq1+i  ) = -a(kkq1+i)
               buf(kq1+i+n) = -a(k+i-1)
               a(kkq1+i) = d(2*q+1)*a(kkq1+i) + d(2*q+2)*a(k+i-1)
               a(k+i-1) = d(2*q+3)*a(k+i-1) - d(2*q+2)*buf(kq1+i)
            end do
         end if
         ! Update columns q+3 to m
         kkj = k + lq + 1
         j = max(q+3,1)
         if(m>=j) then
            call cgemm('n','t',n-j+1,m-j+1,2,cone,a(kkq1+j-q-1), &
               n-((q+1)/nb)*nb,buf(kq1+j-q-1),n,cone,a(kkj),lq)
         end if
         ! Update q and kkq1
         kkq1 = k + lq + 1
         q = q + 2

      else if (pivsiz==-1) then pivsizes
         nzero = nzero + 1 ! count zero pivot
         amaxt_cache = -1 ! mark cache as invalid
         ! Handle a row that is zero
         if (q+1/=m) call swap_cols(matrix_type,n,nb,q,1,rm,a,buf,perm,q+1,m)
         !print *, "fail piv"
         d(2*q+1) = czero
         d(2*q+2) = czero
         ! Store L in A and -LD in buf
         kq1 = q+1+(q+0_long)*n
         buf(kq1) = czero
         a(kkq1) = cone
         do i = 1, n-q-1
            buf(kq1+i) = czero
            a(kkq1+i) = czero
         end do
         kkq1 = kkq1 + lq + 1
         q = q+1
      end if pivsizes

      if (q>=qlast) then
         ! Inner block is complete
         ! Update columns (m+1:min(p,m1+nb-1)) for pivots rm+1:qlast (BLAS3)
         mlast = min(p,nb)
         if (m<mlast) call ma64_update (n,nb,kkq,a,buf,m+1, &
            mlast,rm+1,qlast-nzero)
         rm = qlast
         qlast = min (p,qlast+nbi)
         nzero = 0 ! reset recent zero count
      end if

      if (q == p .or. pivsiz==0) exit pivot
   end do pivot

   ! Update tstats
   if(q<p .and. tstats%num_nothresh==0) &
      tstats%usmall = min(tstats%usmall, max(ubest1,ubest2))
   if (matrix_type == HSL_MATRIX_CPLX_HERM_INDEF) then
      if(tstats%num_zero_pivots/=0) then
         tstats%detsign = 0
         tstats%detlog = rzero
      else
         tstats%detsign = -1
         if (mod(tstats%num_neg,2)==0)tstats%detsign = 1
         tstats%detlog = log(detr)+deti*log(1.0d0*radix(1.0d0))
      end if
   else
      if(tstats%num_zero_pivots/=0) then
         tstats%detlog = rzero
      else
         tstats%detlog = log(detr)+deti*log(1.0d0*radix(1.0d0))
      endif
   end if
   !print *, "exit d = "
   !print "(2es12.4)", d
   !print *, "exit a = "
   !print "(5es12.4)", a
   !print *, "delayed ", p-q, "pivots"

end subroutine factor_solve_block

!*************************************************

!
! Copy trapezoidal block col. held in a from rowwise to colwise, using
! buf as temporary workarray.
!
subroutine row_to_col(n, nb, a_in, a_out)
   integer, intent(in) :: n
   integer, intent(in) :: nb
   complex(wp), dimension(n*(nb+0_long)), intent(inout) :: a_in
   complex(wp), dimension(n*(nb+0_long)), intent(out) :: a_out

   integer, parameter :: xts = 4 ! number of entries to copy in one go
      ! theoretically this should be best if equal to number of entries
      ! in a cache line (I think)

   integer(long) :: i, j, k
   integer :: x, y

   ! copy a_in to a_out, rowwise
   i = 1 ! rowwise
   j = 1 ! colwise
   do x = 1, xts*(nb/xts), xts
      i = 1 + (x-1)*(nb+1) ! diagonal entry for first column of block in a_in
      j = 1 + (x-1)*(n+1)  ! diagonal entry for first column of block in a_out
      do y = x, n
         !call ccopy(xts, a(i), 1, buf(j), n)
         a_out(j:j+(xts-1)*n:n) = a_in(i:i+(xts-1))
         j = j + 1
         i = i + nb
      end do
   end do
   do x = xts*(nb/xts)+1, nb
      i = x
      j = 1+(x-1)*n
      call ccopy(n, a_in(i), nb, a_out(j), 1)
   end do

end subroutine row_to_col

subroutine row_to_col_herm(n, nb, a_in, a_out)
   integer, intent(in) :: n
   integer, intent(in) :: nb
   complex(wp), dimension(n*(nb+0_long)), intent(inout) :: a_in
   complex(wp), dimension(n*(nb+0_long)), intent(out) :: a_out

   integer, parameter :: xts = 4 ! number of entries to copy in one go
      ! theoretically this should be best if equal to number of entries
      ! in a cache line (I think)

   integer(long) :: i, j, k
   integer :: x, y

   ! copy a to buf, rowwise
   i = 1 ! rowwise
   j = 1 ! colwise
   do x = 1, xts*(nb/xts), xts
      i = 1 + (x-1)*(nb+1) ! diagonal entry for first column of block in a
      j = 1 + (x-1)*(n+1)  ! diagonal entry for first column of block in buf
      do y = x, n
         !call ccopy(xts, a(i), 1, buf(j), n)
         a_out(j:j+(xts-1)*n:n) = conjg(a_in(i:i+(xts-1)))
         j = j + 1
         i = i + nb
      end do
   end do
   do x = xts*(nb/xts)+1, nb
      i = x
      j = 1+(x-1)*n
      a_out(j:j+n-1) = conjg(a_in(i:i+(n-1)*nb:nb))
   end do

end subroutine row_to_col_herm

!*************************************************

!
! The eliminated columns (first nb columns) are held rowwise, followed
! by any delayed columns held colwise
!
subroutine col_to_row(n, nb, lfact, colwork, st)
   integer, intent(in) :: n
   integer, intent(in) :: nb
   type(lfactor), intent(inout) :: lfact
   complex(wp), dimension(*), intent(in) :: colwork
   integer, intent(out) :: st

   integer(long) :: i, k
   integer :: j

   st = 0

   ! Reallocate lcol if it is not sufficiently large (ie if there have been
   ! delayed columns inherited by this node)
   if(allocated(lfact%lcol)) then
      if(size(lfact%lcol) < n*lfact%blkn_new) then
         deallocate(lfact%lcol, stat=st)
         allocate(lfact%lcol( n*lfact%blkn_new ), stat=st)
         if(st.ne.0) return
      endif
   else
      allocate(lfact%lcol( n*lfact%blkn_new ), stat=st)
      if(st.ne.0) return
   endif

   ! Copy colwork to lcol, rowwise
   i = 1
   k = 1
   do j = 1, nb
      call ccopy(n-j+1, colwork(k), 1, lfact%lcol(i), nb)
      i = i + nb + 1
      k = k + n + 1
   end do

   ! Copy delays to end of lcol
   lfact%lcol(n*nb+1:n*lfact%blkn_new) = &
      colwork(n*nb+1:n*lfact%blkn_new)

end subroutine col_to_row

!*************************************************

!
! Hermitian case.
! Take trapezoidal block col. held in a, take
! its conjugate of entries in first nelim
! cols and then copy from colwise to rowwise, using
! buf as temporary workarray.
! This means conjg(L) is held by rows ... equivalently, L^H held by cols.
!
subroutine col_to_row_herm(n, nb, lfact, colwork, st)
   integer, intent(in) :: n
   integer, intent(in) :: nb
   type(lfactor), intent(inout) :: lfact
   complex(wp), dimension(*), intent(in) :: colwork
   integer, intent(out) :: st

   integer(long) :: i, k
   integer :: j

   st = 0

   if(allocated(lfact%lcol)) then
      if(size(lfact%lcol) < n*lfact%blkn_new) then
         deallocate(lfact%lcol, stat=st)
         allocate(lfact%lcol( n*lfact%blkn_new ), stat=st)
         if(st.ne.0) return
      endif
   else
      allocate(lfact%lcol( n*lfact%blkn_new ), stat=st)
      if(st.ne.0) return
   endif

   ! Copy colwork to lfact, rowwise
   i = 1
   k = 1
   do j = 1, nb
      lfact%lcol(i:i+(n-j)*nb:nb) = conjg(colwork(k:k+n-j))
      i = i + nb + 1
      k = k + n + 1
   end do

   ! Copy delays to end of lcol
   lfact%lcol(n*nb+1:n*lfact%blkn_new) = &
      colwork(n*nb+1:n*lfact%blkn_new)
end subroutine col_to_row_herm

!*************************************************

!
! TASK_UPDATE_INTERNAL
! A_ik <- A_ik - A_ij A_kj^*
! dest <- dest - src1^* src2
! Remember that the blocks are stored by rows so that
! src1 holds A_ij by rows (or Hermitian case conjg(A_ij) by rows)
! src2 holds A_kj^T by cols (or Hermitian case A_kj^H by cols)
! dest, src1 and src2 all belong to the same node.
!
subroutine update_block_block(transa, m, n, nelim, dest, blk, n1, src1, src2, &
      srcd, ldwork, control)
   character(len=1), intent(in) :: transa ! set to 'C' for Hermitian case
      ! and to 'T' for complex symmetric  
   integer, intent(in) :: m ! number of rows in dest
   integer, intent(in) :: n ! number of columns in dest
   integer, intent(in) :: nelim ! number of columns used for updating
   complex(wp), dimension(*), intent(inout) :: dest ! holds block in L
     ! that is to be updated. 
   type(block_type), intent(inout) :: blk ! destination block  
   integer, intent(in) :: n1 ! leading dim. of src1 and src2
   complex(wp), dimension(*), intent(in) :: src1
   complex(wp), dimension(*), intent(in) :: src2
   complex(wp), dimension(*), intent(in) :: srcd ! inverse of diagonal
   complex(wp), dimension(*), intent(out) :: ldwork ! work array for  
      ! calculating LD
   type(MA86_control), intent(in) :: control

!%%%integer :: t_start, t_end, this_thread
   complex(wp) :: coeff

   call calc_ld(transa, n, nelim, n1, src1, srcd, ldwork)

!%%%if(control%unit_log.gt.0) call system_clock(t_start)

!$ call omp_set_lock(blk%alock)

   coeff = cone
   if(.not.blk%touched) then
      coeff = czero
      blk%touched = .true.
   endif
   call cgemm(transa, 'N', n, m, nelim, -cone, ldwork, nelim, src2, n1,  &
     coeff, dest, n)

!$ call omp_unset_lock(blk%alock)

!%%%if(control%unit_log.gt.0) then
!%%%   call system_clock(t_end)
!%%%   this_thread = 0
!%%%!$    this_thread = omp_get_thread_num()
!%%%   call log_task(control, this_thread, t_start, t_end, "UW")
!%%%endif

end subroutine update_block_block

!*************************************************

!
! If transa = 'T' (symmetric case)
!   takes L (stored by rows) and D^-1, calculates LD (stored by rows)
! If transa = 'C' (Hermitian case)
!   takes conjg(L) (stored by rows) and D^-1, 
!   calculates conjg(L)D (stored by rows)
!
subroutine calc_ld(transa, m, n, ldl, l_mat, d_vec, ld_mat)
   character(len=1), intent(in) :: transa ! set to 'C' for Hermitian case
                                     ! and to 'T' for complex symmetric 
   integer, intent(in) :: m ! number of rows in in block
   integer, intent(in) :: n ! number of cols in in block (=nelim)
   integer, intent(in) :: ldl ! leading dimension of block
   complex(wp), dimension(ldl*m), intent(in) :: l_mat
   complex(wp), dimension(2*n), intent(in) :: d_vec
   complex(wp), dimension(m*n), intent(out) :: ld_mat

   integer :: row, col
   integer :: lptr      ! index into l_mat
   integer :: ldptr     ! index into ld_mat
   complex(wp) :: d1, d2, d3
   complex(wp) :: det
   complex(wp), dimension(2*n) :: dalt ! temporary array storing calculated
      ! values of D

   ! Determine entries of D
   col = 1
   do while(col.le.n)
      if(d_vec(2*col).eq.czero) then 
         ! 1x1 pivot
         dalt(2*col) = czero
         if(d_vec(2*col-1).eq.czero) then
            dalt(2*col-1) = czero
         else
            dalt(2*col-1) = cone/d_vec(2*col-1)
         endif
         col   = col   + 1
      else 
         ! 2x2 pivot
         d1 = d_vec(2*col - 1)
         d2 = d_vec(2*col)
         d3 = d_vec(2*col + 1)
         if (transa == 'C') then
            det = d1*d3 - conjg(d2)*d2
         else
            det = d1*d3 - d2*d2
         end if
         
         dalt(2*col - 1) = d3/det
         dalt(2*col    ) = d2/det
         dalt(2*col + 1) = d1/det

         col = col + 2
      endif
   end do ! cols

   ! Loop over rows
   do row = 1, m
      lptr  = 1 + (row-1)*ldl
      ldptr = 1 + (row-1)*n
      ! Loop over columns
      col = 1
      do while(col.le.n)
         if(dalt(2*col).eq.czero) then 
            ! 1x1 pivot
            ld_mat(ldptr) = l_mat(lptr) * dalt(2*col-1)
            ldptr = ldptr + 1
            lptr  = lptr  + 1
            col   = col   + 1
         else 
            ! 2x2 pivot
            d1 = dalt(2*col - 1)
            d2 = dalt(2*col)
            d3 = dalt(2*col + 1)

            if (transa == 'C') then
               ld_mat(ldptr  ) = d1*l_mat(lptr) - conjg(d2)*l_mat(lptr+1)
            else
               ld_mat(ldptr  ) = d1*l_mat(lptr) - d2*l_mat(lptr+1)
            end if
            ld_mat(ldptr+1) = -d2*l_mat(lptr) + d3*l_mat(lptr+1)

            ldptr = ldptr + 2
            lptr  = lptr  + 2
            col   = col   + 2
         endif
      end do ! cols
   end do ! rows
end subroutine calc_ld

!*************************************************

!
!   Given a destination block dest, update_between performs the update
!                     L_dest <-- L_dest - L_rsrc (L_csrc)^*
!   where L_dest is a submatrix of the block dest of an ancestor
!   of the node snode and L_rsrc and L_csrc are submatrices of contiguous
!   rows that belong to the same block column of snode as the block src
!   (this block col. has local index scol).
!   The first row of L_rsrc is the first row
!   of the block column that corresponds to a row in the block dest and the 
!   last row of L_rsrc is the last row of the block column that corresponds to 
!   a row in the block dest. Similarly, the first/last row of L_csrc is the
!   first/last row of the block column that corresponds to a column in the 
!   block dest. The set of rows and columns of dest thus
!   determine which two sets of contiguous rows in scol are involved.
!   Unless the number of entries updated is very small, use BLAS 3 kernel 
!   gemm by placing its result in a buffer from which we add the 
!   update into the appropriate entries of the
!   destination block dest.
!
subroutine update_between(transa, dest, dnode, n, nelim, sindex, lcol,  &
      lcol_src, d_src, blocks, col_list, row_list, buffer, control, info, &
      st, scol, snode, ldwork)
   character(len=1) :: transa ! set to 'C' for Hermitian and 'T' for symmetric
   integer(long), intent(in) :: dest      ! Destination block
   type(node_type), intent(in) :: dnode   ! Destination node
   integer, intent(in) :: n ! leading dim. the source block column.
   integer, intent(in) :: nelim ! number of eliminations performed
      ! in the source block column (nelim = n if no delayed pivots).
   integer, dimension(:), intent(in) :: sindex ! permuted variable
      ! list for source block column
   type(block_type), dimension(*), intent(inout) :: blocks
   complex(wp), dimension(*), intent(inout) :: lcol ! holds  block
      ! column of factor L that dest belongs to
   complex(wp), dimension(*), intent(in) :: lcol_src ! holds  block
      ! column of factor L that source block belongs to
   complex(wp), dimension(*), intent(in) :: d_src ! holds inverse of
      ! diagonal factor D that source block belongs to
   integer, dimension(*) :: col_list ! size blocks(dest)%blkn
      ! (= number of columns in dest)
   integer, dimension(*) :: row_list ! size blocks(dest)%blkm
      ! (= number of rows in dest)
   complex(wp), dimension(*), intent(out) :: buffer ! work array
   type(MA86_control), intent(in) :: control
   integer, intent(inout) :: info   
   integer, intent(inout) :: st
   !!! following 2 arguments only used in log_task
   integer, intent(in) :: scol ! local index of block column that src 
      ! belongs to within snode
   type(node_type), intent(in) :: snode   ! Source node   
   complex(wp), dimension(:), allocatable, intent(inout) :: ldwork ! work array
      ! for calculating LD

   ! Local scalars
   integer :: cptr ! used in determining the rows that belong to the
      ! source block csrc
   integer :: col_list_sz ! initialise to 0. then increment while
      ! rows involed in csrc are recorded, until
      ! holds the number of columns in dest (= number
      ! of rows in csrc)
   integer :: dcen ! index of end column in dcol
   integer :: dcol ! index of block column that dest belongs to in dnode
   integer :: dcsa ! index of first column in dcol
   integer :: dptr
   integer :: dptr_sa
   integer :: drsa, dren ! point to first and last rows of destination
      ! block dest
   integer :: i
   integer :: j
   integer :: ndiag ! set to int(s1en-s1sa+1) if dest is a
      ! block on diagonal and 0 ow. so is number of triangular rows of update
   integer :: row_list_sz ! initialise to 0. then increment while
      ! rows involed in rsrc are recorded, until
      ! holds the number of rows in dest (= number of rows in rsrc)
   integer :: rptr ! used in determining the rows that belong to the
      ! source block rsrc
   integer :: s1sa, s1en ! point to the first and last rows of
      ! the block csrc within source block col.
   integer :: s2sa, s2en ! point to the first and last rows of
      ! the block rsrc within source block col.
   integer :: size_dnode ! number of rows in dnode
   integer :: size_scol ! size(sindex) ie. no. of rows in source block col.
!%%%integer :: t_start, t_end, this_thread

!%%%if(control%unit_log.gt.0) call system_clock(t_start)

   !
   ! Make a list of incident csrc rows (ie. columns of dest)
   !

   col_list_sz = 0
   row_list_sz = 0

   size_dnode = size(dnode%index)
   size_scol = size(sindex)

   ! Find block column dcol of dnode that dest belongs to. The block
   ! cols are numbered locally within dnode as 1,2,3,...
   dcol = blocks(dest)%bcol - blocks(dnode%blk_sa)%bcol + 1

   ! Set dcsa and dcen to hold indices
   ! of start and end columns in dcol (global column indices)
   dcsa = dnode%sa + (dcol-1)*dnode%nb                
   dcen = min(dnode%sa + dcol*dnode%nb-1, dnode%en)

   ! Compute cptr so that it points to the first row in csrc.
   ! loop while row index within source block column is less than the index
   ! of the first column in dest
   cptr = 1
   do 
      if(sindex(cptr).ge.dcsa) exit
      cptr = cptr + 1
      if(cptr.gt.size_scol) return ! No incident columns
   end do

   ! Set s1sa to point to first row in csrc
   s1sa = cptr 

   ! Now record the rows in csrc. Local row indices are used and
   ! are held in col_list(1:slen-slsa+1)
   do 
      if(sindex(cptr).gt.dcen) exit
      col_list_sz = col_list_sz + 1
      col_list(col_list_sz) = sindex(cptr) - dcsa + 1
      cptr = cptr + 1
      if(cptr.gt.size_scol) exit ! No more rows
   end do

   ! Set slen to point to last row in csrc
   s1en = cptr - 1 
   
   ! Loop over rsrc rows, building row list. Identify required data, form
   ! outer product of it into buffer.

   ! Find first and last rows of destination block
   i = dcol + blocks(dest)%id - blocks(dest)%dblk ! block in snode
   drsa = dnode%index(1 + (i-1)*dnode%nb)
   dren = dnode%index(min(1 + i*dnode%nb - 1, size_dnode))

   ! Find first row in rsrc
   rptr = s1sa
   do 
      if(sindex(rptr).ge.drsa) exit
      rptr = rptr + 1
      if(rptr.gt.size_scol) return ! No incident row! Shouldn't happen.
   end do
   s2sa = rptr ! Points to first row in rsrc

   ! Find the first and last rows of the destination block
   i = dest - blocks(dest)%dblk + 1 ! row block of dest column
   dptr_sa = 1 + (dcol-1 + i-1)*dnode%nb
   !dptr_en = min(dptr_sa + dnode%nb-1, size(dnode%index))

   ! Now record the rows in rsrc. Local row numbers
   ! are held in row_list(1:s2en-s2sa+1)
   dptr = dptr_sa ! Pointer for destination block

   do rptr = s2sa, size_scol
      if(sindex(rptr).gt.dren) exit
      do 
         if(dnode%index(dptr).ge.sindex(rptr)) exit
         dptr = dptr + 1
      end do
      row_list_sz = row_list_sz + 1
      row_list(row_list_sz) = dptr - dptr_sa + 1
   end do
   s2en = rptr - 1 ! Points to last row in rsrc

   ! Store starts of blocks rsrc and csrc for BLAS calls (remember blocks
   ! are stored by rows)
   i = 1 + nelim*(s2sa - 1) ! Start of rsrc
   j = 1 + nelim*(s1sa - 1) ! Start of csrc

   if(size(ldwork).lt.(s1en-s1sa+1)*nelim) then
      deallocate(ldwork,stat=st)
      allocate(ldwork((s1en-s1sa+1)*nelim),stat=st)
      if (st.ne.0) then
         info = MA86_ERROR_ALLOCATION
         call ma86_print_flag(info, control, context='MA86_factor',st=st)
         return
      end if
   endif

   call calc_ld(transa, s1en-s1sa+1, nelim, nelim, lcol_src(j), d_src, ldwork)

   if(n.ge.control%min_width_blas) then
      ! High flop/buffer sz ratio => perform operations into buffer with BLAS

      ndiag = 0
      call cgemm(transa, 'N', int(s1en-s1sa+1), int(s2en-s2sa+1), nelim,    &
         -cone, ldwork, nelim, lcol_src(i), nelim, czero, buffer, col_list_sz)

      !
      ! Apply update
      !

      ! Acquire lock on destination block
!$    call omp_set_lock(blocks(dest)%alock)

      if(.not.blocks(dest)%touched) then
         ! Zero block
         lcol( blocks(dest)%sa : &
            blocks(dest)%sa+blocks(dest)%blkn*blocks(dest)%blkm-1 ) = czero
         blocks(dest)%touched = .true.
      endif

      ndiag = 0
      if(blocks(dest)%id .eq. blocks(dest)%dblk) ndiag = int(s1en-s1sa+1)
      call expand_buffer(lcol(blocks(dest)%sa), blocks(dest)%blkn, row_list, &
         row_list_sz, col_list, col_list_sz, ndiag, buffer)

      ! Release lock
!$    call omp_unset_lock(blocks(dest)%alock)

   else
      ! Low flop/buffer ratio => perform update operation directly
      ! set ndiag if dest is a diagonal block
      ndiag = 0
      if(blocks(dest)%id .eq. blocks(dest)%dblk) ndiag = int(s1en-s1sa+1)

!$    call omp_set_lock(blocks(dest)%alock)

      if(.not.blocks(dest)%touched) then
         ! Zero block
         lcol( blocks(dest)%sa : &
            blocks(dest)%sa+blocks(dest)%blkn*blocks(dest)%blkm-1 ) = czero
         blocks(dest)%touched = .true.
      endif

      call update_direct(transa, lcol, lcol_src, ldwork, blocks(dest)%sa,    &
         blocks(dest)%blkn, 1, i, nelim, nelim, row_list, row_list_sz,       &
         col_list, col_list_sz, ndiag)

!$    call omp_unset_lock(blocks(dest)%alock)

   endif

!%%%if(control%unit_log.gt.0) then
!%%%   call system_clock(t_end)
!%%%   this_thread = 0
!%%%!$    this_thread = omp_get_thread_num()
!%%%   call log_task(control, this_thread, t_start, t_end, "UB", dest, &
!%%%      int(snode%sa,long), int(scol,long))
!%%%endif

end subroutine update_between

!*************************************************

!
! Performs an update_between task directly rather than via a buffer.
!
subroutine update_direct(transa, lcol, lcol_src, lcold_src, &
      dsa, blkn, s1sa, s2sa, &
      srcn, nelim, row_list, rls, col_list, cls, ndiag)
   character(len=1) :: transa ! set to 'C' for Hermitian and 'T' for symmetric
   complex(wp), dimension(*), intent(inout) :: lcol ! holds block
      ! column of factor L that destination block belongs to
   complex(wp), dimension(*), intent(in) :: lcol_src ! holds block
      ! column of factor L that the source blocks belong to
   complex(wp), dimension(*), intent(in) :: lcold_src ! holds block
      ! column of factor LD that the source blocks belong to
   integer, intent(in) :: dsa ! posn of the first entry of the
      ! destination block within lcol
   integer, intent(in) :: blkn ! number of columns in lcol
   integer, intent(in) :: s1sa ! points to first row in source block csrc
   integer, intent(in) :: s2sa ! points to first row in source block rsrc
   integer, intent(in) :: srcn ! number of columns in lcol_src
   integer, intent(in) :: nelim ! number of eliminations performed
      ! in the source block column (nelim = srcn if no delayed pivots).
   integer, intent(in) :: rls  ! size of row_list
   integer, intent(in) :: row_list(rls) ! local row indices for
      ! rows in rsrc (= local row indices for dest)
   integer, intent(in) :: cls  ! size of col_list
   integer, intent(in) :: col_list(cls) ! local row indices for
      ! rows in csrc (= local column indices for dest)
   integer, intent(in) :: ndiag ! Number of triangular rows of update

   integer :: cptr
   integer :: i
   integer :: j
   integer :: k1, k2, k3, k4, l
   integer :: rptr
   complex(wp) :: work1, work2, work3, work4

   ! Note: ndiag is either equal to 0 or, if dest is a diagonal block, it
   ! is equal to the number of cols cls in dest block (and this is
   ! equal to the number of rows rls in dest).
   ! first treat the case when dest is on the diagonal.
   ! loop over the rows of dest
   if (transa == 'C') then
      do j = 1, ndiag
         cptr = dsa + (row_list(j)-1) * blkn - 1
         rptr = s2sa + (j-1) * srcn
         do i = 1, j
            k1 = s1sa + (i-1) * nelim
            work1 = czero
            do l = rptr, rptr + nelim - 1
               work1 = work1 + conjg(lcold_src(k1))*lcol_src(l)
               k1 = k1 + 1
            end do
            k1 = cptr + col_list(i)
            lcol(k1) = lcol(k1) - work1
         end do
      end do

      ! Now consider the case when dest is not a diagonal block
      do j = ndiag+1, rls
         cptr = dsa + (row_list(j)-1) * blkn - 1
         rptr = s2sa + (j-1) * srcn
         i = 4*int(cls/4)
         do i = 1, i, 4
            k1 = s1sa + (i+0-1) * nelim
            k2 = s1sa + (i+1-1) * nelim
            k3 = s1sa + (i+2-1) * nelim
            k4 = s1sa + (i+3-1) * nelim
            work1 = czero; work2 = czero; work3 = czero; work4 = czero
            do l = rptr, rptr + nelim - 1
               work1 = work1 + conjg(lcold_src(k1))*lcol_src(l); k1 = k1 + 1
               work2 = work2 + conjg(lcold_src(k2))*lcol_src(l); k2 = k2 + 1
               work3 = work3 + conjg(lcold_src(k3))*lcol_src(l); k3 = k3 + 1
               work4 = work4 + conjg(lcold_src(k4))*lcol_src(l); k4 = k4 + 1
            end do
            k1 = cptr + col_list(i+0); lcol(k1) = lcol(k1) - work1
            k2 = cptr + col_list(i+1); lcol(k2) = lcol(k2) - work2
            k3 = cptr + col_list(i+2); lcol(k3) = lcol(k3) - work3
            k4 = cptr + col_list(i+3); lcol(k4) = lcol(k4) - work4
         end do
         i = 4*int(cls/4) + 1
         do i = i, cls
            k1 = s1sa + (i-1) * nelim
            work1 = 0
            do l = rptr, rptr + nelim - 1
               work1 = work1 + conjg(lcold_src(k1))*lcol_src(l)
               k1 = k1 + 1
            end do
            k1 = cptr + col_list(i)
            lcol(k1) = lcol(k1) - work1
         end do
      end do 

   else
      ! symmetric case
      do j = 1, ndiag
         cptr = dsa + (row_list(j)-1) * blkn - 1
         rptr = s2sa + (j-1) * srcn
         do i = 1, j
            k1 = s1sa + (i-1) * nelim
            work1 = czero
            do l = rptr, rptr + nelim - 1
               work1 = work1 + lcold_src(k1)*lcol_src(l)
               k1 = k1 + 1
            end do
            k1 = cptr + col_list(i)
            lcol(k1) = lcol(k1) - work1
         end do
      end do

      ! Now consider the case when dest is not a diagonal block
      do j = ndiag+1, rls
         cptr = dsa + (row_list(j)-1) * blkn - 1
         rptr = s2sa + (j-1) * srcn
         i = 4*int(cls/4)
         do i = 1, i, 4
            k1 = s1sa + (i+0-1) * nelim
            k2 = s1sa + (i+1-1) * nelim
            k3 = s1sa + (i+2-1) * nelim
            k4 = s1sa + (i+3-1) * nelim
            work1 = 0; work2 = 0; work3 = 0; work4 = 0
            do l = rptr, rptr + nelim - 1
               work1 = work1 + lcold_src(k1)*lcol_src(l); k1 = k1 + 1
               work2 = work2 + lcold_src(k2)*lcol_src(l); k2 = k2 + 1
               work3 = work3 + lcold_src(k3)*lcol_src(l); k3 = k3 + 1
               work4 = work4 + lcold_src(k4)*lcol_src(l); k4 = k4 + 1
            end do
            k1 = cptr + col_list(i+0); lcol(k1) = lcol(k1) - work1
            k2 = cptr + col_list(i+1); lcol(k2) = lcol(k2) - work2
            k3 = cptr + col_list(i+2); lcol(k3) = lcol(k3) - work3
            k4 = cptr + col_list(i+3); lcol(k4) = lcol(k4) - work4
         end do
         i = 4*int(cls/4) + 1
         do i = i, cls
            k1 = s1sa + (i-1) * nelim
            work1 = 0
            do l = rptr, rptr + nelim - 1
               work1 = work1 + lcold_src(k1)*lcol_src(l)
               k1 = k1 + 1
            end do
            k1 = cptr + col_list(i)
            lcol(k1) = lcol(k1) - work1
         end do
      end do
   end if
end subroutine update_direct

!*************************************************

!
! Optimised sparse expansion by lists
! Entry buffer(i,j) gets added to a(row_list(i),col_list(j))
! Note: by uncommenting the i2 and mm_prefetch lines we can speed this up
!       for the Intel compiler by explicitly prefetching the next row.
!
subroutine expand_buffer(a, blkn, row_list, rls, col_list, cls, ndiag, buffer)
   complex(wp), dimension(*), intent(inout) :: a ! holds L
   integer, intent(in) :: blkn ! number of cols in destination block
   integer, intent(in) :: rls ! size of row_list
   integer, intent(in) :: row_list(rls)
   integer, intent(in) :: cls ! size of col_list
   integer, intent(in) :: col_list(cls)
   integer, intent(in) :: ndiag ! Number of triangular rows of update
   complex(wp), intent(in) :: buffer(rls*cls)

   integer :: i, j, k, rptr, cptr
   !integer :: i2

   rptr = 1
   do j = 1, ndiag
      cptr = 1 + (row_list(j)-1)*blkn - 1
      do i = 1, j
         k = cptr + col_list(i)
         a(k) = a(k) + buffer(rptr)
         rptr = rptr + 1
      end do
      rptr = rptr + (cls-j)
   end do
   do j = ndiag+1, rls-1
      cptr = (row_list(j)-1) * blkn
      !i2 = (row_list(j+1)-row_list(j)) * blkn
      do i = 1, cls
         k = cptr + col_list(i)
         a(k) = a(k) + buffer(rptr)
         !call mm_prefetch(a(k+i2),1)
         rptr = rptr + 1
      end do
   end do

   if(ndiag.lt.rls) then
      cptr = (row_list(rls)-1) * blkn
      do i = 1, cls
         k = cptr + col_list(i)
         a(k) = a(k) + buffer(rptr)
         rptr = rptr + 1
      end do
   endif

end subroutine expand_buffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Update columns jl:jr for pivots pl:pr
! This was based on a routine within ma64.
!
subroutine ma64_update(n,nb,kkq,a,buf,jl,jr,pl,pr)
   integer, intent(in) :: n
   integer, intent(in) :: nb
   integer(long), intent(in) :: kkq
   complex(wp), intent(inout) :: a(*)
   complex(wp), intent(in) :: buf(*)
   integer, intent(in) :: jl,jr,pl,pr

   ! Local variables
   integer i
   integer j ! column index
   integer j2
   integer p
   integer(long) kkj ! Position of diagonal j

   integer, parameter :: bs=64 ! block size for update

   j = max(jl,1)
   j2 = min(nb,jr)

   if(n.lt.j .or. j2.lt.j .or. pl.gt.pr) return

   do i = 1, j2-j+1, bs
      kkj = 1 + (j-1+i-1)*(n+1_long)
      p = min(bs, j2-j+1-i+1)
      call cgemm ('n', 't', n-j+1-i+1, p, pr-pl+1, cone, &
         a(kkq+j-1+n*(pl-1)+i-1), n, buf(j+n*(pl-1)+i-1), n, &
         cone, a(kkj), n)
   end do

end subroutine ma64_update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Swap columns. It may be assumed that 1<=j1<j2<=p.
! This was based on a routine within ma64.
!
subroutine swap_cols(matrix_type,n,nb,q,q1,rm,a,buf,perm,j1,j2)
   integer, intent(in) :: matrix_type
   integer, intent(in) :: n
   integer, intent(in) :: nb
   integer, intent(in) :: q
   integer, intent(in) :: q1
   integer, intent(in) :: rm
   complex(wp), intent(inout) :: a(*)
   complex(wp), intent(inout) :: buf(*)
   integer, intent(inout) :: perm(*)
   integer, intent(in) :: j1,j2

   integer d1,d2 ! Positions of the diagonal entries for columns j1,j2
   integer j ! Index of the leading column of the block
   integer jb ! Size of block
   integer(long) k1 ! Start of current part of column j1
   integer(long) k2 ! Start of current part of column j2
   integer l !
   complex(wp) temp

   !print *, "swapping ", j1, j2

   j = perm(j1)
   perm(j1) = perm(j2)
   perm(j2) = j

   ! Swap columns of buffer
   call cswap(q-q1+1-rm,buf(j1+rm*n),n,buf(j2+rm*n),n)
   !call cswap(q-q1+1,buf(j1),n,buf(j2),n)

   k1 = j1
   k2 = j2

   ! Swap rows
   ! i.e. a(j1, 1:j1-1) <-> a(j2, 1:j1-1)
   l = j1-1
   if (l>0) call cswap(l,a(k1),n,a(k2),n)

   ! calulate values for next phase
   d1 = k1 + l*int(n,long)
   k1 = d1 + 1
   k2 = k2 + l*int(n,long) + n

   ! Swap columns with rows
   ! i.e. a(j1+1:j2-1, j1) <-> a(j2, j1+1:j2-1)
   jb = j2-j1-1
   l = min(jb,nb-j1,j2-1)
   if (l>0.and.j2.ge.1) call cswap(l,a(k1),1,a(k2),n)
   if(matrix_type==HSL_MATRIX_CPLX_HERM_INDEF) then
      a(k1:k1+l-1) = conjg(a(k1:k1+l-1))
      a(k2:k2+(l-1)*n:n) = conjg(a(k2:k2+(l-1)*n:n))
      a(k1+l) = conjg(a(k1+l)) ! No idea what this is for?
   end if

   ! set up values for next phase
   d2 = k2 + l*int(n,long)
   k2 = d2 + 1
   k1 = k1 + l + 1

   ! Swap the diagonals
   ! i.e. a(j1,j1) <-> a(j2,j2)
   temp = a(d1)
   a(d1) = a(d2)
   a(d2) = temp

   ! Swap columns
   ! i.e. a(j2+1:n,j1) <-> a(j2+1:n,j2)
   if (n>j2) call cswap(n-j2,a(k1),1,a(k2),1)

end subroutine swap_cols

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Solves                ( D   ) x = b
!                       (   I )
!
! where D is block diagonal of size q with blocks of size 1 and 2.
! This is taken from the ma64 package.
!
subroutine solveD(n,q,b,d,matrix_type)
   integer, intent (in) :: n ! Matrix order
   integer, intent (in) :: q ! Order of L (and D)
   complex(wp), intent (inout) :: b(n) ! b(n) holds the right-hand
      ! sides on entry and is overwritten by the solution.
   complex(wp), intent (in) :: d(2*q) !  Holds the inverse of D,
      ! except that zero diagonal blocks (pivots) are not inverted.
      ! Diagonal entries are in d(1:2*q-1:2) and entries to the right of
      ! the diagonal are in d(2:2*q-2:2).
   integer :: matrix_type ! -4 for Hermitian or -5 for symmetric

   ! .. Locals ..
   integer :: j ! Column index
   complex(wp) :: temp

   ! Apply operations associated with D
   j = 1
   do
      if (j>q) exit
      if (d(2*j)==czero) then
         ! 1x1 pivot
         b(j) = b(j)*d(2*j-1)
         j = j + 1
      else
         ! 2x2 pivot
         if(matrix_type==HSL_MATRIX_CPLX_HERM_INDEF) then
            temp = b(j)*d(2*j-1) + b(j+1)*conjg(d(2*j))
         else
            temp = b(j)*d(2*j-1) + b(j+1)*d(2*j)
         end if
         b(j+1) = b(j)*d(2*j) + b(j+1)*d(2*j+1)
         b(j) = temp
         j = j + 2
      end if
   end do
end subroutine solveD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Adds D^-1 rhs(index(1:q)) to rhs_local(1:q)
!
subroutine solve_add_D(q,d,rhs,rhs_local,index,matrix_type)
   integer, intent (in) :: q ! Order of L (and D)
   complex(wp), intent (in) :: d(2*q) !  Holds the inverse of D,
      ! except that zero diagonal blocks (pivots) are not inverted.
      ! Diagonal entries are in d(1:2*q-1:2) and entries to the right of
      ! the diagonal are in d(2:2*q-2:2).
   complex(wp), intent(in) :: rhs(*)
   complex(wp), intent(inout) :: rhs_local(*)
   integer, intent(in) :: index(*)
   integer :: matrix_type ! -4 for Hermitian or -5 for symmetric

   ! .. Locals ..
   integer i ! Column index into rhs
   integer i2 ! Column index into rhs for (j+1)
   integer j ! Column index into rhs_local

   ! Apply operations associated with D
   j = 1
   do while(j.le.q)
      i = index(j)
      if (d(2*j)==czero) then
         ! 1x1 pivot
         rhs_local(j) = rhs_local(j) + rhs(i)*d(2*j-1)
         j = j + 1
      else
         ! 2x2 pivot
         i2 = index(j+1)
         if(matrix_type==HSL_MATRIX_CPLX_HERM_INDEF) then
            rhs_local(j)   = rhs_local(j)   + rhs(i)*d(2*j-1) + &
               rhs(i2)*conjg(d(2*j))
         else
            rhs_local(j)   = rhs_local(j)   + rhs(i)*d(2*j-1) + &
               rhs(i2)*d(2*j)
         end if
         rhs_local(j+1) = rhs_local(j+1) + rhs(i)*d(2*j) + rhs(i2)*d(2*j+1)
         j = j + 2
      end if
   end do
end subroutine solve_add_D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine updates the dependency count for blocks in the ancestor
! anode of snode.
!
subroutine calc_dep(cptr, snode, anode, nodes, blocks, swidth, map)
   integer, intent(inout) :: cptr 
     ! pointer into the row indices for snode. On entry points
     ! to the first row corresponding to the column set for anode.
     ! On exit, updated ready for the next ancestor.
     ! It could feasibly point beyond the column set for anode if
     ! there are no entries for anode's columns in the row set of snode.
   integer, intent(in) :: snode ! node
   integer, intent(in) :: anode ! ancestor of snode
   type(node_type), dimension(-1:), intent(in) :: nodes ! node info
   type(block_type), dimension(:), intent(inout) :: blocks ! block info. On 
     ! exit, dependency count updated
   integer, intent(in) :: swidth ! number of block cols in snode
   integer, dimension(:), intent(in) :: map ! For each row k in j-th 
     ! block column of anode, map(k) is set to j

   integer :: cb ! index of block column in anode
   integer(long) :: dblk ! id of diagonal block in anode
   integer :: jlast ! last column of block in anode
   integer :: nb ! set to nodes(anode)%nb (block size for anode)
   integer :: i, jb, k, k1
   integer :: size_snode ! number of rows in snode

   nb = nodes(anode)%nb
   size_snode = size(nodes(snode)%index)
   do 
      if(cptr.gt.size_snode) exit
      if (nodes(snode)%index(cptr).gt.nodes(anode)%en) exit

      ! Find block column of anode
      cb = 1 + (nodes(snode)%index(cptr) - nodes(anode)%sa) / nb

      ! Find diagonal block in block column cb

      ! loop over block columns. blocks(dblk)%last_blk is the last
      ! block in the block column to which dblk belongs and
      ! so blocks(dblk)%last_blk + 1 is first block in the next
      ! block column, which is the diagonal block in that block column

      dblk = nodes(anode)%blk_sa ! first block in anode
      do i = 1, cb-1
         dblk = blocks(dblk)%last_blk + 1
      end do

      ! Increment dep count for each block involved.
      ! loop over rows in snode
      jb = -1 ! Last block
      do i = cptr, size_snode
         k1 = nodes(snode)%index(i)
         k = map(k1)
         ! k is local block number. When we reach a new block,
         ! we increment dep by the number of block columns in snode
         ! (since each block col. in snode will have to update this block)
         if(k.ne.jb) blocks(dblk+k-cb)%dep_initial = &
            blocks(dblk+k-cb)%dep_initial + swidth
         jb = k
      end do

      ! Move cptr to first row in another block of anode
      jlast = min(nodes(anode)%sa + cb*nb - 1, nodes(anode)%en)
      do cptr = cptr, size_snode
         if(nodes(snode)%index(cptr) > jlast) exit
      end do

   end do

end subroutine calc_dep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Initialize the task pool, which is held using a variable of derived
! type taskstack
!
subroutine init_stack(stack, pool_size, control, info, st)
   type(taskstack), intent(out) :: stack ! see description of derived type
   integer, intent(in) :: pool_size
   type(MA86_control), intent(in) :: control  ! see description of derived type
   integer, intent(out) :: info ! error flag. Only possible error is
      ! an allocation error
   integer, intent(out) :: st ! stat parameter

   ! local variables
   integer :: i
   integer :: ncache
   integer :: total_threads


   ! Initialise
   info = 0; st = 0
   total_threads = 1
!$ total_threads = omp_get_max_threads()

   stack%pool_size = pool_size
   stack%total = 0
   stack%active = 0
   stack%freehead = 1
   stack%abort = .false.

   ncache = calc_cache(total_threads-1, control)

   deallocate(stack%ctasks,stat=st)
   deallocate(stack%cheads,stat=st)
!$ deallocate(stack%clocks,stat=st)
   deallocate(stack%tasks,stat=st)
!**deallocate(stack%waiting,stat=st)
   deallocate(stack%next,stat=st)

   allocate(stack%ctasks(control%cache_tq_sz, ncache),  &
            stack%cheads(ncache),   &
!$          stack%clocks(ncache),   &
            stack%tasks(pool_size), &
            stack%next(pool_size),  &
!**         stack%waiting(0:total_threads-1), &
            stat=st)

   if (st.ne.0) then 
     info = MA86_ERROR_ALLOCATION
     return
   end if

!$ call omp_init_lock(stack%lock)

   ! Initialise stack
   do i = 1, stack%pool_size-1
      stack%next(i) = i + 1
   end do
   stack%next(stack%pool_size) = -1
!**stack%waiting(:) = 0.0
   stack%cheads(:) = 0
   stack%prihead(1:4) = -1 ! empty

!$ do i = 1, ncache
!$    call omp_init_lock(stack%clocks(i))
!$ end do

end subroutine init_stack

!*************************************************

!
! Free any resources involved in task pool
!
subroutine cleanup_stack(stack)
   type(taskstack), intent(inout) :: stack ! see description of derived type

   ! local variables
!$ integer :: i ! temporary variable
   integer :: st ! stat parameter

!$ if(allocated(stack%tasks)) call omp_destroy_lock(stack%lock)

   deallocate(stack%tasks,stat=st)
   deallocate(stack%ctasks,stat=st)
   deallocate(stack%cheads,stat=st)
!**deallocate(stack%waiting,stat=st)
   deallocate(stack%next,stat=st)

!$  if(allocated(stack%clocks)) then
!$    do i = 1, size(stack%clocks)
!$       call omp_destroy_lock(stack%clocks(i))
!$    end do
!$    deallocate(stack%clocks,stat=st)
!$  endif

end subroutine cleanup_stack

!*************************************************

!
! This subroutine ensures all components of task are defined
!
subroutine zero_task(task)
   type(dagtask), intent(out) :: task

   task%task_type = 0
   task%dest = 0
   task%src1 = 0
   task%src2 = 0
   task%csrc(:) = 0
   task%rsrc(:) = 0
end subroutine zero_task

!*************************************************

!
! Add a task to the local stack or the task pool.
! In fact we add it to the cache's local stack, if this is full then it
! gets thrown to the task pool. 
!
subroutine add_task(stack, task, control, info, st)
   type(taskstack), intent(inout) :: stack ! holds local stacks and task pool
   type(dagtask), intent(in) :: task ! task to be added
   type(MA86_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(inout) :: st

   ! Possible error returns are info = MA86_ERROR_ALLOCATION
   ! Possible warnings are      info = MA86_WARNING_POOL_SMALL,
   !                                   MA86_WARNING_POOL_SING
   
   integer :: cache, i, j, this_thread

   this_thread = 0
!$ this_thread = omp_get_thread_num()

   cache = calc_cache(this_thread, control)

!$ call omp_set_lock(stack%clocks(cache))

   if(stack%cheads(cache).eq.control%cache_tq_sz) then
      ! Local stack is full.
      ! Clear lower half of stack to task pool.
!$    call omp_set_lock(stack%lock)
      do i = 1, control%cache_tq_sz / 2
         call add_task_g(stack, stack%ctasks(i,cache), &
            control, info, st, locked=.true.)
         if(info.lt.0) then
!$          call omp_unset_lock(stack%lock)
!$          call omp_unset_lock(stack%clocks(cache))
            return
         endif
      end do
!$    call omp_unset_lock(stack%lock)
      j = 1
      do i = control%cache_tq_sz / 2 + 1, control%cache_tq_sz
         stack%ctasks(j,cache) = stack%ctasks(i,cache)
         j = j + 1
      end do
      stack%cheads(cache) = stack%cheads(cache) - control%cache_tq_sz / 2
   endif

   ! Add to top of local stack
   stack%cheads(cache) = stack%cheads(cache) + 1
   stack%ctasks(stack%cheads(cache), cache) = task

!$ call omp_unset_lock(stack%clocks(cache))

end subroutine add_task

!*************************************************

!
! Add a task to the task pool at a given priority.
! The tasks have different priorities and those with the same
! priority are held within the pool using a linked list
!
subroutine add_task_g(stack, task, control, info, st, locked)
   type(taskstack), intent(inout) :: stack ! holds local stacks and task pool
   type(dagtask), intent(in) :: task ! task to be added to task pool.
   type(MA86_control), intent(in) :: control
   integer, intent(inout) :: info ! error flag (only possible error is
      ! an allocation error)
   integer, intent(inout) :: st ! stat parameter
   logical, optional, intent(in) :: locked ! indicates whether stack%lock
      ! is set. If not present, stack%lock will be set and then unset before
      ! return

   ! Possible error return is info = MA86_ERROR_ALLOCATION
   ! Possible warning is info = MA86_WARNING_POOL_SMALL,
   !                            MA86_WARNING_POOL_SING

   integer :: i
   integer :: insert ! holds free location in task pool
   integer :: priority ! priority of task
   type(dagtask), dimension(:), allocatable :: temp_tasks ! only allocated if
      ! task pool is full and has to be reallocated, in which case
      ! temp_task is used to hold a temporary copy.
   integer, dimension(:), allocatable :: temp_next ! used only if
      ! task pool has to be increased in size.

   ! set priority according to type of task
   select case(task%task_type)
   case(TASK_FACTORIZE_COLUMN)
      priority = 1
   case(TASK_UPDATE_INTERNAL)
      priority = 3
   case(TASK_UPDATE_BETWEEN)
      priority = 4
   case(TASK_SLV_FSLV, TASK_SLV_BSLV)
      priority = 1
   !case default
   !   priority = 5
   end select

!$ if(.not.present(locked)) call omp_set_lock(stack%lock)

   !
   ! Find a free spot
   !
   insert = stack%freehead
   if(insert.eq.-1) then
      ! We have overflowed the task pool, we need to reallocate it

      ! Copy tasks
      allocate(temp_tasks(stack%pool_size),stat=st)
      if(st.ne.0) go to 10
      temp_tasks(:) = stack%tasks(:)

      deallocate(stack%tasks,stat=st)
      allocate(stack%tasks(stack%pool_size*2),stat=st)
      if(st.ne.0) go to 10

      stack%tasks(1:stack%pool_size) = temp_tasks(:)
      deallocate(temp_tasks,stat=st)

      ! Copy next
      allocate(temp_next(stack%pool_size),stat=st)
      if(st.ne.0) go to 10
      temp_next(:) = stack%next(:)

      deallocate(stack%next,stat=st)
      allocate(stack%next(stack%pool_size*2),stat=st)
      if(st.ne.0) go to 10

      stack%next(1:stack%pool_size) = temp_next(:)
      deallocate(temp_next,stat=st)

      ! Extend linked list
      stack%freehead = stack%pool_size + 1
      do i = stack%pool_size+1, stack%pool_size*2-1
         stack%next(i) = i + 1
      end do
      stack%next(stack%pool_size*2) = -1

      ! Increase stored pool_size
      stack%pool_size = stack%pool_size*2
      insert = stack%freehead

      select case(info)
      case(MA86_WARNING_POOL_SING)
         ! do nothing, warning already flagged
      case(MA86_WARNING_SINGULAR)
         ! set double warning
         info = MA86_WARNING_POOL_SING
      case default
         ! set normal warning
         info = MA86_WARNING_POOL_SMALL
      end select
      if(control%diagnostics_level.ge.0 .and. control%unit_warning.ge.0) &
         write(control%unit_warning, "(a,i8)")&
            " Task pool size increased to = ", stack%pool_size
   end if

   stack%freehead = stack%next(insert)

   !
   ! Place task in pool in position insert and add into the linked list of 
   ! tasks with same priority
   !
   stack%tasks(insert) = task
   stack%next(insert) = stack%prihead(priority)
   stack%prihead(priority) = insert
   stack%total = stack%total + 1

   !
   ! Update highest priority (the task with the highest priority
   ! is the one with the lowest priority value)
   !
   stack%lowest_priority_value = min(stack%lowest_priority_value, priority)

!$ if(.not.present(locked)) call omp_unset_lock(stack%lock)
   return ! Successful exit

   ! Error handling in case of allocation failure
   10 continue
   info = MA86_ERROR_ALLOCATION
   call MA86_print_flag(info, control, context='MA86_factor',st=st)
   stack%abort = .true.
!$ if(.not.present(locked)) call omp_unset_lock(stack%lock)
   return

end subroutine add_task_g

!*************************************************

!
! Figure out which cache this thread belongs to - allows easy changes between
! machines by requiring only one function be altered.
!
integer function calc_cache(thread, control)
   integer, intent(in) :: thread
   type(MA86_control), intent(in) :: control

   integer :: total_threads

   select case (control%cache_layout)
   case(CACHE_COMPACT)
      calc_cache = thread / control%cache_cores + 1
   case(CACHE_SCATTER)
      total_threads = 1
!$    total_threads = omp_get_max_threads()
      total_threads = max(1, total_threads/control%cache_cores)
         ! (Avoid div by zero)
      calc_cache = mod(thread, total_threads) + 1
   case default ! CACHE_IDENITY
      calc_cache = thread + 1
   end select

end function calc_cache

!*************************************************

!
! Get a task; if none remain then end.
! If we can't find any work in the local task stack then we first try the
! task pool, and if this doesn't contain any, we steal tasks from
! other caches.
!
subroutine get_task(stack, task, control, info, st)
   type(taskstack), intent(inout) :: stack ! holds local stacks and task pool
   type(dagtask), intent(inout) :: task
   type(MA86_control), intent(in) :: control
   integer, intent(inout) :: info ! error flag
   integer, intent(inout) :: st ! stat parameter

   integer :: this_thread, cache

   this_thread = 0
!$ this_thread = omp_get_thread_num()
   cache = calc_cache(this_thread, control)

!$OMP FLUSH(stack)
   !
   ! Check for error termination (no lock required)
   !
   if(stack%abort) then
      ! Error abort
      task%task_type = TASK_DONE
      if(task%task_type.ne.TASK_NONE) then
!$       call omp_set_lock(stack%lock)
         ! decrease count of acive threads
         stack%active = stack%active - 1
!$       call omp_unset_lock(stack%lock)
      endif
!$OMP FLUSH
      return
   endif

!$ call omp_set_lock(stack%clocks(cache))
   if(stack%cheads(cache).ne.0) then
      ! local task stack contains tasks ... take the one off the top and return
      task = stack%ctasks(stack%cheads(cache),cache)
      stack%cheads(cache) = stack%cheads(cache) - 1
!$    call omp_unset_lock(stack%clocks(cache))
      return
   endif
!$ call omp_unset_lock(stack%clocks(cache))

!$ call omp_set_lock(stack%lock)

   ! reduce count of number of active threads
   stack%active = stack%active - 1

   ! If no task in local task stack then we /must/ get one.
   ! First attempt to take a task from the task pool.
   ! If this pool is empty, search for largest local stack
   ! belonging to another cache. If found, move tasks in bottom half 
   ! from this local stack to the task pool (workstealing).
   ! Then take task of highest priority from the pool as next task.
 
   task%task_type = TASK_NONE
   do while(task%task_type.eq.TASK_NONE)
      call get_task_from_pool(stack, task)

      if(info.lt.0) then
!$       call omp_unset_lock(stack%lock)
         return
      endif
      if(task%task_type.eq.TASK_NONE) then
         ! Check if we've finished
         if(stack%active.eq.0) then
!$          call omp_unset_lock(stack%lock)
            task%task_type = TASK_DONE
            return
         endif

!$       call omp_unset_lock(stack%lock)

         ! Spin until a task become available
         call await_task(stack,control)
         ! immediate return if we have to abort
         if(stack%abort) return

         if(stack%cheads(cache).ne.0) then
            ! tasks available in local task stack. Take one from the top.
!$          call omp_set_lock(stack%clocks(cache))
            if(stack%cheads(cache).ne.0) then
               task = stack%ctasks(stack%cheads(cache),cache)
               stack%cheads(cache) = stack%cheads(cache) - 1
!$             call omp_unset_lock(stack%clocks(cache))
!$             call omp_set_lock(stack%lock)
               exit
            endif
!$          call omp_unset_lock(stack%clocks(cache))
         endif

         if(stack%total.le.0) then
            ! nothing left in task pool so look to steal tasks
            ! from another thread. 
            call worksteal(stack, control, info, st)
            if(info.lt.0) return
         endif
!$       call omp_set_lock(stack%lock)
         cycle
      endif
   enddo

   if(task%task_type.ne.TASK_DONE) stack%active = stack%active + 1

!$ call omp_unset_lock(stack%lock)

end subroutine get_task

!*************************************************

!
! Look to other caches to steal work from.
! Steals from the largest local stack. Moves bottom half
! of this stack into the task pool and then moves the
! remaining tasks down the stack
!
subroutine worksteal(stack, control, info, st)
   type(taskstack), intent(inout) :: stack ! holds local stacks and task pool
   type(MA86_control), intent(in) :: control
   integer, intent(inout) :: info ! error flag
   integer, intent(inout) :: st ! stat parameter

   integer :: i, j
   integer :: mi ! index oflargest local stack (-1 if all are empty)
   integer :: msz ! size of largest locak stack
   integer :: cache, total_threads, this_thread

   total_threads = 1
!$ total_threads = omp_get_num_threads()
   this_thread = 0
!$ this_thread = omp_get_thread_num()
   cache = calc_cache(this_thread-1, control)

   ! Find biggest stack to steal from
   msz = -1
   mi = -1
   do i = 1, calc_cache(total_threads-1, control)
      if(cache.eq.i) cycle
      if(stack%cheads(i).gt.msz) then
!$       call omp_set_lock(stack%clocks(i))
         if(stack%cheads(i).gt.msz) then
!$          if(mi.ne.-1) call omp_unset_lock(stack%clocks(mi))
            mi = i
            msz = stack%cheads(i)
         else ! Its no longer bigger, release lock
!$          call omp_unset_lock(stack%clocks(i))
         endif
      endif
   end do
   if(mi.eq.-1) return ! No other caches

   msz = stack%cheads(mi)

   ! Steal half from bottom of the stack mi
!$ call omp_set_lock(stack%lock)
   do i = 1, msz / 2
      call add_task_g(stack, stack%ctasks(i,mi), &
         control, info, st, locked=.true.)
      if(info.lt.0) then
!$       call omp_unset_lock(stack%lock)
         return
      endif
   end do
!$ call omp_unset_lock(stack%lock)
   ! move the remaining tasks down to occupied freed up space
   j = 1
   do i = msz / 2 + 1, msz
      stack%ctasks(j,mi) = stack%ctasks(i,mi)
      j = j + 1
   end do
   stack%cheads(mi) = stack%cheads(mi) - msz / 2
!$ call omp_unset_lock(stack%clocks(mi))
end subroutine worksteal

!*************************************************
!
! Get a task from the task pool. Want a task with as low a priority
! value as possible.
!
subroutine get_task_from_pool(stack, task)
   type(taskstack), intent(inout) :: stack ! holds local stacks and task pool
   type(dagtask), intent(inout) :: task

   integer :: priority ! priority of task
   integer :: ptr ! pointer into stack of tasks
   !
   ! Check there exists a task for us to get
   !
   if(stack%total.le.0) then

      task%task_type = TASK_NONE
      return
   else
      !
      ! Find a task with the lowest priority value (that is, the
      ! task which has the highest priority)
      ! start with the current lowest priority value and increase priority
      ! until a non empty linked list of tasks is found
      !
      priority = stack%lowest_priority_value
      do while(stack%prihead(priority).eq.-1)
         priority = priority + 1
      end do
      stack%lowest_priority_value = priority

      stack%max_pool_size = max(stack%max_pool_size,stack%total)

      !
      ! Grab a task from the stack of tasks
      !
      ptr = stack%prihead(priority)
      task = stack%tasks(ptr)
      stack%prihead(priority) = stack%next(ptr)
      stack%next(ptr) = stack%freehead
      stack%freehead = ptr

      ! reduce count of tasks in stack
      stack%total = stack%total - 1

   endif
end subroutine get_task_from_pool

!*************************************************

!
! Spin until either some work crops up or notification of abort
! received or all work has been performed so nothing to wait for
!
subroutine await_task(stack, control)
   type(taskstack), intent(inout) :: stack ! holds local stacks and task pool
   type(MA86_control), intent(in) :: control

   integer :: this_thread, cache
!**integer :: t_start, t_end, t_rate

   this_thread = 0
!$ this_thread = omp_get_thread_num()
   cache = calc_cache(this_thread, control)

 !**  call system_clock(t_start, t_rate)

!$OMP CRITICAL (await_task_idle)
   do
!$OMP FLUSH(stack)
      if(stack%abort .or. stack%total.gt.0 .or. &
         (stack%total.eq.0 .and. stack%active.eq.0) .or. &
          any(stack%cheads(:) .gt. 4) .or. stack%cheads(cache).ne.0) exit
   end do
!$OMP END CRITICAL (await_task_idle)

!**if (control%time_out.ge.0) then
!**  call system_clock(t_end)
!**  stack%waiting(this_thread) = &
!**     stack%waiting(this_thread) + (t_end-t_start)/real(t_rate)
!**end if

end subroutine await_task

!*************************************************

!
! Sets the abort flag, called on error
!
subroutine set_abort(stack)
   type(taskstack), intent(inout) :: stack ! holds local stacks and task pool

   ! Set the state
!$ call omp_set_lock(stack%lock)
   stack%abort = .true.
!$ call omp_unset_lock(stack%lock)
end subroutine set_abort

!*************************************************

!
! TASK_FSOLV
! B_j <- L_jj^-1 B_j
!
! Hermitian case:
! conjg(L_jj) is stored by rows so we have to take 
! conjugate transpose in blas calls
! Symmetric case: take transpose in blas calls
!
! Note: While diagonal blocks may be trapezoidal, this is handled at the
! level calling this subroutine through a call to slv_fwd_update or
! slv_bwd_update
!
subroutine slv_solve_fwd(transa, n, nelim, col, dest, nrhs, rhs, ldr, &
      control, id)
   character(len=1) :: transa ! set to 'C' for Hermitian and 'T' for symmetric
   integer, intent(in) :: n ! leading dimension of L_jj
   integer, intent(in) :: nelim ! number eliminations (immediate return if =0)
   integer, intent(in) :: col ! start of block column variables in rhs
   complex(wp), dimension(n*nelim), intent(in) :: dest ! holds L_jj
   integer, intent(in) :: nrhs ! number of right-hand sides
   integer, intent(in) :: ldr ! leading extent of rhs
   complex(wp), intent(inout) :: rhs(ldr*nrhs)
   type(MA86_control), intent(in) :: control
   integer(long), intent(in) :: id

!%%%integer :: this_thread, t_start, t_end

   if(nelim.eq.0) return

!%%%if(control%unit_log.gt.0) call system_clock(t_start)

   if(nrhs.eq.1) then
      call ctrsv('Upper', transa, 'Unit', nelim, dest, n, rhs(col), 1)
   else
      call ctrsm('Left', 'Upper', transa, 'Unit', nelim, nrhs, cone, &
         dest, n, rhs(col), ldr)
   endif

!%%%if(control%unit_log.gt.0) then
!%%%   call system_clock(t_end)
!%%%   this_thread = 0
!%%%!$    this_thread = omp_get_thread_num()
!%%%   call log_task(control, this_thread, t_start, t_end, "FS", id)
!%%%endif

end subroutine slv_solve_fwd

!*************************************************

!
! TASK_BSOLV
! B_j <- L_jj^-* B_j
!
! Hermitian case:
! conjg(L_jj) is stored by rows (that is, we hold L_jj^H by cols so no
! transpose needed in blas calls)
! Symmetric case:
! L_jj is stored by rows (that is, we hold L_jj^T by cols so no
! transpose needed in blas calls)
!
! Note: While diagonal blocks may be trapezoidal, this is handled at the
! level calling this subroutine through a call to slv_fwd_update or
! slv_bwd_update
!
subroutine slv_solve_bwd(n, nelim, col, dest, nrhs, rhs, ldr, control, id)
   integer, intent(in) :: n ! leading dimension of L_jj
   integer, intent(in) :: nelim ! number eliminations (immediate return if =0)
   integer, intent(in) :: col ! start of block column variables in rhs
   complex(wp), dimension(n*nelim), intent(in) :: dest ! holds L_jj
   integer, intent(in) :: nrhs ! number of right-hand sides
   integer, intent(in) :: ldr ! leading extent of rhs
   complex(wp), intent(inout) :: rhs(ldr*nrhs)
   type(MA86_control), intent(in) :: control
   integer(long), intent(in) :: id

!%%%integer :: this_thread, t_start, t_end

   if(nelim.eq.0) return

!%%%if(control%unit_log.gt.0) call system_clock(t_start)

   if(nrhs.eq.1) then
      call ctrsv('Upper', 'N', 'Unit', nelim, dest, n, rhs(col), 1)
   else
      call ctrsm('Left', 'Upper', 'N', 'Unit', nelim, nrhs, cone, &
         dest, n, rhs(col), ldr)
   endif

!%%%if(control%unit_log.gt.0) then
!%%%   call system_clock(t_end)
!%%%   this_thread = 0
!%%%!$    this_thread = omp_get_thread_num()
!%%%   call log_task(control, this_thread, t_start, t_end, "BS", id)
!%%%endif

end subroutine slv_solve_bwd

!*************************************************

!
! TASK_FUPD
! B_j <- B_j - L_ij B_i
!
! Hermitian case:
! conjg(L_ij) is stored by rows (ie L_ij^H is stored by cols)
! so we have to take conjugate transpose in blas calls
! Symmetric case:
! L_ij is stored by rows (ie L_ij^T is stored by cols)
! so we have to take transpose in blas calls
!
subroutine slv_fwd_update(transa, m, nelim, col, offset, index, dest, ldd,  &
      nrhs, upd, ldu, rhs, ldr, xlocal, control, id)
   character(len=1) :: transa ! set to 'C' for Hermitian and 'T' for symmetric
   integer, intent(in) :: m ! number of rows in block
   integer, intent(in) :: nelim ! number eliminations (immediate return if =0)
   integer, intent(in) :: col ! start of block column variables in rhs
   integer, intent(in) :: offset ! offset into index we start at
   integer, dimension(*), intent(in) :: index
   integer, intent(in) :: ldd ! leading dimension of dest
   complex(wp), dimension(ldd*m), intent(in) :: dest ! holds L_ij
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldu  ! leading extent of upd
   complex(wp), intent(inout) :: upd(ldu*nrhs) ! vector to update
   integer, intent(in) :: ldr  ! leading extent of rhs
   complex(wp), intent(in) :: rhs(ldr*nrhs) ! rhs vector
   complex(wp), dimension(*), intent(out) :: xlocal
   type(MA86_control), intent(in) :: control
   integer(long), intent(in) :: id

   integer :: i
   integer :: j
   integer :: k
   integer :: r ! right hand side loop variable
   complex(wp) :: w ! temporary work value
!%%%integer :: t_start, t_end, this_thread

   if(nelim.eq.0) return

!%%%if(control%unit_log.gt.0) call system_clock(t_start)

   ! forward substitution
   if(nrhs.eq.1) then
      if(m-nelim.gt.10 .and. nelim.gt.4) then
         !!! Single rhs, BLAS 2

         call cgemv(transa, nelim, m, -cone, dest, ldd, rhs(col), 1, czero,  &
            xlocal, 1)
   
         ! Copy xlocal out
         j = 1
         do i = offset, offset+m-1
            upd(index(i)) = upd(index(i)) + xlocal(j)
            j = j + 1
         end do
      else
         !!! Single rhs, direct update
         j = 1
         if (transa == 'C') then
            do i = offset, offset+m-1
               w = czero
               do k = col, col+nelim-1
                  w = w - conjg(dest(j))*rhs(k)
                  j = j + 1
               end do
               j = j + (ldd-nelim)
               upd(index(i)) = upd(index(i)) + w
            end do
         else
            do i = offset, offset+m-1
               w = czero
               do k = col, col+nelim-1
                  w = w - dest(j)*rhs(k)
                  j = j + 1
               end do
               j = j + (ldd-nelim)
               upd(index(i)) = upd(index(i)) + w
            end do
         end if
      endif
   else
      !!! Multiple rhs, BLAS 3
      call cgemm(transa, 'N', m, nrhs, nelim, -cone, dest, ldd, rhs(col), ldr, &
         czero, xlocal, m)
   
      ! Copy xlocal out
      j = 1
      do i = offset, offset+m-1
         do r = 0, nrhs-1
            upd(index(i)+r*ldu) = upd(index(i)+r*ldu) + xlocal(j+r*m)
         end do
         j = j + 1
      end do
   endif

!%%%if(control%unit_log.gt.0) then
!%%%   this_thread = 0
!%%%!$    this_thread = omp_get_thread_num()
!%%%   call system_clock(t_end)
!%%%   call log_task(control, this_thread, t_start, t_end, "FU", id)
!%%%endif
end subroutine slv_fwd_update

!*************************************************

!
! TASK_BUPD
! B_i <- B_i - L_ij^-* B_j
!
! Hermitian case:
! conjg(L_ij) is stored by rows (that is, we hold L_ij^H by cols so no
! conjugate transpose needed in blas calls)
! Symmetric case:
! L_ij is stored by rows (that is, we hold L_ij^T by cols so no
! transpose needed in blas calls)
!
subroutine slv_bwd_update(m, nelim, col, offset, index, dest, ldd, nrhs, rhs, &
      upd, ldr, xlocal, control, id)
   integer, intent(in) :: m ! number of rows in block
   integer, intent(in) :: nelim ! number eliminations (immediate return if =0)
   integer, intent(in) :: col ! start of block column variables in rhs
   integer, intent(in) :: offset ! offset into index we start at
   integer, dimension(*), intent(in) :: index
   integer, intent(in) :: ldd ! leading dimension of dest
   complex(wp), dimension(ldd*m), intent(in) :: dest ! holds L_ij
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldr  ! leading extent of rhs
   complex(wp), intent(inout) :: rhs(ldr*nrhs)
   complex(wp), intent(inout) :: upd(ldr*nrhs)
   complex(wp), dimension(*), intent(out) :: xlocal
   type(MA86_control), intent(in) :: control
   integer(long), intent(in) :: id

   integer :: i
   integer :: j
   integer :: k
   integer :: r ! right hand side loop variable
   complex(wp) :: w ! temporary work variable
!%%%integer :: t_start, t_end,this_thread

   if(nelim.eq.0) return

!%%%if(control%unit_log.gt.0) call system_clock(t_start)

   ! backward substitution
   if(nrhs.eq.1) then
      if(m-nelim.gt.10 .and. nelim.gt.4) then
         !!! Single right-hand side, BLAS 2

         ! Copy xlocal in
         j = 1
         do i = offset, offset+m-1
            xlocal(j) = rhs(index(i))
            j = j + 1
         end do

         call cgemv('N', nelim, m, -cone, dest, ldd, xlocal, 1, cone, &
            upd(col), 1)
      else
         !!! Single right-hand side, direct update
         j = 1
         do i = offset, offset+m-1
            w = rhs(index(i))
            do k = col, col + nelim - 1
               upd(k) = upd(k) - dest(j)*w
               j = j + 1
            end do
            j = j + (ldd-nelim)
         end do
      endif
   else
      !!! Multiple RHS, BLAS 3

      ! Copy xlocal in
      j = 1
      do i = offset, offset+m-1
         do r = 0, nrhs-1
            xlocal(j+r*m) = rhs(index(i)+r*ldr)
         end do
         j = j + 1
      end do

      call cgemm('N', 'N', nelim, nrhs, m, -cone, dest, ldd, xlocal, m, &
         cone, upd(col), ldr)
   endif

!%%%if(control%unit_log.gt.0) then
!%%%   this_thread = 0
!%%%!$    this_thread = omp_get_thread_num()
!%%%   call system_clock(t_end)
!%%%   call log_task(control, this_thread, t_start, t_end, "BU", id)
!%%%endif
end subroutine slv_bwd_update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debugging, logging and error handling routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Log a task
!
!%%%subroutine log_task(control, thread, start, finish, task, opt1, opt2, opt3)
!%%%   type(ma86_control), intent(in) :: control
!%%%   integer, intent(in) :: thread
!%%%   integer, intent(in) :: start
!%%%   integer, intent(in) :: finish
!%%%   character(len=2), intent(in) :: task
!%%%   integer(long), optional, intent(in) :: opt1
!%%%   integer(long), optional, intent(in) :: opt2
!%%%   integer(long), optional, intent(in) :: opt3
!%%%
!%%%   integer :: arg
!%%%
!%%%   arg = 0
!%%%   if(present(opt1)) arg = arg + 1
!%%%   if(present(opt2)) arg = arg + 1
!%%%   if(present(opt3)) arg = arg + 1
!%%%
!%%%   select case(arg)
!%%%   case(0)
!%%%      write(control%unit_log, "(i4, 2i12, 1x, a2)") &
!%%%         thread, start, finish, task
!%%%   case(1)
!%%%      write(control%unit_log, "(i4, 2i12, 1x, a2, 1x, 1i12)") &
!%%%         thread, start, finish, task, opt1
!%%%   case(2)
!%%%      write(control%unit_log, "(i4, 2i12, 1x, a2, 1x, 2i12)") &
!%%%         thread, start, finish, task, opt1, opt2
!%%%   case(3)
!%%%      write(control%unit_log, "(i4, 2i12, 1x, a2, 1x, 3i12)") &
!%%%         thread, start, finish, task, opt1, opt2, opt3
!%%%   end select
!%%%end subroutine log_task

!*************************************************

!
! Converts a task into a character representation thereof
! Beware not to cause nested i/O !!!
!
!function print_job(task)
!   character(len=26) :: print_job
!   type(dagtask), intent(in) :: task
!
!   select case(task%task_type)
!   case(TASK_FACTORIZE_COLUMN)
!      write(print_job, "(a,'(',i5,')')") "factorize_column", task%dest
!   case(TASK_SOLVE_BLOCK)
!      write(print_job, "(a,'(',i5,',',i5,')')") "solve_block", &
!      task%dest, task%src1
!   case(TASK_UPDATE_INTERNAL)
!      write(print_job, "(a,'(',i5,2(',',i5),')')") "within", task%dest, &
!         task%src1, task%src2
!   case(TASK_UPDATE_BETWEEN)
!      write(print_job, "(a,'(',i5,',',i5,')')") "between", &
!         task%dest, task%src1
!   case(TASK_SLV_FSLV)
!      write(print_job, "(a,'(',i5,')')") "slv_fslv", &
!         task%dest
!   case(TASK_SLV_FUPD)
!      write(print_job, "(a,'(',i5,')')") "slv_fupd", &
!         task%dest
!   case(TASK_SLV_BSLV)
!      write(print_job, "(a,'(',i5,')')") "slv_bslv", &
!         task%dest
!   case(TASK_SLV_BUPD)
!      write(print_job, "(a,'(',i5,')')") "slv_bupd", &
!         task%dest
!   case(TASK_DONE)
!      write(print_job, "(a)") "FINISH"
!   case(TASK_NONE)
!      write(print_job, "(a)") "WAIT"
!   case default
!      write(print_job, "(a)") "UNKNOWN TASK"
!   end select
!end function print_job

!*************************************************

subroutine ma86_print_flag(iflag, control, context, st)
   integer, intent(in) :: iflag
   type(ma86_control), intent(in) :: control
   integer, intent(in), optional :: st
   ! context: is an optional assumed size character array of intent(in).
   ! It describes the context under which the error occured
   character (len=*), optional, intent(in) :: context

   integer :: nout

   if (iflag < 0) then
      nout = control%unit_error
      if (control%diagnostics_level < 0) nout = -1
      if (nout < 0) return
      write (nout,'(/3a,i3)') ' Error return from ',trim(context),&
         '. Error flag = ', iflag
   else
      nout = control%unit_warning
      if (control%diagnostics_level < 0) nout = -1
      if (nout < 0) return
      write (nout,'(/3a,i3)') ' Warning from ',trim(context),&
         '. Warning flag = ', iflag
   end if

   ! Errors
   select case(iflag)
   case(MA86_ERROR_ALLOCATION)
      if (present(st)) write (nout,'(a,i8)') &
         ' Allocation error. stat parameter = ', st
   case(MA86_ERROR_JOB_OOR)
      write (nout,'(a)') ' job out of range.'
   case(MA86_ERROR_ORDER)
      write (nout,'(a)') ' Error in user-supplied elimination order'
   case(MA86_ERROR_X_SIZE)
      write (nout,'(a)') ' Error in size of x. lx or nrhs too small'
   case(MA86_ERROR_SINGULAR)
      write (nout,'(a)') &
    ' Error matrix is singular and control%action=.false'
   case(MA86_ERROR_INFINITY)
      write (nout,'(a)') ' IEEE infinities found in factorization'
   case(MA86_ERROR_STATIC_SMALL)
      write (nout,'(a)') ' Error in control%static'
   case(MA86_ERROR_MATRIX_TYPE)
      write (nout,'(a)') ' matrix_type is invalid'

   ! Warnings

   case(MA86_WARNING_SINGULAR)
      write (nout,'(a)') ' Matrix found to be singular'
   case(MA86_WARNING_POOL_SING)
      write (nout,'(a)') &
    ' Matrix found to be singular and task pool too small'
   case(MA86_WARNING_POOL_SMALL)
      write (nout,'(a)') ' Task pool too small'

   ! Unknown error
   case default
      write (nout,'(a)') ' Unexpected Error. Please report.'
   end select

end subroutine MA86_print_flag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! routines to help interfaces get at private data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure integer function ma86_get_n_complex(keep)
   type(ma86_keep), intent(in) :: keep
   ma86_get_n_complex = keep%n
end function ma86_get_n_complex

end module hsl_MA86_complex
