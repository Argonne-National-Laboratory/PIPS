function varargout = hsl_ma86_expert(varargin)
% HSL_MA86  Sparse Symmetric Indefinite Linear Solver.
%     hsl_ma86_expert is a direct Fortran mex interface. The first parameter
%     specifies the action to be performed. Inputs and outputs depend on the
%     action.
%
%     Complex matrices are assumed to be symmetric unless they are specified as
%     Hermitian by setting the value control.hermitian=true.
%
%     [handle, info] = hsl_ma86_expert('factor', A, control, P)
%        Performs the factorization of a symmetric matrix A and returns an
%        integer handle for the factorization. This is equivalent to calling
%        the Fortran routines ma86_analyse and ma86_factor.
%        The argument control is optional and is described below.
%        The argument P is optional and is a vector as returned by e.g.
%           symamd(A). If it is not present hsl_ma86 will find its own
%           fill-reducing permutation.
%        The argument info is optional and is described below.
%
%     [X, info] = hsl_ma86_expert('solve', handle, B, control)
%        Solves the equation AX = B for X using a factorization previously
%        computed by either the 'factor' or 'backslash' actions.
%        The arguments control and info are optional and are described below.
%
%     [X, info, handle] = hsl_ma86_expert('backslash', A, B, control, P)
%        Combines the 'factor' and 'solve' actions in a more efficient fashion.
%        It solves the equation AX = B for X using a factorization it computes.
%        The arguments control and info are optional and are described below.
%        The argument P is optional and as described for the 'factor' option
%        above.
%        The argument handle is optional, and if supplied then may be used to
%        reference the factorization in future calls to 'solve'.
%
%     hsl_ma86_expert('destroy', handle)
%        Destroys the factorization referecenced by handle, freeing the memory
%        associated with it.
%
%     The optional argument CONTROL may have the following components set. If
%     they are not set then the stated default is used.
%     control.hermitian    - True or false. Determines if a complex matrix is
%                            treated as Hermitian (true) or symmetric (false).
%                            Default is false.
%     control.nb           - Block size to be used. Default is 256.
%     control.nemin        - Maximum number of columns in candidates for
%                            supernode amalgamation. Default is 32.
%     control.num_threads  - Number of threads on which to run. Default is the
%                            maximum available.
%     control.scaling      - Determines if scaling is to be used with values:
%                                   1 : MC77 in the one norm
%                            otherwise:  no scaling
%                            Default is 1.
%     control.small        - Pivots of modulus less than this are treated as
%                            zero. Default is 1e-20.
%     control.static       - If greater than zero static pivoting is used.
%                            Default is 0.0.
%     control.u            - Initial relative pivot tolerance threshold. Default
%                            is 0.01.
%     control.umin         - Relaxed relative pivot tolerance threshold. Default
%                            is 1.0.
%
%     The optional return value INFO will have some of the following components
%     set on exit.
%     info.matrix_rank        - Number of non-zero pivots.
%     info.num_delay          - Number of delayed pivots.
%     info.num_factor         - Number of entries in the factors (after
%                               supernode amalgamation and pivoting).
%     info.num_flops          - Number of floating point operations to form
%                               factors (after supernode amalgamation and
%                               pivoting).
%     info.num_neg            - Number of negative pivots in factors.
%     info.num_perturbed      - Number of perturbed pivots when static pivoting
%                               is used.
%     info.num_two            - Number of 2x2 pivots used in factorization.
%     info.order              - Ordering used. One of 'MeTiS', 'AMD' or 'user'.
%                               If a user has provided a permutation that is
%                               used, otherwise a MeTiS ordering is used.
%                               Should MeTiS not be available then AMD is used
%                               instead.
%     info.usmall             - Threshold parameter actually used. This will
%                               only differ from control.u if control.umin <
%                               control.u.
%     info.order_time         - Wall clock time for Fortran ordering routine
%     info.analyse_time       - Wall clock time for Fortran ma86_analyse call
%     info.factor_time        - Wall clock time for Fortran ma86_factor call
%     info.factor_solve_time  - Wall clock time for Fortran ma86_factor_solve
%                               call
%     info.solve_time         - Wall clock time for Fortran ma86_solve call

error('hsl_ma86 must be compiled before use')
