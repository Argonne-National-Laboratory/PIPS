function [X, varargout] = ma86_backslash(A, B, varargin)
% HSL_MA86_BACKSLASH Sparse Symmetric Indefinite Linear Solve.
%     X = hsl_ma86_backslash(A, B) solves the equation AX=B for X by means of
%     a symmetric indefinite factorization PAP' = LDL'. P is generated
%     automatically to reduce fill-in. A is assumed to be symmetric and only
%     the lower triangular part is referenced. A must be square.
%     Complex matrices are assumed to be symmetric unless they are specified as
%     Hermitian by setting the value control.hermitian=true.
%
%     Usage: X = hsl_ma86_backslash(A, B)
%            [X, info, handle] = ma86_backslash(A, B, control, P)
%
%     control is a structure described below. P is a permuaton such as that
%     output from symamd(A). info is a structure described below. handle
%     provides a reference to the factorization so it may be reused at a later
%     date.
% 
%     control may have the following components set. If they are not set then
%     the stated default is used.
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
%     On return, info will have the following components set.
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
%     info.factor_solve_time  - Wall clock time for Fortran ma86_factor_solve
%                               call
%
%     Note that if the optional return value handle is used, then memory will
%     be tied up storing the factorization. This may be recovered by calling
%     ma86_destroy.
%
%     Please cite HSL as:
%     [1] HSL, a collection of Fortran codes for large-scale scientific
%         computation. See http://www.hsl.rl.ac.uk/.
%
%     This code is described in
%     [2] An indefinite sparse direct solver for large problems on multicore
%         machines. J.D. Hogg and J.A. Scott. Technical Report RAL-TR-2010-011.
%     [3] A modern analyse phase for sparse tree-based direct method.
%         J.D. Hogg and J.A. Scott. Technical Report RAL-TR-2010-031.
%     [4] Design of a multicore sparse Cholesky Factorization using DAGs.
%         J.D. Hogg, J.K. Reid and J.A. Scott.
%         Siam J. Scientific Computing 32(6) pp 3627--3649 (2010)
%
%     See also: ma86_destroy, ma86_factor, ma86_solve

optargin = size(varargin,2);
optargout = max(nargout,1)-1;
if(optargin>2)
   error ('Too many arguments')
end
if(optargout == 0)
   if(optargin == 0)
      X = hsl_ma86_expert('backslash', A, B);
   elseif(optargin == 1)
      X = hsl_ma86_expert('backslash', A, B, varargin{1});
   elseif(optargin == 2)
      X = hsl_ma86_expert('backslash', A, B, varargin{1}, varargin{2});
   end
elseif(optargout == 1)
   if(optargin == 0)
      [X, info] = hsl_ma86_expert('backslash', A, B);
   elseif(optargin == 1)
      [X, info] = hsl_ma86_expert('backslash', A, B, varargin{1});
   elseif(optargin == 2)
      [X, info] = ...
         hsl_ma86_expert('backslash', A, B, varargin{1}, varargin{2});
   end
   varargout(1) = {info};
elseif(optargout == 2)
   if(optargin == 0)
      [X, info, handle] = hsl_ma86_expert('backslash', A, B);
   elseif(optargin == 1)
      [X, info, handle] = ...
         hsl_ma86_expert('backslash', A, B, varargin{1});
   elseif(optargin == 2)
      [X, info, handle] = ...
         hsl_ma86_expert('backslash', A, B, varargin{1}, varargin{2});
   end
   varargout(1) = {info};
   varargout(2) = {handle};
else
    error ('Too many output arguments')
end
