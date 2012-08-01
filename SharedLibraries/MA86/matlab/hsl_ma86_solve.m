function [X, info] = hsl_ma86_solve(handle, B, varargin)
% HSL_MA86_SOLVE  Sparse Symmetric Indefinite Solve.
%     X = hsl_ma86_solve(handle, B) solves the equation AX=B for X given
%     precomputed factors associated with handle. The handle must have been
%     obtained by a prior call to hsl_ma86_factor of hsl_ma86_backslash.
%
%     Usage: X = hsl_ma86_solve(handle, B)
%            [X, info] = hsl_ma86_solve(handle, B, control)
%
%     The optional argument CONTROL may have the following components set. If
%     they are not set then the stated default is used.
%     control.num_threads  - Number of threads on which to run. Default is the
%                            maximum available.
%
%     The optional return value INFO will have some of the following components
%     set on exit.
%     info.solve_time         - Wall clock time for Fortran ma86_solve call
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
%     See also: ma86_backslash, ma86_destroy, ma86_factor

optargin = size(varargin,2);
if(optargin == 0)
   [X, info] = hsl_ma86_expert('solve', handle, B);
elseif(optargin == 1)
   [X, info] = hsl_ma86_expert('solve', handle, B, varargin{1});
else
   error ('Too many arguments')
end
