function hsl_ma86_expert(handle)
% HSL_MA86_DESTROY  Free memory associated with factorization.
%     hsl_ma86_destroy(handle) will free all memory associated with handle.
%     The factorization may not be reused again.
%
%     Usage: hsl_ma86_destroy(handle)
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
%     See also: ma86_backslash, ma86_factor, ma86_solve

hsl_ma86_expert('destroy', handle)
