/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    exprinterpret_cppad.cpp
 * @brief   methods to interpret (evaluate) an expression tree "fast" using CppAD
 * @ingroup EXPRINTS
 * @author  Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "nlpi/pub_expr.h"
#include "nlpi/exprinterpret.h"

#include <cmath>
#include <vector>
using std::vector;

/* CppAD is not thread-safe by itself, but uses some static datastructures
 * To run it in a multithreading environment, a special CppAD memory allocator that is aware of the multiple threads has to be used.
 * This allocator requires to know the number of threads and a thread number for each thread.
 * Conveniently, SCIPs message handler has currently still the same issue, so we can use its routines here.
 */
#ifndef NPARASCIP
#include "scip/message.h"

/** CppAD needs to know a fixed upper bound on the number of threads at compile time. */
#define CPPAD_MAX_NUM_THREADS 48
#endif

/** sign of a value (-1 or +1)
 * 
 * 0.0 has sign +1
 */
#define SIGN(x) ((x) >= 0.0 ? 1.0 : -1.0)

/* in order to use intervals as operands in CppAD,
 * we need to include the intervalarith.hpp very early and require the interval operations to be in the CppAD namespace */
#define SCIPInterval_NAMESPACE CppAD
#include "nlpi/intervalarith.h"

SCIP_Real CppAD::SCIPInterval::infinity = SCIP_DEFAULT_INFINITY;
using CppAD::SCIPInterval;

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wshadow"
#endif

#include <cppad/cppad.hpp>
#ifndef CPPAD_PACKAGE_STRING
#include <cppad/config.h>
#define CPPAD_PACKAGE_STRING PACKAGE_STRING
#endif
#include <cppad/error_handler.hpp>

#ifdef __GNUC__
#pragma GCC diagnostic warning "-Wshadow"
#endif

#ifndef NPARASCIP

/** CppAD callback function that indicates whether we are running in parallel mode */
static bool in_parallel(void)
{
   return SCIPmessagehdlrGetNThreads() > 0;
}

/** CppAD callback function that returns the number of the current thread */
size_t thread_num(void)
{
   return SCIPmessagehdlrGetThreadNum();
}

/** sets up CppAD's datastructures for running in multithreading mode
 * it must be called once before multithreading is started
 */
static char init_parallel(void)
{
   CppAD::thread_alloc::parallel_setup(CPPAD_MAX_NUM_THREADS, in_parallel, thread_num);
   CppAD::parallel_ad<double>();
   CppAD::parallel_ad<SCIPInterval>();

   return 0;
}

/** a dummy variable that can is initialized to the result of init_parallel
 * the purpose is to make sure that init_parallel() is called before any multithreading is started
 */
static char init_parallel_return = init_parallel();

#endif // NPARASCIP

/* Brad recomends using the discrete function feature of CppAD for sign, since it avoids the need for retaping
 * It can be used since it's derivative is almost everywhere 0.0 */

/* sign as function for double */
double sign(const double &x)
{
   return SIGN(x);
}
/* discrete CppAD function sign(double) for use in eval */
CPPAD_DISCRETE_FUNCTION(double, sign)

/* sign as function for SCIPInterval
 * this time outside of the CppAD namespace
 */
SCIPInterval sign(const SCIPInterval& x)
{
   SCIPInterval resultant;

   SCIPintervalSign(&resultant, x);

   return resultant;
}
/* discrete CppAD function sign(SCIPInterval) for use in eval */
CPPAD_DISCRETE_FUNCTION(SCIPInterval, sign)

/** defintion of CondExpOp for SCIPInterval (required by CppAD) */
inline
SCIPInterval CondExpOp(
   enum CppAD::CompareOp cop,
   const SCIPInterval&   left,
   const SCIPInterval&   right,
   const SCIPInterval&   trueCase,
   const SCIPInterval&   falseCase)
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "SCIPInterval CondExpOp(...)",
      "Error: cannot use CondExp with an interval type"
      );

   return SCIPInterval();
}

/** another function that returns whether two intervals are the same (required by CppAD) */
inline
bool EqualOpSeq(
   const SCIPInterval&   x,                  /**< first operand */
   const SCIPInterval&   y                   /**< second operand */
   )
{
   return x == y;
}

/** another function required by CppAD */
inline
bool IdenticalPar(
   const SCIPInterval&   x                   /**< operand */
   )
{
   return true;
}

/** returns whether the interval equals [0,0] */
inline
bool IdenticalZero(
   const SCIPInterval&   x                   /**< operand */
   )
{
   return (x == 0.0);
}

/** returns whether the interval equals [1,1] */
inline
bool IdenticalOne(
   const SCIPInterval&   x                   /**< operand */
   )
{
   return (x == 1.0);
}

/** yet another function that checks whether two intervals are equal */
inline
bool IdenticalEqualPar(
   const SCIPInterval&   x,                  /**< first operand */
   const SCIPInterval&   y                   /**< second operand */
   )
{
   return (x == y);
}

/** greater than zero not defined for intervals */
inline
bool GreaterThanZero(
   const SCIPInterval&   x                   /**< operand */
   )
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "GreaterThanZero(x)",
      "Error: cannot use GreaterThanZero with interval"
      );

   return false;
}

/** greater than or equal zero not defined for intervals */
inline
bool GreaterThanOrZero(
   const SCIPInterval&   x                   /**< operand */
   )
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__ ,
      "GreaterThanOrZero(x)",
      "Error: cannot use GreaterThanOrZero with interval"
      );

   return false;
}

/** less than not defined for intervals */
inline
bool LessThanZero(
   const SCIPInterval&   x                   /**< operand */
   )
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "LessThanZero(x)",
      "Error: cannot use LessThanZero with interval"
      );

   return false;
}

/** less than or equal not defined for intervals */
inline
bool LessThanOrZero(
   const SCIPInterval&   x                   /**< operand */
   )
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "LessThanOrZero(x)",
      "Error: cannot use LessThanOrZero with interval"
      );

   return false;
}

/** conversion to integers not defined for intervals */
inline
int Integer(
   const SCIPInterval&   x                   /**< operand */
   )
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "Integer(x)",
      "Error: cannot use Integer with interval"
      );

   return 0;
}

/** printing of an interval (required by CppAD) */
inline
std::ostream& operator<<(std::ostream& out, const SCIP_INTERVAL& x)
{
   out << '[' << x.inf << ',' << x.sup << ']';
   return out;
}

using CppAD::AD;

/** expression interpreter */
struct SCIP_ExprInt
{
   BMS_BLKMEM*           blkmem;             /**< block memory data structure */
};

/** expression specific interpreter data */
class SCIP_ExprIntData
{
public:
   /* constructor */
   SCIP_ExprIntData()
      : need_retape(true), int_need_retape(true), need_retape_always(false), blkmem(NULL), root(NULL)
   { }

   /* destructor */
   ~SCIP_ExprIntData()
   { }

   vector< AD<double> >  X;                  /**< vector of dependent variables */
   vector< AD<double> >  Y;                  /**< result vector */ 
   CppAD::ADFun<double>  f;                  /**< the function to evaluate as CppAD object */

   vector<double>        x;                  /**< current values of dependent variables */
   double                val;                /**< current function value */
   bool                  need_retape;        /**< will retaping be required for the next point evaluation? */

   vector< AD<SCIPInterval> > int_X;         /**< interval vector of dependent variables */
   vector< AD<SCIPInterval> > int_Y;         /**< interval result vector */
   CppAD::ADFun<SCIPInterval> int_f;         /**< the function to evaluate on intervals as CppAD object */

   vector<SCIPInterval>  int_x;              /**< current interval values of dependent variables */
   SCIPInterval          int_val;            /**< current interval function value */
   bool                  int_need_retape;    /**< will retaping be required for the next interval evaluation? */

   bool                  need_retape_always; /**< will retaping be always required? */

   BMS_BLKMEM*           blkmem;             /**< block memory used to allocate expresstion tree */
   SCIP_EXPR*            root;               /**< copy of expression tree; @todo we should not need to make a copy */
};

#ifdef CPPAD_USER_ATOMIC
/** forward sweep of positive integer power
 * Given the taylor coefficients for x, we have to compute the taylor coefficients for f(x),
 * that is, given tx = (x, x', x'', ...), we compute the coefficients ty = (y, y', y'', ...)
 * in the taylor expansion of f(x) = x^p.
 * Thus, y   = x^p
 *           = tx[0]^p,
 *       y'  = p * x^(p-1) * x'
 *           = p * tx[0]^(p-1) * tx[1],
 *       y'' = p * (p-1) * x^(p-2) * x'^2 + p * x^(p-1) * x''
 *           = p * (p-1) * tx[0]^(p-2) * tx[1]^2 + p * tx[0]^(p-1) * tx[2]
 */
template<class Type>
bool forward_posintpower(
   size_t                      id,           /**< user data identifier, we use it to store exponent p */
   size_t                      k,            /**< derivative order that shall be computed */
   size_t                      n,            /**< number of variables, should be 1 */
   size_t                      m,            /**< dimension of function, should be 1 */
   const CppAD::vector<bool>&  vx,           /**< indicates whether argument is a variable, or empty vector */
   CppAD::vector<bool>&        vy,           /**< vector to store which function values depend on variables, or empty vector */
   const CppAD::vector<Type>&  tx,           /**< values for taylor coefficients of x */
   CppAD::vector<Type>&        ty            /**< vector to store taylor coefficients of y */
   )
{
   assert(id > 1);
   assert(n == 1);
   assert(m == 1);
   assert(tx.size() >= k);
   assert(ty.size() >= k);

   if( vx.size() > 0 )
   {
      assert(vx.size() == 1);
      assert(vy.size() == 1);
      assert(k == 0);

      vy[0] = vx[0];
   }

   switch( k )
   {
   case 0:
      ty[0] = pow(tx[0], id);
      break;

   case 1:
      ty[1] = pow(tx[0], id-1) * tx[1];
      ty[1] *= double(id);
      break;

   case 2:
      if( id > 2 )
      {
         // ty[2] = id * (id-1) * pow(tx[0], id-2) * tx[1] * tx[1] + id * pow(tx[0], id-1) * tx[2];
         ty[2]  = pow(tx[0], id-2) * tx[1] * tx[1];
         ty[2] *= id-1;
         ty[2] += pow(tx[0], id-1) * tx[2];
         ty[2] *= id;
      }
      else
      {
         assert(id == 2);
         // ty[2] = id * tx[1] * tx[1] + id * tx[0] * tx[2];
         ty[2]  = tx[1] * tx[1] + tx[0] * tx[2];
         ty[2] *= id;
      }
      break;

   default:
      return false;
   }

   return true;
}

/** reverse sweep of positive integer power
 * Assume y(x) is a function of the taylor coefficients of f(x) = x^p for x, i.e.,
 *   y(x) = [ x^p, p * x^(p-1) * x', p * (p-1) * x^(p-2) * x'^2 + p * x^(p-1) * x'', ... ].
 * Then in the reverse sweep we have to compute the elements of \partial h / \partial x^[l], l = 0, ..., k,
 * where x^[l] is the l'th taylor coefficient (x, x', x'', ...) and h(x) = g(y(x)) for some function g:R^k -> R.
 * That is, we have to compute
 * px[l] = \partial h / \partial x^[l] = (\partial g / \partial y) * (\partial y / \partial x^[l])
 *       = \sum_{i=0}^k (\partial g / \partial y_i) * (\partial y_i / \partial x^[l])
 *       = \sum_{i=0}^k py[i] * (\partial y_i / \partial x^[l])
 *
 * For k = 0, this means
 * px[0] = py[0] * (\partial y_0 / \partial x^[0])
 *       = py[0] * (\partial x^p / \partial x)
 *       = py[0] * p * tx[0]^(p-1)
 *
 * For k = 1, this means
 * px[0] = py[0] * (\partial y_0 / \partial x^[0]) + py[1] * (\partial y_1 / \partial x^[0])
 *       = py[0] * (\partial x^p / \partial x)     + py[1] * (\partial (p * x^(p-1) * x') / \partial x)
 *       = py[0] * p * tx[0]^(p-1)                 + py[1] * p * (p-1) * tx[0]^(p-2) * tx[1]
 * px[1] = py[0] * (\partial y_0 / \partial x^[1]) + py[1] * (\partial y_1 / \partial x^[1])
 *       = py[0] * (\partial x^p / \partial x')    + py[1] * (\partial (p * x^(p-1) x') / \partial x')
 *       = py[0] * 0                               + py[1] * p * tx[0]^(p-1)
 */
template<class Type>
bool reverse_posintpower(
   size_t                      id,           /**< user data identifier, we use it to store exponent p */
   size_t                      k,            /**< derivative order that shall be computed */
   size_t                      n,            /**< number of variables, should be 1 */
   size_t                      m,            /**< dimension of function, should be 1 */
   const CppAD::vector<Type>&  tx,           /**< values for taylor coefficients of x */
   const CppAD::vector<Type>&  ty,           /**< values for taylor coefficients of y */
   CppAD::vector<Type>&        px,           /**< vector to store partial derivatives of h(x) = g(y(x)) w.r.t. x */
   const CppAD::vector<Type>&  py            /**< values for partial derivatives of g(x) w.r.t. y */
   )
{
   assert(id > 1);
   assert(n == 1);
   assert(m == 1);
   assert(px.size() >= k);
   assert(py.size() >= k);
   assert(tx.size() >= k);

   switch( k )
   {
   case 0:
      // px[0] = py[0] * id * pow(tx[0], id-1);
      px[0]  = py[0] * pow(tx[0], id-1);
      px[0] *= id;
      break;

   case 1:
      // px[0] = py[0] * id * pow(tx[0], id-1) + py[1] * id * (id-1) * pow(tx[0], id-2) * tx[1];
      px[0]  = py[1] * tx[1] * pow(tx[0], id-2);
      px[0] *= id-1;
      px[0] += py[0] * pow(tx[0], id-1);
      px[0] *= id;
      // px[1] = py[1] * id * pow(tx[0], id-1);
      px[1]  = py[1] * pow(tx[0], id-1);
      px[1] *= id;
      break;

   default:
      return false;
   }

   return true;
}

/** computes sparsity of jacobian during a forward sweep
 * For a 1 x q matrix R, we have to return the sparsity pattern of the 1 x q matrix S(x) = f'(x) * R.
 * Since f'(x) is dense, the sparsity of S will be the sparsity of R.
 */
static
bool for_jac_sparse_posintpower(
   size_t                                  id, /**< user data identifier, we use it to store exponent p */
   size_t                                  n,  /**< number of variables, should be 1 */
   size_t                                  m,  /**< dimension of function, should be 1 */
   size_t                                  q,  /**< number of columns in R */
   const CppAD::vector<std::set<size_t> >& r,  /**< sparsity of R, columnwise */
   CppAD::vector<std::set<size_t> >&       s   /**< vector to store sparsity of S, columnwise */
   )
{
   assert(n == 1);
   assert(m == 1);
   assert(id > 1);
   assert(r.size() == 1);
   assert(s.size() == 1);

   s[0] = r[0];

   return true;
}

/** computes sparsity of jacobian during a reverse sweep
 * For a q x 1 matrix S, we have to return the sparsity pattern of the q x 1 matrix R(x) = S * f'(x).
 * Since f'(x) is dense, the sparsity of R will be the sparsity of S.
 */
static
bool rev_jac_sparse_posintpower(
   size_t                                  id, /**< user data identifier, we use it to store exponent p */
   size_t                                  n,  /**< number of variables, should be 1 */
   size_t                                  m,  /**< dimension of function, should be 1 */
   size_t                                  q,  /**< number of rows in R */
   CppAD::vector<std::set<size_t> >&       r,  /**< sparsity of R, rowwise */
   const CppAD::vector<std::set<size_t> >& s   /**< vector to store sparsity of S, rowwise */
   )
{
   assert(n == 1);
   assert(m == 1);
   assert(id > 1);
   assert(r.size() == 1);
   assert(s.size() == 1);

   r[0] = s[0];

   return true;
}

/** computes sparsity of hessian during a reverse sweep
 * Assume V(x) = (g(f(x)))'' R  with f(x) = x^p for a function g:R->R and a matrix R.
 * we have to specify the sparsity pattern of V(x) and T(x) = (g(f(x)))'.
 */
static
bool rev_hes_sparse_posintpower(
   size_t                                  id, /**< user data identifier, we use it to store exponent p */
   size_t                                  n,  /**< number of variables, should be 1 */
   size_t                                  m,  /**< dimension of function, should be 1 */
   size_t                                  q,  /**< number of columns in S and R */
   const CppAD::vector<std::set<size_t> >& r,  /**< sparsity pattern of R */
   const CppAD::vector<bool>&              s,  /**< sparsity pattern of S = g'(y) */
   CppAD::vector<bool>&                    t,  /**< vector to store sparsity pattern of T(x) = (g(f(x)))' */
   const CppAD::vector<std::set<size_t> >& u,  /**< sparsity pattern of U(x) = g''(f(x)) f'(x) R */
   CppAD::vector< std::set<size_t> >&      v   /**< vector to store sparsity pattern of V(x) = (g(f(x)))'' R */
   )
{
   assert(n == 1);
   assert(m == 1);
   assert(id > 0);
   assert(r.size() == 1);
   assert(s.size() == 1);
   assert(t.size() == 1);
   assert(u.size() == 1);
   assert(v.size() == 1);

   // T(x) = g'(f(x)) * f'(x) = S * f'(x), and f' is not identically 0
   t[0] = s[0];

   // V(x) = g''(f(x)) f'(x) f'(x) R + g'(f(x)) f''(x) R
   //      = f'(x) U + S f''(x) R, with f'(x) and f''(x) not identically 0
   v[0] = u[0];
   if( s[0] )
      v[0].insert(r[0].begin(), r[0].end());

   return true;
}

/** tell CppAD about our implementation for x^p, p>=2 integer, for x double-valued */
CPPAD_USER_ATOMIC(
   posintpower               ,
   vector                    ,
   double                    ,
   forward_posintpower       ,
   reverse_posintpower       ,
   for_jac_sparse_posintpower,
   rev_jac_sparse_posintpower,
   rev_hes_sparse_posintpower
   )

/** tell CppAD about our implementation for x^p, p>=2 integer, for x interval-valued */
CPPAD_USER_ATOMIC(
   posintpower               ,
   vector                    ,
   SCIPInterval              ,
   forward_posintpower       ,
   reverse_posintpower       ,
   for_jac_sparse_posintpower,
   rev_jac_sparse_posintpower,
   rev_hes_sparse_posintpower
   )
#else
template<class Type>
void posintpower(size_t exp, vector<Type>& in, vector<Type>& out)
{
   out[0] = pow(in[0], (int)exp);
}
#endif

#ifdef CPPAD_USER_ATOMIC
/** forward sweep of signpower
 * Given the taylor coefficients for x, we have to compute the taylor coefficients for f(x),
 * that is, given tx = (x, x', x'', ...), we compute the coefficients ty = (y, y', y'', ...)
 * in the taylor expansion of f(x) = sign(x)abs(x)^p.
 * Thus, y   = sign(x)abs(x)^p
 *           = sign(tx[0])abs(tx[0])^p,
 *       y'  = p * abs(x)^(p-1) * x'
 *           = p * abs(tx[0])^(p-1) * tx[1],
 *       y'' = p * (p-1) * sign(x) * abs(x)^(p-2) * x'^2 + p * abs(x)^(p-1) * x''
 *           = p * (p-1) * sign(tx[0]) * abs(tx[0])^(p-2) * tx[1]^2 + p * abs(tx[0])^(p-1) * tx[2]
 */
template<class Type>
bool forward_signpower(
   size_t                      id,           /**< user data identifier, we use it to store the pointer to the expression that holds the exponent p */
   size_t                      k,            /**< derivative order that shall be computed */
   size_t                      n,            /**< number of variables, should be 1 */
   size_t                      m,            /**< dimension of function, should be 1 */
   const CppAD::vector<bool>&  vx,           /**< indicates whether argument is a variable, or empty vector */
   CppAD::vector<bool>&        vy,           /**< vector to store which function values depend on variables, or empty vector */
   const CppAD::vector<Type>&  tx,           /**< values for taylor coefficients of x */
   CppAD::vector<Type>&        ty            /**< vector to store taylor coefficients of y */
   )
{
   SCIP_Real p;

   assert(id != 0);
   assert(n == 1);
   assert(m == 1);
   assert(tx.size() >= k);
   assert(ty.size() >= k);

   if( vx.size() > 0 )
   {
      assert(vx.size() == 1);
      assert(vy.size() == 1);
      assert(k == 0);

      vy[0] = vx[0];
   }

   p = SCIPexprGetSignPowerExponent((SCIP_EXPR*)(void*)id);
   assert(p > 1.0);

   switch( k )
   {
   case 0:
      ty[0] = SIGN(tx[0]) * pow(REALABS(tx[0]), p);
      break;

   case 1:
      ty[1] = pow(REALABS(tx[0]), p - 1.0) * tx[1];
      ty[1] *= p;
      break;

   case 2:
      if( p != 2.0 )
      {
         ty[2]  = SIGN(tx[0]) * pow(REALABS(tx[0]), p - 2.0) * tx[1] * tx[1];
         ty[2] *= p - 1.0;
         ty[2] += pow(REALABS(tx[0]), p - 1.0) * tx[2];
         ty[2] *= p;
      }
      else
      {
         // y'' = 2 (sign(x) * x'^2 + |x|*x'') = 2 (sign(tx[0]) * tx[1]^2 + abs(tx[0]) * tx[2])
         ty[2]  = SIGN(tx[0]) * tx[1] * tx[1];
         ty[2] += REALABS(tx[0]) * tx[2];
         ty[2] *= p;
      }
      break;

   default:
      return false;
   }

   return true;
}

/** specialization of forward_signpower template for SCIPinterval
 * @todo try to compute tighter resultants
 */
template<>
bool forward_signpower(
   size_t                             id, /**< user data identifier, we use it to store the pointer to the expression that holds the exponent p */
   size_t                             k,  /**< derivative order that shall be computed */
   size_t                             n,  /**< number of variables, should be 1 */
   size_t                             m,  /**< dimension of function, should be 1 */
   const CppAD::vector<bool>&         vx, /**< indicates whether argument is a variable, or empty vector */
   CppAD::vector<bool>&               vy, /**< vector to store which function values depend on variables, or empty vector */
   const CppAD::vector<SCIPInterval>& tx, /**< values for taylor coefficients of x */
   CppAD::vector<SCIPInterval>&       ty  /**< vector to store taylor coefficients of y */
   )
{
   SCIP_Real p;

   assert(id != 0);
   assert(n == 1);
   assert(m == 1);
   assert(tx.size() >= k);
   assert(ty.size() >= k);

   if( vx.size() > 0 )
   {
      assert(vx.size() == 1);
      assert(vy.size() == 1);
      assert(k == 0);

      vy[0] = vx[0];
   }

   p = SCIPexprGetSignPowerExponent((SCIP_EXPR*)(void*)id);
   assert(p > 1.0);

   switch( k )
   {
   case 0:
      ty[0] = signpow(tx[0], p);
      break;

   case 1:
      ty[1] = pow(abs(tx[0]), p - 1.0) * tx[1];
      ty[1] *= p;
      break;

   case 2:
      if( p != 2.0 )
      {
         ty[2]  = signpow(tx[0], p - 2.0) * square(tx[1]);
         ty[2] *= p - 1.0;
         ty[2] += CppAD::pow(abs(tx[0]), p - 1.0) * tx[2];
         ty[2] *= p;
      }
      else
      {
         // y'' = 2 (sign(x) * x'^2 + |x|*x'') = 2 (sign(tx[0]) * tx[1]^2 + abs(tx[0]) * tx[2])
         ty[2]  = CppAD::sign(tx[0]) * square(tx[1]);
         ty[2] += abs(tx[0]) * tx[2];
         ty[2] *= p;
      }
      break;

   default:
      return false;
   }

   return true;
}

/** reverse sweep of signpower
 * Assume y(x) is a function of the taylor coefficients of f(x) = sign(x)|x|^p for x, i.e.,
 *   y(x) = [ f(x), f'(x), f''(x), ... ].
 * Then in the reverse sweep we have to compute the elements of \partial h / \partial x^[l], l = 0, ..., k,
 * where x^[l] is the l'th taylor coefficient (x, x', x'', ...) and h(x) = g(y(x)) for some function g:R^k -> R.
 * That is, we have to compute
 * px[l] = \partial h / \partial x^[l] = (\partial g / \partial y) * (\partial y / \partial x^[l])
 *       = \sum_{i=0}^k (\partial g / \partial y_i) * (\partial y_i / \partial x^[l])
 *       = \sum_{i=0}^k py[i] * (\partial y_i / \partial x^[l])
 *
 * For k = 0, this means
 * px[0] = py[0] * (\partial y_0 / \partial x^[0])
 *       = py[0] * (\partial f(x) / \partial x)
 *       = py[0] * p * abs(tx[0])^(p-1)
 *
 * For k = 1, this means
 * px[0] = py[0] * (\partial y_0  / \partial x^[0]) + py[1] * (\partial y_1   / \partial x^[0])
 *       = py[0] * (\partial f(x) / \partial x)     + py[1] * (\partial f'(x) / \partial x)
 *       = py[0] * p * abs(tx[0])^(p-1)             + py[1] * p * (p-1) * abs(tx[0])^(p-2) * sign(tx[0]) * tx[1]
 * px[1] = py[0] * (\partial y_0  / \partial x^[1]) + py[1] * (\partial y_1 / \partial x^[1])
 *       = py[0] * (\partial f(x) / \partial x')    + py[1] * (\partial f'(x) / \partial x')
 *       = py[0] * 0                                + py[1] * p * abs(tx[0])^(p-1)
 */
template<class Type>
bool reverse_signpower(
   size_t                      id,           /**< user data identifier, we use it to store the pointer to the expression that holds the exponent p */
   size_t                      k,            /**< derivative order that shall be computed */
   size_t                      n,            /**< number of variables, should be 1 */
   size_t                      m,            /**< dimension of function, should be 1 */
   const CppAD::vector<Type>&  tx,           /**< values for taylor coefficients of x */
   const CppAD::vector<Type>&  ty,           /**< values for taylor coefficients of y */
   CppAD::vector<Type>&        px,           /**< vector to store partial derivatives of h(x) = g(y(x)) w.r.t. x */
   const CppAD::vector<Type>&  py            /**< values for partial derivatives of g(x) w.r.t. y */
   )
{
   SCIP_Real p;

   assert(id != 0);
   assert(n == 1);
   assert(m == 1);
   assert(px.size() >= k);
   assert(py.size() >= k);
   assert(tx.size() >= k);

   p = SCIPexprGetSignPowerExponent((SCIP_EXPR*)(void*)id);
   assert(p > 1.0);

   switch( k )
   {
   case 0:
      // px[0] = py[0] * p * pow(abs(tx[0]), p-1);
      px[0]  = py[0] * pow(REALABS(tx[0]), p - 1.0);
      px[0] *= p;
      break;

   case 1:
      if( p != 2.0 )
      {
         // px[0] = py[0] * p * abs(tx[0])^(p-1) + py[1] * p * (p-1) * abs(tx[0])^(p-2) * sign(tx[0]) * tx[1]
         px[0]  = py[1] * tx[1] * pow(REALABS(tx[0]), p - 2.0) * SIGN(tx[0]);
         px[0] *= p - 1.0;
         px[0] += py[0] * pow(REALABS(tx[0]), p - 1.0);
         px[0] *= p;
         // px[1] = py[1] * p * abs(tx[0])^(p-1)
         px[1]  = py[1] * pow(REALABS(tx[0]), p - 1.0);
         px[1] *= p;
      }
      else
      {
         // px[0] = py[0] * 2.0 * abs(tx[0]) + py[1] * 2.0 * sign(tx[0]) * tx[1]
         px[0]  = py[1] * tx[1] * SIGN(tx[0]);
         px[0] += py[0] * REALABS(tx[0]);
         px[0] *= 2.0;
         // px[1] = py[1] * 2.0 * abs(tx[0])
         px[1]  = py[1] * REALABS(tx[0]);
         px[1] *= 2.0;
      }
      break;

   default:
      return false;
   }

   return true;
}

/** specialization of reverse_signpower for SCIPinterval
 * @todo try to compute tighter resultants
 */
template<>
bool reverse_signpower(
   size_t                              id,           /**< user data identifier, we use it to store the pointer to the expression that holds the exponent p */
   size_t                              k,            /**< derivative order that shall be computed */
   size_t                              n,            /**< number of variables, should be 1 */
   size_t                              m,            /**< dimension of function, should be 1 */
   const CppAD::vector<SCIPInterval>&  tx,           /**< values for taylor coefficients of x */
   const CppAD::vector<SCIPInterval>&  ty,           /**< values for taylor coefficients of y */
   CppAD::vector<SCIPInterval>&        px,           /**< vector to store partial derivatives of h(x) = g(y(x)) w.r.t. x */
   const CppAD::vector<SCIPInterval>&  py            /**< values for partial derivatives of g(x) w.r.t. y */
   )
{
   SCIP_Real p;

   assert(id != 0);
   assert(n == 1);
   assert(m == 1);
   assert(px.size() >= k);
   assert(py.size() >= k);
   assert(tx.size() >= k);

   p = SCIPexprGetSignPowerExponent((SCIP_EXPR*)(void*)id);
   assert(p > 1.0);

   switch( k )
   {
   case 0:
      // px[0] = py[0] * p * pow(abs(tx[0]), p-1);
      px[0]  = py[0] * pow(abs(tx[0]), p - 1.0);
      px[0] *= p;
      break;

   case 1:
      if( p != 2.0 )
      {
         // px[0] = py[0] * p * abs(tx[0])^(p-1) + py[1] * p * (p-1) * abs(tx[0])^(p-2) * sign(tx[0]) * tx[1]
         px[0]  = py[1] * tx[1] * signpow(tx[0], p - 2.0);
         px[0] *= p - 1.0;
         px[0] += py[0] * pow(abs(tx[0]), p - 1.0);
         px[0] *= p;
         // px[1] = py[1] * p * abs(tx[0])^(p-1)
         px[1]  = py[1] * pow(abs(tx[0]), p - 1.0);
         px[1] *= p;
      }
      else
      {
         // px[0] = py[0] * 2.0 * abs(tx[0]) + py[1] * 2.0 * sign(tx[0]) * tx[1]
         px[0]  = py[1] * tx[1] * CppAD::sign(tx[0]);
         px[0] += py[0] * abs(tx[0]);
         px[0] *= 2.0;
         // px[1] = py[1] * 2.0 * abs(tx[0])
         px[1]  = py[1] * abs(tx[0]);
         px[1] *= 2.0;
      }
      break;

   default:
      return false;
   }

   return true;
}

/** computes sparsity of jacobian during a forward sweep
 * For a 1 x q matrix R, we have to return the sparsity pattern of the 1 x q matrix S(x) = f'(x) * R.
 * Since f'(x) is dense, the sparsity of S will be the sparsity of R.
 */
static
bool for_jac_sparse_signpower(
   size_t                                 id,           /**< user data identifier, we use it to store the pointer to the expression that holds the exponent p */
   size_t                                  n,  /**< number of variables, should be 1 */
   size_t                                  m,  /**< dimension of function, should be 1 */
   size_t                                  q,  /**< number of columns in R */
   const CppAD::vector<std::set<size_t> >& r,  /**< sparsity of R, columnwise */
   CppAD::vector<std::set<size_t> >&       s   /**< vector to store sparsity of S, columnwise */
   )
{
   assert(n == 1);
   assert(m == 1);
   assert(id != 0);
   assert(r.size() == 1);
   assert(s.size() == 1);

   s[0] = r[0];

   return true;
}

/** computes sparsity of jacobian during a reverse sweep
 * For a q x 1 matrix S, we have to return the sparsity pattern of the q x 1 matrix R(x) = S * f'(x).
 * Since f'(x) is dense, the sparsity of R will be the sparsity of S.
 */
static
bool rev_jac_sparse_signpower(
   size_t                                 id,  /**< user data identifier, we use it to store the pointer to the expression that holds the exponent p */
   size_t                                  n,  /**< number of variables, should be 1 */
   size_t                                  m,  /**< dimension of function, should be 1 */
   size_t                                  q,  /**< number of rows in R */
   CppAD::vector<std::set<size_t> >&       r,  /**< sparsity of R, rowwise */
   const CppAD::vector<std::set<size_t> >& s   /**< vector to store sparsity of S, rowwise */
   )
{
   assert(n == 1);
   assert(m == 1);
   assert(id != 0);
   assert(r.size() == 1);
   assert(s.size() == 1);

   r[0] = s[0];

   return true;
}

/** computes sparsity of hessian during a reverse sweep
 * Assume V(x) = (g(f(x)))'' R  with f(x) = sign(x)abs(x)^p for a function g:R->R and a matrix R.
 * we have to specify the sparsity pattern of V(x) and T(x) = (g(f(x)))'.
 */
static
bool rev_hes_sparse_signpower(
   size_t                                 id,  /**< user data identifier, we use it to store the pointer to the expression that holds the exponent p */
   size_t                                  n,  /**< number of variables, should be 1 */
   size_t                                  m,  /**< dimension of function, should be 1 */
   size_t                                  q,  /**< number of columns in S and R */
   const CppAD::vector<std::set<size_t> >& r,  /**< sparsity pattern of R */
   const CppAD::vector<bool>&              s,  /**< sparsity pattern of S = g'(y) */
   CppAD::vector<bool>&                    t,  /**< vector to store sparsity pattern of T(x) = (g(f(x)))' */
   const CppAD::vector<std::set<size_t> >& u,  /**< sparsity pattern of U(x) = g''(f(x)) f'(x) R */
   CppAD::vector< std::set<size_t> >&      v   /**< vector to store sparsity pattern of V(x) = (g(f(x)))'' R */
   )
{
   assert(n == 1);
   assert(m == 1);
   assert(id != 0);
   assert(r.size() == 1);
   assert(s.size() == 1);
   assert(t.size() == 1);
   assert(u.size() == 1);
   assert(v.size() == 1);

   // T(x) = g'(f(x)) * f'(x) = S * f'(x), and f' is not identically 0
   t[0] = s[0];

   // V(x) = g''(f(x)) f'(x) f'(x) R + g'(f(x)) f''(x) R
   //      = f'(x) U + S f''(x) R, with f'(x) and f''(x) not identically 0
   v[0] = u[0];
   if( s[0] )
      v[0].insert(r[0].begin(), r[0].end());

   return true;
}

/** tell CppAD about our implementation for signpower for x double-valued */
CPPAD_USER_ATOMIC(
   signpower               ,
   vector                  ,
   double                  ,
   forward_signpower       ,
   reverse_signpower       ,
   for_jac_sparse_signpower,
   rev_jac_sparse_signpower,
   rev_hes_sparse_signpower
   )

/** tell CppAD about our implementation for signpower for x interval-valued */
CPPAD_USER_ATOMIC(
   signpower               ,
   vector                  ,
   SCIPInterval            ,
   forward_signpower       ,
   reverse_signpower       ,
   for_jac_sparse_signpower,
   rev_jac_sparse_signpower,
   rev_hes_sparse_signpower
   )

template<class Type>
void evalSignPower(
   Type&                 resultant,          /**< resultant */
   Type&                 arg,                /**< operand */
   SCIP_EXPR*            expr                /**< expression that holds the exponent */
   )
{
   vector<Type> in(1, arg);
   vector<Type> out(1);

   signpower((size_t)(void*)expr, in, out);

   resultant = out[0];
   return;
}

#else
/** template for evaluation for signpower operator
 * only implemented for real numbers, thus gives error by default
 */
template<class Type>
void evalSignPower(
   Type&                 resultant,          /**< resultant */
   Type&                 arg,                /**< operand */
   SCIP_EXPR*            expr                /**< expression that holds the exponent */
   )
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "evalSignPower()",
      "Error: SignPower not implemented for this value type"
      );
}

/** specialization of signpower evaluation for real numbers
 */
template<>
void evalSignPower(
   CppAD::AD<double>&    resultant,          /**< resultant */
   CppAD::AD<double>&    arg,                /**< operand */
   SCIP_EXPR*            expr                /**< expression that holds the exponent */
   )
{
   SCIP_Real exponent;

   exponent = SCIPexprGetSignPowerExponent(expr);

   if( arg == 0.0 )
      resultant = 0.0;
   else if( arg > 0.0 )
      resultant =  pow( arg, exponent);
   else
      resultant = -pow(-arg, exponent);
}
#endif

/** template for evaluation for minimum operator
 * only implemented for real numbers, thus gives error by default
 * @todo implement own userad function
 */
template<class Type>
void evalMin(
   Type&                 resultant,          /**< resultant */
   Type&                 arg1,               /**< first operand */
   Type&                 arg2                /**< second operand */
   )
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "evalMin()",
      "Error: Min not implemented for this value type"
      );
}

/** specialization of minimum evaluation for real numbers
 */
template<>
void evalMin(
   CppAD::AD<double>&    resultant,          /**< resultant */
   CppAD::AD<double>&    arg1,               /**< first operand */
   CppAD::AD<double>&    arg2                /**< second operand */
   )
{
   resultant = MIN(arg1, arg2);
}

/** template for evaluation for maximum operator
 * only implemented for real numbers, thus gives error by default
 * @todo implement own userad function
 */
template<class Type>
void evalMax(
   Type&                 resultant,          /**< resultant */
   Type&                 arg1,               /**< first operand */
   Type&                 arg2                /**< second operand */
   )
{
   CppAD::ErrorHandler::Call(true, __LINE__, __FILE__,
      "evalMax()",
      "Error: Max not implemented for this value type"
      );
}

/** specialization of maximum evaluation for real numbers
 */
template<>
void evalMax(
   CppAD::AD<double>&    resultant,          /**< resultant */
   CppAD::AD<double>&    arg1,               /**< first operand */
   CppAD::AD<double>&    arg2                /**< second operand */
   )
{
   resultant = MAX(arg1, arg2);
}

/** template for evaluation for square-root operator
 * default is to use the standard sqrt-function
 */
template<class Type>
void evalSqrt(
   Type&                 resultant,          /**< resultant */
   Type&                 arg                 /**< operand */
   )
{
   resultant = sqrt(arg);
}

/** specialization of square-root operator for numbers
 * we perturb the function a little bit so that it's derivatives are defined in 0.0
 */
template<>
void evalSqrt(
   CppAD::AD<double>&    resultant,          /**< resultant */
   CppAD::AD<double>&    arg                 /**< operand */
   )
{
   resultant = sqrt(arg + 1e-20) - 1e-10;
}

/** template for evaluation for absolute value operator
 */
template<class Type>
void evalAbs(
   Type&                 resultant,          /**< resultant */
   Type&                 arg                 /**< operand */
   )
{
   resultant = abs(arg);
}

/** specialization of absolute value evaluation for intervals
 * use sqrt(x^2) for now @todo implement own userad function
 */
template<>
void evalAbs(
   CppAD::AD<SCIPInterval>& resultant,       /**< resultant */
   CppAD::AD<SCIPInterval>& arg              /**< operand */
   )
{
   vector<CppAD::AD<SCIPInterval> > in(1, arg);
   vector<CppAD::AD<SCIPInterval> > out(1);

   posintpower(2, in, out);

   resultant = sqrt(out[0]);
}

/** integer power operation for arbitrary integer exponents */
template<class Type>
void evalIntPower(
   Type&                 resultant,          /**< resultant */
   Type&                 arg,                /**< operand */
   int                   exponent            /**< exponent */
   )
{
   if( exponent > 1 )
   {
      vector<Type> in(1, arg);
      vector<Type> out(1);

      posintpower(exponent, in, out);

      resultant = out[0];
      return;
   }

   if( exponent < -1 )
   {
      vector<Type> in(1, arg);
      vector<Type> out(1);

      posintpower(-exponent, in, out);

      resultant = Type(1.0)/out[0];
      return;
   }

   if( exponent == 1 )
   {
      resultant = arg;
      return;
   }

   if( exponent == 0 )
   {
      resultant = Type(1.0);
      return;
   }

   assert(exponent == -1);
   resultant = Type(1.0)/arg;
}

/** CppAD compatible evaluation of an expression for given arguments and parameters */
template<class Type>
SCIP_RETCODE eval(
   SCIP_EXPR*            expr,               /**< expression */
   const vector<Type>&   x,                  /**< values of variables */
   SCIP_Real*            param,              /**< values of parameters */
   Type&                 val                 /**< buffer to store expression value */
   )
{
   Type* buf;

   assert(expr != NULL);

   /* todo use SCIP_MAXCHILD_ESTIMATE as in expression.c */

   buf = NULL;
   if( SCIPexprGetNChildren(expr) )
   {
      if( BMSallocMemoryArray(&buf, SCIPexprGetNChildren(expr)) == NULL )
         return SCIP_NOMEMORY;

      for( int i = 0; i < SCIPexprGetNChildren(expr); ++i )
         SCIP_CALL( eval(SCIPexprGetChildren(expr)[i], x, param, buf[i]) );
   }

   switch(SCIPexprGetOperator(expr))
   {
   case SCIP_EXPR_VARIDX:
      assert(SCIPexprGetOpIndex(expr) < (int)x.size());
      val = x[SCIPexprGetOpIndex(expr)];
      break;

   case SCIP_EXPR_CONST:
      val = SCIPexprGetOpReal(expr);
      break;

   case SCIP_EXPR_PARAM:
      assert(param != NULL);
      val = param[SCIPexprGetOpIndex(expr)];
      break;

   case SCIP_EXPR_PLUS:
      val = buf[0] + buf[1];
      break;

   case SCIP_EXPR_MINUS:
      val = buf[0] - buf[1];
      break;

   case SCIP_EXPR_MUL:
      val = buf[0] * buf[1];
      break;

   case SCIP_EXPR_DIV:
      val = buf[0] / buf[1];
      break;

   case SCIP_EXPR_SQUARE:
      evalIntPower(val, buf[0], 2);
      break;

   case SCIP_EXPR_SQRT:
      evalSqrt(val, buf[0]);
      break;

   case SCIP_EXPR_REALPOWER:
      val = pow(buf[0], SCIPexprGetRealPowerExponent(expr));
      break;

   case SCIP_EXPR_INTPOWER:
      evalIntPower(val, buf[0], SCIPexprGetIntPowerExponent(expr));
      break;

   case SCIP_EXPR_SIGNPOWER:
      evalSignPower(val, buf[0], expr);
      break;

   case SCIP_EXPR_EXP:
      val = exp(buf[0]);
      break;

   case SCIP_EXPR_LOG:
      val = log(buf[0]);
      break;

   case SCIP_EXPR_SIN:
      val = sin(buf[0]);
      break;

   case SCIP_EXPR_COS:
      val = cos(buf[0]);
      break;

   case SCIP_EXPR_TAN:
      val = tan(buf[0]);
      break;
#if 0 /* these operators are currently disabled */
   case SCIP_EXPR_ERF:
      val = erf(buf[0]);
      break;

   case SCIP_EXPR_ERFI:
      return SCIP_ERROR;
#endif
   case SCIP_EXPR_MIN:
      evalMin(val, buf[0], buf[1]);
      break;

   case SCIP_EXPR_MAX:
      evalMax(val, buf[0], buf[1]);
      break;

   case SCIP_EXPR_ABS:
      evalAbs(val, buf[0]);
      break;

   case SCIP_EXPR_SIGN:
      val = sign(buf[0]);
      break;

   case SCIP_EXPR_SUM:
      val = 0.0;
      for (int i = 0; i < SCIPexprGetNChildren(expr); ++i)
         val += buf[i];
      break;

   case SCIP_EXPR_PRODUCT:
      val = 1.0;
      for (int i = 0; i < SCIPexprGetNChildren(expr); ++i)
         val *= buf[i];
      break;

   case SCIP_EXPR_LINEAR:
   {
      SCIP_Real* coefs;

      coefs = SCIPexprGetLinearCoefs(expr);
      assert(coefs != NULL || SCIPexprGetNChildren(expr) == 0);

      val = SCIPexprGetLinearConstant(expr);
      for (int i = 0; i < SCIPexprGetNChildren(expr); ++i)
         val += coefs[i] * buf[i];
      break;
   }

   case SCIP_EXPR_QUADRATIC:
   {
      SCIP_Real* lincoefs;
      SCIP_QUADELEM* quadelems;
      int nquadelems;
      SCIP_Real sqrcoef;
      Type lincoef;
      vector<Type> in(1);
      vector<Type> out(1);

      lincoefs   = SCIPexprGetQuadLinearCoefs(expr);
      nquadelems = SCIPexprGetNQuadElements(expr);
      quadelems  = SCIPexprGetQuadElements(expr);
      assert(quadelems != NULL || nquadelems == 0);

      SCIPexprSortQuadElems(expr);

      val = SCIPexprGetQuadConstant(expr);

      /* for each argument, we collect it's linear index from lincoefs, it's square coefficients and all factors from bilinear terms
       * then we compute the interval sqrcoef*x^2 + lincoef*x and add it to result */
      int i = 0;
      for( int argidx = 0; argidx < SCIPexprGetNChildren(expr); ++argidx )
      {
         if( i == nquadelems || quadelems[i].idx1 > argidx )
         {
            /* there are no quadratic terms with argidx in its first argument, that should be easy to handle */
            if( lincoefs != NULL )
               val += lincoefs[argidx] * buf[argidx];
            continue;
         }

         sqrcoef = 0.0;
         lincoef = lincoefs != NULL ? lincoefs[argidx] : 0.0;

         assert(i < nquadelems && quadelems[i].idx1 == argidx);
         do
         {
            if( quadelems[i].idx2 == argidx )
               sqrcoef += quadelems[i].coef;
            else
               lincoef += quadelems[i].coef * buf[quadelems[i].idx2];
            ++i;
         } while( i < nquadelems && quadelems[i].idx1 == argidx );
         assert(i == nquadelems || quadelems[i].idx1 > argidx);

         /* this is not as good as what we can get from SCIPintervalQuad, but easy to implement */
         if( sqrcoef != 0.0 )
         {
            in[0] = buf[argidx];
            posintpower(2, in, out);
            val += sqrcoef * out[0];
         }

         val += lincoef * buf[argidx];
      }
      assert(i == nquadelems);

      break;
   }

   case SCIP_EXPR_POLYNOMIAL:
   {
      SCIP_EXPRDATA_MONOMIAL** monomials;
      Type childval;
      Type monomialval;
      SCIP_Real exponent;
      int nmonomials;
      int nfactors;
      int* childidxs;
      SCIP_Real* exponents;
      int i;
      int j;

      val = SCIPexprGetPolynomialConstant(expr);

      nmonomials = SCIPexprGetNMonomials(expr);
      monomials  = SCIPexprGetMonomials(expr);

      for( i = 0; i < nmonomials; ++i )
      {
         nfactors  = SCIPexprGetMonomialNFactors(monomials[i]);
         childidxs = SCIPexprGetMonomialChildIndices(monomials[i]);
         exponents = SCIPexprGetMonomialExponents(monomials[i]);
         monomialval  = SCIPexprGetMonomialCoef(monomials[i]);

         for( j = 0; j < nfactors; ++j )
         {
            assert(childidxs[j] >= 0);
            assert(childidxs[j] <  SCIPexprGetNChildren(expr));

            childval = buf[childidxs[j]];
            exponent = exponents[j];

            /* cover some special exponents separately to avoid calling expensive pow function */
            if( exponent == 0.0 )
               continue;
            if( exponent == 1.0 )
            {
               monomialval *= childval;
               continue;
            }
            if( (int)exponent == exponent )
            {
               Type tmp;
               evalIntPower(tmp, childval, (int)exponent);
               monomialval *= tmp;
               continue;
            }
            if( exponent == 0.5 )
            {
               Type tmp;
               evalSqrt(tmp, childval);
               monomialval *= tmp;
               continue;
            }
            monomialval *= pow(childval, exponent);
         }

         val += monomialval;
      }

      break;
   }

   default:
      return SCIP_ERROR;
   }

   BMSfreeMemoryArrayNull(&buf);

   return SCIP_OKAY;
}

/** analysis an expression tree whether it requires retaping on every evaluation
 * this may be the case if the evaluation sequence depends on values of operands (e.g., in case of abs, sign, signpower, ...)
 */
bool needAlwaysRetape(SCIP_EXPR* expr)
{
   assert(expr != NULL);
   assert(SCIPexprGetChildren(expr) != NULL || SCIPexprGetNChildren(expr) == 0);

   for( int i = 0; i < SCIPexprGetNChildren(expr); ++i )
   {
      if( needAlwaysRetape(SCIPexprGetChildren(expr)[i]) )
         return true;
   }

   switch( SCIPexprGetOperator(expr) )
   {
   case SCIP_EXPR_MIN:
   case SCIP_EXPR_MAX:
   case SCIP_EXPR_ABS:
#ifndef CPPAD_USER_ATOMIC
   case SCIP_EXPR_SIGNPOWER:
#endif
      return true;

   default: ;
   }

   return false;
}

/** replacement for CppAD's default error handler
 * in debug mode, CppAD gives an error when an evaluation contains a nan
 * we do not want to stop execution in such a case, since the calling routine should check for nan's and decide what to do
 * since we cannot ignore this particular error, we ignore all
 * @todo find a way to check whether the error corresponds to a nan and communicate this back
 */
static
void cppaderrorcallback(
   bool               known,                 /**< is the error from a known source? */
   int                line,                  /**< line where error occured */
   const char*        file,                  /**< file where error occured */
   const char*        exp,                   /**< error condition */
   const char*        msg                    /**< error message */
   )
{
   SCIPdebugMessage("ignore CppAD error from %sknown source %s:%d: msg: %s exp: %s\n", known ? "" : "un", file, line, msg, exp);
}

/* install our error handler */
static CppAD::ErrorHandler errorhandler(cppaderrorcallback);

/** gets name and version of expression interpreter */
const char* SCIPexprintGetName(void)
{
   return CPPAD_PACKAGE_STRING;
}

/** gets descriptive text of expression interpreter */
const char* SCIPexprintGetDesc(void)
{
   return "Algorithmic Differentiation of C++ algorithms developed by B. Bell (www.coin-or.org/CppAD)";
}

/** gets capabilities of expression interpreter (using bitflags) */
SCIP_EXPRINTCAPABILITY SCIPexprintGetCapability(
   void
   )
{
   return SCIP_EXPRINTCAPABILITY_FUNCVALUE | SCIP_EXPRINTCAPABILITY_INTFUNCVALUE |
      SCIP_EXPRINTCAPABILITY_GRADIENT | SCIP_EXPRINTCAPABILITY_INTGRADIENT |
      SCIP_EXPRINTCAPABILITY_HESSIAN;
}

/** creates an expression interpreter object */
SCIP_RETCODE SCIPexprintCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPRINT**        exprint             /**< buffer to store pointer to expression interpreter */
   )
{
   assert(blkmem  != NULL);
   assert(exprint != NULL);

   if( BMSallocMemory(exprint) == NULL )
      return SCIP_NOMEMORY;

   (*exprint)->blkmem = blkmem;

   return SCIP_OKAY;
}

/** frees an expression interpreter object */
SCIP_RETCODE SCIPexprintFree(
   SCIP_EXPRINT**        exprint             /**< expression interpreter that should be freed */
   )
{
   assert( exprint != NULL);
   assert(*exprint != NULL);

   BMSfreeMemory(exprint);

   return SCIP_OKAY;
}

/** compiles an expression tree and stores compiled data in expression tree */
SCIP_RETCODE SCIPexprintCompile(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{
   assert(tree    != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   if (!data)
   {
      data = new SCIP_EXPRINTDATA();
      assert( data != NULL );
      SCIPexprtreeSetInterpreterData(tree, data);
      SCIPdebugMessage("set interpreter data in tree %p to %p\n", (void*)tree, (void*)data);
   }
   else
   {
      data->need_retape     = true;
      data->int_need_retape = true;
   }

   int n = SCIPexprtreeGetNVars(tree);

   data->X.resize(n);
   data->x.resize(n);
   data->Y.resize(1);

   data->int_X.resize(n);
   data->int_x.resize(n);
   data->int_Y.resize(1);

   if( data->root != NULL )
   {
      SCIPexprFreeDeep(exprint->blkmem, &data->root);
   }

   SCIP_EXPR* root = SCIPexprtreeGetRoot(tree);

   SCIP_CALL( SCIPexprCopyDeep(exprint->blkmem, &data->root, root) );

   data->need_retape_always = needAlwaysRetape(SCIPexprtreeGetRoot(tree));

   data->blkmem = exprint->blkmem;

   return SCIP_OKAY;
}

/** frees interpreter data */
SCIP_RETCODE SCIPexprintFreeData(
   SCIP_EXPRINTDATA**    interpreterdata     /**< interpreter data that should freed */
   )
{
   assert( interpreterdata != NULL);
   assert(*interpreterdata != NULL);

   if( (*interpreterdata)->root != NULL )
      SCIPexprFreeDeep((*interpreterdata)->blkmem, &(*interpreterdata)->root);   

   delete *interpreterdata;
   *interpreterdata = NULL; 

   return SCIP_OKAY;
}

/** notify expression interpreter that a new parameterization is used
 * this probably causes retaping by AD algorithms
 */
SCIP_RETCODE SCIPexprintNewParametrization(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{
   assert(exprint != NULL);
   assert(tree    != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   if( data != NULL )
   {
      data->need_retape     = true;
      data->int_need_retape = true;
   }

   return SCIP_OKAY;
}

/** evaluates an expression tree */
SCIP_RETCODE SCIPexprintEval(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables */
   SCIP_Real*            val                 /**< buffer to store value */
   )
{
   SCIP_EXPRINTDATA* data;

   assert(exprint != NULL);
   assert(tree    != NULL);
   assert(varvals != NULL);
   assert(val     != NULL);

   data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);
   assert(SCIPexprtreeGetNVars(tree) == (int)data->X.size());
   assert(SCIPexprtreeGetRoot(tree)  != NULL);

   int n = SCIPexprtreeGetNVars(tree);

   if( data->need_retape_always || data->need_retape )
   {
      for( int i = 0; i < n; ++i )
      {
         data->X[i] = varvals[i];
         data->x[i] = varvals[i];
      }

      CppAD::Independent(data->X);

      if( data->root != NULL )
         SCIP_CALL( eval(data->root, data->X, SCIPexprtreeGetParamVals(tree), data->Y[0]) );
      else
         data->Y[0] = 0.0;

      data->f.Dependent(data->X, data->Y);

      data->val = Value(data->Y[0]);
      SCIPdebugMessage("Eval retaped and computed value %g\n", data->val);

      // the following is required if the gradient shall be computed by a reverse sweep later
      // data->val = data->f.Forward(0, data->x)[0];

      data->need_retape = false;
   }
   else
   {
      assert((int)data->x.size() >= n);
      for( int i = 0; i < n; ++i )
         data->x[i] = varvals[i];

      data->val = data->f.Forward(0, data->x)[0];
      SCIPdebugMessage("Eval used foward sweep to compute value %g\n", data->val);
   }

   *val = data->val;

   return SCIP_OKAY;
}

/** evaluates an expression tree on intervals */
extern
SCIP_RETCODE SCIPexprintEvalInt(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        varvals,            /**< interval values of variables */
   SCIP_INTERVAL*        val                 /**< buffer to store interval value of expression */
   )
{
   SCIP_EXPRINTDATA* data;

   assert(exprint != NULL);
   assert(tree    != NULL);
   assert(varvals != NULL);
   assert(val     != NULL);

   data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);
   assert(SCIPexprtreeGetNVars(tree) == (int)data->int_X.size());
   assert(SCIPexprtreeGetRoot(tree)  != NULL);

   int n = SCIPexprtreeGetNVars(tree);

   SCIPInterval::infinity = infinity;

   if( data->int_need_retape || data->need_retape_always )
   {
      for( int i = 0; i < n; ++i )
      {
         data->int_X[i] = varvals[i];
         data->int_x[i] = varvals[i];
      }

      CppAD::Independent(data->int_X);

      if( data->root != NULL )
         SCIP_CALL( eval(data->root, data->int_X, SCIPexprtreeGetParamVals(tree), data->int_Y[0]) );
      else
         data->int_Y[0] = 0.0;

      data->int_f.Dependent(data->int_X, data->int_Y);

      data->int_val = Value(data->int_Y[0]);

      data->int_need_retape = false;
   }
   else
   {
      assert((int)data->int_x.size() >= n);
      for( int i = 0; i < n; ++i )
         data->int_x[i] = varvals[i];

      data->int_val = data->int_f.Forward(0, data->int_x)[0];
   }

   *val = data->int_val;

   return SCIP_OKAY;
}

/** computes value and gradient of an expression tree */
SCIP_RETCODE SCIPexprintGrad(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to a point evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store expression value */
   SCIP_Real*            gradient            /**< buffer to store expression gradient, need to have length at least SCIPexprtreeGetNVars(tree) */
   )
{
   assert(exprint  != NULL);
   assert(tree     != NULL);
   assert(varvals  != NULL || new_varvals == FALSE);
   assert(val      != NULL);
   assert(gradient != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);

   if( new_varvals )
   {
      SCIP_CALL( SCIPexprintEval(exprint, tree, varvals, val) );
   }
   else
      *val = data->val;

   int n = SCIPexprtreeGetNVars(tree);

   vector<double> jac(data->f.Jacobian(data->x));

   for( int i = 0; i < n; ++i )
      gradient[i] = jac[i];

#ifdef SCIP_DEBUG
   SCIPdebugMessage("Grad for "); SCIPexprtreePrint(tree, NULL, NULL, NULL); printf("\n");
   SCIPdebugMessage("x    ="); for (int i = 0; i < n; ++i) printf("\t %g", data->x[i]); printf("\n");
   SCIPdebugMessage("grad ="); for (int i = 0; i < n; ++i) printf("\t %g", gradient[i]); printf("\n");
#endif

   return SCIP_OKAY;
}

/** computes interval value and interval gradient of an expression tree */
SCIP_RETCODE SCIPexprintGradInt(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        varvals,            /**< interval values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable interval values changed since last call to an interval evaluation routine? */
   SCIP_INTERVAL*        val,                /**< buffer to store expression interval value */
   SCIP_INTERVAL*        gradient            /**< buffer to store expression interval gradient, need to have length at least SCIPexprtreeGetNVars(tree) */
   )
{
   assert(exprint  != NULL);
   assert(tree     != NULL);
   assert(varvals  != NULL || new_varvals == false);
   assert(val      != NULL);
   assert(gradient != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);

   if (new_varvals)
      SCIP_CALL( SCIPexprintEvalInt(exprint, tree, infinity, varvals, val) );
   else
      *val = data->int_val;

   int n = SCIPexprtreeGetNVars(tree);

   vector<SCIPInterval> jac(data->int_f.Jacobian(data->int_x));

   for (int i = 0; i < n; ++i)
      gradient[i] = jac[i];

#ifdef SCIP_DEBUG
   SCIPdebugMessage("GradInt for "); SCIPexprtreePrint(tree, NULL, NULL, NULL); printf("\n");
   SCIPdebugMessage("x    ="); for (int i = 0; i < n; ++i) printf("\t [%g,%g]", SCIPintervalGetInf(data->int_x[i]), SCIPintervalGetSup(data->int_x[i])); printf("\n");
   SCIPdebugMessage("grad ="); for (int i = 0; i < n; ++i) printf("\t [%g,%g]", SCIPintervalGetInf(gradient[i]), SCIPintervalGetSup(gradient[i])); printf("\n");
#endif

   return SCIP_OKAY;
}

/** gives sparsity pattern of hessian
 * NOTE: this function might be replaced later by something nicer 
 * Since the AD code might need to do a forward sweep, you should pass variable values in here.
 */
SCIP_RETCODE SCIPexprintHessianSparsityDense(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables */
   SCIP_Bool*            sparsity            /**< buffer to store sparsity pattern of Hessian, sparsity[i+n*j] indicates whether entry (i,j) is nonzero in the hessian */
   )
{
   assert(exprint  != NULL);
   assert(tree     != NULL);
   assert(varvals  != NULL);
   assert(sparsity != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);

   int n = SCIPexprtreeGetNVars(tree);
   int nn = n*n;

   if( data->need_retape_always )
   {
      // @todo can we do something better here, e.g., by looking at the expression tree by ourself?

      for( int i = 0; i < nn; ++i )
         sparsity[i] = TRUE;

#ifdef SCIP_DEBUG
      SCIPdebugMessage("HessianSparsityDense for "); SCIPexprtreePrint(tree, NULL, NULL, NULL); printf("\n");
      SCIPdebugMessage("sparsity = all elements, due to discontinuouities\n");
#endif

      return SCIP_OKAY;
   }

   if( data->need_retape )
   {
      SCIP_Real val;
      SCIP_CALL( SCIPexprintEval(exprint, tree, varvals, &val) );
   }

   SCIPdebugMessage("calling ForSparseJac\n");

   vector<bool> r(nn, false);
   for (int i = 0; i < n; ++i)
      r[i*n+i] = true;
   data->f.ForSparseJac(n, r); // need to compute sparsity for Jacobian first

   SCIPdebugMessage("calling RevSparseHes\n");

   vector<bool> s(1, true);
   vector<bool> sparsehes(data->f.RevSparseHes(n, s));

   for( int i = 0; i < nn; ++i )
      sparsity[i] = sparsehes[i];

#ifdef SCIP_DEBUG
   SCIPdebugMessage("HessianSparsityDense for "); SCIPexprtreePrint(tree, NULL, NULL, NULL); printf("\n");
   SCIPdebugMessage("sparsity ="); for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if (sparsity[i*n+j]) printf(" (%d,%d)", i, j); printf("\n");
#endif

   return SCIP_OKAY;
}

/** computes value and dense hessian of an expression tree
 * the full hessian is computed (lower left and upper right triangle)
 */
SCIP_RETCODE SCIPexprintHessianDense(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to an evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store function value */
   SCIP_Real*            hessian             /**< buffer to store hessian values, need to have size at least n*n */
   )
{
   assert(exprint != NULL);
   assert(tree    != NULL);
   assert(varvals != NULL || new_varvals == FALSE);
   assert(val     != NULL);
   assert(hessian != NULL);

   SCIP_EXPRINTDATA* data = SCIPexprtreeGetInterpreterData(tree);
   assert(data != NULL);

   if( new_varvals )
   {
      SCIP_CALL( SCIPexprintEval(exprint, tree, varvals, val) );
   }
   else
      *val = data->val;

   int n = SCIPexprtreeGetNVars(tree);

   vector<double> hess(data->f.Hessian(data->x, 0));

   int nn = n*n;
   for (int i = 0; i < nn; ++i)
      hessian[i] = hess[i];

#ifdef SCIP_DEBUG
   SCIPdebugMessage("HessianDense for "); SCIPexprtreePrint(tree, NULL, NULL, NULL); printf("\n");
   SCIPdebugMessage("x    ="); for (int i = 0; i < n; ++i) printf("\t %g", data->x[i]); printf("\n");
   SCIPdebugMessage("hess ="); for (int i = 0; i < n*n; ++i) printf("\t %g", hessian[i]); printf("\n");
#endif

   return SCIP_OKAY;
}
