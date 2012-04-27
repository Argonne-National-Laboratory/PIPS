/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/eigval.cxx

    This file is part of ConciBundle, a C/C++ library for convex optimization.

    ConicBundle is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ConicBundle is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

***************************************************************************** */



/* The routines tred2 and imtql2 have been translated by f2c 
   from Eispack and have been adapted to the package so 
   that no additional libraries besides -lm have to be linked.
*/

#include <math.h>
#include "mymath.hxx"
#include "symmat.hxx"

 namespace CH_Matrix_Classes {

Integer Symmatrix::eig(Matrix& A,Matrix& d,bool sort_non_decreasingly) const
{
  chk_init(*this);
  if (sort_non_decreasingly){
    A.init(*this);
  }
  else {
    A.init(*this,-1.);
  }
  d.newsize(nr,1);
  if (nr==0)
    return 0;
  Matrix e(nr,1);
  if (tred2(nr,nr,A.get_store(),d.get_store(),e.get_store(),A.get_store())){
    MEmessage(MatrixError(ME_unspec,"Symmatrix::eig(Matrix&,Matrix&,bool) tred2_ failed",MTsymmetric));
  }
  Integer ret_val=imtql2(nr,nr,d.get_store(),e.get_store(),A.get_store());
  chk_set_init(d,1);
  if (ret_val!=0){
    MEmessage(MatrixError(ME_warning,"Symmatrix::eig(Matrix&,Matrix&,bool) imtql2_ failed",MTsymmetric));
    Indexmatrix sind;
    sortindex(d,sind);
    d=d(sind);
    A=A.cols(sind);
  }
  if (!sort_non_decreasingly) d*=-1;
  return ret_val;
}


/* Subroutine */ Integer Symmatrix::tred2(Integer nm, Integer n, Real *a, 
	Real *d, Real *e, Real *z) const
{
    /* System generated locals */
    Integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;

    /* Builtin functions */
//    Real sqrt(Real), d_sign(Real, Real);

    /* Local variables */
    static Real f, g, h;
    static Integer i, j, k, l;
    static Real scale, hh;
    static Integer ii;



/*     this subroutine is a translation of the algol procedure tred2, */
/*     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson. */
/*     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971). */

/*     this subroutine reduces a real symmetric matrix to a */
/*     symmetric tridiagonal matrix using and accumulating */
/*     orthogonal similarity transformations. */

/*     on input */

/*        nm must be set to the row dimension of two-dimensional */
/*          array parameters as declared in the calling program */
/*          dimension statement. */

/*        n is the order of the matrix. */

/*        a contains the real symmetric input matrix.  only the */
/*          lower triangle of the matrix need be supplied. */

/*     on output */

/*        d contains the diagonal elements of the tridiagonal matrix. */

/*        e contains the subdiagonal elements of the tridiagonal */
/*          matrix in its last n-1 positions.  e(1) is set to zero. */

/*        z contains the orthogonal transformation matrix */
/*          produced in the reduction. */

/*        a and z may coincide.  if distinct, a is unaltered. */

/*     questions and comments should be directed to burton s. garbow, */
/*     mathematics and computer science div, argonne national laboratory 
*/

/*     this version dated august 1983. */

/*     ------------------------------------------------------------------ 
*/

    /* Parameter adjustments */
    z_dim1 = nm;
    z_offset = z_dim1 + 1;
    z -= z_offset;
    --e;
    --d;
    a_dim1 = nm;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    for (i = 1; i <= n; ++i) {
	for (j = i; j <= n; ++j) {
	    z[j + i * z_dim1] = a[j + i * a_dim1];
	}
	d[i] = a[n + i * a_dim1];
    }

    if (n == 1) {
	goto L510;
    }

/*     .......... for i=n step -1 until 2 do -- .......... */
    i__1 = n;
    for (ii = 2; ii <= i__1; ++ii) {
	i = n + 2 - ii;
	l = i - 1;
	h = 0.;
	scale = 0.;
	if (l < 2) {
	    goto L130;
	}
/*     .......... scale row (algol tol then not needed) .......... */
	for (k = 1; k <= l; ++k) {
	    scale += abs( d[k]);
	}

	if (scale != 0.) {
	    goto L140;
	}
L130:
	e[i] = d[l];

	for (j = 1; j <= l; ++j) {
	    d[j] = z[l + j * z_dim1];
	    z[i + j * z_dim1] = 0.;
	    z[j + i * z_dim1] = 0.;
	}

	goto L290;

L140:
	for (k = 1; k <= l; ++k) {
	    d[k] /= scale;
	    h += d[k] * d[k];
	}

	f = d[l];
	g = -d_sign(::sqrt(h),f);
	e[i] = scale * g;
	h -= f * g;
	d[l] = f - g;

/*     .......... form a*u .......... */
	for (j = 1; j <= l; ++j) {
	    e[j] = 0.;
	}

	for (j = 1; j <= l; ++j) {
	    f = d[j];
	    z[j + i * z_dim1] = f;
	    g = e[j] + z[j + j * z_dim1] * f;

	    for (k = j+1; k <= l; ++k) {
		g += z[k + j * z_dim1] * d[k];
		e[k] += z[k + j * z_dim1] * f;
	    }

	    e[j] = g;
	}

/*     .......... form p .......... */
	f = 0.;

	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    e[j] /= h;
	    f += e[j] * d[j];
	}

	hh = f / (h + h);

/*     .......... form q .......... */
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    e[j] -= hh * d[j];
	}

/*     .......... form reduced a .......... */
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    f = d[j];
	    g = e[j];
            Real *zpj=z+j+j * z_dim1;
            Real *ep=e+j;
            Real *dp=d+j;
	    for (k = j; k <= l; ++k) {
		//z[k + j * z_dim1] -= f * e[k] + g * d[k];
                (*zpj++)-= f*(*ep++)+g*(*dp++);
	    }

	    d[j] = z[l + j * z_dim1];
	    z[i + j * z_dim1] = 0.;
	}

L290:
	d[i] = h;
    }

/*     .......... accumulation of transformation matrices .......... */
    i__1 = n;
    for (i = 2; i <= i__1; ++i) {
	l = i - 1;
	z[n + l * z_dim1] = z[l + l * z_dim1];
	z[l + l * z_dim1] = 1.;
	h = d[i];
	if (h == 0.) {
	    goto L380;
	}

	//for (k = 1; k <= l; ++k) {
	//    d[k] = z[k + i * z_dim1] / h;
	//}
        mat_xeya(l,d+1,z+1+i*z_dim1,1./h);

	for (j = 1; j <= l; ++j) {
	    //g = 0.;
	    //for (k = 1; k <= l; ++k) {
	    //    g += z[k + i * z_dim1] * z[k + j * z_dim1];
	    //}
            g=mat_ip(l,z+1+i * z_dim1,z+1+j * z_dim1);

	    //for (k = 1; k <= l; ++k) {
	    //    z[k + j * z_dim1] -= g * d[k];
	    //}
            mat_xpeya(l,z+1+j*z_dim1,d+1,-g);
	}

L380:
	//for (k = 1; k <= l; ++k) {
	//    z[k + i * z_dim1] = 0.;
	//}
        mat_xea(l,z+1+i * z_dim1,0.);

    }

L510:
    for (i = 1; i <= n; ++i) {
	d[i] = z[n + i * z_dim1];
	z[n + i * z_dim1] = 0.;
    }

    z[n + n * z_dim1] = 1.;
    e[1] = 0.;
    return 0;
} /* tred2_ */

/* Subroutine */ Integer Symmatrix::imtql2(Integer nm, Integer n, Real *d, 
	Real *e, Real *z) const
{
    /* System generated locals */
    Integer z_dim1, z_offset, i__1, i__2, i__3;

    /* Builtin functions */
//    Real sqrt(Real), d_sign(Real *, Real *);

    /* Local variables */
    static Real b, c, f, g;
    static Integer i, j, k, l, mi;
    static Real p, r, s;
    static Integer ii, mml;
    static Real tst1, tst2;
    static Integer ierr;



/*     this subroutine is a translation of the algol procedure imtql2, */
/*     num. math. 12, 377-383(1968) by martin and wilkinson, */
/*     as modified in num. math. 15, 450(1970) by dubrulle. */
/*     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971). */

/*     this subroutine finds the eigenvalues and eigenvectors */
/*     of a symmetric tridiagonal matrix by the implicit ql method. */
/*     the eigenvectors of a full symmetric matrix can also */
/*     be found if  tred2  has been used to reduce this */
/*     full matrix to tridiagonal form. */

/*     on input */

/*        nm must be set to the row dimension of two-dimensional */
/*          array parameters as declared in the calling program */
/*          dimension statement. */

/*        n is the order of the matrix. */

/*        d contains the diagonal elements of the input matrix. */

/*        e contains the subdiagonal elements of the input matrix */
/*          in its last n-1 positions.  e(1) is arbitrary. */

/*        z contains the transformation matrix produced in the */
/*          reduction by  tred2, if performed.  if the eigenvectors */
/*          of the tridiagonal matrix are desired, z must contain */
/*          the identity matrix. */

/*      on output */

/*        d contains the eigenvalues in ascending order.  if an */
/*          error exit is made, the eigenvalues are correct but */
/*          unordered for indices 1,2,...,ierr-1. */

/*        e has been destroyed. */

/*        z contains orthonormal eigenvectors of the symmetric */
/*          tridiagonal (or full) matrix.  if an error exit is made, */
/*          z contains the eigenvectors associated with the stored */
/*          eigenvalues. */

/*        ierr is set to */ /*C. Helmberg: ierr is now return value*/
/*          zero       for normal return, */
/*          j          if the j-th eigenvalue has not been */
/*                     determined after 30 iterations. */

/*     calls pythag for  dsqrt(a*a + b*b) . */

/*     questions and comments should be directed to burton s. garbow, */
/*     mathematics and computer science div, argonne national laboratory 
*/

/*     this version dated august 1983. */

/*     ------------------------------------------------------------------ 
*/

    /* Parameter adjustments */
    z_dim1 = nm;
    z_offset = z_dim1 + 1;
    z -= z_offset;
    --e;
    --d;

    /* Function Body */
    ierr = 0;
    if (n == 1) {
	goto L1001;
    }

    i__1 = n;
    for (i = 2; i <= i__1; ++i) {
/* L100: */
	e[i - 1] = e[i];
    }

    e[n] = 0.;

    i__1 = n;
    for (l = 1; l <= i__1; ++l) {
	j = 0;
/*     .......... look for small sub-diagonal element .......... */
L105:
	i__2 = n;
	for (mi = l; mi <= i__2; ++mi) {
	    if (mi == n) {
		goto L120;
	    }
	    tst1 = abs(d[mi])+abs(d[mi+1]);
	    tst2 = tst1 + abs(e[mi]);
	    if (tst2 == tst1) {
		goto L120;
	    }
/* L110: */
	}

L120:
	p = d[l];
	if (mi == l) {
	    goto L240;
	}
	if (j == 30) {
	    goto L1000;
	}
	++j;
/*     .......... form shift .......... */
	g = (d[l + 1] - p) / (e[l] * 2.);
	r = ::sqrt(g * g + 1.);
	g = d[mi] - p + e[l] / (g + d_sign(r, g));
	s = 1.;
	c = 1.;
	p = 0.;
	mml = mi - l;
/*     .......... for i=mi-1 step -1 until l do -- .......... */
	i__2 = mml;
	for (ii = 1; ii <= i__2; ++ii) {
	    i = mi - ii;
	    f = s * e[i];
	    b = c * e[i];
	    r = ::sqrt(f * f + g * g);
	    e[i + 1] = r;
	    if (r == 0.) {
		goto L210;
	    }
	    s = f / r;
	    c = g / r;
	    g = d[i + 1] - p;
	    r = (d[i] - g) * s + c * 2. * b;
	    p = s * r;
	    d[i + 1] = g + p;
	    g = c * r - b;
/*     .......... form vector .......... */
	    i__3 = n;
            Real* zpi=z+i*z_dim1;
            Real* zpi1=zpi+z_dim1;
	    for (k = 1; k <= i__3; ++k) {
		//f = z[k + (i + 1) * z_dim1];
                f=*(++zpi1);
		//z[k + (i + 1) * z_dim1] = s * z[k + i * z_dim1] + c * f;
                *zpi1=s*(*(++zpi))+c*f;
		//z[k + i * z_dim1] = c * z[k + i * z_dim1] - s * f;
                *zpi*=c; *zpi-=s*f;
/* L180: */
	    }

/* L200: */
	}

	d[l] -= p;
	e[l] = g;
	e[mi] = 0.;
	goto L105;
/*     .......... recover from underflow .......... */
L210:
	d[i + 1] -= p;
	e[mi] = 0.;
	goto L105;
L240:
	;
    }
/*     .......... order eigenvalues and eigenvectors .......... */
    i__1 = n;
    for (ii = 2; ii <= i__1; ++ii) {
	i = ii - 1;
	k = i;
	p = d[i];

	i__2 = n;
	for (j = ii; j <= i__2; ++j) {
	    if (d[j] >= p) {
		goto L260;
	    }
	    k = j;
	    p = d[j];
L260:
	    ;
	}

	if (k == i) {
	    goto L300;
	}
	d[k] = d[i];
	d[i] = p;

	i__2 = n;
	//for (j = 1; j <= i__2; ++j) {
	    //p = z[j + i * z_dim1];
	    //z[j + i * z_dim1] = z[j + k * z_dim1];
	    //z[j + k * z_dim1] = p;
            
/* L280: */
	//}
        mat_swap(i__2,z+1+i*z_dim1,z+1+k*z_dim1);

L300:
	;
    }

    goto L1001;
/*     .......... set error -- no convergence to an */
/*                eigenvalue after 30 iterations .......... */
L1000:
    ierr = l;
L1001:
    return ierr;
} /* imtql2_ */

}

