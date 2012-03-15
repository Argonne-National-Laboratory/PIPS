/****************************************************************
Copyright (C) 1997 Lucent Technologies
All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the name of Lucent or any of its entities
not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior
permission.

LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
****************************************************************/

/* Sample mex function (in MATLAB 4.x format) for getting functions,
   gradients, and dense Hessians from an AMPL .nl file.  Start with

	[x,bl,bu,v,cl,cu] = amplfunc('stub')

   to read in a problem (discarding the previous problem, if any).
   The return values are:

	x = primal initial guess
	v = dual initial guess
	bl, bu = lower and upper bounds on x
	cl, cu = lower and upper bounds on c (the constraint bodies).

   Then

	[f,c] = amplfunc(x,0)

   gives the function (f) and constraint bodies (c) at x;

	[g,Jac] = amplfunc(x,1)

   gives the gradient g of f, the Jacobian matrix J of c, and

	W = amplfunc(v)

   gives the Hessian W of the Lagrangian function L = f + v*c
   (at the last x at which amplfunc(x,1) was called).

   After finding optimal values for x and v,

	amplfunc('solution message',x,v)

   to write a stub.sol file.
*/

#include "mex.h"
#undef printf
#include "asl_pfgh.h"

static char msgbuf[256];

 static real*
sizechk(Matrix *mp, char *who, fint m)
{
	int m1, n1;
	m1 = mxGetM(mp);
	n1 = mxGetN(mp);
	if (m1 != m || n1 != 1 && m1) {
		sprintf(msgbuf,
			"Expected %s to be %d x 1 rather than %d x %d\n",
			who, m, m1, n1);
		mexErrMsgTxt(msgbuf);
		}
	return mxGetPr(mp);
	}

 static void
at_end(void)
{
	if (cur_ASL)
		ASL_free(&cur_ASL);
	}

 static void
usage(void)
{
	mexErrMsgTxt("amplfunc usage:\n\n\
	[x,bl,bu,v,cl,cu] = amplfunc('stub')\nor\n\
	[f,c] = amplfunc(x,0)\nor\n\
	[g,Jac] = amplfunc(x,1)\nor\n\
	W = amplfunc(v)\nor\n\
	amplfunc('solution message',x,v)\nwith\n\
	x = primal, v = dual variables");
	}

 void
mexFunction(int nlhs, Matrix **plhs, int nrhs, Matrix **prhs)
{
	FILE *nl;
	char *buf1, buf[512], *what;
	static fint n, nc, nz;
	fint nerror;
	real *J1, *W, *c, *f, *g, *v, *t, *x;
	static real *J;
	cgrad *cg, **cgp, **cgpe;
	static size_t Jsize;
	Jmp_buf err_jmp0;
	ASL *asl = cur_ASL;

	if (nrhs == 1 && mxIsString(prhs[0])) {
		if (nlhs != 6)
			usage();
		if (mxGetString(prhs[0], buf1 = buf, sizeof(buf)))
			mexErrMsgTxt("Expected 'stub' as argument\n");
		at_end();
		mexAtExit(at_end);
		asl = ASL_alloc(ASL_read_pfgh);
		return_nofile = 1;
		if (!(nl = jac0dim(buf1,strlen(buf)))) {
			sprintf(msgbuf, "Can't open %.*s\n",
				sizeof(msgbuf)-20, buf);
			mexErrMsgTxt(msgbuf);
			}
		if (n_obj <= 0)
			printf("Warning: objectve == 0\n");
		n = n_var;
		nc = n_con;
		nz = nzc;
		J = (real *)M1alloc(nz*sizeof(real));
		X0 = mxGetPr(plhs[0] = mxCreateFull(n, 1, REAL));
		LUv = mxGetPr(plhs[1] = mxCreateFull(n, 1, REAL));
		Uvx = mxGetPr(plhs[2] = mxCreateFull(n, 1, REAL));
		pi0 = mxGetPr(plhs[3] = mxCreateFull(nc, 1, REAL));
		LUrhs = mxGetPr(plhs[4] = mxCreateFull(nc, 1, REAL));
		Urhsx = mxGetPr(plhs[5] = mxCreateFull(nc, 1, REAL));
		pfgh_read(nl, ASL_findgroups);
		Jsize = nc*n*sizeof(real);
		return;
		}

	if (!filename)
		mexErrMsgTxt("amplfunc(\"stub\") has not been called\n");
	nerror = -1;
	err_jmp1 = &err_jmp0;
	if (nlhs == 2) {
		if (nrhs != 2)
			usage();
		x = sizechk(prhs[0],"x",n);
		t = sizechk(prhs[1],"0 or 1", 1);
		if (t[0] == 0.) {
			f = mxGetPr(plhs[0] = mxCreateFull(1, 1, REAL));
			c = mxGetPr(plhs[1] = mxCreateFull(nc, 1, REAL));
			if (setjmp(err_jmp0.jb)) {
				sprintf(msgbuf, "Trouble evaluating %s\n",
					what);
				mexErrMsgTxt(msgbuf);
				}
			what = "f";
			*f = objval(0, x, &nerror);
			what = "c";
			conval(x, c, &nerror);
			return;
			}
		g = mxGetPr(plhs[0] = mxCreateFull(n, 1, REAL));
		J1 = mxGetPr(plhs[1] = mxCreateFull(nc, n, REAL));
		what = "g";
		objgrd(0, x, g, &nerror);
		if (nc) {
			memset(J1, 0, Jsize);
			what = "J";
			jacval(x, J, &nerror);
			cgp = Cgrad;
			for(cgpe = cgp + nc; cgp < cgpe; J1++)
				for(cg = *cgp++; cg; cg = cg->next)
					J1[nc*cg->varno] = J[cg->goff];
			}
		return;
		}
	if (nlhs == 0 && nrhs == 3) {
		/* eval2('solution message', x, v): x = primal, v = dual */
		if (!mxIsString(prhs[0]))
			usage();
		x = sizechk(prhs[1],"x",n);
		v = sizechk(prhs[2],"v",nc);
		if (mxGetString(prhs[0], buf, sizeof(buf)))
			mexErrMsgTxt(
			 "Expected 'solution message' as first argument\n");
		write_sol(buf, x, v, 0);
		return;
		}
	if (nlhs != 1 || nrhs != 1)
		usage();
	v = sizechk(prhs[0],"v",nc);
	W = mxGetPr(plhs[0] = mxCreateFull(n, n, REAL));

	what = "W";
	fullhes(W, n, 0, 0, v);
	}
