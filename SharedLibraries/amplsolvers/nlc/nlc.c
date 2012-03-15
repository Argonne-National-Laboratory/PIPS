/****************************************************************
Copyright (C) 1997, 2001 Lucent Technologies
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

#include "nlp.h"
#undef f_OPNUM
#include "opcode.hd"
#include "nlc.h"
#define asl ((ASL_fg*)cur_ASL)

#define Egulp 400

typedef struct maxinfo {
	int ncond;
	int ndv;
	int nvt;
	int needT1;
	/* !! plus goo for function evaluations */
	} maxinfo;

static maxinfo bmax, cmax, omax;

 static int cwant = 2, derkind = 2, djoff = -1, owant = 2;
 static int needx0check, nvseen, output_time, *vseen;
 static list **c_ifset, **com_ifset, **o_ifset;
 int branches;
 static int Fortran1, condlevel, fwalk, iflevel, maxa, ncond, ndv, ndvmax,
	npd, nplterm, nv1, nv2, nvt, nvtmax;
 static int *cvmap;
 static fint nJ;
 static expr *Plterms;
 static char *intsk;
 efunc *r_ops[N_OPS];
 extern efunc *r_op[];
#include "r_ops1.hd"

 static real *bounds;
 extern char *progname;
 char *Half = "0.5", *One = "1.", *Negone = "-1.", *Zero = "0.";
 char *opEQ = "==", *opGE = ">=", *opGT = ">", *opLE = "<=",
	*opLT = "<", *opNE = "!=";
 char *opAND = "&&", *opOR = "||";
 char *T = "t1", *T1 = "t2";	/* temporary variable names */
 char *offlfmt1 = "&%s", *offlfmt2 = "&%s, &%s";
 static char *assign_fmt = "\t%s = %s;\n";
 static char *binop_fmt = "\t%s = %s %s %s;\n";
 static char *goto_fmt = "\tgoto L%d;\n";
 static char *ifgo_fmt = "\tif (%s %s %s) goto L%d;\n";
 static char *label_fmt = " L%d:\n";
 static char *ifstart_fmt = "\tif (%s %s %s) {\n";
 static char *elseif_fmt = "\t} else if (%s %s %s) {\n";
 static char *else_fmt = "\t} else {\n";
 static char *endif_fmt = "\t}\n";
 static char *case_fmt = "   case %d:\n";
 char *cond_fmt = "cond[%d]";
 static char *break_fmt = "\tbreak;\n";
 static char *endswitch_fmt = "\t}\n";
 static char *zerdiv_fmt = "\tzerdiv_(&%s);";
 static char *call_fmt = "\t%s(&%s);\n";
 static char *dv_fmt = "dv[%d]";
 static char *Fortstar = "";
 char *pd_fmt = "pd[%d]";
 static char *tv_fmt = "v[%d]";
 static char *x_fmt = "x[%d]";
 static char seen[N_OPS], *declare[N_OPS];
 static char *grad_fmt = "g[%d]", *jac_fmt = "J[%d]";
 static char *g0_fmt = "\tg[%d] = %s;\n", *j0_fmt = "\tJ[%d] = %s;\n";
 static char *eos = ";\n";
 static char *star = "";
 static char *Void = "void";
 static char *xcheck = "\tfint wantfg = *needfg;\n\
	if (xcheck(x) && wantfg == 2)\n\t\twantfg = 3;\n";
 static char *xcheck0 = "\txcheck(x);\n";
 static char *xcheckdcl = "(real *x)";
 static char *xkind = "\tif (!(xkind & %d)) {\n\t\txkind |= %d;\n";
 static real negone = -1.0, one = 1.0;
 extern real edagread_one;	/* &edagread_one = special derp->c value */

 static expr_nx *nums;
 static ograd **stog;

 static Adjoint *dwalk_one;
 static real *rdwalk_one;
 static Adjoint A0;
 static int *acount, *dvtfree, *rownos;
 static int densejac, krflag, ndvfree, ndvtmax, needT1, nvtfree,
	    pdts, stor_grad;

 typedef struct
pdl {
	struct pdl *next;
	int i, ts;
	} pdl;

 static pdl *pdlbusy, *pdlfree;

 static int
new_pd(VOID)
{
	pdl *p;
	if (!iflevel)
		return npd++;
	if (p = pdlfree)
		pdlfree = p->next;
	else {
		p = (pdl *)mem(sizeof(pdl));
		p->i = npd++;
		}
	p->next = pdlbusy;
	pdlbusy = p;
	p->ts = pdts;
	return p->i;
	}

 static int
#ifdef KR_headers
pdlsave(opdb)
	pdl **opdb;
#else
pdlsave(pdl **opdb)
#endif
{
	condlevel++;
	if (iflevel++)
		*opdb = pdlbusy;
	else {
		pdlfree = 0;
		*opdb = 0;
		}
	pdlbusy = 0;
	return ++pdts;
	}

 static void
pdlreset(VOID)
{
	pdl *p, *pnext;
	for(p = pdlbusy; p; p = pnext) {
		pnext = p->next;
		p->next = pdlfree;
		pdlfree = p;
		}
	pdlbusy = 0;
	}

 static void
#ifdef KR_headers
pdlrestore(opdb, mts)
	pdl *opdb;
	int mts;
#else
pdlrestore(pdl *opdb, int mts)
#endif
{
	pdl *p, *pnext;
	condlevel--;
	if (--iflevel) {
		pnext = pdlbusy;
		while(p = pnext) {
			pnext = p->next;
			p->next = opdb;
			opdb = p;
			}
		pnext = pdlfree;
		while((p = pnext) && p->ts >= mts) {
			pnext = p->next;
			p->next = opdb;
			opdb = p;
			}
		pdlfree = pnext;
		pdlbusy = opdb;
		}
	}

#define NDVTGULP 1000

 static void
dvtfree_inc(VOID)
{
	int i, j;
	i = ndvtmax;
	j = ndvtmax += NDVTGULP;
	dvtfree = (int *)Realloc(dvtfree, j*sizeof(int));
	while(i > nvtfree)
		dvtfree[--j] = dvtfree[--i];
	nvtfree += NDVTGULP;
	}

 static int
new_dv(VOID)
{ return ndvfree ? dvtfree[--ndvfree] : ndv++; }

 static void
#ifdef KR_headers
dv_free(k)
	int k;
#else
dv_free(int k)
#endif
{
	if (ndvfree >= nvtfree)
		dvtfree_inc();
	dvtfree[ndvfree++] = k;
	}

 static int
#ifdef KR_headers
new_vt(deriv)
	int deriv;
#else
new_vt(int deriv)
#endif
{
	if (deriv & 2)
		return -(new_pd() + 1);
	return nvtfree < ndvtmax ? dvtfree[nvtfree++] : nvt++;
	}

 static void
#ifdef KR_headers
vt_free(k)
	int k;
#else
vt_free(int k)
#endif
{
	if (ndvfree >= nvtfree)
		dvtfree_inc();
	dvtfree[--nvtfree] = k;
	}

 static list *list_freelist;

 static list *
#ifdef KR_headers
new_list(nxt)
	list *nxt;
#else
new_list(list *nxt)
#endif
{
	list *rv;
	if (rv = list_freelist)
		list_freelist = rv->next;
	else
		rv = (list *)mem(sizeof(list));
	rv->next = nxt;
	return rv;
	}

 static v_i *vi_freelist;

 static v_i *
#ifdef KR_headers
new_vi(v, i)
	real v;
	int i;
#else
new_vi(real v, int i)
#endif
{
	v_i *rv;
	if (rv = vi_freelist)
		vi_freelist = rv->u.next;
	else
		rv = (v_i *)mem(sizeof(v_i));
	rv->u.v = v;
	rv->i = i;
	return rv;
	}

 static void
#ifdef KR_headers
free_vi(vi)
	v_i *vi;
#else
free_vi(v_i *vi)
#endif
{
	vi->u.next = vi_freelist;
	vi_freelist = vi;
	}

 static char *
#ifdef KR_headers
fpval(r)
	real r;
#else
fpval(real r)
#endif
{
	static char buf[32];
	if (r >= Infinity)
		return "1.7e308";
	else if (r <= negInfinity)
		return "-1.7e308";
	g_fmt(buf, r);
	return buf;
	}

 void
#ifdef KR_headers
assign(a, b)
	char *a;
	char *b;
#else
assign(char *a, char *b)
#endif
{
	if (a)
		printf(assign_fmt, a, b);
	}

 void
#ifdef KR_headers
binop(a, b, op, c)
	char *a;
	char *b;
	char *op;
	char *c;
#else
binop(char *a, char *b, char *op, char *c)
#endif
{
	if (a) {
		if (a == b && !Fortran)
			printf("\t%s %s= %s;\n", a, op, c);
		else
			printf(binop_fmt, a, b, op, c);
		}
	}

 void
#ifdef KR_headers
Goto(n)
	int n;
#else
Goto(int n)
#endif
{ printf(goto_fmt, n); }

 void
#ifdef KR_headers
ifgo(a, op, b, n)
	char *a;
	char *op;
	char *b;
	int n;
#else
ifgo(char *a, char *op, char *b, int n)
#endif
{ printf(ifgo_fmt, a, op, b, n); }

 void
#ifdef KR_headers
label(n)
	int n;
#else
label(int n)
#endif
{ printf(label_fmt, n); }

 void
#ifdef KR_headers
ifstart(a, op, b)
	char *a;
	char *op;
	char *b;
#else
ifstart(char *a, char *op, char *b)
#endif
{ printf(ifstart_fmt, a, op, b); }

 void
elsestart(VOID)
{ printf(else_fmt); }

 void
#ifdef KR_headers
elseif(a, b, c)
	char *a;
	char *b;
	char *c;
#else
elseif(char *a, char *b, char *c)
#endif
{ printf(elseif_fmt, a, b, c); }

 void
endif(VOID)
{ printf(endif_fmt); }

 int
#ifdef KR_headers
Switch(v, n)
	char *v;
	int n;
#else
Switch(char *v, int n)
#endif
{
	int i, j, rv;
	if (Fortran) {
		rv = ++branches;
		branches += n;
		printf("\tgoto(");
		for(i = j = 0, --n; i < n; i++) {
			printf("%d,", rv+i);
			if (++j == 9) {
				printf("\n     x\t\t");
				j = 0;
				}
			}
		printf("%d),%s\n", rv+n, v);
		return rv;
		}
	printf("\tswitch(%s) {\n", v);
	return 0;
	}

 void
#ifdef KR_headers
Case(k)
	int k;
#else
Case(int k)
#endif
{ printf(case_fmt, k); }

 void
#ifdef KR_headers
Break(k)
	int k;
#else
Break(int k)
#endif
{ printf(break_fmt, k); }

 void
#ifdef KR_headers
endswitch(lbl)
	int lbl;
#else
endswitch(int lbl)
#endif
{ printf(endswitch_fmt, lbl); }

 void
#ifdef KR_headers
domain(s, x)
	char *s;
	char *x;
#else
domain(char *s, char *x)
#endif
{
	if (Fortran)
		printf("\tcall domain(\"%s\", %s)\n", s, x);
	else
		printf("\tdomain_(\"%s\", &%s, %dL);\n", s, x, strlen(s));
	}

 void
#ifdef KR_headers
zerdiv(s)
	char *s;
#else
zerdiv(char *s)
#endif
{ printf(zerdiv_fmt, s); }

 void
#ifdef KR_headers
call(what, arg)
	char *what;
	char *arg;
#else
call(char *what, char *arg)
#endif
{ printf(call_fmt, what, arg); }

 char *
#ifdef KR_headers
call1(what, a)
	char *what;
	char *a;
#else
call1(char *what, char *a)
#endif
{
	static char buf[48];
	sprintf(buf, "%s(%s)", what, a);
	return buf;
	}

 char *
#ifdef KR_headers
call2(what, a, b)
	char *what;
	char *a;
	char *b;
#else
call2(char *what, char *a, char *b)
#endif
{
	static char buf[64];
	sprintf(buf, "%s(%s, %s)", what, a, b);
	return buf;
	}

 void
#ifdef KR_headers
introuble(who, x)
	char *who;
	char *x;
#else
introuble(char *who, char *x)
#endif
{
	if (!Fortran)
		printf("\tif (errno) in_trouble(\"%s\",%s);\n",
			who, x);
	}

 void
#ifdef KR_headers
introuble2(who, x, y)
	char *who;
	char *x;
	char *y;
#else
introuble2(char *who, char *x, char *y)
#endif
{
	if (!Fortran)
		printf("\tif (errno) in_trouble2(\"%s\",%s,%s);\n",
			who, x, y);
	}

 char *
#ifdef KR_headers
num(x)
	int x;
#else
num(int x)
#endif
{
	static char buf[16];
	sprintf(buf, "%d", x);
	return buf;
	}

 static char op_type[] = {
#include "op_type.hd"
	};

 static void
#ifdef KR_headers
dvset(dp, deriv)
	register real *dp;
	int deriv;
#else
dvset(register real *dp, int deriv)
#endif
{
	register dLR *d;
	real t;

	t = *dp;
	d = Make_dLR(dp);
	if (!deriv) {
		d->kind = dLR_UNUSED;
		return;
		}
	if (t == -1.) {
		d->kind = dLR_negone;
		d->o.vp = &negone;
		}
	else if (t == 1.) {
		d->kind = dLR_one;
		d->o.vp = &one;
		}
	else {
		d->kind = dLR_PD;
		d->o.i = new_pd();
		}
	}

 static void
#ifdef KR_headers
Lset(L, nlin)
	linpart *L;
	int nlin;
#else
Lset(linpart *L, int nlin)
#endif
{
	linpart *Le;
	double t;
	register dLR *d;

	for(Le = L + nlin; L < Le; L++) {
		t = L->fac;
		d = Make_dLR(&L->fac);
		if (t == 1.)
			d->kind = dLR_one;
		else if (t == -1.)
			d->kind = dLR_negone;
		else {
			d->kind = dLR_VP;
			*(d->o.vp = (real *)mem(sizeof(real))) = t;
			}
		}
	}

 static void
#ifdef KR_headers
mpd(e, dp, deriv)
	expr *e;
	real *dp;
	int deriv;
#else
mpd(expr *e, real *dp, int deriv)
#endif
{
	int i, j;
	register dLR *d = Make_dLR(dp);

	if (!deriv) {
		d->kind = dLR_UNUSED;
		if (e->op != f_OPNUM && e->op != f_OPVARVAL && (i = e->a) > 0)
			vt_free(i);
		}
	else if (e->op == f_OPNUM) {
		d->kind = dLR_VP;
		d->o.vp = &((expr_nx *)e)->v;
		i = 0;
		}
	else {
		i = e->a;
		if (e->op == f_OPVARVAL)
			if (i >= nv1) {
				i = (expr_v *)e - var_e;
				if (!(j = cvmap[i -= nv1]))
					j = cvmap[i] = -(++npd);
				d->o.i = -(1 + j);
				d->kind = dLR_PD;
				}
			else {
				d->kind = dLR_VARVAL;
				d->o.i = i;
				}
		else if (i < 0) {
			d->kind = dLR_PD;
			d->o.i = -(1 + i);
			}
		else if (deriv & 2) {
			d->kind = dLR_PD;
			e->a = -((d->o.i = new_pd()) + 1);
			}
		else /* DEBUG */ {
			fprintf(Stderr, "Bug in mpd!\n");
			exit(13);
			}
		}
	}

 static int ewalk(expr*, int);

 static void
#ifdef KR_headers
lwalk(ep, neg)
	expr **ep;
	int neg;
#else
lwalk(expr **ep, int neg)
#endif
{
	int i, j;
	expr *e = *ep;
 top:
	switch(Intcast e->op) {
		case OPNOT:
			neg = 1 - neg;
			e = *ep = e->L.e;
			goto top;
		case OPOR:
			e->op = neg ? f_OPAND : f_OPOR;
			goto andor;
		case OPAND:
			e->op = neg ? f_OPOR : f_OPAND;
 andor:
			lwalk(&e->L.e, neg);
			ep = &e->R.e;
			e = *ep;
			goto top;
		case LT:
			e->op = neg ? f_GE : f_LT;
			goto compare;
		case LE:
			e->op = neg ? f_GT : f_LE;
			goto compare;
		case EQ:
			e->op = neg ? f_NE : f_EQ;
			goto compare;
		case GE:
			e->op = neg ? f_LT : f_GE;
			goto compare;
		case GT:
			e->op = neg ? f_LE : f_GT;
			goto compare;
		case NE:
			e->op = neg ? f_EQ : f_NE;
 compare:
			i = ewalk(e->L.e,0);
			j = ewalk(e->R.e,0);
			if (i > 0)
				vt_free(i);
			if (j > 0)
				vt_free(j);
		}
	}

 static char *
#ifdef KR_headers
f_OPVARVAL1(e, buf)
	expr *e;
	char *buf;
#else
f_OPVARVAL1(expr *e, char *buf)
#endif
{
	int k;
	Adjoint *A;
	char *fmt;

	if ((k = e->a) < nv1) {
		if (k < 0) {
			fmt = pd_fmt;
			k = Fortran1 - k;
			}
		else {
			A = Adjp(&adjoints[k]);
			if (!A->seen) {
				A->seen = 1;
				vseen[nvseen++] = e->a;
				}
			fmt = x_fmt;
			k += Fortran;
			}
		}
	else {
		k = cvmap[(expr_v *)e - var_e - nv1];
		if (k < 0) {
			fmt = pd_fmt;
			k = Fortran1 - k;
			}
		else {
			fmt = tv_fmt;
			k += Fortran1;
			}
		}
	sprintf(buf, fmt, k);
	return buf;
	}

 char *
#ifdef KR_headers
f_OPNUM1(e, rv)
	register expr *e;
	char *rv;
#else
f_OPNUM1(register expr *e, char *rv)
#endif
{
	g_fmt(rv,((expr_nx *)e)->v);
	return rv;
	}

 static int
#ifdef KR_headers
ewalk(e, deriv)
	expr *e;
	int deriv;
#else
ewalk(expr *e, int deriv)
#endif
{
	int deriv1, i, j, k, k1, lderiv, mts, rderiv;
	expr *e1, **ep, **epe;
	expr_if *eif;
	expr_va *eva;
	expr_f *ef;
	argpair *ap, *ape;
	de *d;
	derp *dp;
	dLR *LR, **LRp;
	efunc *op;
	Adjoint *A;
	double t;
	pdl *opdb;
	static int achk[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 0, 0, 10, 11 };

	k = Intcast e->op;
	e->op = r_op[k];
	deriv1 = deriv & 1;
	if (achk[k1 = op_type[k]] && e->a == nv1)
		deriv = deriv1 = 0;
	switch(k1) {
		case 11: /* OP2POW */
			e->dL = 0;
			/* no break */
		case 1: /* unary */
			j = ewalk(e->L.e, deriv1);
			i = new_vt(deriv);
			if (deriv)
				dvset(k == OPCPOW ? &e->dR : &e->dL, deriv);
			if (j > 0)
				vt_free(j);
			seen[k] = 1;
 reti:
			return e->a = i;

		case 2: /* binary */
			k1 = 1;
			switch(k) {
				case OPPLUS:
				case OPMINUS:
				case OPREM:
					e->dL = 1.;
					break;
				case OPMULT:
					k1 = 3;
					/* no break */
				default:
					e->dL = 0.;
				}
			if (lderiv = rderiv = deriv1) {
				lderiv = achk[op_type[Intcast e->L.e->op]]
					&& e->L.e->a != nv1 ? k1 : 0;
				rderiv = achk[op_type[Intcast e->R.e->op]]
					&& e->R.e->a != nv1 ? k1 : 0;
				}
			i = ewalk(e->L.e,lderiv);
			j = ewalk(e->R.e,rderiv);
			if (deriv) {
				switch(k) {
					case OPREM:
					case OPintDIV:
					case OPprecision:
					case OPround:
					case OPtrunc:
						rderiv = 0;
						break;
					case OPMULT:
						mpd(e->R.e, &e->dL, lderiv);
						mpd(e->L.e, &e->dR, rderiv);
						goto retnew;
					case OP_atan2:
						needT1 = 1;
					}
				dvset(&e->dL, lderiv);
				dvset(&e->dR, rderiv);
				k = new_vt(deriv1);
				if (i > 0)
					vt_free(i);
				if (j > 0)
					vt_free(j);
				return e->a = k;
				}
			k = i;
			i = new_vt(0);
			if (j > 0)
				vt_free(j);
			if (k > 0)
				vt_free(k);
			goto reti;

		case 3: /* vararg (min, max) */
			eva = (expr_va *)e;
			d = eva->L.d;
			condlevel++;
			if (!(i = ewalk(d->e,deriv)))
				i = new_vt(deriv);
			while(e1 = (++d)->e) {
				pdlreset();
				if ((j = ewalk(e1,deriv1)) > 0)
					vt_free(j);
				}
			condlevel--;
			if (eva->a != nv1) {
				eva->next = (expr_va *)ncond++;
				/* arrange to find this expr_va */
				/* when walking derps */
				dp = eva->R.D;
				LRp = Make_dLRp(dp->c.rp);
				*LRp = LR = (dLR *)mem(sizeof(dLR));
				LR->kind = dLR_VARARG;
				LR->o.eva = eva;
				}
			goto reti;

		case 4: /* piece-wise linear */
			if (deriv) {
				e->dL = 0;
				dvset(&e->dL, 1);
				}
			LR = Make_dLR(&e->dR);
			LR->kind = nplterm++;
			LR->o.ep = Plterms;
			Plterms = e;
 retnew:
			return e->a = new_vt(deriv);

		case 5: /* if */
			eif = (expr_if *)e;
			eif->next = (expr_if*) -1;
			lwalk(&eif->e, 0);
			mts = pdlsave(&opdb);
			if ((i = ewalk(eif->T,deriv1)) > 0)
				vt_free(i);
			pdlreset();
			if ((i = ewalk(eif->F,deriv1)) > 0)
				vt_free(i);
			pdlrestore(opdb, mts);
			if (eif->a != nv1) {
				eif->next = (expr_if *)ncond++;
				/* arrange to find this expr_if */
				/* when walking derps */
				dp = eif->D;
				LRp = Make_dLRp(dp->c.rp);
				*LRp = LR = (dLR *)mem(sizeof(dLR));
				LR->kind = dLR_IF;
				LR->o.eif = eif;
				}
			goto retnew;

		case 6: /* sumlist */
			ep = e->L.ep;
			epe = e->R.ep;
			i = ewalk(*ep++, deriv);
			if (i < 0)
				deriv = deriv1;
			j = ewalk(*ep++, deriv);
			if (i > 0) {
				if (j > 0)
					vt_free(j);
				}
			else if (!(i = j))
				i = new_vt(deriv);
			do {
				if ((j = ewalk(*ep++, deriv1)) > 0)
					vt_free(j);
				}
				while(ep < epe);
			goto reti;

		case 7: /* function call */
			ef = (expr_f *)e;
			i = k = 0;
			for(ap = ef->ap, ape = ef->ape; ap < ape; ap++)
				if (j = ewalk(ap->e,deriv))
					if (j < 0) {
						i = j;
						deriv = deriv1;
						}
					else
						k++;
			if (k)
				for(ap = ef->ap; ap < ape; ap++) {
					op = ap->e->op;
					if (op != f_OPNUM
					 && op != f_OPHOL
					 && op != f_OPVARVAL
					 && (j = ap->e->a) > 0)
						vt_free(j);
					}
			if (!i)
				goto retnew;
			goto reti;

		case 8: /* Hollerith */
			break;

		case 9: /* number */
			t = ((expr_n *)e)->v;
			((expr_nx *)e)->v = t;
			((expr_nx *)e)->next = nums;
			nums = (expr_nx *)e;
			break;

		case 10: /* variable value */
			i = (expr_v *)e - var_e;
			A = Adjp(&adjoints[e->a]);
			A->ifset = condlevel ? 1 : 0;
			if (!A->seen) {
				A->seen = 1;
				vseen[nvseen++] = i;
				}
			break;
		/*DEBUG*/default:
		/*DEBUG*/ fprintf(Stderr, "bad opnumber %d in ewalk\n", k);
		/*DEBUG*/ exit(1);
		}
	return 0;
	}

 static void
#ifdef KR_headers
max_save(m)
	maxinfo *m;
#else
max_save(maxinfo *m)
#endif
{
	if (nvtmax < nvt)
		nvtmax = nvt;
	m->nvt = nvtmax - 1;
	m->ndv = ndvmax;
	m->ncond = ncond;
	m->needT1 = needT1;
	}

 static void
#ifdef KR_headers
max_restore(m)
	maxinfo *m;
#else
max_restore(maxinfo *m)
#endif
{
	nvtmax = m->nvt + 1;
	ndvmax = m->ndv;
	ncond = m->ncond;
	needT1 = m->needT1;
	}

 static void
#ifdef KR_headers
vreset(ifset)
	list **ifset;
#else
vreset(list **ifset)
#endif
{
	Adjoint *A;
	list *lastif = 0;
	int k;

	while(nvseen) {
		k = vseen[--nvseen];
		A = Adjp(&adjoints[var_e[k].a]);
		if (A->ifset)
			(lastif = new_list(lastif))->item.i = k;
		*A = A0;
		}
	*ifset = lastif;
	}

 static void
#ifdef KR_headers
ewalkvt(e, ifset)
	expr *e;
	list **ifset;
#else
ewalkvt(expr *e, list **ifset)
#endif
{
	nvtfree = ndvtmax;
	ewalk(e,want_derivs);
	vreset(ifset);
	}

 static void
#ifdef KR_headers
zset(z)
	register int *z;
#else
zset(register int *z)
#endif
{
	register int *ze;
	if (z)
		for(ze = z + 2**z; z < ze; z += 2)
			var_e[z[1]].a = z[2];
	}

 static void
#ifdef KR_headers
ndvtreset(z)
	int *z;
#else
ndvtreset(int *z)
#endif
{
	nvt = 1;
	nvtfree = ndvtmax;
	ndv = Fortran;
	ndvfree = 0;
	zset(z);
	}

 static void
#ifdef KR_headers
comwalk(i, n)
	int i;
	int n;
#else
comwalk(int i, int n)
#endif
{
	cexp *c, *ce;
	list **ifset;

	if (i >= n)
		return;
	zset(zaC[i]);
	ifset = com_ifset + i;
	c = cexps + i;
	for(ce = cexps + n; c < ce; c++) {
		Lset(c->L, c->nlin);
		ndvtreset(zaC[i]);
		ewalkvt(c->e, ifset++);
		if (nvtmax < nvt)
			nvtmax = nvt;
		}
	}

 static void dwalk0(derp *, derp*);
 static void dwalk1(derp *, derp*);

 static void
#ifdef KR_headers
dwalk0(d, d0)
	derp *d;
	derp *d0;
#else
dwalk0(derp *d, derp *d0)
#endif
{
	int i;

	if (d == d0)
		return;
	do {
		if ((i = d->b.rp - adjoints - nv1) > 0)
			acount[i]++;
		}
		while((d = d->next) != d0);
	}

 static void
#ifdef KR_headers
cond_magic(LR)
	dLR *LR;
#else
cond_magic(dLR *LR)
#endif
{
	expr_va *eva;
	expr_if *eif;
	derp *D, *D0, Dsave;
	de *d;
	int iftype;

	if (LR->kind == dLR_IF) {
		eif = LR->o.eif;
		D = eif->D;
		iftype = 1;
		}
	else {
		eva = LR->o.eva;
		D = eva->R.D;
		iftype = 0;
		}
	Dsave = *D;
	D->c.rp = &edagread_one;
	if (iftype) {
		D->a.rp = eif->Tv.rp;
		D->next = eif->dT;
		dwalk1(D, eif->d0);
		D->a.rp = eif->Fv.rp;
		D->next = eif->dF;
		dwalk1(D, eif->d0);
		}
	else {
		D0 = eva->d0;
		for(d = eva->L.d; d->e; d++) {
			D->a.rp = d->dv.rp;
			D->next = d->d;
			dwalk1(D, D0);
			}
		*D = Dsave;
		}
	*D = Dsave;
	}

#ifdef DEBUG
 /*debug*/static int awrite, derps_in, derps_out, stats, vwrite;
#endif

 static void
#ifdef KR_headers
dwalk1(d, d0)
	derp *d;
	derp *d0;
#else
dwalk1(derp *d,  derp *d0)
#endif
{
	Adjoint *a, *b;
	dLR *c;
	int i, j, j0, neg;

	if (d == d0)
		return;
	dwalk0(d, d0);
	do {
		j = -1;
		if ((j0 = d->b.rp - adjoints - nv1) > 0)
			j = --acount[j0];
		b = Adjp(d->b.rp);
		if (!b->storage)	/* defined var used only in if */
			continue;
		i = d->a.rp - adjoints;
		if (maxa < i)
			maxa = i;
		c = dLRp(*d->c.rp);
		if (c->kind >= dLR_VARARG) {
			if (!output_time)
				cond_magic(c);
			goto jfree;
			}
		if (i < nv1 && !fwalk)
			continue;
		a = Adjp(d->a.rp);
		switch(a->storage) {
			case STOR_VP:
				a->storage = STOR_DV;
				a->o.i = new_dv();
			default:
				goto jfree;
			case 0:
				break;
			}
		if (b->storage == STOR_IMPLICIT) {
			a->neg = b->neg;
			switch(c->kind) {
			  case dLR_negone:
				a->neg = 1 - b->neg;
				/* no break */
			  case dLR_one:
				a->storage = STOR_IMPLICIT;
				break;
			  case dLR_VP:
				a->storage = STOR_VP;
				a->o.vp = c->o.vp;
				break;
			  case dLR_PD:
				a->storage = STOR_PD;
				a->o.i = c->o.i;
				break;
			  case dLR_VARVAL:
				a->storage = STOR_VARVAL;
				a->o.i = c->o.i;
				break;
			  default:
				/*DEBUG*/ fprintf(Stderr,
				"\nBad c->kind = %d with STOR_IMPLICIT\n",
						c->kind);
				/*DEBUG*/ exit(10);
			  }
 jfree:
			if (!j && b->storage == STOR_DV)
				dv_free(b->o.i);
			continue;
			}
		neg = b->neg;
		switch(c->kind) {
		  case dLR_negone:
			neg = 1 - neg;
			/* no break */
		  case dLR_one:
			if (j || fwalk && b->storage == STOR_PD)
				break;
			a->storage = b->storage;
			a->o = b->o;
			a->neg = neg;
			continue;
		  }
		a->storage = STOR_DV;
		a->o.i = !j && b->storage == STOR_DV ? b->o.i : new_dv();
		}
		while((d = d->next) != d0);
	if (b && b->storage == STOR_DV)
		dv_free(b->o.i);
	}

 static void
#ifdef KR_headers
areset(cg)
	cgrad *cg;
#else
areset(cgrad *cg)
#endif
{
	Adjoint *A;
	int k;
	while(nvseen) {
		A = Adjp(&adjoints[var_e[vseen[--nvseen]].a]);
		*A = A0;
		}

	for(k = nv1; k <= maxa; k++) {
		A = Adjp(&adjoints[k]);
		*A = A0;
		}
	maxa = 0;
	for(; cg; cg = cg->next) {
		A = Adjp(&adjoints[cg->varno]);
		*A = A0;
		}
	}

 static void
#ifdef KR_headers
dwalk(d, cg, z, L)
	register derp *d;
	cgrad *cg;
	int *z;
	list *L;
#else
dwalk(register derp *d, cgrad *cg, int *z, list *L)
#endif
{
	register real *rp;
	register derp *d1;
	int i;
	Adjoint *A;
	static cgrad *cg0;

	areset(cg0);
	cg0 = cg;
	ndv = Fortran;
	ndvfree = 0;
	zset(z);
	for(; L; L = L->next) {
		i = L->item.i;
		A = Adjp(&adjoints[var_e[i].a]);
		A->seen = 1;
		vseen[nvseen++] = i;
		A->storage = STOR_DEFV;
		A->o.i = ndv++;
		}
	if (d1 = d) {
		dwalk_one->storage = STOR_IMPLICIT;
		rp = d->b.rp;
		do {
			if (d1->b.rp == rp)
				d1->b.rp = rdwalk_one;
			}
			while(d1 = d1->next);
		dwalk1(d, 0);
		}
	if (ndvmax < ndv)
		ndvmax = ndv;
	}

 static void
#ifdef KR_headers
funnelwalk(f)
	funnel *f;
#else
funnelwalk(funnel *f)
#endif
{
	cplist *cl;
	dLR *d;
	int k;

	if (!want_derivs)
		return;
	for(; f; f = f->next) {
		k = f->ce - cexps;
		dwalk(f->fcde.d, 0, zaC[k], com_ifset[k]);
		for(cl = f->cl; cl; cl = cl->next) {
			d = Make_dLR(cl->cfa);
			d->kind = dLR_PD;
			d->o.i = npd++;
			}
		}
	}

 static void
#ifdef KR_headers
cde_walk(d, n, m, ifset, c1st, z)
	cde *d;
	int n;
	maxinfo *m;
	list **ifset;
	int *c1st;
	int **z;
#else
cde_walk(cde *d, int n, maxinfo *m, list **ifset, int *c1st, int **z)
#endif
{
	cde *De;
	int i, i1, j, k;
	cexp1 *c1;

	j = *c1st++;
	for(De = d + n; d < De; d++, ifset++) {
		i = j;
		k = *c1st++;
		ndvtreset(*z++);
		while(j < k) {
			c1 = cexps1 + j++;
			Lset(c1->L, c1->nlin);
			ewalk(c1->e, want_derivs);
			}
		ewalkvt(d->e, ifset);
		while(i < k)
			if (!cvmap[i1 = i++ + ncom0])
				cvmap[i1] = nvt++;
		if (nvtmax < nvt)
			nvtmax = nvt;
		if (want_derivs)
			dwalk(d->d, 0, 0, *ifset);
		}
	max_save(m);
	}

 static void
#ifdef KR_headers
cde_dwalk(d, n, m, ifset, z, cgp)
	cde *d;
	int n;
	maxinfo *m;
	list **ifset;
	int **z;
	cgrad **cgp;
#else
cde_dwalk(cde *d, int n, maxinfo *m, list **ifset, int **z, cgrad **cgp)
#endif
{
	cde *De;

	if (!want_derivs)
		return;
	for(De = d + n; d < De; d++, ifset++) {
		dwalk(d->d, *cgp++, *z++, *ifset);
		if (m->ndv < ndv)
			m->ndv = ndv;
		}
	}

 static void derivs(derp *, derp *);
 static list **dstored;

 static void
#ifdef KR_headers
dstore(a, d)
	Adjoint *a;
	derp *d;
#else
dstore(Adjoint *a, derp *d)
#endif
{
	a->stored = 1;
	if (dstored)
		(*dstored = new_list(*dstored))->item.D = d;
	}

 static void
derprestore(VOID)
{
	Adjoint *a;
	list *L, *Lnext;

	if (L = *dstored) {
		*dstored = 0;
		do {
			Lnext = L->next;
			a = Adjp(L->item.D->a.rp);
			a->stored = 0;
			L->next = list_freelist;
			list_freelist = L;
			}
			while(L = Lnext);
		}
	}

 static void
#ifdef KR_headers
cond_dmagic(LR)
	dLR *LR;
#else
cond_dmagic(dLR *LR)
#endif
{
	expr_va *eva;
	expr_if *eif;
	derp *D, *D0;
	de *d;
	char cbuf[16];
	int i, j, n;
	list **ds0, *dstore1;

	ds0 = dstored;
	dstored = &dstore1;
	dstore1 = 0;
	if (LR->kind == dLR_IF) {
		eif = LR->o.eif;
		D = eif->D;
		D->c.rp = &edagread_one;
		sprintf(cbuf, cond_fmt, Intcast eif->next);
		if (eif->dT != eif->d0) {
			printf(ifstart_fmt, cbuf, opNE, "0");
			D->a.rp = eif->Tv.rp;
			D->next = eif->dT;
			derivs(D, eif->d0);
			derprestore();
			if (eif->dF != eif->d0) {
				printf(else_fmt);
 just_else:
				D->a.rp = eif->Fv.rp;
				D->next = eif->dF;
				derivs(D, eif->d0);
				derprestore();
				}
			}
		else {
			printf(ifstart_fmt, cbuf, opEQ, "0");
			goto just_else;
			}
#ifdef DEBUG
		--derps_in;
#endif
		printf(endif_fmt);
		D->next = eif->d0;
		}
	else {
		eva = LR->o.eva;
		D = eva->R.D;
		D->c.rp = &edagread_one;
		sprintf(cbuf, cond_fmt, Intcast eva->next);
		for(n = 0, d = eva->L.d; d->e; d++)
			n++;
		i = Switch(cbuf, n);
		j = i + n;
		D0 = eva->d0;
#ifdef DEBUG
		++derps_in;
#endif
		for(d = eva->L.d; ;) {
			Case(i++);
			D->a.rp = d->dv.rp;
			D->next = d->d;
#ifdef DEBUG
			--derps_in;
#endif
			derivs(D, D0);
			derprestore();
			if (!(++d)->e)
				break;
			Break(j);
			}
		endswitch(j);
		}
	dstored = ds0;
	}

 static void
#ifdef KR_headers
derivs(d, d0)
	derp *d;
	derp *d0;
#else
derivs(derp *d, derp *d0)
#endif
{
	Adjoint *a, *b;
	dLR *c;
	char *aop, *fmt, *ginc, *mult;
	int havenum, i, neg, num2chk, stored;
	derp *dnext;
	char bufb[32], bufc[32], bufg[32], bufi[32];
	double t, t1;
	v_i *vi;

	for(; d != d0; d = dnext) {
		dnext = d->next;
		c = dLRp(*d->c.rp);
		if (c->kind >= dLR_VARARG) {
			cond_magic(c);
			cond_dmagic(c);
			continue;
			}
#ifdef DEBUG
		++derps_in;
#endif
		b = Adjp(d->b.rp);
		if (!b->storage)	/* defined var used only in if */
			continue;
		neg = b->neg;
		a = Adjp(d->a.rp);
		if (c->kind == dLR_one
		 && a->storage == b->storage && a->o.i == b->o.i) {
			switch(a->storage) {
		  		case STOR_GRAD:
		 		case STOR_JAC:
				case STOR_DEFV:
					break;
				default:
					dstore(a,d);
					a->neg = neg;
					continue;
				}
			}
		num2chk = 0;
		stored = a->stored;
		switch(a->storage) {
		  case STOR_VI:
			i = a->o.vi->i;
			fmt = dv_fmt;
			break;
		  case STOR_PD:
			if (b->storage == STOR_IMPLICIT
			 && (stor_grad != STOR_DEFV
				|| a->storage == STOR_PD
					&& c->kind == dLR_PD
					&& a->o.i == c->o.i))
				continue;
			stored = a->stored = 1;
			if (a->neg)
				neg = 1 - neg;
			fmt = pd_fmt;
			i = a->o.i + Fortran;
			break;
		  case STOR_DV:
		  case STOR_DEFV:
			num2chk = 1;
			i = a->o.i;
			fmt = dv_fmt;
			break;
		  case STOR_GRAD:
			i = a->o.og->varno + Fortran;
			fmt = grad_fmt;
			break;
		  case STOR_JAC:
			i = (djoff >= 0
				? djoff + a->o.cg->varno*n_con
				: a->o.cg->goff) + Fortran;
			fmt = jac_fmt;
			break;
		  default:
			continue;
			}
		mult = "*";
		havenum = 0;
		switch(b->storage) {
			case STOR_PD:
				sprintf(bufb, pd_fmt, b->o.i + Fortran);
				break;
			case STOR_DV:
			case STOR_DEFV:
				sprintf(bufb, dv_fmt, b->o.i);
				break;
			case STOR_VARVAL:
				sprintf(bufb, x_fmt, b->o.i + Fortran);
				break;
			case STOR_VI:
				t = b->o.vi->u.v;
				goto have_t;
			case STOR_VP:
				t = *b->o.vp;
 have_t:
				havenum = 1;
				if (t < 0) {
					t = -t;
					neg = 1 - neg;
					}
				g_fmt(bufb, t);
				break;
			case STOR_IMPLICIT:
				havenum = 1;
				t = 1.;
				mult = "";
				bufb[0] = 0;
				break;
			default:
				/*DEBUG*/ fprintf(Stderr,
				"\nBad b->storage = %d\n", b->storage);
				/*DEBUG*/ exit(12);
				}
		switch(c->kind) {
			case dLR_VP:
				if (havenum) {
					t *= *c->o.vp;
					if (!stored && num2chk) {
						havenum = 2;
						break;
						}
					mult = "";
					bufb[0] = 0;
					}
				else
					t = *c->o.vp;
 num_fmt:
				if (t < 0) {
					t = -t;
					neg = 1 - neg;
					}
				g_fmt(bufc, t);
				break;
			case dLR_PD:
				sprintf(bufc, pd_fmt, c->o.i + Fortran);
				break;
			case dLR_VARVAL:
				sprintf(bufc, x_fmt, c->o.i + Fortran);
				break;
			case dLR_negone:
				neg = 1 - neg;
			case dLR_one:
				if (havenum) {
					if (!stored) {
						havenum = 2;
						if (*mult)
							break;
						}
					goto num_fmt;
					}
				if (*mult) {
					mult = "";
					bufc[0] = 0;
					}
				else
					strcpy(bufc, One);
				break;
			default:
				/*DEBUG*/fprintf(Stderr,
				"\nBad c->kind = %d\n", c->kind);
				/*DEBUG*/ exit(11);
			}
		ginc = 0;
		sprintf(bufg, fmt, i);
		if (stored) {
			aop = neg ? "-" : "+";
			if (Fortran)
				printf("\t%s = %s %s ", bufg, bufg, aop);
			else
				printf("\t%s %s= ", bufg, aop);
			}
		else {
			if (neg && havenum)
				t = -t;
			switch(a->storage) {
				case STOR_DV:
				case STOR_DEFV:
					if (havenum == 2) {
						a->storage = STOR_VI;
						a->o.vi = new_vi(t, a->o.i);
						a->neg = 0;
						continue;
						}
					break;
				case STOR_VI:
					if (havenum == 2) {
						a->o.vi->u.v += t;
						continue;
						}
					a->storage = STOR_DV;
					a->o.i = (vi = a->o.vi)->i;
					if (vi->u.v)
						g_fmt(ginc = bufi, vi->u.v);
					sprintf(bufg, dv_fmt, vi->i);
					free_vi(vi);
					break;
				case STOR_GRAD:
					t1 = a->o.og->coef;
					goto t1_test;
				case STOR_JAC:
					t1 = a->o.cg->coef;
 t1_test:
					if (t1)
						g_fmt(ginc = bufi, t1);
				}
			dstore(a,d);
			if (!*mult && !neg
			 && (!bufb[0] && !strcmp(bufc,bufg))
			  || !bufc[0] && !strcmp(bufb,bufg))
				continue;
			if (!strcmp(bufg,bufb) && !Fortran
			 && *mult == '*' && !neg && !ginc) {
				printf("\t%s *= %s;\n", bufg, bufc);
#ifdef DEBUG
				++derps_out;
				if (awrite)
					goto a_write;
#endif
				continue;
				}
			printf("\t%s %s", bufg, neg ? "= -" : "= ");
			}
		printf("%s%s%s", bufb, mult, bufc);
		if (ginc)
			printf(" + %s", ginc);
		printf(eos);
#ifdef DEBUG
		++derps_out;
		if (awrite) {
 a_write:
		if (!*bufb)
			strcpy(bufb, One);
		if (!*bufc)
			strcpy(bufc, One);
		  if (Fortran) {
			printf("\twrite(*,*) '%d: %s', %s, ' * ', %s%s",
				derps_in, neg ? "-(" : "", bufb, bufc,
				neg ? ",')'" : "");
			if (ginc)
				printf(", ' +  %s'", ginc);
			printf(", ': a[%d] = ', %s\n",
				d->a.rp - adjoints, bufg);
		  } else {
			printf("fprintf(Stderr, \"%d: %s", derps_in,
				neg ? "-(%g * %g)" : "%g * %g");
			if (ginc)
				printf(" + %s", ginc);
			printf(": a[%d] = %%g\\n\", %s, %s, %s);\n",
				d->a.rp - adjoints, bufb, bufc, bufg);
		  }
			}
#endif
		}
	}

 static void
#ifdef KR_headers
co_derivs(d, L, gj0_fmt)
	derp *d;
	list *L;
	char *gj0_fmt;
#else
co_derivs(derp *d, list *L, char *gj0_fmt)
#endif
{
	char buf[32];
	Adjoint *a;
	int j, k;
	ograd *og;

	for(j = Fortran; L; L = L->next, j++) {
		a = Adjp(&adjoints[var_e[k = L->item.i].a]);
		if (k < nv1)
			if (stor_grad == STOR_DEFV) {
				sprintf(buf, dv_fmt, j + Fortran);
				assign(buf, Zero);
				}
			else {
				og = stog[k];
				g_fmt(buf, stor_grad == STOR_GRAD ? og->coef
						: ((cgrad *)og)->coef);
				printf(gj0_fmt, k + Fortran, buf);
				}
		else {
			a->storage = STOR_DEFV;
			a->o.i = j;
			sprintf(buf, dv_fmt, j);
			assign(buf, Zero);
			}
		a->stored = 1;
		}
	derivs(d,0);
	}

 static char *
#ifdef KR_headers
vprod(t, k)
	real t;
	int k;
#else
vprod(real t, int k)
#endif
{
	static char buf[64];
	int i;

	k += Fortran;
	if (t == 1.)
		sprintf(buf, x_fmt, k);
	else if (t == -1.) {
		buf[0] = '-';
		sprintf(buf+1, x_fmt, k);
		}
	else {
		i = g_fmt(buf, t);
		buf[i++] = '*';
		sprintf(buf+i, x_fmt, k);
		}
	return buf;
	}

 static void
#ifdef KR_headers
plcommon(deriv)
	int deriv;
#else
plcommon(int deriv)
#endif
{
	register plterm *p;
	register expr *e;
	register dLR *LR;

	for(e = Plterms; e; e = LR->o.ep) {
		p = e->L.p;
		LR = dLRp(e->dR);
		printf("\tcommon /pltcom/ bs%d\n\
	double precision bs%d(%d)\n", LR->kind, LR->kind, 2*p->n - 1);
		}
	printf(
	 "\tcommon /xkindc/ xkind\n\tinteger xkind\n\tsave /xkindc/\n");

	if (npd)
		printf(
	"\tcommon /pdcomn/ pd\n\tdouble precision pd(%d)\n\tsave /pdcomn/\n",
			npd);
	if (!deriv && needx0check)
		printf("\tlogical xchk\n");
	}

 static char *
#ifdef KR_headers
rv_output(rv, eval, og)
	char *rv;
	char *eval;
	ograd *og;
#else
rv_output(char *rv, char *eval, ograd *og)
#endif
{
	char *s;

	for(; og; og = og->next)
		if (og->coef)
			break;
	if (og) {
		s = vprod(og->coef, og->varno);
		if (strcmp(eval, Zero))
			binop(rv, eval, "+", s);
		else
			assign(rv, s);
		while(og = og->next)
			if (og->coef)
				binop(rv, rv, "+",
					vprod(og->coef, og->varno));
		return rv;
		}
	if (Fortran)
		assign(rv, eval);
	return eval;
	}

 static void
#ifdef KR_headers
obj_derivs(i)
	int i;
#else
obj_derivs(int i)
#endif
{
	ograd *og, *og0;
	Adjoint *a;
	cde *c;
	char buf[32];

	if (!want_derivs)
		return;
	c = obj_de + i;
	og0 = Ograd[i];
	dwalk(c->d, (cgrad*)og0, zao[i], o_ifset[i]);	/* added 20010312 */
	for(og = og0; og; og = og->next) {
		a = Adjp(&adjoints[og->varno]);
		a->storage = STOR_GRAD;
		a->o.og = og;
		a->stored = 0;
		stog[og->varno] = og;
		}
	stor_grad = STOR_GRAD;
	co_derivs(c->d, o_ifset[i], g0_fmt);
	for(og = og0; og; og = og->next) {
		a = Adjp(&adjoints[og->varno]);
		if (!a->stored && (!intsk || !intsk[og->varno])) {
			g_fmt(buf, og->coef);
			printf(g0_fmt, og->varno + Fortran, buf);
			}
		}
	}

 static char *
#ifdef KR_headers
con_linadd(i, s)
	int i;
	char *s;
#else
con_linadd(int i, char *s)
#endif
{
	cgrad *cg;
	char *s1;

	for(cg = Cgrad[i]; cg; cg = cg->next)
		if (cg->coef) {
			s1 = vprod(cg->coef, cg->varno);
			if (strcmp(s,Zero))
				binop(T, s, "+", s1);
			else
				assign(T, s1);
			s = T;
			while(cg = cg->next)
				if (cg->coef)
					binop(s, s, "+",
						vprod(cg->coef, cg->varno));
			break;
			}
	return s;
	}

 static void
#ifdef KR_headers
con_derivs(i)
	int i;
#else
con_derivs(int i)
#endif
{
	cgrad *cg, *cg0;
	Adjoint *a;
	char buf[32];
	int i1, j;
	cde *c = con_de + i;

	if (densejac)
		djoff = i;
	printf("\n%s   /*** derivatives for constraint %d ***/\n\n",
		Fortstar, i+1);
	dwalk(c->d, cg0 = Cgrad[i], zac[i], c_ifset[i]);
	for(cg = cg0; cg; cg = cg->next) {
		a = Adjp(&adjoints[cg->varno]);
		a->storage = STOR_JAC;
		a->stored = 0;
		a->o.cg = cg;
		stog[cg->varno] = (ograd *)cg;
		}
	stor_grad = STOR_JAC;
	co_derivs(con_de[i].d, c_ifset[i], j0_fmt);
	for(cg = cg0; cg; cg = cg->next) {
		a = Adjp(&adjoints[cg->varno]);
		if (!a->stored && (!intsk || !intsk[cg->varno])) {
			g_fmt(buf, cg->coef);
			j = densejac ? i + cg->varno*n_con : cg->goff;
			printf(j0_fmt, j + Fortran, buf);
			}
		}
	if (densejac) {
		i1 = i + 1;
		for(cg = cg0; cg; cg = cg->next)
			rownos[cg->varno] = i1;
		for(j = 0; j < nv1; j++)
			if (rownos[j] != i1)
				printf(j0_fmt, i + j*n_con + Fortran, Zero);
		}
	}

 static int
nzcgrad(VOID)
{
	cgrad *cg;
	int i;
	for(i = 0; i < n_con; i++)
		for(cg = Cgrad[i]; cg; cg = cg->next)
			if (cg->coef)
				return 1;
	return 0;
	}

 static char *
#ifdef KR_headers
cv_name(L, buf)
	linpart *L;
	char *buf;
#else
cv_name(linpart *L, char *buf)
#endif
{
	expr_v ev;
	expr *ep;

	ep = (expr *)((char *)L->v.rp - ((char *)&ev.v - (char *)&ev));
	return f_OPVARVAL1(ep, buf);
	}

 static void
#ifdef KR_headers
com_out(e, L, nlin, k0)
	expr *e;
	linpart *L;
	int nlin;
	int k0;
#else
com_out(expr *e, linpart *L, int nlin, int k0)
#endif
{
	char buf[32], bufg[32], res[32], vn[32];
	char *s;
	double t;
	linpart *Le;
	int asg, j, k, op;
	efuncb *eb;
	dLR *d;

	printf("\n%s\t/*** defined variable %d ***/\n\n", star, k0+1);
	if ((j = cvmap[k0]) < 0) {
		j = Fortran1 - j;
		s = pd_fmt;
		}
	else {
		eb = (efuncb *)e->op;
		if (eb != (efuncb *)OPNUM && eb != (efuncb *)OPVARVAL)
			e->a = j;
		j += Fortran1;
		s = tv_fmt;
		}
	sprintf(res, s, j);
	s = callb(e,buf);
	if (!L) {
		if (strcmp(res,s))
			assign(res, s);
#ifdef DEBUG
		goto done;
#else
		return;
#endif
		}
	asg = !strcmp(s, Zero);
	for(Le = L + nlin; L < Le; L++, asg = 0) {
		d = dLRp(L->fac);
		op = '+';
		switch(k = d->kind) {
			case dLR_negone:
				op = '-';
			case dLR_one:
				break;
			case dLR_VP:
				t = *d->o.vp;
				if (t < 0. && !asg) {
					t = -t;
					op = '-';
					}
				g_fmt(bufg, t);
				break;
			default:/*DEBUG*/ fprintf(Stderr,
				"Bad d->kind = %d in com_walk\n", d->kind);
				/*DEBUG*/ exit(14);
			}
		if (asg)
			printf(op == '-' ? "\t%s = -" : "\t%s = ", res);
		else if (Fortran || res != s)
			printf("\t%s = %s %c ", res, s, op);
		else
			printf("\t%s %c= ", res, op);
		if (k == dLR_VP)
			printf("%s*", bufg);
		printf("%s%s", cv_name(L,vn), eos);
		s = res;
		}
#ifdef DEBUG
 done:
	if (vwrite)
	printf("fprintf(Stderr, \"defvar[%d] = %%g\\n\", %s);\n", k0, res);
#endif
	}

 static char *
#ifdef KR_headers
putout(e, i, j, what, k, z)
	expr *e;
	int i;
	int j;
	char *what;
	int k;
	int *z;
#else
putout(expr *e, int i, int j, char *what, int k, int *z)
#endif
{
	static char buf[32];
	cexp1 *c1;

	++k;
	ndvtreset(z);
	if (i < j) {
		if (what)
		printf("\n\n%s\t/*** defined variables for %s %d ***/\n",
			star, what, k);
		for(c1 = cexps1 + i; i < j; i++, c1++)
			com_out(c1->e, c1->L, c1->nlin, i + ncom0);
		}
	printf(what ? "\n%s  /***  %s %d  ***/\n\n"
		    : "\n%s  /***  objective ***/\n\n", star, what, k);
	return callb(e,buf);
	}

 static void
#ifdef KR_headers
x0check(i, j, k, deriv)
	int i;
	int j;
	int k;
	int deriv;
#else
x0check(int i, int j, int k, int deriv)
#endif
{
	cexp *c;

	if (!needx0check)
		return;
	printf(deriv ? xcheck : xcheck0);
	if (i < j) {
		printf(xkind, k, k, k);
		for(c = cexps + i; i < j; i++, c++)
			com_out(c->e, c->L, c->nlin, i);
		printf("\t%s", endif_fmt);
		}
	}

 static void
#ifdef KR_headers
funnel_set(f, gj0_fmt)
	funnel *f; char *gj0_fmt;
#else
funnel_set(funnel *f, char *gj0_fmt)
#endif
{
	cplist *cl;
	dLR *d;
	Adjoint *a;
	int k;
	list *Lp;
	char *b, *b1, buf[32], bufb[32];

	stor_grad = STOR_DEFV;
	fwalk = 1;
	do {
		printf("\n%s\t/*** funnel ***/\n\n", Fortstar);
		k = f->ce - cexps;
		dwalk(f->fcde.d, 0, zaC[k], Lp = com_ifset[k]);
		co_derivs(f->fcde.d, Lp, gj0_fmt);
		cl = f->cl;
		do {
			d = dLRp(*cl->cfa); /*!! was Make_dLR(cl->cfa) */
			sprintf(buf, pd_fmt, d->o.i + Fortran);
			a = Adjp(cl->ca.rp);
			b = bufb;
			switch(a->storage) {
			  case STOR_IMPLICIT:
				b = a->neg ? "-1" : "1";
				break;
			  case STOR_VP:
				g_fmt(b, *a->o.vp);
				break;
			  case STOR_PD:
				b1 = b;
				if (a->neg)
					*b1++ = '-';
				sprintf(b1, pd_fmt, a->o.i + Fortran);
				break;
			  case STOR_DV:
			  case STOR_DEFV:
				sprintf(b, dv_fmt, a->o.i);
				break;
			  default:/*DEBUG*/ fprintf(Stderr,
				"Bad a->storage in funnel_set\n");
				/*DEBUG*/exit(15);
				}
			*a = A0;
			assign(buf, b);
			}
			while(cl = cl->next);
		}
		while(f = f->next);
	fwalk = 0;
	}

 static void
end(VOID)
{
	printf(Fortran ? "\tend\n" : "}\n");
	branches = 0;
	}

 void
#ifdef KR_headers
fobj_output(deriv)
	int deriv;
#else
fobj_output(int deriv)
#endif
{
	int *c1, i, lbl;
	char buf[24], *eval;
	ograd *og;
	static char *rv[2] = { "feval0", "feval" };

	sprintf(buf, n_obj ? "%d" : "*", o_vars);
	printf(deriv ?
	"\n\tdouble precision function feval(nobj, needfg, x, g)\n\
	integer nobj, needfg\n\
	double precision x(%s), g(%s)\n"
	: "\n\tdouble precision function feval0(nobj, x)\n\
	integer nobj\n\
	double precision x(%s)\n", buf, buf);

	if (n_obj > 1)
		printf("\tinteger nobj1\n");
	else if (!n_obj) {
		printf("\tend\n");
		return;
		}

	want_derivs = deriv;
	for(og = 0, i = 0; i < n_obj; i++) {
		for(og = Ograd[i]; og; og = og->next)
			if (og->coef)
				goto break2;
		}
 break2:
	cde_dwalk(obj_de, n_obj, &omax, o_ifset, zao, (cgrad**)Ograd);
	if (needx0check)
		printf("%s\tlogical xcheck\n\texternal xcheck\n",
			deriv ? "\tinteger wantfg\n" : "");
	if (omax.nvt) {
		printf("\n\tdouble precision v(%d)", omax.nvt);
		if (omax.ndv > 1 && deriv)
			printf(", dv(%d)", omax.ndv-1);
		printf(og ? ", rv\n" : "\n");
		if (omax.ncond > 1 && want_derivs)
			printf("\tinteger cond(%d)\n", omax.ncond-1);
		}
	else if (og)
		printf("\n\tdouble precision rv\n");

	if (omax.needT1)
		printf("\tdouble precision t1, t2\n");

	plcommon(deriv);

	x0check(combc, ncom0, 2, deriv);

	if ((f_b || f_o) && deriv) {
		printf("\tif (wantfg .ge. 2) then\n");
		if (f_b)
			printf("\t\tcall funelb(x)\n");
		if (f_o)
			funnel_set(f_o, g0_fmt);
		printf("\t\tendif\n");
		}

	if (n_obj > 1) {
		printf("\tnobj1 = nobj + 1\n");
		lbl = Switch("nobj1", n_obj);
		}
	c1 = o_cexp1st;
	for(i = 0;; c1++) {
		if (n_obj > 1)
			Case(lbl+i);
		if (deriv)
			printf("\tif (mod(wantfg,2) .eq. 1) then\n");
		eval = putout(obj_de[i].e, c1[0], c1[1],
				n_obj > 1 ? "objective" : NULL, i, zao[i]);
		rv_output(rv[deriv], eval, Ograd[i]);
		if (deriv) {
			printf("\n\tendif\n\n\tif (wantfg .gt. 1) then\n");
			obj_derivs(i);
			printf("\n\tendif\n");
			}
		if (++i >= n_obj)
			break;
		printf("\treturn\n");
		}
	end();
	}

 void
#ifdef KR_headers
obj_output(deriv)
	int deriv;
#else
obj_output(int deriv)
#endif
{
	int *c1, i;
	ograd *og;
	static char rv[] = "rv";
	static char *header[4] = {
		"0_(fint *nobj, real *x)",
		"_(fint *nobj, fint *needfg, real *x, real *g)",
		"0_(nobj, x)\n\tfint *nobj;\n real *x;",
		"_(nobj, needfg, x, g)\n fint *nobj, *needfg;\n real *x, *g;"
		};
	char *eval, *s;

	printf(" real\nfeval%s\n", header[deriv + krflag]);

	if (!n_obj) {
		/*{*/ printf("{ return 0.; }\n");
		return;
		}
	want_derivs = deriv;
	printf("{");
	for(og = 0, i = 0; i < n_obj; i++) {
		for(og = Ograd[i]; og; og = og->next)
			if (og->coef)
				goto break2;
		}
 break2:
	cde_dwalk(obj_de, n_obj, &omax, o_ifset, zao, (cgrad**)Ograd);
	if (omax.nvt) {
		printf("\n\treal v[%d]", omax.nvt);
		if (omax.ndv)
			printf(", dv[%d]", omax.ndv);
		printf("%s;\n", og ? ", rv" : "");
		}
	else if (og)
		printf("\n\treal rv;\n");

	if (omax.ncond)
		printf("\tstatic int cond[%d];\n", omax.ncond);

	if (omax.needT1)
		printf("\treal t1, t2;\n");

	x0check(combc, ncom0, 2, deriv);

	if ((f_b || f_o) && deriv) {
		printf("\tif (wantfg & 2) {\n");
		if (f_b)
			printf("\t\tfunnelb(x);\n");
		if (f_o)
			funnel_set(f_o, g0_fmt);
		printf("\t\t}\n");
		}

	printf(n_obj > 1 ? "\n\tswitch(*nobj) {\n" : "\n");
	c1 = o_cexp1st;
	for(i = 0; i < n_obj; i++, c1++) {
		if (n_obj > 1)
			printf("\n  case %d:\n", i);
		if (deriv)
			printf("\tif (wantfg & 1) {\n");
		eval = putout(obj_de[i].e, c1[0], c1[1],
			n_obj > 1 ? "objective" : NULL, i, zao[i]);
		s = rv_output(rv, eval, Ograd[i]);
		if (deriv) {
			printf("\t;}\n\n\tif (wantfg & 2) {\n");
			obj_derivs(i);
			printf("\t}\n");
			}
		printf("\n\treturn %s;\n", s);
		}
	if (n_obj > 1)
		printf("\n\t}\n");
	end();
	}

 void
#ifdef KR_headers
fcon_output(deriv)
	int deriv;
#else
fcon_output(int deriv)
#endif
{
	int *c1, i;
	char *s;

	printf(deriv ? "\n\tsubroutine ceval(needfg, x, c, J)\n\
	integer needfg\n\tdouble precision x(%d), c(%d), J(%ld)\n"
	: "\n\tsubroutine ceval0(x, c)\n\
	double precision x(%d), c(%d)\n",
		c_vars ? c_vars : 1, n_con ? n_con : 1, nJ ? nJ : 1);

	if (!n_con) {
		printf("\tend\n");
		return;
		}

	want_derivs = deriv;
	cde_dwalk(con_de, n_con, &cmax, c_ifset, zac, Cgrad);
	if (needx0check)
		printf("%s\tlogical xcheck\n\texternal xcheck\n",
			deriv ? "\tinteger wantfg\n" : "");

	if (cmax.nvt) {
		printf("\n\tdouble precision v(%d)", cmax.nvt);
		printf(cmax.ndv > 1  && deriv ? ", dv(%d)\n" : "\n",
			cmax.ndv-1);
		}
	if (cmax.ncond > 1)
		printf("\tinteger cond(%d)\n", cmax.ncond-1);

	if (cmax.needT1)
		printf("\tdouble precision t1, t2\n");
	else if (nzcgrad())
		printf("\tdouble precision t1\n");

	plcommon(deriv);

	x0check(comb, combc, 1, deriv);

	if (deriv)
		printf("\n\tif (mod(wantfg,2) .gt. 0) then\n");
	c1 = c_cexp1st;
	for(i = 0; i < n_con; i++, c1++) {
		s = putout(con_de[i].e, c1[0], c1[1], "constraint", i, zac[i]);
		printf("\tc(%d) = %s\n", i+1, con_linadd(i,s));
		}
	if (deriv) {
		printf("\tendif\n\n\tif (wantfg .gt. 1) then\n");
		if (f_b)
			printf("\t\tcall funelb(x)\n");
		if (f_c)
			funnel_set(f_c, j0_fmt);
		for(i = 0; i < n_con; i++)
			con_derivs(i);
		printf("\tendif\n");
		djoff = -1;
		}
	end();
	}

 void
#ifdef KR_headers
con_output(deriv)
	int deriv;
#else
con_output(int deriv)
#endif
{
	int *c1, i;
	char *s;
	static char *header[4] = {
		"0_(real *x, real *c)",
		"_(fint *needfg, real *x, real *c, real *J)",
		"0_(x, c)\n real *x, *c;",
		"_(needfg, x, c, J)\n fint *needfg;\n real *x, *c, *J;"
		};
	printf("\n void\nceval%s\n{", header[deriv + krflag]);

	if (!n_con) {
		printf("}\n");
		return;
		} /*}*/

	want_derivs = deriv;
	cde_dwalk(con_de, n_con, &cmax, c_ifset, zac, Cgrad);
	if (cmax.nvt) {
		printf("\n\treal v[%d]", cmax.nvt);
		printf(cmax.ndv && deriv ? ", dv[%d];\n" : ";\n", cmax.ndv);
		}
	if (cmax.ncond)
		printf("\tstatic int cond[%d];\n", cmax.ncond);

	if (cmax.needT1)
		printf("\treal t1, t2;\n");
	else if (nzcgrad())
		printf("\treal t1;\n");

	x0check(comb, combc, 1, deriv);
	if (deriv)
		printf("\n\tif (wantfg & 1) {\n");
	c1 = c_cexp1st;
	for(i = 0; i < n_con; i++, c1++) {
		s = putout(con_de[i].e, c1[0], c1[1], "constraint", i, zac[i]);
		printf("\tc[%d] = %s;\n", i, con_linadd(i,s));
		}
	if (deriv) {
		printf("\t;}\n\tif (wantfg & 2) {\n");
		if (f_b)
			printf("\tfunnelb(x);\n");
		if (f_c)
			funnel_set(f_c, j0_fmt);
		for(i = 0; i < n_con; i++)
			con_derivs(i);
		printf("\t}\n");
		djoff = -1;
		}
	end();
	}

 void
output(VOID)
{
	int i, j;
	plterm *p;
	expr *e;
	real *b, *be;
	dLR *LR;
	char buf[32], *x0;
	cexp *c;
	static char *fhead[2] = { "(real *x)", "(x) real *x;" };

#ifdef DEBUG
	if (awrite | vwrite)
		printf("#include \"stdio.h\"\n");
#endif
	printf("#include \"math.h\"\n#include \"errno.h\"\n\
#ifndef fint\n\
#ifndef Long\n\
#include \"arith.h\"	/* for Long */\n\
#ifndef Long\n\
#define Long long\n\
#endif\n\
#endif\n\
#define fint Long\n\
#endif\n\
#ifndef real\n\
#define real double\n\
#endif\n");
	if (!krflag)
		printf("#ifdef __cplusplus\nextern \"C\" {\n#endif\n");
	if (krflag < 2)
		printf(" %s\n %s\n %s\n %s\n %s\n %s\n %s\n %s\n",
			"real acosh_(real *);",
			"real asinh_(real *);",
			"real acoshd_(real *, real *);",
			"real asinhd_(real *, real *);",
			"void in_trouble(char *, real);",
			"void in_trouble2(char *, real, real);",
			"void domain_(char *, real *, fint);",
			"void zerdiv_(real *);");
	printf(" fint auxcom_[1] = { %d /* nlc */ };\n", nlc);
	printf(" fint funcom_[%d] = {\n\
	%d /* nvar */,\n\
	%d /* nobj */,\n\
	%d /* ncon */,\n\
	%d /* nzc */,\n\
	%d /* densejac */",
		n_con && !densejac ? nzc + nv1 + n_obj + 6 : n_obj + 5,
		nv1, n_obj, n_con, densejac ? nv1*n_con : nzc, densejac);

	if (n_obj) {
		printf(",\n\n\t/* objtype (0 = minimize, 1 = maximize) */\n");
		for(i = 0; i < n_obj; i++)
			printf("%s\n\t%d", i ? "," : "", objtype[i]);
		}

	if (n_con && !densejac) {
		printf(",\n\n\t/* colstarts */\n");
		for(i = 0; i <= n_var; i++)
			printf("\t%d,\n", A_colstarts[i]);
		printf("\n\t/* rownos */\n\t1");
		for(i = 1; i < nzc; i++)
			printf(",\n\t%d", rownos[i]);
		}
	printf(" };\n\n");

	for(i = j = 0; i < N_OPS; i++)
		if (seen[i] && declare[i])
			printf(j++ ? ", %s_()" : " extern real %s_()",
				declare[i]);
	if (j)
		printf(";\n");
	for(e = Plterms; e; e = LR->o.ep) {
		LR = dLRp(e->dR);
		p = e->L.p;
		i = 2*p->n - 1;
		printf(" real bs%d[%d] = {\n", LR->kind, i);
		for(b = p->bs, be = b + i; b < be; b++) {
			g_fmt(buf, *b);
			printf("\t%s%s\n", buf,
				b + 1 < be ? "," : "};\n");
			}
		}

	printf(
" real boundc_[1+%d+%d] /* Infinity, variable bounds, constraint bounds */ = {\n\t\t1.7e308",
		2*nv1, 2*n_con);
	b = bounds;
	be = b + 2*(n_con + nv1);
	while(b < be)
		printf(",\n\t\t%s", fpval(*b++));
	printf("};\n\n");

	printf(" real x0comn_[%d] = {\n", nv1);
	for(i = 0; i < nv1; i++)
		printf("%s\t\t%s", i ? ",\n" : "", fpval(X0[i]));
	printf(" };\n\n");

	if (npd)
		printf(" static real pd[%d];\n", npd);

	if (f_b && derkind & 2) {
		printf("\n static void\nfunnelb%s\n{\n", fhead[krflag>>1]);
		if ((i = bmax.ndv) > 0)
			printf("\treal dv[%d];\n", i);
		funnel_set(f_b, "Botch<%d> = %s;");
		printf("\t}\n");
		}

	if (needx0check) {
		printf("static real old_x[%d];\nstatic int xkind = -1;\n\n\
 static int\nxcheck%s\n{\n\treal", nv1, xcheckdcl);
		if (comb > 0) {
			printf(" *x0 = x,");
			x0 = "x0";
			}
		else
			x0 = "x";
		printf(" *x1 = old_x, *xe = x + %d;\n", nv1);
		if ((i = bmax.nvt) > 0)
			printf("\treal v[%d];\n", i);
		printf("\terrno = 0;\n\
	if (xkind >= 0) {\n\t\twhile(*%s++ == *x1++)\n\
		\tif (%s == xe)\n\t\t\t\treturn 0;\n\t\t--%s, --x1;\n\t\t}\n\
	do *x1++ = *%s++;\n\t\twhile(%s < xe);\n\txkind = 0;\n",
			x0,x0,x0,x0,x0);
		for(i = 0, c = cexps; i < comb; c++, i++)
			com_out(c->e, c->L, c->nlin, i);
		printf("\treturn 1;\n\t}\n");
		}
	for(i = 1; i < 3; i++) {
		if (owant & i)
			obj_output(i-1);
		if (cwant & i)
			con_output(i-1);
		}
	if (!krflag)
		printf("#ifdef __cplusplus\n\t}\n#endif\n");
	}

 void
foutput(VOID)
{
	int colrow, i, j, nb, ncr;
	plterm *p;
	real *b, *be;
	expr *e;
	dLR *LR;
	char buf[32];
	cexp *c;

	ncr = n_obj;
	if (colrow = n_con && !densejac)
		ncr += nzc + n_var + 1;
	printf(
	 "\tblock data\n\tcommon /funcom/ nvar, nobj, ncon, nzc, densej%s\n",
		ncr ? ", colrow" : "");
	printf("\tinteger nvar, nobj, ncon, nzc, densej");
	if (ncr)
		printf(", colrow(%d)", ncr);

	nb = 2*(n_con + o_vars);
	printf("\n\n\tcommon /boundc/ bounds\n\
	double precision bounds(%d)\n\
	common /x0comn/ x0\n\tdouble precision x0(%d)\n",
		nb + 1, nv1);
	printf("\tcommon /auxcom/ nlc\n\tinteger nlc\n");
	plcommon(1);

	printf("\tdata nvar/%d/, nobj/%d/, ncon/%d/, nzc/%d/, densej/%d/\n",
		nv1, n_obj, n_con, colrow ? nzc : n_con*nv1, densejac);
	if (n_obj) {
		printf("\n*\t*** objtype (0 = minimize, 1 = maximize) ***\n");
		for(i = 0; i < n_obj; i++)
			printf("\tdata colrow(%d)/%d/\n", i+1, objtype[i]);
		}
	if (colrow) {
		printf("\n*\t*** colstarts ***\n");
		j = n_obj + 1;
		for(i = 0; i <= n_var; i++)
			printf("\tdata colrow(%d)/%d/\n", i+j, A_colstarts[i]+1);
		j += n_var + 1;
		printf("\n*\t*** rownos ***\n\tdata colrow(%d)/1/\n", j);
		for(i = 1; i < nzc; i++)
			printf("\tdata colrow(%d)/%d/\n", i+j, rownos[i]);
		}
	for(i = 0; i < nv1; i++)
		printf("\tdata x0(%d)/%s/\n", i+1, fpval(X0[i]));
	for(e = Plterms; e; e = LR->o.ep) {
		LR = dLRp(e->dR);
		p = e->L.p;
		i = 2*p->n - 1;
		j = 0;
		for(b = p->bs, be = b + i; b < be; b++) {
			g_fmt(buf, *b);
			printf("\tdata bs%d(%d)/%s/\n", LR->kind, ++j, buf);
			}
		}

	b = bounds;
	be = b + nb;
	i = 2;
	printf("\tdata bounds(1)/1.7e+308/\n");
	while(b < be)
		printf("\tdata bounds(%d)/%s/\n", i++, fpval(*b++));
	printf("\tdata nlc/%d/\n", nlc);
	printf("\tdata xkind/-1/\n\tend\n\n");

	if (f_b && derkind & 2) {
		printf("\tsubroutine funelb(x)\n");
		printf("\tdouble precision x(%d)\n", n_var);
		plcommon(1);
		if ((i = bmax.ndv) > 0)
			printf("\tdouble precision dv(%d)\n", i);
		funnel_set(f_b, "Botch<%d> = %s");
		printf("\tend\n");
		}

	if (needx0check) {
		printf("\tlogical function xcheck(x)\n\
	double precision x(%d)\n\
	double precision oldx(%d)\n", nv1, nv1);
		if ((i = bmax.nvt) > 0)
			printf("\tdouble precision v(%d)\n", i);
		if (comb > 0)
			plcommon(1);
		else
			printf(
		  "\tcommon /xkindc/ xkind\n\tinteger xkind\n\tsave /xkindc/\n");
		printf("\tif (xkind .lt. 0) then\n\t\ti = 1\n\t\tgoto 20\n\t\tendif\n\
	do 10 i = 1, %d\n\t\tif (x(i) .ne. oldx(i)) goto 20\n 10\t\tcontinue\n\
	xcheck = .false.\n\treturn\n 20 \tdo 30 i = i, %d\n\
 30\t\toldx(i) = x(i)\n\txkind = 0\n", nv1, nv1);
		for(i = 0, c = cexps; i < comb; c++, i++)
			com_out(c->e, c->L, c->nlin, i);
		printf("\txcheck = .true.\n\tend\n");
		}

	for(i = 1; i < 3; i++) {
		if (owant & i)
			fobj_output(i-1);
		if (cwant & i)
			fcon_output(i-1);
		}
	}

 static void
#ifdef KR_headers
usage(rc)
	int rc;
#else
usage(int rc)
#endif
{
	static char *op[] = {
		"-1 { produce feval0, ceval0 (function values only) }",
		"-2 { produce feval, ceval (functions and gradients, default) }",
		"-3 { produce feval0, ceval0, feval, ceval }",
		"-4 { produce feval, ceval0 }",
		"-5 { produce feval0, ceval }",
#ifdef DEBUG
		"-a { write adjoint operations on Stderr (DEBUG option) }",
#endif
		"-d { dense Jacobian }",
		"-f { Fortran output }",
		"-i { no derivatives for integer variables }",
		"-k { K&R C output (default = ANSI C) }",
#ifdef DEBUG
		"-s { show statistics on Stderr (DEBUG option) }",
		"-v { write def var assignments on Stderr (DEBUG option) }",
#endif
		0};
	char **o;

	fprintf(Stderr, "usage: %s [options] file[.nl]\noptions:\n", progname);
	for(o = op; *o; o++)
		fprintf(Stderr, "\t%s\n", *o);
	exit(rc);
	}

 static void
#ifdef KR_headers
cant(s1, s2)
	char *s1;
	char *s2;
#else
cant(char *s1, char *s2)
#endif
{
	fprintf(Stderr, "Can't open %s", s1);
	if (s2)
		fprintf(Stderr, " or %s.nl", s2);
	fputc('\n', Stderr);
	exit(1);
	}

 static void
get_rownos(VOID)
{
	int i = n_con, i1, j, j1;
	cgrad *cg;
	memset((char *)rownos, 0, nzc*sizeof(int));
	while(i1 = i)
		for(cg = Cgrad[--i]; cg; cg = cg->next)
			rownos[cg->goff] = i1;
	if (!Fortran)
		for(i = 0; i <= n_var; i++)
			A_colstarts[i]++;
	if (intsk) {
		i1 = j = 0;
		--intsk;
		for(i = 1; i <= n_var; i++) {
			j1 = A_colstarts[i] - 1;
			if (!intsk[i])
				for(; j < j1; j++)
					rownos[i1++] = rownos[j];
			A_colstarts[i] = i1 + 1;
			j = j1;
			}
		++intsk;
		nzc = i1;
		if (nlvbi && (nlvci < nlvc || nlvoi < nlvo)
		 || nlvci && nlvoi < nlvo) {
			for(i = 0; i < n_var; i++)
				A_colstarts[i] = A_colstarts[i+1];
			i = n_con;
			while(--i >= 0)
			  for(cg = Cgrad[i]; cg; cg = cg->next)
			    if (!intsk[cg->varno])
				cg->goff = --A_colstarts[cg->varno] - 1;
			}
		}
	}

 static void
#ifdef KR_headers
nlvzap(i, j)
	int i;
	int j;
#else
nlvzap(int i, int j)
#endif
{
	memset(intsk + i - j, 1, j);
	}

 int
#ifdef KR_headers
main(argc, argv)
	int argc;
	char **argv;
#else
main(int argc, char **argv)
#endif
{
	FILE *nl;
	char *s, *s0;
	int i, ncom, nv;
	fint L;
	expr_nx *enx;

	ASL_alloc(ASL_read_fg);
	progname = *argv;
	g_fmt_decpt = 1;
 nextarg:
	while ((s = *++argv) && *s == '-')
		for(s0 = s, --argc;;)
		    switch(*++s) {
			default:
				fprintf(Stderr, "invalid option %s\n",s0);
				usage(1);
			case 0:
				goto nextarg;
			case '1':
				want_derivs = 0;
				cwant = owant = derkind = 1;
				break;
			case '2':
			case '3':
				want_derivs = 1;
				cwant = owant = derkind = *s - '0';
				break;
			case '4':
				want_derivs = cwant = 1;
				derkind = owant = 2;
				break;
			case '5':
				want_derivs = owant = 1;
				derkind = cwant = 2;
				break;
#ifdef DEBUG
			case 'a':		/*debug*/
				awrite = 1;	/*debug*/
				break;
			case 's':		/*debug*/
				stats = 1;	/*debug*/
				break;
			case 'v':		/*debug*/
				vwrite = 1;	/*debug*/
				break;
#endif
			case 'd':
				densejac = 1;
				break;
			case 'f':
				Fortran = 1;
				g_fmt_decpt = 2;
				g_fmt_E = 'd';
				Half = ".5d0";
				One = "1d0";
				Zero = "0.d+00";
				Negone = "-1d0";
				opEQ = ".eq.";
				opGE = ".ge.";
				opGT = ".gt.";
				opLE = ".le.";
				opLT = ".lt.";
				opNE = ".ne.";
				opAND = ".and.";
				opOR = ".or.";
				offlfmt1 = "%s";
				offlfmt2 = "%s, %s";
				assign_fmt = "\t%s = %s\n";
				binop_fmt = "\t%s = %s %s %s\n";
				goto_fmt = "\tgo to %d\n";
				ifgo_fmt = "\tif (%s %s %s) go to %d\n";
				label_fmt = "%5d continue\n";
				ifstart_fmt = "\tif (%s %s %s) then\n";
				elseif_fmt = "\telse if (%s %s %s) then\n";
				else_fmt = "\telse\n";
				endif_fmt = "\tendif\n";
				call_fmt = "\t%s(%s)\n";
				case_fmt = label_fmt;
				cond_fmt = "cond(%d)";
				break_fmt = goto_fmt;
				endswitch_fmt = label_fmt;
				zerdiv_fmt = "\tcall zerdiv(%s)\n";
				dv_fmt = "dv(%d)";
				Fortstar = "*";
				pd_fmt = "pd(%d)";
				tv_fmt = "v(%d)";
				x_fmt = "x(%d)";
				grad_fmt = "g(%d)";
				jac_fmt = "J(%d)";
				g0_fmt = "\tg(%d) = %s\n";
				j0_fmt = "\tJ(%d) = %s\n";
				star = "*";
				eos = "\n";
				xcheck = "\twantfg = needfg\n\
	if (xcheck(x) .and. wantfg .eq. 2) wantfg = 3\n";
				xcheck0 = "\txchk = xcheck(x)\n";
				xkind =
		"\tif (mod(xkind,2*%d) .lt. %d) then\n\t\txkind = xkind + %d\n";
				break;
			case 'i':
				skip_int_derivs = 1;
				break;
			case 'k':
				krflag = 2;
				Void = "";
				xcheckdcl = "(x) real *x;";
				break;
			case '?':
				usage(0);
			}
	if (argc > 2 || !s)
		usage(1);
	return_nofile = 1;
	nl = jacdim0(s, L = strlen(s));
	if (!nl) {
		if ((L -= 3) > 0 && !strcmp(s + L, ".nl")) {
			s[L] = 0;
			nl = jacdim0(s, L);
			if (!nl)
				cant(s, s);
			}
		else
			cant(s, 0);
		}
	for(i = 0; i < N_OPS; i++)
		r_ops[i] = (efunc *)i;

	nv1 = c_vars > o_vars ? c_vars : o_vars;
	ncom = (i = comb + comc + como) + comc1 + como1;
	nv2 = nv1 + ncom;

	c_cexp1st = (int *)Malloc((n_con + n_obj + 2)*sizeof(int));
	o_cexp1st = c_cexp1st + n_con + 1;
	zac = (int **)Malloc((n_con + n_obj + i)*sizeof(int*));
	zao = zac + n_con;
	zaC = zao + n_obj;

	if (n_con)
		if (densejac) {
			rownos = (int *)Malloc(nv1*sizeof(int));
			memset((char *)rownos, 0, nv1*sizeof(int));
			}
		else {
			rownos = (int *)Malloc((nzc + nv1 + 1)*sizeof(int));
			A_colstarts = rownos + nzc;
			}

	LUv = bounds = (real *)Malloc((3*nv1+2*n_con)*sizeof(real));
	LUrhs = LUv + 2*nv1;
	X0 = LUrhs + 2*n_con;

	size_expr_n = sizeof(expr_nx);
	fg_read(nl,0);

	needx0check = comb > 0 || derkind & 2;
	c_ifset = (list **)Malloc((n_con + n_obj + ncom0)*sizeof(list *));
	o_ifset = c_ifset + n_con;
	com_ifset = o_ifset + n_obj;
	op_type[OP1POW] = 2;
	op_type[OP2POW] = 11;

	stog = (ograd **)Malloc(nv1*sizeof(ograd *));

	declare[OP_asinh] = "asinh";
	declare[OP_asinh+1] = "asinhd";
	declare[OP_acosh] = "acosh";
	declare[OP_acosh+1] = "acoshd";
	declare[OPPLTERM] = "plterm";

	dvtfree = (int *)Malloc(NDVTGULP*sizeof(int));
	ndvtmax = nvtfree = NDVTGULP;

	Fortran1 = Fortran - 1;

	if (skip_int_derivs) {
		intsk = (char *)Malloc(nv1);
		memset(intsk, 0, nv1);
		if (nlvbi)
			nlvzap(nlvb, nlvbi);
		if (nlvci)
			nlvzap(nlvb+nlvc, nlvci);
		if (nlvoi)
			nlvzap(nlvb+nlvc+nlvo, nlvoi);
		}
	if (n_con && !densejac)
		get_rownos();

	vseen = (int *)Malloc((nv2 + ncom)*sizeof(int));
	cvmap = vseen + nv2;
	npd = 0;
	for(i = 0; i < ncom0; i++)
		cvmap[i] = -(++npd);
	if (ncom1)
		memset((char *)&cvmap[ncom0], 0, ncom1*sizeof(int));

	ndv = ncond = Fortran;

	dvset(&edagread_one, 1);
#ifdef X64_bit_pointers
	{Adjoint *Ap = (Adjoint *)Malloc(amax*sizeof(Adjoint));
	memset((char *)Ap, 0, amax*sizeof(Adjoint));
	for(i = 0; i < amax; i++)
		*(Adjoint **)&adjoints[i] = Ap++;
	}
#else
	memset((char *)adjoints, 0, amax*sizeof(real));
#endif
	rdwalk_one = &adjoints[nv1];
	dwalk_one = Adjp(rdwalk_one);
	if ((i = amax - nv1) > 0) {
		acount = (int *)Malloc(i*sizeof(int));
		memset((char *)acount, 0, i*sizeof(int));
		}

	comwalk(0,comb);
	if (f_b)
		funnelwalk(f_b);
	max_save(&bmax);
	comwalk(comb, combc);
	if (f_c)
		funnelwalk(f_c);
	cde_walk(con_de, n_con, &cmax, c_ifset, c_cexp1st, zac);
	max_restore(&bmax);
	comwalk(combc, ncom0);
	if (f_o)
		funnelwalk(f_o);
	cde_walk(obj_de, n_obj, &omax, o_ifset, o_cexp1st, zao);

	nv = nv1 + ncom;
	for(i = 0; i < nv; i++)
		var_e[i].op = (efunc *)f_OPVARVAL1;

	for(enx = nums; enx; enx = enx->next)
		enx->op = f_OPNUM1;

	if (n_obj || n_con) {
		output_time = 1;
		if (Fortran)
			foutput();
		else
			output();
		}
#ifdef DEBUG
	if (stats) {
		*stub_end = 0;
		fprintf(Stderr, "%s: derps_in %d derps_out %d nderps %d\n",
			filename, derps_in, derps_out, nderps);
		}
#endif
	return 0;
	}

 char *
#ifdef KR_headers
e_val(e, buf)
	expr *e;
	char *buf;
#else
e_val(expr *e, char *buf)
#endif
{
	int i;
	if ((i = e->a) >= 0)
		sprintf(buf, tv_fmt, Fortran1 + i);
	else
		sprintf(buf, pd_fmt, Fortran1 - i);
	return buf;
	}
