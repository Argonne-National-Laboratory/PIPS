/****************************************************************
Copyright (C) 2009 AMPL Optimization LLC
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

AMPL Optimization LLC DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS.  IN NO EVENT SHALL AMPL Optimization LLC OR ANY OF ITS
ENTITIES BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES
OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
SOFTWARE.
****************************************************************/
#include "gurobi_c.h"
#include "getstub.h"
#include "signal.h"
#include "time.h"

#ifndef Sig_ret_type
#define Sig_ret_type void
#define SigRet /*nothing*/
#endif

#ifndef SigRet
#define SigRet /*nothing*/
#endif

#ifndef Sig_func_type
typedef void sig_func_type(int);
#endif

 typedef struct
Dims {
	double	*c;
	double	*x;
	double	*y;
	double	*y0;
	int	*cstat;
	int	*rstat;
	SufDesc	*csd;
	SufDesc *rsd;
	char	*mb, *mbend;
	int	kiv;
	int	missing;
	int	nc0, nv0;
	int	objsense;
	} Dims;

 typedef struct
Filename {
	struct Filename *next;
	char *name;
	} Filename;

 typedef struct
Ext_info {
	char *ext;
	int aftersol;
	} Ext_info;

 typedef struct
mint_values {
	int L;
	int U;
	int val;
	} mint_values;

 enum { /* sf_mint f values */
	set_iis		= 0,
	set_relax	= 1,
	set_mipstval	= 2,
	set_objno	= 3,
	set_sos		= 4,
	set_sos2	= 5,
	set_timing	= 6,
	set_basis	= 7,
	set_intstart	= 8,
	set_outlev	= 9,
	set_bestbound	= 10,
	set_solnsens	= 11,
	set_retmipgap	= 12,
	set_rays	= 13
	};

 static mint_values
mint_val[14] = {
	/* set_iis */		{0, 1, 0},
	/* set_relax */		{0, 1, 0},
	/* set_mipstval */	{0, 1, 1},
	/* set_objno */		{0, 0/*n_obj*/,	1},
	/* set_sos */		{0, 1, 1},
	/* set_sos2 */		{0, 1, 1},
	/* set_timing */	{0, 3, 0},
	/* set_basis */		{0, 3, 3},
	/* set_intstart */	{0, 1, 1},
	/* set_outlev */	{0, 1, 0},
	/* set_bestbound */	{0, 1, 0},
	/* set_solnsens */	{0, 1, 0},
	/* set_retmipgap */	{0, 7, 0},
	/* set_rays */		{0, 3, 3}
	};

#define want_iis	mint_val[0].val
#define relax		mint_val[1].val
#define mipstval	mint_val[2].val
#define nobjno		mint_val[3].U
#define objno		mint_val[3].val
#define sos		mint_val[4].val
#define sos2		mint_val[5].val
#define time_flag	mint_val[6].val
#define basis		mint_val[7].val
#define intstart	mint_val[8].val
#define outlev		mint_val[9].val
#define bestbound	mint_val[10].val
#define solnsens	mint_val[11].val
#define retmipgap	mint_val[12].val
#define rays		mint_val[13].val

 static Filename *Wflist, *Wflist1;
 static GRBmodel *grbmodel;
 static char *logfile, verbuf[64];
 static double Times[5];
 static int breaking, wantlog;
 static jmp_buf Jb;

#if GRB_VERSION_MAJOR >= 3
 static double ams_eps, ams_epsabs;
 static int ams_limit;
 static char *ams_stub;
#endif

 static void
badretfmt(int rc, char *fmt, ...)
{
	ASL *asl = cur_ASL;
	va_list ap;
	char buf[8192], *s;
	int k;

	va_start(ap, fmt);
	k = Vsnprintf(buf, sizeof(buf)-1, fmt, ap) + 1;
	if (rc) {
		solve_result_num = rc;
		memcpy(s = (char*)M1alloc(k), buf, k);
		asl->i.uinfo = s;
		}
	if (!rc)
		fprintf(Stderr, "%s\n", buf);
	va_end(ap);
	}

 static void
failed(GRBenv *env, const char *what)
{
	const char *s = GRBgeterrormsg(env);
	if (s)
		badretfmt(501, "%s failed:\n\t%s.\n", what, s);
	else
		badretfmt(501, "%s failed.\n", what);
	longjmp(Jb,1);
	}

 static void enamefailed(GRBenv *env, const char *what, const char *name);

 static void
namefailed(const char *what, const char *name)
{
	badretfmt(506, "%s(\"%s\") failed.", what, name);
	longjmp(Jb,1);
	}

 static void
badival(Option_Info *oi, keyword *kw, int t, int L, int U)
{
	printf("rejecting %s %d; must be between %d and %d\n",
		kw->name, t, L, U);
	badopt_ASL(oi);
	}

 static char *
sf_mint(Option_Info *oi, keyword *kw, char *v)
{
	int t;
	char *rv;
	int i = Intcast kw->info;
	mint_values *m = mint_val + i;

	if (*v == '?' && v[1] <= ' ') {
		printf("%s=%d\n", kw->name, m->val);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	t = (int)strtol(v, &rv, 10);
	if (rv == v) {
		printf("Expected an integer value for %s, not \"%s\"\n",
			kw->name, v);
		badopt_ASL(oi);
		return v;
		}
	if (t < m->L || t > m->U) {
		badival(oi,kw,t,m->L,m->U);
		return rv;
		}
	m->val = t;
	return rv;
	}

 static char *
sf_dpar(Option_Info *oi, keyword *kw, char *v)
{
	GRBenv *env;
	double p[4], t;
	char *parname, *rv;

	env = (GRBenv*)oi->uinfo;
	parname = (char*)kw->info;

	if (*v == '?' && v[1] <= ' ') {
		if (GRBgetdblparam(env, parname, &t))
			namefailed("GRBgetdblparam", parname);
		printf("%s=%.g\n", kw->name, t);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	t = strtod(v, &rv);
	if (rv == v) {
		printf("Expected a numeric value for %s, not \"%s\"\n",
			kw->name, v);
		badopt_ASL(oi);
		return v;
		}
	if (GRBsetdblparam(env, parname, t)) {
		if (GRBgetdblparaminfo(env, parname, p, p+1, p+2, p+3))
			namefailed("GRBsetdblparam", parname);
		badretfmt(506, "%s must be >= %.g and <= %.g.", kw->name, p[1], p[2]);
		badopt_ASL(oi);
		}
	return rv;
	}

 static void
int_rangerr(Option_Info *oi, keyword *kw)
{
	GRBenv *env;
	char *fmt, *parname;
	int p[4];

	env = (GRBenv*)oi->uinfo;
	parname = (char*)kw->info;
	if (GRBgetintparaminfo(env, parname, p, p+1, p+2, p+3))
		namefailed("GRBsetintparam", parname);
	fmt = p[2] == p[1] + 1
		? "%s must be %d or %d."
		: "%s must be >= %d and <= %d.";
	badretfmt(506, fmt, kw->name, p[1], p[2]);
	badopt_ASL(oi);
	}

 static char *
sf_ipar(Option_Info *oi, keyword *kw, char *v)
{
	GRBenv *env;
	int t;
	char *parname, *rv;

	env = (GRBenv*)oi->uinfo;
	parname = (char*)kw->info;

	if (*v == '?' && v[1] <= ' ') {
		if (GRBgetintparam(env, parname, &t))
			namefailed("GRBgetintparam", parname);
		printf("%s=%d\n", kw->name, t);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	t = (int)strtol(v, &rv, 10);
	if (rv == v) {
		printf("Expected an integer value for %s, not \"%s\"\n",
			kw->name, v);
		badopt_ASL(oi);
		return v;
		}
	if (GRBsetintparam(env, parname, t))
		int_rangerr(oi, kw);
	return rv;
	}

 static char*
sf_iparlog(Option_Info *oi, keyword *kw, char *v)
{
	GRBenv *env;
	int t;
	char *parname, *rv;

	env = (GRBenv*)oi->uinfo;
	parname = (char*)kw->info;

	if (*v == '?' && v[1] <= ' ') {
		if (GRBgetintparam(env, parname, &t))
			namefailed("GRBgetintparam", parname);
		printf("%s=%d\n", kw->name, t);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	t = (int)strtol(v, &rv, 10);
	if (rv == v) {
		printf("Expected an integer value for %s, not \"%s\"\n",
			kw->name, v);
		badopt_ASL(oi);
		return v;
		}
	if (GRBsetintparam(env, parname, t))
		int_rangerr(oi, kw);
	else if (t)
		++wantlog;
	return rv;
	}

 static Filename *
fn_value(char **pv, const char *what, keyword *kw)
{
	ASL *asl = cur_ASL;
	Filename *f;
	char *s, *t, *v;
	int c, q;
	size_t L;

	v = *pv;
	q = *v;
	if (q == '"' || q == '\'') {
		s = ++v;
		for(;;) {
			if (!(c = *s)) {
				printf("Bad %s \"%s\" for %s\n", what, v, kw->name);
				*pv = v;
				return 0;
				}
			++s;
			if (c == q && *s != q)
				break;
			}
		}
	else {
		q = 0;
		for(s = v; *s > ' '; ++s);
		if (s == v) {
			printf("Missing %s for %s\n", what, kw->name);
			return 0;
			}
		}
	L = s - v;
	f = M1alloc(sizeof(Filename) + L + 1);
	f->name = t = (char*)(f + 1);
	if (q) {
		for(s = v;; ++s, ++t) {
			if ((*t = *s) == q) {
				if (*++s == q)
					continue;
				break;
				}
			}
		}
	else {
		memcpy(t, v, L);
		t += L;
		}
	*t = 0;
	*pv = s;
	return f;
	}

#if GRB_VERSION_MAJOR > 1 /*{*/
#define GRB_MAJ2(x) x

 static char *
sf_spar(Option_Info *oi, keyword *kw, char *v)
{
	Filename *f;
	GRBenv *env;
	char *parname, tbuf[GRB_MAX_STRLEN + 8];

	env = (GRBenv*)oi->uinfo;
	parname = (char*)kw->info;

	if (*v == '?' && v[1] <= ' ') {
		memset(tbuf, 0, sizeof(tbuf));
		if (GRBgetstrparam(env, parname, tbuf))
			namefailed("GRBgetstrparam", parname);
		printf("%s=\"%s\"\n", kw->name, tbuf);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	if ((f = fn_value(&v, "value", kw))
	 && GRBsetstrparam(env, parname, f->name))
		enamefailed(env, "GRBsetstrparam", parname);
	return v;
	}
#else
#define GRB_MAJ2(x) /*nothing*/
#endif /*}*/

 static char *
sf_wfile(Option_Info *oi, keyword *kw, char *v)
{
	Ext_info *e;
	Filename *f, **pf;
	GRBenv *env;
	char *dot, *t;
	int q;
	static Ext_info W_ext[] = {
		{"bas",1},
		{"lp",0},
		{"mps",0},
		{"prm",0},
		{"sol",1},
		{0,0}};

	env = (GRBenv*)oi->uinfo;

	q = *v;
	if (q == '?' && v[1] <= ' ') {
		for(f = Wflist; f; f = f->next)
			printf("%s=\"%s\"\n", kw->name, f->name);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	if (!(f = fn_value(&v, "file name", kw))) {
		printf("Bad file name \"%s\" for %s\n", v, kw->name);
		badopt_ASL(oi);
		return v;
		}
	dot = 0;
	for(t = f->name; *t; )
		if (*t++ == '.')
			dot = t;
	if (dot)
		for(e = W_ext; e->ext; ++e)
			if (!strcmp(e->ext, dot))
				goto good_ext;
	printf("File name for %s must end in one of\n", kw->name);
	for(e = W_ext; e->ext; ++e)
		printf("\t.%s\n", e->ext);
	badopt_ASL(oi);
	goto ret;
 good_ext:
	pf = e->aftersol ? &Wflist1 : &Wflist;
	f->next = *pf;
	*pf = f;
 ret:
	return v;
	}

#if GRB_VERSION_MAJOR >= 3 /*{*/

 static char *
sf_pf(Option_Info *oi, keyword *kw, char *v)
{
	FILE *f;
	Filename *fn;
	GRBenv *env;
	char buf[512], *fname, *s, *s1, *se;
	int lineno;
	static char extra[]  = "Line %d of paramfile \"%s\":\n\t"
		"expected a name and value, but got \"%s\".";
	static char failed[] = "Line %d of paramfile \"%s\":\n\t"
		"GRBsetstrparam(\"Dummy\", \"%s\") failed:\n\t%s.";
	static char missing[] = "Missing value in line %d of paramfile \"%s\".";

	env = (GRBenv*)oi->uinfo;

	if ((fn = fn_value(&v, "value", kw))) {
		fname = fn->name;
		if (!(f = fopen(fname, "r"))) {
			badretfmt(511, "Cannot open paramfile \"%s\".", fname);
			longjmp(Jb,1);
			}
		lineno = 0;
 nextline:
		while(fgets(buf, sizeof(buf), f)) {
			++lineno;
			for(s = buf; *s <= ' '; ++s)
				if (!*s)
					goto nextline;
			if (*s == '#')
				goto nextline;
			for(s1 = s; *++s1 > ' '; );
			while(*s1 <= ' ')
				if (!*s1++) {
					badretfmt(512, missing, lineno, fname);
					longjmp(Jb,1);
					}
			for(se = s1; *++se; );
			while(--se > s1 && *se <= ' ');
			se[1] = 0;
			while(*++s1 > ' ');
			if (*s1) {
				for(se = s1; *++se; ) {
					if (*se > ' ') {
						badretfmt(513, extra, lineno, fname, s);
						longjmp(Jb,1);
						}
					}
				*s1 = 0;
				}
			if (GRBsetstrparam(env, "Dummy", s)) {
				badretfmt(514, failed, lineno, fname, s, GRBgeterrormsg(env));
				longjmp(Jb,1);
				}
			}
		fclose(f);
		}
	return v;
	}

 static char aggfill_desc[] = "amount of fill allowed during aggregation during\n\
			gurobi's presolve (default 10)";
#endif /*}*/

 static char aggregate_desc[] = "whether to use aggregation during Gurobi presolve:\n\
			0 = no (sometimes reduces numerical errors)\n\
			1 = yes (default)";

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char ams_eps_desc[] = "relative tolerance for reporting alternate MIP solutions\n\
			(default = no limit)";

 static char ams_epsabs_desc[] = "absolute tolerance for reporting alternate MIP solutions\n\
			(default = no limit)";

 static char ams_limit_desc[] = "limit on number of alternate MIP solutions written\n\
			(default = number of available alternate solutions)";

 static char ams_stub_desc[] = "stub for alternate MIP solutions.  The number of\n\
			alternative MIP solution files written is determined\n\
			by three keywords:\n\
			  ams_limit gives the maximum number of files written;\n\
			  ams_eps gives a relative tolerance on the objective\n\
				values of alternative solutions; and\n\
			  ams_epsabs gives an absolute tolerance on how much\n\
				worse the objectives can be.";

 static char barconvtol_desc[] = "tolerance on the relative difference between the\n\
			primal and dual objectives for stopping the barrier\n\
			algorithm (default 1e-8)";

 static char barcorrectors_desc[] = "Limit on the number of central corrections done in\n\
			each barrier iteration (default -1 = automatic choice)";

 static char bariterlimit_desc[] = "Limit on the number of barrier iterations (default none)";

 static char barorder_desc[] = "Ordering used to reduce fill in sparse-matrix factorizations\n\
				during the barrier algorithm:\n\
		       -1 = automatic choice\n\
			0 = approximate minimum degree\n\
			1 = nested dissection";
#endif /*}*/

 static char basis_desc[] = "whether to use or return a basis:\n\
			0 = no\n\
			1 = use incoming basis (if provided)\n\
			2 = return final basis\n\
			3 = both (1 + 2 = default)";

 static char bestbound_desc[] = "whether to return suffix .bestbound for the\n\
		best known bound on the objective value:\n\
			0 = no (default)\n\
			1 = yes";

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char crossover_desc[] = "how to transform a barrier solution to a basic one:\n\
		       -1 = automatic choice (default)\n\
			0 = none: return an interior solution\n\
			1 = push dual vars first, finish with primal simplex\n\
			2 = push dual vars first, finish with dual simplex\n\
			3 = push primal vars first, finish with primal simplex\n\
			4 = push primal vars first, finish with dual simplex";
 static char crossoverbasis_desc[] = "strategy for initial basis construction during crossover:\n\
			0 = favor speed (default)\n\
			1 = favor numerical stability";
#endif /*}*/

#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
 static char branchdir_desc[] = "Which child node to explore first when branching:\n\
			-1 = explore \"down\" branch first\n\
			 0 = explore \"most promising\" branch first (default)\n\
			 1 = explore \"up\" branch first";
#endif /*}*/
 static char cutagg_desc[] = "maximum number of constraint aggregation passes\n\
		during cut generation (-1 = default = no limit);\n\
		overrides \"cuts\"";

 static char cutoff_desc[] = "target objective value:  stop searching once an objective\n\
		value better than the target is found; default:\n\
		-Infinity for minimizing, +Infinity for maximizing";

#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
 static char cutpasses_desc[] = "maximum number of cutting-plane passes to do\n\
		during root-cut generation; default = -1 ==> automatic choice";
#endif /*}*/

 static char cuts_desc[] = "global cut generation control, valid unless overridden\n\
		by individual cut-type controls:\n\
		       -1 = automatic choice (default)\n\
			0 = no cuts\n\
			1 = conservative cut generation\n\
			2 = aggressive cut generation"
			GRB_MAJ2("\n\t\t\t3 = very aggressive cut generation")
		;

 static char feastol_desc[] = "primal feasibility tolerance (default 1e-6)";

 static char gomory_desc[] = "maximum number of Gomory cut passes during cut generation\n\
		(-1 = default = no limit); overrides \"cuts\"";

 static char heurfrac_desc[] = "fraction of time to spend in MIP heuristics (default 0.05)";

 static char iisfind_desc[] = "whether to return an IIS (via suffix .iis) when\n\
		the problem is infeasible:\n\
			0 = no (default)\n\
			1 ==> yes";

#if GRB_VERSION_MAJOR > 1 /*{*/
 static char iismethod_desc[] = "which method to use when finding an IIS (irreducible\n\
		infeasible set of constraints, including variable bounds):\n\
		       -1 = automatic choice (default)\n\
			0 = often faster than method 1\n\
			1 = can find a smaller IIS than method 0";
#endif /*}*/

#if GRB_VERSION_MAJOR >= 4 /*{*/
 static char isg_desc[] = "optimality gap below which the MIP solver switches from\n\
		trying to improve the best bound to trying to find better\n\
		feasible solutions (default 0)";

 static char ist_desc[] = "execution seconds after which the MIP solver switches from\n\
		trying to improve the best bound to trying to find better\n\
		feasible solutions (default Infinity)";
#endif /*}*/

 static char intfeastol_desc[] = "feasibility tolerance for integer variables (default 1e-05)";

 static char intstart_desc[] = "when there are integer variables, whether to use\n\
		an initial guess (if available):\n\
			0 = no\n\
			1 = yes (default)";

 static char logfile_desc[] = "name of file to receive log lines (default: none)"
		GRB_MAJ2(";\n\t\t\timplies outlev = 1")
		;

 static char logfreq_desc[] = "interval in seconds between log lines (default 5)";

 static char maxmipsub_desc[] = "maximum number of nodes for RIMS heuristic to explore\n\
		on MIP problems (default 500)";

#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
 static char minrelnodes_desc[] = "number of nodes for the Minimum Relaxation heuristic\n\
		to explore at the MIP root node when a feasible solution has\n\
		not been found by any other heuristic (default 0)";
#endif /*}*/
#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char mipfocus_desc[]  = "MIP solution strategy:\n\
			0 = balance finding good feasible solutions and\n\
			    proving optimality (default)\n\
			1 = favor finding feasible solutions\n\
			2 = favor proving optimality\n\
			3 = focus on improving the best objective bound";
#endif /*}*/

 static char mipstart_desc[] = "whether to use initial guesses in problems with\n\
		integer variables:\n\
			0 = no\n\
			1 = yes (default)";

#if GRB_VERSION_MAJOR >= 4 /*{*/
 static char nodemethod_desc[] = "algorithm used to solve relaxed MIP node problems:\n\
			0 = primal simplex\n\
			1 = dual simplex (default)\n\
			2 = barrier";
#endif /*}*/


 static char objno_desc[] = "objective to optimize:\n\
			0 = none\n\
			1 = first (default, if available),\n\
			2 = second (if available), etc.";

#if GRB_VERSION_MAJOR > 1 /*{*/
 static char objscale_desc[] = "how to scale the objective:\n\
			0 ==> automatic choice (default)\n\
			negative >= -1 ==> divide by max abs. coefficient\n\
					   raised to this power\n\
			positive ==> divide by this value";
#endif /*}*/

 static char opttol_desc[] = "optimality tolerance on reduced costs (default 1e-6)";

 static char outlev_desc[] = "whether to write Gurobi log lines (chatter) to stdout:\n\
			0 = no (default)\n\
			1 = yes (see logfreq)";

#define Overrides_cuts "overrides \"cuts\"; choices as for \"cuts\""
 static char overrides_cuts[] = Overrides_cuts;

 static char perturb_desc[] = "magnitude of simplex perturbation (when needed; default 2e-4)";

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char param_desc[] = "general way to specify values of both documented and\n\
		undocumented Gurobi parameters; value should be a quoted string\n\
		(delimited by ' or \") containing a parameter name, a space, and\n\
		the value to be assigned to the parameter.  Can appear more\n\
		than once.  Cannot be used to query current parameter values.";

 static char paramfile_desc[] = "name of file (surrounded by 'single' or \"double\" quotes if the\n\
		name contains blanks) of parameter names and values for them.\n\
		Lines that start with # are ignored.  Otherwise, each nonempty\n\
		line should contain a name and a value, separated by a space.";

 static char predeprow_desc[] = "whether Gurobi's presolve should remove linearly\n\
		dependent constraint-matrix rows:\n\
		       -1 = only for continuous models\n\
			0 = never\n\
			1 = for all models";

 static char predual_desc[] = "whether gurobi's presolve should form the dual of a\n\
				continuous model:\n\
		       -1 = automatic choice (default)\n\
			0 = no\n\
			1 = yes\n\
			2 = form both primal and dual and use two threads to\n\
			     choose heuristically between them";
#endif /*}*/

 static char presolve_desc[] = "whether to use Gurobi's presolve:\n\
		       -1 (default) = automatic choice\n\
			0 = no\n\
			1 = conservative presolve\n\
			2 = aggressive presolve";

 static char pricing_desc[] = "pricing strategy:\n\
		       -1 = automatic choice (default)\n\
			0 = partial pricing\n\
			1 = steepest edge\n\
			2 = Devex\n\
			3 = quick-start steepest edge";

#if GRB_VERSION_MAJOR >= 4 /*{*/
 static char psdtol_desc[] = "maximum diagonal perturbation to correct indefiniteness\n\
		in quadratic objectives (default 1e-6)";
#endif /*}*/

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char pumppasses_desc[] = "number of feasibility-pump passes to do after the\n\
			MIP root when no other root heuristoc found a\n\
			feasible solution (default 0)";
#endif /*}*/

#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
 static char rays_desc[] = "Whether to return suffix .unbdd when the objective is unbounded\n\
		or suffix .dunbdd when the constraints are infeasible:\n\
			0 = neither\n\
			1 = just .unbdd\n\
			2 = just .dunbdd\n\
			3 = both (default)";
#endif /*}*/

 static char relax_desc[] = "whether to enforce integrality:\n\
			0 = yes (default)\n\
			1 = no: treat integer and binary variables\n\
				as continuous";

 static char return_mipgap_desc[] =
		"Whether to return mipgap suffixes or include mipgap values\n\
		(|objectve - best_bound|) in the solve_message:  sum of\n\
			1 = return relmipgap suffix (relative to |obj|);\n\
			2 = return absmipgap suffix (absolute mipgap);\n\
			4 = suppress mipgap values in solve_message.\n\
		Default = 0.  The suffixes are on the objective and problem.\n\
		Returned suffix values are +Infinity if no integer-feasible\n\
		solution has been found, in which case no mipgap values are\n\
		reported in the solve_message.";

 static char scale_desc[] = "whether to scale the problem:\n\
			0 = no\n\
			1 = yes (default)";

 static char simplex_desc[] = "which algorithm to use:"
#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
			"\n\
			-1 automatic (default): 3 for LP, 2 for QP, 1 for MIP\n\
				root node"
#endif /*}*/
			"\n\
			0 = primal simplex\n\
			1 = dual simplex (default)"
#if GRB_VERSION_MAJOR >= 3 /*{*/
			"\n\
			2 = barrier"
#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5
			"\n\
			3 = nondeterministic concurrent (several solves in\n\
				parallel)\n\
			4 = deterministic concurrent"
#endif
#endif /*}*/
			;

 static char solnsens_desc[] = "whether to return suffixes for solution sensitivities, i.e.,\n\
		ranges of values for which the optimal basis remains optimal:\n\
			0 = no (default)\n\
			1 = yes:  suffixes return on variables are\n\
				.sensobjlo = smallest objective coefficient\n\
				.sensobjhi = greatest objective coefficient\n\
				.senslblo = smallest variable lower bound\n\
				.senslbhi = greatest variable lower bound\n\
				.sensublo = smallest variable upper bound\n\
				.sensubhi = greatest variable upper bound\n\
			suffixes for constraints are\n\
				.sensrhslo = smallest right-hand side value\n\
				.sensrhshi = greatest right-hand side value";

 static char sos_desc[] = "whether to honor declared suffixes .sosno and .ref describing\n\
		SOS sets:\n\
			0 = no\n\
			1 = yes (default):  each distinct nonzero .sosno\n\
				value designates an SOS set, of type 1 for\n\
				positive .sosno values and of type 2 for\n\
				negative values.  The .ref suffix contains\n\
				corresponding reference values.";

 static char sos2_desc[] = "whether to tell GUROBI about SOS2 constraints for nonconvex\n\
		piecewise-linear terms:\n\
			1 = no\n\
			2 = yes (default), using suffixes .sos and .sosref\n\
				provided by AMPL.";

 static char threads_desc[] = "maximum threads to use on MIP problems\n\
		(default 0 ==> max possible)";

 static char timing_desc[] = "whether to report timing:\n\
			0 (default) = no\n\
			1 = report times on stdout\n\
			2 = report times on stderr";

 static char varbranch_desc[] = "MIP branch variable selection strategy:\n\
		       -1 = automatic choice (default)\n\
			0 = pseudo reduced-cost branching\n\
			1 = pseudo shadow-price branching\n\
			2 = maximum infeasibility branching\n\
			3 = strong branching";

 static char writeprob_desc[] = "name of a GUROBI-format file to be written (for debugging);\n\
		must end in one of \".bas\", \".lp\", \".mps\", \".prm\", or \".sol\";\n\
		can appear more than once (with different filenames).";

 /* WS_desc_ASL = modified solvers/ws_desc.c: extra initial tab; must not be static. */
 char WS_desc_ASL[] = "=... solution report without -AMPL: sum of\n\
			1 ==> write .sol file\n\
			2 ==> print primal variable values\n\
			4 ==> print dual variable values\n\
			8 ==> do not print solution message";

#if GRB_VERSION_MAJOR > 1 /*{*/

 static char multprice_norm_desc[] = "choice of norm used in multiple pricing:\n\
		       -1 = automatic choice (default)\n\
			0, 1, 2, 3 = specific choices:  hard to describe,\n\
				but sometimes a specific choice will perform\n\
				much better than the automatic choice.";

 static char nodefiledir_desc[] = "directory where MIP tree nodes are written after memory\n\
		for them exceeds nodefilestart; default \".\"";

 static char nodefilestart_desc[] = "gigabytes of memory to use for MIP tree nodes;\n\
		default = Infinity (no limit, i.e., no node files written)";

#if GRB_VERSION_MAJOR >= 4 /*{*/
 static char premiqpmethod_desc[] = "how Gurobi's presolve should treat MIQP problems:\n\
			-1 = automatic choice (default)\n\
			 0 = leave the problem as an MIQP\n\
			 1 = try to transform an MIQP to an MILP";
#endif /*}*/

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char prepasses_desc[] = "limit on the number of Gurobi presolve passes:\n\
		       -1 = automatic choice (automatic)\n\
			n >= 0: at most n passes";
#endif /*}*/

 static char quad_desc[] = "whether simplex should use quad-precision:\n\
		       -1 = automatic choice (default)\n\
			0 = no\n\
			1 = yes";

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char resultfile_desc[] = "name of a file of extra information written after\n\
				completion of optimization.  The name's suffix\n\
				determines what is written:\n\
			.sol	solution vector\n\
			.bas	simplex basis\n\
			.mst	integer variable solution vector";

 static char rins_desc[] = "how often to apply the RINS heuristic for MIP problems:\n\
		       -1 = automatic choice (default)\n\
			0 = never\n\
			n > 0: every n-th node";
#endif /*}*/

 static char rootmethod_desc[] = "algorithm for MIP root relaxation:\n\
			0 = primal simplex\n\
			1 = dual simplex (default)"
#if GRB_VERSION_MAJOR >= 3
			"\n\
			2 = barrier"
#endif
			;
#endif /*}*/

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char symmetry_desc[] = "MIP symmetry detection:\n\
		       -1 = automatic choice (default)\n\
			0 = none\n\
			1 = conservative\n\
			2 = agressive";
#endif /*}*/

#define VP (void*)

#if GRB_VERSION_MAJOR >= 4
#define Method "Method"
#else
#define Method "LPMethod"
#endif

 static keyword
keywds[] = {	/* must be in alphabetical order */

#if GRB_VERSION_MAJOR >= 3
	{ "aggfill", sf_ipar, "AggFill", aggfill_desc },
#endif
	{ "aggregate", sf_ipar, "Aggregate", aggregate_desc },
#if GRB_VERSION_MAJOR >= 3 /*{*/
	{ "ams_eps", D_val, &ams_eps, ams_eps_desc },
	{ "ams_epsabs", D_val, &ams_epsabs, ams_epsabs_desc },
	{ "ams_limit", I_val, &ams_limit, ams_limit_desc },
	{ "ams_stub", C_val, &ams_stub, ams_stub_desc },
	{ "barconvtol", sf_dpar, "BarConvTol", barconvtol_desc },
	{ "barcorrectors", sf_ipar, "BarCorrectors", barcorrectors_desc },
	{ "bariterlim",  sf_ipar, "BarIterLimit", bariterlimit_desc },
	{ "barorder", sf_ipar, "BarOrder", barorder_desc },
#endif /*}*/
	{ "basis", sf_mint, VP set_basis, basis_desc },
	{ "bestbound", sf_mint, VP set_bestbound, bestbound_desc },
#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
	{ "branchdir", sf_ipar, "BranchDir", branchdir_desc },
#endif /*}*/
	{ "cliquecuts", sf_ipar, "CliqueCuts", overrides_cuts },
	{ "covercuts", sf_ipar, "CoverCuts", overrides_cuts },
#if GRB_VERSION_MAJOR >= 3 /*{*/
	{ "crossover", sf_ipar, "Crossover", crossover_desc },
	{ "crossoverbasis", sf_ipar, "CrossoverBasis", crossoverbasis_desc },
#endif /*}*/
	{ "cutagg", sf_ipar, "CutAggPasses", cutagg_desc },
	{ "cutoff", sf_dpar, "Cutoff", cutoff_desc },
#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
	{ "cutpasses", sf_ipar, "CutPasses", cutpasses_desc },
#endif /*}*/
	{ "cuts", sf_ipar, "Cuts", cuts_desc },
	{ "feastol", sf_dpar, "FeasibilityTol", feastol_desc },
	{ "flowcover", sf_ipar, "FlowCoverCuts", "flowcover cuts:  " Overrides_cuts },
	{ "flowpath", sf_ipar, "FlowPathCuts", "flowpath cuts:  " Overrides_cuts },
	{ "gomory", sf_ipar, "GomoryPasses", gomory_desc },
	{ "gubcover", sf_ipar, "GUBCoverCuts", "gubcover cuts:  " Overrides_cuts },
	{ "heurfrac", sf_dpar, "Heuristics", heurfrac_desc },
	{ "iisfind", sf_mint, VP set_iis, iisfind_desc },
#if GRB_VERSION_MAJOR > 1 /*{*/
	{ "iismethod", sf_ipar, "IISMethod", iismethod_desc },
#endif /*}*/
	{ "implied", sf_ipar, "ImpliedCuts", "implied cuts:  " Overrides_cuts },
#if GRB_VERSION_MAJOR >= 4 /*{*/
	{ "improvegap", sf_dpar, "ImproveStartGap", isg_desc },
	{ "improvetime", sf_dpar, "ImproveStartTime", ist_desc },
#endif /*}*/
	{ "intfeastol", sf_dpar, "IntFeasTol", intfeastol_desc },
	{ "intstart", sf_mint, VP set_intstart, intstart_desc },
	{ "iterlim", sf_dpar, "IterationLimit", "iteration limit (default: no limit)" },
	{ "logfile", C_val, &logfile, logfile_desc },
	{ "logfreq", sf_iparlog, "DisplayInterval", logfreq_desc },
	{ "lpmethod", sf_ipar, Method, "synonym for \"method\"" },
	{ "maxmipsub", sf_ipar, "SubMIPNodes", maxmipsub_desc },
	{ "method", sf_ipar, Method, simplex_desc },
#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
	{ "minrelnodes", sf_ipar, "MinRelNodes", minrelnodes_desc },
#endif /*}*/
#if GRB_VERSION_MAJOR >= 3
	{ "mipfocus", sf_ipar, "MIPFocus", mipfocus_desc },
#endif
	{ "mipgap", sf_dpar, "MipGap", "max. relative MIP optimality gap (default 1e-4)" },
#if GRB_VERSION_MAJOR >= 3
	{ "mipgapabs", sf_dpar, "MipGapAbs", "absolute MIP optimality gap (default 1e-10)" },
#endif
	{ "mipsep", sf_ipar, "MIPSepCuts", "MIPsep cuts:  " Overrides_cuts },
	{ "mipstart", sf_mint, VP set_mipstval, mipstart_desc },
	{ "mircuts", sf_ipar, "MIRCuts", "MIR cuts:  " Overrides_cuts },
#if GRB_VERSION_MAJOR >= 4
	{ "modkcuts", sf_ipar, "ModKCuts", "mod-k cuts:  " Overrides_cuts },
#endif
#if GRB_VERSION_MAJOR > 1 /*{*/
	{"multprice_norm", sf_ipar, "NormAdjust", multprice_norm_desc},
#if GRB_VERSION_MAJOR >= 3
	{ "networkcuts", sf_ipar, "NetworkCuts", "Network cuts:  " Overrides_cuts },
#endif
	{"nodefiledir", sf_spar, "NodefileDir", nodefiledir_desc},
	{"nodefilestart", sf_dpar, "NodefileStart", nodefilestart_desc},
#endif /*}*/
	{ "nodelim", sf_dpar, "NodeLimit", "maximum MIP nodes to explore (default: no limit)" },
#if GRB_VERSION_MAJOR >= 4 /*{*/
	{ "nodemethod", sf_ipar, "NodeMethod", nodemethod_desc },
#endif /*}*/
#if GRB_VERSION_MAJOR > 1 /*{*/
	{ "normadjust", sf_ipar, "NormAdjust", "synonym for multprice_norm" },
#endif /*}*/
	{ "objno", sf_mint, VP set_objno, objno_desc },
#if GRB_VERSION_MAJOR > 1 /*{*/
	{ "objscale", sf_dpar, "ObjScale", objscale_desc },
#endif /*}*/
	{ "opttol", sf_dpar, "OptimalityTol", opttol_desc },
	{ "outlev", sf_mint, VP set_outlev, outlev_desc },
#if GRB_VERSION_MAJOR >= 3
	{ "param", sf_spar, "Dummy", param_desc },
	{ "paramfile", sf_pf, 0, paramfile_desc },
#endif
	{ "perturb", sf_dpar, "PerturbValue", perturb_desc },
	{ "pivtol", sf_dpar, "MarkowitzTol", "Markowitz pivot tolerance (default 7.8125e-3)" },
	/*GRB_MAJ2(({"precrush", sf_ipar, "PreCrush", precrush_desc},))*/
#if GRB_VERSION_MAJOR >= 3
	{ "predeprow", sf_ipar, "PreDepRow", predeprow_desc },
	{ "predual", sf_ipar, "PreDual", predual_desc },
#if GRB_VERSION_MAJOR >= 4
	{ "premiqpmethod", sf_ipar, "PreMIQPMethod", premiqpmethod_desc },
#endif
	{ "prepases", sf_ipar, "PrePasses", "deprecated synonym for prepasses" },
	{ "prepasses", sf_ipar, "PrePasses", prepasses_desc },
#endif
	{ "presolve", sf_ipar, "Presolve", presolve_desc },
	{ "pricing", sf_ipar, "SimplexPricing", pricing_desc },
#if GRB_VERSION_MAJOR >= 4
	{ "psdtol", sf_dpar, "PSDTol", psdtol_desc },
#endif
#if GRB_VERSION_MAJOR >= 3
	{ "pumppasses", sf_ipar, "PumpPasses", pumppasses_desc },
#endif
#if GRB_VERSION_MAJOR > 1 /*{*/
	{ "quad", sf_ipar, "Quad", quad_desc},
#endif /*}*/
#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
	{ "rays", sf_mint, VP set_rays, rays_desc },
#endif /*}*/
	{ "relax", sf_mint, VP set_relax, relax_desc },
#if GRB_VERSION_MAJOR >= 3
	{ "resultfile", sf_spar, "ResultFile", resultfile_desc },
#endif
	{ "return_mipgap", sf_mint, VP set_retmipgap, return_mipgap_desc },
#if GRB_VERSION_MAJOR >= 3
	{ "rins", sf_ipar, "RINS", rins_desc },
#endif
#if GRB_VERSION_MAJOR > 1 /*{*/
	{"rootmethod", sf_ipar, "RootMethod", rootmethod_desc},
#endif /*}*/
	{ "scale", sf_ipar, "ScaleFlag", scale_desc },
	{ "simplex", sf_ipar, Method, "synonym for lpmethod" },
	{ "solnlimit", sf_ipar, "SolutionLimit", "maximum MIP solutions to find (default 2e9)" },
	{ "solnsens", sf_mint, VP set_solnsens, solnsens_desc },
	{ "sos", sf_mint, VP set_sos, sos_desc },
	{ "sos2", sf_mint, VP set_sos2, sos2_desc },
#if GRB_VERSION_MAJOR >= 3
	{ "submipcuts", sf_ipar, "SubMIPCuts", "sub-MIP cuts:  " Overrides_cuts },
	{ "symmetry", sf_ipar, "Symmetry", symmetry_desc },
#endif
	{ "threads", sf_ipar, "Threads", threads_desc },
	{ "timelim", sf_dpar, "TimeLimit", "limit on solve time (in seconds; default: no limit)" },
	{ "timing", sf_mint, VP set_timing, timing_desc },
	{ "varbranch", sf_ipar, "VarBranch", varbranch_desc },
	{ "wantsol", WS_val, 0, WS_desc_ASL+5 },
	{ "writeprob", sf_wfile, 0, writeprob_desc }
#if GRB_VERSION_MAJOR > 1 /*{*/
	,{"zerohalfcuts", sf_ipar, "ZeroHalfCuts", "zero-half cuts:  " Overrides_cuts }
#endif /*}*/
	};
#undef Method

 static Option_Info
Oinfo = { "gurobi", verbuf, "gurobi_options", keywds, nkeywds, 0, verbuf,
	   0,0,0,0,0, 20110907 };

 static void
enamefailed(GRBenv *env, const char *what, const char *name)
{
	fprintf(Stderr, "%s(\"%s\") failed:\n\t%s.\n", what, name, GRBgeterrormsg(env));
	++Oinfo.n_badopts;
	}

 static char iis_table[] = "\n\
0	non	not in the iis\n\
1	low	at lower bound\n\
2	fix	fixed\n\
3	upp	at upper bound\n\
4	mem	member\n\
5	pmem	possible member\n\
6	plow	possibly at lower bound\n\
7	pupp	possibly at upper bound\n\
8	bug\n";

 static SufDecl
suftab[] = {
	{ "absmipgap", 0, ASL_Sufkind_obj   | ASL_Sufkind_outonly },
	{ "absmipgap", 0, ASL_Sufkind_prob  | ASL_Sufkind_outonly },
	{ "bestbound", 0, ASL_Sufkind_obj   | ASL_Sufkind_outonly },
	{ "bestbound", 0, ASL_Sufkind_prob  | ASL_Sufkind_outonly },
#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
	{ "dunbdd", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
#endif /*}*/
	{ "iis", iis_table, ASL_Sufkind_var | ASL_Sufkind_outonly },
	{ "iis", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
	{ "ref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "relmipgap", 0, ASL_Sufkind_obj   | ASL_Sufkind_outonly },
	{ "relmipgap", 0, ASL_Sufkind_prob  | ASL_Sufkind_outonly },
	{ "senslbhi",  0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "senslblo",  0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensobjhi", 0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensobjlo", 0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensrhshi", 0, ASL_Sufkind_con | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensrhslo", 0, ASL_Sufkind_con | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensubhi",  0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensublo",  0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sos", 0, ASL_Sufkind_var },
	{ "sos", 0, ASL_Sufkind_con },
	{ "sosno", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "sosref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "sstatus", 0, ASL_Sufkind_var, 1 },
	{ "sstatus", 0, ASL_Sufkind_con, 1 }
#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
	,{ "unbdd", 0, ASL_Sufkind_var | ASL_Sufkind_outonly }
#endif /*}*/
	};
 static void
show_times(void)
{
	FILE *f;
	int i;

	Times[3] = xectim_();
	Times[4] = time(0) - Times[4];
	for(i = 1; i <= 2; i++)
	    if (time_flag & i) {
		f = i == 1 ? stdout : Stderr;
		fprintf(f, "\nTimes (seconds):\nInput =  %g"
			"\nSolve =  %g (summed over threads)"
			"\nOutput = %g\nElapsed ",
			Times[1] - Times[0], Times[2] - Times[1],
			Times[3] - Times[2]);
		fprintf(f, Times[4] < 1. ? "< 1\n" : "= %g\n", Times[4]);
		}
	}

 Sig_ret_type
intcatch(int n)
{
	printf("\n<BREAK> (gurobi)\n", n);
	fflush(stdout);
	if (++breaking > 3)
		longjmp(Jb, 2);
	signal(SIGINT, intcatch);
	if (grbmodel)
		GRBterminate(grbmodel);
	SigRet;
	}

 static char*
retiis(ASL *asl, GRBenv *env, GRBmodel *mdl, Dims *d, const char *what, int *srp)
{
	char buf[128], *rv;
	int *c, i, j, k, kv, m, n, nr, *s, *v;

	m = n_con;
	n = n_var;
	nr = n + nranges;
	if (GRBgetintattrarray(mdl, "IISConstr", 0, m, s = d->rstat))
		failed(env, "GRBgetintattrarray(\"IISConstr\")");
	c = v = 0;
	for(i = k = kv = 0; i < m; ++i) {
		if (s[i]) {
			c = (int*)M1zapalloc(m*sizeof(int));
			for(; i < m; ++i)
				if (s[i]) {
					c[i] = 4;
					++k;
					}
			break;
			}
		}
	if (GRBgetintattrarray(mdl, "IISLB", 0, nr, s = d->cstat))
		failed(env, "GRBgetintattrarray(\"IISLB\")");
	for(i = 0; i < n; ++i) {
		if (s[i]) {
			v = (int*)M1zapalloc(n*sizeof(int));
			for(; i < n; ++i)
				if (s[i]) {
					v[i] = 1;
					++kv;
					}
			break;
			}
		}
	for(i = n; i < nr; ++i) {
		if (s[i]) {
			if (!c)
				c = (int*)M1zapalloc(m*sizeof(int));
			for(; i < nr; ++i)
				if (s[i] && !c[i]) {
					c[i] = 4;
					++k;
					}
			break;
			}
		}
	if (GRBgetintattrarray(mdl, "IISUB", 0, nr, s = d->cstat))
		failed(env, "GRBgetintattrarray(\"IISUB\")");
	for(i = 0; i < n; ++i) {
		if (s[i]) {
			if (!v)
				v = (int*)M1zapalloc(n*sizeof(int));
			for(; i < n; ++i)
				if (s[i] && !v[i]) {
					v[i] = 3;
					++kv;
					}
			break;
			}
		}
	for(i = n; i < nr; ++i) {
		if (s[i]) {
			if (!c)
				c = (int*)M1zapalloc(m*sizeof(int));
			for(; i < nr; ++i)
				if (s[i] && !c[i]) {
					c[i] = 4;
					++k;
					}
			break;
			}
		}
	if (c)
		suf_iput("iis", ASL_Sufkind_con, c);
	if (v)
		suf_iput("iis", ASL_Sufkind_var, v);
	*srp = 201;
	if (k) {
		j = Snprintf(buf, sizeof(buf), "%s\nReturning an IIS of %d constraints",
			what, k);
		j += Snprintf(buf+j, sizeof(buf)-j, kv ? " and %d variables." : ".", kv);
		}
	else if (kv)
		j = Snprintf(buf, sizeof(buf), "%s\nReturning an IIS of %d variables.",
			what, kv);
	else {
		j = Snprintf(buf, sizeof(buf), "%s; empty IIS!");
		*srp = 202;
		}
	rv = (char*)M1alloc(++j);
	memcpy(rv, buf, j);
	return rv;
	}

 static void
dpf(Dims  *d, const char *fmt, ...)
{
	size_t L;
	va_list ap;

	if ((L = d->mbend - d->mb) > 0) {
		va_start(ap, fmt);
		d->mb += Vsnprintf(d->mb, L, fmt, ap);
		va_end(ap);
		}
	}

 static void
missing_msg(Dims *d)
{
	static const char *missing[3] = {
		"primal",
		"dual",
		"primal or dual"
		};
	dpf(d, "\nNo %s variables returned.", missing[d->missing - 1]);
	}

#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
#ifdef RAYDEBUG
#define Debug(x) x
#else
#define Debug(x) /*nothing*/
#endif

 static int
send_ray(ASL *asl, Dims *d, GRBenv *env, GRBmodel *mdl)
{
	int i, n, n0;
	real *c, t, *y;

	n = n_var;
	n0 = d->nv0;
	y = (real *)M1zapalloc(n0 * 2*sizeof(real));
	c = y + n;
	if (n > n0)
		n = n0;
	if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_UNBDRAY, 0, n, y)) {
		Debug(printf("Get UnbdRay failed: %s\n", GRBgeterrormsg(env)));
		return 0;
		}
	if (!GRBgetdblattrarray(mdl, GRB_DBL_ATTR_OBJ, 0, n, c)) {
		t = 0.;
		for(i = 0; i < n; ++i)
			t += c[i]*y[i];
		if (d->objsense * t > 0.)
			for(i = 0; i < n; ++i)
				y[i] = -y[i];
		}
	suf_rput("unbdd", ASL_Sufkind_var, y);
	return 1;
	}

 static int
send_dray(ASL *asl, Dims *d, GRBenv *env, GRBmodel *mdl)
{
	char *sense;
	int i, n, n0;
	real *rhs, t, t1, *y;

	n = n_con;
	n0 = d->nc0;
	y = (real *)M1zapalloc(n*sizeof(real) + n0*(sizeof(real) + 1));
	rhs = y + n;
	sense = (char*)(rhs + n0);
	if (n > n0)
		n = n0;
	if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_FARKASDUAL, 0, n, y)) {
		Debug(printf("Get FarkasDual failed: %s\n", GRBgeterrormsg(env)));
		return 0;
		}
	if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_RHS, 0, n, rhs))
		Debug(printf("Get RHS failed: %s\n", GRBgeterrormsg(env)));
	else if (GRBgetcharattrarray(mdl, GRB_CHAR_ATTR_SENSE, 0, n, sense))
		Debug(printf("Get SENSE failed: %s\n", GRBgeterrormsg(env)));
	else {
		t = 0.;
		for(i = 0; i < n; ++i) {
			t1 = y[i]*rhs[i];
			if (sense[i] == '>')
				t1 = -t1;
			t += t1;
			}
		if (d->objsense * t < 0.)
			for (i = 0; i < n; ++i)
				y[i] = -y[i];
		}
	suf_rput("dunbdd", ASL_Sufkind_con, y);
	return 1;
	}

#undef Debug
#endif /*}*/

 static char*
statmsg(ASL *asl, GRBenv *env, GRBmodel *mdl, int i, Dims *d, int *wantobj)
{
	char buf[64], *rv, *rv1;
	int m, n, nc, nv, nvr, objwant, sr, srd, srp;
	real *x, *y;
	size_t L;

	*wantobj = 0;
	if (i) {
		solve_result_num = 502;
		d->x = d->y = 0;
		switch(i) {
		  case GRB_ERROR_OUT_OF_MEMORY:
			rv = "ran out of memory";
			break;
		  case GRB_ERROR_NO_LICENSE:
			rv = "invalid license";
			break;
		  case GRB_ERROR_SIZE_LIMIT_EXCEEDED:
			rv = "problem size limit exceeded";
			break;
		  case GRB_ERROR_IIS_NOT_INFEASIBLE:
			rv = "bug: IIS problem is infeasible";
			break;
#if GRB_VERSION_MAJOR >= 4
		  case GRB_ERROR_Q_NOT_PSD:
			rv = "quadratic objective is not positive definite";
			solve_result_num = 524;
			break;
#endif
		  default:
			Snprintf(buf, sizeof(buf), "surprise return %d from GRBoptimize", i);
			rv = M1alloc(strlen(buf)+1);
			strcpy(rv, buf);
		  }
		return rv;
		}
	if (GRBgetintattr(mdl, GRB_INT_ATTR_STATUS, &i))
		failed(env, "GRBgetintattr(STATUS)");
	nc = n_con;
	nv = n_var;
	nvr = nv + nranges;
	objwant = 1;
	sr = 0;
	switch(i) {
	  case GRB_OPTIMAL:
		rv = "optimal solution";
		break;
	  case GRB_INFEASIBLE:
		objwant = 0;
		nc = srd = 0;
#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
		if (rays & 2)
			srd = send_dray(asl, d, env, mdl);
#endif /*}*/
		if (!want_iis) {
			sr = 200;
			rv = "infeasible";
			}
		else if (GRBcomputeIIS(mdl)) {
			sr = 202;
			rv = "infeasible; no IIS found";
			}
		else
			rv = retiis(asl, env, mdl, d, "infeasible", &sr);
#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
		if (srd) {
 have_srd:
			sr += 3;
			rv1 = (char*)M1alloc(L = strlen(rv) + 40);
			snprintf(rv1, L, "%s; constraint.dunbdd returned.", rv);
			rv = rv1;
			}
#endif /*}*/
		break;
	  case GRB_INF_OR_UNBD:
		objwant = 0;
		nc = nv = 0;
#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
		srd = srp = 0;
		if (rays & 1)
			srp = send_ray(asl, d, env, mdl);
		if (rays & 2 && srp != 1)
			srd = send_dray(asl, d, env, mdl);
#endif /*}*/
		if (!want_iis) {
			sr = 300;
			rv = "infeasible or unbounded";
			}
		else if (GRBcomputeIIS(mdl)) {
			sr = 301;
			rv = "infeasible or unbounded; no IIS";
			srp = srd = 0;
			}
		else
			rv = retiis(asl, env, mdl, d, "infeasible or unbounded", &sr);
#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
		if (srd)
			goto have_srd;
		if (srp)
			goto have_srp;
#endif /*}*/
		break;
	  case GRB_UNBOUNDED:
		rv = "unbounded";
		objwant = 0;
		nc = nv = srp = 0;
		sr = 300;
#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
		if (rays & 1) {
			srp = send_ray(asl, d, env, mdl);
 have_srp:
			sr += 2;
			rv1 = (char*)M1alloc(L = strlen(rv) + 30);
			snprintf(rv1, L, "%s; variable.unbdd returned.", rv);
			rv = rv1;
			}
#endif /*}*/
		break;
	  case GRB_CUTOFF:
		rv = "cutoff reached";
		sr = 400;
		break;
	  case GRB_ITERATION_LIMIT:
		rv = "iteration limit";
		sr = 401;
		break;
	  case GRB_NODE_LIMIT:
		rv = "node limit";
		sr = 402;
		break;
	  case GRB_TIME_LIMIT:
		rv = "time limit";
		sr = 403;
		break;
	  case GRB_SOLUTION_LIMIT:
		rv = "solution limit";
		sr = 404;
		break;
	  case GRB_INTERRUPTED:
		rv = "interrupted";
		sr = 600;
		break;
	  case GRB_NUMERIC:
		rv = "numeric error";
		sr = 520;
		break;
#ifdef GRB_SUBOPTIMAL
	  case GRB_SUBOPTIMAL:
		rv = "suboptimal";
		sr = 100;
		break;
#endif
	  default:
		Snprintf(buf, sizeof(buf), "surprise status %d after GRBoptimize", i);
		rv = (char*)M1alloc(strlen(buf)+1);
		sr = 530;
		strcpy(rv, buf);
	  }
	solve_result_num = sr;
	x = y = 0;
	if ((n = nvr + nc) > 0) {
		x = (real*)M1alloc(n*sizeof(real));
		d->y0 = y = x + nvr;
		}
	m = 0;
	if (!nv)
		x = 0;
	else if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_X, 0, nvr, x)) {
		x = 0;
		m = 1;
		}
	if (!nc)
		y = 0;
	else if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_PI, 0, nc, y)) {
		y = 0;
		m += 2;
		}
	d->missing = m;
	d->x = x;
	d->y = y;
	if (!x)
		objwant = 0;
	*wantobj = objwant;
	return rv;
	}

 static void
stat_map(int *stat, int n, int *map, int mx, char *what)
{
	int bad, i, i1, j, j1;
	static char badfmt[] = "gurobi driver: %s[%d] = %d\n";

	bad = i1 = j1 = 0;
	for(i = 0; i < n; i++) {
		if ((j = stat[i]) >= 0 && j <= mx)
			stat[i] = map[j];
		else {
			stat[i] = 0;
			i1 = i;
			j1 = j;
			if (!bad++)
				fprintf(Stderr, badfmt, what, i, j);
			}
		}
	if (bad > 1) {
		if (bad == 2)
			fprintf(Stderr, badfmt, what, i1, j1);
		else
			fprintf(Stderr,
		"gurobi driver: %d messages about bad %s values suppressed.\n",
				bad-1, what);
		}
	}

 static void
get_input_statuses(ASL *asl, GRBenv *env, GRBmodel *mdl, Dims *d)
{
	SufDesc *sd;
	int i, m, n, nvr, *rs, *rsta;
	real *lu;
	static int vmap[] = {-3, 0, -3, -1, -2, -3, -3};
	static int cmap[] = {-1, 0, -1, -1, -1, -1, -1};

	sd = d->csd;
	n = n_var;
	if (!(sd->kind & ASL_Sufkind_input))
		return;
	sd = d->rsd;
	m = n_con;
	if (!(sd->kind & ASL_Sufkind_input))
		return;
	stat_map(d->cstat, n, vmap, 6, "incoming cstat");
	stat_map(d->rstat, m, cmap, 6, "incoming rstat");
	nvr = n + nranges;
	if (nvr > n) {
		rs = d->rstat;
		rsta = d->cstat + n;
		lu = LUrhs;
		for(i = 0; i < m; ++i, lu += 2) {
			if (lu[0] > negInfinity && lu[0] < lu[1] && lu[1] < Infinity) {
				if (rs[i] == 0)
					*rsta = -1;
				++rsta;
				}
			}
		}
	if (GRBsetintattrarray(mdl, GRB_INT_ATTR_VBASIS, 0, nvr, d->cstat))
		failed(env, "GRBsetintattrarray(\"VBasis\")");
	if (GRBsetintattrarray(mdl, GRB_INT_ATTR_CBASIS, 0, m, d->rstat))
		failed(env, "GRBsetintattrarray(\"CBasis\")");
	return;
	}

 static void
intbasis_fail(Dims *d, const char *call)
{
	dpf(d, "\nintbasis trouble: GRB%s failed.", call);
	}

 static GRBmodel*
fixed_model(ASL *asl, GRBmodel *mdl0, Dims *d)
{
	GRBenv *env;
	GRBmodel *mdl;
	double f, *y;
	int i;
	static char *statusname[] = {
		"infeasible",
		"infeasible or unbounded",
		"unbounded",
		"cutoff",
		"iteration limit",
		"node limit",
		"time limit",
		"solution limit",
		"interrupted",
		"numeric difficulty"
		};

	if (!(mdl = GRBfixedmodel(mdl0)))
		return 0;
	if (!(env = GRBgetenv(mdl))) {
		dpf(d, "\nGRBgetenv failed in fixed_model().");
 badret:
		GRBfreemodel(mdl);
		return 0;
		}
	if (GRBsetintparam(env, "Presolve", 0)) {
		intbasis_fail(d, "setintparam(\"Presolve\")");
		goto badret;
		}
	if (GRBoptimize(mdl)) {
		intbasis_fail(d, "optimize()");
		goto badret;
		}
	if (GRBgetintattr(mdl, GRB_INT_ATTR_STATUS, &i)) {
		intbasis_fail(d, "getintattr()");
		goto badret;
		}
	if (i != GRB_OPTIMAL) {
		if (i >= GRB_INFEASIBLE && i <= GRB_NUMERIC)
			dpf(d, "\nGRBoptimize of fixed model: %s.",
				statusname[i-GRB_INFEASIBLE]);
		else
			dpf(d, "\nSurprise status %d after GRBoptimize of fixed model.",
				i);
		goto badret;
		}
	if (d->missing & 2 && (y = d->y0)
	 && !GRBgetdblattrarray(mdl, GRB_DBL_ATTR_PI, 0, n_con, y)) {
		d->y = y;
		d->missing &= ~2;
		}
	if (!GRBgetdblattr(mdl, GRB_DBL_ATTR_ITERCOUNT, &f)) {
		if (f > 0.)
			dpf(d, "\nplus %.0f simplex iteration%s for intbasis",
				f, "s" + (f == 1.));
		}
	return mdl;
	}

 static void
get_output_statuses(ASL *asl, GRBmodel *mdl, Dims *d)
{
	int i, j, m, n, nr, rv, *s;
	static int vmap[4] = {2, 4, 3, 1};

	m = n_con;
	n = n_var;
	nr = n + nranges;
	rv = 1;
	if (GRBgetintattrarray(mdl, GRB_INT_ATTR_VBASIS, 0, nr, d->cstat)
	 || GRBgetintattrarray(mdl, GRB_INT_ATTR_CBASIS, 0, m,  d->rstat)) {
		/*failed(env, "GRBgetintattrarray(\"VBasis\")");*/
		d->csd->kind &= ~ASL_Sufkind_output;
		d->rsd->kind &= ~ASL_Sufkind_output;
		goto ret;
		}
	rv = 0;
	s = d->cstat;
	for(i = 0; i < n; ++i) {
		if ((j = s[i] + 3) >= 0 && j <= 3)
			s[i] = vmap[j];
		else {
			badretfmt(504, "Surprise VBasis[%d] = %d.", i, j-3);
			goto ret;
			}
		}
	s = d->rstat;
	for(i = 0; i < m; ++i) {
		j = ++s[i];
		if (j < 0 || j > 1) {
			badretfmt(505, "Surprise CBasis[%d] = %d.", i, j-1);
			goto ret;
			}
		}
 ret:
	if (rv)
		dpf(d, "\nNo basis.");
	}

 static void
nl_iv_adj(ASL *asl, int j, int k, char *vtype, real *x)
{
	/* This will be needed once gurobi can handle nonlinear discrete variables, */
	/* e.g., in QPs. */

	int i, i0;
	real *L, *U;

	L = LUv;
	U = Uvx;
	i0 = k - j;
	if (vtype)
		for(i = i0; i < k; ++i)
			vtype[i] = L[i] == 0. && U[i] == 1. ? 'B' : 'I';
	if (x)
		for(i = i0; i < k; ++i) {
			if (x[i] < L[i])
				x[i] = L[i];
			else if (x[i] > U[i])
				x[i] = U[i];
			else
				x[i] = floor(x[i] + .5);
			}
	}

 typedef struct
Sensname {
	char *aname, *gname;
	int iscon;
	} Sensname;

 static void
put_sens(ASL *asl, GRBmodel *mdl)
{
	Sensname *sn, *sne;
	int len[2], nc, nv;
	real *a;
	static Sensname Snames[] = {
		{ "senslbhi",  "SALBUp", 0 },
		{ "senslblo",  "SALBLow", 0 },
		{ "sensobjhi", "SAObjUp", 0 },
		{ "sensobjlo", "SAObjLow", 0 },
		{ "sensrhshi", "SARHSUp", 1 },
		{ "sensrhslo", "SARHSLow", 1 },
		{ "sensubhi",  "SAUBUp", 0 },
		{ "sensublo",  "SAUBLow", 0 }};
	static int ak[2] = {ASL_Sufkind_var, ASL_Sufkind_con};

	len[0] = nv = n_var;
	len[1] = nc = n_con;
	a = (real*)M1alloc((6*nv + 2*nc)*sizeof(real));
	for(sn = Snames, sne = sn + sizeof(Snames)/sizeof(Sensname); sn < sne; ++sn) {
		if (GRBgetdblattrarray(mdl, sn->gname, 0, len[sn->iscon], a))
  	  		namefailed("GRBgetdblattrarray", sn->gname);
		suf_rput(sn->aname, ak[sn->iscon], a);
		a += len[sn->iscon];
		}
	}

#if 0
 static void
make_int(int j, int k, real *x)
{
	int i;
	for(i = k - j; i < k; ++i)
		x[i] = floor(x[i] + .5);
	}
#endif

#if 0
 static int __stdcall
my_querycb(GRBmodel *model, void *qcbdata, int where, void *usrdata)
{
	static time_t last_time;
	time_t now = time((time_t*)0);
	if (last_time != now) {
		fflush(stdout);
		last_time = now;
		}
	return 0;
	}
#endif

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static int
ams_write(ASL *asl, GRBenv *env, GRBmodel *mdl, Dims *dims, int nsols,
	  real bestobj, int objprec)
{
	Option_Info oi;
	char *fname, *fname_end, msg[64];
	enum {fname_endlen = 32};
	int havetol, i, j, nvr;
	real ba, *c, obj, t, ta, *x;
	size_t L;

	memset(&oi, 0, sizeof(oi));
	oi.wantsol = 9;
	if (--nsols > ams_limit && ams_limit > 0)
		nsols = ams_limit;
	nvr = n_var + nranges;
	L = strlen(ams_stub);
	x = (real*)Malloc(nvr*sizeof(real) + L + fname_endlen);
	fname = (char*)(x + nvr);
	fname_end = fname + L;
	memcpy(fname, ams_stub, L);
	havetol = 0;
	if (ams_eps > 0. || ams_epsabs > 0.) {
		havetol = 1;
		if ((ba = bestobj) < 0.)
			ba = -ba;
		}
	c = dims->c;
	for(i = 1; i <= nsols; ++i) {
		if (GRBsetintparam(env, "SolutionNumber", i))
			namefailed("GRBsetintparam", "SolutionNumber");
		if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_Xn, 0, nvr, x))
			namefailed("GRBgetdblattrarray", GRB_DBL_ATTR_Xn);
		if (!GRBgetdblattr(mdl, GRB_DBL_ATTR_OBJVAL, &obj)) {
			obj = 0.;
			for(j = 0; j < nvr; ++j)
				obj += c[j]*x[j];
			}
		if (havetol) {
			t = dims->objsense*(obj - bestobj);
			if (ams_epsabs > 0. && t > ams_epsabs)
				break;
			if (ams_eps > 0. && t > 0.) {
				if ((ta = obj) < 0.)
					ta = -ta;
				if (ta < ba)
					ta = ba;
				if (t/ta > ams_eps)
					break;
				}
			}
		Snprintf(msg, sizeof(msg), "Alternative MIP solution %d, objective = %.*g",
			i, objprec, obj);
		Snprintf(fname_end, fname_endlen, "%d.sol", i);
		if (write_solf_ASL(asl, msg, x, 0, &oi, fname))
			break;
		dpf(dims, "\n%s", msg);
		}
	free(x);
	if (GRBsetintparam(env, "SolutionNumber", 0)) /* restore */
		namefailed("GRBsetintparam", "SolutionNumber");
	return i - 1;
	}
#endif /*}*/

#define MBL 8192

 static GRBenv *env0;

 FILE *GRB_nl_ASL;	/* Foil gcc optimization bug: with -O, the "nl = 0" */
			/* assignment below does not happen, and after longjmp */
			/* we get an erroneous attempt to fclose(nl). */

 int
main(int argc, char **argv)
{
	ASL *asl;
	Dims dims;
	/*FILE *nl;*/ #define nl GRB_nl_ASL
	Filename *fn;
	GRBenv *env;
	GRBmodel *mdl, *mdl1;
	char mbuf[MBL];
	char *hx0, *sense, *solmsg, *sostype, *stub, *vtype;
	int i, j, k, lvi, nc, nfree, nlvi, nsosnz, nrange, nsos;
	int nv, nvr, nz, nzcr, objprec, rc, wantobj;
	int *cs, *csr, *rnr, *rsta, *sosbeg, *sosind, *sostypes;
	int vinfo[3], *vlen, *vlenr;
	ograd *og;
	real absmipgap, f, obj, oc, relmipgap, t;
	real *A, *Ar, *bb, *lu, *lxr, *rhs, *sosref, *uxr, *x, *y;
	sig_func_type *oic;
	size_t L;
	static int sos_types[2] = { GRB_SOS_TYPE1, GRB_SOS_TYPE2 };
#if GRB_VERSION_MAJOR >= 4
	fint *colqf, nelqf, *rowqf;
	int *colq, i1, j1, nelq, *rowq;
	real *qmat;
	nelqf = 0;
#endif

	Times[0] = xectim_();
	Times[4] = time(0);
#ifdef LICENSE_FILE
	if (!(stub = getenv("GRB_LICENSE_FILE")) || !*stub)
		putenv("GRB_LICENSE_FILE=" LICENSE_FILE);
#endif

	oic = 0;
	env0 = 0;
	mdl = 0;
	memset(&dims, 0, sizeof(dims));
	vinfo[0] = vinfo[1] = vinfo[2] = 0;
	GRBversion(&vinfo[0], &vinfo[1], &vinfo[2]);
	Snprintf(verbuf, sizeof(verbuf), "Gurobi %d.%d.%d", vinfo[0], vinfo[1], vinfo[2]);
	Lic_info_add_ASL = "Portions Copyright Gurobi Optimization, Inc., 2008.";

	asl = ASL_alloc(ASL_read_fg);
	if (!(stub = getstub(&argv, &Oinfo)))
		usage_ASL(&Oinfo, 1);
	nl = jac0dim(stub, 0);
	if (nlc > 0) {
		asl->i.uinfo = "Gurobi can't handle nonlinear constraints.";
		solve_result_num = 522;
		rc = 1;
		goto bailout;
		}
	if (n_cc > 0) {
		asl->i.uinfo = "Gurobi can't handle complementarity constraints.";
		solve_result_num = 567;
		rc = 1;
		goto bailout;
		}
	if (!(nobjno = n_obj))
		objno = 0;
	rc = 1;
	/* Passing 0 file logfile; would make logfile an option if we could	*/
	/* specify it after processing $gurobi_options, but we cannot do so,	*/
	/* and we need to have env before reading $gurobi_options, as values	*/
	/* in env may be changed by $gurobi_options. */
	if ((!env0 && GRBloadenv(&env0,0)) || !env0) {
		solve_result_num = 500;
		solmsg = "Could not create the gurobi environment.";
		goto ws_now;
		}
	Oinfo.uinfo = (char*)env0;
	nlvi = nlvbi + nlvci + nlvoi;
	lvi = nbv + niv;
	dims.kiv = nlvi + lvi;
#if 0
	if (dims.kiv)
		basis = 0;
#endif
	GRBsetintparam(env0, "OutputFlag", 0); /* possibly changed by getopts... */
	rc = setjmp(Jb);
	if (rc) {
 bailout:
		if (nl)
			fclose(nl);
		--rc;
		if (solve_result_num > 0 && asl->i.uinfo){
			if (amplflag | (Oinfo.wantsol & 1))
				rc = 0;
			L = strlen(Oinfo.bsname) + strlen(asl->i.uinfo);
			sprintf(solmsg = (char*)M1alloc(L+3), "%s: %s",
				Oinfo.bsname, asl->i.uinfo);
 ws_now:
			write_sol(solmsg, 0, 0, &Oinfo);
			}
		goto done;
		}
	if (getopts(argv, &Oinfo)) {
		if ((solmsg = asl->i.uinfo))
			goto ws_now;
		solve_result_num = 503;
		solmsg = "Bad $gurobi_options.";
		goto ws_now;
		}
	if (relax)
		dims.kiv = 0;
	breaking = 3;
	oic = signal(SIGINT, intcatch);
	nrange = nranges;
	dims.nv0 = nv = n_var;
	nvr = nv + nrange;
	dims.nc0 = nc = n_con;
	nz = nzc;
	nzcr = nz + nrange;
	A_vals = A = (real*)Malloc((2*(nvr+nc)+nzcr+nrange) * sizeof(real)
				+ (2*nvr + nzcr + 1)*sizeof(int) + nc + nvr + nv);
	Ar = A + nz;
	LUv = lxr = A + nzcr;
	Uvx = uxr = lxr + nvr;
	pi0 = y = uxr + nvr;
	rhs = y + nc;
	lxr += nv;
	uxr += nv;
	vlen = vlenr = (int*)(rhs + nc);
	A_colstarts = cs = vlen + nvr;
	A_rownos = rnr = cs + nvr + 1;
	sense = (char*)(rnr + nzcr);
	rnr += nz;
	vtype = sense + nc;
	havex0 = hx0 = vtype + nvr;
	suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));
	dims.cstat = (int*)M1zapalloc((nvr+nc+2)*sizeof(int));
	dims.rstat = dims.cstat + nvr + 1;
	rsta = dims.cstat + nv;
	dims.csd = suf_iput("sstatus", ASL_Sufkind_var, dims.cstat);
	dims.rsd = suf_iput("sstatus", ASL_Sufkind_con, dims.rstat);
	if (dims.kiv)
		want_xpi0 = 1;
#if GRB_VERSION_MAJOR >= 4
	qp_read(nl,0);
#else
	fg_read(nl,0);
#endif
	nl = 0;	/* was closed by qp_read */
	nfree = nsos = 0;
	dims.c = 0;
	dims.objsense = 1;
	oc = 0.;
	if (objno > 0) {
		i = objno - 1;
#if GRB_VERSION_MAJOR < 4
		if (i >= n_obj - nlo) {
			asl->i.uinfo = "Gurobi cannot handle nonlinear objectives.";
			solve_result_num = ASL_readerr_nonlin;
			rc = 1;
			goto bailout;
			}
#else
		/* must call mqpcheck() after qp_read so objconst() will work right */
		nelqf = mqpcheck(i, &rowqf, &colqf, &qmat);
		if (nelqf < 0) {
			if (nelqf == -2) {
				solve_result_num = 523;
				asl->i.uinfo = "a quadratic objective involving division by 0";
				}
			else {
				solve_result_num = 521;
				asl->i.uinfo = "a nonlinear objective";
				}
			rc = 1;
				goto bailout;
			}
#endif
		if (objtype[i])
			dims.objsense = -1;
		oc = objconst(i);
		dims.c = M1zapalloc(nvr * sizeof(real));
		for(og = Ograd[i]; og; og = og->next)
			dims.c[og->varno] = og->coef;
		}
	if (dims.kiv || (sos
	    && suf_get("sosno", ASL_Sufkind_var | ASL_Sufkind_input)
	    && suf_get("ref", ASL_Sufkind_var | ASL_Sufkind_input))) {
		i = ASL_suf_sos_explict_free;
		if (!sos)
			i |= ASL_suf_sos_ignore_sosno;
		if (!sos2)
			i |= ASL_suf_sos_ignore_amplsos;
		if ((nsos = suf_sos(i, &nsosnz, &sostype, 0, 0,
				&sosbeg, &sosind, &sosref))) {
			nv = n_var;
			nc = n_con;
			nvr = nv + nrange;
			lvi = nbv + niv;
			}
		}
	for(i = 0; i < nv; ++i)
		*vlenr++ = cs[i+1] - cs[i];
	csr = cs + nv;
	for(lu = LUrhs, i = 0; i < nc; ++i, lu += 2) {
		rhs[i] = lu[0];
		if (lu[0] <= negInfinity) {
			if (lu[1] >= Infinity)
				++nfree;
			else {
				sense[i] = '<';
				rhs[i] = lu[1];
				}
			}
		else if (lu[1] >= Infinity)
			sense[i] = '>';
		else {
			sense[i] = '=';
			if (lu[1] > lu[0]) {
				*rnr++ = i;
				*csr++ = Ar - A;
				*Ar++ = -1.;
				*vlenr++ = 1;
				*lxr++ = 0.;
				*uxr++ = lu[1] - lu[0];
				*rsta++ = y[i] > 0. ? -1 : y[i] < 0. ? -2 : 0;
				}
			}
		}
	if (nfree) {
		fprintf(stderr, "Botch: gurobi cannot handle %d free rows.\n", nfree);
		return 1;
		}

	memset(vtype, 'C', nvr);
	if (dims.kiv) {
		if (nlvi) {
			k = nlvb;
			if ((j = nlvbi))
				nl_iv_adj(asl, j, k, vtype, X0);
			k = nlvc;
			if ((j = nlvci))
				nl_iv_adj(asl, j, k, vtype, X0);
			k += nlvo - nlvc;
			if ((j = nlvoi))
				nl_iv_adj(asl, j, k, vtype, X0);
			}
		k = nv - lvi;
		if ((j = nbv)) {
			memset(vtype+k, 'B', j);
			k += j;
			if (X0)
				nl_iv_adj(asl, j, k, 0, X0);
			}
		if ((j = niv)) {
			memset(vtype+k, 'I', j);
			k += j;
			if (X0)
				nl_iv_adj(asl, j, k, 0, X0);
			}
		}
	if (GRBloadmodel(env0, &mdl, "foo",
			nvr, nc, dims.objsense, oc, dims.c,
			sense, rhs, cs, vlen, A_rownos, A_vals,
			LUv, Uvx, vtype,  0,0) || !mdl)
		failed(env0, "GRBloadmodel");
	x = 0;
	if (X0 && (!dims.kiv || mipstval)) {
		x = X0;
		for(i = 0; i < nv; ++i)
			if (!hx0[i])
				x[i] = GRB_UNDEFINED;
		}
	free(A);

	if (!(env = GRBgetenv(mdl)))
		failed(env0, "GRBgetenv");
	Oinfo.uinfo = (char*)env;

#if GRB_VERSION_MAJOR >= 4
	if (nelqf) {
		for(i = j = nelq = 0; i < nv; ++i) {
			for(k = colqf[i+1]; j < k; ++j)
				if (rowqf[j] <= i)
					++nelq;
			}
		colq = (int*)Malloc(nelq*sizeof(int));
		for(i = j = j1 = 0; i < nv; ++i) {
			for(k = colqf[i+1]; j < k; ++j) {
				if ((i1 = rowqf[j]) <= i) {
					colq[j1] = i;
					rowqf[j1] = i1;
					qmat[j1] = qmat[j];
					if (i1 == i)
						qmat[j1] *= 0.5;
					++j1;
					}
				}
			}
		if (sizeof(fint) == sizeof(int))
			rowq = (int*)rowqf;
		else {
			rowq = (int*)Malloc(nelq*sizeof(int));
			for(i = 0; i < nelq; ++i)
				rowq[i] = rowqf[i];
			free(rowqf);
			}
		free(colqf);
		if (GRBaddqpterms(mdl, nelq, rowq, colq, qmat))
			failed(env, "GRBaddqpterms");
		free(rowq);
		free(qmat);
		}
#endif
	if (nsos) {
		sostypes = (int*)Malloc(nsos*sizeof(int));
		for(i = 0; i < nsos; ++i)
			sostypes[i] = sos_types[sostype[i] - '1'];
		if (GRBaddsos(mdl, nsos, nsosnz, sostypes, sosbeg, sosind, sosref))
			failed(env, "GRBaddsos");
		free(sostypes);
  	  	}
	if (basis & 1)
		get_input_statuses(asl, env, mdl, &dims);
	if (x && GRBsetdblattrarray(mdl, GRB_DBL_ATTR_START, 0, nv, x))
		failed(env, "GRBsetintattrarray(START)");
	if (logfile) {
		if (GRBsetintparam(env, "OutputFlag", 1))
			namefailed("GRBsetintparam", "OutputFlag");
#if GRB_VERSION_MAJOR == 1
		if (GRBsetlogfilename(env, logfile))
			failed(env, "GRBsetlogfilename");
#elif GRB_VERSION_MAJOR == 2
		if (GRBsetstrparam(env, "LogfileName", logfile))
			namefailed("GRBsetstrparam", "LogfileName");
#else
		if (GRBsetstrparam(env, GRB_STR_PAR_LOGFILE, logfile))
			namefailed("GRBsetstrparam", GRB_STR_PAR_LOGFILE);
#endif
		}
	else if (outlev && GRBsetintparam(env, "OutputFlag", 1))
		namefailed("GRBsetintparam", "OutputFlag");
	if ((fn = Wflist)) {
		GRBupdatemodel(mdl);
		do {
			if (GRBwrite(mdl, fn->name))
				enamefailed(env, "GRBwrite", fn->name);
			} while((fn = fn->next));
		}
#if GRB_VERSION_MAJOR >= 4 && GRB_VERSION_MINOR >= 5 /*{*/
	if (rays)
		GRBsetintparam(env, "InfUnbdInfo", 1);
#endif /*}*/
	breaking = 1;
	Times[1] = xectim_();
	grbmodel = mdl;
	i = GRBoptimize(mdl);
	grbmodel = 0;
	Times[2] = xectim_();
	solmsg = statmsg(asl, env, mdl, i, &dims, &wantobj);
	dims.mb = mbuf;
	dims.mbend = mbuf + sizeof(mbuf);
	dpf(&dims, "%s: %s", verbuf, solmsg);
	absmipgap = relmipgap = Infinity;
	objprec = 0;
	if (wantobj && !GRBgetdblattr(mdl, GRB_DBL_ATTR_OBJVAL, &obj)) {
		dpf(&dims, "; objective %.*g", objprec = obj_prec(), obj);
		if ((bestbound | (retmipgap^4)) && objno > 0) {
			if (GRBgetdblattr(mdl, GRB_DBL_ATTR_OBJBOUND, &f))
				f = Infinity * dims.objsense;
			else {
				if ((absmipgap = obj - f) < 0.)
					absmipgap = -absmipgap;
				if ((t = obj) < 0.)
					t = -t;
				relmipgap = absmipgap / (1e-10 + t);
				}
			if (retmipgap & 1) {
				bb = (real*)M1zapalloc(nobjno*sizeof(real));
				bb[objno - 1] = relmipgap;
				suf_rput("relmipgap", ASL_Sufkind_obj, bb);
				suf_rput("relmipgap", ASL_Sufkind_prob, bb);
				}
			if (retmipgap & 2) {
				bb = (real*)M1zapalloc(nobjno*sizeof(real));
				bb[objno - 1] = absmipgap;
				suf_rput("absmipgap", ASL_Sufkind_obj, bb);
				suf_rput("absmipgap", ASL_Sufkind_prob, bb);
				}
			if (bestbound) {
				bb = (real*)M1zapalloc(nobjno*sizeof(real));
				bb[objno - 1] = f;
				suf_rput("bestbound", ASL_Sufkind_obj, bb);
				suf_rput("bestbound", ASL_Sufkind_prob, bb);
				}
			}
		}
	else
		solnsens = 0;
#if GRB_VERSION_MAJOR >= 3
	if (!GRBgetintattr(mdl, GRB_INT_ATTR_BARITERCOUNT, &i) && i > 0)
		dpf(&dims, "\n%d barrier iterations", i);
#endif
	if (!GRBgetdblattr(mdl, GRB_DBL_ATTR_ITERCOUNT, &f) && f > 0.)
		dpf(&dims, "\n%.0f simplex iterations", f);
	if (dims.kiv && !GRBgetdblattr(mdl, GRB_DBL_ATTR_NODECOUNT, &f) && f > 0.)
		dpf(&dims, "\n%.0f branch-and-cut nodes", f);
	mdl1 = mdl;
	if (dims.kiv && ((basis & 2) | solnsens))
		mdl1 = fixed_model(asl, mdl, &dims);
	if (absmipgap > 0. && absmipgap < Infinity && !(retmipgap & 4))
		dpf(&dims, "\nabsmipgap = %.3g, relmipgap = %.3g", absmipgap, relmipgap);
	if (mdl1) {
		if (solnsens)
			put_sens(asl, mdl1);
		if (basis & 2)
			get_output_statuses(asl, mdl1, &dims);
		if (mdl1 != mdl)
			GRBfreemodel(mdl1);
		}
	if (dims.missing)
		missing_msg(&dims);
#if GRB_VERSION_MAJOR >= 3
	if (dims.kiv && ams_stub > 0 && wantobj
	 && !GRBgetintattr(mdl, GRB_INT_ATTR_SOLCOUNT, &i) && i > 1
	 && (k = ams_write(asl, env, mdl, &dims, i, obj, objprec))) {
		dpf(&dims, "\n%d alternative MIP solution%s written to %s1.sol",
			k, "s" + (k == 1), ams_stub);
		if (k >= 2)
			dpf(&dims, "\n%s %s2.sol", k == 2 ? "and" : "...", ams_stub);
		if (k < --i)
			dpf(&dims, "\nIgnoring %d other inferior alternative MIP solutions.",
				i - k);
		}
#endif
#if 0
	if (dims.kiv && dims.x) {
		/* return integer values for integer variables */
		if (nlvi) {
			k = nlvb;
			if ((j = nlvbi))
				make_int(j, k, dims.x);
			k = nlvc;
			if ((j = nlvci))
				make_int(j, k, dims.x);
			k += nlvo - nlvc;
			if ((j = nlvoi))
				make_int(j, k, dims.x);
			}
		k = nv - lvi;
		if ((j = nbv)) {make_int(j, k, dims.x);
			k += j;
			make_int(j, k, dims.x);
			}
		if ((j = niv)) {
			k += j;
			make_int(j, k, dims.x);
			}
		}
#endif
	write_sol(mbuf, dims.x, dims.y, &Oinfo);
	for(fn = Wflist1; fn; fn = fn->next)
		if (GRBwrite(mdl, fn->name))
			enamefailed(env, "GRBwrite", fn->name);
 done:
	if (oic)
		signal(SIGINT, oic);
	if (mdl)
		GRBfreemodel(mdl);
	if (env0)
		GRBfreeenv(env0);
	ASL_free(&asl);
	show_times();
	return rc;
	}
