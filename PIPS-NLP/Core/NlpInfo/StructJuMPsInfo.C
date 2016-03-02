

#include "StructJuMPsInfo.h"
#include "StructJuMPInput.h"

#include "sData.h"
#include "sVars.h"
#include "sTree.h"

#include "../../PIPS-NLP/par_macro.h"

StructJuMPsInfo::StructJuMPsInfo()
{
	assert(false);
}
StructJuMPsInfo::StructJuMPsInfo(sData *data_in):sInfo(data_in)
{
	assert(false);
}

StructJuMPsInfo::StructJuMPsInfo(sData *data_in, stochasticInput& in)
	:sInfo(data_in)
{
	PAR_DEBUG("StructJuMPsInfo ( data_in , stochasticInput)  - id "<< nodeId());
	parent = NULL;
	stochInput = &(dynamic_cast<StructJuMPInput&>(in));

	PAR_DEBUG("  in StructJuMPsInfo constr comm"<<data_in->stochNode->commWrkrs<<" "<<mpiComm<<"  "<<MPI_COMM_NULL);

	//	iAmDistrib=0;
	//	  if( MPI_COMM_NULL!=mpiComm) {
	int size;
	MPI_Comm_size(mpiComm, &size);
	PAR_DEBUG("size of parallel procs "<<size);
	//	    iAmDistrib = size==1?0:1;
	//	  }
	assert(MPI_COMM_NULL!=mpiComm);
	assert(size == gnprocs);
	createChildren(data_in,*stochInput);

//	data_in->inputNlp = this;
}

StructJuMPsInfo::StructJuMPsInfo(sData *data_in, stochasticInput& in, const int idx)
	:sInfo(data_in)
{
	PAR_DEBUG("StructJuMPsInfo ( data_in , structJuMPInput, "<<idx<<") id ("<<nodeId()<<")");

	stochInput = &(dynamic_cast<StructJuMPInput&>(in));
//	data_in->inputNlp = this;

//	iAmDistrib = 0;
//	if (MPI_COMM_NULL != mpiComm) {
		int size;
		MPI_Comm_size(mpiComm, &size);
//		iAmDistrib = size == 1 ? 0 : 1;
//	}
	PAR_DEBUG("number of children "<<data_in->children.size());
	createChildren(data_in,*stochInput);
}


StructJuMPsInfo::~StructJuMPsInfo()
{

}

int StructJuMPsInfo::nodeId()
{
	return stochNode->id();
}

void StructJuMPsInfo::createChildren(sData *data_in, stochasticInput& in){
	PAR_DEBUG("createChildren");
//	int mype_;
//		MPI_Comm_rank(in.prob->comm/* MPI_COMM_WORLD*/, &mype_);

	for (size_t it = 0; it < data_in->children.size(); it++) {
		if (stochNode->children[it]->commWrkrs != MPI_COMM_NULL) {
			AddChild(new StructJuMPsInfo(data_in->children[it], in, it));
		}
		else {
			PAR_DEBUG("comm null "<<MPI_COMM_NULL<<" commwrk "<<stochNode->children[it]->commWrkrs);
			AddChild(new sInfoDummy());
		}
		children[it]->parent = this;
		//	this->iAmDistrib = 0;
	}
}

double StructJuMPsInfo::ObjValue(NlpGenVars * vars){
	PAR_DEBUG("StructJuMPsInfo - ObjValue - id "<<nodeId()<<" nlocalVar "<<locNx);
	sVars* svars = dynamic_cast<sVars*>(vars);
	StochVector& vars_X = dynamic_cast<StochVector&>(*svars->x);
	OoqpVector& local_X = *(dynamic_cast<StochVector&>(*svars->x).vec);

	double local_var[locNx];
	local_X.copyIntoArray(local_var);

	double robj = 0.0;
	if(parent==NULL)
	{

		double objv = 0.0;
		if(gmyid == 0)
		{
			PAR_DEBUG("StructJuMPsInfo - ObjValue - parent null "<<" nodeid"<<nodeId());
			assert(nodeId() == 0);
			CallBackData cbd = {stochInput->prob->userdata,nodeId(),nodeId()};
			double obj;
			stochInput->prob->eval_f(local_var,local_var,&obj,&cbd);
			objv += obj;
		}
		for(size_t it=0;it<children.size();it++) {
			objv += children[it]->ObjValue(svars->children[it]);
		}
		MPI_Allreduce(&objv, &robj, 1, MPI_DOUBLE, MPI_SUM, mpiComm);
		PAR_DEBUG(" local obj "<<objv<< "  global obj"<<robj);
	}
	else
	{
		PAR_DEBUG("StructJuMPsInfo - ObjValue - parent "<<parent<<" nodeid "<<nodeId());
		int parid = parent->stochNode->id();
		assert(parid == 0);
		double parent_var[parent->locNx];
		OoqpVector* parent_X = (vars_X.parent->vec);
		parent_X->copyIntoArray(parent_var);
		assert(nodeId() != 0);
		PAR_DEBUG("StructJuMPsInfo - ObjValue - "<<nodeId());
		CallBackData cbd = {stochInput->prob->userdata,nodeId(),nodeId()};
		stochInput->prob->eval_f(parent_var,local_var,&robj,&cbd);
	}

	PAR_DEBUG("return ObjValue "<<robj);
	return robj;
}

int StructJuMPsInfo::ObjGrad(NlpGenVars * vars, OoqpVector *grad){
	PAR_DEBUG("ObjGrad");
	sVars * svars = dynamic_cast<sVars*>(vars);
	StochVector& vars_X = dynamic_cast<StochVector&>(*svars->x);
	OoqpVector& local_X = *(dynamic_cast<StochVector&>(*svars->x).vec);
	StochVector* sGrad = dynamic_cast<StochVector*>(grad);
	PAR_DEBUG("sGrad size "<<sGrad->vec->n);

	assert(parent == NULL);
	std::vector<double> local_grad(locNx,0.0);
	for(size_t it=0; it<children.size(); it++){
		(children[it])->ObjGrad_FromSon(svars->children[it],sGrad->children[it], &local_grad[0]);
	}

	double local_var[locNx];
	local_X.copyIntoArray(local_var);

	if(gmyid == 0) {
//		PAR_DEBUG("S1 -- ----");
		double grad_temp[locNx];
		CallBackData cbd = {stochInput->prob->userdata, nodeId(), nodeId()};
//		PAR_DEBUG("S2 -- ----");
		stochInput->prob->eval_grad_f(local_var,local_var,grad_temp,&cbd);
		for(int i = 0;i<locNx;i++)
			local_grad[i] += grad_temp[i];
//		PAR_DEBUG("S3 -- ----");
	}

	double rgrad[locNx];
	MPI_Allreduce(&local_grad[0], rgrad, locNx, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	sGrad->vec->copyFromArray(rgrad);

	return 1;
}

void StructJuMPsInfo::ObjGrad_FromSon(NlpGenVars* vars, OoqpVector* grad, double* pgrad)
{
	PAR_DEBUG(parent);
	PAR_DEBUG("ObjGrad -- with parent vector size = "<<parent->locNx);
	assert(parent!=NULL);
//	assert(parent->locNx == pgrad.size());
//	PAR_DEBUG("S1 -- ----");
	sVars * svars = dynamic_cast<sVars*>(vars);
	StochVector& vars_X = dynamic_cast<StochVector&>(*svars->x);
	OoqpVector& local_X = *(dynamic_cast<StochVector&>(*svars->x).vec);
//	PAR_DEBUG("S2 -- ----");

	double local_var[locNx];
	local_X.copyIntoArray(local_var);
	PAR_DEBUG("local_var set size "<<locNx);

	assert(parent->stochNode->id()==0);
	double parent_var[parent->locNx];
	OoqpVector* parent_X = (vars_X.parent->vec);
	parent_X->copyIntoArray(parent_var);


	double parent_part[parent->locNx];
	double this_part[locNx];
	CallBackData cbd_parent = {stochInput->prob->userdata, nodeId(), parent->stochNode->id()};
	stochInput->prob->eval_grad_f(parent_var,local_var,parent_part,&cbd_parent);
	CallBackData cbd_this = {stochInput->prob->userdata, nodeId(), nodeId()};
	stochInput->prob->eval_grad_f(parent_var,local_var,this_part,&cbd_this);

	StochVector* sGrad = dynamic_cast<StochVector*>(grad);
	PAR_DEBUG("sGrad size "<<sGrad->vec->n);
	sGrad->vec->copyFromArray(this_part);
	for(int i = 0;i<parent->locNx;i++)
		pgrad[i] += parent_part[i];
}


void StructJuMPsInfo::ConstraintBody(NlpGenVars * vars, OoqpVector *conEq,OoqpVector *conIneq){
	PAR_DEBUG("ConstraintBody");
//	PAR_DEBUG("S1 -- ----");
	sVars * svars = dynamic_cast<sVars*>(vars);
	StochVector& vars_X = dynamic_cast<StochVector&>(*svars->x);
	OoqpVector& local_X = *(dynamic_cast<StochVector&>(*svars->x).vec);
//	PAR_DEBUG("S2 -- ----");

	double local_var[locNx];
	local_X.copyIntoArray(local_var);
	PAR_DEBUG("local_var set size "<<locNx);


	double coneq[locMy];
	double coninq[locMz];
	if(parent!=NULL) {
		assert(parent->stochNode->id()==0);
		double parent_var[parent->locNx];
		OoqpVector* parent_X = (vars_X.parent->vec);
		parent_X->copyIntoArray(parent_var);

		CallBackData cbd = {stochInput->prob->userdata,nodeId(),nodeId()};
		stochInput->prob->eval_g(parent_var,local_var,coneq,coninq,&cbd);
	}
	else
	{
		assert(nodeId()==0);
		CallBackData cbd = {stochInput->prob->userdata,nodeId(),nodeId()};
		stochInput->prob->eval_g(local_var,local_var,coneq,coninq,&cbd);
	}

	StochVector* sconeq = dynamic_cast<StochVector*>(conEq);
	PAR_DEBUG("sconeq size "<<sconeq->vec->n);
	assert(sconeq->vec->n == locMy);
	sconeq->vec->copyFromArray(coneq);

	StochVector* sconinq = dynamic_cast<StochVector*>(conIneq);
	PAR_DEBUG("sconinq size "<<sconinq->vec->n);
	assert(sconinq->vec->n == locMz);
	sconinq->vec->copyFromArray(coninq);

	for(size_t it=0; it<children.size(); it++){
		(children[it])->ConstraintBody(svars->children[it],sconeq->children[it],sconinq->children[it]);
	}

}

void StructJuMPsInfo::JacFull(NlpGenVars* vars, GenMatrix* JacA, GenMatrix* JaC)
{
//	note: no linking constraint handling
	PAR_DEBUG("JacFull");

	long long mA, nA, mC, nC, mB,nB,mD,nD;
	Amat->getSize(mA,nA);
	Cmat->getSize(mC,nC);
	Bmat->getSize(mB,nB);
	Dmat->getSize(mD,nD);
	PAR_DEBUG(" Amat "<<mA<<"  "<<nA<<"  nz"<<Amat->numberOfNonZeros());
	PAR_DEBUG(" Cmat "<<mC<<"  "<<nC<<"  nz"<<Cmat->numberOfNonZeros());
	PAR_DEBUG(" Bmat "<<mB<<"  "<<nB<<"  nz"<<Bmat->numberOfNonZeros());
	PAR_DEBUG(" Dmat "<<mD<<"  "<<nD<<"  nz"<<Dmat->numberOfNonZeros());

	sVars * svars = dynamic_cast<sVars*>(vars);
	StochVector& vars_X = dynamic_cast<StochVector&>(*svars->x);
	OoqpVector* local_X = (vars_X.vec);
	double local_var[locNx];
	local_X->copyIntoArray(local_var);
	//update A , B , C, D matrix
	if(parent == NULL){
		//only B D matrix
		assert(nodeId() == 0);
		PAR_DEBUG("JacFull -- parent is NULL");
		int e_nz = Bmat->numberOfNonZeros();
		int i_nz = Dmat->numberOfNonZeros();
		PAR_DEBUG("Bmat nz "<<e_nz<<" Dmat nz "<<i_nz);

		CallBackData cbd = {stochInput->prob->userdata,nodeId(),nodeId()};
////		have to request the structure again
//		stochInput->prob->eval_jac_g(local_var,local_var,
//				&e_nz,NULL,NULL,NULL,
//				&i_nz,NULL,NULL,NULL,&cbd);
//		PAR_DEBUG("Bmat nz "<<e_nz<<" Dmat nz "<<i_nz);  //should fixed structure from stochcasticInput interface

		std::vector<int> e_rowidx(e_nz);
		std::vector<int> e_colptr(locNx+1,0);
		std::vector<double> e_elts(e_nz);

		std::vector<int> i_rowidx(i_nz);
		std::vector<int> i_colptr(locNx+1,0);
		std::vector<double> i_elts(i_nz);

		stochInput->prob->eval_jac_g(local_var,local_var,
					&e_nz,&e_elts[0],&e_rowidx[0],&e_colptr[0],
					&i_nz,&i_elts[0],&i_rowidx[0],&i_colptr[0],&cbd);

		Bmat->copyMtxFromDouble(Bmat->numberOfNonZeros(),&e_elts[0]);
		Dmat->copyMtxFromDouble(Dmat->numberOfNonZeros(),&i_elts[0]);
	}
	else{
		//all A B C D
		PAR_DEBUG("JacFull -- with parent - nodeid "<<nodeId());
		double parent_var[parent->locNx];
		OoqpVector* parent_X = (vars_X.parent->vec);
		parent_X->copyIntoArray(parent_var);

		int e_nz_Amat = Amat->numberOfNonZeros();
		int i_nz_Cmat = Cmat->numberOfNonZeros();

		CallBackData cbd_link = {stochInput->prob->userdata,nodeId(),parent->stochNode->id()};
//		PAR_DEBUG("nz amat "<<e_nz_Amat<<"  cmat "<<i_nz_Cmat);
//		stochInput->prob->eval_jac_g(parent_var,local_var,
//					&e_nz_Amat,NULL,NULL,NULL,
//					&i_nz_Cmat,NULL,NULL,NULL,&cbd_link);
//		PAR_DEBUG("nz amat "<<e_nz_Amat<<"  cmat "<<i_nz_Cmat);  //should fixed structure from stochcasticInput interface

		int e_amat_rowidx[e_nz_Amat];
		int e_amat_colptr[parent->locNx+1];
		double e_amat_elts[e_nz_Amat];

		int i_cmat_rowidx[i_nz_Cmat];
		int i_cmat_colptr[parent->locNx+1];
		double i_cmat_elts[i_nz_Cmat];

		stochInput->prob->eval_jac_g(parent_var,local_var,
				&e_nz_Amat,e_amat_elts,e_amat_rowidx,e_amat_colptr,
				&i_nz_Cmat,i_cmat_elts,i_cmat_rowidx,i_cmat_colptr, &cbd_link);

		int e_nz_Bmat = Bmat->numberOfNonZeros();
		int i_nz_Dmat = Dmat->numberOfNonZeros();

		CallBackData cbd_diag = {stochInput->prob->userdata,nodeId(),nodeId()};
//		stochInput->prob->eval_jac_g(parent_var,local_var,
//						&e_nz_Bmat,NULL,NULL,NULL,
//						&i_nz_Dmat,NULL,NULL,NULL,&cbd_diag);

		int e_bmat_rowidx[e_nz_Bmat];
		int e_bmat_colptr[locNx+1];
		double e_bmat_elts[e_nz_Bmat];

		int i_dmat_rowidx[i_nz_Dmat];
		int i_dmat_colptr[locNx+1];
		double i_dmat_elts[i_nz_Dmat];

		stochInput->prob->eval_jac_g(parent_var,local_var,
				&e_nz_Bmat,e_bmat_elts,e_bmat_rowidx,e_bmat_colptr,
				&i_nz_Dmat,i_dmat_elts,i_dmat_rowidx,i_dmat_colptr, &cbd_diag);

		Amat->copyMtxFromDouble(Amat->numberOfNonZeros(),e_amat_elts);
		Bmat->copyMtxFromDouble(Bmat->numberOfNonZeros(),e_bmat_elts);
		Cmat->copyMtxFromDouble(Cmat->numberOfNonZeros(),i_cmat_elts);
		Dmat->copyMtxFromDouble(Dmat->numberOfNonZeros(),i_dmat_elts);
	}

	for(size_t it=0; it<children.size(); it++)
		children[it]->JacFull(svars->children[it], NULL,NULL);
}


void StructJuMPsInfo::Hessian(NlpGenVars * nlpvars, SymMatrix *Hess)
{
	PAR_DEBUG("Hessian -- at root");
	//update Qdiag and Qborder
	long long mqi, nqi, mqb, nqb;
	Qdiag->getSize(mqi,nqi);
	Qborder->getSize(mqb,nqb);
	PAR_DEBUG(" Qdiag "<<mqi<<"  "<<nqi<<"  nz"<<Qdiag->numberOfNonZeros());
	PAR_DEBUG(" Qborder "<<mqb<<"  "<<nqb<<"  nz"<<Qborder->numberOfNonZeros());

	sVars * vars = dynamic_cast<sVars*>(nlpvars);
	StochVector& vars_X = dynamic_cast<StochVector&>(*vars->x);
	StochVector& vars_Y = dynamic_cast<StochVector&>(*vars->y);
	StochVector& vars_Z = dynamic_cast<StochVector&>(*vars->z);
	OoqpVector* local_X = vars_X.vec;
	OoqpVector* local_Y = vars_Y.vec; //eq con
	OoqpVector* local_Z = vars_Z.vec; //ieq con

	double local_var[locNx];
	local_X->copyIntoArray(local_var);
	double local_y[locMy];
	double local_z[locMz];
	local_Y->copyIntoArray(local_y);
	local_Z->copyIntoArray(local_z);
	std::vector<double> lam(locMy+locMz,0.0);
	int i=0;
	for(i=0;i<locMy;i++) lam[i] = -local_y[i];
	for(;i<locMy+locMz;i++) lam[i] = -local_z[i];

	int nzqd = Qdiag->numberOfNonZeros();
	std::vector<double> elts(nzqd,0.0);

	if(gmyid == 0) {
		int rowidx[nzqd];
		int colptr[locNx+1];
		CallBackData cbd = {stochInput->prob->userdata,0,0};
		stochInput->prob->eval_h(local_var,local_var,&lam[0],&nzqd,&elts[0],rowidx,colptr,&cbd);
	}

	for(size_t it=0; it<children.size(); it++){
		children[it]->Hessian_FromSon(vars->children[it],&elts[0]);
	}

	//MPI ALL REDUCE
	double g_elts[nzqd];
	MPI_Allreduce(&elts[0], g_elts, nzqd, MPI_DOUBLE, MPI_SUM, mpiComm);

	Qdiag->copyMtxFromDouble(nzqd,g_elts);
}

void StructJuMPsInfo::Hessian_FromSon(NlpGenVars* nlpvars, double *parent_hess){
	PAR_DEBUG("Hessian_FromSon");
	long long mqi, nqi, mqb, nqb;
	Qdiag->getSize(mqi,nqi);
	Qborder->getSize(mqb,nqb);
	PAR_DEBUG(" Qdiag "<<mqi<<"  "<<nqi<<"  nz"<<Qdiag->numberOfNonZeros());
	PAR_DEBUG(" Qborder "<<mqb<<"  "<<nqb<<"  nz"<<Qborder->numberOfNonZeros());

	sVars * vars = dynamic_cast<sVars*>(nlpvars);
	StochVector& vars_X = dynamic_cast<StochVector&>(*vars->x);
	StochVector& vars_Y = dynamic_cast<StochVector&>(*vars->y);
	StochVector& vars_Z = dynamic_cast<StochVector&>(*vars->z);
	OoqpVector* local_X = vars_X.vec;
	OoqpVector* local_Y = vars_Y.vec; //eq con
	OoqpVector* local_Z = vars_Z.vec; //ieq con

	double local_var[locNx];
	local_X->copyIntoArray(local_var);
	double local_y[locMy];
	double local_z[locMz];
	local_Y->copyIntoArray(local_y);
	local_Z->copyIntoArray(local_z);
	std::vector<double> lam(locMy+locMz,0.0);
	int i=0;
	for(i=0;i<locMy;i++) lam[i] = -local_y[i];
	for(;i<locMy+locMz;i++) lam[i] = -local_z[i];

	double parent_var[parent->locNx];
	OoqpVector* parent_X = (vars_X.parent->vec);
	parent_X->copyIntoArray(parent_var);


	//pnzqd
	{
		int pnzqd = parent->Qdiag->numberOfNonZeros();
		PAR_DEBUG("  -- Qdiag nz parent "<<pnzqd);
		double elts[pnzqd];
		int rowidx[pnzqd];
		int colptr[parent->locNx+1];
		CallBackData cbd_pnzqd = {stochInput->prob->userdata,nodeId(),0};
		stochInput->prob->eval_h(parent_var,local_var,&lam[0],&pnzqd,elts,rowidx,colptr,&cbd_pnzqd);
		for(int i=0;i<pnzqd;i++) parent_hess[i] += elts[i];
	}

	//nzqd
	{
		int nzqd = Qdiag->numberOfNonZeros();
		double elts[nzqd];
		int rowidx[nzqd];
		int colptr[locNx+1];
		CallBackData cbd_nzqd = {stochInput->prob->userdata,nodeId(),nodeId()};
		stochInput->prob->eval_h(parent_var,local_var,&lam[0],&nzqd,elts,rowidx,colptr,&cbd_nzqd);
		Qdiag->copyMtxFromDouble(nzqd,elts);
	}

	//nzqb
	{
		int nzqb = Qborder->numberOfNonZeros();
		double elts[nzqb];
		int rowidx[nzqb];
		int colptr[parent->locNx+1];
		CallBackData cbd_nzqb = {stochInput->prob->userdata,0,nodeId()};
		stochInput->prob->eval_h(parent_var,local_var,&lam[0],&nzqb,elts,rowidx,colptr,&cbd_nzqb);
		Qborder->copyMtxFromDouble(nzqb,elts);
	}
}

void StructJuMPsInfo::get_InitX0(OoqpVector* vX){
	PAR_DEBUG("get_InitX0 - "<<nodeId());
	StochVector* vars_X = dynamic_cast<StochVector*>(vX);
	OoqpVector* local_X = vars_X->vec;
	PAR_DEBUG("size of x "<<locNx<<" - vX->n "<<vX->n);
//	assert(locNx == vX->n);
	assert(children.size() == vars_X->children.size());

	double temp_var[locNx];
	CallBackData cbd = {stochInput->prob->userdata,nodeId(),nodeId()};
	stochInput->prob->init_x0(temp_var,&cbd);
//	local_X->print();
	local_X->copyFromArray(temp_var);
//	local_X->print();
	for(size_t it=0; it<children.size(); it++)
	    children[it]->get_InitX0(vars_X->children[it]);
}

void StructJuMPsInfo::writeSolution(NlpGenVars* vars)
{
	PAR_DEBUG("writeSolution");
}
