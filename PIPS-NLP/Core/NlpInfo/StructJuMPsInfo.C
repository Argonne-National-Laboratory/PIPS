

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
	data_in->inputNlp = this;

	PAR_DEBUG("  in StructJuMPsInfo constr comm"<<data_in->stochNode->commWrkrs<<" "<<mpiComm<<"  "<<MPI_COMM_NULL);

	iAmDistrib=0;
	  if( MPI_COMM_NULL!=mpiComm) {
	    int size;
	    MPI_Comm_size(mpiComm, &size);
	    PAR_DEBUG("size of parallel procs "<<size);
	    iAmDistrib = size==1?0:1;
	  }
	createChildren(data_in,*stochInput);

}

StructJuMPsInfo::StructJuMPsInfo(sData *data_in, StructJuMPInput& in, const int idx)
	:sInfo(data_in)
{
	PAR_DEBUG("StructJuMPsInfo ( data_in , structJuMPInput, "<<idx<<") id ("<<nodeId()<<")");

	stochInput = &(dynamic_cast<StructJuMPInput&>(in));
	data_in->inputNlp = this;

	iAmDistrib = 0;
	if (MPI_COMM_NULL != mpiComm) {
		int size;
		MPI_Comm_size(mpiComm, &size);
		iAmDistrib = size == 1 ? 0 : 1;
	}
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

void StructJuMPsInfo::createChildren(sData *data_in, StructJuMPInput& in){
	PAR_DEBUG("createChildren");
//	int mype_;
//		MPI_Comm_rank(in.prob->comm/* MPI_COMM_WORLD*/, &mype_);

	for (size_t it = 0; it < data_in->children.size(); it++) {
		if (stochNode->children[it]->commWrkrs != MPI_COMM_NULL) {
			AddChild(new StructJuMPsInfo(data_in->children[it], in, it));
		} else {
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

	double objv = 0.0;
	double local_var[locNx];
	local_X.copyIntoArray(local_var);

	if(parent!=NULL)
	{
		PAR_DEBUG("StructJuMPsInfo - ObjValue - parent "<<parent);
		int parid = parent->stochNode->id();
		assert(parid == 0);
		double parent_var[parent->locNx];
		OoqpVector* parent_X = (vars_X.parent->vec);
		parent_X->copyIntoArray(parent_var);
		double obj;
		assert(nodeId() != 0);
		PAR_DEBUG("StructJuMPsInfo - ObjValue - "<<nodeId());
		CallBackData cbd = {stochInput->prob->userdata,nodeId(),nodeId()};
		stochInput->prob->eval_f(parent_var,local_var,&obj,&cbd);
		objv += obj;
	}
	else
	{
		if(gmyid == 0)
		{
			PAR_DEBUG("StructJuMPsInfo - ObjValue - parent null ");
			PAR_DEBUG("StructJuMPsInfo - ObjValue - "<<nodeId());
			assert(nodeId() == 0);
			CallBackData cbd = {stochInput->prob->userdata,nodeId(),nodeId()};
			double obj;
			stochInput->prob->eval_f(local_var,local_var,&obj,&cbd);
			objv += obj;
		}
	}

	for(size_t it=0;it<children.size();it++) {
		objv += children[it]->ObjValue(svars->children[it]);
	}
	PAR_DEBUG("ObjValue "<<objv);
	if(iAmDistrib){
		double robj;
		MPI_Allreduce(&objv, &robj, 1, MPI_DOUBLE, MPI_SUM, mpiComm);
		objv = robj;
	}
	PAR_DEBUG("return ObjValue "<<objv);
	return objv;
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
	sconeq->vec->copyFromArray(coninq);

	for(size_t it=0; it<children.size(); it++){
		(children[it])->ConstraintBody(svars->children[it],sconeq->children[it],sconinq->children[it]);
	}

}

void StructJuMPsInfo::Hessian(NlpGenVars * vars, SymMatrix *Hess)
{
	PAR_DEBUG("Hessian");
}

void StructJuMPsInfo::JacFull(NlpGenVars* vars, GenMatrix* JacA, GenMatrix* JaC)
{
	PAR_DEBUG("JacFull");

	//update A , B , C, D matrix
}

void StructJuMPsInfo::get_InitX0(OoqpVector* vX){
	PAR_DEBUG("get_InitX0"<<vX);
}

void StructJuMPsInfo::Hessian_FromSon(NlpGenVars * vars, double *tempFromParH){
	PAR_DEBUG("Hessian_FromSon");
}

void StructJuMPsInfo::writeSolution(NlpGenVars* vars)
{
	PAR_DEBUG("writeSolution");
}
