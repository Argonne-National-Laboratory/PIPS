

#include "StructJuMPsInfo.h"

StructJuMPsInfo::StructJuMPsInfo()
{

}

StructJuMPsInfo::StructJuMPsInfo(sData *data_in){
	std::cout<<"StructJuMPsInfo with data_in only"<<std::endl;
}
StructJuMPsInfo::StructJuMPsInfo(sData *data_in, stochasticInput& in){
	std::cout<<"StructJuMPsInfo with data_in & stochinput"<<std::endl;
}

//StructJuMPsInfo::StructJuMPsInfo(sData *data_in, StructJuMPInput& in){
//	std::cout<<"StructJuMPsInfo with data_in & structJuMPInput"<<std::endl;
//}


StructJuMPsInfo::~StructJuMPsInfo()
{

}


void StructJuMPsInfo::createChildren(sData *data_in, StructJuMPInput& in){

}

double StructJuMPsInfo::ObjValue(NlpGenVars * vars){
	return 0;
}

void StructJuMPsInfo::ConstraintBody(NlpGenVars * vars, OoqpVector *conEq,OoqpVector *conIneq){

}

int StructJuMPsInfo::ObjGrad(NlpGenVars * vars, OoqpVector *grad){
	return 0;
}

void StructJuMPsInfo::Hessian(NlpGenVars * vars, SymMatrix *Hess){

}

void StructJuMPsInfo::JacFull(NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC){

}

void StructJuMPsInfo::get_InitX0(OoqpVector* vX){

}

void StructJuMPsInfo::Hessian_FromSon(NlpGenVars * vars, double *tempFromParH){

}

void StructJuMPsInfo::ObjGrad_FromSon( NlpGenVars * vars, OoqpVector *grad,double *tempFromParGrad ){

}
