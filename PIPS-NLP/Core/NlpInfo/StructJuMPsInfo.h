#ifndef STRUCTJUMPSINFO
#define STRUCTJUMPSINFO

#include "sInfo.h"

#include "StructJuMPInput.h"

class StructJuMPsInfo : public sInfo
{
protected:

public:
	StructJuMPsInfo();
	StructJuMPsInfo(sData *data_in);
	StructJuMPsInfo(sData *data_in, stochasticInput& in);
	~StructJuMPsInfo();

	void createChildren(sData *data_in, StructJuMPInput& in);

	virtual double ObjValue(NlpGenVars * vars);

	virtual void ConstraintBody(NlpGenVars * vars, OoqpVector *conEq,
			OoqpVector *conIneq);

	virtual int ObjGrad(NlpGenVars * vars, OoqpVector *grad);

	virtual void Hessian(NlpGenVars * vars, SymMatrix *Hess);

	virtual void JacFull(NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC);

	virtual void get_InitX0(OoqpVector* vX);

	virtual void Hessian_FromSon(NlpGenVars * vars, double *tempFromParH);

	virtual void ObjGrad_FromSon( NlpGenVars * vars, OoqpVector *grad,double *tempFromParGrad );
};


#endif
