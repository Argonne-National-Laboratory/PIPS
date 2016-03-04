#include "StructJuMPInput.h"

#include "../PIPS-NLP/par_macro.h"


#include <string>
#include <iostream>
#include <sstream>

StructJuMPInput::StructJuMPInput(PipsNlpProblemStruct* p)
{
	PAR_DEBUG("enter constructor StructJuMPInput - "<<p);
	this->prob = p;
	useInputDate = 1;
	datarootname = "StructJuMP";
	PAR_DEBUG("exit constructor StructJuMPInput - ");
}
StructJuMPInput::~StructJuMPInput() {

}

int StructJuMPInput::nScenarios() {
	return prob->nscen;
}

void StructJuMPInput::get_prob_info(int nodeid)
{
	PAR_DEBUG("get_prob_info  - nodeid "<<nodeid);
	int nv = 0;
	int mc = 0;

	CallBackData data={prob->userdata,nodeid,nodeid};
//	PAR_DEBUG("data -"<<&data);
	PAR_DEBUG("get_prob_info callback data ptr - "<<data.row_node_id);
	PAR_DEBUG("get_prob_info userdata ptr - "<<prob->userdata);
	prob->prob_info(&nv,NULL,NULL,&mc,NULL,NULL,&data);
	nvar_map[nodeid] = nv;
	ncon_map[nodeid] = mc;
	PAR_DEBUG("ncon,nvar"<<mc<<", "<<nv);

	std::vector<double> collb(nv);
	std::vector<double> colub(nv);
	std::vector<double> rowlb(mc);
	std::vector<double> rowub(mc);
	prob->prob_info(&nv, &collb[0], &colub[0], &mc, &rowlb[0], &rowub[0], &data);

	collb_map[nodeid] = collb;
	colub_map[nodeid] = colub;
	rowlb_map[nodeid] = rowlb;
	rowub_map[nodeid] = rowub;

	int e_mc=0;
	int i_mc=0;
	bool e=true; //equality constraint must be at front of list
	for(int i=0;i<rowlb.size();i++)
	{
		if(rowlb[i] == rowub[i]){
			e_mc++;
			assert(e); //equality constraint must be at front of list
		}
		else{
			e = false;
			assert(rowlb[i]<rowub[i]);
			i_mc++;
		}
	}
	assert(i_mc+e_mc == mc);
	e_ncon_map[nodeid] = e_mc;
	i_ncon_map[nodeid] = i_mc;
	PAR_DEBUG("end get_prob_info  - nodeid "<<nodeid<< "(ncon,nvar)=("<<mc<<","<<nv<<") -- e_ncon "<<e_mc<< " i_ncon "<<i_mc);
}

int StructJuMPInput::nFirstStageVars()
{
	PAR_DEBUG("nFirstStageVars");
	std::map<int,int>::iterator it = nvar_map.find(0);
	if(it!=nvar_map.end())
		return  it->second;
	get_prob_info(0);
	return nvar_map[0];
}
int StructJuMPInput::nFirstStageCons(){
	PAR_DEBUG("nFirstStageCons");
	std::map<int,int>::iterator it = ncon_map.find(0);
	if(it!=ncon_map.end())
		return  it->second;
	get_prob_info(0);
	return ncon_map[0];
}
int StructJuMPInput::nSecondStageVars(int scen){
	int nodeid = scen + 1;
	PAR_DEBUG("nSecondStageVars - "<<nodeid);
	std::map<int,int>::iterator it = nvar_map.find(nodeid);
	if(it!=nvar_map.end())
		return  it->second;
	get_prob_info(nodeid);
	return nvar_map[nodeid];
}
int StructJuMPInput::nSecondStageCons(int scen){
	int nodeid = scen + 1;
	PAR_DEBUG("nSecondStageCons - "<<nodeid);
	std::map<int,int>::iterator it = ncon_map.find(nodeid);
	if(it!=ncon_map.end())
		return  it->second;
	get_prob_info(nodeid);
	return ncon_map[nodeid];
}

std::vector<double> StructJuMPInput::getFirstStageColLB(){
	std::map<int, std::vector<double> >::iterator it = collb_map.find(0);
	if(it!=collb_map.end())
		return  it->second;
	get_prob_info(0);
	return collb_map[0];
}
std::vector<double> StructJuMPInput::getFirstStageColUB(){
	std::map<int, std::vector<double> >::iterator it = colub_map.find(0);
	if(it!=colub_map.end())
		return  it->second;
	get_prob_info(0);
	return colub_map[0];
}

std::vector<double> StructJuMPInput::getFirstStageObj(){
	PAR_DEBUG("getFirstStageObj - 0");
	assert(nvar_map.find(0)!=nvar_map.end());
	int nvar = nvar_map[0];
	double x0[nvar];
	std::vector<double> grad(nvar);
	CallBackData data = {prob->userdata,0,0};
	prob->eval_grad_f(x0,x0,&grad[0],&data);
	PAR_DEBUG("end getFirstStageObj - 0");
	return grad;
}
std::vector<std::string> StructJuMPInput::getFirstStageColNames(){
	assert(nvar_map.find(0)!=nvar_map.end());
	int nvar = nvar_map[0];
	std::vector<std::string> cnames(nvar);
	for(int i=0;i<nvar;i++)
	{
		std::ostringstream oss;
		oss<<"x"<<i;
		cnames[i] = oss.str();
	}
	return cnames;
}
std::vector<double> StructJuMPInput::getFirstStageRowLB(){
	std::map<int, std::vector<double> >::iterator it = rowlb_map.find(0);
	if(it!=rowlb_map.end())
		return  it->second;
	get_prob_info(0);
	return rowlb_map[0];
}
std::vector<double> StructJuMPInput::getFirstStageRowUB(){
	std::map<int, std::vector<double> >::iterator it = rowub_map.find(0);
	if(it!=rowub_map.end())
		return  it->second;
	get_prob_info(0);
	return rowub_map[0];
}
std::vector<std::string> StructJuMPInput::getFirstStageRowNames(){
	assert(ncon_map.find(0)!=ncon_map.end());
	int ncon = ncon_map[0];
	std::vector<std::string> cnames(ncon);
	for(int i=0;i<ncon;i++)
	{
		std::ostringstream oss;
		oss<<"c"<<i;
		cnames[i] = oss.str();
	}
	return cnames;
}
bool StructJuMPInput::isFirstStageColInteger(int col){
	return false;
}

std::vector<double> StructJuMPInput::getSecondStageColLB(int scen){
	int nodeid = scen + 1;
	std::map<int, std::vector<double> >::iterator it = collb_map.find(nodeid);
	if(it!=collb_map.end())
		return  it->second;
	get_prob_info(nodeid);
	return collb_map[nodeid];
}
std::vector<double> StructJuMPInput::getSecondStageColUB(int scen){
	int nodeid = scen + 1;
	std::map<int, std::vector<double> >::iterator it = colub_map.find(nodeid);
	if(it!=colub_map.end())
		return it->second;
	get_prob_info(nodeid);
	return colub_map[nodeid];
}
// objective vector, already multiplied by probability
std::vector<double> StructJuMPInput::getSecondStageObj(int scen){
	PAR_DEBUG("getSecondStageObj -  nodeid "<<scen+1);
	int nodeid = scen + 1;
	assert(nvar_map.find(nodeid)!=nvar_map.end());
	assert(nvar_map.find(0)!=nvar_map.end());
	int n0 = nvar_map[0];
	int n1 = nvar_map[nodeid];
	double x0[n0];
	double x1[n1];
	std::vector<double> grad(n1);
	CallBackData data = {prob->userdata,nodeid,nodeid};
	prob->eval_grad_f(x0,x1,&grad[0],&data);
	PAR_DEBUG("end getSecondStageObj - nodeid"<<scen+1);
	return grad;
}
std::vector<std::string> StructJuMPInput::getSecondStageColNames(int scen){
	int nodeid = scen + 1;
	PAR_DEBUG("getSecondStageColNames - nodeid"<<nodeid);
	assert(nvar_map.find(nodeid)!=nvar_map.end());
	assert(nvar_map.find(0)!=nvar_map.end());
	int i0 = nvar_map[0];
	int nvar = nvar_map[nodeid];
	std::vector<std::string> cnames(nvar);
	for(int i=0;i<nvar;i++)
	{
		std::ostringstream oss;
		oss<<"x"<<(i+i0);
		cnames[i] = oss.str();
	}
	return cnames;
}
std::vector<double> StructJuMPInput::getSecondStageRowUB(int scen){
	int nodeid = scen + 1;
	PAR_DEBUG("getSecondStageRowUB - nodeid"<<nodeid);
	std::map<int, std::vector<double> >::iterator it = rowub_map.find(nodeid);
	if(it!=rowub_map.end())
		return  it->second;
	get_prob_info(nodeid);
	return rowub_map[nodeid];
}
std::vector<double> StructJuMPInput::getSecondStageRowLB(int scen){
	int nodeid = scen + 1;
	PAR_DEBUG("getSecondStageRowLB - nodeid"<<nodeid);
	std::map<int, std::vector<double> >::iterator it = rowlb_map.find(nodeid);
	if(it!=rowlb_map.end())
		return  it->second;
	get_prob_info(nodeid);
	return rowlb_map[nodeid];
}
std::vector<std::string> StructJuMPInput::getSecondStageRowNames(int scen){
	int nodeid = scen + 1;
	PAR_DEBUG("getSecondStageRowNames - nodeid"<<nodeid);
	assert(ncon_map.find(0)!=ncon_map.end());
	assert(ncon_map.find(nodeid)!=ncon_map.end());
	int ncon = ncon_map[nodeid];
	int i0 = ncon_map[0];
	std::vector<std::string> cnames(ncon);
	for(int i=0;i<ncon;i++)
	{
		std::ostringstream oss;
		oss<<"c"<<(i+i0);
		cnames[i] = oss.str();
	}
	return cnames;
}
double StructJuMPInput::scenarioProbability(int scen){
	return 1.0/scen;
}
bool StructJuMPInput::isSecondStageColInteger(int scen, int col){
	return false;
}

// returns the column-oriented first-stage constraint matrix (A matrix)
CoinPackedMatrix StructJuMPInput::getFirstStageConstraints(){
	PAR_DEBUG("getFirstStageConstraints - Amat -  ");
	int nvar = nvar_map[0];
	CallBackData cbd = {prob->userdata,0,0};
	std::vector<double> x0(nvar,1.0);
	int e_nz, i_nz;
	prob->eval_jac_g(&x0[0],&x0[0],
			&e_nz,NULL,NULL,NULL,
			&i_nz,NULL,NULL,NULL,&cbd);

	std::vector<int> e_rowidx(e_nz);
	std::vector<int> e_colptr(nvar+1,0);
	std::vector<double> e_elts(e_nz);

	std::vector<int> i_rowidx(i_nz);
	std::vector<int> i_colptr(nvar+1,0);
	std::vector<double> i_elts(i_nz);

	prob->eval_jac_g(&x0[0],&x0[0],
				&e_nz,&e_elts[0],&e_rowidx[0],&e_colptr[0],
				&i_nz,&i_elts[0],&i_rowidx[0],&i_colptr[0],&cbd);


	int e_ncon = e_ncon_map[0];
	int i_ncon = i_ncon_map[0];

	amat.copyOf(true,e_ncon,nvar,e_nz,&e_elts[0],&e_rowidx[0],&e_colptr[0],0);
	CoinPackedMatrix i_amat;
	i_amat.copyOf(true,i_ncon,nvar,i_nz,&i_elts[0],&i_rowidx[0],&i_colptr[0],0);

	amat.bottomAppendPackedMatrix(i_amat);

	assert(amat.getNumCols()==nvar);
	assert(amat.getNumRows()==ncon_map[0]);
	PAR_DEBUG("return getFirstStageConstraints "<<amat.getNumRows()<<" x "<<amat.getNumCols()<<" nz "<<amat.getNumElements());

//	CoinPackedMatrix arow;
//	arow.reverseOrderedCopyOf(amat);
//	PAR_DEBUG("arow "<<arow->getNumRows()<<" "<<arow->getNumCols()<<" "<<arow->getNumElements());

//	CoinPackedMatrix Arow, testmat;
//	int m = 2, n = 3;
//	int colptr[5] ={1,2,4,5,5};
//	int rowidx[4] ={1,1,2,1};
//	double vals[4] = {1.1,1.2,1.3,1.4};
//	testmat.copyOf(true,m,n,4,vals,rowidx,colptr,0);
//	PAR_DEBUG("testmat "<<testmat.getNumRows()<<" "<<testmat.getNumCols()<<" "<<testmat.getNumElements());
//	Arow.reverseOrderedCopyOf( testmat );
//	PAR_DEBUG("testmat "<<testmat.getNumRows()<<" "<<testmat.getNumCols()<<" "<<testmat.getNumElements()<<
//			" Arow -"<<Arow.getNumRows()<<" "<<Arow.getNumCols()<<" "<<Arow.getNumElements());

	return amat;
}
// returns the column-oriented second-stage constraint matrix (W matrix)
CoinPackedMatrix StructJuMPInput::getSecondStageConstraints(int scen){
	int nodeid = scen + 1;
	PAR_DEBUG("getSecondStageConstraints - Wmat -  "<<(scen+1));
	//Bmat + Dmat  = Wmat

	std::map<int,CoinPackedMatrix>::iterator it = wmat_map.find(nodeid);
	if(it!=wmat_map.end()){
		PAR_DEBUG("return (quick) getSecondStageConstraints - Wmat -  "<<it->second.getNumRows()<<" x "<<it->second.getNumCols()
				<<" nz "<<it->second.getNumElements());
		return it->second;
	}

	int nvar = nvar_map[nodeid];

	CallBackData cbd = {prob->userdata,nodeid,nodeid};
	std::vector<double> x0(nvar_map[0],1.0);
	std::vector<double> x1(nvar,1.0);

	int e_nz, i_nz;
	prob->eval_jac_g(&x0[0],&x1[0],
			&e_nz,NULL,NULL,NULL,
			&i_nz,NULL,NULL,NULL,&cbd);

	std::vector<int> e_rowidx(e_nz);
	std::vector<int> e_colptr(nvar+1,0);
	std::vector<double> e_elts(e_nz);

	std::vector<int> i_rowidx(i_nz);
	std::vector<int> i_colptr(nvar+1,0);
	std::vector<double> i_elts(i_nz);

	prob->eval_jac_g(&x0[0],&x1[0],
			&e_nz,&e_elts[0],&e_rowidx[0],&e_colptr[0],
			&i_nz,&i_elts[0],&i_rowidx[0],&i_colptr[0],&cbd);

	int e_ncon = e_ncon_map[nodeid];
	int i_ncon = i_ncon_map[nodeid];

	CoinPackedMatrix wmat;
	CoinPackedMatrix i_wmat;
	wmat.copyOf(true,e_ncon,nvar,e_nz,&e_elts[0],&e_rowidx[0],&e_colptr[0],0);
	i_wmat.copyOf(true,i_ncon,nvar,i_nz,&i_elts[0],&i_rowidx[0],&i_colptr[0],0);

	wmat.bottomAppendPackedMatrix(i_wmat);

	wmat_map[scen] = wmat;
	assert(wmat.getNumCols()==nvar);
	assert(wmat.getNumRows()==ncon_map[nodeid]);
	PAR_DEBUG("return getSecondStageConstraints - Wmat -  "<<wmat.getNumRows()<<" x "<<wmat.getNumCols()<<" nz "<<wmat.getNumElements());
	return wmat;
}
// returns the column-oriented matrix linking the first-stage to the second (T matrix)
CoinPackedMatrix StructJuMPInput::getLinkingConstraints(int scen){
	PAR_DEBUG("getLinkingConstraints - Tmat -  "<<(scen+1));
	int nodeid = scen + 1;

	//Amat + Cmat  = Tmat

	std::map<int,CoinPackedMatrix>::iterator it = tmat_map.find(nodeid);
	if(it!=tmat_map.end()) {
		PAR_DEBUG("return (quick) getLinkingConstraints - Tmat - "<<it->second.getNumRows()<<" x "<<it->second.getNumCols()
				<<" nz "<<it->second.getNumElements());
		return it->second;
	}

	int nvar = nvar_map[0];
	PAR_DEBUG(" nvar "<<nvar);
	CallBackData cbd = {prob->userdata,nodeid,0};
	std::vector<double> x0(nvar_map[0],1.0);
	std::vector<double> x1(nvar_map[nodeid],1.0);

	int e_nz, i_nz;
	prob->eval_jac_g(&x0[0],&x1[0],
			&e_nz,NULL,NULL,NULL,
			&i_nz,NULL,NULL,NULL,&cbd);

	std::vector<int> e_rowidx(e_nz);
	std::vector<int> e_colptr(nvar+1,0);
	std::vector<double> e_elts(e_nz);

	std::vector<int> i_rowidx(i_nz);
	std::vector<int> i_colptr(nvar+1,0);
	std::vector<double> i_elts(i_nz);


	prob->eval_jac_g(&x0[0],&x1[0],
			&e_nz,&e_elts[0],&e_rowidx[0],&e_colptr[0],
			&i_nz,&i_elts[0],&i_rowidx[0],&i_colptr[0],&cbd);

	int e_ncon = e_ncon_map[nodeid];
	int i_ncon = i_ncon_map[nodeid];

	CoinPackedMatrix tmat;
	CoinPackedMatrix i_tmat;
	tmat.copyOf(true,e_ncon,nvar,e_nz,&e_elts[0],&e_rowidx[0],&e_colptr[0],0);
	PAR_DEBUG("e_tmat "<<tmat.getNumElements()<< " exp:"<<e_elts.size());
	assert(e_elts.size() == tmat.getNumElements());
	i_tmat.copyOf(true,i_ncon,nvar,i_nz,&i_elts[0],&i_rowidx[0],&i_colptr[0],0);
	PAR_DEBUG("i_tmat "<<i_tmat.getNumElements());
	assert(i_elts.size() == i_tmat.getNumElements());

	tmat.bottomAppendPackedMatrix(i_tmat);
	PAR_DEBUG("tmat "<<tmat.getNumElements());
	assert(tmat.getNumElements() == (e_nz+i_nz));

	tmat_map[nodeid] = tmat;
	assert(tmat.getNumCols()==nvar);
	assert(tmat.getNumRows()==ncon_map[nodeid]);
	PAR_DEBUG("return getLinkingConstraints - Tmat - "<<tmat.getNumRows()<<" x "<<tmat.getNumCols()<<" nz "<<tmat.getNumElements());
	return tmat;
}

CoinPackedMatrix StructJuMPInput::getFirstStageHessian(){
	PAR_DEBUG("getFirstStageHessian - ");
	assert(nvar_map.find(0)!=nvar_map.end());

	int nvar = nvar_map[0];
	int ncon = ncon_map[0];
	std::vector<double> x0(nvar,1.0);
	std::vector<double> lam(ncon,1.0);

	int nz;
	CallBackData cbd = {prob->userdata,0,0};
	prob->eval_h(&x0[0],&x0[0],&lam[0],&nz, NULL,NULL,NULL,&cbd);
	PAR_DEBUG("getFirstStageHessian - nz "<<nz);

	std::vector<int> rowidx(nz);
	std::vector<int> colptr(nvar+1,0);
	std::vector<double> elts(nz);

	PAR_DEBUG("getFirstStageHessian - before value");
	prob->eval_h(&x0[0],&x0[0],&lam[0],&nz,&elts[0],&rowidx[0],&colptr[0],&cbd);
	PAR_DEBUG("getFirstStageHessian - after value");

	qamat.copyOf(true,nvar,nvar,nz,&elts[0],&rowidx[0],&colptr[0],0);

	PAR_DEBUG("return getFirstStageHessian - Qamat - "<<qamat.getNumRows()<<" x "<<qamat.getNumCols()
			<<" nz "<<qamat.getNumElements());
	return qamat;
}
// Q_i
CoinPackedMatrix StructJuMPInput::getSecondStageHessian(int scen){
	PAR_DEBUG("getSecondStageHessian - "<<scen+1);
	int nodeid = scen + 1;
	assert(nvar_map.find(nodeid)!=nvar_map.end());
	std::map<int,CoinPackedMatrix>::iterator it = qwmat_map.find(nodeid);
	if(it!=qwmat_map.end()){
		PAR_DEBUG("return (quick) getSecondStageHessian - Qwmat - "<<it->second.getNumRows()<<" x "<<it->second.getNumCols()
				<<" nz "<<it->second.getNumElements());
		return it->second;
	}

	int n0 = nvar_map[0];
	int nvar = nvar_map[nodeid];
	int ncon = ncon_map[nodeid];
	std::vector<double> x0(n0,1.0);
	std::vector<double> x1(nvar,1.0);
	std::vector<double> lam(ncon,1.0);

	int nz;
	CallBackData cbd = {prob->userdata,nodeid,nodeid};
	prob->eval_h(&x0[0],&x1[0],&lam[0],&nz, NULL,NULL,NULL,&cbd);

	std::vector<int> rowidx(nz);
	std::vector<int> colptr(nvar+1,0);
	std::vector<double> elts(nz);

	prob->eval_h(&x0[0],&x1[0],&lam[0],&nz,&elts[0],&rowidx[0],&colptr[0],&cbd);

	CoinPackedMatrix qwmat;
	qwmat.copyOf(true,nvar,nvar,nz,&elts[0],&rowidx[0],&colptr[0],0);

	qwmat_map[nodeid] = qwmat;
	PAR_DEBUG("return getSecondStageHessian - Qwmat - "<<qwmat.getNumRows()<<" x "<<qwmat.getNumCols()
					<<" nz "<<qwmat.getNumElements());
	return qwmat;
}

// column-oriented, \hat Q_i
// Note: this has the second-stage variables on the rows and first-stage on the columns
CoinPackedMatrix StructJuMPInput::getSecondStageCrossHessian(int scen){
	PAR_DEBUG("getSecondStageCrossHessian - "<<scen+1);
	int nodeid = scen + 1;
	assert(nvar_map.find(nodeid)!=nvar_map.end());
	assert(nvar_map.find(0)!=nvar_map.end());
	std::map<int,CoinPackedMatrix>::iterator it = qtmat_map.find(nodeid);
	if(it!=qtmat_map.end()){
		PAR_DEBUG("return (quick) getSecondStageCrossHessian - Qtmat - "<<it->second.getNumRows()<<" x "<<it->second.getNumCols()
				<<" nz "<<it->second.getNumElements());
		return it->second;
	}

	int n0 = nvar_map[0];
	int n1 = nvar_map[nodeid];
	int ncon = ncon_map[nodeid];
	std::vector<double> x0(n0,1.0);
	std::vector<double> x1(n1,1.0);
	std::vector<double> lam(ncon,1.0);


	int nz;
	CallBackData cbd = {prob->userdata,0,nodeid};
	prob->eval_h(&x0[0],&x1[0],&lam[0],&nz, NULL,NULL,NULL,&cbd);
	PAR_DEBUG(" nz "<<nz<< " - (0,"<<nodeid<<") - "<<n0<<" - "<<ncon);
	std::vector<int> rowidx(nz,0);
	std::vector<int> colptr(n0+1,0);
	std::vector<double> elts(nz,0.0);

	prob->eval_h(&x0[0],&x1[0],&lam[0],&nz,&elts[0],&rowidx[0],&colptr[0],&cbd);

	CoinPackedMatrix qtmat(true,n1,n0,nz,&elts[0],&rowidx[0],&colptr[0],0);

	qtmat_map[nodeid] = qtmat;
	PAR_DEBUG("return getSecondStageCrossHessian - Qtmat - "<<qtmat.getNumRows()<<" x "<<qtmat.getNumCols()
						<<" nz "<<qtmat.getNumElements());
	return qtmat;
}


// some problem characteristics that could be helpful to know

// all scenarios have the same number of variables and constraints
bool StructJuMPInput::scenarioDimensionsEqual(){
	return false;
}
// constraint (and hessian) matrices are identical for each scenario,
// column and row bounds and objective are allowed to vary
bool StructJuMPInput::onlyBoundsVary(){
	return false;
}
// all scenarios equally likely
bool StructJuMPInput::allProbabilitiesEqual(){
	return true;
}
// all second-stage variables continuous
bool StructJuMPInput::continuousRecourse(){
	return true;
}
