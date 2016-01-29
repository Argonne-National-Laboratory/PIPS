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
	return prob->nnodes;
}

void StructJuMPInput::get_prob_info(int nodeid)
{
	PAR_DEBUG("get_prob_info - prob_info - "<<nodeid);
	int nv = 0;
	int mc = 0;

	CallBackData data={prob->userdata,nodeid,nodeid};
//	PAR_DEBUG("data -"<<&data);
	PAR_DEBUG("get_prob_info callback data ptr - "<<&data);
	PAR_DEBUG("get_prob_info userdata ptr - "<<prob->userdata);
	prob->prob_info(&nv,NULL,NULL,&mc,NULL,NULL,&data);
	nvar_map[nodeid] = nv;
	ncon_map[nodeid] = mc;

	std::vector<double> collb(nv);
	std::vector<double> colub(nv);
	std::vector<double> rowlb(mc);
	std::vector<double> rowub(mc);
	prob->prob_info(&nv, &collb[0], &colub[0], &mc, &rowlb[0], &rowub[0], &data);

	collb_map[nodeid] = collb;
	colub_map[nodeid] = colub;
	rowlb_map[nodeid] = rowlb;
	rowub_map[nodeid] = rowub;
	PAR_DEBUG("get_prob_info - insert - "<<nodeid);
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
	PAR_DEBUG("nSecondStageVars");
	std::map<int,int>::iterator it = nvar_map.find(nodeid);
	if(it!=nvar_map.end())
		return  it->second;
	get_prob_info(nodeid);
	return nvar_map[nodeid];
}
int StructJuMPInput::nSecondStageCons(int scen){
	int nodeid = scen + 1;
	PAR_DEBUG("nSecondStageCons");
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
	std::map<int, std::vector<double> >::iterator it = rowub_map.find(nodeid);
	if(it!=rowub_map.end())
		return  it->second;
	get_prob_info(nodeid);
	return rowub_map[nodeid];
}
std::vector<double> StructJuMPInput::getSecondStageRowLB(int scen){
	int nodeid = scen + 1;
	std::map<int, std::vector<double> >::iterator it = rowlb_map.find(nodeid);
	if(it!=rowlb_map.end())
		return  it->second;
	get_prob_info(nodeid);
	return rowlb_map[nodeid];
}
std::vector<std::string> StructJuMPInput::getSecondStageRowNames(int scen){
	int nodeid = scen + 1;
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
	int nvar = nvar_map[0];
	int ncon = ncon_map[0];
	std::vector<int> rowindx(0);
	std::vector<int> colptr(nvar+1,0);
	std::vector<double> elts(0);
	amat.copyOf(true,ncon,nvar,0,&elts[0],&rowindx[0],&colptr[0],0);
	assert(amat.getNumCols()==nvar);
	assert(amat.getNumRows()==ncon);
	return amat;
}
// returns the column-oriented second-stage constraint matrix (W matrix)
CoinPackedMatrix StructJuMPInput::getSecondStageConstraints(int scen){
	int nodeid = scen + 1;
	std::map<int,CoinPackedMatrix>::iterator it = wmat_map.find(nodeid);
	if(it!=wmat_map.end())
		return it->second;

	int nvar = nvar_map[nodeid];
	int ncon = ncon_map[nodeid];
	std::vector<int> rowidx(0);
	std::vector<int> colptr(nvar+1,0);
	std::vector<double> elts(0);
	CoinPackedMatrix wmat;
	wmat.copyOf(true,ncon,nvar,0,&elts[0],&rowidx[0],&colptr[0],0);
	wmat_map[scen] = wmat;
	assert(wmat.getNumCols()==nvar);
	assert(wmat.getNumRows()==ncon);
	return wmat;
}
// returns the column-oriented matrix linking the first-stage to the second (T matrix)
CoinPackedMatrix StructJuMPInput::getLinkingConstraints(int scen){
	int nodeid = scen + 1;
	std::map<int,CoinPackedMatrix>::iterator it = tmat_map.find(nodeid);
	if(it!=tmat_map.end())
		return it->second;

	int nvar = nvar_map[0];
	int ncon = ncon_map[nodeid];
	std::vector<int> rowidx(0);
	std::vector<int> colptr(nvar+1,0);
	std::vector<double> elts(0);
	CoinPackedMatrix tmat;
	tmat.copyOf(true,ncon,nvar,0,&elts[0],&rowidx[0],&colptr[0],0);
	tmat_map[nodeid] = tmat;
	assert(tmat.getNumCols()==nvar);
	assert(tmat.getNumRows()==ncon);
	return tmat;
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
