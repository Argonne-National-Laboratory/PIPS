#include "StructJuMPInput.h"

#include <string>
#include <iostream>
#include <sstream>

StructJuMPInput::StructJuMPInput(PipsNlpProblemStruct* p)
{
	this->prob = p;
	useInputDate = 1;
	datarootname = "StructJuMP";
}
StructJuMPInput::~StructJuMPInput() {

}

int StructJuMPInput::nScenarios() {
	return prob->nnodes;
}

void StructJuMPInput::get_prob_info(int nodeid)
{
	std::cout<<"call jump - prob_info - "<<nodeid<<std::endl;
	int nv;
	int mc;
	double* col_ub;
	double* col_lb;
	double* row_lb;
	double* row_ub;

	CallBackData data={nodeid,nodeid};
	prob->prob_info(&nv, col_lb, col_ub, &mc, row_lb, row_ub, &data);
	nvar_map[nodeid] = nv;
	ncon_map[nodeid] = mc;
	std::vector<double> collb; collb.assign(col_lb,col_lb+nv);
	std::vector<double> colub; colub.assign(col_ub,col_ub+nv);
	std::vector<double> rowlb; rowlb.assign(row_lb,row_lb+mc);
	std::vector<double> rowub; colub.assign(col_ub,col_ub+mc);
	collb_map[nodeid] = collb;
	colub_map[nodeid] = colub;
	rowlb_map[nodeid] = rowlb;
	rowub_map[nodeid] = rowub;
}

int StructJuMPInput::nFirstStageVars()
{
	std::cout<<"nFirstStageVars"<<std::endl;
	std::map<int,int>::iterator it = nvar_map.find(0);
	if(it!=nvar_map.end())
		return  it->second;
	get_prob_info(0);
	return nvar_map[0];
}
int StructJuMPInput::nFirstStageCons(){
	std::cout<<"nFirstStageCons"<<std::endl;
	std::map<int,int>::iterator it = ncon_map.find(0);
	if(it!=ncon_map.end())
		return  it->second;
	get_prob_info(0);
	return ncon_map[0];
}
int StructJuMPInput::nSecondStageVars(int scen){
	std::cout<<"nSecondStageVars"<<std::endl;
	std::map<int,int>::iterator it = nvar_map.find(scen);
	if(it!=nvar_map.end())
		return  it->second;
	get_prob_info(scen);
	return nvar_map[scen];
}
int StructJuMPInput::nSecondStageCons(int scen){
	std::cout<<"nSecondStageCons"<<std::endl;
	std::map<int,int>::iterator it = ncon_map.find(scen);
	if(it!=ncon_map.end())
		return  it->second;
	get_prob_info(scen);
	return ncon_map[scen];
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
	assert(nvar_map.find(0)!=nvar_map.end());
	int nvar = nvar_map[0];
	double x0[nvar];
	std::vector<double> grad(nvar);
	CallBackData data = {0,0};
	prob->eval_grad_f(x0,NULL,&grad[0],&data);
	return grad;
}
std::vector<std::string> StructJuMPInput::getFirstStageColNames(){
	assert(nvar_map.find(0)!=nvar_map.end());
	int nvar = nvar_map[0];
	std::vector<std::string> cnames(nvar);
	for(int i=0;i<nvar;i++)
	{
		std::ostringstream oss;
		oss<<"x"<<i<<std::endl;
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
		oss<<"c"<<i<<std::endl;
		cnames[i] = oss.str();
	}
	return cnames;
}
bool StructJuMPInput::isFirstStageColInteger(int col){
	return false;
}

std::vector<double> StructJuMPInput::getSecondStageColLB(int scen){
	std::map<int, std::vector<double> >::iterator it = collb_map.find(scen);
	if(it!=collb_map.end())
		return  it->second;
	get_prob_info(scen);
	return collb_map[scen];
}
std::vector<double> StructJuMPInput::getSecondStageColUB(int scen){
	std::map<int, std::vector<double> >::iterator it = colub_map.find(scen);
	if(it!=colub_map.end())
		return it->second;
	get_prob_info(scen);
	return colub_map[scen];
}
// objective vector, already multiplied by probability
std::vector<double> StructJuMPInput::getSecondStageObj(int scen){
	assert(nvar_map.find(scen)!=nvar_map.end());
	assert(nvar_map.find(0)!=nvar_map.end());
	int n0 = nvar_map[0];
	int n1 = nvar_map[scen];
	double x0[n0];
	double x1[n1];
	std::vector<double> grad(n1);
	CallBackData data = {0,scen};
	prob->eval_grad_f(x0,x1,&grad[0],&data);
	return grad;
}
std::vector<std::string> StructJuMPInput::getSecondStageColNames(int scen){
	assert(nvar_map.find(scen)!=nvar_map.end());
	assert(nvar_map.find(0)!=nvar_map.end());
	int i0 = nvar_map[0];
	int nvar = nvar_map[scen];
	std::vector<std::string> cnames(nvar);
	for(int i=0;i<nvar;i++)
	{
		std::ostringstream oss;
		oss<<"x"<<(i+i0)<<std::endl;
		cnames[i] = oss.str();
	}
	return cnames;
}
std::vector<double> StructJuMPInput::getSecondStageRowUB(int scen){
	std::map<int, std::vector<double> >::iterator it = rowub_map.find(scen);
	if(it!=rowub_map.end())
		return  it->second;
	get_prob_info(scen);
	return rowub_map[scen];
}
std::vector<double> StructJuMPInput::getSecondStageRowLB(int scen){
	std::map<int, std::vector<double> >::iterator it = rowlb_map.find(scen);
	if(it!=rowlb_map.end())
		return  it->second;
	get_prob_info(scen);
	return rowlb_map[scen];
}
std::vector<std::string> StructJuMPInput::getSecondStageRowNames(int scen){
	assert(ncon_map.find(0)!=ncon_map.end());
	assert(ncon_map.find(scen)!=ncon_map.end());
	int ncon = ncon_map[scen];
	int i0 = ncon_map[0];
	std::vector<std::string> cnames(ncon);
	for(int i=0;i<ncon;i++)
	{
		std::ostringstream oss;
		oss<<"c"<<(i+i0)<<std::endl;
		cnames[i] = oss.str();
	}
	return cnames;
}
double StructJuMPInput::scenarioProbability(int scen){
	return 1.0/prob->nnodes;
}
bool StructJuMPInput::isSecondStageColInteger(int scen, int col){
	return false;
}

// returns the column-oriented first-stage constraint matrix (A matrix)
CoinPackedMatrix StructJuMPInput::getFirstStageConstraints(){
	int nvar = nvar_map[0];
	int ncon = ncon_map[0];
	std::vector<int> elements(nvar+1,0);
	amat.copyOf(true,ncon,nvar,0,0,0,&elements[0],0);
	assert(amat.getNumCols()==nvar);
	assert(amat.getNumRows()==ncon);
	return amat;
}
// returns the column-oriented second-stage constraint matrix (W matrix)
CoinPackedMatrix StructJuMPInput::getSecondStageConstraints(int scen){
	int nvar = nvar_map[scen];
	int ncon = ncon_map[scen];
	std::vector<int> elements(nvar+1,0);
	wmat.copyOf(true,ncon,nvar,0,0,0,&elements[0],0);
	assert(wmat.getNumCols()==nvar);
	assert(wmat.getNumRows()==ncon);
	return wmat;
}
// returns the column-oriented matrix linking the first-stage to the second (T matrix)
CoinPackedMatrix StructJuMPInput::getLinkingConstraints(int scen){
	int nvar = nvar_map[0];
	int ncon = ncon_map[scen];
	std::vector<int> elements(nvar+1,0);
	tmat.copyOf(true,ncon,nvar,0,0,0,&elements[0],0);
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
