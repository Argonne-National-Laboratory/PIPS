/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef AMPLNLINTERFACE_HPP
#define AMPLNLINTERFACE_HPP

#include <map>
#include <vector>
#include <string>


struct ASL_pfgh;
struct SufDecl;
struct SufDesc;
class  AmplSuffix;

enum { kMinimize = 0, kMaximize = 1 } ;

enum Suffix_NumType
{
  Suffix_Int,
  Suffix_Double
};

enum Suffix_Type
{
  Suffix_Var,
  Suffix_Con
};


class AmplData_NL {

public:
	std::vector<double> collb, colub, rowlb, rowub, objGrad;
	std::vector<int> VarsStatus, ConsStatus;
	std::vector<std::string> rownames, colnames;
	bool didLoad;
	ASL_pfgh *locASL;
	int nzH;
	
	AmplData_NL() : didLoad(false),locASL(NULL), nzH(-1) {}
	ASL_pfgh* initASL(char *nlFileName[], AmplSuffix* amplSuffix=NULL);
	void initialize(int nvar, int ncons);

	int Ampl_nnz_Hessian_Tri();

};



class AmplSuffix
{
public:
  AmplSuffix();
  ~AmplSuffix();

  void DefineSuffix(std::string suffix_name, Suffix_Type suffix_type, Suffix_NumType suffix_numtype)
  {
	Suf_Idx.push_back(suffix_name);
	Suf_NumType.push_back(suffix_numtype);
	Suf_Type.push_back(suffix_type);
  }

  int* GetSuffixVal_Int(ASL_pfgh* asl, std::string suffix_string, Suffix_Type source);

  double* GetSuffixVal_Double(ASL_pfgh* asl,std::string suffix_string, Suffix_Type source);

  void InitSpaceForSuffixes(ASL_pfgh* asl);

private:
  /**@name Default Compiler Generated Methods
   * (Hidden to avoid implicit creation/calling).
   * These methods are not implemented and 
   * we do not want the compiler to implement
   * them for us, so we declare them private
   * and do not define them. This ensures that
   * they will not be implicitly created/called. */
  //@{
  /** Default Constructor */
  //AmplSuffixHandler();

  /** Copy Constructor */
  AmplSuffix(const AmplSuffix&);

  /** Overloaded Equals Operator */
  void operator=(const AmplSuffix&);
  //@}

  SufDecl* suftab_;

  std::vector<std::string> 		Suf_Idx;
  std::vector<Suffix_NumType> 	Suf_NumType;
  std::vector<Suffix_Type> 		Suf_Type;

};



#endif
