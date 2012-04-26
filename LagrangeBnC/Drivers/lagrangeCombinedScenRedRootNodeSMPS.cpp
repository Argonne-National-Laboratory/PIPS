
#include "SMPSInput.hpp"
#include "CbcLagrangeSolver.hpp"
#include "CbcRecourseSolver.hpp"
#include "lagrangeCombinedScenRedRootNode.hpp"
#include <boost/scoped_ptr.hpp>

using namespace std;
using boost::scoped_ptr;

int main(int argc, char **argv) {

	
	MPI_Init(&argc, &argv);

	int mype;
	MPI_Comm_rank(MPI_COMM_WORLD,&mype);

	if (argc != 3) {
		if (mype == 0) printf("Usage: %s [SMPS root name] [number per subproblem]\n",argv[0]);
		return 1;
	}

	string smpsrootname(argv[1]);
	int nper = atoi(argv[2]);

	SMPSInput input(smpsrootname+".cor",smpsrootname+".tim",smpsrootname+".sto");

	//lagrangeCombinedScenRedRootNode<CbcLagrangeSolver,CbcRecourseSolver>(input,nper,false);
	lagrangeCombinedScenRedRootNode<CbcLagrangeSolver,CbcRecourseSolver>(input,nper,true);

	MPI_Finalize();

	return 0;

}
