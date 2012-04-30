#include "ucRollingModel.hpp"
#include "CbcLagrangeSolver.hpp"
#include "ClpRecourseSolver.hpp"
#include "conicBundleDriver.hpp"
#include <boost/scoped_ptr.hpp>

using namespace std;
using boost::scoped_ptr;

int main(int argc, char **argv) {

        
        MPI_Init(&argc, &argv);

        int mype;
        MPI_Comm_rank(MPI_COMM_WORLD,&mype);

        if (argc != 3) {
		printf("Usage: %s [num scenarios] [time horizon]\n",argv[0]);
                return 1;
        }
	
	int nscen = atoi(argv[1]);
	int T = atoi(argv[2]);
	ucRollingModel model("../../../apps/unitcommitment_rolling/Illinois",nscen,0,T);

        conicBundleDriver<CbcLagrangeSolver,ClpRecourseSolver>(model);

        MPI_Finalize();

        return 0;

}
