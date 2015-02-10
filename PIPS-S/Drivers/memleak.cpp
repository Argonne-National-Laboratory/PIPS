#include "BAData.hpp"
#include "rawInput.hpp"
#include "PIPSSInterface.hpp"
#include <boost/scoped_ptr.hpp>
#include <cstdlib>
#include <vector>

//#include "BAVector.hpp"

using boost::scoped_ptr; // replace with unique_ptr for C++11
using namespace std;

int main(int argc, char **argv) {
	
  cout << "Start" << endl;
  MPI_Init(&argc, &argv);
  
  int mype;
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  
  
  string datarootname("/sandbox/petra/work/temp/rawinput/ssndata/problemdata");
  int nscen = 4;
  
  scoped_ptr<rawInput> s(new rawInput(datarootname,nscen));

  BAContext ctx(MPI_COMM_WORLD);
  ctx.initializeAssignment(nscen);

  BADimensions dims(*s, ctx);

  cout << " ctx and dims created" << endl;

  vector<denseBAVector> v1,v2;

  denseBAVector baVec1Copy;

  {   
    denseBAVector baVec1,baVec2;
    baVec1.allocate(dims, ctx, PrimalVector);

    v1.push_back(baVec1);

    baVec2.allocate(dims,ctx,PrimalVector);
    v1.push_back(baVec2);
    
    
    cout << "vector1 has " << v1.size() << " elements\n";
    denseBAVector baVec3 = v1[0];

    //testing assignement operator
    baVec1Copy = v1[1];

    baVec2 = baVec1;

    v2.push_back(baVec3);
    v2.push_back(baVec1Copy);
    v2.push_back(baVec2);
  }
  denseBAVector baVec;
  baVec.allocate(dims, ctx, PrimalVector);  
  v2.push_back(baVec);
  
  cout << "vector2 has " << v2.size() << " elements\n";

  MPI_Finalize();
  cout << "driver exiting ..." << endl;
  return 0;
}

