/* PIPS
   Authors: Miles Lubin and Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysRootAugEmtlSym.h"
#include "QpGenStochData.h"
#include "EmtlDenSymMatrix.h"
#include "EmtlStochSymIndefSolver.h"
#include "EmtlSymPSDSolver.h"
#include "EmtlVector.h"
#include <cmath>

#ifdef STOCH_TESTING
extern double g_iterNumber;
#endif

sLinsysRootAugEmtlSym::sLinsysRootAugEmtlSym(sFactory * factory_, QpGenStochData * prob_, const EmtlContext &ctx_)
  : sLinsysRootAugEmtl(factory_, prob_, ctx_)
{
  // for some reason createSolver called from sLinsysRootAugEmtl
  // constructor doesn't map to this class's version
  delete solver;
  solver = createSolver(prob_, kkt);
}

sLinsysRootAugEmtlSym::sLinsysRootAugEmtlSym(sFactory* factory_,
			       QpGenStochData* prob_,
			       OoqpVector* dd_, 
			       OoqpVector* dq_,
			       OoqpVector* nomegaInv_,
			       OoqpVector* rhs_,
			       const EmtlContext &ctx_)
  : sLinsysRootAugEmtl(factory_, prob_, dd_, dq_, nomegaInv_, rhs_, ctx_)
{ 

}

sLinsysRootAugEmtlSym::~sLinsysRootAugEmtlSym()
{
}


DoubleLinearSolver*
sLinsysRootAugEmtlSym::createSolver(QpGenStochData* prob, SymMatrix* kktmat_)
{

  EmtlDenSymMatrix* kktmat = dynamic_cast<EmtlDenSymMatrix*>(kktmat_);
  return new EmtlStochSymIndefSolver(kktmat, locnx);
  //assert(locnx == kktmat->size());
  //return new EmtlSymPSDSolver(kktmat);
}


// given col in processor grid, and starting column of buffer,
// see which is the first column in the buffer that belongs to this processor
// also used here for rows
static inline int firstcol(const int mycol,const int startcol,const int npcol)
{
  return (mycol - startcol + (startcol/npcol+1)*npcol)%npcol;
}

// figure out the largest number of columns we can calculate
// that will fit in the buffers
// we calculate from startcol inclusive to endcol exclusive
// this involves finding the zero of a quadratic equation:
// (elements in lower diagonal from startcol to endcol) - max_elts
static int calcendcol(int startcol, int n, int max_elts)
{
  double discrim = pow(n+0.5,2.)+2.*(.5*pow((double)startcol,2.)-(n+.5)*startcol-max_elts);
  double ans;
  if (discrim > 0) {
    ans = static_cast<int>((n+.5)-sqrt(discrim));
  } else {
    ans = n;
  }
  ans = MIN(ans,n);
  assert(ans > startcol);
  return ans;
}
 
static int countcols(int startcol, int endcol, int n)
{
  return -.5*pow((double)endcol,2)+(n+.5)*endcol+.5*pow((double)startcol,2)-(n+.5)*startcol;
}

// upper bound on how elements will be sent to each processor
static inline int recvUB(int startcol, int endcol, int nprow, int npcol, int n) {
  int count = 0;
  for (int i = startcol; i < endcol; i += npcol) {
    count += utilities::MaxLocalLength(n - i, nprow);
  }
  return count;
}

// figure out the largest number of columns we can calculate
// that will fit in a buffer that is padded such that each
// processor recvs the same number of elements
static inline int calcendcol2(int startcol, int n, int max_elts, int nprow, int npcol,
                              int &UB) {
  int guess = calcendcol(startcol, n, max_elts);
  // shouldn't need more than one iteration here
  while ((UB=recvUB(startcol, guess, nprow, npcol, n))*nprow*npcol > max_elts) {
    guess--;
  }
  return guess;

}
#ifdef __bg__
/// BGP memory -- delete after done debugging
#include <spi/kernel_interface.h>
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>
#include <malloc.h>

void
getMemSize(long long *mem)
{
  long long total;
  int node_config;
  _BGP_Personality_t personality;

  Kernel_GetPersonality(&personality, sizeof(personality));
  total = BGP_Personality_DDRSizeMB(&personality);

  node_config  = BGP_Personality_processConfig(&personality);
  if (node_config == _BGP_PERS_PROCESSCONFIG_VNM) total /= 4;
  else if (node_config == _BGP_PERS_PROCESSCONFIG_2x2) total /= 2;
  total *= 1024*1024;

  *mem = total;
}

/* gives total memory used 
 *  NOTE: this does not account for static data, only heap
 */
void
getUsedMem(long long *mem)
{
  long long alloc;
  struct mallinfo m;

  m = mallinfo();
  alloc = m.hblkhd + m.uordblks;

  *mem = alloc;
}

/* gives available memory 
 *  NOTE: this does not account for static data
 */
void
getFreeMem(long long *mem)
{
  long long total, alloc;

  getMemSize(&total);
  getUsedMem(&alloc);

  *mem = total - alloc;
}
void memstat() {
  long long m1,m2,m3;
  getMemSize(&m1);
  getUsedMem(&m2);
  getFreeMem(&m3);
  cout << m1 << " " << m2 << " " << m3 << endl;
}
///
#endif


const double MAX_MB_FOR_SEND_BUFFERS = 200;

void sLinsysRootAugEmtlSym::factor2(QpGenStochData *prob, Variables *vars)
{
  EmtlDenSymMatrix& kktd = dynamic_cast<EmtlDenSymMatrix&>(*kkt);
  int nxP = locnx;
  
  const int ELTS_PER_BUFFER = 1048576*MAX_MB_FOR_SEND_BUFFERS/
                                (2*sizeof(double));

  initializeKKT(prob, vars);
  
  // initial guess at max number of elements we'll need to receive
  int recvsize = static_cast<int>(1.05*ELTS_PER_BUFFER/(ctx.nprow()*ctx.npcol()));


  double *recvbuffer  = new double[recvsize];
  double *sendbuffer  = new double[ELTS_PER_BUFFER];
  double *schurbuffer = new double[ELTS_PER_BUFFER];
  int *recvcounts_real = new int[ctx.nprocs()];
  int *recvcounts = new int[ctx.nprocs()]; // with padding

  memset(recvcounts, 0, ctx.nprocs()*sizeof(int));

  bool memreport = true;

  int *nr_counts = new int[ctx.nprow()];
  for (int i = 0; i < ctx.nprow(); i++) {
    nr_counts[i] = utilities::LocalLength(nxP, i, ctx.nprow());
  }
  
  // First tell children to factorize.
  for(unsigned int c=0; c<children.size(); c++) {
    children[c]->factor2(prob->children[c], vars);
  }
 
  int endcol;
  for (int startcol = 0; startcol < nxP; startcol = endcol) {
    int UB;
    endcol = calcendcol2(startcol, nxP, ELTS_PER_BUFFER, ctx.nprow(), ctx.npcol(), UB);
    int numcols = endcol-startcol;
    /*if (ctx.mype() == 0) {
      printf("startcol: %d endcol: %d UB: %d\n", startcol, endcol, UB);
    }*/
    memset(schurbuffer, 0, ELTS_PER_BUFFER*sizeof(double));
    memset(recvcounts_real, 0, ctx.nprocs()*sizeof(int));
    for (int i = 0; i < ctx.npcol()*ctx.nprow(); i++) {
      recvcounts[i] = UB;
    }

    for (unsigned int c=0; c<children.size(); c++) {
      if(children[c]->mpiComm == MPI_COMM_NULL)
      	continue;
    
      children[c]->stochNode->resMon.recFactTmChildren_start();    
      //---------------------------------------------
      children[c]->symAddColsToDenseSchurCompl(prob->children[c], schurbuffer, startcol, endcol);
      //---------------------------------------------
      children[c]->stochNode->resMon.recFactTmChildren_stop();    
    
    }

    // only to improve timing of reduce 
    //MPI_Barrier(ctx.comm());

    //printf("pe %d got columns\n", ctx.mype());
    stochNode->resMon.recReduceTmLocal_start(); 
    // now the fun part, first rearrange the elements so that
    // we can call reducescatter
    // that is, all elements belonging to the first proc go first, etc
    
    // create a map from columns to location of first element in schurbuffer
    int *start = new int[numcols];
    start[0] = 0;
    for (int i = 1; i < numcols; i++) {
      start[i] = start[i-1] + nxP-(startcol+i-1);
    }


    int pcol, prow, pe, offset, j, i;
    const int npcol = ctx.npcol(), nprow = ctx.nprow();
    #pragma omp parallel for schedule(dynamic) private(pcol, prow, pe, offset, j, i) COLLAPSE(2)
    for (pcol = 0; pcol < npcol; pcol++) {
      for (prow = 0; prow < nprow; prow++) {
        pe = ctx.get_pnum(prow, pcol);
        offset = pe*UB;
        for(j = firstcol(pcol,startcol,ctx.npcol()); 
              j < numcols; j+= ctx.npcol()) {
          //printf("pe %d for %d %d at col %d\n",ctx.mype(),prow,pcol,j);
          for (i = firstcol(prow,j+startcol,ctx.nprow()); 
                      i < nxP-(j+startcol); i += ctx.nprow()) {
            sendbuffer[offset++] = schurbuffer[start[j]+i];
            recvcounts_real[pe]++;
          }
        }
        //printf("pe %d loaded buffer for %d %d\n",ctx.mype(),prow,pcol); 
      }
    }
    //assert(desti == countcols(startcol, endcol, nxP));
    //printf("pe %d loaded send buffer\n", ctx.mype());
    if (UB > recvsize) {
      // this should happen very rarely
      printf("proc %d needed to resize from %d to %d\n", ctx.mype(), recvsize, recvcounts[ctx.mype()]);
      delete [] recvbuffer;
      recvbuffer = new double[UB];
      recvsize = UB;
    }

    double rtime = MPI_Wtime();

    stochNode->resMon.recReduceScatterTmLocal_start();
    MPI_Reduce_scatter(sendbuffer, recvbuffer, recvcounts, MPI_DOUBLE, 
      MPI_SUM, ctx.comm());
    stochNode->resMon.recReduceScatterTmLocal_stop();
   
		#ifdef __bg__
		if (memreport && ctx.mype() == 0) {
      printf("memory report after redscat: ");
      memstat();
		  //memreport = false;
    }
    #endif

    rtime = MPI_Wtime()-rtime;
    //printf("Reduce-scatter took %f on proc %d, which got %d\n", rtime, ctx.mype(), recvcounts[ctx.mype()]);


    // now unpack on the local processor
    // each column is already in continuous memory
    int myrecv = recvcounts_real[ctx.mype()];
    if ( myrecv > 0 ) {
      assert(myrecv <= UB);
      int fcol = firstcol(ctx.mycol(),startcol,ctx.npcol());
      int lcol = (startcol+fcol)
                    /ctx.npcol();
      int buffi = 0, j, endj, nlrows;
      // precalculate map from columns to position in buffer 
      // so we can parallelize
      /*
      std::vector<int> buffis, lrows;
      buffis.reserve(numcols/(ctx.nprow()*ctx.npcol()));
      lrows.reserve(numcols/(ctx.nprow()*ctx.npcol()));
      buffis.push_back(0);
      for (j = 0; buffi < myrecv; j++) {
        int realcol = startcol+fcol+j*ctx.npcol();
        // first global row we're getting in this column
        int realrow = realcol+firstcol(ctx.myrow(),realcol,ctx.nprow());
        int lrow = realrow/ctx.nprow();
        int nlrows = nr_counts[ctx.myrow()]-lrow;
        buffi += nlrows;
        buffis.push_back(buffi);
        lrows.push_back(lrow);
      }
      endj = j;
      assert(buffi == myrecv);
      
      #pragma omp parallel for schedule(dynamic) private(nlrows)
      for (j = 0; j < endj; j++) {
        nlrows = nr_counts[ctx.myrow()]-lrows[j];

        memcpy(kktd.mat->data+(j+lcol)*kktd.getNR()+lrows[j],
                recvbuffer+buffis[j], nlrows*sizeof(double));
      }*/
        
      // SERIAL LOOP
      for (int j = 0; buffi < myrecv; j++) {
        int realcol = startcol+fcol+j*ctx.npcol();
        // first global row we're getting in this column
        int realrow = realcol+firstcol(ctx.myrow(),realcol,ctx.nprow());
        int lrow = realrow/ctx.nprow();
        int nlrows = nr_counts[ctx.myrow()]-lrow;
        //printf("%d %d\n", ctx.mype(), j);
        memcpy(kktd.mat->data+(j+lcol)*kktd.getNR()+lrow,
                recvbuffer+buffi, nlrows*sizeof(double));
        buffi += nlrows;
      }
    }

    delete [] start;
    stochNode->resMon.recReduceTmLocal_stop(); 

    
  }
  //printf("done reducing\n");

  delete [] recvcounts;
  delete [] recvcounts_real;
  delete [] recvbuffer;
  delete [] sendbuffer;
  delete [] schurbuffer;
  delete [] nr_counts;
  
  finalizeKKT(prob, vars);
  
  //double val = kktd.getVal(PROW,PCOL);
  //if (cinfo.mype == 0) {
  //  printf("(%d,%d) --- %f\n", PROW, PCOL, val);
  //}

  factorizeKKT();

#ifdef TIMING
  afterFactor();
#endif

}



