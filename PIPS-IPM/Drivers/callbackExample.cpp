#include "StochInputTree.h"
#include "PIPSIpmInterface.h"
#include "sFactoryAug.h"
#include "sFactoryAugSchurLeaf.h"
#include "sFactoryAugMumpsLeaf.h"
//#include "MehrotraStochSolver.h"
#include "GondzioStochSolver.h"


#include "mpi.h"

#define LINKING_CONS 1

extern int gOuterSolve;
extern int gInnerSCsolve;

extern "C" typedef int (*FNNZ)(void* user_data, int id, int* nnz);

/* Row-major format */
extern "C" typedef int (*FMAT)(void* user_data, int id, int* krowM, int* jcolM, double* M);

extern "C" typedef int (*FVEC)(void* user_data, int id, double* vec, int len);


/** Problem parameters and data */
class ProbData
{
public: //data
  int nScenarios;

public: //methods
  ProbData(int nScenarios)
  {
     this->nScenarios = nScenarios;
  };

//  ~ProbData();
};

extern "C" {

// m \times n - submatrix (equality and inequality is the same)
// DONE (ERROR: submatrices do not have the same number of columns
// for equality and inequality equations)
int nSize(void* user_data, int id, int* nnz)
{
   if ( id == 0 )
     *nnz = 9;
   else if ( id == 3 || id == 9 || id == 21)
     *nnz = 6;
   else
    *nnz = 5;
   return 0;
}

//y=equality m x n - submatrix
//DONE (CHECKED)
int mySize(void* user_data, int id, int* nnz)
{
   if ( id == 21 )
     *nnz = 8;
   else if (id == 10)
     *nnz = 4;
   else if ( id == 0 )
     *nnz = 3;
   else if( 10 < id && id <= 23 )
      *nnz = 7;
   else
      *nnz = 0;
   return 0;
}

// z=inequality m x n - submatrix
// DONE
int mzSize(void* user_data, int id, int* nnz)
{
   if( id == 3 || id == 9 )
     *nnz = 9;
   else if (id == 10)
     *nnz = 3;
   else if ( 0 <= id && id <= 9)
     *nnz = 7;
   else
      *nnz = 0;

   return 0;
}

// y=equality, l=linking, m \times n
// DONE (CHECKED)
int mylSize(void* user_data, int id, int* nnz)
{

   *nnz = 0;

   return 0;
}

// z=inequality, l=linking, m \times n
// DONE (CHECKED)
int mzlSize(void* user_data, int id, int* nnz)
{
   *nnz = 0;

   return 0;
}

// DONE (CHECKED)
int nnzMatEqStage1(void* user_data, int id, int* nnz)
{
   if( id == 0 )
     *nnz = 6;
   else if ( id > 0 && id < 11)
     *nnz = 0;
   else if ( id == 21)
     *nnz = 4;
   else
     *nnz = 3;
   return 0;
}

// DONE (CHECKED)
int nnzMatIneqStage1(void* user_data, int id, int* nnz)
{
   if(( id == 3 ) || (id == 9 ))
      *nnz = 4;
   else if ( id == 0)
     *nnz = 9;
   else if (id <= 10)
     *nnz = 3;
   else
     *nnz = 0;
   return 0;
}

// DONE (CHECKED)
int nnzMatEqStage2(void* user_data, int id, int* nnz)
{
   if (id < 10)
      *nnz = 0;
   else if (id == 10)
      *nnz = 9;
   else if( id == 21 )
      *nnz = 13;
   else
     *nnz = 12;
   return 0;
}

// DONE (CHECKED)
int nnzMatIneqStage2(void* user_data, int id, int* nnz)
{
   if ( id == 0)
     *nnz = 0;
   else if( id == 3 || id == 9)
     *nnz = 13;
   else if (id == 10)
     *nnz = 3;
   else if (id < 10)
     *nnz = 12;
   else
     *nnz = 0;
   return 0;
}

// DONE (CHECKED)
int nnzMatEqLink(void* user_data, int id, int* nnz)
{
   *nnz = 0;

   return 0;
}

// DONE (CHECKED)
int nnzMatIneqLink(void* user_data, int id, int* nnz)
{
   *nnz = 0;

   return 0;
}

// DONE (CHECKED)
int nnzAllZero(void* user_data, int id, int* nnz)
{
   *nnz = 0;
   return 0;
}

// DONE (CHECKED)
int vecAllZero(void* user_data, int id, double* vec, int len)
{
   int i;
   for( i = 0; i < len; i++ )
      vec[i] = 0.0;

   return 0;
}

// DONE (CHECKED)
int vecEqRhs(void* user_data, int id, double* vec, int len)
{
   if (id == 10)
   {
    vec[0] = 0;
     vec[1] = 0;
    vec[2] = 0;
    vec[3] = 3.1150345460011648;
//     vec[4] = 2.8682471601853763;
//    vec[5] = 2.7578793581978331;
//    vec[6] = 2.6466753561204408;
   }
   else if (id == 11)
   {
    vec[0] = 2.6469978527163391;
    vec[1] = 2.6217500008263914;
    vec[2] = 2.5362489258990175;
     vec[3] = 3.0285483627333685;
     vec[4] = 3.5490520058344135;
    vec[5] = 3.8335485972911396;
    vec[6] = 3.7931765573895335;
   }
   else if (id == 12)
   {
    vec[0] = 3.9967640325717175;
    vec[1] = 4.1835481024983876;
    vec[2] = 4.1632030743291679;
     vec[3] = 4.0946213641646949;
     vec[4] = 3.9588509809393955;
    vec[5] = 3.9833774731136717;
    vec[6] = 3.8997352078425225;
   }
   else if (id == 13)
   {
    vec[0] = 3.9473236153946694;
    vec[1] = 3.9294681937168408;
    vec[2] = 3.9700395879722423;
    vec[3] = 3.9028760909027946;
     vec[4] = 3.6953135778080046;
     vec[5] = 3.6938516318291503;
    vec[6] = 0.1310966731219193;
   }
   else if (id == 14)
   {
    vec[0] = 0.080921960266879867;
    vec[1] = 0.05469347604221346;
    vec[2] = 0.058115332111554147;
    vec[3] = 0.088905491069479869;
    vec[4] = 0.13565900242179263;
     vec[5] = 0.265660295747978;
    vec[6] = 0.5336551407578024;
   }
   else if (id == 15)
   {
    vec[0] = 0.79479230954117974;
    vec[1] = 0.75259458182259953;
    vec[2] = 0.68759611704822032;
    vec[3] = 0.67618483907675331;
    vec[4] = 0.72067354994448041;
     vec[5] = 0.76970058933623264;
    vec[6] = 0.77083517146723923;
   }
   else if (id == 16)
   {
    vec[0] = 0.84384116782006646;
    vec[1] = 0.96471780254650508;
    vec[2] = 1.1152244860019158;
    vec[3] = 1.1095297564597497;
    vec[4] = 0.91224337898746466;
     vec[5] = 0.64655479034778152;
    vec[6] = 0.41618443431309615;
   }
   else if (id == 17)
   {
    vec[0] = 0.29645110927251767;
    vec[1] = 0.22916820701513163;
     vec[2] = 0.0;
    vec[3] = 0.0;
    vec[4] = 0.0;
    vec[5] = 0.0;
    vec[6] = 0.0;
   }
   else if ( id == 21 )
   {
     vec[0] = 0.0;
    vec[1] = 0.0;
    vec[2] = 0.0;
    vec[3] = 0.0;
    vec[4] = 0.0;
    vec[5] = 0.0;
    vec[6] = 0.0;
    vec[7] = 0.0;
   }
   else if (17 < id && id <= 23)
   {
     vec[0] = 0.0;
    vec[1] = 0.0;
    vec[2] = 0.0;
    vec[3] = 0.0;
    vec[4] = 0.0;
    vec[5] = 0.0;
    vec[6] = 0.0;
   }
   return 0;
}

// DONE (CHECKED)
int vecIneqRhs(void* user_data, int id, double* vec, int len)
{

   int i;

      for( i = 0; i < len; i++ )
         vec[i] = 0.0;

   return 0;
}

// DONE (CHECKED)
int vecIneqRhsLink(void* user_data, int id, double* vec, int len)
{
   //vec[0] = 4.0;

   return 0;
}

// DONE (CHECKED)
int vecIneqRhsActive(void* user_data, int id, double* vec, int len)
{
   int i;
   for( i = 0; i < len; i++ )
      vec[i] = 1.0;

   return 0;
}

// DONE (CHECKED)
int vecIneqRhsActiveLink(void* user_data, int id, double* vec, int len)
{
   //vec[0] = 1.0;

   return 0;
}

// DONE (CHECKED)
int vecObj(void* user_data, int id, double* vec, int len)
{
   vec[0] = 0.0870147443485377;
   vec[1] = 0.0817961426809888;
   vec[2] = 0.134037429705466;

   return 0;
}

// DONE (CHECKED)
int vecXlb(void* user_data, int id, double* vec, int len)
{
   int i;
   for( i = 0; i < len; i++ )
      vec[i] = 0.0;

   return 0;
}

// DONE (CHECKED)
int vecXlbActive(void* user_data, int id, double* vec, int len)
{
   int i;
   for( i = 0; i < len; i++ )
      vec[i] = 1.0;

   return 0;
}

// DONE (CHECKED)
int vecLinkRhs(void* user_data, int id, double* vec, int len)
{
   //
   return 0;
}

// DONE (CHECKED)
int matAllZero(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
    return 0;
}

// DONE (CHECKED)
int matEqStage1(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
   if( id == 0 )
    {
     M[0] = -1.0;
     M[1] = 1.0;
     M[2] = -1.0;
     M[3] = 1.0;
     M[4] = -0.2791456700785919;
     M[5] = 1.0;

      krowM[0] = 0;
     krowM[1] = 2;
     krowM[2] = 4;
     krowM[3] = 6;

     jcolM[0] = 1;
     jcolM[1] = 4;
     jcolM[2] = 2;
     jcolM[3] = 5;
     jcolM[4] = 3;
     jcolM[5] = 8;
     }
    else if( id == 11 )
    {
     M[0] = -1.0;
     M[1] = -1.0;
     M[2] = -0.16350071001847705;

     krowM[0] = 0;
     krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
    }
   else if (id == 12)
   {
     M[0] = -1.0;
     M[1] = -1.0;
     M[2] = -0.19685471646215402;

     krowM[0] = 0;
     krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
   }
   else if (id == 13)
   {
     M[0] = -1.0;
     M[1] = -1.0;
     M[2] = -0.202695406649129;

     krowM[0] = 0;
     krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
   }
   else if (id == 14)
   {
     M[0] = -1.0;
     M[1] = -1.0;
     M[2] = -0.45153311537522034;

     krowM[0] = 0;
     krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
   }
   else if (id == 15)
   {
     M[0] = -1.0;
     M[1] = -1.0;
     M[2] = -0.34119787231863596;

     krowM[0] = 0;
     krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
   }
   else if (id == 16)
   {
     M[0] = -1.0;
     M[1] = -1.0;
     M[2] = -0.1080766877949475;

     krowM[0] = 0;
     krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
   }
   else if (id == 17)
   {
     M[0] = -1.0;
     M[1] = -1.0;
     M[2] = -0.384276132690269;

     krowM[0] = 0;
     krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
   }
   else if (id == 18)
   {
     M[0] = -1.0;
     M[1] = -1.0;
     M[2] = -0.452150528236392;

     krowM[0] = 0;
     krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
   }
   else if (id == 19)
   {
     M[0] = -1.0;
     M[1] = -1.0;
     M[2] = -0.3612012261832522;

     krowM[0] = 0;
     krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
   }
   else if (id == 20)
   {
     M[0] = -1.0;
     M[1] = -1.0;
     M[2] = -0.09114665517232735;

     krowM[0] = 0;
     krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
   }
   else if (id == 21)
   {
     M[0] = 1.0;
     M[1] = -1.0;
     M[2] = -1.0;
     M[3] = -0.16278100269280252;

     krowM[0] = 0;
     krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;
     krowM[4] = 4;

     jcolM[0] = 3;
     jcolM[1] = 1;
     jcolM[2] = 2;
     jcolM[3] = 3;
   }
   else if (id == 22)
   {
     M[0] = -1.0;
     M[1] = -1.0;
     M[2] = -0.3005857732544753;

     krowM[0] = 0;
     krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[0] = 2;
     jcolM[0] = 3;
   }
   else if (id == 23)
   {
     M[0] = -1.0;
     M[1] = -1.0;
     M[2] = -0.11703452488540654;

     krowM[0] = 0;
     krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[0] = 2;
     jcolM[0] = 3;
   }
   return 0;
}



// DONE (CHECKED)
int matIneqLink(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
    // krowM[0] = 0;
    // krowM[1] = 0;

    return 0;
}

// DONE (CHECKED)
int matIneqStage1(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
   if( id == 0 )
    {
     krowM[0] = 0;
      krowM[1] = 0;
     krowM[2] = 0;
     krowM[3] = 0;
     krowM[4] = 1;
     krowM[5] = 2;
     krowM[6] = 6;
     krowM[7] = 9;

     jcolM[0] = 6;
     jcolM[1] = 7;
     jcolM[2] = 4;
     jcolM[3] = 5;
     jcolM[4] = 6;
     jcolM[5] = 8;
     jcolM[6] = 4;
     jcolM[7] = 5;
     jcolM[8] = 7;

     M[0] = 1.0;
     M[1] = 1.0;
     M[2] = -1.0;
     M[3] = 1.0;
     M[4] = -1.0;
     M[5] = 1.0;
     M[6] = 0.7;
     M[7] = -1.6666666666666667;
     M[8] = -1.0;
    }
    else if( id == 1 )
    {
     M[0]=-1.0;
      M[1]=-1.0;
      M[2]=-0.3125616473473225;

     krowM[0] = 0;
      krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
    }
    else if( id == 2 )
    {
     M[0]=-1.0;
      M[1]=-1.0;
      M[2]=-0.42299791324855057;

     krowM[0] = 0;
      krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
    }
   else if( id == 3 )
    {
     M[0]=1.0;
      M[1]=-1.0;
     M[2]=-1.0;
     M[3]=-0.24027346743022576;

     krowM[0] = 0;
      krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;
     krowM[4] = 4;

     jcolM[0] = 1;
     jcolM[1] = 1;
     jcolM[2] = 2;
     jcolM[3] = 3;
    }
    else if( id == 4 )
    {
     M[0]=-1.0;
      M[1]=-1.0;
      M[2]=-0.30056068941739067;

     krowM[0] = 0;
      krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
    }
    else if( id == 5 )
    {
     M[0]=-1.0;
      M[1]=-1.0;
      M[2]=-0.2889374044533182;

     krowM[0] = 0;
      krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
    }
    else if( id == 6 )
    {
     M[0]=-1.0;
      M[1]=-1.0;
      M[2]=-0.4240748903345609;

     krowM[0] = 0;
      krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
    }
    else if( id == 7 )
    {
     M[0]=-1.0;
      M[1]=-1.0;
      M[2]=-0.08264330469163483;

     krowM[0] = 0;
      krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
    }
    else if( id == 8 )
    {
     M[0]=-1.0;
      M[1]=-1.0;
      M[2]=-0.311081709650903;

     krowM[0] = 0;
      krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
    }
    else if( id == 9 )
    {
     M[0]= 1.0;
      M[1]=-1.0;
     M[2]=-1.0;
      M[3]=-0.12498970392247788;

     krowM[0] = 0;
      krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;
     krowM[4] = 4;

     jcolM[0] = 1;
     jcolM[1] = 1;
     jcolM[2] = 2;
     jcolM[3] = 3;
    }
    else if( id == 10 )
    {
     M[0]=-1.0;
      M[1]=-1.0;
      M[2]=-0.18218392485489793;

     krowM[0] = 0;
      krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;

     jcolM[0] = 1;
     jcolM[1] = 2;
     jcolM[2] = 3;
    }
   return 0;
}

// DONE (CHECKED)
int matEqStage2(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
    if( id == 10 )
    {
     M[0] = 1.0;
     M[1] = 1.0;
     M[2] = -1.0;
     M[3] = 1.0;
     M[4] = -1.0;
     M[5] = 1.0;
     M[6] = 0.7;
     M[7] = -1.6666666666666667;
     M[8] = -1.0;

      krowM[0] = 0;
     krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 6;
     krowM[4] = 9;

     jcolM[0] = 2;
     jcolM[1] = 3;
     jcolM[2] = 0;
     jcolM[3] = 1;
     jcolM[4] = 2;
     jcolM[5] = 4;
     jcolM[6] = 0;
     jcolM[7] = 1;
     jcolM[8] = 3;
   }
   else if( id == 21 )
    {
     M[0] = -1.0;
     M[1] = 1.0;
     M[2] = 1.0;
     M[3] = 1.0;
     M[4] = 1.0;
     M[5] = 1.0;
     M[6] = -1.0;
     M[7] = 1.0;
     M[8] = -1.0;
     M[9] = 1.0;
     M[10] = 0.7;
     M[11] = -1.6666666666666667;
     M[12] = -1.0;

     krowM[0] = 0;
      krowM[1] = 1;
     krowM[2] = 2;
      krowM[3] = 3;
     krowM[4] = 4;
      krowM[5] = 5;
     krowM[6] = 6;
      krowM[7] = 10;
     krowM[8] = 13;

     jcolM[0] = 0;
     jcolM[1] = 1;
     jcolM[2] = 2;
     jcolM[3] = 5;
     jcolM[4] = 3;
     jcolM[5] = 4;
     jcolM[6] = 1;
     jcolM[7] = 2;
     jcolM[8] = 3;
     jcolM[9] = 5;
     jcolM[10] = 1;
     jcolM[11] = 2;
     jcolM[12] = 4;
   }
    else if (id > 10 && id <= 23)
    {
      M[0] = 1.0;
      M[1] = 1.0;
     M[2] = 1.0;
      M[3] = 1.0;
     M[4] = 1.0;
      M[5] = -1.0;
     M[6] = 1.0;
      M[7] = -1.0;
     M[8] = 1.0;
      M[9] = 0.7;
     M[10] = -1.6666666666666667;
      M[11] = -1.0;

      krowM[0] = 0;
      krowM[1] = 1;
      krowM[2] = 2;
     krowM[3] = 3;
      krowM[4] = 4;
      krowM[5] = 5;
     krowM[6] = 9;
     krowM[7] = 12;

      jcolM[0] = 0;
      jcolM[1] = 1;
     jcolM[2] = 4;
      jcolM[3] = 2;
     jcolM[4] = 3;
      jcolM[5] = 0;
     jcolM[6] = 1;
      jcolM[7] = 2;
     jcolM[8] = 4;
      jcolM[9] = 0;
     jcolM[10] = 1;
      jcolM[11] = 3;
    }
    return 0;
}

// DONE (CHECKED)
int matIneqStage2(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
   if ( id == 3 || id == 9)
   {
     M[0]=-1.0;
     M[1]=1.0;
     M[2]=1.0;
     M[3]=1.0;
     M[4]=1.0;
     M[5]=1.0;
     M[6]=-1.0;
     M[7]=1.0;
     M[8]=-1.0;
     M[9]=1.0;
     M[10]=0.7;
     M[11]=-1.6666666666666667;
     M[12]=-1.0;

      krowM[0] = 0;
      krowM[1] = 1;
     krowM[2] = 2;
     krowM[3] = 3;
     krowM[4] = 4;
     krowM[5] = 5;
     krowM[6] = 6;
     krowM[7] = 7;
     krowM[8] = 10;
     krowM[9] = 13;

      jcolM[0] = 0;
     jcolM[1] = 1;
     jcolM[2] = 2;
     jcolM[3] = 5;
     jcolM[4] = 3;
     jcolM[5] = 4;
     jcolM[6] = 1;
     jcolM[7] = 2;
     jcolM[8] = 3;
     jcolM[9] = 5;
     jcolM[10] = 1;
     jcolM[11] = 2;
     jcolM[12] = 4;
   }
   else if ( id == 10 )
   {
     M[0] = 1.0;
    M[1] = 1.0;
    M[2] = 1.0;

    krowM[0] = 0;
     krowM[1] = 1;
    krowM[2] = 2;
    krowM[3] = 3;

    jcolM[0] = 0;
    jcolM[1] = 1;
    jcolM[2] = 4;

   }
   else if (id >=1 && id < 10)
   {
     M[0] = 1.0;
    M[1] = 1.0;
    M[2] = 1.0;
    M[3] = 1.0;
    M[4] = 1.0;
     M[5] = -1.0;
    M[6] = 1.0;
    M[7] = -1.0;
    M[8] = 1.0;
    M[9] = 0.7;
     M[10] = -1.66666666666666666666666666666;
    M[11] = -1.0;

    krowM[0] = 0;
    krowM[1] = 1;
    krowM[2] = 2;
    krowM[3] = 3;
    krowM[4] = 4;
    krowM[5] = 5;
    krowM[6] = 9;
    krowM[7] = 12;

    jcolM[0] = 0;
    jcolM[1] = 1;
    jcolM[2] = 4;
    jcolM[3] = 2;
    jcolM[4] = 3;
    jcolM[5] = 0;
    jcolM[6] = 1;
    jcolM[7] = 2;
    jcolM[8] = 4;
    jcolM[9] = 0;
    jcolM[10] = 1;
    jcolM[11] = 3;
   }
   return 0;
}

// DONE (CHECKED)
int matEqLink(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
   return 0;
}

} /* extern C */


int main(int argc, char ** argv) {

  MPI_Init(&argc, &argv);

  int nScenarios = 20;

  // set callbacks
  FNNZ nCall = &nSize;
  FNNZ myCall = &mySize;
  FNNZ mzCall = &mzSize;
  FNNZ mylCall = &mylSize;
  FNNZ mzlCall = &mzlSize;
  FNNZ fnnzQ = &nnzAllZero;
  FNNZ fnnzA = &nnzMatEqStage1;
  FNNZ fnnzB = &nnzMatEqStage2;

  FNNZ fnnzC = &nnzMatIneqStage1;//FNNZ fnnzC = &nnzAllZero;
  FNNZ fnnzD = &nnzMatIneqStage2;//FNNZ fnnzD = &nnzAllZero;

#if LINKING_CONS
  FNNZ fnnzBl = &nnzMatEqLink;
  FVEC fbl = &vecLinkRhs;
  FMAT fBl = &matEqLink;
  FMAT fDl = &matIneqLink;
  FNNZ fnnzDl = &nnzMatIneqLink;

  FVEC fdlupp = &vecIneqRhsLink;
  FVEC fdllow = &vecAllZero;
  FVEC fidlupp = &vecIneqRhsActiveLink;
  FVEC fidllow = &vecAllZero;

#else
  FNNZ fnnzBl = &nnzAllZero;
  FVEC fbl = &vecAllZero;
  FMAT fBl = &matAllZero;

  FMAT fDl = &matAllZero;
  FNNZ fnnzDl = &nnzAllZero;

  FVEC fdlupp = &vecAllZero;
  FVEC fdllow = &vecAllZero;
  FVEC fidlupp = &vecAllZero;
  FVEC fidllow = &vecAllZero;
#endif

  FVEC fc = &vecObj;
  FVEC fb = &vecEqRhs;

  FVEC fclow = &vecAllZero;
  FVEC fcupp = &vecIneqRhs;//FVEC fcupp = &vecAllZero;
  FVEC fxlow = &vecXlb;
  FVEC fxupp = &vecAllZero;
  FVEC ficlow = &vecAllZero;
  FVEC fixlow = &vecXlbActive;
  FVEC ficupp = &vecIneqRhsActive;//FVEC ficupp = vecAllZero;
  FVEC fixupp = &vecAllZero;


  FMAT fQ = &matAllZero;
  FMAT fA = &matEqStage1;
  FMAT fB = &matEqStage2;
  FMAT fC = &matIneqStage1;//FMAT fC = &matAllZero;
  FMAT fD = &matIneqStage2;//FMAT fD = &matAllZero;

  ProbData probData(nScenarios);

  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);


#if LINKING_CONS
  //build the problem tree
  StochInputTree::StochInputNode dataLinkCons(&probData, 0,
                  nCall, myCall, mylCall, mzCall, mzlCall,
                  fQ, fnnzQ, fc,
                  fA, fnnzA,
                  fB, fnnzB,
                  fBl, fnnzBl,
                  fb, fbl,
                  fC, fnnzC,
                  fD, fnnzD,
                  fDl, fnnzDl,
                  fclow, ficlow, fcupp, ficupp,
                  fdllow, fidllow, fdlupp, fidlupp,
                  fxlow, fixlow, fxupp, fixupp, false );

  StochInputTree* root = new StochInputTree(dataLinkCons);
#else

//  int nx0 = 2;
//  int my0 = 2;
//  int mz0 = 1;

  //build the problem tree
  StochInputTree::StochInputNode data(&probData, 0,
                  nCall, myCall, mzCall, //myl0, mzl0
                  fQ, fnnzQ, fc,
                  fA, fnnzA,
                  fB, fnnzB,
                 //fBL, fnnzBL,
                  fb, //fbl
                  fC, fnnzC,
                  fD, fnnzD,
                 //fDL, fnnzDl,
                  fclow, ficlow, fcupp, ficupp,
                 //fdllow, fidllow, fdlupp, fidlupp,
                  fxlow, fixlow, fxupp, fixupp, false );

  StochInputTree* root = new StochInputTree(data);
#endif


  for( int id = 1; id <= nScenarios; id++ ) {
#if LINKING_CONS
     StochInputTree::StochInputNode dataLinkConsChild(&probData, id,
               nCall, myCall, mylCall, mzCall, mzlCall,
               fQ, fnnzQ, fc,
               fA, fnnzA,
               fB, fnnzB,
               fBl, fnnzBl,
               fb, fbl,
               fC, fnnzC,
               fD, fnnzD,
               fDl, fnnzDl,
               fclow, ficlow, fcupp, ficupp,
               fdllow, fidllow, fdlupp, fidlupp,
               fxlow, fixlow, fxupp, fixupp, false);

     root->AddChild(new StochInputTree(dataLinkConsChild));
#else
//     int nx = 2;
//     int my = 2;
//     int mz = 1;

     StochInputTree::StochInputNode data(&probData, id,
               nCall, myCall, mzCall, //myl, mzl
               fQ, fnnzQ, fc,
               fA, fnnzA,
               fB, fnnzB,
               //fBL, fnnzBL,
               fb, //fbl
               fC, fnnzC,
               fD, fnnzD,
               //fDL, fnnzDL,
               fclow, ficlow, fcupp, ficupp,
               //fdllow, fidllow, fdlupp, fidlupp,
               fxlow, fixlow, fxupp, fixupp, false);

     root->AddChild(new StochInputTree(data));
#endif

  }

  if( rank == 0 )
     cout << "Using a total of " << size << " MPI processes." << endl;

  /* use BiCGStab for outer solve */
  gOuterSolve = 2;
  gInnerSCsolve = 0;
#if defined(WITH_MUMPS_LEAF)
  PIPSIpmInterface<sFactoryAugMumpsLeaf, GondzioStochSolver> pipsIpm(root, MPI_COMM_WORLD, SCALER_EQUI_STOCH, PRESOLVER_NONE);
#elif defined(WITH_PARDISO) && !defined(PARDISO_BLOCKSC)
  PIPSIpmInterface<sFactoryAugSchurLeaf, GondzioStochSolver> pipsIpm(root, MPI_COMM_WORLD, SCALER_EQUI_STOCH, PRESOLVER_NONE);
#else
  PIPSIpmInterface<sFactoryAug, GondzioStochSolver> pipsIpm(root, MPI_COMM_WORLD, SCALER_EQUI_STOCH, PRESOLVER_NONE);
#endif

  if( rank == 0 )
     cout << "PIPSIpmInterface created" << endl;

  if( rank == 0 )
     cout << "solving..." << endl;

  pipsIpm.go();

  const double objective = pipsIpm.getObjective();
  if( rank == 0 )
     cout << "solving finished ... objective value: " << objective << endl;

  delete root;

  MPI_Finalize();

  return 0;
}
