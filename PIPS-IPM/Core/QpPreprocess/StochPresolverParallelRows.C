/*
 * StochPresolverParallelRows.C
 *
 *  Created on: 02.05.2018
 *      Author: bzfuslus
 */

#include "StochPresolverParallelRows.h"

namespace rowlib
{
   bool operator==(rowWithColInd const& row1, rowWithColInd const& row2)
   {
      return row1.id == row2.id;
   }

   std::size_t hash_value(rowWithColInd const& row)
   {
      std::size_t seed = 0;
      boost::hash_combine(seed, (row.lengthA+row.lengthB) );
      for( int i = 0; i < row.lengthA; i++ )
         boost::hash_combine(seed, row.colIndicesA[i]);
      for( int i = 0; i < row.lengthB; i++ )
         boost::hash_combine(seed, (row.colIndicesB[i] + row.offset_nA) );
      return seed;
   }

   bool operator==(rowWithEntries const& row1, rowWithEntries const& row2)
   {
      return row1.id == row2.id;
   }

   std::size_t hash_value(rowWithEntries const& row)
   {
      std::size_t seed = 0;
      for( int i = 0; i < row.lengthA; i++ )
      {
         // Instead of hashing the normalized double coefficient, use an integer representation:
         int value_to_hash;
         double mantisse = std::frexp( (row.norm_entriesA[i] + offset_hash_double), &value_to_hash);
         value_to_hash += 10*( (int)trunc(mantisse*10000) );
         boost::hash_combine(seed, value_to_hash);
      }
      for( int i = 0; i < row.lengthB; i++ )
      {
         int value_to_hash;
         double mantisse = std::frexp( (row.norm_entriesB[i] + offset_hash_double), &value_to_hash);
         value_to_hash += 10*( (int)trunc(mantisse*10000) );
         boost::hash_combine(seed, value_to_hash);
      }
      return seed;
   }
}

StochPresolverParallelRows::StochPresolverParallelRows(PresolveData& presData)
: StochPresolverBase(presData)
{
   currCmat = NULL;
   currCmatTrans = NULL;
   currDmat = NULL;
   currDmatTrans = NULL;
   currNnzRowC = NULL;
   norm_Amat = NULL;
   norm_Bmat= NULL;
   norm_Cmat = NULL;
   norm_Dmat = NULL;
   norm_b = NULL;
   norm_cupp = NULL;
   norm_clow = NULL;
   norm_icupp = NULL;
   norm_iclow = NULL;

   mA = 0;
   nA = 0;
}

StochPresolverParallelRows::~StochPresolverParallelRows()
{
 // todo
   rowsFirstHashTable = boost::unordered_set<rowlib::rowWithColInd, boost::hash<rowlib::rowWithColInd> >();
   rowsSecondHashTable  = boost::unordered_set<rowlib::rowWithEntries, boost::hash<rowlib::rowWithEntries> >();
}

bool StochPresolverParallelRows::applyPresolving(int& nelims)
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   if( myRank == 0 )
      cout<<"Before Parallel Row Presolving:"<<endl;
   countRowsCols();

   if( myRank == 0 ) cout<<"Start Parallel Row Presolving..."<<endl;

   StochGenMatrix& matrixA = dynamic_cast<StochGenMatrix&>(*(presProb->A));
   StochGenMatrix& matrixC = dynamic_cast<StochGenMatrix&>(*(presProb->C));
   presData.resetRedCounters();
   int nRowElims = 0;

   // for children:
   for( size_t child_it = 0; child_it< matrixA.children.size(); child_it++)
   {
      // copy and normalize A,B,C,D and b,c,d:
      if( setNormalizedPointers((int)child_it, matrixA, matrixC ) )
      {
         cout<<"Normalized pointers are set for child #"<<(int)child_it<<endl;
         // Prepare unordered set 'rowsFirstHashTable':
         rowsFirstHashTable.clear();

         // Per row, add row to the set 'rowsFirstHashTable':
         if( norm_Amat )
            insertRowsIntoHashtable( rowsFirstHashTable, norm_Amat, norm_Bmat, EQUALITY_SYSTEM );
         assert( (int)rowsFirstHashTable.size() <= mA );
         if( norm_Cmat )
         {
            insertRowsIntoHashtable( rowsFirstHashTable, norm_Cmat, norm_Dmat, INEQUALITY_SYSTEM );
            assert( (int)rowsFirstHashTable.size() <= mA + norm_Cmat->m);
         }

         // Prints for visualization:
         /*std::cout << "unordered_set rowsFirstHashTable has size " << rowsFirstHashTable.size() << '\n';
         for (size_t i=0; i<rowsFirstHashTable.bucket_count(); ++i)
         {
             std::cout << "bucket #" << i << " contains: ";
             for (boost::unordered_set<rowlib::rowWithColInd>::local_iterator it = rowsFirstHashTable.begin(i); it!=rowsFirstHashTable.end(i); ++it)
                 std::cout << " row id#"<<it->id;
             std::cout << "\n";
         }*/

         // Second Hashing: Per bucket, do Second Hashing:
         rowsSecondHashTable.clear();
         for( size_t i = 0; i < rowsFirstHashTable.bucket_count(); ++i )
         {
            for( boost::unordered_set<rowlib::rowWithColInd>::local_iterator it =
                  rowsFirstHashTable.begin(i); it != rowsFirstHashTable.end(i); ++it )
            {
               rowsSecondHashTable.emplace(it->id, it->offset_nA, it->lengthA, it->colIndicesA, it->norm_entriesA,
                     it->lengthB, it->colIndicesB, it->norm_entriesB);
            }

            // Prints for visualization:
            /*if( !rowsSecondHashTable.empty() )
            {
               std::cout << "unordered_set rowsSecondHashTable has size " << rowsSecondHashTable.size() << '\n';
               for (size_t i=0; i<rowsSecondHashTable.bucket_count(); ++i)
               {
                  if( rowsSecondHashTable.bucket_size(i) > 0 )
                  {
                     std::cout << "bucket #" << i << " contains: "<<rowsSecondHashTable.bucket_size(i)<<" entries.";
                     for (boost::unordered_set<rowlib::rowWithEntries>::local_iterator it = rowsSecondHashTable.begin(i); it!=rowsSecondHashTable.end(i); ++it)
                     {
                        std::cout << " row id#"<<it->id<<", entries: ";
                        for(int k=0; k<it->lengthA; k++)
                           std::cout <<it->norm_entriesA[k]<<", ";
                        for(int k=0; k<it->lengthB; k++)
                           std::cout <<it->norm_entriesB[k]<<", ";
                     }
                     std::cout << "\n";
                  }
               }
            }*/
            // Compare the rows in the final (from second hash) bin:
            compareRowsInSecondHashTable(nRowElims);

            rowsSecondHashTable.clear();
         }

         // Objects created with new have to be deleted at the end of each child (norm_Amat etc)
         deleteNormalizedPointers((int)child_it, matrixA, matrixC);

         /*std::cout << "unordered_set rows has size " << rows.size() <<" after deleting"<< '\n';
         for (size_t i=0; i<rows.bucket_count(); ++i)
         {
             std::cout << "bucket #" << i << " contains: ";
             for (boost::unordered_set<rowlib::row>::local_iterator it = rows.begin(i); it!=rows.end(i); ++it)
                 std::cout << " colIndices:"<<it->colIndices[0];
             std::cout << "\n";
         }*/
      }
   }

   updateNnzColParent(MPI_COMM_WORLD);

   //SparseStorageDynamic norm_storage = matrixA.Bmat->getStorageDynamicRef();

   rowsFirstHashTable.clear();
   rowsSecondHashTable.clear();

   MPI_Barrier( MPI_COMM_WORLD);
   synchronize(nRowElims);
   if( myRank == 0 )
      cout<<"Removed "<<nRowElims<<" Rows in Parallel Row Presolving."<<endl;
   //assert(0);

   if( myRank == 0 )
      cout<<"After Parallel Row Presolving:"<<endl;
   countRowsCols();

   //countDuplicateRows(matrixC, INEQUALITY_SYSTEM);
   //countDuplicateRows(matrixA, EQUALITY_SYSTEM);

   return true;
}

/** If it is no dummy child, sets normalized pointers:
 * For the matrix blocks A, B, C and D and their transposed matrices and for the lhs/rhs.
 * Sets the pointers to currNnzRow, currNnzRowC, currNnzColChild, currRedColParent.
 * Sets mA and nA correctly.
 */
bool StochPresolverParallelRows::setNormalizedPointers(int it, StochGenMatrix& matrixA, StochGenMatrix& matrixC)
{
   // check if it is no dummy child and copy the matrices:

   if( !childIsDummy(matrixA, it, EQUALITY_SYSTEM) )
   {
      norm_Amat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.children[it]->Amat)->getStorageDynamicRef());
      norm_Bmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.children[it]->Bmat)->getStorageDynamicRef());
      norm_b = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).children[it]->vec)->cloneFull();
      setCPAmatsChild( presProb->A, it, EQUALITY_SYSTEM);
      setCPBmatsChild( presProb->A, it, EQUALITY_SYSTEM);
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->children[it]->vec);
   }
   else
   {
      norm_Amat = NULL;
      norm_Bmat = NULL;
      norm_b = NULL;

      currAmat = NULL;
      currAmatTrans = NULL;
      currBmat = NULL;
      currBmatTrans = NULL;
      currNnzRow = NULL;
   }
   if( !childIsDummy(matrixC, it, INEQUALITY_SYSTEM) )
   {
      norm_Cmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.children[it]->Amat)->getStorageDynamicRef());
      norm_Dmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.children[it]->Bmat)->getStorageDynamicRef());
      norm_cupp = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).children[it]->vec)->cloneFull();
      norm_clow = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).children[it]->vec)->cloneFull();
      norm_icupp = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).children[it]->vec)->cloneFull();
      norm_iclow = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).children[it]->vec)->cloneFull();

      currCmat = dynamic_cast<SparseGenMatrix*>(matrixC.children[it]->Amat)->getStorageDynamic();
      currCmatTrans = dynamic_cast<SparseGenMatrix*>(matrixC.children[it]->Amat)->getStorageDynamicTransposed();
      currDmat = dynamic_cast<SparseGenMatrix*>(matrixC.children[it]->Bmat)->getStorageDynamic();
      currDmatTrans = dynamic_cast<SparseGenMatrix*>(matrixC.children[it]->Bmat)->getStorageDynamicTransposed();
      currNnzRowC = dynamic_cast<SimpleVector*>(presData.nRowElemsC->children[it]->vec);
   }
   else
   {
      norm_Cmat = NULL;
      norm_Dmat = NULL;
      norm_cupp = NULL;
      norm_clow = NULL;
      norm_icupp = NULL;
      norm_iclow = NULL;

      currCmat = NULL;
      currCmatTrans = NULL;
      currDmat = NULL;
      currDmatTrans = NULL;
      currNnzRowC = NULL;
   }
   if( !norm_Amat && !norm_Cmat )   // case no child exists
   {
      currRedColParent = NULL;
      currNnzColChild = NULL;
      return false;
   }
   else  // set mA, nA correctly
   {
      if( norm_Amat )
      {
         mA = norm_Amat->m;
         nA = norm_Amat->n;
      }
      else
      {
         mA = 0;
         nA = 0;
      }
      currRedColParent = dynamic_cast<SimpleVector*>(presData.redCol->vec);
      currNnzColChild = dynamic_cast<SimpleVector*>(presData.nColElems->children[it]->vec);
   }
   if(norm_Amat && norm_Cmat) // only for assertions
   {
      assert( norm_Amat->n == norm_Cmat->n );
      assert( norm_Bmat->n == norm_Dmat->n );
      assert( norm_Amat->m == norm_Bmat->m );
      assert( norm_Cmat->m == norm_Dmat->m );
   }

   // normalize the rows:
   if( norm_Amat )
      normalizeBLocksRowwise( EQUALITY_SYSTEM, norm_Amat, norm_Bmat, norm_b, NULL, NULL, NULL);
   if( norm_Cmat )
      normalizeBLocksRowwise( INEQUALITY_SYSTEM, norm_Cmat, norm_Dmat, norm_cupp, norm_clow, norm_icupp, norm_iclow);

   return true;
}

void StochPresolverParallelRows::deleteNormalizedPointers(int it, StochGenMatrix& matrixA, StochGenMatrix& matrixC)
{
   if( !childIsDummy(matrixA, it, EQUALITY_SYSTEM) )
   {
      assert( norm_Amat && norm_Bmat && norm_b );
      delete norm_Amat;
      delete norm_Bmat;
      delete norm_b;
   }
   if( !childIsDummy(matrixC, it, INEQUALITY_SYSTEM) )
   {
      assert( norm_Cmat && norm_Dmat && norm_cupp && norm_clow && norm_icupp && norm_iclow );
      delete norm_Cmat;
      delete norm_Dmat;
      delete norm_cupp;
      delete norm_clow;
      delete norm_icupp;
      delete norm_iclow;
   }
}

void StochPresolverParallelRows::normalizeBLocksRowwise( SystemType system_type, SparseStorageDynamic* Ablock, SparseStorageDynamic* Bblock,
      SimpleVector* Rhs, SimpleVector* Lhs, SimpleVector* iRhs, SimpleVector* iLhs)
{
   assert( Ablock && Bblock && Rhs );
   int nRows = Ablock->m;
   assert( nRows == Bblock->m);
   if( nRows != Rhs->n )
      cout<<"nRows="<<nRows<<", Rhs->n="<<Rhs->n<<endl;
   assert( nRows == Rhs->n );

   if( system_type == INEQUALITY_SYSTEM )
   {
      assert( Lhs && iRhs && iLhs );
      assert( Lhs->n == nRows && iRhs->n == nRows && iLhs->n == nRows );
   }

   for( int i=0; i<nRows; i++)   // row i
   {
      double maxValue = 0.0;
      bool negateRow = false;
      const int rowStartA = Ablock->rowptr[i].start;
      const int rowEndA = Ablock->rowptr[i].end;
      if( rowStartA < rowEndA )
      {
         if( Ablock->M[rowStartA] < 0)
            negateRow = true;
         for(int k=rowStartA; k<rowEndA; k++)
         {
            if( fabs(Ablock->M[k]) > maxValue )
               maxValue = fabs(Ablock->M[k]);
         }
      }
      const int rowStartB = Bblock->rowptr[i].start;
      const int rowEndB = Bblock->rowptr[i].end;
      if( rowStartB < rowEndB )
      {
         if( Bblock->M[rowStartB] < 0)
            negateRow = true;
         for(int k=rowStartB; k<rowEndB; k++)
         {
            if( fabs(Bblock->M[k]) > maxValue )
               maxValue = fabs(Bblock->M[k]);
         }
      }
      // normalize the row by dividing all entries by maxValue and possible by -1, if negateRow.
      if(negateRow)
         maxValue *= -1.0;
      if( rowStartA < rowEndA )
      {
         for(int k=rowStartA; k<rowEndA; k++)
            Ablock->M[k] /= maxValue;
      }
      if( rowStartB < rowEndB )
      {
         for(int k=rowStartB; k<rowEndB; k++)
            Bblock->M[k] /= maxValue;
      }
      if( system_type == EQUALITY_SYSTEM )
         Rhs->elements()[i] /= maxValue;
      else
      {
         if( iLhs->elements()[i] != 0.0 )
            Lhs->elements()[i] /= maxValue;
         if( iRhs->elements()[i] != 0.0 )
            Rhs->elements()[i] /= maxValue;
         // if multiplied by a negative value, lhs and rhs have to be swapped:
         if(negateRow)
         {
            std::swap(Lhs->elements()[i], Rhs->elements()[i]);
            std::swap(iLhs->elements()[i], iRhs->elements()[i]);
         }
      }
   }
}

/** Inserts all non-empty rows into the given unordered set 'rows'.
 * Ablock and Bblock are supposed to be the two matrix blocks of a child.
 */
void StochPresolverParallelRows::insertRowsIntoHashtable( boost::unordered_set<rowlib::rowWithColInd, boost::hash<rowlib::rowWithColInd> > &rows,
      SparseStorageDynamic* Ablock, SparseStorageDynamic* Bblock, SystemType system_type )
{
   if( Ablock )
   {
      assert( Bblock );
      assert( Ablock->m == Bblock->m );
      if( system_type == EQUALITY_SYSTEM )
         assert( mA == Ablock->m );
      for(int i=0; i<Ablock->m; i++)
      {
         // calculate rowId:
         int rowId = i;
         if( system_type == INEQUALITY_SYSTEM )
            rowId += mA;

         // calculate rowlength of A- and B-block
         const int rowStartA = Ablock->rowptr[i].start;
         const int rowlengthA = Ablock->rowptr[i].end - rowStartA;
         const int rowStartB = Bblock->rowptr[i].start;
         const int rowlengthB = Bblock->rowptr[i].end - rowStartB;

         // colIndices and normalized entries are set as pointers to the original data.
         // create and insert the new element:
         if( rowlengthA != 0 || rowlengthB != 0 )
            rows.emplace(rowId, nA, rowlengthA, &(Ablock->jcolM[rowStartA]), &(Ablock->M[rowStartA]),
                           rowlengthB, &(Bblock->jcolM[rowStartB]), &(Bblock->M[rowStartB]));

         else  // only for assertions
         {
            if( system_type == EQUALITY_SYSTEM )
               assert( currNnzRow->elements()[rowId] == 0.0 );
            else
               assert( currNnzRowC->elements()[rowId] == 0.0 );
         }
      }
   }
}

/*
 * Per bucket in rowsSecondHashTable, compare the containing rows and check if they are parallel.
 * If so, consider the different possible cases.
 */
bool StochPresolverParallelRows::compareRowsInSecondHashTable(int& nRowElims)
{
   // todo The action has to be applied to the original matrices (not the normalized copies)
   // todo make sure that more than 2 parallel rows are treated correctly

   if( rowsSecondHashTable.empty() )
      return true;

   for (size_t i=0; i<rowsSecondHashTable.bucket_count(); ++i)
   {
      for (boost::unordered_set<rowlib::rowWithEntries>::local_iterator it1 = rowsSecondHashTable.begin(i);
            it1!=rowsSecondHashTable.end(i); ++it1)
      {
         // either pairwise comparison OR lexicographical sorting and then compare only neighbors.
         // Here: parwise comparison:
         for (boost::unordered_set<rowlib::rowWithEntries>::local_iterator it2 = it1;
                                    it2!=rowsSecondHashTable.end(i); ++it2)
         {
            // When two parallel rows are found, check if they are both =, both <=, or = and and <=
            if( checkRowsAreParallel( *it1, *it2) )
            {
               cout<<"Found two rows with parallel coefficients: Row "<<it1->id<<" and "<<it2->id<<endl;
               if( it1->id < mA && it2->id < mA )
               {
                  // Case both constraints are equalities
                  if( norm_b->elements()[it1->id] != norm_b->elements()[it2->id])
                     return false;
                  // delete row2 in the original system:
                  eliminateOriginalRow((int)it2->id, nRowElims);

               }
               else if( it1->id >= mA && it2->id >= mA )
               {
                  // Case both constraints are inequalities
                  const int rowId1 = it1->id - mA;
                  const int rowId2 = it2->id - mA;
                  // Set newxlow/newxupp to the bounds of the second inequality:
                  double newxlow = -std::numeric_limits<double>::max();
                  double newxupp = std::numeric_limits<double>::max();
                  if( norm_iclow->elements()[rowId2] != 0.0 )
                     newxlow = norm_clow->elements()[rowId2];
                  if( norm_icupp->elements()[rowId2] != 0.0 )
                     newxupp = norm_cupp->elements()[rowId2];
                  // test if new bounds are tightening:
                  if( newBoundsTightenOldBounds(newxlow, newxupp, rowId1, norm_iclow->elements(),
                        norm_icupp->elements(), norm_clow->elements(), norm_cupp->elements()) )
                  {
                     // tighten the constraint bounds:
                     setNewBounds(rowId1, newxlow, newxupp, norm_iclow->elements(), norm_clow->elements(),
                           norm_icupp->elements(), norm_cupp->elements());
                     // todo: apply the lhs/rhs tightening to the original data (it1)
                  }
                  // delete row2 in the original system:
                  eliminateOriginalRow(rowId2, nRowElims);

               }
               else
               {  // Case one constraint is an equality, one an inequality
                  int id1 = it1->id;
                  int id2 = it2->id;
                  bool swappedId1Id2 = false;
                  if( id1 >= mA )   // swap ids so that id2 is the inequality constraint.
                  {
                     std::swap(id1, id2);
                     swappedId1Id2 = true;
                  }

                  const int ineqRowId = id2 - mA;
                  if( norm_iclow->elements()[ineqRowId] != 0.0 )
                  {
                     if( norm_clow->elements()[ineqRowId] > norm_b->elements()[id1] )
                        return false;
                  }
                  if( norm_icupp->elements()[ineqRowId] != 0.0 )
                  {
                     if( norm_cupp->elements()[ineqRowId] < norm_b->elements()[id1] )
                        return false;
                  }

                  // delete the inequality constraint id2 (which could be either it1 or it2 now):
                  if(swappedId1Id2)
                  {
                     // delete row2 in the original system:
                     eliminateOriginalRow(it1->id, nRowElims);
                  }
                  else
                  {
                     // delete row2 in the original system:
                     eliminateOriginalRow(it2->id, nRowElims);
                  }
               }
            }
         }
      }
   }
   return true;
}

/**
 * Compare two rowWithEntries rows if the normalized coefficients are the same, at the same
 * columns indices. If yes, return true. Else, return false.
 */
bool StochPresolverParallelRows::checkRowsAreParallel( rowlib::rowWithEntries row1, rowlib::rowWithEntries row2)
{
   assert( row1.id >= 0 && row2.id >= 0 );
   if( row1.id == row2.id )
      return false;

   if( row1.lengthA != row2.lengthA || row1.lengthB != row2.lengthB )
      return false;
   for( int i=0; i<row1.lengthA; i++)
   {
      if( row1.colIndicesA[i] != row2.colIndicesA[i] )
         return false;
      if( !PIPSisEQ_withTolerance(row1.norm_entriesA[i], row2.norm_entriesA[i], tol_compare_double) )
         return false;
   }
   for( int i=0; i<row1.lengthB; i++)
   {
      if( row1.colIndicesB[i] != row2.colIndicesB[i] )
         return false;
      if( !PIPSisEQ_withTolerance(row1.norm_entriesB[i], row2.norm_entriesB[i], tol_compare_double) )
         return false;
   }
   return true;
}

/** Eliminate the row in the original system. */
void StochPresolverParallelRows::eliminateOriginalRow(int rowId, int& nRowElims)
{
   assert(rowId>=0);
   if(rowId < mA) // equality row
   {
      assert( norm_Amat );
      if(currNnzRow->elements()[rowId] != 0.0) // check if row was already removed
      {
         cout<<"Delete row in A with id# "<<rowId<<endl;
         nRowElims++;
         removeRow(rowId, currAmat, currAmatTrans, currBmat, currBmatTrans,
               currNnzRow, currRedColParent, currNnzColChild);
      }
   }
   else  // inequality row
   {
      assert( norm_Cmat );
      assert(rowId < mA + norm_Cmat->m );
      const int rowIdC = rowId - mA;
      if(currNnzRowC->elements()[rowIdC] != 0.0) // check if row was already removed
      {
         cout<<"Delete row in C with id# "<<rowIdC<<endl;
         nRowElims++;
         removeRow(rowIdC, currCmat, currCmatTrans, currDmat, currDmatTrans,
               currNnzRowC, currRedColParent, currNnzColChild);
      }
   }
}

/**
 * Remove row rowIdx in Ablock and Bblock. Removes the corresponding column in
 * AblockTrans and BblockTrans. Additionally, sets nnzRow[rowIdx] to 0.0.
 * Increments redColParent by one at each column index the row had an entry.
 * Decrements nnzColChild by one at each column index the row had an entry.
 */
void StochPresolverParallelRows::removeRow(int rowIdx, SparseStorageDynamic* Ablock, SparseStorageDynamic* AblockTrans,
      SparseStorageDynamic* Bblock, SparseStorageDynamic* BblockTrans, SimpleVector* nnzRow,
      SimpleVector* redColParent, SimpleVector* nnzColChild)
{
   assert( Ablock && Bblock);
   assert( AblockTrans && BblockTrans);
   assert( nnzRow && nnzColChild && redColParent );
   assert( rowIdx>=0 && rowIdx<Ablock->m );
   assert( Ablock->m == Bblock->m );
   assert( Ablock->m == nnzRow->n );
   assert( Ablock->n == redColParent->n );
   assert( Bblock->n == nnzColChild->n );

   const int rowStartA = Ablock->rowptr[rowIdx].start;
   const int rowEndA = Ablock->rowptr[rowIdx].end;
   // delete row in AblockTrans:
   for(int k=rowStartA; k<rowEndA; k++)
   {
      const int colIdx = Ablock->jcolM[k];
      double tmp = 0.0;
      removeEntryInDynamicStorage(*AblockTrans, colIdx, rowIdx, tmp);
      // increment redColParent[colIdx]:
      redColParent->elements()[colIdx]++;
   }
   // delete row in Ablock:
   clearRow(*Ablock, rowIdx);

   const int rowStartB = Bblock->rowptr[rowIdx].start;
   const int rowEndB = Bblock->rowptr[rowIdx].end;
   // delete row in BblockTrans:
   for(int k=rowStartB; k<rowEndB; k++)
   {
      const int colIdx = Bblock->jcolM[k];
      double tmp = 0.0;
      removeEntryInDynamicStorage(*BblockTrans, colIdx, rowIdx, tmp);
      // decrement nnzColChild[colIdx]:
      nnzColChild->elements()[colIdx]--;
   }
   // delete row in Bblock:
   clearRow(*Bblock, rowIdx);

   // set nnzRow[rowIdx] to 0.0:
   nnzRow->elements()[rowIdx] = 0.0;

}

void StochPresolverParallelRows::countDuplicateRows(StochGenMatrix& matrix, SystemType system_type)
{
   //TODO: several duplicate rows are counted too often because each duplicate pair is counted as one
   // ie, if rows i,j,k,l are duplicate, they are counted (i,j) (i,k) (i,l) (j,k) (j,l)(k,l) as 6 instead of 4
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   int duplicRow = 0;

   // the B_0 block:
   if( myRank == 0)
   {
      currAmat = dynamic_cast<SparseGenMatrix*>(matrix.Bmat)->getStorageDynamic();
      if(system_type==EQUALITY_SYSTEM)
         currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vec);
      else
         currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vec);
      for( int i=0; i<currNnzRow->n; i++)
      {
         for( int j=i+1; j<currNnzRow->n; j++)
         {
            double* nRow = currNnzRow->elements();
            if( nRow[i] != 0.0 && nRow[i] == nRow[j] &&
                  compareCoefficients(*currAmat, i, j))
            {
               duplicRow++;
               //cout<<"Row "<<i<<" and "<<j<<" of B_0 are duplicate rows."<<endl;
            }
         }
      }
      if(system_type==EQUALITY_SYSTEM)
         cout<<"There are "<<duplicRow<<" duplicate rows in the A_0 block of A."<<endl;
      else
         cout<<"There are "<<duplicRow<<" duplicate rows in the A_0 block of C."<<endl;
   }

   // the children:
   for( size_t it = 0; it< matrix.children.size(); it++)
   {
      if( (system_type==EQUALITY_SYSTEM && !childIsDummy( matrix, (int)it, EQUALITY_SYSTEM)) ||
            (system_type == INEQUALITY_SYSTEM && !childIsDummy( matrix, (int)it, INEQUALITY_SYSTEM)) )
      {
         currAmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamic();
         currBmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamic();
         if(system_type==EQUALITY_SYSTEM)
            currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->children[it]->vec);
         else
            currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->children[it]->vec);
         for( int i=0; i<currNnzRow->n; i++)
         {
            for( int j=i+1; j<currNnzRow->n; j++)
            {
               double* nRow = currNnzRow->elements();
               if( nRow[i] != 0.0 && nRow[i] == nRow[j] &&
                     compareCoefficients(*currAmat, i, j) && compareCoefficients(*currBmat, i, j) )
               {
                  duplicRow++;
                  //cout<<"Row "<<i<<" and "<<j<<" of child "<<it<<" are duplicate rows."<<endl;
               }
            }
         }
      }
   }

   // the linking constraints:
   //TODO: more efficient communication, maybe vector containing info about all rows?
   if( (system_type==EQUALITY_SYSTEM && hasLinking(EQUALITY_SYSTEM)) ||
         (system_type==INEQUALITY_SYSTEM && hasLinking(INEQUALITY_SYSTEM)))
   {
      int nDuplicLinkRow = 0;
      if(system_type==EQUALITY_SYSTEM)
         currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
      else
         currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl);
      for(int i=0; i<currNnzRow->n; i++)
      {
         for( int j=i+1; j<currNnzRow->n; j++)
         {
            int notParallel = 0;
            double* nRow = currNnzRow->elements();
            currBlmat = dynamic_cast<SparseGenMatrix*>(matrix.Blmat)->getStorageDynamic();
            if( nRow[i] != 0.0 && nRow[i] == nRow[j] &&
                  compareCoefficients(*currBlmat, i, j) )
            {
               // check other F_i blocks as well:
               for( size_t it = 0; it< matrix.children.size(); it++)
               {
                  if( (system_type==EQUALITY_SYSTEM && !childIsDummy( matrix, (int)it, EQUALITY_SYSTEM)) ||
                        (system_type==INEQUALITY_SYSTEM && !childIsDummy( matrix, (int)it, INEQUALITY_SYSTEM)) )
                  {
                     currBlmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Blmat)->getStorageDynamic();
                     if( !compareCoefficients(*currBlmat, i, j))
                        notParallel++;
                  }
                  if(system_type==EQUALITY_SYSTEM)
                     currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
                  else
                     currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl);
               }
               if( iAmDistrib )
                  MPI_Allreduce(MPI_IN_PLACE, &notParallel, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
               if( myRank == 0 && notParallel == 0 )
               {
                  duplicRow++;
                  nDuplicLinkRow++;
                  //cout<<"Row "<<i<<" and "<<j<<" of Linking Constraints are duplicate rows."<<endl;
               }
            }
         }
      }
      if(myRank == 0)
      {
         assert(nDuplicLinkRow <= currNnzRow->n);
         //duplicRow += nDuplicLinkRow;
         if(system_type==EQUALITY_SYSTEM)
            cout<<"There are "<<nDuplicLinkRow<<" duplicate rows in the linking rows of A."<<endl;
         else
            cout<<"There are "<<nDuplicLinkRow<<" duplicate rows in the linking rows of C."<<endl;
      }
   }

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &duplicRow, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   if(myRank == 0)
   {
      if(system_type==EQUALITY_SYSTEM)
         cout<<"There are "<<duplicRow<<" duplicate rows in A."<<endl;
      else
         cout<<"There are "<<duplicRow<<" duplicate rows in C."<<endl;
   }
}

/* Compares the matrix coeffients in row i and j.
 * Returns true if they are the same or false if not.
 * The method assumes ordered rows, so it should only be called when it is
 * ensured that the order in each row is intact.
 * Empty rows can be considered as duplicate rows.
 */
bool StochPresolverParallelRows::compareCoefficients(SparseStorageDynamic& matrix, int i, int j) const
{
   assert( currNnzRow->elements()[i] == currNnzRow->elements()[j] );
   assert( i>=0 && i<matrix.m );
   assert( j>=0 && j<matrix.m );

   int rowLen = matrix.rowptr[i].end - matrix.rowptr[i].start;

   if( rowLen != matrix.rowptr[j].end - matrix.rowptr[j].start )
      return false;  // row entries can be in Amat or Bmat blocks
   if( rowLen == 0.0 )
      return true;  // empty rows

   bool allNegatedCoeffs = false;
   // compare first entry to determine the sign:
   int k_i = matrix.rowptr[i].start;
   int k_j = matrix.rowptr[j].start;
   if( matrix.jcolM[k_i] != matrix.jcolM[k_j] )
      return false;
   if( matrix.M[k_i] == - matrix.M[k_j] )
      allNegatedCoeffs = true;
   else if ( matrix.M[k_i] != matrix.M[k_j] )
      return false;

   for( int k=1; k<rowLen; k++ )
   {
      k_i = matrix.rowptr[i].start + k;
      k_j = matrix.rowptr[j].start + k;

      if( matrix.jcolM[k_i] != matrix.jcolM[k_j] )
         return false;
      if( (allNegatedCoeffs && matrix.M[k_i] != -matrix.M[k_j]) ||
            (!allNegatedCoeffs && matrix.M[k_i] != matrix.M[k_j]) )
         return false;
   }
   return true;
}
