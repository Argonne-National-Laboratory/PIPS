// $Id: CglTreeInfo.hpp 1033 2011-06-19 16:49:13Z stefan $
// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglTreeInfo_H
#define CglTreeInfo_H

#include "OsiCuts.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinHelperFunctions.hpp"
class CglStored;
/** Information about where the cut generator is invoked from. */

class CglTreeInfo {
public:
  /// The level of the search tree node 
  int level;
  /** How many times the cut generator was already invoked in this search tree
      node */
  int pass;
  /** The number of rows in the original formulation. Some generators may not
      want to consider already generated rows when generating new ones. */
  int formulation_rows;
  /** Options 
      1 - treat costed integers as important
      2 - switch off some stuff as variables semi-integer
      4 - set global cut flag if at root node
      8 - set global cut flag if at root node and first pass
      16 - set global cut flag and make cuts globally valid
      32 - last round of cuts did nothing - maybe be more aggressive
      64 - in preprocessing stage
      128 - looks like solution
      256 - want alternate cuts
      512 - in sub tree (i.e. parent model)
      1024 - in must call again mode or after everything mode
  */
  int options;
  /// Set true if in tree (to avoid ambiguity at first branch)
  bool inTree;
  /** Replacement array.  Before Branch and Cut it may be beneficial to strengthen rows
      rather than adding cuts.  If this array is not NULL then the cut generator can
      place a pointer to the stronger cut in this array which is number of rows in size.

      A null (i.e. zero elements and free rhs) cut indicates that the row is useless 
      and can be removed.

      The calling function can then replace those rows.
  */
  OsiRowCut ** strengthenRow;
  /// Optional pointer to thread specific random number generator
  CoinThreadRandom * randomNumberGenerator;
  /// Default constructor 
  CglTreeInfo ();
 
  /// Copy constructor 
  CglTreeInfo (
    const CglTreeInfo &);
  /// Clone
  virtual CglTreeInfo * clone() const;

  /// Assignment operator 
  CglTreeInfo &
    operator=(
    const CglTreeInfo& rhs);
  
  /// Destructor 
  virtual
    ~CglTreeInfo ();
  /// Take action if cut generator can fix a variable (toValue -1 for down, +1 for up)
  virtual bool fixes(int , int , int ,bool) {return false;}
  /** Initalizes fixing arrays etc - returns >0 if we want to save info
      0 if we don't and -1 if is to be used */
  virtual int initializeFixing(const OsiSolverInterface * ) {return 0;}
  
};

/** Derived class to pick up probing info. */
typedef struct {
  //unsigned int oneFixed:1; //  nonzero if variable to 1 fixes all
  //unsigned int sequence:31; //  variable (in matrix) (but also see cliqueRow_)
  unsigned int fixes;
} cliqueEntry;

class CglTreeProbingInfo : public CglTreeInfo {
public:
  /// Default constructor 
  CglTreeProbingInfo ();
  /// Constructor from model
  CglTreeProbingInfo (const OsiSolverInterface * model);
 
  /// Copy constructor 
  CglTreeProbingInfo (
    const CglTreeProbingInfo &);
  /// Clone
  virtual CglTreeInfo * clone() const;

  /// Assignment operator 
  CglTreeProbingInfo &
    operator=(
    const CglTreeProbingInfo& rhs);
  
  /// Destructor 
  virtual
    ~CglTreeProbingInfo ();
  OsiSolverInterface * analyze(const OsiSolverInterface & si, int createSolver=0);
  /** Take action if cut generator can fix a variable 
      (toValue -1 for down, +1 for up)
      Returns true if still room, false if not  */
  virtual bool fixes(int variable, int toValue, int fixedVariable,bool fixedToLower);
  /** Initalizes fixing arrays etc - returns >0 if we want to save info
      0 if we don't and -1 if is to be used */
  virtual int initializeFixing(const OsiSolverInterface * model) ;
  /// Fix entries in a solver using implications
  int fixColumns(OsiSolverInterface & si) const;
  /// Fix entries in a solver using implications for one variable
  int fixColumns(int iColumn, int value, OsiSolverInterface & si) const;
  /// Packs down entries
  int packDown();
  /// Generate cuts from implications
  void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
		    const CglTreeInfo info) const;
  /// Entries for fixing variables
  inline cliqueEntry * fixEntries() const
  { convert(); return fixEntry_;}
  /// Starts of integer variable going to zero
  inline int * toZero() const
  { convert(); return toZero_;}
  /// Starts of integer variable going to one
  inline int * toOne() const
  { convert(); return toOne_;}
  /// List of 0-1 integer variables
  inline int * integerVariable() const
  { return integerVariable_;}
  /// Backward look up
  inline int * backward() const
  { return backward_;}
  /// Number of variables
  inline int numberVariables() const
  { return numberVariables_;}
  /// Number of 0-1 variables
  inline int numberIntegers() const
  { return numberIntegers_;}
private:
  /// Converts to ordered
  void convert() const;
protected:
  /// Entries for fixing variables
  mutable cliqueEntry * fixEntry_;
  /// Starts of integer variable going to zero
  mutable int * toZero_;
  /// Starts of integer variable going to one
  mutable int * toOne_;
  /// List of 0-1 integer variables
  int * integerVariable_;
  /// Backward look up
  int * backward_;
  /// Entries for fixing variable when collecting
  mutable int * fixingEntry_;
  /// Number of variables
  int numberVariables_;
  /// Number of 0-1 variables
  int numberIntegers_;
  /// Maximum number in fixEntry_
  int maximumEntries_;
  /// Number entries in fixingEntry_ (and fixEntry_) or -2 if correct style
  mutable int numberEntries_;
};
inline int sequenceInCliqueEntry(const cliqueEntry & cEntry)
{ return cEntry.fixes&0x7fffffff;}
inline void setSequenceInCliqueEntry(cliqueEntry & cEntry,int sequence)
{ cEntry.fixes = sequence|(cEntry.fixes&0x80000000);}
inline bool oneFixesInCliqueEntry(const cliqueEntry & cEntry)
{ return (cEntry.fixes&0x80000000)!=0;}
inline void setOneFixesInCliqueEntry(cliqueEntry & cEntry,bool oneFixes)
{ cEntry.fixes = (oneFixes ? 0x80000000 : 0)|(cEntry.fixes&0x7fffffff);}

#endif
