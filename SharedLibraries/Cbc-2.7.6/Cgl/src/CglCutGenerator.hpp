// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglCutGenerator_H
#define CglCutGenerator_H

#include "OsiCuts.hpp"
#include "OsiSolverInterface.hpp"
#include "CglTreeInfo.hpp"

//-------------------------------------------------------------------
//
// Abstract base class for generating cuts.
//
//-------------------------------------------------------------------
///
/** Cut Generator Base Class

This is an abstract base class for generating cuts.  A specific cut 
generator will inherit from this class.
*/
class CglCutGenerator  {
  
public:
    
  /**@name Generate Cuts */
  //@{
  /** Generate cuts for the model data contained in si.
  The generated cuts are inserted into and returned in the
  collection of cuts cs.
  */
  virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
			     const CglTreeInfo info = CglTreeInfo()) const=0; 
  //@}

    
  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglCutGenerator (); 
 
  /// Copy constructor 
  CglCutGenerator ( const CglCutGenerator &);

  /// Clone
  virtual CglCutGenerator * clone() const = 0;

  /// Assignment operator 
  CglCutGenerator & operator=(const CglCutGenerator& rhs);

  /// Destructor 
  virtual ~CglCutGenerator ();

  /** Create C++ lines to set the generator in the current state.
      The output must be parsed by the calling code, as each line
      starts with a key indicating the following:<BR>
      0: must be kept (for #includes etc)<BR>
      3: Set to changed (not default) values<BR>
      4: Set to default values (redundant)<BR>

      Keys 1, 2, 5, 6, 7, 8 are defined, but not applicable to 
      cut generators.
  */
  virtual std::string generateCpp( FILE * ) {return "";}

  /// This can be used to refresh any information
  virtual void refreshSolver(OsiSolverInterface * ) {}
  //@}
  
  /**@name Gets and Sets */
  //@{
  /**
     Get Aggressiveness - 0 = neutral, 100 is normal root node.
     Really just a hint to cut generator
  */
  inline int getAggressiveness() const
  { return aggressive_;}

  /**
     Set Aggressiveness - 0 = neutral, 100 is normal root node.
     Really just a hint to cut generator
  */
  inline void setAggressiveness(int value)
  { aggressive_=value;}

  /**
     Returns true if may generate Row cuts in tree (rather than root node).
     Used so know if matrix will change in tree.  Really
     meant so column cut generators can still be active
     without worrying code.
     Default is true
  */
  virtual bool mayGenerateRowCutsInTree() const;
  /// Return true if needs optimal basis to do cuts
  virtual bool needsOptimalBasis() const;
  /// Return maximum length of cut in tree
  virtual int maximumLengthOfCutInTree() const
  { return COIN_INT_MAX;}
  //@}
  
  // test this class
  //static void unitTest();
  
// private:
  
  /**
     Aggressiveness - 0 = neutral, 100 is normal root node.
     Really just a hint to cut generator
  */
  int aggressive_;
};

#endif
