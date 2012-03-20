// $Id: CglMessage.hpp 948 2011-01-04 23:28:13Z lou $
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglMessage_H
#define CglMessage_H


#include "CoinPragma.hpp"

// This deals with Cgl messages (as against Osi messages etc)

#include "CoinMessageHandler.hpp"
enum CGL_Message
{
  CGL_INFEASIBLE,
  CGL_CLIQUES,
  CGL_FIXED,
  CGL_PROCESS_STATS,
  CGL_SLACKS,
  CGL_PROCESS_STATS2,
  CGL_PROCESS_SOS1,
  CGL_PROCESS_SOS2,
  CGL_UNBOUNDED,
  CGL_ELEMENTS_CHANGED1,
  CGL_ELEMENTS_CHANGED2,
  CGL_MADE_INTEGER,
  CGL_ADDED_INTEGERS,
  CGL_POST_INFEASIBLE,
  CGL_POST_CHANGED,
  CGL_DUMMY_END
};

/** This deals with Cgl messages (as against Osi messages etc)
 */
class CglMessage : public CoinMessages {

public:

  /**@name Constructors etc */
  //@{
  /** Constructor */
  CglMessage(Language language=us_en);
  //@}

};

#endif
