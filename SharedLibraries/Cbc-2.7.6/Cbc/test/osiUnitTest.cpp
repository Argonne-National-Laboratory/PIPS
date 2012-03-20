// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).
// $Id$

#include "CoinPragma.hpp"

#include <iostream>

#include "OsiUnitTests.hpp"
#include "OsiCbcSolverInterface.hpp"

using namespace OsiUnitTest;

//----------------------------------------------------------------
// to see parameter list, call unitTest -usage
//----------------------------------------------------------------

int main (int argc, const char *argv[])
{
/*
  Start off with various bits of initialisation that don't really belong
  anywhere else.

  Synchronise C++ stream i/o with C stdio. This makes debugging
  output a bit more comprehensible. It still suffers from interleave of cout
  (stdout) and cerr (stderr), but -nobuf deals with that.
 */
  std::ios::sync_with_stdio() ;
/*
  Suppress an popup window that Windows shows in response to a crash. See
  note at head of file.
 */
  WindowsErrorPopupBlocker();

/*
  Process command line parameters.
 */
  std::map<std::string,std::string> parms;
  if (processParameters(argc,argv,parms) == false)
    return 1;

  std::string mpsDir = parms["-mpsDir"] ;
  std::string netlibDir = parms["-netlibDir"] ;

  /*
    Test Osi{Row,Col}Cut routines.
   */
  {
    OsiCbcSolverInterface cbcSi;
    testingMessage( "Testing OsiRowCut with OsiCbcSolverInterface\n" );
    OSIUNITTEST_CATCH_ERROR(OsiRowCutUnitTest(&cbcSi,mpsDir), {}, cbcSi, "rowcut unittest");
  }
  {
    OsiCbcSolverInterface cbcSi;
    testingMessage( "Testing OsiColCut with OsiCbcSolverInterface\n" );
    OSIUNITTEST_CATCH_ERROR(OsiColCutUnitTest(&cbcSi,mpsDir), {}, cbcSi, "colcut unittest");
  }
  {
    OsiCbcSolverInterface cbcSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiCbcSolverInterface\n" );
    OSIUNITTEST_CATCH_ERROR(OsiRowCutDebuggerUnitTest(&cbcSi,mpsDir), {}, cbcSi, "rowcut debugger unittest");
  }

  /*
    Run the OsiCbc class test. This will also call OsiSolverInterfaceCommonUnitTest.
   */
  testingMessage( "Testing OsiCbcSolverInterface\n" );
  OSIUNITTEST_CATCH_ERROR(OsiCbcSolverInterfaceUnitTest(mpsDir,netlibDir), {}, "cbc", "osicbc unittest");

  /*
    We have run the specialised unit test.
    Check now to see if we need to run through the Netlib problems.
   */
  if (parms.find("-testOsiSolverInterface") != parms.end())
  {
    // Create vector of solver interfaces
    std::vector<OsiSolverInterface*> vecSi(1, new OsiCbcSolverInterface);

    testingMessage( "Testing OsiSolverInterface on Netlib problems.\n" );
    OSIUNITTEST_CATCH_ERROR(OsiSolverInterfaceMpsUnitTest(vecSi,netlibDir), {}, "cbc", "netlib unittest");

    delete vecSi[0];
  }
  else
    testingMessage( "***Skipped Testing of OsiCbcSolverInterface on Netlib problems, use -testOsiSolverInterface to run them.***\n" );

  /*
    We're done. Report on the results.
   */
  std::cout.flush();
  outcomes.print();

  int nerrors;
  int nerrors_expected;
  outcomes.getCountBySeverity(TestOutcome::ERROR, nerrors, nerrors_expected);

  if (nerrors > nerrors_expected)
    std::cerr << "Tests completed with " << nerrors - nerrors_expected << " unexpected errors." << std::endl ;
  else
    std::cerr << "All tests completed successfully\n";

  return nerrors - nerrors_expected;
}
