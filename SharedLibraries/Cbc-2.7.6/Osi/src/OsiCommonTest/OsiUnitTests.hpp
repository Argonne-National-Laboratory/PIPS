// Copyright (C) 2010
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef OSISOLVERINTERFACETEST_HPP_
#define OSISOLVERINTERFACETEST_HPP_

#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <map>

class OsiSolverInterface;
class CoinPackedVectorBase;

/** A function that tests that a lot of problems given in MPS files (mostly the NETLIB problems) solve properly with all the specified solvers.
 *
 * The routine creates a vector of NetLib problems (problem name, objective,
 * various other characteristics), and a vector of solvers to be tested.
 *
 * Each solver is run on each problem. The run is deemed successful if the
 * solver reports the correct problem size after loading and returns the
 * correct objective value after optimization.

 * If multiple solvers are available, the results are compared pairwise against
 * the results reported by adjacent solvers in the solver vector. Due to
 * limitations of the volume solver, it must be the last solver in vecEmptySiP.
 */
void OsiSolverInterfaceMpsUnitTest
  (const std::vector<OsiSolverInterface*> & vecEmptySiP,
   const std::string& mpsDir);

/** A function that tests the methods in the OsiSolverInterface class.
 * Some time ago, if this method is compiled with optimization,
 * the compilation took 10-15 minutes and the machine pages (has 256M core memory!)...
 */
void OsiSolverInterfaceCommonUnitTest
  (const OsiSolverInterface* emptySi,
   const std::string& mpsDir,
   const std::string& netlibDir);

/** A function that tests the methods in the OsiColCut class. */
void OsiColCutUnitTest
  (const OsiSolverInterface * baseSiP,
   const std::string & mpsDir);

/** A function that tests the methods in the OsiRowCut class. */
void OsiRowCutUnitTest
  (const OsiSolverInterface * baseSiP,
   const std::string & mpsDir);

/** A function that tests the methods in the OsiRowCutDebugger class. */
void OsiRowCutDebuggerUnitTest
  (const OsiSolverInterface * siP,
   const std::string & mpsDir);

/** A function that tests the methods in the OsiCuts class. */
void OsiCutsUnitTest();

namespace OsiUnitTest {

class TestOutcomes;

/** verbosity level of unit tests
 * 0 (= default) for minimal output */
extern unsigned int verbosity;

/** behaviour on failing test
 * 0 (= default) continue
 * 1 press any key to continue
 * 2 stop with abort()
 */
extern unsigned int haltonerror;

/** a global TestOutcomes object to store test outcomes during a run of the Osi unittest
 */
extern TestOutcomes outcomes;

void failureMessage(const std::string &solverName,
		    const std::string &message) ;
void failureMessage(const OsiSolverInterface &si,
		    const std::string &message) ;
void failureMessage(const std::string &solverName,
		    const std::string &testname, const std::string &testcond) ;
void failureMessage(const OsiSolverInterface &si,
		    const std::string &testname, const std::string &testcond) ;
void testingMessage(const char *const msg) ;

bool equivalentVectors(const OsiSolverInterface * si1,
		       const OsiSolverInterface * si2,
		       double tol,
		       const double * v1,
		       const double * v2,
		       int size) ;

bool compareProblems(OsiSolverInterface *osi1, OsiSolverInterface *osi2) ;

bool isEquivalent(const CoinPackedVectorBase &pv, int n, const double *fv) ;

/** Utility routine to process Osi unittest command line parameters.
 * An unrecognised parameter which keyword is not in the ignorekeywords set will trigger the help message and a return value of false.
 * The ignorekeywords mapping specifies for each keyword to ignore the number of following parameters that should be ignored.
 * This should be replaced with the one of the standard CoinUtils parameter mechanisms.
 */
bool processParameters (int argc, const char **argv, std::map<std::string,std::string>& parms, const std::map<std::string,int>& ignorekeywords=std::map<std::string,int>());

class TestOutcome {
public:
	typedef enum {
		NOTE     = 0,
		PASSED   = 1,
		WARNING  = 2,
		ERROR    = 3,
		LAST     = 4
	} SeverityLevel;

	static std::string SeverityLevelName[LAST];

	std::string     component;
	std::string     testname;
	std::string     testcond;
	SeverityLevel   severity;
	bool            expected;

	std::string     filename;
	int             linenumber;

	TestOutcome(const std::string& comp, const std::string& tst, const char* cond, SeverityLevel sev, const char* file, int line, bool exp = false)
	: component(comp), testname(tst), testcond(cond), severity(sev), expected(exp), filename(file), linenumber(line)
	{ }

	void print() const;
};

class TestOutcomes : public std::list<TestOutcome> {
public:
	void add(std::string comp, std::string tst, const char* cond, TestOutcome::SeverityLevel sev, const char* file, int line, bool exp = false)
	{
		push_back(TestOutcome(comp, tst, cond, sev, file, line, exp));
	}

	void add(const OsiSolverInterface& si, std::string tst, const char* cond, TestOutcome::SeverityLevel sev, const char* file, int line, bool exp = false);

	void print() const;

	/** given a severity level, computes the total and the expected number of outcomes for this severity level
	 * this is implemented as O(1) operation
	 */
	void getCountBySeverity(TestOutcome::SeverityLevel sev, int& total, int& expected) const;
};

#define OSIUNITTEST_QUOTEME_(x) #x
#define OSIUNITTEST_QUOTEME(x) OSIUNITTEST_QUOTEME_(x)

#define OSIUNITTEST_ADD_OUTCOME(component, testname, testcondition, severity, expected) \
	OsiUnitTest::outcomes.add(component, testname, testcondition, severity, __FILE__, __LINE__, expected)

#define OSIUNITTEST_ASSERT_SEVERITY_EXPECTED(condition, failurecode, component, testname, severity, expected) \
	{ \
		if( condition ) { \
			OSIUNITTEST_ADD_OUTCOME(component, testname, #condition, OsiUnitTest::TestOutcome::PASSED, false); \
		  if (OsiUnitTest::verbosity >= 2) { \
		  	std::string successmsg(__FILE__ ":" OSIUNITTEST_QUOTEME(__LINE__) ": "); \
			  successmsg = successmsg + testname; \
			  successmsg = successmsg + " (condition \'" #condition "\') passed.\n"; \
				OsiUnitTest::testingMessage(successmsg.c_str()); \
		  } \
		} else { \
			OSIUNITTEST_ADD_OUTCOME(component, testname, #condition, severity, expected); \
			OsiUnitTest::failureMessage(component, testname, #condition); \
      switch( OsiUnitTest::haltonerror ) { \
        case 2: if( severity >= OsiUnitTest::TestOutcome::ERROR ) abort(); break; \
        case 1: std::cout << std::endl << "press any key to continue..." << std::endl; getchar(); \
        default: ; } \
			failurecode; \
    } \
	}

#define OSIUNITTEST_ASSERT_ERROR(condition, failurecode, component, testname) \
	OSIUNITTEST_ASSERT_SEVERITY_EXPECTED(condition, failurecode, component, testname, OsiUnitTest::TestOutcome::ERROR, false)

#define OSIUNITTEST_ASSERT_WARNING(condition, failurecode, component, testname) \
	OSIUNITTEST_ASSERT_SEVERITY_EXPECTED(condition, failurecode, component, testname, OsiUnitTest::TestOutcome::WARNING, false)

#define OSIUNITTEST_CATCH_SEVERITY_EXPECTED(trycode, catchcode, component, testname, severity, expected) \
	{\
		try {\
			trycode; \
		  OSIUNITTEST_ADD_OUTCOME(component, testname, #trycode " did not threw exception", OsiUnitTest::TestOutcome::PASSED, false); \
		  if (OsiUnitTest::verbosity >= 2) { \
		  	std::string successmsg( __FILE__ ":" OSIUNITTEST_QUOTEME(__LINE__) ": "); \
			  successmsg = successmsg + testname; \
			  successmsg = successmsg + " (code \'" #trycode "\') did not throw exception.\n"; \
				OsiUnitTest::testingMessage(successmsg.c_str()); \
		  } \
		} catch (CoinError& e) { \
			std::stringstream errmsg; \
			errmsg << #trycode " threw CoinError: " << e.message(); \
			if (e.className().length() > 0) \
			  errmsg << " in " << e.className(); \
			if (e.methodName().length() > 0) \
			  errmsg << " in " << e.methodName(); \
			if (e.lineNumber() >= 0) \
			  errmsg << " at " << e.fileName() << ":" << e.lineNumber(); \
			OSIUNITTEST_ADD_OUTCOME(component, testname, errmsg.str().c_str(), severity, expected); \
			OsiUnitTest::failureMessage(component, testname, errmsg.str().c_str()); \
      switch( OsiUnitTest::haltonerror ) { \
        case 2: if( severity >= OsiUnitTest::TestOutcome::ERROR ) abort(); break; \
        case 1: std::cout << std::endl << "press any key to continue..." << std::endl; getchar(); \
        default: ; } \
			catchcode; \
		} catch (...) { \
		  std::string errmsg; \
		  errmsg = #trycode; \
		  errmsg = errmsg + " threw unknown exception"; \
		  OSIUNITTEST_ADD_OUTCOME(component, testname, errmsg.c_str(), severity, false); \
		  OsiUnitTest::failureMessage(component, testname, errmsg.c_str()); \
		  catchcode; \
		} \
	}

#define OSIUNITTEST_CATCH_ERROR(trycode, catchcode, component, testname) \
	OSIUNITTEST_CATCH_SEVERITY_EXPECTED(trycode, catchcode, component, testname, OsiUnitTest::TestOutcome::ERROR, false)

#define OSIUNITTEST_CATCH_WARNING(trycode, catchcode, component, testname) \
	OSIUNITTEST_CATCH_SEVERITY_EXPECTED(trycode, catchcode, component, testname, OsiUnitTest::TestOutcome::WARNING, false)

} // end namespace OsiUnitTest

#endif /*OSISOLVERINTERFACETEST_HPP_*/
