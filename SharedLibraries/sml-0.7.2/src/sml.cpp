/* (c) 2008,2009 Jonathan Hogg and Andreas Grothey, University of Edinburgh
 *
 * This file is part of SML.
 *
 * SML is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, using version 3 of the License.
 *
 * SML is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 */

#include "sml.h"
#include "AmplModel.h"
#include "backend.h"
#include "ExpandedModel.h"
#include "GlobalVariables.h"
#include <iostream>
#include <sys/stat.h>

#ifdef HAVE_DIRECT_H
#include <direct.h> // for mkdir() under MinGW
#endif

using namespace std;

/* global variables */
string GlobalVariables::datafilename = "";
const string GlobalVariables::amplcommand = "ampl";

bool GlobalVariables::logParseModel = false;
PrintLevelValues GlobalVariables::prtLvl = PRINT_NONE;

extern int yydebug;
void parse_data(AmplModel*, const string& datafilename);
int parse_model(const string& modelfilename);

string sml_version() {
  return PACKAGE_NAME " " PACKAGE_VERSION;
}

static void writeCopyright(ostream &out) {
   out << "Structure-conveying Modelling Language, " << sml_version() << endl;
   out << "(c) 2008,2009 Jonathan Hogg and Andreas Grothey, "
      "University of Edinburgh." << endl;
   out << "Released under LGPL v3" << endl;
}

static int createTmpDirIfNotPresent() {

   int err = 0;
   bool fl_exists;
   struct stat sbuf;
   const char *tmpdirname = "tmp";

   fl_exists = !(stat(tmpdirname, &sbuf) == -1);

   // tmp doesn't exist
   if (!fl_exists){
#ifdef HAVE_DIRECT_H
      err = mkdir(tmpdirname);
#else
      err = mkdir(tmpdirname, S_IRWXU);
#endif
      if (err)
        cerr << "ERROR: Failed to create temporary directory '"
             << tmpdirname << "'.\n";
   }

   // tmp is present
   else {

     // check that it's a directory with RWX permissions
     bool isusable = S_ISDIR(sbuf.st_mode) &&
                     ((sbuf.st_mode & S_IRWXU) == S_IRWXU);

     if (!isusable) {
       err = 1;
       cerr << "ERROR: Cannot use '" << tmpdirname
            << "' as temporary directory.\n";
     }
   }

   return err;
}

/* ----------------------------------------------------------------------------
main
---------------------------------------------------------------------------- */
ExpandedModelInterface* sml_generate(const string& modelfilename,
                                     const string& datafilename, bool debug) {

   int errcode;

   if (debug)
     GlobalVariables::prtLvl = PRINT_INFO;

   // ensure that the directory for temporary files exists and can be used
   errcode = createTmpDirIfNotPresent();
   if (errcode)
     return NULL;

   GlobalVariables::datafilename = datafilename;

   cout << "Reading model file '" << modelfilename << "'..." << endl;
   int rv = parse_model(modelfilename);
   if (rv)
     return NULL;

#ifdef PARSE_DATA
   /** @todo At the moment we are using AMPL to process the data file (see
    *  process_model() in backend.cpp or expandSet() expandStagesOfComp()
    *  in StochModel.cpp). Eventually, we will do it ourselves.
    */
   cout << "Reading data file '" << datafilename << "'..." << endl;
   parse_data(AmplModel::root, datafilename);
#endif

   // change working directory back to tmp/ for processing model
   errcode = chdir("tmp");
   if (errcode){
      cerr << "ERROR: Failed to change working directory to 'tmp/'\n";
      return NULL;
   }

   AmplModel::root->addDummyObjective();
   AmplModel::root->dump("logModel.dat");

   /* Write out and run ampl on ExpandedModels */
   errcode = process_model(AmplModel::root, GlobalVariables::datafilename);
   if (errcode)
     return NULL;

   if (GlobalVariables::prtLvl >= PRINT_LOG)
      cout << "------------- Generate ExpandedModel tree ------------ \n";
   ExpandedModel *em = AmplModel::root->createExpandedModel("root", "");

   if (GlobalVariables::prtLvl >= PRINT_LOG)
      em->print();

   return em;
}

