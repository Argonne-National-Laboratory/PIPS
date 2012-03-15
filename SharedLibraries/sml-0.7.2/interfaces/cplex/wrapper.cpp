/* (c) 2009 Marco Colombo, University of Edinburgh
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

#include <iostream>
#include <fstream>
#include <cstring> // for strcmp()
#include "sml.h"
#include "sml-cplex.h"

using namespace std;

const string progname = "smlcplex";

void writeHelp(ostream& out, const string& programname) {
  out << "Syntax:\n";
  out << "   " << programname
      << " [OPTIONS] modelfile datafile\n\n";
  out << "Option summary:\n";
  out << " -d            Enables debug information when reading model file.\n";
  out << " --help        Displays this help information.\n";
  out << " modelfile     File containing SML model.\n";
  out << " datafile      File containing SML data." << endl;
}

int analyseOptions(int argc, char **argv, string &modelfilename,
                   string &datafilename, bool &debug) {

  int found = 0;

  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-d") == 0) {
      debug = true;
    }
    else if (strcmp(argv[i], "--help") == 0) {
      writeHelp(cout, progname);
      return 1;
    }
    else {
      switch (found) {
      case 0:
        // first proper argument is the model file to read
        modelfilename = argv[i];
        break;
      case 1:
        // next one is data file
        datafilename = argv[i];
        break;
      default:
        break;
      }
      found++;
    }
  }

  if (modelfilename == "" || datafilename == "") {
    cerr << "ERROR: all of modelfile and datafile must be supplied."
         << endl << endl;
    return 1;
  }

  return 0;
}

int main(int argc, char **argv) {

   string modelfilename = "";
   string datafilename = "";
   bool debug = false;

   cout << "SML CPLEX, " << sml_version() << endl;
   cout << "(c) 2009 Marco Colombo, University of Edinburgh" << endl << endl;

   int rv = analyseOptions(argc, argv, modelfilename, datafilename, debug);
   if (rv)
     return rv;

   ExpandedModelInterface *em = sml_generate(modelfilename, datafilename, debug);
   if (!em)
     return 1;

   rv = SML_CPLEX_driver(em);
   if (rv)
     return rv;

   // solution
   string solnfilename = "model.sol";
   ofstream outfile(("../" + solnfilename).c_str());
   em->outputSolution(outfile);
   outfile.close();
   cout << "Solution written to file '" << solnfilename << "'" << endl;

   // clean up
   delete em;

   return 0;
}
