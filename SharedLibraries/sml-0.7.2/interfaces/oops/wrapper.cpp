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

#include <iostream>
#include <fstream>
#include <cstring> // for strcmp() and strncmp()
#include <cstdlib> // for exit()
#include "sml.h"
#include "sml-oops.h"

using namespace std;

const string progname = "smloops";

void writeHelp(ostream& out, const string& programname) {
   out << "Syntax:" << endl;
   out << "   " << programname
       << " [OPTIONS] modelfile datafile\n\n";
   out << "Option summary:" << endl;
   out << " -d                  Enables debug information when reading model "
      "file." << endl;
   out << " --help              Displays this help information." << endl;
   out << " --output=outfile," << endl;
   out << "   -o outfile        Write solution to file outfile." << endl;
   out << " modelfile           File containing SML model." << endl;
   out << " datafile            File containing SML data." << endl;
}

int analyseOptions(int argc, char **argv, string &modelfilename,
                   string &datafilename, string &outfilename, bool &debug) {
   int found = 0;
   for (int i=1;i<argc;i++){
      if (strcmp(argv[i], "-d")==0){
         debug = true;
      }
      else if(strcmp(argv[i], "--help")==0) {
         writeHelp(cout, progname);
         exit(0);
      }
      else if(strcmp(argv[i], "-o")==0) {
         if(i+1==argc) {
            cerr << "-o supplied without filename" << endl;
            return 1;
         }
         outfilename = argv[++i];
         if(outfilename.at(0)=='-') {
            cerr << "-o supplied without filename" << endl;
            return 1;
         }
      }
      else if(strncmp(argv[i], "--output=", 9)==0) {
         outfilename = (argv[i]+9);
      }
      else if(*(argv[i]) == '-') {
         cerr << "Unrecognised option '" << argv[i] << "'" << endl;
         return 1;
      }
      else {
         if (found==0){
            // first proper argument is the model file to read
            modelfilename = argv[i];
            found++;
         }else if(found==1){
            // next one is data file
            datafilename = argv[i];
            found++;
         }else{
            cerr << "ERROR: too many filenames." << endl;
            return 1;
         }
      }
   }

   if(modelfilename=="" || datafilename=="") {
      cerr << "ERROR: both modelfile and datafile " 
         "must be supplied." << endl << endl;
      return 1;
   }

   return 0;
}

/* ----------------------------------------------------------------------------
main
---------------------------------------------------------------------------- */
int main(int argc, char **argv) {
   string modelfilename = "";
   string datafilename = "";
   string outfilename = "";
   bool debug = false;

   cout << "SML OOPS, " << sml_version() << endl;
   cout << "(c) 2009 Andreas Grothey and Jonathan Hogg, University of Edinburgh"
        << endl << endl;

   int rv = analyseOptions(argc, argv,
                           modelfilename, datafilename, outfilename, debug);
   if (rv)
     return rv;

   if(debug) {
      cout << "======================================================" << endl;
      cout << "----------------- Call OOPS generator ----------------" << endl;
   }

   ExpandedModelInterface *em = sml_generate(modelfilename,datafilename,debug);
   if (!em)
     return 1;

   cout << "Calling OOPS..." << endl;
   SML_OOPS_driver(em);

   if(outfilename!="") {
      ofstream outfile(("../" + outfilename).c_str());
      em->outputSolution(outfile);
      outfile.close();
      cout << "Solution written to file '" << outfilename << "'" << endl;
   }

   // clean up
   delete em;

   return 0;
}
