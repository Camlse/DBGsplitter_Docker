//
// Created by Leandro Ishi Soares de Lima on 19/03/2019.
//

#include "differential_analysis.h"
#include <fstream>
#include "Utils.h"
#include <vector>
#include <string>
#include "global.h"

using namespace std;

namespace differential_analysis {

    void differential_analysis::populateParser() {
      getParser()->push_front (new OptionOneParam (STR_PAR_GRAPH_PREFIX, "Prefix of the name of the built files related to the graph",  false, "graph"));
    }

    differential_analysis::differential_analysis() : Tool("differential_analysis") //give a name to our tool
    {
      populateParser();
    }

    void differential_analysis::execute() {
      //get pars
      string dirWhereToolIsInstalled = getDirWhereToolIsInstalled();
      string prefix = getInput()->getStr(STR_PAR_GRAPH_PREFIX);
      int nbCores = getInput()->getInt(STR_PAR_NB_CORES);

      //read nbReplicates
      {
        std::ifstream nbReplicatesIF;
        openFileForReading(prefix + ".nb_replicates", nbReplicatesIF);
        nbReplicatesIF >> nbReplicates;
        nbReplicatesIF.close();
      }


      //read the files where we have to perform differential analysis
      vector <string> allDiffAnalysisFiles = readFileAsVectorString(prefix + string(".diff_analysis_files"));

      //execute all diff analyses
      for (const auto &file : allDiffAnalysisFiles) {
        string commandLine;
        {
          stringstream commandLineSS;
          commandLineSS << "Rscript " << dirWhereToolIsInstalled << "/run_diff_ana.R " << file << " " << nbReplicates;
          commandLine = commandLineSS.str();
        }
        executeCommand(commandLine);
      }
    }
}