#include "identify_nodes_ov_repeats.hpp"
#include "Utils.h"
#include "global.h"
using namespace std;

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
// The Tool constructor allows to give a name to our tool.
// This name appears when one gets help in the command line or in the final output
IdentifyNodesOverlappingRepeats::IdentifyNodesOverlappingRepeats ()  : Tool ("IdentifyNodesOverlappingRepeats") //give a name to our tool
{
  populateParser();
}

void IdentifyNodesOverlappingRepeats::populateParser() {
  getParser()->push_front (new OptionOneParam (STR_PAR_GRAPH_PREFIX, "Prefix of the name of the built files related to the graph",  false, "graph"));
  getParser()->push_front (new OptionOneParam (STR_PAR_USE_REF_TO_IDENTIFY_REPEAT_NODES, "Repeat identification mode (reference or de-novo)",  false, "reference"));
  getParser()->push_front (new OptionOneParam (STR_PAR_STARLONG_PATH, "Path to STARlong binary",  false, "/data2/leandro/repeats_on_transcriptome/STAR_2.5.3a/STAR-2.5.3a/bin/Linux_x86_64_static/STARlong"));
  getParser()->push_front (new OptionOneParam (STR_PAR_STARGENOME_DIR, "Path to STAR genome dir",  false, "/data2/leandro/repeats_on_transcriptome/ref_genome/STAR_index"));
  getParser()->push_front (new OptionOneParam (STR_PAR_BEDTOOLS_PATH, "Path to bedtools",  false, "/data2/leandro/repeats_on_transcriptome/bedtools/bedtools2/bin/bedtools"));
  getParser()->push_front (new OptionOneParam (STR_PAR_REPEAT_MASKER_BED_TRACK_PATH, "Path to Repeat Masker bed track",  false, "/data2/leandro/repeats_on_transcriptome/ref_repeats_track/hg38_repeatmasker_tracks.bed"));
  getParser()->push_front (new OptionOneParam (STR_PAR_REPEAT_INTERSECTION_THRESHOLD, "The minimum size of an intersection of a node and a repeat to consider the node as a repeat",  false, "20"));
}


/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void IdentifyNodesOverlappingRepeats::execute ()
{
  //get pars
  string dirWhereToolIsInstalled = getDirWhereToolIsInstalled();
  string prefix = getInput()->getStr(STR_PAR_GRAPH_PREFIX);
  string referenceOrDenovo = getInput()->getStr(STR_PAR_USE_REF_TO_IDENTIFY_REPEAT_NODES);
  string starlong = getInput()->getStr(STR_PAR_STARLONG_PATH);
  string stargenome = getInput()->getStr(STR_PAR_STARGENOME_DIR);
  string bedtools = getInput()->getStr(STR_PAR_BEDTOOLS_PATH);
  string rmBedTrack = getInput()->getStr(STR_PAR_REPEAT_MASKER_BED_TRACK_PATH);
  string repeatIntersectionThreshold = getInput()->getStr(STR_PAR_REPEAT_INTERSECTION_THRESHOLD);
  int nbCores = getInput()->getInt(STR_PAR_NB_CORES);


  if (referenceOrDenovo=="reference")  {
    //builds the command line
    string commandLine;
    {
      stringstream commandLineSS;
      commandLineSS << "bash " << dirWhereToolIsInstalled << "/pipeline.sh ";
      commandLineSS << prefix << " ";
      commandLineSS << starlong << " ";
      commandLineSS << stargenome << " ";
      commandLineSS << bedtools << " ";
      commandLineSS << rmBedTrack << " ";
      commandLineSS << repeatIntersectionThreshold << " ";
      commandLineSS << nbCores << " ";
      commandLine = commandLineSS.str();
    }
    executeCommand(commandLine);
  }else if (referenceOrDenovo=="de-novo") {
    fatalError("De-novo repeat identification mode not yet implemented");
  }else {
    fatalError("Unknown value for parameter -rep-id-mode");
  }
}