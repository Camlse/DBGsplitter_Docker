//! [snippet1]

#include "build_dbg.hpp"



#ifdef EYTA_DEBUG
#include "debug.h"
#endif

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
build_dbg::build_dbg ()  : Tool ("build_dbg"), //give a name to our tool
                           nbConditions(-1), nbReplicates(-1)
{
    populateParser();
}


void build_dbg::populateParser () {
    //getParser()->push_front (new OptionNoParam (STR_PAR_SIMPLIFY, "Simplify the graph using GATB's error correction procedures (tested in genomics, not in transcriptomics). See https://github.com/GATB/gatb-core/blob/v1.4.1/gatb-core/src/gatb/debruijn/impl/Simplifications.hpp.",  false));
    getParser()->push_front (new OptionOneParam (STR_PAR_STEP_1_FOLDER, "Which is the workdir?",  false, "step1/"));
    getParser()->push_front (new OptionOneParam (STR_PAR_GRAPH_PREFIX, "Prefix of the name of the built files related to the graph",  false, "graph"));
    getParser()->push_front (new OptionOneParam (STR_PAR_BETA_MAX, "Maximum beta layer",  false, "100"));
    getParser()->push_front (new OptionOneParam (STR_PAR_SHORT_READS_MIN_ABUNDANCE, "An integer, k-mers present strictly less than this"
        "number of times in the dataset will be discarded. This is applied only on the short reads graph.",  false, "2"));
    getParser()->push_front (new OptionOneParam (STR_PAR_SHORT_READS_RELATIVE_CUTOFF, "A threshold. Edges which counts are relatively smaller than this threshold are removed. This is applied only on the short reads graph.",  false, "0.02"));
    getParser()->push_front (new OptionOneParam (STR_PAR_K, "K-mer size",  false, "41"));
    getParser()->push_front (new OptionOneParam (STR_PAR_SHORT_READS, "A text file containing the input reads. Each input read is described in a line containing 3 columns: 1) path to a fasta, fastq, or .gz file containing the reads; 2) condition of the input reads (int); 3) replicate number of the input reads (int). This file needs a header. Files can have the same condition and replicate (in case of paired-end reads, for example). Check https://gitlab.inria.fr/lishisoa/DBGSplitter/blob/master/tests/test_map2k2/MAP2K2 for an example.",  true));
}

void build_dbg::checkReadFiles() {
    string readsFile = getInput()->getStr(STR_PAR_SHORT_READS);
    bool header=true;
    ifstream input;
    openFileForReading(readsFile, input);

    for(string line; getline( input, line ); )
    {
        //parse header
        if (header) {
            header=false;
            continue;
        }

        //ignore empty lines
        if (line.size()==0)
            continue;

        //read the info
        stringstream ss;
        ss << line;

        //check if the path is ok
        string path;
        ss >> path;
        ifstream file;
        openFileForReading(path, file);
        if (!file.is_open()) {
            stringstream ss;
            ss << "Error opening file " << path << " in " << readsFile << endl;
            fatalError(ss.str());
        }
        file.close();

        //check if condition is fine
        int cond;
        if (!(ss >> cond)) {
            stringstream error;
            error << "Conditions must be an integer. Offending line in " << readsFile << ": " << line;
            fatalError(error.str());
        }
        if (cond > nbConditions) nbConditions = cond;

        //check if replicate is fine
        int repl;
        if (!(ss >> repl)) {
            stringstream error;
            error << "Replicates must be an integer. Offending line in " << readsFile << ": " << line;
            fatalError(error.str());
        }

        if (repl > nbReplicates) nbReplicates = repl;

        //add it to condRepl2Reads
        string condRepl;
        {stringstream condReplSS; condReplSS << "C" << cond << "_" << "R" << repl; condRepl = condReplSS.str(); }
        condRepl2Reads[condRepl].push_back(path);
    }
    input.close();

    cout << "Summary of input conditions, replicates and files: " << endl;
    cout << "Conditions: " << nbConditions << endl;
    cout << "Replicates: " << nbReplicates << endl;
    int index=0;
    for (const auto &pair : condRepl2Reads) {
        cout << pair.first << ":";
        for (const auto &file : pair.second) {
            cout << file << ", ";
            reads2CondReplIndex[file]=index;
        }
        cout << endl;
        index++;
    }
}


//create the reads file given a condRepl2Reads
void build_dbg::createReadsFile(const string &readsFile) {
    ofstream fout;
    openFileForWriting(readsFile, fout);

    for (const auto &pair : condRepl2Reads) {
        for (const auto &file : pair.second) {
            auto boostPath(boost::filesystem::canonical(file));
            fout << boost::filesystem::canonical(file).string() << endl; //print the canonical path
        }
    }
    fout.close();
}

void buildSequence (
        const Graph& graph,
        const Node& startingNode,
        size_t length,
        size_t nbContigs,
        const string& consensusRight,
        const string& consensusLeft,
        Sequence& seq
)
{
    /** Shortcuts. */
    Data&  data     = seq.getData();

    /** We set the sequence comment. */
    stringstream ss1;
    ss1 << nbContigs << "__len__" << length << " ";
    seq._comment = ss1.str();

    /** We set the data length. */
    seq.getData().resize (length);

    //We fill the data
    string finalSequence = consensusLeft + graph.toString (startingNode) + consensusRight;
    for (size_t i=0;i<finalSequence.size();i++)
        data[i] = finalSequence[i];
}

//TODO: there is a bug with queryAbundance here!!
double computeUnitigCoverage(const gatb::core::debruijn::impl::Graph& graph, Node &startingNode, const Path& consensusRight, const Path& consensusLeft) {
    vector<CountNumber> kmerCoverages;
    kmerCoverages.push_back(graph.queryAbundance(startingNode));
    //kmerCoverages.push_back(startingNode.abundance);
    auto currentNode = startingNode;

    for_each(consensusRight.path.begin(), consensusRight.path.end(), [&](const Nucleotide &nucleotide) {
        currentNode = graph.successor(currentNode, nucleotide);
        kmerCoverages.push_back(graph.queryAbundance(currentNode));
        //kmerCoverages.push_back(currentNode.abundance);
    });
    currentNode = graph.reverse(startingNode);
    for_each(consensusLeft.path.begin(), consensusLeft.path.end(), [&](const Nucleotide &nucleotide) {
        currentNode = graph.successor(currentNode, nucleotide);
        kmerCoverages.push_back(graph.queryAbundance(currentNode));
        //kmerCoverages.push_back(currentNode.abundance);
    });

    //compute the average
    double coverage=0;
    for(const auto &kmerCoverage : kmerCoverages)
        coverage+=kmerCoverage;
    coverage=coverage/kmerCoverages.size();

    return coverage;
}




//Represents a node and its level - if a node is visited in another BFS level, then it is a different node for us here
struct NodeLevel {
    Node node;
    int level;
    NodeLevel(const Node& node, int level):node(node),level(level){}

    bool operator==(const NodeLevel &that)const {
        return this->node==that.node && this->level==that.level;
    }
    bool operator<(const NodeLevel &that)const {
        if (this->node!=that.node)
            return this->node<that.node;
        return this->level<that.level;
    }
};

vector<size_t> computeBeta(const gatb::core::debruijn::impl::Graph& graph, const Node &s, int maxBetaLevel, bool inBeta) {
    //the explored nodes and the queue
    set<NodeLevel> exploredNodes;
    queue<NodeLevel> q;

    //add the source
    NodeLevel source(s, 0);
    q.push(source);
    exploredNodes.insert(source);

    while (q.empty()==false) {
        //get the node to be explored
        NodeLevel currentNode = q.front();
        q.pop();

        //compute next beta
        int betaLevel = currentNode.level+1;

        if (betaLevel>maxBetaLevel)
            break; //no need to continue - we reached the max

        //get the neighbours to be explored
        auto neighbours = (inBeta ? graph.predecessors(currentNode.node) : graph.successors(currentNode.node));
        neighbours.iterate([&](const Node &possibleNeighbour) {
            NodeLevel possibleNeighbourWithLevel(possibleNeighbour, betaLevel);
            if (exploredNodes.count(possibleNeighbourWithLevel)==0) {
                //unexplored neighbour for the moment - add it to the queue
                q.push(possibleNeighbourWithLevel);
                exploredNodes.insert(possibleNeighbourWithLevel);
            }
        });
    }

    //count the # of nodes per BFS layer
    vector<size_t> nbNodesPerBFSLayer(maxBetaLevel+1, 0);
    for (const auto &nodeLevel : exploredNodes)
        nbNodesPerBFSLayer[nodeLevel.level]++;


    //DEBUG
    //for (size_t level=0; level<nbNodesPerBFSLayer.size(); level++)
    //    cout << "Level: " << level << ", # nodes: " << nbNodesPerBFSLayer[level] << endl;

    return nbNodesPerBFSLayer;
}

vector<size_t> computeInBeta(const gatb::core::debruijn::impl::Graph& graph, const Node &s, int maxBetaLevel) {
    return computeBeta(graph, s, maxBetaLevel, true);
}
vector<size_t> computeOutBeta(const gatb::core::debruijn::impl::Graph& graph, const Node &s, int maxBetaLevel) {
    return computeBeta(graph, s, maxBetaLevel, false);
}

//computes the shannon entropy of infoVector
double computeEntropyCore(const vector<string> &infoVector) {
    //compute the sums
    map<string, int> infoToSums;
    for (const auto &info : infoVector)
        infoToSums[info]++;

    //compute the probs
    map<string, double> infoToProb;
    for (const auto &pair : infoToSums)
        infoToProb[pair.first]= ((double)(pair.second))/infoVector.size();

    //compute the entropy
    double entropy = 0;
    for (const auto &pair : infoToProb)
        entropy-=pair.second*log2(pair.second);

    return entropy;
}


/*
 * Here we compute the lowest entropy between all entropies taking as the basic info:
 * 1 nt (in the only possible frame) - can capture stuff like AAAAAAAAAAAAAAAAAA
 * 2 nts (in the two possible frames) - can capture stuff like ATATATATATATATATATATATAT
 * 3 nts (in the three possible frames) - can capture stuff like ATCATCATCATCATCATCATCATCATCATCATC
 * ...
 * 6 nts (in the six possible frames) - can capture stuff like ATCGGAATCGGAATCGGAATCGGAATCGGAATCGGAATCGGA
 * We do from 1 to 6 since microsatellites have usually this size: https://en.wikipedia.org/wiki/Microsatellite
*/
double computeMicrosatelliteEntropy(const string &seq) {
    double minEntropy = std::numeric_limits<double>::max();
    for (int size=1;size<=6;size++) { //for each size of microsatellite
        for (int frame=0;frame<size;frame++) { //for each frame
            //builds infoVector
            vector<string> infoVector;
            int i=frame;
            while (i+size<=seq.size()) {
                infoVector.push_back(seq.substr(i, size));
                i+=size;
            }
            double entropy = computeEntropyCore(infoVector);

            minEntropy=min(minEntropy, entropy);
        }
    }
    return minEntropy;
}

//function used to build the unitigs
void construct_linear_seqs (const gatb::core::debruijn::impl::Graph& graph, const string& prefix, bool computeStats, int maxBetaLevel, SolidKmerOracle* solidKmerOracle, bool verbose)
{
    using namespace gatb::core::debruijn::impl;
    using namespace gatb::core::tools::misc::impl;

    string linear_seqs_name = prefix+".unitigs";
    string unitigSizeFilename = linear_seqs_name+".size";
    string unitigsFeaturesFilename = linear_seqs_name+".features";



    ofstream unitigsFeaturesFile;
    if (computeStats) {
        openFileForWriting(unitigsFeaturesFilename, unitigsFeaturesFile);
        unitigsFeaturesFile << "unitig_id\tcompression_level\tcoverage\tentropy";
        for (size_t i = 0; i < maxBetaLevel + 1; i++)
            unitigsFeaturesFile << "\tOutBetaLevel_" << i;
        unitigsFeaturesFile << "\tMaxOutBeta";
        for (size_t i = 0; i < maxBetaLevel + 1; i++)
            unitigsFeaturesFile << "\tInBetaLevel_" << i;
        unitigsFeaturesFile << "\tMaxInBeta";
        unitigsFeaturesFile << "\tMaxBeta" << endl;
    }

    IBank* outputBank = new BankFasta (linear_seqs_name.c_str());
    LOCAL (outputBank);

    // We create a Terminator object - this will mark the nodes that are already built
    MPHFTerminator terminator (graph);

    // We create a BranchingTerminator object - this will mark the nodes where to stop the traversal
    BranchingTerminator branchingTerminator(graph);

    // We create a Traversal instance to traverse unitigs
    Traversal* traversal = Traversal::create (TRAVERSAL_UNITIG, graph, branchingTerminator);
    LOCAL (traversal);

    Path consensusRight;
    Path consensusLeft;
    Sequence seq (Data::ASCII);
    u_int64_t nbContigs = 0;
    BankFasta::setDataLineSize(0);

    //We loop through the nodes and build the unitigs
    GraphIterator<Node> graphIt = graph.iterator();
    ProgressGraphIterator<Node, ProgressTimerAndSystem> progressIt (graph.iterator(), "Graph: building unitigs");
    Iterator<Node> *it = NULL;
    if (verbose)
        it = &progressIt;
    else
        it = &graphIt;


    for (it->first(); !it->isDone(); it->next()) {
        auto &startingNode = it->item();

        if (terminator.is_marked(startingNode))
            continue;

        if (graph.isNodeDeleted(startingNode))
            continue;

        auto reversedNode = graph.reverse(startingNode);
        int lenRight = traversal->traverse (startingNode, DIR_OUTCOMING, consensusRight);
        int lenLeft = traversal->traverse (reversedNode, DIR_OUTCOMING, consensusLeft);
        int lenTotal = graph.getKmerSize() + lenRight + lenLeft;

        //mark the traversed nodes
        terminator.mark(startingNode);
        auto currentNode = startingNode;
        for_each(consensusRight.path.begin(), consensusRight.path.end(), [&](const Nucleotide &nucleotide) {
            currentNode = graph.successor(currentNode, nucleotide);
            terminator.mark(currentNode);
        });
        auto rightestNode = currentNode;
        currentNode = reversedNode;
        for_each(consensusLeft.path.begin(), consensusLeft.path.end(), [&](const Nucleotide &nucleotide) {
            currentNode = graph.successor(currentNode, nucleotide);
            terminator.mark(currentNode);
        });
        auto leftestNode = graph.reverse(currentNode);

        // We get the unitig strings
        string consensusLeftStr;
        {
            stringstream ss;
            ss << consensusLeft;
            consensusLeftStr=reverse_complement(ss.str());
        }

        string consensusRightStr;
        {
            stringstream ss;
            ss << consensusRight;
            consensusRightStr=ss.str();
        }


        /** We create the contig sequence. */
        buildSequence(graph, startingNode, lenTotal, nbContigs, consensusRightStr, consensusLeftStr, seq);

        if (solidKmerOracle) { //if it was given, use it
            //we add it to the oracle
            solidKmerOracle->addUnitigSequence(seq.toString(), leftestNode);
        }

        /** We add the sequence into the output bank. */
        outputBank->insert (seq);

        if (computeStats) {
            //compute the compression level (#nb of kmers inside this unitig)
            int compressionLevel=lenTotal-graph.getKmerSize()+1;

            //We calculate the coverage of this unitig, i.e. the average k-mer count
            double coverage = computeUnitigCoverage(graph, startingNode, consensusRight, consensusLeft);

            //compute the in-beta (i.e. the in-beta of the leftest node) and the out-beta (i.e. the out-beta of the rightest node)
            size_t maxBeta=0;
            size_t maxInBeta=0;
            size_t maxOutBeta=0;
            vector<size_t> inBetas;
            vector<size_t> outBetas;
            {
                boost::timer::cpu_timer timer;

                //compute all in and out betas
                inBetas = computeInBeta(graph, leftestNode, maxBetaLevel);
                outBetas = computeOutBeta(graph, rightestNode, maxBetaLevel);

                //and now the max in and out betas
                maxInBeta = *max_element(inBetas.begin(), inBetas.end());
                maxOutBeta = *max_element(outBetas.begin(), outBetas.end());

                //and now the max beta
                maxBeta=max(maxInBeta, maxOutBeta);
                auto nanoseconds = boost::chrono::nanoseconds(timer.elapsed().user + timer.elapsed().system);
                auto seconds = boost::chrono::duration_cast<boost::chrono::seconds>(nanoseconds);
                //cout << "[TIMER_BETA] " << seconds.count() << endl;
            }


            //compute the sequence entropy
            double entropy=0;
            {
                boost::timer::cpu_timer timer;
                entropy = computeMicrosatelliteEntropy(seq.toString());
                auto nanoseconds = boost::chrono::nanoseconds(timer.elapsed().user + timer.elapsed().system);
                auto seconds = boost::chrono::duration_cast<boost::chrono::seconds>(nanoseconds);
                //cout << "[TIMER_ENTROPY] " << seconds.count() << endl;
            }

            //print everything to unitigsFeaturesFile
            unitigsFeaturesFile << nbContigs << "\t" << compressionLevel << "\t" << coverage << "\t" << entropy;
            for (size_t i=0; i<maxBetaLevel+1; i++)
                unitigsFeaturesFile << "\t" << outBetas[i];
            unitigsFeaturesFile << "\t" << maxOutBeta;
            for (size_t i=0; i<maxBetaLevel+1; i++)
                unitigsFeaturesFile << "\t" << inBetas[i];
            unitigsFeaturesFile << "\t" << maxInBeta;
            unitigsFeaturesFile << "\t" << maxBeta << endl;
        }



        //increase the number of contigs
        nbContigs += 1;
    }

    //print the nb of contigs to unitigSizeFilename file
    ofstream unitigSizeFile;
    openFileForWriting(unitigSizeFilename, unitigSizeFile);
    unitigSizeFile << nbContigs;

    //close all files
    outputBank->flush ();
    unitigSizeFile.close();

    if (computeStats)
        unitigsFeaturesFile.close();
}






void build_dbg::buildUnitigs (gatb::core::debruijn::impl::Graph& graph, //the graph where to build the unitigs
                    const string& prefix, //the prefix to use for the files
                    SolidKmerOracle *solidKmerOracle, //the solid kmer oracle - will store the unitigs and allow us for query for solid kmers, if != NULL
                    /* parameters for computing the edge counts
                     * the edge count is the count of the (k+1)-mer in the read sets
                     * it is only computed if these 3 parameters are not null */
                    gatb::core::debruijn::impl::Graph *graphKP1, //the graph with k+1-mers
                    SolidKmerOracle *solidKmerOracleKP1, //the solid kmer oracle of the k+1-mer graph - allowing us to remove the FPs in the k+1-graph
                    vector<vector<int>> *edgeIndexKP1_2_CountVector, //the edge index to count vector
                    /* other parameters */
                    bool computeStats, //compute the stats to use ML to infer which nodes are due to repeats
                    int maxBetaLevel, //the maximum beta level
                    bool verbose) //verbosity level
{
    //compute the unitigs of the graph
    int kmerSize = graph.getKmerSize();
    string linear_seqs_name = prefix+".unitigs";
    construct_linear_seqs(graph, prefix, computeStats, maxBetaLevel, solidKmerOracle, verbose);

    //builds and outputs .nodes and .edges.dbg files
    //TODO: go until GraphOutput<KMER_SPAN(7)>?
    typedef boost::variant <
        GraphOutput<KMER_SPAN(0)>,
        GraphOutput<KMER_SPAN(1)>,
        GraphOutput<KMER_SPAN(2)>,
        GraphOutput<KMER_SPAN(3)>
    >  GraphOutputVariant;

    GraphOutputVariant graphOutput;
    //TODO: check this
    if (kmerSize < KMER_SPAN(0))  {       graphOutput = GraphOutput<KMER_SPAN(0)>(&graph, solidKmerOracle, prefix, graphKP1, solidKmerOracleKP1, edgeIndexKP1_2_CountVector, condRepl2Reads.size()); }
    else if (kmerSize < KMER_SPAN(1))  {  graphOutput = GraphOutput<KMER_SPAN(1)>(&graph, solidKmerOracle, prefix, graphKP1, solidKmerOracleKP1, edgeIndexKP1_2_CountVector, condRepl2Reads.size()); }
    else if (kmerSize < KMER_SPAN(2))  {  graphOutput = GraphOutput<KMER_SPAN(2)>(&graph, solidKmerOracle, prefix, graphKP1, solidKmerOracleKP1, edgeIndexKP1_2_CountVector, condRepl2Reads.size()); }
    else if (kmerSize < KMER_SPAN(3))  {  graphOutput = GraphOutput<KMER_SPAN(3)>(&graph, solidKmerOracle, prefix, graphKP1, solidKmerOracleKP1, edgeIndexKP1_2_CountVector, condRepl2Reads.size()); }
    else { throw gatb::core::system::Exception ("Graph failure because of unhandled kmer size %d", kmerSize); }
    boost::apply_visitor (EdgeConstructionVisitor(linear_seqs_name, true),  graphOutput);

    #ifdef EYTA_DEBUG
    //EYTATester::testSolidKmerOracle(*solidKmerOracle);
    #endif

    //save disk space
    boost::filesystem::remove(linear_seqs_name);
}


/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void build_dbg::execute ()
{
    //check read files
    checkReadFiles();

    //read parameters
    int kmerSize = getInput()->getInt(STR_PAR_K);
    int maxBeta = getInput()->getInt(STR_PAR_BETA_MAX);
    int minAbundance = getInput()->getInt(STR_PAR_SHORT_READS_MIN_ABUNDANCE);
    string step1Folder = getInput()->getStr(STR_PAR_STEP_1_FOLDER);
    string prefix = getInput()->getStr(STR_PAR_GRAPH_PREFIX);
    prefix = step1Folder + "/" + prefix;
    int nbCores = getInput()->getInt(STR_PAR_NB_CORES);
    double relativeCutoff = 0;
    if (getInput()->get(STR_PAR_SHORT_READS_RELATIVE_CUTOFF))
        relativeCutoff = getInput()->getDouble(STR_PAR_SHORT_READS_RELATIVE_CUTOFF);

    /*
     * Removing the simplify parameter
    bool simplify = getInput()->get(STR_PAR_SIMPLIFY) != 0;
     */

    //create the reads file
    string readsFile(prefix+string(".readsFile"));
    createReadsFile(readsFile);


    /*
    [graph options]

   [kmer count options]
          -in                             (1 arg) :    reads file
          -kmer-size                      (1 arg) :    size of a kmer  [default '31']
          -abundance-min                  (1 arg) :    min abundance threshold for solid kmers  [default '2']
          -abundance-max                  (1 arg) :    max abundance threshold for solid kmers  [default '2147483647']
          -abundance-min-threshold        (1 arg) :    min abundance hard threshold (only used when min abundance is "auto")  [default '2']
          -histo-max                      (1 arg) :    max number of values in kmers histogram  [default '10000']
          -solidity-kind                  (1 arg) :    way to compute counts of several files (sum, min, max, one, all, custom)  [default 'sum']
          -solidity-custom                (1 arg) :    when solidity-kind is cutom, specifies list of files where kmer must be present  [default '']
          -max-memory                     (1 arg) :    max memory (in MBytes)  [default '5000']
          -max-disk                       (1 arg) :    max disk   (in MBytes)  [default '0']
          -solid-kmers-out                (1 arg) :    output file for solid kmers (only when constructing a graph)  [default '']
          -out                            (1 arg) :    output file  [default '']
          -out-dir                        (1 arg) :    output directory  [default '.']
          -out-tmp                        (1 arg) :    output directory for temporary files  [default '.']
          -out-compress                   (1 arg) :    h5 compression level (0:none, 9:best)  [default '0']
          -storage-type                   (1 arg) :    storage type of kmer counts ('hdf5' or 'file')  [default 'hdf5']

      [kmer count, algorithmic options options]
             -minimizer-type   (1 arg) :    minimizer type (0=lexi, 1=freq)  [default '0']
             -minimizer-size   (1 arg) :    size of a minimizer  [default '10']
             -repartition-type (1 arg) :    minimizer repartition (0=unordered, 1=ordered)  [default '0']

   [bloom options]
          -bloom        (1 arg) :    bloom type ('basic', 'cache', 'neighbor')  [default 'neighbor']
          -debloom      (1 arg) :    debloom type ('none', 'original' or 'cascading')  [default 'cascading']
          -debloom-impl (1 arg) :    debloom impl ('basic', 'minimizer')  [default 'minimizer']

   [branching options]
          -branching-nodes (1 arg) :    branching type ('none' or 'stored')  [default 'stored']
          -topology-stats  (1 arg) :    topological information level (0 for none)  [default '0']

   [general options]
          -config-only       (0 arg) :    dump config only
          -nb-cores          (1 arg) :    number of cores  [default '0']
          -verbose           (1 arg) :    verbosity level  [default '1']
          -integer-precision (1 arg) :    integers precision (0 for optimized value)  [default '0']
     */

    //Builds the SR DBG using GATB
    //TODO: by using create() and assigning to a Graph object, the copy constructor does a shallow or deep copy??
    //TODO: bug - if you run with -abundance-min 100, where we have 0 solid kmers, this bugs
    cerr << "Building Short Reads graph..." << endl;
    Graph* graph = gatb::core::debruijn::impl::Graph::createAsPointer(
        "-in %s -kmer-size %d -abundance-min %d -out %s -nb-cores %d",
        readsFile.c_str(), kmerSize, minAbundance, prefix.c_str(), nbCores);

    //computes the edge list for each node
    graph->precomputeAdjacency(nbCores, true);

    /*
     * Removing the simplify parameter
    //simplify the graph if the parameter is set
    if (simplify) {
        //TODO: finish this implementation, so that we have a better simplification for RNA-seq data
        //Remove tips/bubbles that are SNPs or indels of length < 3 NT
        //graph->simplifySNPsAndIndelsLessThan3NT(nbCores, true);
        cerr << "Simplifying the graph..." << endl;
        graph->simplify(nbCores, true);
        cerr << "Simplifying the graph - Done!" << endl;
    }
      */

    //do the error removal by removing edges from the precomputed adjacency list
    if (relativeCutoff > 0)
        graph->relativeErrorRemoval(relativeCutoff, prefix, nbCores);


    //get the counts for each edge (k+1-mer)
    cerr << "Getting edges count..." << endl;

    //TODO: makes it better - like not building the k+1 graph?
    cout << "Getting all k+1-mers..." << endl;
    Graph* graphKPlus1 = gatb::core::debruijn::impl::Graph::createAsPointer(
        "-in %s -kmer-size %d -abundance-min 0 -out %s_k_plus_1 -nb-cores %d -verbose 0",
        readsFile.c_str(), kmerSize+1, prefix.c_str(), nbCores);
    //build the unitigs of the graph
    //TODO: in fact I just need to fill in the SolidKmerOracle so that I know which edges correspond to true (k+1)-mers
    //TODO: so we are doing more than we need here...
    string prefixKPlus1 = prefix + "_kp1";
    SolidKmerOracle solidKmerOracleKPlus1(graphKPlus1, prefixKPlus1);
    buildUnitigs(*graphKPlus1, prefixKPlus1, &solidKmerOracleKPlus1);
    cout << "Getting all k+1-mers... - Done!" << endl;


    vector<string> allReadFilesNames = readFileAsVectorString(readsFile);
    vector<vector<int>> edgeIndexKP1_2_CountVector(
        (size_t)(graphKPlus1->getInfo()["kmers_nb_solid"]->getInt()),
        vector<int>(condRepl2Reads.size(), 0)); //TODO: memory issues: change vector<int> for a int[]?


    int readFileIndex=0;
    for (const string &readFileName : allReadFilesNames) {
        //TODO: CHANGE THIS - WE JUST NEED TO COMPUTE THE K+1 KMER COUNTS (DSK), NOT THE WHOLE GRAPH...
        cerr << "Getting edge count for read file " << readFileName << "..." << endl;
        Graph graphEdgesCount = gatb::core::debruijn::impl::Graph::create(
            "-in %s -kmer-size %d -abundance-min 0 -out %s_edge_count_%d -nb-cores %d -verbose 0",
            readFileName.c_str(), kmerSize+1, prefix.c_str(), readFileIndex, nbCores);

        //we now go through all kmers of the graph and store their count
        graphEdgesCount.iterator().iterate([&](Node &node){
            string nodeSeq = graphEdgesCount.toString(node);
            auto nodeIngraphKPlus1 = graphKPlus1->buildNode(nodeSeq.c_str());
            auto index = graphKPlus1->nodeMPHFIndex(nodeIngraphKPlus1);

            //we update the count
            edgeIndexKP1_2_CountVector[index][ reads2CondReplIndex[readFileName] ] += node.abundance;
        });

        ++readFileIndex;
    }

    //we output the edge counts
    //TODO: REMOVE ME AFTER EVERYTHING SEEMS FINE...
    //TODO: REMOVE ME AFTER EVERYTHING SEEMS FINE...
    /*
    string edgeCountFilename=prefix+string(".edges.counts");
    ofstream edgeCountFile;
    openFileForWriting(edgeCountFilename, edgeCountFile);
    for (const auto& pair : edgeIndexKP1_2_CountVector) {
        edgeCountFile << pair.first << "\t";
        for (const auto &count : pair.second)
            edgeCountFile << count << "\t";
        edgeCountFile << endl;
    }
    edgeCountFile.close();
     */
    //TODO: REMOVE ME AFTER EVERYTHING SEEMS FINE...
    //TODO: REMOVE ME AFTER EVERYTHING SEEMS FINE...
    cerr << "Getting edges count - Done!" << endl;


    //build the unitigs of the graph
    SolidKmerOracle solidKmerOracle(graph, prefix);
    buildUnitigs(*graph, prefix, &solidKmerOracle, graphKPlus1, &solidKmerOracleKPlus1, &edgeIndexKP1_2_CountVector, true, maxBeta, true);
    //buildUnitigs(*graph, prefix, NULL, graphKPlus1, &solidKmerOracleKPlus1, &edgeIndexKP1_2_CountVector, true, maxBeta, true);

    //save solidKmerOracle to the disk
    //TODO: delete the graph, since the SKO is responsible to doing so
    cerr << "Saving " << prefix << ".sko ..." << endl;
    std::ofstream solidKmerOracleOF(prefix+".sko");
    boost::archive::text_oarchive textOArchive(solidKmerOracleOF);
    textOArchive << solidKmerOracle;
    cerr << "Saving " << prefix << ".sko ... - Done!" << endl;


    cerr << "Saving additional files..." << endl;
    //save condRepl
    {
        std::ofstream condRepl2ReadsOF;
        openFileForWriting(prefix + ".cond_repl", condRepl2ReadsOF);
        for (const auto &condRepl2ReadsPair : condRepl2Reads)
            condRepl2ReadsOF << condRepl2ReadsPair.first << endl;
        condRepl2ReadsOF.close();
    }
    //save nbReplicates
    {
        std::ofstream nbReplicatesOF;
        openFileForWriting(prefix + ".nb_replicates", nbReplicatesOF);
        nbReplicatesOF << nbReplicates;
        nbReplicatesOF.close();
    }
    cerr << "Saving additional files - Done!" << endl;

    cerr << "Building Short Reads graph - Done!" << endl;
}
