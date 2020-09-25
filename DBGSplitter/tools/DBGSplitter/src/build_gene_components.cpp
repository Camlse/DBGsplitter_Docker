//
// Created by Leandro Ishi Soares de Lima on 16/03/18.
//

#include "build_gene_components.h"
using namespace std;

namespace build_gene_components {
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
    build_gene_components::build_gene_components ()  : Tool ("build_gene_components") //give a name to our tool
    {
      populateParser();
    }

    void build_gene_components::populateParser() {
      getParser()->push_front (new OptionOneParam (STR_PAR_STEP_1_FOLDER, "Step1 workdir",  false, "step1/"));
      getParser()->push_front (new OptionOneParam (STR_PAR_STEP_2_FOLDER, "Step2 workdir",  false, "step2/"));
      getParser()->push_front (new OptionOneParam (STR_PAR_STEP_3_FOLDER, "Step3 workdir",  false, "step3/"));
      getParser()->push_front (new OptionOneParam (STR_PAR_GRAPH_PREFIX, "Prefix of the name of the built files related to the graph",  false, "graph"));
      getParser()->push_front (new OptionOneParam (STR_PAR_K, "K-mer size",  false, "41"));
    }

    string nodeIdToStr(const graph_t &graph, const Vertex &node) {
      stringstream ss;
      ss << graph[node].id << "_" << graph[node].strand;
      return ss.str();
    }

    string nodeIdToStr(const graph_t &graph, const Vertex &node, char inOrOutEdges) {
      stringstream ss;
      ss << nodeIdToStr(graph, node) << "_" << inOrOutEdges;
      return ss.str();
    }

    //for each component, output 4 files:
    //.edges (cytoscape)
    //.nodes (cytoscape)
    //.inputR (in the mail)
    //.nodesToTest (39_F (source node) 12_R 10_F (target nodes)...)
    void build_gene_components::printComponent(const string &prefix, const graph_t &graph, const set<Vertex> &nodesInComponent,
                        int nbOfCondRepl, vector<string> &filesToDiffAnalysis) {
      if (nodesInComponent.size()<4) {
        cout << "Component " << prefix << " has less than 4 nodes and will not be output" << endl;
        return;
      }

      //create the .nodes file
      {
        ofstream cytoscapeNodesFile;
        openFileForWriting(prefix+string(".nodes"), cytoscapeNodesFile);
        cytoscapeNodesFile << "id\tseq\tlength" << endl;
        for (const auto &node : nodesInComponent)
          cytoscapeNodesFile << nodeIdToStr(graph, node) << "\t" << graph[node].name << "\t" << graph[node].weight << endl;
        cytoscapeNodesFile.close();
      }

      //create the .edges, .inputR and .nodesToTest files
      {
        //create the headers for each file
        ofstream cytoscapeEdgesFile;
        openFileForWriting(prefix+string(".edges"), cytoscapeEdgesFile);
        cytoscapeEdgesFile << "from\tto";
        for (const auto &condRepl : condRepls)
          cytoscapeEdgesFile << "\t" << condRepl;
        cytoscapeEdgesFile << endl;

        ofstream inputRFile;
        openFileForWriting(prefix+string(".inputR"), inputRFile);
        inputRFile << "source";
        for (int degree=1; degree<=4; degree++) {
          for (const auto &condRepl : condRepls)
            inputRFile << "\tquantif_tg_" << degree << "_" << condRepl;
        }
        inputRFile << "\tnb_targets" << endl;

        ofstream nodesToTestFile;
        openFileForWriting(prefix+string(".nodesToTest"), nodesToTestFile);
        nodesToTestFile << "source\ttargets" << endl;


        //define a function to print the nodes
        auto printNodes = [&](const Vertex &node, bool outDirection) {
            //get the neighbours of the node that are in the component
            vector<Edge> neighbourEdges;
            vector<Vertex> neighbourVertices;

            if (outDirection) {
              graph_t::out_edge_iterator ei, ei_end;
              for (boost::tie(ei, ei_end) = out_edges(node, graph); ei != ei_end; ++ei) {
                auto uvEdge = *ei;
                auto neighbour = target(uvEdge, graph);
                if (nodesInComponent.count(neighbour)){ //is the neighbour in the component?
                  //yes
                  neighbourEdges.push_back(uvEdge);
                  neighbourVertices.push_back(neighbour);
                }
              }
            }else {
              graph_t::in_edge_iterator ei, ei_end;
              for (boost::tie(ei, ei_end) = in_edges(node, graph); ei != ei_end; ++ei) {
                auto uvEdge = *ei;
                auto neighbour = source(uvEdge, graph);
                if (nodesInComponent.count(neighbour)){ //is the neighbour in the component?
                  //yes
                  neighbourEdges.push_back(uvEdge);
                  neighbourVertices.push_back(neighbour);
                }
              }
            }

            int degree = neighbourVertices.size();
            char inOrOutEdges = (outDirection ? 'O' : 'I');
            bool printToRFiles = (degree >= 2);

            //print the source to inputRFile and nodesToTestFile
            if (printToRFiles) {
              inputRFile << nodeIdToStr(graph, node, inOrOutEdges);
              nodesToTestFile << nodeIdToStr(graph, node, inOrOutEdges);
            }

            //main print of the 3 files
            for (int neighbourI=0;neighbourI<neighbourEdges.size();neighbourI++) {
              auto uvEdge = neighbourEdges[neighbourI];
              auto neighbour = neighbourVertices[neighbourI];

              //print the edge to cytoscapeEdgesFile
              //there is no need to print this two times, thus we will only print this when outDirection is true
              if (outDirection) {
                cytoscapeEdgesFile << nodeIdToStr(graph, node) << "\t" << nodeIdToStr(graph, neighbour);
                for (int count : graph[uvEdge].countVector)
                  cytoscapeEdgesFile << "\t" << count;
                cytoscapeEdgesFile << endl;
              }


              //print the edge counts to inputRFile
              if (printToRFiles) {
                for (int count : graph[uvEdge].countVector)
                  inputRFile << "\t" << count;
              }

              //print the target to nodesToTestFile
              if (printToRFiles)
                nodesToTestFile << "\t" << nodeIdToStr(graph, neighbour);
            }

            //print the rest of the inputRFile
            if (printToRFiles) {
              for (int degreeIndex=degree;degreeIndex<4;degreeIndex++) {
                for (int i=0; i<nbOfCondRepl; i++)
                  inputRFile << "\t-";
              }
              inputRFile << "\t" << degree << endl;
            }

            //print the rest of the nodesToTestFile
            if (printToRFiles) {
              nodesToTestFile << endl;
            }
        };

        for (const auto &node : nodesInComponent) {
          printNodes(node, true); //process the out edges of the node
          printNodes(node, false);//process the in edges of the node
        }

        //add to filesToDiffAnalysis
        filesToDiffAnalysis.push_back(prefix+string(".inputR"));

        cytoscapeEdgesFile.close();
        inputRFile.close();
        nodesToTestFile.close();
      }
    }



    /*
     * TODO: I THINK THESE ARE NOT NEEDED ANYMORE
    //recursively finds all paths between s and t
    void findAllSTPathsCore(const FilteredGraph &graph, Vertex s, Vertex t, vector< vector<Vertex> > &allPaths, vector<Vertex> &currentPath, set<Vertex> &explored) {
      //add s to the explored nodes
      explored.insert(s);

      //add this path to allPaths if t is the last node
      if (s==t)
        allPaths.push_back(currentPath);
      else {
        FGAdjacencyIterator vi, vi_end;
        for (boost::tie(vi, vi_end) = adjacent_vertices(s, graph); vi != vi_end; ++vi) {
          Vertex v = *vi;

          //check if v should be explored
          if (explored.count(v)==0) { //if v is not already explored
            //yes, v should be explored
            //Configure the currentPath accordingly
            currentPath.push_back(v);

            //call DFS
            findAllSTPathsCore(graph, v, t, allPaths, currentPath, explored);

            //DesConfigure the currentPath accordingly
            currentPath.pop_back();
          }
        }
      }

      //already explored all paths with this prefix
      explored.erase(s);
    }

    //function that finds all paths between s and t
    //sets some variables and call findAllSTPathsCore
    vector< vector<Vertex> > findAllSTPaths(const FilteredGraph &graph, Vertex s, Vertex t) {
      vector< vector<Vertex> > allPaths; //will contain all paths from s to t
      vector<Vertex> currentPath;
      set<Vertex> explored; //will keep track of the nodes that were already explored on currentPath

      currentPath.push_back(s);
      findAllSTPathsCore(graph, s, t, allPaths, currentPath, explored);
      return allPaths;
    }
     */


    //recursively finds all paths between s and and any target t
    void findAllSTPathsCore(const FilteredGraph &graph, Vertex s, vector< vector<Vertex> > &allPaths, vector<Vertex> &currentPath, set<Vertex> &explored) {
      //add s to the explored nodes
      explored.insert(s);

      //add this path to allPaths if s is a target
      if (out_degree(s, graph)==0)
        allPaths.push_back(currentPath);
      else {
        FGAdjacencyIterator vi, vi_end;
        for (boost::tie(vi, vi_end) = adjacent_vertices(s, graph); vi != vi_end; ++vi) {
          Vertex v = *vi;

          //check if v should be explored
          if (explored.count(v)==0) { //if v is not already explored
            //yes, v should be explored
            //Configure the currentPath accordingly
            currentPath.push_back(v);

            //call DFS
            findAllSTPathsCore(graph, v, allPaths, currentPath, explored);

            //DesConfigure the currentPath accordingly
            currentPath.pop_back();
          }
        }
      }

      //already explored all paths with this prefix
      explored.erase(s);
    }


    //finds all paths in a graph
    vector< vector<Vertex> > findAllSTPaths(const FilteredGraph &graph) {
      vector< vector<Vertex> > allPaths; //will contain all paths in the graph

      //for each source s, call FindAllSTPaths(graph, s)
      FilteredGraph::vertex_iterator vi, vi_end;
      for (boost::tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
        auto s = *vi;

        if (in_degree(s, graph)==0) {
          //find all paths starting with s
          vector<Vertex> currentPath; //current path being built
          currentPath.push_back(s);
          set<Vertex> explored; //will keep track of the nodes that were already explored on currentPath
          findAllSTPathsCore(graph, s, allPaths, currentPath, explored);
        }
      }

      return allPaths;
    }


    //cluster similar sequences
    vector< vector < int > >  clusterSequences(const vector <string> &allSTPathsAsSequences, double snpRate) {
      //init clusters
      vector< vector < int > > clusters;
      map<int, int> stringId2ClusterId; //define the dict stringId2ClusterId

      for (int i=0; i<allSTPathsAsSequences.size(); i++) {
        vector< int> clusterDeS;
        clusterDeS.push_back(i);
        clusters.push_back(clusterDeS);

        stringId2ClusterId[i]=i;
      }



      //make a union-find algorithm where the comparison operator is the edit distance
      for (int i=0; i<allSTPathsAsSequences.size()-1; i++) {
        for (int j=i+1; j<allSTPathsAsSequences.size(); j++) {
          const string &s1 = allSTPathsAsSequences[i];
          const string &s2 = allSTPathsAsSequences[j];
          int editDistance = computeEditDistance(s1, s2);
          int lengthDifference = abs( ((int)s1.size()) - ((int)s2.size()));

          //check if edit distance is good (if both sequences are similar enough)
          if ( lengthDifference <= 2 && //we allow at most 2 indels of one base each, because 3 could be an alternative acceptor or donor
               editDistance <= 1 + snpRate * max(s1.size(), s2.size()) + lengthDifference) {
            int lowerIndex=min(stringId2ClusterId[i], stringId2ClusterId[j]);
            int upperIndex=max(stringId2ClusterId[i], stringId2ClusterId[j]);


            clusters[lowerIndex].insert(clusters[lowerIndex].end(), clusters[upperIndex].begin(), clusters[upperIndex].end());
            clusters[upperIndex].clear();
          }
        }
      }


      //clean clusters
      vector< vector < int > > cleanedClusters;
      for (const vector < int > &cluster : clusters) {
        if (cluster.empty() == false)
          cleanedClusters.push_back(cluster);
      }

      return cleanedClusters;
    }

    //gets a path and the graph as input and returns the average of the coverage of the edges in the path
    //Note: the edges are not all the k+1-mers in the path - the edges are the edges in the unitig graph (or the k+1-mers between the unitigs)
    //TODO: the most covered path is defined by the edge counts, not k-mer or unitig count
    //TODO: we did it like this because we just need to get a representative sequence and what interests us later is the edge coverage, not kmer/unitig coverage
    double getCoverage(const vector<Vertex> &path, const graph_t &graph) {
      double sum=0;
      int nbOfEdges=0;

      for (int i=0; i<path.size()-1; i++) {
        Vertex edgeSource = path[i];
        Vertex edgeTarget = path[i+1];
        Edge pathEdge = edge(edgeSource, edgeTarget, graph).first;
        sum += graph[pathEdge].totalCount;
        nbOfEdges++;
      }

      return sum/nbOfEdges;
    }


    //returns, for each cluster, the id of the most covered path
    //We chose to keep only the heaviest path in the cluster to represent all the sequences in the cluster (to remove the SNPs and indels)
    //This will simplify the graph and it will give more power to the statistical analysis
    //and will also allow us to focus on the splice sites, and remove the small variations
    vector <int> getMostCoveredPathFromClusters(const vector< vector < int > > &clusters, const vector< vector<Vertex> > &allSTPaths, const graph_t &graph) {
      vector<int> clusterMostCoveredPath;

      for (const vector < int > &cluster : clusters) {
        //we assume that the first path is the most covered one
        int indexMostCoveredPath=cluster[0];
        double coverageMostCoveredPath = getCoverage( allSTPaths[cluster[0]], graph );

        //we search the most covered path
        for (int j=1; j<cluster.size(); j++) {
          double coverageOfPathJ = getCoverage(allSTPaths[cluster[j]], graph);
          if (coverageOfPathJ > coverageMostCoveredPath) {
            coverageMostCoveredPath = coverageOfPathJ;
            indexMostCoveredPath = cluster[j];
          }
        }

        //we add the index of the most covered path to clusterMostCoveredPath
        clusterMostCoveredPath.push_back(indexMostCoveredPath);
      }

      return clusterMostCoveredPath;
    }

    //Output the bubble found
    //Complexity O(BFS) = O(n lg n) - because we use set and they are ordered data structures
    void build_gene_components::outputComplexBubble (const string &prefix, const Vertex &s, const Vertex &t, int complexBubbleIndex, int componentId, const graph_t &graph,
                              SubgraphVertexFilterAddSet* subgraphVertexFilter, FilteredGraph* filteredGraph, int nbOfCondRepl, int k, vector<string> &filesToDiffAnalysis) {
      //get the nodes of the graph that are reachable from s
      set<Vertex> verticesReachableFromS;
      {
        bfs_get_visited_nodes_visitor bfsVisitor(verticesReachableFromS);
        boost::breadth_first_search(graph, s, boost::visitor(bfsVisitor));
      }

      //build a graph only with the nodes reachable from s
      subgraphVertexFilter->removeAllVertices();
      for (const Vertex &u : verticesReachableFromS)
        subgraphVertexFilter->add(u);

      //get the nodes that are now reachable from t in the reverse BFS
      boost::reverse_graph<FilteredGraph> reversedGraph = boost::make_reverse_graph(*filteredGraph);
      set<Vertex> verticesInTheComplexBubble; //these are the vertices that are reachable from s and that can reach t
      {
        bfs_get_visited_nodes_visitor bfsVisitor(verticesInTheComplexBubble);
        boost::breadth_first_search(reversedGraph, t, boost::visitor(bfsVisitor));
      }


      //build the graph only with the nodes that are reachable from s and that can reach t
      subgraphVertexFilter->removeAllVertices();
      for (const Vertex &u : verticesInTheComplexBubble)
        subgraphVertexFilter->add(u);

      //output the graph
      printComponent(prefix+string(".complex_bubble_")+to_string(complexBubbleIndex)+".comp_uncompressed_"+to_string(componentId), graph, verticesInTheComplexBubble, nbOfCondRepl, filesToDiffAnalysis);

      /*
       * compressing bubble code - which was moved to compress the components (graph)

      //get all the paths between s and t
      vector< vector<Vertex> > allSTPaths = findAllSTPaths(*filteredGraph, s, t);

      //transform the paths into sequences
      vector <string> allSTPathsAsSequences;
      for (const vector<Vertex> &stPath : allSTPaths)
        allSTPathsAsSequences.push_back(buildSequenceUsingDBG(graph, stPath, k));

      //TODO: this is not needed, we should remove it!!!
      {
        ofstream complexBubbleSequenceFile;
        openFileForWriting(prefix+string(".complex_bubble_")+to_string(complexBubbleIndex)+string(".sequences.fa"), complexBubbleSequenceFile);
        for (int i=0; i<allSTPathsAsSequences.size(); i++)
          complexBubbleSequenceFile << ">seq" << i << endl << allSTPathsAsSequences[i] << endl;
        complexBubbleSequenceFile.close();
      }
      //TODO: this is not needed, we should remove it!!!

      //cluster the sequences
      //TODO: put 0.01 as a parameter
      vector< vector < int > >  clusters = clusterSequences(allSTPathsAsSequences, 0.01);


      //TODO: this is not needed, we should remove it!!!
      {
        ofstream clusterFile;
        openFileForWriting(prefix+string(".complex_bubble_")+to_string(complexBubbleIndex)+string(".sequences.clusters.fa"), clusterFile);
        int i=0;
        for (const vector < int > &cluster : clusters) {
          clusterFile << "Cluster " << i << ":" << endl;
          for (int idString : cluster) {
            clusterFile << allSTPathsAsSequences[idString] << endl;
          }
          clusterFile << endl << endl;
          i++;
        }
        clusterFile.close();
      }
      //TODO: this is not needed, we should remove it!!!

      //get the heaviest sequence as the representative of the cluster
      vector <int> clusterMostCoveredPath = getMostCoveredPathFromClusters(clusters, allSTPaths, graph);

      //output the heaviest sequences to a file
      string heaviestSequenceFilename = prefix+string(".complex_bubble_")+to_string(complexBubbleIndex)+string(".sequences.heaviest.fa");
      {
        ofstream heaviestSequencesFile;
        openFileForWriting(heaviestSequenceFilename, heaviestSequencesFile);
        for (int i=0; i<clusterMostCoveredPath.size(); i++) {
          heaviestSequencesFile << ">cluster_" << i << "_heaviest_path" << endl;
          heaviestSequencesFile << allSTPathsAsSequences[clusterMostCoveredPath[i]] << endl;
        }
        heaviestSequencesFile.close();
      }


      //build the unitig graph from the heaviest sequences using GATB
      gatb::core::debruijn::impl::Graph graphFromHeaviestSequences = gatb::core::debruijn::impl::Graph::create(
          "-in %s -kmer-size %d -abundance-min 0 -out %s -nb-cores %d -verbose 0",
          heaviestSequenceFilename.c_str(), k, heaviestSequenceFilename.c_str(), 1); //TODO: increase nb-cores?
      buildUnitigs(graphFromHeaviestSequences, heaviestSequenceFilename);

      //remove some unneeded files
      boost::filesystem::remove(heaviestSequenceFilename+".unitigs");
      boost::filesystem::remove(heaviestSequenceFilename+".h5");
      boost::filesystem::remove(heaviestSequenceFilename+".unitigs.features");
      boost::filesystem::remove(heaviestSequenceFilename+".unitigs.size");

      //TODO: requantify the edges of the new graph using the previous edge quantification
        //TODO: to quantify one k+1-mer (edge) of the new graph, we go through each k+1-mer (edge) of the previous graph and check if it is similar enough (with Camille's function)
        //TODO: if it is, we add the quantification to the new edge
      */
    }

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
    void build_gene_components::execute ()
    {
      string step1Folder=getInput()->getStr(STR_PAR_STEP_1_FOLDER);
      string step2Folder=getInput()->getStr(STR_PAR_STEP_2_FOLDER);
      string step3Folder=getInput()->getStr(STR_PAR_STEP_3_FOLDER);
      string prefix = getInput()->getStr(STR_PAR_GRAPH_PREFIX);
      string step1Prefix = step1Folder+string("/")+prefix;
      string step2Prefix = step2Folder+string("/")+prefix;
      string step3Prefix = step3Folder+string("/")+prefix;
      int k = getInput()->getInt(STR_PAR_K);
      int nbCores = getInput()->getInt(STR_PAR_NB_CORES);

      //read condRepls
      {
        std::ifstream condRepl2ReadsIF;
        openFileForReading(step1Prefix + ".cond_repl", condRepl2ReadsIF);
        string condRepl;
        while (condRepl2ReadsIF >> condRepl)
          condRepls.push_back(condRepl);
        condRepl2ReadsIF.close();
      }

      //get the nodes that are identified as repeats
      cerr << "Getting nodes identified as repeats..." << endl;
      set<int> repeatNodes;
      {
        string repeatNodesFilename(step2Prefix + ".unitigs.fasta_mapped.unitigs_containing_repeats_or_unmapped");
        ifstream repeatNodesFileReader;
        openFileForReading(repeatNodesFilename, repeatNodesFileReader);

        int nodeId;
        while (repeatNodesFileReader >> nodeId)
          repeatNodes.insert(nodeId);
        repeatNodesFileReader.close();
      }
      cerr << "Getting nodes identified as repeats... - Done!" << endl;

      //create the nodes of the graph
      int nbOfNodes = (countNbNodes(step1Prefix + ".unitigs.size")-repeatNodes.size()) * 2; //*2 because of the RC
      map<pair<int, char>, int> labelToIndex; //label is an entry in the .nodes file and index is the index of that node in the Boost Graph. E.g.: <14, 'R'> and 500
      map<int, pair<int, char> > indexToLabel;
      graph_t graph(nbOfNodes);
      cerr << "Creating the repeat-free Boost graph..." << endl;
      {
        string nodesFilename(step1Prefix + ".nodes");
        ifstream nodesFileReader(nodesFilename);
        int id;
        string seq;
        int index = 0;
        while (nodesFileReader >> id >> seq) {
          if (repeatNodes.count(id)>0) //I should not add this node - this is a repeat node
            continue;

          //add the FW node
          Vertex vF = vertex(index, graph);
          labelToIndex[make_pair(id, 'F')] = index;
          indexToLabel[index] = make_pair(id, 'F');

          //TODO: sequence should be stored as a 2-bit array
          graph[vF].id = id;
          graph[vF].name = seq;
          graph[vF].strand = 'F';
          graph[vF].weight = seq.length() - k + 1;
          index++;

          //add the RC node
          Vertex vR = vertex(index, graph);
          labelToIndex[make_pair(id, 'R')] = index;
          indexToLabel[index] = make_pair(id, 'R');

          //TODO: sequence should be stored as a 2-bit array
          graph[vR].id = id;
          graph[vR].name = reverse_complement(seq);
          graph[vR].strand = 'R';
          graph[vR].weight = seq.length() - k + 1;
          index++;
        }
        nodesFileReader.close();
      }
      cerr << "Creating the repeat-free Boost graph... - Done!" << endl;

      //create the edges of the graph
      cerr << "Creating the edges of the Boost graph..." << endl;
      int nbOfCondRepl=0;
      {
        string edgesFilename(step1Prefix + ".edges.dbg");
        ifstream edgesFileReader;
        openFileForReading(edgesFilename, edgesFileReader);
        string line;
        int index = 0;

        while (std::getline(edgesFileReader, line)) {
          //declare the local vars
          int from, to;
          string label;
          int count;
          int totalCount=0;
          vector<int> countVector;

          //read all the vars
          stringstream ss;
          ss << line;
          ss >> from >> to >> label;
          while (ss >> count) {
            countVector.push_back(count);
            totalCount+=count;
          }
          nbOfCondRepl=countVector.size();

          //check if the source or the target is a repeat node - if yes, we should not add this edge
          if (repeatNodes.count(from)>0 || repeatNodes.count(to)>0)
            continue;

          //add the edge
          pair<Edge, bool> return_from_add_edge = boost::add_edge(vertex(labelToIndex[make_pair(from, label[0])], graph),
                                                                  vertex(labelToIndex[make_pair(to, label[1])], graph),
                                                                  graph);
          if (return_from_add_edge.second) {
            //configure the edge
            graph[return_from_add_edge.first].id = index;
            graph[return_from_add_edge.first].from = from;
            graph[return_from_add_edge.first].to = to;
            graph[return_from_add_edge.first].label = label;
            graph[return_from_add_edge.first].countVector = countVector;
            graph[return_from_add_edge.first].totalCount = totalCount;
            graph[return_from_add_edge.first].targetWeight = graph[target(return_from_add_edge.first,
                                                                          graph)].weight;
            graph[return_from_add_edge.first].sourceWeight = graph[source(return_from_add_edge.first,
                                                                          graph)].weight;

            index++;
          }
        }
        edgesFileReader.close();
      }
      cerr << "Creating the edges of the Boost graph... - Done!" << endl;

      //these are the args to the diff analysis script
      //it is easier to do it here
      vector<string> filesToDiffAnalysis;


      //print this graph, just for viewing
      cerr << "Printing the repeat_free graph..." << endl;
      set<Vertex> allNodes;
      for (auto vp = vertices(graph); vp.first != vp.second; ++vp.first) {
        Vertex v = *vp.first;
        allNodes.insert(v);
      }
      printComponent(step3Prefix+string(".repeat_free"), graph, allNodes, nbOfCondRepl, filesToDiffAnalysis);
      cerr << "Printing the repeat_free graph... - Done!" << endl;

      //get the components of this graph
      //to do this, we have to use the undirected version of the graph
      //builds an undirected vesion of the graph
      cerr << "Getting components..." << endl;
      vector<int> componentOfThisNode;
      vector<set<Vertex>> nodesInComponent;
      int num;
      {
        //make a copy of the graph
        graph_t udCopyOfGraph;
        boost::copy_graph(graph, udCopyOfGraph);

        //goes through all the edges and add the reverse edge to the graph, creating an undirected graph
        graph_t::edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei) {
          auto uvEdge = *ei;
          add_edge(target(uvEdge, graph), source(uvEdge, graph), udCopyOfGraph);
        }

        //get the components on the udCopyOfGraph
        componentOfThisNode = vector<int>(num_vertices(udCopyOfGraph));
        num = connected_components(udCopyOfGraph, &componentOfThisNode[0]);
        nodesInComponent = vector<set<Vertex>>(num);
        for (int i = 0; i != componentOfThisNode.size(); ++i)
          nodesInComponent[componentOfThisNode[i]].insert(vertex(i, graph));

      }
      cerr << "Getting components... - Done!" << endl;

      //output the components on the graph now
      cerr << "Printing all components..." << endl;
      for (int i = 0; i < num; i++)
        printComponent(step3Prefix+string(".comp_uncompressed_")+to_string(i), graph, nodesInComponent[i], nbOfCondRepl, filesToDiffAnalysis);
      cerr << "Printing all components... - Done!" << endl;


      /*
      cerr << "Removing small errors from the graph..." << endl;
      list<graph_t> compressedGraphs;
      for (int i = 0; i < num; i++) {
        //check if this component has at least 4 nodes
        if (nodesInComponent[i].size() < 4) {
          cout << "Component " << i << " has less than 4 nodes and will not be considered for further analysis" << endl;
          continue;
        }

        //create a filtered graph that will represent this component
        SubgraphEdgeFilterAddSet* subgraphEdgeFilter = new SubgraphEdgeFilterAddSet(&graph);
        SubgraphVertexFilterAddSet* subgraphVertexFilter = new SubgraphVertexFilterAddSet(&graph, subgraphEdgeFilter);
        FilteredGraph* componentGraph = new FilteredGraph(graph,
                                                         SubgraphEdgeFilterForFG<SubgraphEdgeFilterAddSet>(subgraphEdgeFilter),
                                                         SubgraphVertexFilterForFG<SubgraphVertexFilterAddSet>(subgraphVertexFilter));
        for (const auto &node : nodesInComponent[i])
          subgraphVertexFilter->add(node);

        //main step of removing small errors
        //TODO: do this iteratively until convergence?
        cerr << "Compressing graph " << i << "..." << endl;
        graph_t compressedGraph = compressGraph(*componentGraph, k);
        cerr << "Compressing graph " << i << " - Done!" << endl;


        //debug code - checking level 1 compression
        {
          auto verticesIterator = vertices(compressedGraph);
          set<Vertex> nodesInThisComponent(verticesIterator.first, verticesIterator.second); //TODO: should probably fix up printComponent() function
          printComponent(step3Prefix+string(".comp_level_1_compression_")+to_string(i), compressedGraph, nodesInThisComponent, nbOfCondRepl, filesToDiffAnalysis);
        }



        //parition the nodes with the same single source and target
        map<tuple<Vertex,Vertex>, set<Vertex>> pairSourceTarget2Nodes;

        //TODO: change the way we iterate, this is really verbose
        auto beginAndEndVI = vertices(compressedGraph);
        typedef decltype(beginAndEndVI.first) VertexIteratorType;
        VertexIteratorType beginVI, endVI;
        //goes through all nodes
        for (tie(beginVI, endVI)=beginAndEndVI; beginVI != endVI; ++beginVI) {
          Vertex node = *beginVI;

          //check if the in_degree and out_degree of the node is 1
          if (in_degree(node, compressedGraph)==1 && out_degree(node, compressedGraph)==1) {
            //get the source and target vertices
            Edge outEdge = *(out_edges(node, compressedGraph).first);
            Vertex target = boost::target(outEdge, compressedGraph);
            Edge inEdge = *(in_edges(node, compressedGraph).first);
            Vertex source = boost::source(inEdge, compressedGraph);

            //add it to the partitioned list
            pairSourceTarget2Nodes[make_tuple(source, target)].insert(node);
          }
        }

        //for each partition having >= 2 nodes, we check if the difference between the nodes is <= 1 nt
        for (const auto& mapPair : pairSourceTarget2Nodes) {
          Vertex source = std::get<0>(mapPair.first);
          Vertex target = std::get<1>(mapPair.first);
          const set<Vertex> &nodes = mapPair.second;
          if (nodes.size()>=2) {
            //check if we should compress or not
            bool weShouldCompress=true;
            auto nodeIt = nodes.begin();
            const string &firstNodeSeq = compressedGraph[*nodeIt].name;
            for(++nodeIt; nodeIt != nodes.end(); ++nodeIt) {
              if (computeEditDistance(firstNodeSeq, compressedGraph[*nodeIt].name) > 1) {
                weShouldCompress=false;
                break;
              }
            }

            if (weShouldCompress) {
              //compress the graph
              //in this case, we just choose one node to be kept, increase its coverage with the other nodes counts, and remove the other nodes

              //TODO: the node to be kept sould be the most abundant one. We don't do this for now since we don't have the node abundances, but this must be changed
              //TODO: for now, we choose the node with the highest edge abundances, should be enough, right?
              Vertex nodeToBeKept = *(nodes.begin());
              int maxEdgeAbundance = -1;

              for (const auto& node : nodes) {
                //computes the edge abundance of node
                Edge outEdge = *(out_edges(node, compressedGraph).first);
                Edge inEdge = *(in_edges(node, compressedGraph).first);
                int edgeAbundance = compressedGraph[outEdge].totalCount + compressedGraph[inEdge].totalCount;
                if (edgeAbundance > maxEdgeAbundance) {
                  //update
                  maxEdgeAbundance = edgeAbundance;
                  nodeToBeKept = node;
                }
              }

              //remove the other nodes
              for (const auto& node : nodes) {
                if (node!=nodeToBeKept) {
                  clear_vertex(node, compressedGraph);
                  remove_vertex(node, compressedGraph);
                }
              }
            }
          }
        }

        //recompress the graph
        compressedGraphs.push_back(compressGraph(compressedGraph, k));


        /*
        //get all paths from sources to targets in this graph
        //TODO: if the graph has no sources (i.e. a cycle, we won't find any path)
        cout << "# of nodes: " << nodesInComponent[i].size() << endl;
        vector< vector<Vertex> > allSTPaths = findAllSTPaths(*componentGraph);
        cout << "# of paths: " << allSTPaths.size() << endl;

        //transform the paths into sequences
        vector <string> allSTPathsAsSequences;
        for (const vector<Vertex> &stPath : allSTPaths)
          allSTPathsAsSequences.push_back(buildSequenceUsingDBG(graph, stPath, k));


        //TODO: this is not needed, we should remove it!!!
        {
          ofstream compUncompressedSequenceFile;
          openFileForWriting(step3Prefix+string(".comp_uncompressed_")+to_string(i)+string(".sequences.fa"), compUncompressedSequenceFile);
          for (int i=0; i<allSTPathsAsSequences.size(); i++)
            compUncompressedSequenceFile << ">seq" << i << endl << allSTPathsAsSequences[i] << endl;
          compUncompressedSequenceFile.close();
        }
        //TODO: this is not needed, we should remove it!!!


        //cluster the sequences
        //TODO: put 0.01 as a parameter
        vector< vector < int > >  clusters = clusterSequences(allSTPathsAsSequences, 0.01);

        //some checks for continuity after getting the clusters
        if (clusters.size()==0) {
          stringstream ss;
          ss << "I did not find any sequence for component " << i << "... is this an error?";
          fatalError(ss.str());
        }
        if (clusters.size()==1) {
          cout << "Component " << i << " has only one cluster of sequences and will not be considered for further analysis (probably a component representing only SNPs/indels)" << endl;
          continue;
        }


        //TODO: this is not needed, we should remove it!!!
        {
          ofstream clusterFile;
          openFileForWriting(step3Prefix+string(".comp_uncompressed_")+to_string(i)+string(".sequences.clusters.fa"), clusterFile);
          int i=0;
          for (const vector < int > &cluster : clusters) {
            clusterFile << "Cluster " << i << ":" << endl;
            for (int idString : cluster) {
              clusterFile << allSTPathsAsSequences[idString] << endl;
            }
            clusterFile << endl << endl;
            i++;
          }
          clusterFile.close();
        }
        //TODO: this is not needed, we should remove it!!!

        //get the heaviest sequence as the representative of the cluster
        vector <int> clusterMostCoveredPath = getMostCoveredPathFromClusters(clusters, allSTPaths, graph);

        //output the heaviest sequences to a file
        string heaviestSequenceFilename = step3Prefix+string(".comp_uncompressed_")+to_string(i)+string(".sequences.heaviest.fa");
        {
          ofstream heaviestSequencesFile;
          openFileForWriting(heaviestSequenceFilename, heaviestSequencesFile);
          for (int i=0; i<clusterMostCoveredPath.size(); i++) {
            heaviestSequencesFile << ">cluster_" << i << "_heaviest_path" << endl;
            heaviestSequencesFile << allSTPathsAsSequences[clusterMostCoveredPath[i]] << endl;
          }
          heaviestSequencesFile.close();
        }


        //build the unitig graph from the heaviest sequences using GATB
        gatb::core::debruijn::impl::Graph graphFromHeaviestSequences = gatb::core::debruijn::impl::Graph::create(
            "-in %s -kmer-size %d -abundance-min 0 -out %s -nb-cores %d -verbose 0",
            heaviestSequenceFilename.c_str(), k, heaviestSequenceFilename.c_str(), 1); //TODO: increase nb-cores?
        buildUnitigs(graphFromHeaviestSequences, heaviestSequenceFilename);

        //remove some unneeded files
        boost::filesystem::remove(heaviestSequenceFilename+".h5");
        boost::filesystem::remove(heaviestSequenceFilename+".unitigs.features");
        boost::filesystem::remove(heaviestSequenceFilename+".unitigs.size");

        //TODO:
        //requantify the edges of the new graph using the previous edge quantification
        //to quantify one k+1-mer (edge) of the new graph, we go through each k+1-mer (edge) of the previous graph and check if it is similar enough, i.e. if the edit distance between the two edges is < snpRate
        //if it is, we add the quantification to the new edge


      }
      cerr << "Removing small errors from the graph - Done!" << endl;

      //print the compressed components
      //output the components on the graph now
      cerr << "Printing all compressed components..." << endl;
      int index=0;
      for (const auto &compressedComponent : compressedGraphs) {
        auto verticesIterator = vertices(compressedComponent);
        set<Vertex> nodesInThisComponent(verticesIterator.first, verticesIterator.second); //TODO: should probably fix up printComponent() function
        printComponent(step3Prefix+string(".comp_compressed_")+to_string(index++), compressedComponent, nodesInThisComponent, nbOfCondRepl, filesToDiffAnalysis);
      }

      cerr << "Printing all compressed components... - Done!" << endl;
      */

      cerr << "Getting the complex bubbles..." << endl;
      {
        //create a filtered graph that will represent all the flow graphs rooted at each vertex s
        SubgraphEdgeFilterAddSet* subgraphEdgeFilter = new SubgraphEdgeFilterAddSet(&graph);
        SubgraphVertexFilterAddSet* subgraphVertexFilter = new SubgraphVertexFilterAddSet(&graph, subgraphEdgeFilter);
        FilteredGraph* filteredGraph = new FilteredGraph(graph,
                                                         SubgraphEdgeFilterForFG<SubgraphEdgeFilterAddSet>(subgraphEdgeFilter),
                                                         SubgraphVertexFilterForFG<SubgraphVertexFilterAddSet>(subgraphVertexFilter));


        //For each node v of the graph, compute the 1st-layer dominators of v
        //and if a node s is such that immDom(s)==v, then s is an endpoint of a (v,s)-multiwise bubble
        cerr << "Finding endpoints of the complex bubbles..." << endl;
        map<Vertex, set<Vertex> > nodeToComplexBubbleEndpoints;
        for (auto vp = vertices(graph); vp.first != vp.second; ++vp.first) {
          Vertex v = *vp.first;
          if (out_degree(v, graph) >= 2) { //pruning
            //get the flowgraph of  of v
            {
              //get the vertices in the flowgraph of v
              set<Vertex> verticesInTheFlowgraphOfV;
              bfs_get_visited_nodes_visitor bfsVisitor(verticesInTheFlowgraphOfV);
              boost::breadth_first_search(graph, v, boost::visitor(bfsVisitor));

              //build the flowgraph
              subgraphVertexFilter->removeAllVertices();
              for (const Vertex &u : verticesInTheFlowgraphOfV)
                subgraphVertexFilter->add(u);
            }


            //find all endpoints of bubbles starting in v
            {
              //get the type of vertex_index (probably size_t)
              typedef boost::property_map<graph_t, boost::vertex_index_t>::type VerticeIndexType;

              //this will map the index of the vertices to a Vertex
              //it will represent, for a vertex v, <Immediate dom of v (Vertex), index of v (VertexIndexType)>
              typedef boost::iterator_property_map<vector<Vertex>::iterator, VerticeIndexType> PredMap;

              //let u be a vertex, and imm. dom (u) = v, then domTreePredVector[u]=v
              vector<Vertex> domTreePredVector(num_vertices(*filteredGraph), -1);
              PredMap domTreePredMap =
                  make_iterator_property_map(domTreePredVector.begin(), boost::get(boost::vertex_index, *filteredGraph));

              // call the Lengauer-Tarjan dominator tree algorithm
              lengauer_tarjan_dominator_tree(*filteredGraph, v, domTreePredMap);



              //debug
              /*
              cerr << "Computing dominator tree with root = " << nodeIdToStr(*filteredGraph, v) << endl;
              for (int u = 0; u < domTreePredVector.size(); u++) {
                if (domTreePredVector[u]!=-1)
                  cerr << "imm_dom[" << nodeIdToStr(*filteredGraph, u) << "] = " << nodeIdToStr(*filteredGraph, domTreePredVector[u]) << endl;
              }
               */

              //find which nodes are the endpoints of a complex bubble starting with v
              //if a node u is such that immDom(u)==v, then u is an endpoint of a (v,u)-complex bubble
              for (int u = 0; u < domTreePredVector.size(); u++) {
                //check if the immediate dominator of u is v
                if (domTreePredVector[u] == v) { //yeah

                  //the following conditions treat the outneighbours of the source:
                  //the immediate dominator of an outneighbour of the source will always be the source, so this is a special case
                  //for now, i just forbid an outneighbour of the source to be an endpoint of a (v,u)-complex bubble
                  if (edge(v,u,*filteredGraph).second==false) //if there is no edge v->u (i.e., if u is not an out neighbour of v)
                    nodeToComplexBubbleEndpoints[v].insert(u);
                }
              }
            }
          }
        }
        cerr << "Finding endpoints of the complex bubbles... - Done!" << endl;

        cerr << "Removing non-maximal complex bubbles..." << endl;

        cerr << "[DEBUG] MULTI-WISE BUBBLES FOUND: " << endl;
        for (const auto &pairS_setOfTs : nodeToComplexBubbleEndpoints) {
          Vertex s = pairS_setOfTs.first;
          for (Vertex t : pairS_setOfTs.second)
            cerr << "(" << nodeIdToStr(graph, s) << ", " << nodeIdToStr(graph, t) << ")" << endl;
        }
        cerr << endl << endl;

        //let b=(s,t) and b'=(s',t') be two complex bubbles
        //b' is non-maximal if b' is contained in b
        //this can be verified if there exists:
        //1) a path from s to s';
        //2) a path from t' to t;
        //3) no path from s' to s;
        //4) no path from t to t';

        //compute the transitive closure graph
        TC_graph_t TC;
        transitive_closure(graph, TC);

        //get the bubbles that are not maximal
        set<pair<Vertex, Vertex>> endpointsOfNonMaximalComplexBubbles;
        for (const auto &pairS_setOfTs : nodeToComplexBubbleEndpoints) {
          Vertex s = pairS_setOfTs.first;
          for (Vertex t : pairS_setOfTs.second) {
            //here, we check the (s,t)-complex bubble

            //we get sPrime and tPrime
            for (const auto &pairSPrime_setOfTPrimes : nodeToComplexBubbleEndpoints) {
              Vertex sPrime = pairSPrime_setOfTPrimes.first;
              for (Vertex tPrime : pairSPrime_setOfTPrimes.second) {
                if (s!=sPrime || t!=tPrime) { //at least one of these must differ
                  /*
                  cerr << "Checking if ComplexBubble(" << nodeIdToStr(graph, s) << "," << nodeIdToStr(graph, t) << ") contains ComplexBubble("
                  << nodeIdToStr(graph, sPrime) << "," << nodeIdToStr(graph, tPrime) <<")..." << endl;
                   */

                  //check if there is:
                  if ((s==sPrime || edge(s, sPrime, TC).second==true) && //1) a path from s to s';
                      (t==tPrime || edge(tPrime, t, TC).second==true) && //2) a path from t' to t;
                      (s==sPrime || edge(sPrime, s, TC).second==false) && //3) no path from s' to s;
                      (t==tPrime || edge(t, tPrime, TC).second==false )) //4) no path from t to t';
                  {
                    //yes, this ComplexBubble (sPrime, tPrime) is not maximal
                    cerr << "[DEBUG] ComplexBubble(" << nodeIdToStr(graph, s) << "," << nodeIdToStr(graph, t) << ") contains MultiwiseBubble("
                    << nodeIdToStr(graph, sPrime) << "," << nodeIdToStr(graph, tPrime) <<")" << endl;
                    endpointsOfNonMaximalComplexBubbles.insert(make_pair(sPrime, tPrime));
                  }
                }
              }
            }
          }
        }

        //get the bubbles that are maximal now
        cerr << "[DEBUG] COMPLEX BUBBLES FOUND: " << endl;
        set<pair<Vertex, Vertex>> endpointsOfMaximalComplexBubbles;
        for (const auto &pairS_setOfTs : nodeToComplexBubbleEndpoints) {
          Vertex s = pairS_setOfTs.first;
          for (Vertex t : pairS_setOfTs.second) {
            //here, we check if the (s,t)-complex bubble is not in endpointsOfNonMaximalComplexBubbles
            pair<Vertex, Vertex> stPair = make_pair(s, t);
            if (endpointsOfNonMaximalComplexBubbles.count(stPair) == 0) {
              cerr << "(" << nodeIdToStr(graph, s) << "," << nodeIdToStr(graph, t) << ")" << endl;
              //it is not, so it is a maximal complex bubble
              endpointsOfMaximalComplexBubbles.insert(stPair);
            }
          }
        }
        cerr << "Removing non-maximal complex bubbles - Done!" << endl;

        cerr << "Outputting the complex bubbles..." << endl;
        int complexBubbleIndex=0;
        for (const auto &stPair : endpointsOfMaximalComplexBubbles) {
          outputComplexBubble(step3Prefix, stPair.first, stPair.second, complexBubbleIndex, componentOfThisNode[stPair.first], graph, subgraphVertexFilter, filteredGraph, nbOfCondRepl, k, filesToDiffAnalysis);
          complexBubbleIndex++;
        }
        cerr << "Outputting the complex bubbles - Done!" << endl;

        cerr << "Output " << complexBubbleIndex << " complex bubbles!" << endl;

        delete subgraphEdgeFilter;
        delete subgraphVertexFilter;
        delete filteredGraph;
      }
      cerr << "Getting the complex bubbles - Done!" << endl;


      //print filesToDiffAnalysis
      ofstream filesToDiffAnalysisFile;
      openFileForWriting(step3Prefix+string(".diff_analysis_files"), filesToDiffAnalysisFile);
      for (const auto &file : filesToDiffAnalysis)
        filesToDiffAnalysisFile << file << endl;
      filesToDiffAnalysisFile.close();
    }


}
