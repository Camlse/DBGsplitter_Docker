//
// Created by Leandro Ishi Soares de Lima on 16/03/18.
//

#ifndef DBGSPLITTER_BUILD_GENE_COMPONENTS_H
#define DBGSPLITTER_BUILD_GENE_COMPONENTS_H

#include <gatb/gatb_core.hpp>
#include <string>
#include <sstream>
#include "Utils.h"
#include "global.h"
#include "BoostGraphDefs.h"
#include "BoostGraphFilters.h"
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/dominator_tree.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/transitive_closure.hpp>
#include <cstdlib>
#include "BoostGraphFilters.h"
#include <boost/graph/reverse_graph.hpp>
#include "GraphOutput.h"
#include "build_dbg.hpp"
#include <tuple>



using namespace std;

namespace build_gene_components {
    class bfs_get_visited_nodes_visitor : public boost::default_bfs_visitor {
    public:
        bfs_get_visited_nodes_visitor(set<Vertex>& vertices):vertices(vertices){ }

        template <class VertexType, class GraphType>
        void discover_vertex(VertexType u, const GraphType &g) {
          vertices.insert(u);
        }
        set<Vertex>& vertices; //vertices found so far
    };

    //given a path of vertex described by two iterators, build the sequence of this path
    template <class InputIterator, class GraphType>
    string buildSequenceUsingDBG(InputIterator begin, InputIterator end, const GraphType &graph, int k) {
      string sequence="";
      sequence = graph[*begin].name.substr(graph[*begin].name.size()-k);
      ++begin;
      --end;
      for_each(begin, end, [&](const Vertex &v) {
          sequence += graph[v].name.substr(k - 1);
      });
      sequence += string(1, graph[*end].name[k-1]);
      return sequence;
    }

    //compress linear paths into one node
    //declare a structure that will represent a compressed node
    struct CompressedNode {
        int id;
        string seq;
        Vertex oldIdOfTheRightestNode;
        list<Vertex> outNeighboursWithOldId;
    };

    //compress the graph
    template <class GraphType>
    graph_t compressGraph(const GraphType &graph, int k) {
      auto beginAndEndVI = vertices(graph);
      typedef decltype(beginAndEndVI.first) VertexIteratorType;
      VertexIteratorType beginVI, endVI;
      map<Vertex, bool> node2Marked; //will mark all nodes that were already visited
      map<Vertex, int> oldId2NewUnitigId; //will map the old vertex id in the graph to the new unitig id

      list<CompressedNode> compressedNodes;

      //goes through all nodes
      int id=0;
      for (tie(beginVI, endVI)=beginAndEndVI; beginVI != endVI; ++beginVI) {
        Vertex node = *beginVI;
        if (node2Marked[node]) //node is already marked
          continue;

        node2Marked[node] = true; //mark the node

        //path that we have to compact
        list<Vertex> pathToCompact{node};
        oldId2NewUnitigId[node] = id;

        //try to extend to the right
        Vertex currentNode = node;
        while (out_degree(currentNode, graph)==1) {
          //update current node
          Edge edge = *(out_edges(currentNode, graph).first);
          currentNode = target(edge, graph);

          //add the node to the path
          pathToCompact.push_back(currentNode);
          oldId2NewUnitigId[currentNode] = id;

          //mark the node
          node2Marked[currentNode]=true;
        }
        auto oldIdOfTheRightestNode = currentNode;

        //try to extend to the left
        currentNode = node;
        while (in_degree(currentNode, graph)==1) {
          //update current node
          Edge edge = *(in_edges(currentNode, graph).first);
          currentNode = source(edge, graph);

          //add the node to the path
          pathToCompact.push_front(currentNode);
          oldId2NewUnitigId[currentNode] = id;

          //mark the node
          node2Marked[currentNode]=true;
        }

        //create the compressed node
        CompressedNode compressedNode;
        compressedNode.id = id++;
        compressedNode.oldIdOfTheRightestNode = oldIdOfTheRightestNode;
        compressedNode.seq = buildSequenceUsingDBG(pathToCompact.begin(), pathToCompact.end(), graph, k);

        //fill in outNeighboursWithOldId
        auto outEdgesLastNode = out_edges(pathToCompact.back(), graph);
        typedef decltype(outEdgesLastNode.first) EdgeIteratorType;
        EdgeIteratorType beginEI, endEI;
        for(tie(beginEI, endEI)=outEdgesLastNode; beginEI != endEI; ++beginEI) {
          auto edge = *beginEI;
          compressedNode.outNeighboursWithOldId.push_back(target(edge, graph));
        }
      }


      //create the compressed graph
      graph_t compressedGraph(compressedNodes.size());
      for (const auto &compressedNode : compressedNodes) {
        Vertex vF = compressedNode.id;
        compressedGraph[vF].id = compressedNode.id;
        compressedGraph[vF].name = compressedNode.seq;
        compressedGraph[vF].strand = 'F'; //TODO: is this problematic???
        compressedGraph[vF].weight = compressedNode.seq.length() - k + 1;
      }

      //create the edges of the compressed graph
      id = 0;
      for (const auto &compressedNode : compressedNodes) {
        for (const auto& neighbourWithOldId : compressedNode.outNeighboursWithOldId) {
          //add the edge
          Edge edge = boost::add_edge(compressedNode.id, oldId2NewUnitigId[neighbourWithOldId], compressedGraph).first;

          Edge oldEdge = boost::edge(compressedNode.oldIdOfTheRightestNode, neighbourWithOldId, graph).first;

          compressedGraph[edge].id = id++;
          compressedGraph[edge].from = compressedNode.id;
          compressedGraph[edge].to = oldId2NewUnitigId[neighbourWithOldId];
          compressedGraph[edge].label = string("FF"); //TODO: is this problematic???
          compressedGraph[edge].countVector = graph[oldEdge].countVector;
          compressedGraph[edge].totalCount = graph[oldEdge].totalCount;
          compressedGraph[edge].targetWeight = compressedGraph[target(edge,compressedGraph)].weight;
          compressedGraph[edge].sourceWeight = compressedGraph[source(edge,compressedGraph)].weight;
        }
      }

      //we finished constructing the compressed graph, return it
      return compressedGraph;
    }

    class build_gene_components : public Tool {
    private:
        void populateParser();
        void printComponent(const string &prefix, const graph_t &graph, const set<Vertex> &nodesInComponent,
                            int nbOfCondRepl, vector<string> &filesToDiffAnalysis);
        void outputComplexBubble (const string &prefix, const Vertex &s, const Vertex &t, int complexBubbleIndex, int componentId, const graph_t &graph,
                                                         SubgraphVertexFilterAddSet* subgraphVertexFilter, FilteredGraph* filteredGraph, int nbOfCondRepl, int k, vector<string> &filesToDiffAnalysis);
        vector<string> condRepls;

    public:
        // Constructor
        build_gene_components ();

        // Actual job done by the tool is here
        void execute ();

        //overriding this in order to exit the tool when finding a problem with the arguments
        IProperties* run (int argc, char* argv[])
        {
          IProperties* toReturn = Tool::run(argc, argv);
          if (!toReturn)
            std::exit(1);
          return toReturn;
        }
    };
}



#endif //DBGSPLITTER_BUILD_GENE_COMPONENTS_H
