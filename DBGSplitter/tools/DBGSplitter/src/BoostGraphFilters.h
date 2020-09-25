//
// Created by Leandro Ishi Soares de Lima on 18/04/17.
//

#ifndef DBGSPLITTER_BOOSTGRAPHFILTERS_H
#define DBGSPLITTER_BOOSTGRAPHFILTERS_H

#include <vector>
#include <set>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/filtered_graph.hpp>
#include "BoostGraphDefs.h"

using namespace std;

namespace build_gene_components {
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    //This class uses a set<Edge> as predicate, where if the edge is in the set, then the edge IS PRESENT
    //Good to use when you add VERY FEW things from the graph, and remove A LOT of them
    class SubgraphEdgeFilterAddSet {
    private:
        graph_t *g;
        set<Edge> predicate;
    public:
        SubgraphEdgeFilterAddSet() : g(NULL), predicate() { } //default constructor, needed by boost
        SubgraphEdgeFilterAddSet(graph_t *g) : g(g), predicate() { } //at the start, the graph is empty!

        void removeAllEdges() { predicate.clear(); }

        void add(const Edge &e) {
          int id = (*g)[e].id;
          predicate.insert(e);
        }

        void remove(const Edge &e) {
          int id = (*g)[e].id;
          predicate.erase(e);
        }

        //check if the edge is present or not
        bool check(const Edge &e) const {
          int id = (*g)[e].id;
          return predicate.count(e)>0;
        }
    };

    //This class uses a set<Vertex> as predicate, where if the vertex is in the set, then the vertex IS PRESENT
    //Good to use when you add VERY FEW things from the graph, and remove A LOT of them
    class SubgraphVertexFilterAddSet {
    private:
        graph_t *g;
        set<Vertex> predicate;
        SubgraphEdgeFilterAddSet *subgraphEdgeFilter;
    public:
        SubgraphVertexFilterAddSet() : g(NULL), predicate(),
                                 subgraphEdgeFilter(NULL) { } //default constructor, needed by boost
        SubgraphVertexFilterAddSet(graph_t *g, SubgraphEdgeFilterAddSet *subgraphEdgeFilter) : g(g),
                                                                                   predicate(),
                                                                                   subgraphEdgeFilter(subgraphEdgeFilter) { } //at the start, the graph is empty!

        void removeAllVertices() {
          subgraphEdgeFilter->removeAllEdges();
          predicate.clear();
        }

        void add(const Vertex &v) {
          auto outEdgeIt = out_edges(v, *g);
          for_each(outEdgeIt.first, outEdgeIt.second, [&](const Edge &e) {
              if (this->check(target(e, *g)))
                subgraphEdgeFilter->add(e);
          });
          auto inEdgeIt = in_edges(v, *g);
          for_each(inEdgeIt.first, inEdgeIt.second, [&](const Edge &e) {
              if (this->check(source(e, *g)))
                subgraphEdgeFilter->add(e);
          });
          predicate.insert(v);
        }

        void remove(const Vertex &v) {
          auto outEdgeIt = out_edges(v, *g);
          for_each(outEdgeIt.first, outEdgeIt.second, [&](const Edge &e) {
              if (this->check(target(e, *g)))
                subgraphEdgeFilter->remove(e);
          });
          auto inEdgeIt = in_edges(v, *g);
          for_each(inEdgeIt.first, inEdgeIt.second, [&](const Edge &e) {
              if (this->check(source(e, *g)))
                subgraphEdgeFilter->remove(e);
          });
          predicate.erase(v);
        }

        //check if the vertex is present or not
        bool check(const Vertex &v) const {
          return predicate.count(v)>0;
        }
    };
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////



    //This is what has to be given to the Filtered graph
    template <class VertexFilterClass>
    class SubgraphVertexFilterForFG {
    public:
        VertexFilterClass *subgraphVertexFilter; //An object simply does not work here... It needs to be a pointer... I have no goddamn idea why
        SubgraphVertexFilterForFG():subgraphVertexFilter(NULL){}
        SubgraphVertexFilterForFG(VertexFilterClass *subgraphVertexFilter):subgraphVertexFilter(subgraphVertexFilter){}
        bool operator()(const Vertex &v) const { return subgraphVertexFilter->check(v); }
    };
    //This is what has to be given to the Filtered graph
    template <class EdgeFilterClass>
    class SubgraphEdgeFilterForFG {
    public:
        EdgeFilterClass *subgraphEdgeFilter; //An object simply does not work here... It needs to be a pointer... I have no goddamn idea why
        SubgraphEdgeFilterForFG():subgraphEdgeFilter(NULL){}
        SubgraphEdgeFilterForFG(EdgeFilterClass *subgraphEdgeFilter):subgraphEdgeFilter(subgraphEdgeFilter){}
        bool operator()(const Edge &e) const { return subgraphEdgeFilter->check(e); }
    };
    typedef boost::filtered_graph <graph_t, SubgraphEdgeFilterForFG<SubgraphEdgeFilterAddSet>, SubgraphVertexFilterForFG<SubgraphVertexFilterAddSet> > FilteredGraph;
    typedef boost::graph_traits<FilteredGraph>::adjacency_iterator FGAdjacencyIterator;
    typedef boost::graph_traits<FilteredGraph>::vertex_descriptor FilteredGraphVertex;
}

#endif //DBGSPLITTER_BOOSTGRAPHFILTERS_H
