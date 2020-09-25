//
// Created by Leandro Ishi Soares de Lima on 16/03/18.
//

#ifndef DBGSPLITTER_BOOSTGRAPHDEFS_H
#define DBGSPLITTER_BOOSTGRAPHDEFS_H
#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/properties.hpp>
#include <string>
#include <vector>

using namespace std;


namespace build_gene_components {
    //define the structures needed by these classes
    struct VertexInfo {
        int id; //do not use unsigned values //TODO: I DON'T THINK WE REALLY NEED THIS - WE CAN USE BOOST TO RETRIEVE THE ID
        string name; //probably represent this as (unitig id, pos)
        char strand; //strand
        int weight; //how many kmers we have here
    };

    //edge informations
    struct EdgeInfo {
        int id; //do not use unsigned values //TODO: I DON'T THINK WE REALLY NEED THIS - WE CAN USE BOOST TO RETRIEVE THE ID
        int from; //TODO: I DON'T THINK WE REALLY NEED THIS - WE CAN USE BOOST TO RETRIEVE THE SOURCE NODE
        int to; //TODO: I DON'T THINK WE REALLY NEED THIS - WE CAN USE BOOST TO RETRIEVE THE TARGET NODE
        string label;
        vector<int> countVector;
        int totalCount;
        int sourceWeight;
        int targetWeight;
    };

    //some typedefs for making life easier
    typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::bidirectionalS, VertexInfo, EdgeInfo> graph_t;
    typedef boost::graph_traits<graph_t>::vertex_descriptor Vertex;
    typedef boost::graph_traits<graph_t>::edge_descriptor Edge;
    typedef boost::adjacency_list <> TC_graph_t;
}


#endif //DBGSPLITTER_BOOSTGRAPHDEFS_H_H
