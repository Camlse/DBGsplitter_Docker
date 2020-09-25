//
// Created by Leandro Ishi Soares de Lima on 10/10/17.
//

#ifndef KSGATB_SOLIDKMERORACLE_H
#define KSGATB_SOLIDKMERORACLE_H

#include <vector>
#include <string>
#include <sstream>
#include <gatb/gatb_core.hpp>
#include "ExceptionWithFileAndLine.h"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>

using namespace std;


//represents the unitig id a GATB node maps to
//pos is not always defined - it is optional
struct UnitigInfo {
    long unitigId;
    char strand;
    int pos;
    static UnitigInfo InvalidKmer; //represents an invalid kmer

    UnitigInfo (long unitigId, char strand, int pos=-1) : unitigId(unitigId), strand(strand), pos(pos){}
    UnitigInfo(){}
    string toString() const;
};

//exceptions in case of errors
class InexistentKmer : public ExceptionWithFileAndLine {
    using ExceptionWithFileAndLine::ExceptionWithFileAndLine;
};
class InvalidUnitig : public ExceptionWithFileAndLine {
    using ExceptionWithFileAndLine::ExceptionWithFileAndLine;
};
class InvalidKmer : public ExceptionWithFileAndLine {
    using ExceptionWithFileAndLine::ExceptionWithFileAndLine;
};

//Used to store the unitigs sequence and to check if a kmer really exist in the graph
class SolidKmerOracle {
public:
    //constructor
    //Graph: the graph to work on
    SolidKmerOracle(Graph* graph, const string &graphPrefix) : graph(graph), graphPrefix(graphPrefix)
        ,unitigsSequences(), nodeIdToUnitigId(((size_t)graph->getInfo()["kmers_nb_solid"]->getInt())) {}

    //constructor strictly for the serialization
    SolidKmerOracle(){}
    ~SolidKmerOracle(){
        delete graph;
    }

    //@brief: Add an unitig sequence that is present on the graph, by adding the sequence to unitigsSequences and also filling up nodeIdToUnitigId
    //@params:
    //  seq: the unitig sequence - this is intentionally seq and not const string &seq
    //  startingNode: the first node of the unitig sequence in order to traverse the graph
    void addUnitigSequence (string seq, const Node &leftestNode);

    //@brief: checks if the kmer given as input is in the graph or not. If it is, return an UnitigIdStrandPos of the kmer. Otherwise, throws either an invalid kmer or an inexistent kmer exception
    //this is intentionally string kmer and not const string &kmer
    UnitigInfo contains (string kmer) const;

    //get the unitig sequence stored in this oracle
    string getUnitigSequence (int index) const{
        if (index>=unitigsSequences.size())
            throw InvalidUnitig(__FILE__, __LINE__);
        return unitigsSequences[index];
    }

    bool isTrustableUnitig (int index, int minSizeUnitigs) const {
        return getUnitigSequence(index).size()>=getKmerSize()+minSizeUnitigs;
    }

    //get the nb of unitigs
    int getNbOfUnitigs () const {
        return unitigsSequences.size();
    }

    //get the kmer size
    long getKmerSize() const {
        return graph->getKmerSize();
    }
private:
    //for easy testing
    friend class EYTATester;

    //attributes
    Graph* graph; //the graph where this oracle works on
    string graphPrefix; //graph is written to graphPrefix+".h5"

    //TODO: vector<string> is not good here... when it reaches the limit, it is going to double and copy all the elements...
    //TODO: find something more efficient...
    vector<string> unitigsSequences; //the unitig sequences themselves - TODO: we should not use string, but a bit array (look at codeseed of the ModelDirect, Canonical, etc)
    //TODO: find something more efficient...

    vector < long > nodeIdToUnitigId; //maps node to unitig id, strand pos

    //serialization
    friend class boost::serialization::access;
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar  & graphPrefix;
        ar  & unitigsSequences;
        ar  & nodeIdToUnitigId;
    }
    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar  & graphPrefix;
        ar  & unitigsSequences;
        ar  & nodeIdToUnitigId;

        //when loading the object, load also the graph
        graph = gatb::core::debruijn::impl::Graph::loadAsPointer(graphPrefix+".h5");
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};
BOOST_CLASS_VERSION(SolidKmerOracle, 1)

#endif //KSGATB_SOLIDKMERORACLE_H
