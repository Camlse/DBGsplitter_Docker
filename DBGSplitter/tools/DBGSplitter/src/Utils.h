//
// Created by Leandro Ishi Soares de Lima on 05/06/16.
//

#ifndef KISSPLICE_UTILS_H
#define KISSPLICE_UTILS_H

#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <set>
#include <map>
#include <ctime>
#include <sstream>
#include <iostream>
#include <limits>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <whereami.h>
#include <pstream.h>
#include "global.h"

using namespace std;

//TODO: include utils here, not replicate the function
//#include <../../../modules/Utils.h>
char complement(char b);
string reverse_complement(const string &seq);
//Read all strings in the readsFile file and return them as a vector of strings
vector<string> readFileAsVectorString(const string &readsFile);

string toUpperCase (const string &s);

string getTime();

//Read all strings in the readsFile file and return them as a vector of strings
vector<string> getVectorStringFromFile(const string &readsFile);

void openFileForReading(const string &filePath, ifstream &stream);
void openFileForWriting(const string &filePath, ofstream &stream);

void fatalError (const string &message);

long countNbNodes (const string &unitigSizeFilename);

//read removed edges into a set of string
set<string> readRemovedEdgesIntoASetOfString(const string &filename);

void executeCommand(const string &command, bool verbose=true, const string &messageIfItFails="");

//get dir where the tool is installed
string getDirWhereToolIsInstalled();

//compute edit distance between s1 and s2
int computeEditDistance(const string &s1, const string &s2);
#endif //KISSPLICE_UTILS_H
