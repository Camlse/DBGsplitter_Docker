//
// Created by Leandro Ishi Soares de Lima on 05/06/16.
//

#include "Utils.h"

using namespace std;

char complement(char b)
{
  switch(b)
  {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'G': return 'C';
    case 'C': return 'G';

    case 'a': return 't';
    case 't': return 'a';
    case 'g': return 'c';
    case 'c': return 'g';

    case 'N': return 'N';
    case '*': return '*';
  }
  return '?';
}

string reverse_complement(const string &seq)
{
  string s(seq.begin(),seq.end());
  string::iterator pos;

  reverse(s.begin(), s.end());

  for(pos=s.begin();pos!=s.end();++pos)
    *pos=complement(*pos);

  return s;
}


//Read all strings in the readsFile file and return them as a vector of strings
vector<string> readFileAsVectorString(const string &readsFile) {
  vector<string> allReadFilesNames;
  string tempStr;

  ifstream readsFileStream(readsFile, ifstream::in);
  while (getline(readsFileStream, tempStr)) {
    if (tempStr.size() > 0)
      allReadFilesNames.push_back(tempStr);
  }
  readsFileStream.close();

  return allReadFilesNames;
}



string getTime() {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[1024];
  tstruct = *localtime(&now);
  // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
  // for more information about date/time format
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

  return buf;
}



string toUpperCase (const string &s) {
  string toReturn(s);
  //transform the string
  for (int j=0;j<s.size();j++)
    toReturn[j]=toupper(s[j]);
  return toReturn;
}


//Read all strings in the readsFile file and return them as a vector of strings
vector<string> getVectorStringFromFile(const string &readsFile) {
  vector<string> allReadFilesNames;
  string tempStr;

  ifstream readsFileStream;
  openFileForReading(readsFile, readsFileStream);
  while (getline(readsFileStream, tempStr)) {
    if (tempStr.size() > 0)
      allReadFilesNames.push_back(tempStr);
  }
  readsFileStream.close();

  return allReadFilesNames;
}

void openFileForReading(const string &filePath, ifstream &stream) {
  stream.open(filePath);
  if (!stream.is_open()) {
    stringstream ss;
    ss << "Error opening file " << filePath;
    fatalError(ss.str());
  }
}

void openFileForWriting(const string &filePath, ofstream &stream) {
  stream.open(filePath);
  if (!stream.is_open()) {
    stringstream ss;
    ss << "Error opening file " << filePath;
    fatalError(ss.str());
  }
}

void fatalError (const string &message) {
  cerr << endl << endl << "[FATAL ERROR] " << message << endl << endl;
  cerr.flush();
  exit(1);
}


long countNbNodes (const string &unitigSizeFilename) {
  ifstream unitigSizeFile;
  openFileForReading(unitigSizeFilename, unitigSizeFile);
  long nbNodes;
  unitigSizeFile >> nbNodes;
  unitigSizeFile.close();
  return nbNodes;
}


//read removed edges into a set of string
set<string> readRemovedEdgesIntoASetOfString(const string &filename) {
  set<string> removedEdges;
  vector<string> removedEdgesAsVector = getVectorStringFromFile(filename);
  removedEdges.insert(removedEdgesAsVector.begin(), removedEdgesAsVector.end());
  return removedEdges;
}

//execute a command
void executeCommand(const string &command, bool verbose, const string &messageIfItFails) {
  // run a process and create a streambuf that reads its stdout and stderr
  if (verbose)
    cerr << "Executing " << command << "..." << endl;

  //create the process
  redi::ipstream proc(command, redi::pstreams::pstdout | redi::pstreams::pstderr);
  string line;

  // read child's stdout
  while (getline(proc.out(), line)) {
    if (verbose)
      cout << line << endl;
  }
  // read child's stderr
  while (getline(proc.err(), line)) {
    if (verbose)
      cerr << line << endl;
  }

  //check exit status
  proc.close();
  if (proc.rdbuf()->exited()) {
    if (proc.rdbuf()->status() != 0) {
      stringstream ss;
      ss << "Error executing " << command << ". Exit status: " << proc.rdbuf()->status() << endl;
      if (messageIfItFails != "")
        ss << "Message: " << messageIfItFails << endl;
      fatalError(ss.str());
    }
    if (verbose)
      cerr << "Executing " << command << " - Done!" << endl;
  }
  else
    fatalError("On executeCommand()");
}


string getDirWhereToolIsInstalled() {
  char* path = NULL;
  int length, dirname_length;
  string toReturn;

  length = wai_getExecutablePath(NULL, 0, &dirname_length);
  if (length > 0)
  {
    path = (char*)malloc(length + 1);
    if (!path)
      //error, malloc did not work
      fatalError("Error on Utils.cpp::getDirWhereDBGWASIsInstalled()");

    //get the executable path
    wai_getExecutablePath(path, length, &dirname_length);

    //get the executable dir
    path[length] = '\0';
    path[dirname_length] = '\0';

    //save it in toReturn
    toReturn = string(path);

    //free the memory
    free(path);
  }
  else {
    fatalError("Error on Utils.cpp::getDirWhereDBGWASIsInstalled()");
  }

  return toReturn;
}


//compute edit distance between s1 and s2
int computeEditDistance(const string &s1, const string &s2) {
  const std::size_t len1 = s1.size(), len2 = s2.size();
  std::vector<int> col(len2+1), prevCol(len2+1);

  for (int i = 0; i < prevCol.size(); i++)
    prevCol[i] = i;
  for (int i = 0; i < len1; i++) {
    col[0] = i+1;
    for (int j = 0; j < len2; j++)
      // note that std::min({arg1, arg2, arg3}) works only in C++11,
      // for C++98 use std::min(std::min(arg1, arg2), arg3)
      //TODO: treat the evetuel Ns that we can have here?
      col[j+1] = std::min({ prevCol[1 + j] + 1, col[j] + 1, prevCol[j] + (s1[i]==s2[j] ? 0 : 1) });
    col.swap(prevCol);
  }
  return prevCol[len2];
}

