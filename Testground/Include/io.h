#pragma once
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <exception>
#include "MurmurHash.h"
#include "SignedRandomProjections.h"

std::string exec(const char* cmd);

void kmerizeSquiggleSRPSliding(std::vector<double> &squiggle, std::vector<int>& vec, int dim, int K, int L);

void kmerizeSquiggleSRPPartition(std::vector<double> &squiggle, std::vector<int>& vec, int dim, int K, int L);

void vectorFeaturesSquiggleCompiled(std::istream& in, std::vector<int>& vec, std::istream& labelIn, double& label, int& dim, int K, int L);

void parseSignalDirectory(std::istringstream& cliOutput, std::string& path);

void parseSquiggle(std::istringstream& cliOutput, std::vector<double>& vec, std::string& temp);

void getSquiggleVector(std::string fileName, std::vector<double>& vec, std::string& h5lsPath, std::string& h5dumpPath);

void getLabel(std::istream& labelIn, double& label);

void writeCSVResults(std::ostream& out, size_t sketch_size, double preprocessing_time, double query_time, std::vector<double>& estimates);
/*
Writes a CSV file row 
Assumes out starts at the start of the row (does not write leading newline)
Writes the following: 
[sketch_size in cells], [preprocessing time in seconds], [total querying time in seconds], KDE1, KDE2, KDE3, KDE4, (all doubles)

*/

void printList(std::vector< double >& list);

void printList(std::vector< int >& list);

void printList(std::vector< std::chrono::microseconds >& list);

std::string path_to_root();




