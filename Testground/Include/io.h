#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "MurmurHash.h"
#include "SignedRandomProjections.h"

std::string exec(const char* cmd);

void VectorFeatures(std::istream& in, std::vector<double>& vec, int& label, size_t& dimensions);
/*
Reads a CSV-coded set
element0, element1, element2, ... elementN
from stream "in"
Clears vec and and places element0, ... into vec
*/

void VectorFeaturesAdjacent(std::istream& in, std::vector<int>& vec, double& label);

void kmerize(std::string sequence, std::vector<int>& vec, std::string alphabet, int k);

void KmerizeMurmur(std::string sequence, std::vector<int>& vec, int k);

void KmerizeSquiggleSRPSliding(std::vector<double> &squiggle, std::vector<int>& vec, int dim, int K, int L);

void KmerizeSquiggleSRPPartition(std::vector<double> &squiggle, std::vector<int>& vec, int dim, int K, int L);

void VectorFeaturesSquiggleCompiled(std::istream& in, std::vector<int>& vec, std::istream& labelIn, double& label, int& dim, int K, int L);

void VectorFeaturesFastaMurmur(std::istream& in, std::vector<int>& vec, std::istream& labelIn, double& label, int& k);

void VectorFeaturesFasta(std::istream& in, std::vector<int>& vec, std::string& label, std::string& alphabet, int& k);

void WriteCSVResults(std::ostream& out, size_t sketch_size, double preprocessing_time, double query_time, std::vector<double>& estimates); 
/*
Writes a CSV file row 
Assumes out starts at the start of the row (does not write leading newline)
Writes the following: 
[sketch_size in cells], [preprocessing time in seconds], [total querying time in seconds], KDE1, KDE2, KDE3, KDE4, (all doubles)

*/

void printList(std::vector< double >& list);

void printList(std::vector< std::chrono::microseconds >& list);

std::vector<std::pair< std::vector<double>, int>> open_CSV(std::string path, int width);

std::vector<std::pair< std::vector<double>, int>> open_sparse(std::string path, int smallest_id, int largest_id);

std::string path_to_root();




