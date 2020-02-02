#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "MurmurHash.h"
#include "SignedRandomProjections.h"


void kmerize(std::string sequence, std::vector<int>& vec, std::string alphabet, int k);

void VectorFeaturesFastKraken(std::istream& in, std::vector<int>& vec, double& label, std::string& alphabet, int& k, std::string& fastWhat);

bool KrakenSequenceFeatures(std::istream& in, std::string& sequence, int& label, std::string fastWhat); 






