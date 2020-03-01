#pragma once

#include "io.h"
#include "util.h"
#include "reservoir.h"

#include "L2Hash.h"
#include "MinHash.h"
#include "MurmurHash.h"

#include <vector>
#include <iostream>
#include <fstream>

#include <chrono>
#include <map>
#include <algorithm>

void testRACEKrakenInput(size_t n_hashes, int hash_power, size_t hash_range, std::ifstream& data, size_t threshold, std::string& path, size_t dimensions, size_t n_exp, std::string fastWhat, std::string alphabet, int k, bool printCounters);
