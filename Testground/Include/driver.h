#pragma once

#include "io.h"
#include "util.h"
#include "RACE.h"
#include "reservoir.h"

#include "L2Hash.h"
#include "MinHash.h"
#include "MurmurHash.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <dirent.h>

#include <chrono>


/*
Driver for the experiments that we will run

Basically, this just wraps RACE, RS, SKA, HBE
and any other method we might choose to implement
in an easy-to-call function that conducts the experiment
and deposits the results in a format that is easy to deal with
*/


// convenience struct for packaging data we want
// to plot or analyze later in experiments
struct ExperimentResult
{
    size_t sketch_size;
    double preprocessing_time;
    double query_time;
};

std::vector< std::pair<int, int> > clusteringExpMinHashStreamSquiggle(size_t n_hashes, int hash_power, size_t hash_range, std::ifstream& data, std::ifstream& labelIn, size_t n_samples_per_bucket, std::vector<int>& clusters, std::string& path, size_t dimensions, size_t total, int dim, int K);

std::vector< std::pair<int, int> > clusteringExpMinHashStreamSquiggle(size_t nHashes, size_t nPerBucket, size_t total, size_t minHashRange, int minHashPower, size_t minHashDim, int srpPower, int srpDim, const char* fast5Dir, std::string h5dumpPath, std::string h5lsPath, std::ifstream& labelIn, std::vector<int>& clusters, std::string& outputPath);
