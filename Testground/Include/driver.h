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

//ExperimentResult test_RACE(size_t n_hashes, size_t dimensions, int w, int hash_power, size_t hash_range, std::ifstream& data, std::ifstream& queries, std::vector<double>& estimates);

std::vector<int> getWs();

std::vector<int> getPowers();


std::vector< std::pair<int, int> > clusteringExpMinHashStream(size_t n_hashes, int hash_power, size_t hash_range, std::ifstream& data, size_t threshold, std::vector<int>& clusters, std::string& path, size_t dimensions, size_t n_exp);

// In the function below, the numbers printed in the console do not correspond to clusters, but instead to the read number in the fasta file (first read is labeled 0).
std::vector< std::pair<int, int> > clusteringExpMinHashStreamMurmur(size_t n_exp, size_t n_hashes, int hash_power, size_t hash_range, std::ifstream& data, std::ifstream& labelIn, size_t n_samples_per_bucket, std::vector<int>& clusters, std::string& path, size_t dimensions, size_t total, int k);

std::vector< std::pair<int, int> > RSExpStream(size_t n_hashes, std::ifstream& data, size_t n_samples, std::vector<int>& clusters, size_t dimensions, size_t total);
