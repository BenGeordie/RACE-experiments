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

void testRACEKrakenInput(size_t n_hashes, int hash_power, size_t hash_range, std::ifstream& data, size_t threshold, std::string& path, size_t dimensions, size_t n_exp, std::string fastWhat, std::string alphabet, int k, bool printCounters);
