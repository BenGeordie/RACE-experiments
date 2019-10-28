#pragma once

#include <iostream>
#include <iomanip>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <string>

#include <cstring>
#include <random>

#include "Debug.h"



typedef std::vector< double > samples;
typedef int counters;

class Reservoir
{
public:
    Reservoir(size_t R, size_t range, size_t threshold, size_t n_exp);
    ~Reservoir();

    void add(std::vector< double > vec, int *hashes);
    void clear();

    void pprint(std::ostream& out, std::string& path, bool printCounters);

    private:
        size_t _R, _range, _n_exp, _threshold;
        counters* _counters;
        samples* _samples; // One "samples" vector for each trial.
    
        const uint8_t magic_number = 0x4D; // magic number for binary file IO
        const uint8_t file_version_number = 0x01; // file version number
};

