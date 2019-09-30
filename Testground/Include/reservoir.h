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



typedef std::vector< std::vector< double > > reserve;
typedef std::map<int, int> statistics;
typedef int remainders;
// Should I reserve space for vectors?

class Reservoir
{
public:
    Reservoir(size_t R, size_t range, size_t L, std::vector<int>& clusters, size_t total);
    ~Reservoir();

    void add(std::vector< double > vec, int *hashes);
    void clear();

    void pprint(std::ostream& out, int width, bool format, std::vector<int>& clusters);
    void pprint(std::ostream& out, std::string& path);
    int count();
    int countLabels(size_t& r);
    int countElements(size_t& r);


    private:
        size_t _R, _range, _L;
        reserve* _buckets;
        statistics* _stats;
        remainders* _rems;
        const uint8_t magic_number = 0x4D; // magic number for binary file IO
        const uint8_t file_version_number = 0x01; // file version number
};

