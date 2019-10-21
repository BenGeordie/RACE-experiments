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
#include "io.h"


// change reserve
typedef std::vector< std::pair< int, std::string > > reserve;
// Should I reserve space for vectors?

class Reservoir
{
public:
    Reservoir(size_t R, size_t range, size_t L, size_t total);
    ~Reservoir();
    
    void add(int label, int *hashes, std::string fastaEntry, std::ofstream& outputFasta);
    void clear();
    
    void pprint(std::ostream& outPath);
    void pprint(std::ostream& out, std::string& path);
    int count();
    int countLabels();
    int countElements(size_t& r);
    
    
private:
    size_t _R, _range, _L, _total;
    reserve _sample;
    std::set<int> _labelSet;
    int* _race;
    const uint8_t magic_number = 0x4D; // magic number for binary file IO
    const uint8_t file_version_number = 0x01; // file version number
};

