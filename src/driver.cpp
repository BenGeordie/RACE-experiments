#include "driver.h"

void testRACEKrakenInput(size_t n_hashes, int hash_power, size_t hash_range, std::ifstream& data, size_t threshold, std::string& path, size_t dimensions, size_t n_exp, std::string fastWhat, std::string alphabet, int k, bool printCounters) {
    
    // For latency
    // using namespace std::chrono;
    // std::vector<microseconds> latencies;
    
    // reset stream
    data.clear();
    data.seekg(0, std::ios::beg);
    
    // input buffer for vector
    std::vector<int> vec;
    
    // input buffer for label
    double label;
    std::vector<double> labelvec;
    
    // Initialize reservoirs (run n_exp of RACEs with R = n_hashes, hash_range buckets per function, accepts vector if average of counters in hash buckets <= threshold)
    Reservoir reservoirs(n_hashes, hash_range, threshold, n_exp);

    // Initialize hash tables and hash buffer
    MinHash hash = MinHash(dimensions, n_hashes*hash_power*n_exp);
    int* raw_hashes = new int[n_hashes*hash_power*n_exp];
    int* rehashes = new int[n_hashes*n_exp];
    
    int idx = 0;
    do{
        // auto start = high_resolution_clock::now();
        
        // Stream a fasta/q entry and produce a vector of kmer token numbers
        VectorFeaturesFastKraken(data, vec, label, alphabet, k, fastWhat);
        if (vec.size() == 0)continue;
        
        // Hash em
        hash.getHash(vec,raw_hashes);
        rehash(raw_hashes, rehashes, n_hashes*n_exp, hash_power);
        
        // For the moment, we are adding the labels instead of the actual vector for ease of analysis.
        labelvec.clear();
        labelvec.push_back(label);
        labelvec.push_back((double) idx);
        reservoirs.add(labelvec, rehashes);
        if(idx % 1000 == 0)
            std::cout << '.' << std::flush;
        idx++;
        // auto stop = high_resolution_clock::now();
        // auto duration = duration_cast<microseconds>(stop - start);
        // latencies.push_back(duration);
    }
    while(data);
    reservoirs.pprint(std::cout, path, printCounters);
    //    printList(latencies);

}


