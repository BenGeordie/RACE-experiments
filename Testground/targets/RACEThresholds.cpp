#include "Debug.h"
#include "io.h"
#include "driver.h"
#include <chrono>

/*
Writes throughput timing information to clog and writes the output samples to cout
*/

int main(int argc, char **argv){
    if (argc < 8){
        std::clog<<"Usage: "<<std::endl; 
        std::clog<<"RACEThresholds <tau> <race_range> <race_repetitions> <hash_power> <k> <data> <output>"<<std::endl; 
        std::clog<<"tau: (floating point RACE sampling threshold)"<<std::endl; 
        std::clog<<"race_range: (integer range for each ACE)"<<std::endl; 
        std::clog<<"race_repetitions: (integer number of ACE repetitions)"<<std::endl; 
        std::clog<<"hash_power: (integer number of minhashes to compute for each ACE)"<<std::endl; 
        std::clog<<"k: (integer that decides how long the hashable k-mers are)"<<std::endl; 
        std::clog<<"data: (filename of Kraken-classified fastq file)"<<std::endl; 
        std::clog<<"output: (filename of fastq output file)"<<std::endl; 
        std::clog<<std::endl<<"Example usage:"<<std::endl; 
        std::clog<<std::endl<<"RACEThresholds 15.0 50 4 1 6 input.fastq output.fastq"<<std::endl; 
        return -1; 
    }

    double tau = std::stod(argv[1]);
    int race_range = std::stoi(argv[2]);
    int race_repetitions = std::stoi(argv[3]);
    int hash_power = std::stoi(argv[4]);
    int kmer_k = std::stoi(argv[5]);

    std::ifstream datastream(argv[6]);
    std::ofstream samplestream(argv[7]);

    // reset stream
    datastream.clear(); 
    datastream.seekg(0, std::ios::beg); 






}
    
    // // For latency
    // // using namespace std::chrono;
    // // std::vector<microseconds> latencies;
    
    // // reset stream
    // data.clear();
    // data.seekg(0, std::ios::beg);
    
    // // input buffer for vector
    // std::vector<int> vec;
    
    // // input buffer for label
    // double label;
    // std::vector<double> labelvec;
    
    // // Initialize reservoirs (run n_exp of RACEs with R = n_hashes, hash_range buckets per function, accepts vector if average of counters in hash buckets <= threshold)
    // Reservoir reservoirs(n_hashes, hash_range, threshold, n_exp);

    // // Initialize hash tables and hash buffer
    // MinHash hash = MinHash(dimensions, n_hashes*hash_power*n_exp);
    // int* raw_hashes = new int[n_hashes*hash_power*n_exp];
    // int* rehashes = new int[n_hashes*n_exp];
    
    // int idx = 0;
    // do{
    //     // auto start = high_resolution_clock::now();
        
    //     // Stream a fasta/q entry and produce a vector of kmer token numbers
    //     VectorFeaturesFastKraken(data, vec, label, alphabet, k, fastWhat);
    //     if (vec.size() == 0)continue;
        
    //     // Hash em
    //     hash.getHash(vec,raw_hashes);
    //     rehash(raw_hashes, rehashes, n_hashes*n_exp, hash_power);
        
    //     // For the moment, we are adding the labels instead of the actual vector for ease of analysis.
    //     labelvec.clear();
    //     labelvec.push_back(label);
    //     labelvec.push_back((double) idx);
    //     reservoirs.add(labelvec, rehashes);
    //     if(idx % 1000 == 0)
    //         std::cout << '.' << std::flush;
    //     idx++;
    //     // auto stop = high_resolution_clock::now();
    //     // auto duration = duration_cast<microseconds>(stop - start);
    //     // latencies.push_back(duration);
    // }
    // while(data);
    // reservoirs.pprint(std::cout, path, printCounters);
    // //    printList(latencies);
