#include "io.h"
#include "SequenceMinHash.h"
#include "RACE.h"
#include "util.h"

#include <chrono>
#include <sstream>
#include <vector>

/*
Writes throughput timing information to clog and writes the output samples to cout
*/

int main(int argc, char **argv){
    if (argc < 8){
        std::clog<<"Usage: "<<std::endl; 
        std::clog<<"RACEThresholds <tau> <race_range> <race_repetitions> <hash_power> <k> <data> <output>"<<std::endl; 
        std::clog<<"tau: (csv list of one or more floating point RACE sampling thresholds. Note: Don't include spaces between the elements)"<<std::endl; 
        std::clog<<"race_range: (integer range for each ACE)"<<std::endl; 
        std::clog<<"race_repetitions: (integer number of ACE repetitions)"<<std::endl; 
        std::clog<<"hash_power: (integer number of minhashes to compute for each ACE)"<<std::endl; 
        std::clog<<"k: (integer that decides how long the hashable k-mers are)"<<std::endl; 
        std::clog<<"data: (filename of Kraken-classified fastq file)"<<std::endl; 
        std::clog<<"output: (filename of fastq output file)"<<std::endl; 
        std::clog<<std::endl<<"Example usage:"<<std::endl; 
        std::clog<<std::endl<<"RACEThresholds 1.0,2.0,2.5,15.0 50 4 1 6 input.fastq output"<<std::endl; 
        std::clog<<"\t Outputs: output-1.0.fastq"<<std::endl; 
        std::clog<<std::endl<<"RACEThresholds 1.0 50 4 1 6 input.fastq output.fastq"<<std::endl; 
        std::clog<<"\t Outputs: output-1.0.fastq"<<std::endl; 
        return -1; 
    }

    // handle multiple tau's
    // double tau = std::stod(argv[1]);
    std::vector<double> taus; 

    std::stringstream ss(argv[1]); // get set of taus
    double tau = 0;
    while(ss >> tau){
        taus.push_back(tau); 
        if (ss.peek() == ',')
            ss.ignore();
    }
    ss.clear();


    int race_range = std::stoi(argv[2]);
    int race_repetitions = std::stoi(argv[3]);
    int hash_power = std::stoi(argv[4]);
    int kmer_k = std::stoi(argv[5]);

    std::ifstream datastream(argv[6]);

    // Create vector of ofstreams - one for each tau
    std::vector<std::ofstream> samplestreams; 
    std::string baseoutputfilename(argv[7]); 
    for(size_t i = 0; i < taus.size(); i++){
        std::string filename = baseoutputfilename; 
        filename += "-"; 
        filename += std::to_string(taus[i]); 
        filename += ".fastq"; 
        std::ofstream s(filename); 
        samplestreams.push_back(std::move(s));
    }
    // std::ofstream samplestream(argv[7]);

    // reset stream
    datastream.clear();
    datastream.seekg(0, std::ios::beg);

    std::string sequence;
    int label;


    SequenceMinHash hash = SequenceMinHash(race_repetitions*hash_power);
    int* raw_hashes = new int[race_repetitions*hash_power]; 
    int* rehashes = new int[race_repetitions];

    RACE sketch = RACE(race_repetitions,race_range); 

    int i = 0; 
    auto all_start = std::chrono::high_resolution_clock::now();
    do{
        auto start = std::chrono::high_resolution_clock::now();
        bool success = KrakenSequenceFeatures(datastream, sequence, label, "fastq");
        if (!success) continue;
        // std::cout<<label<<":\t"<<sequence<<std::endl;

        // auto start = std::chrono::high_resolution_clock::now();
        // begin timing here: 

        hash.getHash(kmer_k, sequence, raw_hashes); 
        // now that we have the sequence and label
        // feed the sequence into the RACE structure
        // first rehash so that the arrays can fit into RACE
        rehash(raw_hashes, rehashes, race_repetitions, hash_power);
        // then simultaneously query and add 
        double KDE = sketch.query_and_add(rehashes); 
        // note: KDE is on a scale from [0,N] not the normalized interval [0,1]

        for(size_t i = 0; i < taus.size(); i++){
            double tau = taus[i]; 
                if (KDE < tau){
                // then keep this sample 
                samplestreams[i]<<label;
                samplestreams[i]<<'\t';
                samplestreams[i]<<sequence; 
                samplestreams[i]<<std::endl; 
            }
        }
        auto finish = std::chrono::high_resolution_clock::now();
        std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count() << ",";
        i++; 
    }
    while(datastream);

    auto all_finish = std::chrono::high_resolution_clock::now();
    std::clog <<"Processed "<<i<<" sequences in "<< std::chrono::duration_cast<std::chrono::seconds>(all_finish-all_start).count()<<" seconds";

}

