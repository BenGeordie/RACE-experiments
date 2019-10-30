#include "io.h"
#include "SequenceMinHash.h"
#include "RACE.h"
#include "util.h"

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

    std::string sequence;
    int label;


    SequenceMinHash hash = SequenceMinHash(race_repetitions*hash_power);
    int* raw_hashes = new int[race_repetitions*hash_power]; 
    int* rehashes = new int[race_repetitions];

    RACE sketch = RACE(race_repetitions,race_range); 

    do{
        bool success = KrakenSequenceFeatures(datastream, sequence, label, "fastq");
        if (!success) continue;
        // std::cout<<label<<":\t"<<sequence<<std::endl;

        // begin timing here: 

        hash.getHash(kmer_k, sequence, raw_hashes); 
        // now that we have the sequence and label
        // feed the sequence into the RACE structure
        // first rehash so that the arrays can fit into RACE
        rehash(raw_hashes, rehashes, race_repetitions, hash_power);
        // then simultaneously query and add 
        double KDE = sketch.query_and_add(rehashes); 
        // note: KDE is on a scale from [0,N] not the normalized interval [0,1]
        if (KDE < tau){
            // then keep this sample 
            samplestream<<label;
            samplestream<<'\t';
            samplestream<<sequence; 
            samplestream<<std::endl; 
        }

    }
    while(datastream);
}
