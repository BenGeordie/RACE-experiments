#include "io.h"
#include "SequenceMinHash.h"
#include "RACE.h"
#include "util.h"

#include <chrono>
#include <vector>
#include <string>
#include <random>
#include <cstring>

/*
Writes throughput timing information to clog and writes the output samples to cout
*/

int main(int argc, char **argv){
    if (argc < 4){
        std::clog<<"Usage: "<<std::endl; 
        std::clog<<"RSReservoir <n_samples> <data> <output>"<<std::endl; 
        std::clog<<"n_samples: (integer number of samples stored)"<<std::endl; 
        std::clog<<"data: (filename of Kraken-classified fastq file)"<<std::endl; 
        std::clog<<"output: (filename of fastq output file)"<<std::endl; 
        std::clog<<"[--seed <int>]: (integer seed for random number generator. Optional. Defaults to time(NULL))"<<std::endl; 
        std::clog<<std::endl<<"Example usage:"<<std::endl; 
        std::clog<<std::endl<<"RSReservoir 200 input.fastq output.txt --seed 420"<<std::endl; 
        return -1; 
    }
    int n_samples = std::stod(argv[1]);

    std::ifstream datastream(argv[2]);
    std::ofstream samplestream(argv[3]);

    unsigned int seed = time(NULL); 
    for (int i = 0; i < argc; ++i){
        if (std::strcmp("--seed",argv[i]) == 0){
            if ((i+1) < argc){
                seed = std::stoi(argv[i+1]);
            } else {
                std::cerr<<"Invalid seed"<<std::endl; 
                return -1;
            }
        }
    }

    // reset stream
    datastream.clear();
    datastream.seekg(0, std::ios::beg);

    std::string sequence;
    int label;

    std::mt19937 generator(seed);

    std::vector< std::string > sequences;
    std::vector< int > labels;

    auto all_start = std::chrono::high_resolution_clock::now();

    // ridiculously ad-hoc implementation of reservoir sampling
    int i = 0;
    do{
        auto start = std::chrono::high_resolution_clock::now();
        bool success = KrakenSequenceFeatures(datastream, sequence, label, "fastq");
        if (!success) continue;


        if (i < n_samples){
                sequences.push_back(sequence);
                labels.push_back(label);
        } else {
            std::uniform_int_distribution<int> distribution(0,i-1);
            int j = distribution(generator); 
            if (j < n_samples){
                sequences[j] = sequence;
                labels[j] = label; 
            }
        }
        auto finish = std::chrono::high_resolution_clock::now();
        std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count() << ",";

        i++;
    }
    while(datastream);

    auto all_finish = std::chrono::high_resolution_clock::now();
    std::clog <<"Processed "<<i<<" sequences in "<< std::chrono::duration_cast<std::chrono::seconds>(all_finish-all_start).count()<<" seconds";


    for(int i = 0; i < sequences.size(); i++){
        samplestream<<labels[i];
        samplestream<<'\t';
        samplestream<<sequences[i]; 
        samplestream<<std::endl; 
    }

}
