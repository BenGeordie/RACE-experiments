#include "Debug.h"
#include "io.h"
#include "driver.h"
#include <chrono>

/*
Usage:
RACEReservoirs <data> <b> <range> <nhashes>

Writes throughput timing information to clog and writes the output samples to cout

*/

int main(int argc, char **argv){

    if (argc < 3){
        std::clog<<"Usage: "<<std::endl;
        std::clog<<"RACEReservoirs <data> <b> <range> <nhashes> <k>"<<std::endl; 
        std::clog<<"data: (filename of Kraken-classified fastq file)"<<std::endl; 
        std::clog<<"b: (integer number of buckets in each reservoir)"<<std::endl; 
        std::clog<<"range: (integer range of ACE)"<<std::endl; 
        std::clog<<"nhashes: (integer number of minhashes)"<<std::endl; 
        std::clog<<"k: (integer k-mer )"<<std::endl; 



        std::clog<<"RACEReservoirs <data> <b> <range> <k>"<<std::endl; 


        std::clog<<"--distanceID <int> (default 0)"<<std::endl; 
        std::clog<<"--skipD <int> (skips first columns in data vector, default 0)"<<std::endl; 
        std::clog<<"Supported distanceIDs:"<<std::endl;
        std::clog<<"        0: L2 (default)"<<std::endl;
        std::clog<<"        1: L1 "<<std::endl;
        std::clog<<"        2: Angular distance"<<std::endl;
        return 0;
    }

}



typedef int w;
typedef int power;
typedef int clusters;

int main(){
    
    std::ifstream fastqIn("/Users/benitogeordie/bg31/Research/ENA_datasets/SRR1056036/SRR1056036_classified.fastq");
    size_t n_exp = 50;
    std::string pathRACE = "/Users/benitogeordie/bg31/Research/Compression Project/promethion_after_race_minhash.txt";

    for (size_t n_hashes : {3}) {
        for(int hash_power : {4}){
            for(size_t hash_range : {100}){
                std::cout << "RACE" << std::endl;
                std::cout << " R = " << n_hashes << " p = " << hash_power << " range = " << hash_range << '\n';
                for(size_t threshold : {10, 20, 30, 40, 50}){
                    std::cout << "threshold = " << threshold << '\n';
                    std::cout << "Not Pretokenized" << '\n';
                    testRACEKrakenInput(n_hashes, hash_power, hash_range, fastqIn, threshold, pathRACE, 1024, n_exp, "fastq", "ACGTN", 5, true);
                }
            }
        }
    }
}

