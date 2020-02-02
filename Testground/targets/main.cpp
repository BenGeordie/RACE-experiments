#include "Debug.h"
#include "io.h"
#include "driver.h"
#include <chrono>

typedef int w;
typedef int power;
typedef int clusters;

int main(){
    
    std::ifstream fastqIn("/data/SRR1056036_classified.fastq");
    size_t n_exp = 50;
    std::string pathRACE = "/Users/benitogeordie/bg31/Research/Compression Project/promethion_after_race_minhash.txt";

    for (size_t n_hashes : {3}) { // race reps 
        for(int hash_power : {4}){ // number of minhashes that get rehashed 
            for(size_t hash_range : {100}){ // range of each race 
                std::cout << "RACE" << std::endl;
                std::cout << " R = " << n_hashes << " p = " << hash_power << " range = " << hash_range << '\n';
                for(size_t threshold : {10, 20, 30, 40, 50}){ // threshold 
                    std::cout << "threshold = " << threshold << '\n';
                    std::cout << "Not Pretokenized" << '\n';
                    // pathRACE is not current (deprecated) 
                    // minhash_dimensions (deprecated)
                    // n_exp : number of experiments
                    // 
                    testRACEKrakenInput(n_hashes, hash_power, hash_range, fastqIn, threshold, pathRACE, 1024, n_exp, "fastq", "ACGTN", 5, true);
                }
            }
        }
    }
}


