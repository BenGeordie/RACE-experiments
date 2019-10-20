#include "Debug.h"
#include "io.h"
#include "driver.h"
#include <chrono>

typedef int w;
typedef int power;
typedef int clusters;

int main(){
    
    // Vector containing names of ground truth clusters:
    std::vector<int> clusters;
    for(int i = 0; i < 1000; i++)
        clusters.push_back(i);
    
    using namespace std::chrono;
    
    std::ifstream labelIn("/Users/benitogeordie/bg31/Research/Compression\ Project/promethion_kmeans_ground_truth.csv");
    std::ifstream in("/Users/benitogeordie/Downloads/vboza-deepnano-e8a621e17b9f/r9/promethion_basecalled.fasta");
    size_t n_hashes = 1;
    std::string pathRACE = "/Users/benitogeordie/bg31/Research/Compression\ Project/promethion_after_race_minhash.txt";
    std::string pathRS = "/Users/benitogeordie/bg31/Research/Compression\ Project/promethion_after_rs.txt";
    //    std::cout << "RS:" << std::endl;
    for(size_t samples = 800; samples < 900; samples += 100){
        std::cout << "RS" << std::endl;
//        clusteringExpMinHashStreamMurmur(n_hashes, 1, 1, in, labelIn, samples, clusters, pathRS, 1024, 3600, 5);
        clusteringExpMinHashStream(n_hashes, 1, 1, labelIn, samples, clusters, pathRS, 1024, 3600);
    }
    for(int hash_power=4; hash_power<5; ++hash_power){
        for(size_t hash_range : {200}){
            std::cout << "RACE" << std::endl;
            std::cout << " p = " << hash_power << " range = " << hash_range << '\n';
            for(size_t n_samples_per_bucket = 11; n_samples_per_bucket < 12; ++n_samples_per_bucket){
                clusteringExpMinHashStream(n_hashes, hash_power, hash_range, labelIn, n_samples_per_bucket, clusters, pathRACE, 1024, 800);
//                clusteringExpMinHashStreamMurmur(n_hashes, hash_power, hash_range, in, labelIn, n_samples_per_bucket, clusters, pathRACE, 1024, n_samples_per_bucket*100, 6);
            }
        }
    }

}


