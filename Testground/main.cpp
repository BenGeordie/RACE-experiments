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
    
    std::ifstream in("/Users/benitogeordie/bg31/Research/Compression\ Project/promethion_kmeans_ground_truth.csv");
    size_t n_hashes = 1;
    std::string pathRACE = "/Users/benitogeordie/bg31/Research/Compression\ Project/promethion_after_race_minhash.txt";
    std::string pathRS = "/Users/benitogeordie/bg31/Research/Compression\ Project/promethion_after_rs.txt";
    for(int hash_power=4; hash_power<5; ++hash_power){
        for(size_t hash_range : {400}){
            std::cout << " p = " << hash_power << " range = " << hash_range << '\n';
            for(size_t n_samples_per_bucket = 9; n_samples_per_bucket < 10; ++n_samples_per_bucket){
                auto start = high_resolution_clock::now();
                clusteringExpMinHashStream(n_hashes, hash_power, hash_range, in, n_samples_per_bucket, clusters, pathRACE, 1024, n_samples_per_bucket*100);
                auto stop = high_resolution_clock::now();
                auto duration = duration_cast<microseconds>(stop - start);
                std::cout << "Time taken by function: "
                << duration.count() << " microseconds" << std::endl;
            }
        }
    }
    
//    std::cout << "RS:" << std::endl;
//    for(size_t samples = 300; samples < 1000; samples += 100){
//        clusteringExpMinHashStream(n_hashes, 1, 1, in, samples, clusters, pathRS, 1024, 3600);
//    }
}


