#include "Debug.h"
#include "io.h"
#include "driver.h"
#include <dirent.h>
#include <chrono>
//#include "Include/catch.hpp"
//
//#define CATCH_CONFIG_MAIN

typedef int w;
typedef int power;
typedef int clusters;

int main(){
    
//    // TESTING IO
//    std::string sequence = "AAAAA";
//    std::vector<int> vec;
//    int k = 3;
//    KmerizeMurmur(sequence, vec, k);
//    
//    std::string sequence2 = "AAABA";
//    std::vector<int> vec2;
//    KmerizeMurmur(sequence2, vec2, k);
//    
//
//
//    std::string kmer = "AAA";
//    std::vector<int> ref;
//    ref.push_back(MurmurHash(&kmer, (k + 1) * sizeof(char), k));
//    std::cout << "vec contents" << std::endl;
//    std::copy(vec.begin(), vec.end(), std::ostream_iterator<int>(std::cout, " "));
//    std::cout << std::endl;
//    std::cout << "vec2 contents" << std::endl;
//    std::copy(vec2.begin(), vec2.end(), std::ostream_iterator<int>(std::cout, " "));
//    std::cout << std::endl;
//    std::cout << "ref contents" << std::endl;
//    std::copy(ref.begin(), ref.end(), std::ostream_iterator<int>(std::cout, " "));
//    std::cout << std::endl;
    
    
    
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
    for(size_t samples = 400; samples < 500; samples += 100){
//        clusteringExpMinHashStreamMurmur(n_hashes, 1, 1, in, labelIn, samples, clusters, pathRS, 1024, 3600, 5);
        clusteringExpMinHashStream(n_hashes, 1, 1, labelIn, samples, clusters, pathRS, 1024, 3600);
    }
    for(int hash_power=4; hash_power<5; ++hash_power){
        for(size_t hash_range : {200}){
            std::cout << " p = " << hash_power << " range = " << hash_range << '\n';
            for(size_t n_samples_per_bucket = 4; n_samples_per_bucket < 5; ++n_samples_per_bucket){
                clusteringExpMinHashStream(n_hashes, hash_power, hash_range, labelIn, n_samples_per_bucket, clusters, pathRACE, 1024, n_samples_per_bucket*100);
//                clusteringExpMinHashStreamMurmur(n_hashes, hash_power, hash_range, in, labelIn, n_samples_per_bucket, clusters, pathRACE, 1024, n_samples_per_bucket*100, 6);
            }
        }
    }

}


