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
    
//    std::ifstream labelIn("/Users/benitogeordie/bg31/Research/Compression\ Project/promethion_kmeans_ground_truth.csv");
//    std::ifstream labelIn("/Users/benitogeordie/bg31/Research/Compression Project/promethion_classified_token_vectors_k5_taxid.csv");
    std::ifstream labelIn("/Users/benitogeordie/bg31/Research/Compression Project/SRR1056036_token_vectors_k5.csv");
    std::ifstream fastqIn("/Users/benitogeordie/bg31/Research/ENA_datasets/SRR1056036/SRR1056036_classified.fastq");
//    std::ifstream in("/Users/benitogeordie/Downloads/vboza-deepnano-e8a621e17b9f/r9/promethion_basecalled.fasta");
    size_t n_exp = 50;
    std::string pathRACE = "/Users/benitogeordie/bg31/Research/Compression\ Project/promethion_after_race_minhash.txt";
    std::string pathRS = "/Users/benitogeordie/bg31/Research/Compression\ Project/promethion_after_rs.txt";
    //    std::cout << "RS:" << std::endl;
    for(size_t samples = 800; samples < 900; samples += 100){
        std::cout << "RS" << std::endl;
//        clusteringExpMinHashStreamMurmur(n_hashes, 1, 1, in, labelIn, samples, clusters, pathRS, 1024, 3600, 5);
//        clusteringExpMinHashStream(n_hashes, 1, 1, labelIn, samples, clusters, pathRS, 1024, 3600);
    }
    for (size_t n_hashes : {3}) {
        for(int hash_power : {2}){
            for(size_t hash_range : {200}){
                std::cout << "RACE" << std::endl;
                std::cout << " R = " << n_hashes << " p = " << hash_power << " range = " << hash_range << '\n';
                for(size_t threshold : {10, 20, 30, 40, 50}){
                    std::cout << "threshold = " << threshold << '\n';
//                    std::cout << "Pretokenized" << '\n';
//                    clusteringExpMinHashStream(n_hashes, hash_power, hash_range, labelIn, threshold, pathRACE, 1024, n_exp);
                    std::cout << "Not Pretokenized" << '\n';
                    clusteringExpMinHashFastKraken(n_hashes, hash_power, hash_range, fastqIn, threshold, pathRACE, 1024, n_exp, "fastq", "ACGTN", 5);
                    //                clusteringExpMinHashStreamMurmur(n_hashes, hash_power, hash_range, in, labelIn, n_samples_per_bucket, clusters, pathRACE, 1024, n_samples_per_bucket*100, 6);
                }
            }
        }
    }
}


