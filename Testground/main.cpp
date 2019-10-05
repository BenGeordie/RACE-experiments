#include "Debug.h"
#include "io.h"
#include "driver.h"
#include <dirent.h>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

typedef int w;
typedef int power;
typedef int clusters;

int main(){
    
    // CODE TO ITERATE THROUGH FILES IN DIRECTORY
    
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir ("/Users/benitogeordie/bg31/Research/gridion_single")) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL) {
            printf ("%s\n", ent->d_name);
        }
        closedir (dir);
    } else {
        /* could not open directory */
        perror ("");
        return EXIT_FAILURE;
    }
    
    // CODE TO GET RESULT OF CLI
    
    const char* cmd = "/Users/benitogeordie/anaconda3/bin/h5dump -y -d /Raw/Reads/Read_10028/Signal /Users/benitogeordie/bg31/Research/gridion_single/0a3eaef6-5643-486e-b659-3e964adbd38f.fast5";
    std::cout << exec(cmd);

    // RUN RACE
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
//        clusteringExpMinHashStream(n_hashes, 1, 1, labelIn, samples, clusters, pathRS, 1024, 3600);
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


