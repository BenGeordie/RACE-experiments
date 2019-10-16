#include "driver.h"

//std::vector< std::pair<int, int> > clusteringExpMinHashStreamSquiggle(size_t n_hashes, int hash_power, size_t hash_range, std::ifstream& data, std::ifstream& labelIn, size_t n_samples_per_bucket, std::vector<int>& clusters, std::string& path, size_t dimensions, size_t total, int dim, int K){
//
//    using namespace std::chrono;
//    std::vector<microseconds> latencies;
//
//    // reset stream
//    data.clear();
//    data.seekg(0, std::ios::beg);
//
//    labelIn.clear();
//    labelIn.seekg(0, std::ios::beg);
//
//    // input buffer for vector
//    std::vector<int> vec;
//    // input buffer for label
//    double label;
//    std::vector<double> labelvec;
//
//    // output  vector
//    std::vector< std::pair<int, int> > output;
//
//    // Initialize reservoirs (R reservoirs of hash_power buckets,
//    // each containing n_samples samples, each of them a vector)
//    Reservoir reservoirs(n_hashes, hash_range, n_samples_per_bucket, clusters, total);
//
//    // Initialize hash tables
//    MinHash hash = MinHash(dimensions, n_hashes*hash_power);
//    int* raw_hashes = new int[n_hashes*hash_power];
//    int* rehashes = new int[n_hashes];
//
//    int idx = 0;
//    do{
//        auto start = high_resolution_clock::now();
//        vectorFeaturesSquiggleCompiled(data, vec, labelIn, label, dim, K, 1);
//        //        std::cout << label << std::endl;
//        if (vec.size() == 0)continue;
//
//        hash.getHash(vec, raw_hashes);
//        rehash(raw_hashes, rehashes, n_hashes, hash_power);
//
//
//        labelvec.clear();
//        labelvec.push_back(label);
//        reservoirs.add(labelvec, rehashes);
//        if(idx % 1000 == 0)
//            std::cout << '.' << std::flush;
//        idx++;
//        auto stop = high_resolution_clock::now();
//        auto duration = duration_cast<microseconds>(stop - start);
//        latencies.push_back(duration);
//    }
//    while(data);
//
//    //    for(size_t r = 0; r < n_hashes; r++){
//    //        int clus = reservoirs.countLabels(r); // Haven't adapted to MinHash format
//    //        //        std::cout << "clus : " << clus << std::endl << std::flush;
//    //        int elem = reservoirs.countElements(r);
//    //        output.push_back(std::make_pair(clus, elem));
//    //    }
//    reservoirs.pprint(std::cout, path);
//    printList(latencies);
//    return output;
//}

std::vector< std::pair<int, int> > clusteringExpMinHashStreamSquiggle(size_t nHashes, size_t nPerBucket, size_t total, size_t minHashRange, int minHashPower, size_t minHashDim, int srpPower, int srpDim, const char* fast5Dir, std::string& destinationDir, std::string h5dumpPath, std::string h5lsPath, std::ifstream& labelIn, std::vector<int>& clusters, std::string& outputPath){
    
    std::vector< std::pair<int, int> > output;
    std::vector<double> squiggles;
    std::vector<int> vec;
    double label;
    
    // Code to iterate through directory from: https://stackoverflow.com/questions/612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (fast5Dir)) != NULL) {
        Reservoir reservoir(nHashes, minHashRange, nPerBucket, total, fast5Dir, destinationDir);
        MinHash hash = MinHash((int) minHashDim, (int) nHashes*minHashPower);
        int* raw_hashes = new int[nHashes*minHashPower];
        int* rehashes = new int[nHashes];
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL) {
            vec.clear();
            squiggles.clear();
            getSquiggleVector((string) ent->d_name, squiggles, h5lsPath, h5dumpPath);
            kmerizeSquiggleSRPSliding(squiggles, vec, srpDim, srpPower, 1);
            getLabel(labelIn, label);
            hash.getHash(vec, raw_hashes);
            rehash(raw_hashes, rehashes, (int) nHashes, minHashPower);
            reservoir.add(label, rehashes, (string) ent->d_name);
        }
        closedir (dir);
        return output;
    } else {
        /* could not open directory */
        perror ("COULD NOT OPEN DIRECTORY");
        return output;
    }
}


