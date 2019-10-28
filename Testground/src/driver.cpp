#include "driver.h"

std::vector<int> getWs(){
    std::string temp;
    std::cout << "Please enter w values separated by commas" << std::endl;
    std::getline(std::cin, temp);
    std::stringstream ss(temp);
    std::vector<int> Ws;
    int element;
    while (ss >> element)
    {
        Ws.push_back(element);
        if (ss.peek() == ',')
            ss.ignore();
    }
    // courtesy of: https://stackoverflow.com/questions/1894886/parsing-a-comma-delimited-stdstring
    return Ws;
}
std::vector<int> getPowers(){
    std::string temp;
    std::cout << "Please enter power values separated by commas" << std::endl;
    std::getline(std::cin, temp);
    std::stringstream ss(temp);
    std::vector<int> Powers;
    int element;
    while (ss >> element)
    {
        Powers.push_back(element);
        if (ss.peek() == ',')
            ss.ignore();
    }
    // courtesy of: https://stackoverflow.com/questions/1894886/parsing-a-comma-delimited-stdstring
    return Powers;
}

std::vector< std::pair<int, int> > clusteringExpMinHashStream(size_t n_hashes, int hash_power, size_t hash_range, std::ifstream& data, size_t threshold, std::string& path, size_t dimensions, size_t n_exp){
    
    // n_samples_per_bucket -> threshold
    // total -> n_exp
    
    using namespace std::chrono;
    std::vector<microseconds> latencies;
    
    // reset stream
    data.clear();
    data.seekg(0, std::ios::beg);
    
    // input buffer for vector
    std::vector<int> vec;
    // input buffer for label
    double label;
    std::vector<double> labelvec;
    
    // output  vector
    std::vector< std::pair<int, int> > output;
    
    // Initialize reservoirs (R reservoirs of hash_power buckets,
    // each containing n_samples samples, each of them a vector)
    Reservoir reservoirs(n_hashes, hash_range, threshold, n_exp);
//    reservoirs.pprint(std::cout, path);
    // Initialize hash tables
    MinHash hash = MinHash(dimensions, n_hashes*hash_power*n_exp);
    int* raw_hashes = new int[n_hashes*hash_power*n_exp];
    int* rehashes = new int[n_hashes*n_exp];
    
    int idx = 0;
    do{
        auto start = high_resolution_clock::now();
        VectorFeaturesAdjacent(data, vec, label);
//        std::cout << label << std::endl;
        if (vec.size() == 0)continue;
        
        hash.getHash(vec,raw_hashes);
        rehash(raw_hashes, rehashes, n_hashes*n_exp, hash_power);
        
        // For the moment, we are adding the labels instead of the actual vector for ease of analysis.
        labelvec.clear();
        labelvec.push_back(label); 
        labelvec.push_back((double) idx);
        //        labelvec.push_back(label);
        reservoirs.add(labelvec, rehashes);
        if(idx % 1000 == 0)
            std::cout << '.' << std::flush;
        idx++;
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        latencies.push_back(duration);
    }
    while(data);
    
    //    for(size_t r = 0; r < n_hashes; r++){
    //        int clus = reservoirs.countLabels(r); // Haven't adapted to MinHash format
    //        //        std::cout << "clus : " << clus << std::endl << std::flush;
    //        int elem = reservoirs.countElements(r);
    //        output.push_back(std::make_pair(clus, elem));
    //    }
    reservoirs.pprint(std::cout, path);
//    printList(latencies);
    return output;
}

std::vector< std::pair<int, int> > clusteringExpMinHashFastKraken(size_t n_hashes, int hash_power, size_t hash_range, std::ifstream& data, size_t threshold, std::string& path, size_t dimensions, size_t n_exp, std::string fastWhat, std::string alphabet, int k) {
    
    using namespace std::chrono;
    std::vector<microseconds> latencies;
    
    // reset stream
    data.clear();
    data.seekg(0, std::ios::beg);
    
    // input buffer for vector
    std::vector<int> vec;
    // input buffer for label
    double label;
    std::vector<double> labelvec;
    
    // output  vector
    std::vector< std::pair<int, int> > output;
    
    // Initialize reservoirs (R reservoirs of hash_power buckets,
    // each containing n_samples samples, each of them a vector)
    Reservoir reservoirs(n_hashes, hash_range, threshold, n_exp);
    //    reservoirs.pprint(std::cout, path);
    // Initialize hash tables
    MinHash hash = MinHash(dimensions, n_hashes*hash_power*n_exp);
    int* raw_hashes = new int[n_hashes*hash_power*n_exp];
    int* rehashes = new int[n_hashes*n_exp];
    
    int idx = 0;
    do{
        auto start = high_resolution_clock::now();
        VectorFeaturesFastKraken(data, vec, label, alphabet, k, fastWhat);
        //        std::cout << label << std::endl;
        if (vec.size() == 0)continue;
        
        hash.getHash(vec,raw_hashes);
        rehash(raw_hashes, rehashes, n_hashes*n_exp, hash_power);
        
        // For the moment, we are adding the labels instead of the actual vector for ease of analysis.
        labelvec.clear();
        labelvec.push_back(label);
        labelvec.push_back((double) idx);
        //        labelvec.push_back(label);
        reservoirs.add(labelvec, rehashes);
        if(idx % 1000 == 0)
            std::cout << '.' << std::flush;
        idx++;
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        latencies.push_back(duration);
    }
    while(data);
    
    //    for(size_t r = 0; r < n_hashes; r++){
    //        int clus = reservoirs.countLabels(r); // Haven't adapted to MinHash format
    //        //        std::cout << "clus : " << clus << std::endl << std::flush;
    //        int elem = reservoirs.countElements(r);
    //        output.push_back(std::make_pair(clus, elem));
    //    }
    reservoirs.pprint(std::cout, path);
    //    printList(latencies);
    return output;
}

// In the function below, the numbers printed in the console do not correspond to clusters, but instead to the read number in the fasta file (first read is labeled 0).
std::vector< std::pair<int, int> > clusteringExpMinHashStreamMurmur(size_t n_exp, size_t n_hashes, int hash_power, size_t hash_range, std::ifstream& data, std::ifstream& labelIn, size_t n_samples_per_bucket, std::vector<int>& clusters, std::string& path, size_t dimensions, size_t total, int k){
    
    using namespace std::chrono;
    std::vector<microseconds> latencies;
    
    // reset stream
    data.clear();
    data.seekg(0, std::ios::beg);
    
    labelIn.clear();
    labelIn.seekg(0, std::ios::beg);
    
    // input buffer for vector
    std::vector<int> vec;
    // input buffer for label
    double label;
    std::vector<double> labelvec;
    
    // output  vector
    std::vector< std::pair<int, int> > output;
    
    // Initialize reservoirs (R reservoirs of hash_power buckets,
    // each containing n_samples samples, each of them a vector)
    Reservoir reservoirs(n_hashes, hash_range, n_samples_per_bucket, total);
    
    // Initialize hash tables
    MinHash hash = MinHash(dimensions, n_hashes*hash_power);
    int* raw_hashes = new int[n_hashes*hash_power];
    int* rehashes = new int[n_hashes];
    
    int idx = 0;
    do{
        auto start = high_resolution_clock::now();
        VectorFeaturesFastaMurmur(data, vec, labelIn, label, k);
        //        std::cout << label << std::endl;
        if (vec.size() == 0)continue;
        
        hash.getHash(vec, raw_hashes);
        rehash(raw_hashes, rehashes, n_hashes, hash_power);
        
        
        labelvec.clear();
        labelvec.push_back(label);
        reservoirs.add(labelvec, rehashes);
        if(idx % 1000 == 0)
            std::cout << '.' << std::flush;
        idx++;
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        latencies.push_back(duration);
    }
    while(data);
    
    //    for(size_t r = 0; r < n_hashes; r++){
    //        int clus = reservoirs.countLabels(r); // Haven't adapted to MinHash format
    //        //        std::cout << "clus : " << clus << std::endl << std::flush;
    //        int elem = reservoirs.countElements(r);
    //        output.push_back(std::make_pair(clus, elem));
    //    }
    reservoirs.pprint(std::cout, path);
    printList(latencies);
    return output;
}

std::vector< std::pair<int, int> > RSExpStream(size_t n_hashes, std::ifstream& data, size_t n_samples, std::vector<int>& clusters, size_t dimensions, size_t total){
    // reset stream
    data.clear();
    data.seekg(0, std::ios::beg);
    
    // input buffer for vector
    std::vector<double> vec;
    // input buffer for label
    int label;
    std::vector<double> labelvec;
    
    // output  vector
    std::vector< std::pair<int, int> > output;
    
    // Initialize reservoirs (R reservoirs of hash_power buckets,
    // each containing n_samples samples, each of them a vector)
    Reservoir reservoirs(n_hashes, 1, n_samples, total);
    
    int idx = 0;
    do{
        VectorFeatures(data, vec, label, dimensions); //VectorFeatures?
        if (vec.size() == 0)continue;
        
        // We are adding the labels instead of the actual vector for ease of analysis.
        labelvec.clear();
        labelvec.push_back(label);
        int* fake_hash = new int[n_hashes];
        for(size_t r = 0; r < n_hashes; r++)
            fake_hash[r] = 1;
        reservoirs.add(labelvec, fake_hash);
        idx++;
    }
    while(data);
    
//    for(size_t r = 0; r < n_hashes; r++){
//        int clus = reservoirs.countLabels(r);
//        int elem = reservoirs.countElements(r);
//        output.push_back(std::make_pair(clus, elem));
//    }
//    reservoirs.pprint(std::cout, 3, true, clusters);
    return output;
}


