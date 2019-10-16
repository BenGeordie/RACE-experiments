#include "reservoir.h"
#include <ctime>
#include <cstdlib>

// rewrite reserve to counters.

Reservoir::Reservoir(size_t R, size_t range, size_t L, size_t total, const char* oldDir, std::string& newDir){
    // parameters: R = number of ACE repetitions
    // range = size of each ACE array
    // L = size of reservoir bucket
    _R = R;
    _range = range;
    _L = L;
    _total = total;
    _old = oldDir;
    _new = newDir;
    
    
    // Each bucket is initialized with size L+1 because the 0th entry is used
    // to store a counter of vectors that have been passed into the bucket.
    _race = new int[_R*_range];
    for (int r = 0; r < _R; r++){
        for (int ran = 0; ran < _range; ran++){
            _race[r*(_range) + ran] = 0;
        }
    }
}

Reservoir::~Reservoir(){
    delete[] _race;
    _sample.clear();
    _labelSet.clear();
}


void Reservoir::add(int label, int *hashes, std::string fileName){
    /*
     Input: hashes = A set of R integer hash values, one hash value for each ACE
     repetition;
     A vector that is hashed to given hash values.
     */
    if (_sample.size() >= _total) {
        return;
    }
    double averageCount = 0;
    #pragma omp parallel for
    for (size_t r = 0; r < _R; r++){
        size_t bucket = hashes[r] % _range;
        size_t bucketIndex = r*(_range) + bucket;
        _race[bucketIndex]++;
        averageCount += (double) _race[bucketIndex];
    }
    if (averageCount/_R < _L) {
        _sample.push_back(make_pair(label, fileName));
        _labelSet.insert(label);
        
        std::string cmdString = "cp " + _old + fileName + " " + _new + fileName;
        const char* cmd = cmdString.c_str();
        system(cmd);
    }
}

void Reservoir::clear(){
    memset(_race, 0, _R*(_range)*sizeof(*_race));
}

// DONE EDITING TIL HERE

void Reservoir::pprint(std::ostream& outPath){
    // For each r in range _R
    std::cout << "RACE" << std::endl;
    for (size_t r = 0; r < _R; r++){
        // std::cout << "r = " << r << std::endl;
        std::cout << r << ": [" << _race[r*_range];
        for (size_t i = 0; i < (_range); i++){
            std::cout << ", " << _race[r*_range+i];
        }
        std::cout << "]" << std::endl;
    }
    // STILL NEEDA PRINT ACTUAL VECTORS, OR FOR EACH SAMPLE, STORE THE COPIES OF ORIGINAL FILES IN A SEPARATE FOLDER?
    
}

// WHY DO I HAVE ANOTHER PRINT FUNCTION THIS IS STUPID AND WHY IS IT SO FRIKKIN LONG
//void Reservoir::pprint(std::ostream& out, std::string& path){
//    std::ofstream out_data(path);
//    std::vector<int> clusterCountArrayForGraph;
//    std::vector<int> elemCountArrayForGraph;
//    std::set<std::vector<double>> setOfSequences;
//    double average_clusters = 0;
//    double average_elements = 0;
//    for(size_t r = 0; r < _R; ++r){
//        setOfSequences.clear();
////        std::cout << "R = " << r <<'\n';
//        for(int ran = 0; ran < _range + 1; ++ran){
////            std::cout << "\n-------------------------\n";
////            std::cout << ran << '\n';
//            for(int i=1; i < _sample[r*(_range+1) + ran].size(); ++i){
////                std::cout << _buckets[r*(_range+1) + ran][i][0] << ", ";
//                setOfSequences.insert(_buckets[r*(_range+1) + ran][i]);
//            }
//        }
////        std::cout << "\n-------------------------\n";
//        int n_clusters = setOfSequences.size();
//        clusterCountArrayForGraph.push_back(n_clusters);
////        out_data << "number of clusters = " << n_clusters << '\n';
////        std::cout << "number of clusters = " << n_clusters << '\n';
//        average_clusters += n_clusters;
//        int count = countElements(r);
//        elemCountArrayForGraph.push_back(count);
////        out_data << "number of elements = " << count << '\n';
////        std::cout << "number of elements = " << count << '\n';
//        average_elements += count;
//    }
//    average_clusters /= _R;
//    average_elements /= _R;
////    out_data << "average number of clusters = " << average_clusters << '\n';
////    std::cout << "average number of clusters = " << average_clusters << '\n';
////    out_data << "number of elements = " << average_elements << '\n';
////    std::cout << "number of elements = " << average_elements << '\n';
//
//    // Print clusterCountArrayForGraph
//    std::cout << "clusterCountArrayForGraph = [";
//    std::cout << clusterCountArrayForGraph[0];
//    for(size_t i = 1; i < clusterCountArrayForGraph.size(); ++i){
//        std::cout << ", " << clusterCountArrayForGraph[i];
//    }
//    std::cout << "]" << std::endl;
//
//    // Print elemCountArrayForGraph
//    std::cout << "elemCountArrayForGraph = [";
//    std::cout << elemCountArrayForGraph[0];
//    for(size_t i = 1; i < elemCountArrayForGraph.size(); ++i){
//        std::cout << ", " << elemCountArrayForGraph[i];
//    }
//    std::cout << "]" << std::endl;
//}

// MODIFY TO ONLY LOOK THROUGH ONE VECTOR
int Reservoir::countLabels(){
    return (int) _labelSet.size();
}

//int Reservoir::countElements(size_t& r){
//    int elements = 0;
//    // For each bucket (displayed as a row)
//    for (size_t i = 0; i < (_range+1); i++){
//        elements += _sample[r*(_range+1) + i].size() - 1;
//    }
//    return elements;
//}
