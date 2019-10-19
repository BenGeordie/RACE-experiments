#include "reservoir.h"
#include <ctime>
#include <cstdlib>

Reservoir::Reservoir(size_t R, size_t range, size_t L, std::vector<int>& clusters, size_t total){
    // parameters: R = number of ACE repetitions
    // range = size of each ACE array
    // L = size of reservoir bucket
    _R = R;
    _range = range;
    _L = L;
    
    
    // Each bucket is initialized with size L+1 because the 0th entry is used
    // to store a counter of vectors that have been passed into the bucket.
    _buckets = new reserve[_R*(_range + 1)]();
    _stats = new statistics[_R*(_range + 1)]();
    _rems = new remainders[_R]();
    for (int r = 0; r < _R; r++){
        for (int ran = 0;  ran < _range + 1; ran++){
            std::map<int, int> clus_map;
            for (int cluster : clusters)
                clus_map[cluster] = 0; // Counts for each cluster
            _stats[r * (_range+1) + ran] = clus_map;
        }
    }
    for (int r = 0; r < _R; r++){
        for (int ran = 0; ran < _range + 1; ran++){ // _range+1 => the extra bucket is to random sample to fill remaining slots if some buckets are empty. e.g. range = 5 and L = 2. If total samples less than 10 then remaining samples will be randomly sampled from this bucket.
            std::vector<double> v_0;
            v_0.push_back(0);
            _buckets[r*(_range+1) + ran].push_back(v_0);
        }
    }
    for (int r = 0; r < _R; r++){
        _rems[r] = total;
    }
}

Reservoir::~Reservoir(){
    delete[] _buckets;
    delete[] _stats;
}



void Reservoir::add(std::vector< double > vec, int *hashes){
    /*
     Input: hashes = A set of R integer hash values, one hash value for each ACE
     repetition;
     A vector that is hashed to given hash values.
     */
    #pragma omp parallel for
    srand( time( NULL ) );
    for (size_t r = 0; r < _R; r++){
        if (_rems[r] == 0) {
            continue;
        }
        size_t bucket = hashes[r] % _range;
        size_t bucket_index = r*(_range+1) + bucket;
        // If statement below is to make sure all L indices in each reservoir
        // bucket is filled
//        std::cout << vec[0] << " index: " << bucket_index << " contents: "<< _buckets[bucket_index][0][0] << std::endl;
//        std::cout <<"r: " << r << std::endl;
//        std::cout <<bucket << std::endl;
//        std::cout <<bucket_index << std::endl;
//        std::cout <<_buckets[bucket_index][0][0] << std::endl;
        if (_buckets[bucket_index][0][0] < _L){
            _buckets[bucket_index].push_back(vec);
            _rems[r] --;
        }
        // Else condition below is for when there are more than L items; some
        // have to be discarded.
        else{
            /*
             First +1 is because we're adding a new vector;
             P(keep a vector) = 1/(total no. of vectors).
             total no. of vectors = no. of vectors already passed into bucket + 1.
             
             Second +1 is again an offset since 0th entry countains counter.
             */
            int index = rand() % ((int)_buckets[bucket_index][0][0]+1) + 1;
            if(index < (_L+1)){
                _buckets[bucket_index][index] = vec;
            }
//            else {
//                if (_buckets[r*(_range+1)+_range][0][0] < _rems[r]){
//                    _buckets[r*(_range+1)+_range].push_back(vec);
//                }else{
//                    int remindex = rand() % ((int)_buckets[r*(_range+1)+_range][0][0]+1) + 1;
//                    if(remindex < (_rems[r]+1)){
////                        std::cout << remindex << std::endl;
////                        std::cout << _buckets[r*(_range+1)+_range].size() << std::endl;
//                        _buckets[r*(_range+1)+_range][remindex] = vec;
//                    }
//                }
//                _buckets[r*(_range+1)+_range][0][0] ++;
//                // +1 offset to _rems[r] because first entry is counter of elements directed to the bucket.
//                for(size_t i = _buckets[r*(_range+1)+_range].size() - 1; i >= _rems[r] + 1; --i){
//                    _buckets[r*(_range+1)+_range].erase(_buckets[r*(_range+1)+_range].begin()+i);
//                }
//            }
        }
        
        _buckets[bucket_index][0][0]++;
//        std::cout << "label = " << vec[0] << std::endl;
        _stats[bucket_index][int(vec[0])] = _stats[bucket_index][int(vec[0])] + 1;
    }
}

void Reservoir::clear(){
    memset(_buckets, 0, _R*(_range+1)*sizeof(*_buckets));
}

void Reservoir::pprint(std::ostream& out, int width, bool format, std::vector<int>& clusters){
    // For each r in range _R
    for (size_t r = 0; r < _R; r++){
     std::cout << "r = " << r << std::endl;
     // For each bucket (displayed as a row)
     for (size_t i = 0; i < (_range+1); i++){
         if (_buckets[r*(_range+1) + i][0][0] > 0){
             if (format)
                 std::cout << std::string((width+1)*_L + 1, '-') << std::endl;
             // Some buckets have less than L elements. Find bucket size first.
             size_t bucketsize = _buckets[r*(_range+1) + i].size();
             // bucketsize includes the first element, which is just a counter. So below we use bucketsize-1.
             for (size_t l = 0; l < bucketsize - 1; l++){
                 std::cout << '|' << std::setw(width);
                 printVec(_buckets[r*(_range+1) + i][l+1]);
                 
             }
             std::cout << "Res " << i << " size = " << _buckets[r*(_range+1) + i][0][0] << std::endl;
             for(int n : clusters){
                 std::cout << n << " = " << _stats[r*(_range+1) + i][n]*100/_buckets[r*(_range+1) + i][0][0] << "%" << std::endl;
             }
         }
         
     }
     if (format)
         out << std::string((width+1)*_L + 1, '-') << std::endl;
     out << std::endl;
    }
}

void Reservoir::pprint(std::ostream& out, std::string& path){
    std::ofstream out_data(path);
    std::vector<int> clusterCountArrayForGraph;
    std::vector<int> elemCountArrayForGraph;
    std::set<std::vector<double>> setOfSequences;
    double average_clusters = 0;
    double average_elements = 0;
    for(size_t r = 0; r < _R; ++r){
        setOfSequences.clear();
//        std::cout << "R = " << r <<'\n';
        for(int ran = 0; ran < _range + 1; ++ran){
//            std::cout << "\n-------------------------\n";
//            std::cout << ran << '\n';
            for(int i=1; i < _buckets[r*(_range+1) + ran].size(); ++i){
//                std::cout << _buckets[r*(_range+1) + ran][i][0] << ", ";
                setOfSequences.insert(_buckets[r*(_range+1) + ran][i]);
                std::copy(_buckets[r*(_range+1) + ran][i].begin(), _buckets[r*(_range+1) + ran][i].end(),
                          std::ostream_iterator<double>(std::cout, ", "));
            }
        }
//        std::cout << "\n-------------------------\n";
        int n_clusters = setOfSequences.size();
        clusterCountArrayForGraph.push_back(n_clusters);
//        out_data << "number of clusters = " << n_clusters << '\n';
//        std::cout << "number of clusters = " << n_clusters << '\n';
        average_clusters += n_clusters;
        int count = countElements(r);
        elemCountArrayForGraph.push_back(count);
//        out_data << "number of elements = " << count << '\n';
//        std::cout << "number of elements = " << count << '\n';
        average_elements += count;
    }
    average_clusters /= _R;
    average_elements /= _R;
//    out_data << "average number of clusters = " << average_clusters << '\n';
//    std::cout << "average number of clusters = " << average_clusters << '\n';
//    out_data << "number of elements = " << average_elements << '\n';
//    std::cout << "number of elements = " << average_elements << '\n';
    
    // Print clusterCountArrayForGraph
    std::cout << "clusterCountArrayForGraph = [";
    std::cout << clusterCountArrayForGraph[0];
    for(size_t i = 1; i < clusterCountArrayForGraph.size(); ++i){
        std::cout << ", " << clusterCountArrayForGraph[i];
    }
    std::cout << "]" << std::endl;
    
    // Print elemCountArrayForGraph
    std::cout << "elemCountArrayForGraph = [";
    std::cout << elemCountArrayForGraph[0];
    for(size_t i = 1; i < elemCountArrayForGraph.size(); ++i){
        std::cout << ", " << elemCountArrayForGraph[i];
    }
    std::cout << "]" << std::endl;
}

int Reservoir::count(){
    int clusters = 0;
    for (size_t r = 0; r < _R; r++){
        // For each bucket (displayed as a row)
        for (size_t i = 0; i < (_range+1); i++){
            if (_buckets[r*(_range+1) + i][0][0] > 10 /*&& _buckets[r*_range + i][0][0] < 1000*/){
                clusters++;
            }
        }
    }
    return clusters;
}

int Reservoir::countLabels(size_t& r){
//    std::cout << "THE R VALUE IS : " << r << std::endl << std::flush;
    std::vector<int> labels;
    int n;
    // For each bucket (displayed as a row)
    for (size_t i = 0; i < (_range+1); i++){
        size_t lim = _L;
        if(i == _range){
            lim = _rems[r];
        }
        for (size_t j = 0; j < _buckets[r*(_range+1) + i][0][0] && j < lim; j++){ //_buckets[r*_range + i][0][0] shows number of elements in bucket.
            n = _buckets[r * (_range+1) + i][j+1][0]; // +1 for offset
            std::vector<int>::iterator it = std::find(labels.begin(), labels.end(), n);
            if(it == labels.end())
//                std::cout << "n : " << n << std::endl << std::flush;
                labels.push_back(n);
        }
    }
    return labels.size();
}

int Reservoir::countElements(size_t& r){
    int elements = 0;
    // For each bucket (displayed as a row)
    for (size_t i = 0; i < (_range+1); i++){
        elements += _buckets[r*(_range+1) + i].size() - 1;
    }
    return elements;
}
