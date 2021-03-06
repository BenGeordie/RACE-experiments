#include "reservoir.h"
#include <ctime>
#include <cstdlib>

Reservoir::Reservoir(size_t R, size_t range, size_t threshold, size_t n_exp){
    /* parameters:
     n_exp = number of experiments
     R = number of ACE repetitions
     range = size of each ACE array
     threshold = threshold for deciding whether a vector is kept
     */
    _R = R;
    _range = range;
    _n_exp = n_exp;
    _threshold = threshold;
    
    
    _samples = new samples[n_exp];
    _counters = new counters[n_exp*R*range];
    
    // Fill each index in _samples with an empty vector
    for (size_t n = 0; n < _n_exp; n++) {
        _samples[n] = std::vector<double> {};
    }
    
    // Set each counter to 0
    for (size_t n = 0; n < _n_exp; n++) {
        for (size_t r = 0; r < _R; r++) {
            for (size_t ran = 0; ran < _range; ran++){
                _counters[(n * _R * _range) + (r * _range) + ran] = 0;
            }
        }
    }
}

Reservoir::~Reservoir(){
    delete[] _samples;
    delete[] _counters;
}

//DONE TIL HERE

void Reservoir::add(std::vector< double > vec, int *hashes){
    /*
     Input: hashes = A set of R integer hash values, one hash value for each ACE
     repetition;
     A vector that is hashed to given hash values.
     */
    #pragma omp parallel for
    for (size_t n = 0; n < _n_exp; n++) {
        double total = 0;
        for (size_t r = 0; r < _R; r++) {
            size_t bucket = hashes[n*_R+r] % _range;
            size_t bucket_index = n*_R*(_range)+r*(_range) + bucket;
            _counters[bucket_index] += 1;
            total += _counters[bucket_index];
        }
        if (total <= _threshold * _R) {
            // Assume vec's first entry is the label.
            _samples[n].push_back(vec[0]);
        }
    }
}

void Reservoir::clear(){
    memset(_samples, 0, _n_exp*sizeof(*_samples));
    memset(_counters, 0, _n_exp*_R*_range*sizeof(*_counters));
}

void Reservoir::pprint(std::ostream& out, std::string& path, bool printCounters){
    std::ofstream out_data(path);
    std::vector<size_t> clusterCountArrayForGraph;
    std::vector<size_t> elemCountArrayForGraph;
    std::set<double> setOfSequences;
    double average_clusters = 0;
    double average_elements = 0;
    for(size_t n = 0; n < _n_exp; ++n){
        setOfSequences.clear();
        for (auto entry : _samples[n]) {
            setOfSequences.insert(entry);
        }
        size_t n_clusters = setOfSequences.size();
        clusterCountArrayForGraph.push_back(n_clusters);
        average_clusters += n_clusters;
        size_t count = _samples[n].size();
        elemCountArrayForGraph.push_back(count);
        average_elements += count;
    }
    average_clusters /= _R;
    average_elements /= _R;
    
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
    
    
    // Print counters
    if (printCounters) {
        for (size_t n = 0; n < _n_exp; n++) {
            std::cout << "Exp " << n << '\n';
            auto outer = n * _R * _range;
            for (size_t r = 0; r < _R; r++) {
                auto mid = r * _range;
                for (size_t ran = 0; ran < _range; ran++) {
                    auto index = outer + mid + ran;
                    std::cout << _counters[index] << ",";
                }
                std::cout << std::endl;
                std::cout << "----------------------------------------------------";
                std::cout << std::endl;
            }
        }
    }
}
