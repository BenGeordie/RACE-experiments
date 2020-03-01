#include "SequenceMinHash.h"

SequenceMinHash::SequenceMinHash(int number_of_hashes){
    _numhashes = number_of_hashes; 
}



void SequenceMinHash::getHash(size_t k, const std::string& sequence, int* hashes){
    // getHash(size_t k, std::string& sequence, int* hashes)
    // hashes had better be pre-allocated to _numhashes!! 
    // I do this because this is faster in a loop 
    #pragma omp parallel for 
    for (int n=0; n < _numhashes; n++) {

        unsigned int hashed_value;
        unsigned int minhashed_value;

        minhashed_value = std::numeric_limits<unsigned int>::max(); 
        hashes[n] = 0; 

        const char* seq = sequence.c_str(); 
        int len = sequence.length(); 

        // for each kmer in the sequence
        for (int start = 0; (start + k + 1) < len; start++){
            hashed_value = MurmurHash(seq + start, sizeof(char)*k, n); 
            // std::cout<<start<<" - "<<start +k<<" | "<<len<<" : "<<hashed_value<<std::endl; 
            // std::cout<<"Hashing: ";
            // for(int xxx = start; xxx < start+k; xxx++){
            //     std::cout<<seq[xxx]; 
            // }std::cout<<std::endl; 

            if (hashed_value < minhashed_value){
                minhashed_value = hashed_value;
                hashes[n] = MurmurHash(seq + start, sizeof(char)*k, n*3); 
                // proxy for returning the string itself
            }
        }

    }
    return; 
}
