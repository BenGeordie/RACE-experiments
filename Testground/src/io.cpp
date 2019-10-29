#include "io.h"

#include <iostream>
#include <exception>


void kmerize(std::string sequence, std::vector<int>& vec, std::string alphabet, int k){
// Kmerize produces a vector that contains the numbers of the tokens. This means that we need to know what the numbers are. This means we need to  calculate it, and we need to know what k is. Right we need to k to kmerize in the first place. We need to know the size of alphabet used in sequence. E.g. this one is 21??
// Then go through string from idx of 0 to len(sequence) - k. Then look at substring. i 0 to <size. alphabetsize^i + position in the alphabet. Maybe I can pass a string of the alphabet. then use find to get the index and use this for numbering.

    auto alpha_size = alphabet.size();
    #pragma omp parallel for
    for(int i=0; i < sequence.length() - k + 1; ++i){
        std::string kmer;
        int tokenNum = 0;
//        std::cout << sequence.substr(i, k) << std::endl;
        try {
            kmer = sequence.substr(i, k);
        }
        catch(...) {
            std::cout << sequence << '\n';
            std::cout << "length = " << sequence.length() << '\n';
            std::cout << "i = " << i << '\n';
        }
        for(int j=0; j < k; ++j){
            try {
                tokenNum += (int) alphabet.find(kmer.at(j)) * pow(alpha_size, (double) j);
            }
            catch(...) {
                std::cout << kmer << '\n';
            }
            
//            std::cout << "mid " << j << " " << tokenNum << std::endl;
//            std::cout << "idx in alphabet: " << alphabet.find(sequence.substr(i, k).at(j)) << std::endl;
//            std::cout << "idx to the power of " << j << " = " << pow((double) alphabet.find(sequence.substr(i, k).at(j)), (double) j) << std::endl;
        }
//        std::cout << tokenNum << std::endl;
        vec.push_back(tokenNum);
    }
//    for(auto i : vec){
//        std::cout << i << std::endl;
//    }
}

void VectorFeaturesFastKraken(std::istream& in, std::vector<int>& vec, double& label, std::string& alphabet, int& k, std::string& fastWhat) {
    int interval;
    char begin;
    if (fastWhat == "fasta") {
        interval = 2;
        begin = '>';
    } else if (fastWhat == "fastq") {
        interval = 4;
        begin = '@';
    } else {
        interval = NULL;
        begin = NULL;
    }
    int current = 0;
    // Parsing the label line, extract cluster name, assign to the reference to label
    vec.clear();
    std::string labelTemp;
    std::getline(in, labelTemp);
    if (labelTemp.length() > 0 && labelTemp.at(0) == begin) {
        std::size_t start = labelTemp.find("taxid|")+6; // .find method finds index of the space before cluster name in uniref. So add 1 to get index of start of cluster name.
        labelTemp = labelTemp.substr(start, labelTemp.length());
    }
//    std::cout << labelTemp << '\n';
    
//    try {
//        labelTemp = labelTemp.substr(start, labelTemp.length());
//        std::cout << "labelTemp" << labelTemp << '\n';
//    }
//    catch (...) {
//        std::cout << labelTemp << '\n';
//    }
    std::stringstream ss(labelTemp);
    ss >> label;
    current += 1;
    
    // Parsing the sequence, pass to kmerize
    std::string sequence;
    std::string temp;
    // Get current position
    
    while(current != 0 && in) {
        //        std::cout << "peek" << in.peek() << std::endl;
        std::getline(in, temp); // each vector occupies a single line.
        if (current == 1 && temp.size() != 0) {
            sequence += temp;
        }
        current = (current + 1) % interval;
    }
    
    if (sequence.length() > 0){
        kmerize(sequence, vec, alphabet, k);
    }
}
