#include "io.h"

#include <iostream>
#include <exception>

void VectorFeatures(std::istream& in, std::vector<double>& vec, int& label, size_t& dimensions) {

    std::string line;
    std::getline(in, line); // each vector occupies a single line. 
    std::stringstream ss(line); 
    vec.clear(); 
    double element;

    while (ss >> element)
    {
        vec.push_back(element); 
        if (ss.peek() == ',')
            ss.ignore(); 
    }
    const int last_idx = (int)(vec.size() - 1);
//    std::cout<<last_idx;
    if(last_idx == dimensions){
        label = int(vec[last_idx]);
        vec.erase(vec.begin()+last_idx);
    }
    // courtesy of: https://stackoverflow.com/questions/1894886/parsing-a-comma-delimited-stdstring
}

void VectorFeaturesAdjacent(std::istream& in, std::vector<int>& vec, double& label) {
    
    std::string line;
    std::getline(in, line); // each vector occupies a single line.
    std::stringstream ss(line);
    vec.clear();
    int element;
    
    
    ss >> element;
    label = (double) element;
//    std::cout << "label" << label <<'\n';
    if(ss.peek() == ',')
        ss.ignore();
    
    while (ss >> element)
    {
        
        vec.push_back(element);
        if (ss.peek() == ',')
            ss.ignore();
    }
//    for(int ele : vec){
//        std::cout<<ele<<'\n';
//    }
}

void kmerize(std::string sequence, std::vector<int>& vec, std::string alphabet, int k){
// Kmerize produces a vector that contains the numbers of the tokens. This means that we need to know what the numbers are. This means we need to  calculate it, and we need to know what k is. Right we need to k to kmerize in the first place. We need to know the size of alphabet used in sequence. E.g. this one is 21??
// Then go through string from idx of 0 to len(sequence) - k. Then look at substring. i 0 to <size. alphabetsize^i + position in the alphabet. Maybe I can pass a string of the alphabet. then use find to get the index and use this for numbering.
    for(int i=0; i < sequence.length()-k; ++i){
        int tokenNum = 0;
//        std::cout << sequence.substr(i, k) << std::endl;
        for(int j=0; j<k; ++j){
            tokenNum += (int) alphabet.find(sequence.substr(i, k).at(j)) * pow(26, (double) j);
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

void KmerizeMurmur(std::string sequence, std::vector<int>& vec, int k) {
//    std::cout << sequence << std::endl;
    auto len = sequence.length();
    #pragma omp parallel for
    for(int i=0; i < len-k+1; ++i){
        int hash;
        if (len < k) {
            hash = MurmurHash(&sequence, (len + 1) * sizeof(char), len);
        } else {
            std::string kmer = sequence.substr(i, k);
            hash = MurmurHash(&kmer, (k+1) * sizeof(char), k);
        }
//        std::cout<<sequence.substr(i, k)<<hash<<std::endl;
        std::vector<int>::iterator it = std::find(vec.begin(), vec.end(), hash);
        if (it == vec.end()) {
            vec.push_back(hash);
        }
    }
}

/*
 Takes a vector of squiggle signal. For each k-mer, hash it with a hash function (srp?). then add hash to vec.
 */
void KmerizeSquiggleSRPSliding(std::vector<double> &squiggle, std::vector<int>& vec, int dim, int K, int L) {
    auto len = squiggle.size();
    // CHANGES I STILL NEED TO MAKE OR CONSIDER:
    // 1. REPLACE K AND L WITH SOMETHING ELSE?
    // 2. AM I FINE WITH THE VECTOR OF ARRAYS FORMAT FOR VEC, OR WOULD AN ARRAY OF VECTORS BE BETTER?
    // SAME GOES FOR THE NEXT FUNCTION/METHOD
    SignedRandomProjection *proj = new SignedRandomProjection(dim, K * L); // CHANGE – DONT KEEP AS K * L
    for(int i=0; i < len-dim+1; ++i){
        double *toHash = new double[dim];
        int *hashes;
        if (len < dim) {
            for (int j = 0; j < len; ++j) {
                toHash[j] = squiggle[i + j];
            }
            for (int j = len; j < dim; ++j) {
                toHash[j] = 0; // could be changed from 0 to what?
            }
            hashes = proj->getHash(toHash, dim);
        } else {
            for (int j = 0; j < dim; ++j) {
                toHash[j] = squiggle[i + j];
            }
            hashes = proj->getHash(toHash, dim);
        }
        //        std::cout<<sequence.substr(i, k)<<hash<<std::endl;
        // Concatenate binary values in "hashes" array into a base-10 integer.
        int hashConcat = 0;
        for (int j = 0; j < K*L; ++j) {
            hashConcat += hashes[j]*pow(2, j);
        }
        std::vector<int>::iterator it = std::find(vec.begin(), vec.end(), hashConcat);
        if (it == vec.end()) {
            vec.push_back(hashConcat);
        }
    }
}

void KmerizeSquiggleSRPPartition(std::vector<double> &squiggle, std::vector<int>& vec, int dim, int K, int L) {
    auto len = squiggle.size();
    SignedRandomProjection *proj = new SignedRandomProjection(dim, K * L); // CHANGE – DONT KEEP AS K * L
    for(int i=0; i < len-dim+1; ++dim){
        double *toHash = new double[dim];
        int *hashes;
        if (len < dim) {
            for (int j = 0; j < len; ++j) {
                toHash[j] = squiggle[i + j];
            }
            for (int j = len; j < dim; ++j) {
                toHash[j] = 0; // could be changed from 0 to what?
            }
            hashes = proj->getHash(toHash, dim);
        } else {
            for (int j = 0; j < dim; ++j) {
                toHash[j] = squiggle[i + j];
            }
            hashes = proj->getHash(toHash, dim);
        }
        //        std::cout<<sequence.substr(i, k)<<hash<<std::endl;
        // Concatenate binary values in "hashes" array into a base-10 integer.
        int hashConcat = 0;
        for (int j = 0; j < K*L; ++j) {
            hashConcat += hashes[j]*pow(2, j);
        }
        std::vector<int>::iterator it = std::find(vec.begin(), vec.end(), hashConcat);
        if (it == vec.end()) {
            vec.push_back(hashConcat);
        }
    }
}

void VectorFeaturesSquiggleCompiled(std::istream& in, std::vector<int>& vec, std::istream& labelIn, double& label, int& dim, int K, int L) {
    // Parsing the label line, extract cluster name, assign to the reference to label
    vec.clear();
    std::string temp;
    std::getline(in, temp);
    std::stringstream ssSquig(temp);
    std::string labelTemp;
    std::getline(labelIn, labelTemp); // each vector occupies a single line.
    std::stringstream ss(labelTemp);
    
    double squig;
    double element;
    
    // Pass label into "label"
    ss >> element;
    label = element;
    labelTemp.clear();
    
    // Initialize squiggle vector to pass into kmerize
    std::vector<double> squiggle;
    
    // Get current position
    
    while (ssSquig >> squig)
    {
        squiggle.push_back(squig);
        if (ssSquig.peek() == ',')
            ssSquig.ignore();
    }
    
    temp.clear();
    
    KmerizeSquiggleSRPSliding(squiggle, vec, dim, K, L);
}

void VectorFeaturesFastaMurmur(std::istream& in, std::vector<int>& vec, std::istream& labelIn, double& label, int& k) {
    
    // Parsing the label line, extract cluster name, assign to the reference to label
    vec.clear();
    std::string temp;
    std::string labelTemp;
    std::getline(labelIn, labelTemp); // each vector occupies a single line.
    std::stringstream ss(labelTemp);
    int element;
    
    ss >> element;
    label = (double) element;
//    std::cout << "label: " << label << std::endl;
    labelTemp.clear();
    
    // Parsing the sequence, pass to kmerize
    std::string sequence;
    temp.clear();
    // Get current position
    
    while(in.peek() != '>' && in){
        //        std::cout << "peek" << in.peek() << std::endl;
        std::getline(in, temp); // each vector occupies a single line.
        sequence += temp;
    }
    
    //    std::cout << sequence << std::endl;
    
    KmerizeMurmur(sequence, vec, k);
}

void VectorFeaturesFasta(std::istream& in, std::vector<int>& vec, std::string& label, std::string& alphabet, int& k) {
    
    // Parsing the label line, extract cluster name, assign to the reference to label
    vec.clear();
    std::getline(in, label);
    std::size_t start = label.find(" ")+1; // .find method finds index of the space before cluster name in uniref. So add 1 to get index of start of cluster name.
    std::size_t end = label.find(" n=");
    label = label.substr(start, end-start);
//    std::cout << '.' << std::flush;
    
    // Parsing the sequence, pass to kmerize
    std::string sequence;
    std::string temp;
    // Get current position
    
    while(in.peek() != '>' && in){
//        std::cout << "peek" << in.peek() << std::endl;
        std::getline(in, temp); // each vector occupies a single line.
        sequence += temp;
    }
    
//    std::cout << sequence << std::endl;
    
    kmerize(sequence, vec, alphabet, k);
    
}


void WriteCSVResults(std::ostream& out, size_t sketch_size, double preprocessing_time, double query_time, std::vector<double>& estimates){
    // Writes a CSV file row d
    // Assumes out starts at the start of the row (does not write leading newline)
    // Writes the following: 
    // [sketch_size in cells], [preprocessing time in seconds], [total querying time in seconds], KDE1, KDE2, KDE3, KDE4, (all doubles)
    out << sketch_size; 
    out << ',' << preprocessing_time;
    out << ',' << query_time;
    for(size_t i = 0; i < estimates.size(); i++)
        out << ',' << estimates[i];
    out<<std::endl; 
}

void printList(std::vector< double >& list){
    std::cout << '[';
    std::cout << list[0];
    for(int i = 1; i < list.size(); i++){
        std::cout << ", " << list[i];
    }
    std::cout << ']';
}

void printList(std::vector< std::chrono::microseconds >& list){
    std::cout << '[';
    std::cout << list[0].count();
    for(int i = 1; i < list.size(); i++){
        std::cout << ", " << list[i].count();
    }
    std::cout << ']';
}

std::vector<std::pair< std::vector<double>, int>> open_CSV(std::string path, int width){
    // Opens CSV file with specified row width.
    std::vector< std::pair< std::vector<double>, int > > people;
    std::ifstream data(path);
    if(!data.is_open()) std::cout << "ERROR: File Open" << '\n';
    if(data.is_open()){
        int label;
        std::pair< std::vector<double>, int > pair;
        std::vector<double> person;
        int i = 1;
        std::string temp;
        double tempD;
        while(data.good()){
            if(!(i % width)){
                std::getline(data, temp, '\n');
                std::istringstream iss (temp);
                iss >> label;
                pair = std::make_pair(person, label);
                people.push_back(pair);
                person.clear();
            } else {
                std::getline(data, temp, ',');
                //                std::cout << temp << std::endl;
                std::string::size_type sz;     // alias of size_t
                //                tempD = std::stod (temp, &sz);
                try {
                    tempD = std::stod (temp, &sz);
                } catch (std::exception& e)
                {
                    std::cout << temp << std::endl;
                    std::cout << "Standard exception: " << e.what() << std::endl;
                }
                person.push_back(tempD);
            }
            i++;
        }
        data.close();
    }
    return people;
}

std::vector<std::pair< std::vector<double>, int>> open_sparse(std::string path, int smallest_id, int largest_id){
    std::vector< std::pair< std::vector<double>, int > > dataset;
    std::ifstream data(path);
    if(!data.is_open()) std::cout << "ERROR: File Open" << '\n';
    if(data.is_open()){
        int element;
        int label;
        int prevlabel = -1; // prevent duplicate rows
        std::string temp;
        std::pair< std::vector<double>, int > pair;
        std::vector<double> vec;
        while(data.good()){
            temp.clear();
            vec.clear();
            std::getline(data, temp, '\n');
            std::stringstream ss(temp);
            temp.clear();
            std::getline(ss, temp, ',');
            std::istringstream iss (temp);
            iss >> label;
            if(label == prevlabel)
                continue;
            prevlabel = label;
            int previous = -1; // duplicate detector
            //            std::cout << "label: " << label << std::endl;
            int last = smallest_id - 1;
            while(ss){
                temp.clear();
                std::getline(ss, temp, ',');
                std::istringstream iss (temp);
                iss >> element;
                if(element - last > 1){
                    for(int i = last+1; i < element; i++){
                        vec.push_back(0.0);
                        //                        std::cout << "filler: " << i << std::endl;
                        if(i == previous){
                            std::cout << "DUPLICATE DUPLICATE DUPLICATE DUPLICATE DUPLICATE";
                        } else {
                            previous = i;
                        }
                    }
                }
                if(element != previous){
                    vec.push_back(1.0);
                    //                    std::cout << "element: " << element << std::endl;
                    previous = element;
                }
                last = element;
            }
            for(int i = last+1; i <= largest_id; i++){
                vec.push_back(0.0);
                //                std::cout << "filler 2: " << i << std::endl;
                if(i == previous){
                    std::cout << "DUPLICATE DUPLICATE DUPLICATE DUPLICATE DUPLICATE";
                } else {
                    previous = i;
                }
            }
            pair = std::make_pair(vec, label);
            dataset.push_back(pair);
        }
        data.close();
        vec.clear();
    }
    return dataset;
}

std::string path_to_root(){
    std::string file_path = __FILE__;
    std::string dir_path = file_path.substr(0, file_path.rfind("\\"));
    dir_path = dir_path.substr(0, dir_path.find("src"));
    std::cout << dir_path << std::endl;
    return dir_path;
}
