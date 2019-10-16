#include "io.h"

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}


/*
 Takes a vector of squiggle signal. For each k-mer, hash it with a hash function (srp?). then add hash to vec.
 */
void kmerizeSquiggleSRPSliding(std::vector<double> &squiggle, std::vector<int>& vec, int dim, int K, int L) {
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

void kmerizeSquiggleSRPPartition(std::vector<double> &squiggle, std::vector<int>& vec, int dim, int K, int L) {
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

void vectorFeaturesSquiggleCompiled(std::istream& in, std::vector<int>& vec, std::istream& labelIn, double& label, int& dim, int K, int L) {
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
    
    kmerizeSquiggleSRPSliding(squiggle, vec, dim, K, L);
}

void parseSignalDirectory(std::istringstream& cliOutput, std::string& path) {
    while (std::getline(cliOutput, path)) {
        if (path.find("Signal") != std::string::npos) {
            path = path.substr(0, path.find("Signal") + 6);
            break;
        }
    }
}

void parseSquiggle(std::istringstream& cliOutput, std::vector<double>& squiggles, std::string& temp) {
    while (std::getline(cliOutput, temp)) {
        // temp.size() > 3 is a condition because the empty line after the last squiggle signal is 3 characters long (for the tabs?)
        if (temp.find("DATA") == std::string::npos && temp.find("HDF5") == std::string::npos && temp.find("}") == std::string::npos && temp.size() > 3) {
            // std::cout << temp << temp.size() << std::endl;
            // auto commaLoc = temp.find(",");
            std::stringstream ss(temp);
            double squiggle = 0;
            ss >> squiggle;
            squiggles.push_back(squiggle);
        }
    }
}

void getSquiggleVector(std::string fileName, std::vector<double>& squiggles, std::string& h5lsPath, std::string& h5dumpPath) {
    squiggles.clear();
    std::string path;
    
    // First use h5ls to find what the signal directory path. Use parser to isolate this path from the other contents of the fast5 file.
    std::string cmdString = h5lsPath + " -r " + fileName;
    const char* cmd = cmdString.c_str();
    std::istringstream cliOutputDir(exec(cmd));
    parseSignalDirectory(cliOutputDir, path);
    cliOutputDir.clear();
    
    // Then use h5dump to give squiggle signals from the directory path. Use parser to isolate the squiggle signals and put each number as a double in a vector.
    cmdString = h5dumpPath + " -w 1 -y -d " + path + " " + fileName;
    cmd = cmdString.c_str();
    std::istringstream cliOutputPath(exec(cmd));
    parseSquiggle(cliOutputPath, squiggles, path);
    cliOutputPath.clear();
}

void getLabel(std::istream& labelIn, double& label) {
    std::string temp;
    std::getline(labelIn, temp); // each vector occupies a single line.
    std::stringstream ss(temp);
    double element;
    ss >> element;
    label = element;
    temp.clear();
}

void writeCSVResults(std::ostream& out, size_t sketch_size, double preprocessing_time, double query_time, std::vector<double>& estimates){
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
    if (list.size() > 0) {
        std::cout << list[0];
        for(int i = 1; i < list.size(); i++){
            std::cout << ", " << list[i];
        }
    }
    std::cout << ']';
}

void printList(std::vector< int >& list){
    std::cout << '[';
    if (list.size() > 0) {
        std::cout << list[0];
        for(int i = 1; i < list.size(); i++){
            std::cout << ", " << list[i];
        }
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

std::string path_to_root(){
    std::string file_path = __FILE__;
    std::string dir_path = file_path.substr(0, file_path.rfind("\\"));
    dir_path = dir_path.substr(0, dir_path.find("src"));
    std::cout << dir_path << std::endl;
    return dir_path;
}
