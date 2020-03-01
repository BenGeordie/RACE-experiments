#include "io.h"
#include "SequenceMinHash.h"
#include "RACE.h"
#include "util.h"

#include <chrono>
#include <limits>
#include <algorithm>

/*
Writes throughput timing information to clog and writes the output samples to cout

Details of KNN buffer algorithm: 

We find the most diverse sequences through exhaustive search through the buffer

First we store the first n_samples sequences in a buffer. At the end of this process, 
we find the minimum pairwise distance between sequences in the buffer. 
Let D = minimum pairwise distance 
and let i1 and i2 be the indices of the pair that have the minimum pairwise distance

When a new element X is added, we compute the distance to each sequence already in the buffer. 
We compute the distance d(X,Xi) for each element Xi in the buffer. 
We remember the smallest 2 distances and the corresponding indices i and j. 

If the smallest distance D' > D, then
1. Consider the two sequences at i1 and i2 (the previous smallest pair)
2. Find out which sequence has smaller distance to X
3. Delete that sequence and replace it with X
We can't just update D = D' because we might have deleted the element for which d(X,y) = D'
4. if we replaced the index of the smallest distance, D is now the second-smallest distance
5. if not, then D is now the smallest distance 


Note: The distance can be the Jaccard distance, edit distance or string hamming distance
You can pick which distance metric you want to use. 
*/

// ugly ugly ugly its so ugly to declare a function here 
// but I only have like two hours to write this thing
double knnbuffer_distance(const std::string &s1, const std::string &s2, int k, int distance_id){
    double d = 0; 
    switch(distance_id){
        case 0: d = JaccardDistance(s1,s2,k);
        break; 
        case 1: d = HammingDistance(s1,s2); 
        break; 
        case 2: d = EditDistance(s1,s2); 
        break;
    }
    return d; 
}


int main(int argc, char **argv){
    if (argc < 6){
        std::clog<<"Usage: "<<std::endl;
        std::clog<<"Coresets <window_size> <samples_per_window> <distance_id> <k> <data> <output>"<<std::endl;
        std::clog<<"window_size: (size of each coreset window)"<<std::endl;
        std::clog<<"samples_per_window: (number of samples stored in each window)"<<std::endl;
        std::clog<<"Supported IDs:0\t Jaccard distance\n             :1\t Hamming distance";
        std::clog<<"\n             :2\t Levenshtein distance (edit distance)"<<std::endl;
        std::clog<<"k: (integer that decides how long the k-mers are)"<<std::endl;
        // k is only used for the Jaccard distance but has to be specified because its simpler for 
        // me to implement and cruncy cruncy time crunchy AHHHHHH
        std::clog<<std::endl<<"Example: "<<std::endl; 
        std::clog<<std::endl<<"Coresets 200 10 1 6 input.fastq output.txt"<<std::endl; 
        return -1; 
    }

    std::clog<<"Reading parameters"<<std::endl; 

    int window_size = std::stoi(argv[1]);
    int samples_per_window = std::stoi(argv[2]);
    int distance_id = std::stoi(argv[3]);
    int kmer_k = std::stoi(argv[4]);
    std::ifstream datastream(argv[5]);
    std::ofstream samplestream(argv[6]);

    // reset stream
    datastream.clear();
    datastream.seekg(0, std::ios::beg);

    std::string sequence;
    int label;

    std::vector< std::string > window_sequences(window_size);
    std::vector< int > window_labels(window_size);
    double* distances = new double[window_size]();

    // ridiculously ad-hoc implementation of greedy kcenter sampling
    int i = 0;
    auto all_start = std::chrono::high_resolution_clock::now();
    do{
        auto start = std::chrono::high_resolution_clock::now();
        bool success = KrakenSequenceFeatures(datastream, sequence, label, "fastq");
        if (!success) continue;

        int index = i%window_size;

        if ((index == 0) && (i > 0)){
            // then we have to greedy max min (GMM) the window
            std::vector< int > center_indices;
            center_indices.push_back(0); 
            for (int c = 0; c < samples_per_window-1; c++){
                std::clog<<'+'<<std::flush;
                // for each center we need to find, first reset distances
                for (int j = 0; j < window_size; j++){distances[j]= std::numeric_limits<double>::max();}
                // for each center we have already found update min distances
                for (int center_idx : center_indices){

                    #pragma omp parallel for
                    for (int j = 0; j < window_size; j++){
                        double d = knnbuffer_distance(window_sequences[center_idx],window_sequences[j],kmer_k,distance_id);
                        if (d < distances[j]){distances[j] = d;}
                    }
                }
                // now find the maximum of the minimums
                double maxd = 0; 
                int maxid = 0; 
                for (int j = 0; j < window_size; j++){
                    if (distances[j] > maxd){
                        maxd = distances[j]; 
                        maxid = j; 
                    }
                } // and add it to the center_indices
                center_indices.push_back(maxid);
            }
            // now print out the kcenters we found 
            for (int center_idx : center_indices){
                samplestream<<window_labels[center_idx]<<'\t'<<window_sequences[center_idx]<<std::endl; 
            }
            std::clog<<std::endl<<std::flush;

        } else {
            window_labels[index] = label; 
            window_sequences[index] = sequence; 
        }
        auto finish = std::chrono::high_resolution_clock::now();
        std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count() << ",";
   
        i++;
        if (i%1000==0)
            std::clog<<'.'<<std::flush;
    }
    while(datastream);
    auto all_finish = std::chrono::high_resolution_clock::now();
    std::clog <<"Processed "<<i<<" sequences in "<< std::chrono::duration_cast<std::chrono::seconds>(all_finish-all_start).count()<<" seconds";
}
