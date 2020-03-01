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
        std::clog<<"KNNBuffer <n_samples> <distance_id> <k> <data> <output>"<<std::endl;
        std::clog<<"n_samples: (integer number of samples stored by the KNNBuffer)"<<std::endl;
        std::clog<<"distance_id: (integer that species the string distance)"<<std::endl;
        std::clog<<"Supported IDs:0\t Jaccard distance\n             :1\t Hamming distance";
        std::clog<<"\n             :2\t Levenshtein distance (edit distance)"<<std::endl;
        std::clog<<"k: (integer that decides how long the k-mers are)"<<std::endl;
        // k is only used for the Jaccard distance but has to be specified because its simpler for 
        // me to implement and cruncy cruncy time crunchy AHHHHHH
        std::clog<<std::endl<<"Example: "<<std::endl; 
        std::clog<<std::endl<<"KNNBuffer 200 1 6 input.fastq output.txt"<<std::endl; 
        return -1; 
    }

    int n_samples = std::stoi(argv[1]);
    int distance_id = std::stoi(argv[2]);
    int kmer_k = std::stoi(argv[3]);

    std::ifstream datastream(argv[4]);
    std::ofstream samplestream(argv[5]);

    // reset stream
    datastream.clear();
    datastream.seekg(0, std::ios::beg);

    std::string sequence;
    int label;

    std::vector< std::string > sequences;
    std::vector< int > labels;


    /*
    std::string test1 = "TTAGGACTTAG";//ATCATGGCTCAGATTGAACGCTGGCGGCATGCC"; 
    std::string test2 = "TTAGGGCGGCATGCTTAG";//TAACACATGCAAGTCGAACGAGATCTTCGGATC"; 
    std::cout<<EditDistance(test1, test2)<<std::endl; 
    std::cout<<JaccardDistance(test1, test2,4)<<std::endl; 
    std::cout<<HammingDistance(test1, test2)<<std::endl; 
    */

    // ridiculously ad-hoc implementation of greedy kcenter sampling

    // First we store the first n_samples sequences in a buffer. At the end of this process, 
    // we find the minimum pairwise distance between sequences in the buffer. 
    // Let D = minimum pairwise distance 
    // and let i1 and i2 be the indices of the pair that have the minimum pairwise distance

    // When a new element X is added, we compute the distance to each sequence already in the buffer. 
    // We compute the distance d(X,Xi) for each element Xi in the buffer. 
    // We remember the smallest 2 distances and the corresponding indices i and j. 

    // If the smallest distance D' > D, then
    // 1. Consider the two sequences at i1 and i2 (the previous smallest pair)
    // 2. Find out which sequence has smaller distance to X
    // 3. Delete that sequence and replace it with X
    // We can't just update D = D' because we might have deleted the element for which d(X,y) = D'
    // 4. if we replaced the index of the smallest distance, D is now the second-smallest distance
    // 5. if not, then D is now the smallest distance 


    int i = 0;
    int buffer_idx1 = 0; 
    int buffer_idx2 = 0;
    double buffer_D = 0;
    do{
        bool success = KrakenSequenceFeatures(datastream, sequence, label, "fastq");
        if (!success) continue;


        if (i < n_samples){
            // First we store the first n_samples sequences in a buffer
            sequences.push_back(sequence);
            labels.push_back(label);
        } else if (i == n_samples){
            // std::cout<<"finding distances"<<std::endl; 
            // we find the minimum pairwise distance between sequences in the buffer
            double min_d = std::numeric_limits<double>::max();
            for(int j = 0; j < n_samples; j++){
                for(int k = j+1; k < n_samples; k++){
                    double d = knnbuffer_distance(sequences[j],sequences[k],kmer_k,distance_id); 
                    // std::cout<<d<<","; 
                    if (d < min_d){
                        min_d = d;
                        buffer_idx1 = j;
                        buffer_idx2 = k;
                    }
                }
            }
            buffer_D = min_d;
            // std::cout<<std::endl<<"Min d = "<<buffer_D<<std::endl; 
            std::cout<<std::flush; 
        } else {
            // here is where it gets interesting

            // When a new element X is added, we compute the distance to each sequence already in the buffer. 
            // We remember the smallest 2 distances and the corresponding indices i and j. 
            double min_d1 = std::numeric_limits<double>::max();
            double min_d2 = std::numeric_limits<double>::max();
            int sample_idx1 = 0; 
            int sample_idx2 = 0; 

            for (int j = 0; j < n_samples; j++){
                double d = knnbuffer_distance(sequence,sequences[j],kmer_k,distance_id); 
                if (d < min_d1){
                    min_d2 = min_d1; 
                    sample_idx2 = sample_idx1; 
                    min_d1 = d;
                    sample_idx1 = j; 
                } else if (d < min_d2){
                    min_d2 = d; 
                    sample_idx2 = j; 
                }
            }

            // std::cout<<std::endl<<"Min d for sample = "<<min_d1<<" : "<<min_d2<<" | buffer_D = "<<buffer_D<<std::endl; 

            // If the smallest distance D' > D, then
            if (min_d1 > buffer_D){
                // 1. Consider the two sequences at i1 and i2 (the previous smallest pair)
                // 2. Find out which sequence has smaller distance to X
                double d1 = knnbuffer_distance(sequence, sequences[buffer_idx1], kmer_k, distance_id); 
                double d2 = knnbuffer_distance(sequence, sequences[buffer_idx2], kmer_k, distance_id); 

                // 3. Delete that sequence and replace it with X
                if (d1 < d2){
                    sequences[buffer_idx1] = sequence; 
                    labels[buffer_idx1] = label; 
                    // 4. if we replaced the index of the smallest distance, D is now the second-smallest distance
                    // 5. if not, then D is now the smallest distance 
                    if (sample_idx1 != buffer_idx1){
                        buffer_D = min_d1; 
                    } else {
                        buffer_D = min_d2; 
                    }
                } else {
                    sequences[buffer_idx2] = sequence; 
                    labels[buffer_idx2] = label; 
                    // 4. if we replaced the index of the smallest distance, D is now the second-smallest distance
                    // 5. if not, then D is now the smallest distance 
                    if (sample_idx2 != buffer_idx1){
                        buffer_D = min_d1; 
                    } else {
                        buffer_D = min_d2; 
                    }
                }
            }
        }
        i++;
    }
    while(datastream);
    for(int i = 0; i < n_samples; i++){
        samplestream<<labels[i];
        samplestream<<'\t';
        samplestream<<sequences[i]; 
        samplestream<<std::endl; 
    }

}
