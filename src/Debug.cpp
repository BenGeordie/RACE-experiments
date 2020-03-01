#include "Debug.h"

using namespace std;

void printPair(std::pair<std::vector<double>, int> pair){
    // Print label
    std::cout << "Label: " << pair.second << std::endl;
    
    // Print numvec
    unsigned long vecsize = pair.first.size();
    std::cout << "{";
    for(int i = 0; i < vecsize; i++){
        std::cout << pair.first[i];
        if (i == vecsize - 1){
            std::cout  << "}";
            continue;
        }
        std::cout << ",";
    }
}

void printVec(std::vector<double> vec){
    unsigned long vecsize = vec.size();
    std::cout << "{";
    for(int i = 0; i < vecsize; i++){
        std::cout << vec[i];
        if (i == vecsize - 1){
            std::cout << "}";
            continue;
        }
        std::cout << ",";
    }
}

void drawVec(std::vector< double >& vec, int width){
    unsigned long vecsize = vec.size();
    for(int i = 0; i < vecsize; ++i){
        if(i % width == 0)
            std::cout << '\n';
        if(vec[i] < 100){
            std::cout << " ";
            if(vec[i] < 10){
                std::cout << " ";
            }
        }
        std::cout << vec[i];
    }
}

double distance(std::vector< double >& vec1, std::vector< double >& vec2){
    unsigned long vecsize = vec1.size();
    double distance_sq = 0;
    for(int i = 0; i < vecsize; ++i){
        distance_sq += (vec1[i] - vec2[i])*(vec1[i] - vec2[i]);
    }
    double distance = std::sqrt(distance_sq);
    return distance;
}

std::map<int, std::vector<std::vector<double>>> sort(std::vector<std::pair<std::vector<double>, int>>& dataset, int classes){
    std::map<int, std::vector<std::vector<double>>> sorted;
    for(int i = 0; i < classes; ++i){
        std::vector<std::vector<double>> number;
        sorted[i] = number;
    }
    for(numpair pair : dataset){
        sorted[pair.second].push_back(pair.first);
    }
    return sorted;
}

std::vector<double> average(std::vector<std::vector<double>>& vectors){
    std::vector<double> average;
    unsigned long vecsize = vectors[0].size();
    for(int i = 0; i < vecsize; ++i)
        average.push_back(0);
    for(std::vector<double> vector : vectors){
        for(int i = 0; i < vecsize; ++i){
            average[i] += vector[i];
        }
    }
    unsigned long n_vectors = vectors.size();
    for(int i = 0; i < vecsize; ++i)
        average[i] = average[i]/n_vectors;
    return average;
}

double averageDist(std::map<int, std::vector<std::vector<double>>>& sorted, int number){
    double averageDist = 0;
    unsigned long n_vectors = sorted[number].size();
    int counter = 0;
    for(int i = 0; i < n_vectors; ++i){
        std::vector<int> seen;
        seen.push_back(i);
        for(int j = 0; j < n_vectors; ++j){
            std::vector<int>::iterator it = std::find(seen.begin(), seen.end(), j);
            if(it == seen.end()){
                averageDist += distance(sorted[number][i], sorted[number][j]);
                counter++;
            }
        }
    }
    averageDist = averageDist/counter;
    return averageDist;
}

double averageDistDifNum(std::map<int, std::vector<std::vector<double>>>& sorted, int number1, int number2){
    double averageDist = 0;
    unsigned long n_vectors1 = sorted[number1].size();
    unsigned long n_vectors2 = sorted[number2].size();
    int counter = 0;
    for(int i = 0; i < n_vectors1; ++i){
        for(int j = 0; j < n_vectors2; ++j){
            averageDist += distance(sorted[number1][i], sorted[number2][j]);
            counter++;
        }
    }
    averageDist = averageDist/counter;
    return averageDist;
}
