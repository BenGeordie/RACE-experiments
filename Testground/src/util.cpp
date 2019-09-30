#include "util.h"

double cauchy(const std::vector<double>& x, const std::vector<double>& y, double sigma){
	double distance = L2Distance(x,y); 	
	return 1.0/(1.0 + distance*distance/(sigma*sigma) ); 
}

double inverse_multiquadric(const std::vector<double>& x, const std::vector<double>& y, double sigma){
	double distance = L2Distance(x,y); 
	return sigma/sqrt(distance*distance + sigma*sigma); 
}

double rational_quadratic(const std::vector<double>& x, const std::vector<double>& y, double sigma){
	double distance = L2Distance(x,y); 
	return 1.0 - (distance*distance)/(distance*distance + sigma); 
}

double student_t(const std::vector<double>& x, const std::vector<double>& y, int power){
	double distance = L2Distance(x,y);
	return 1.0/(1.0 + pow(distance, power) ); 
}

double L2Distance(const std::vector<double>& x, const std::vector<double>& y){
	// returns L2 distance 
	size_t dimension = x.size(); 
	double distance = 0; 
	for(size_t i = 0; i < dimension; i++)
		distance += pow(x[i] - y[i], 2.0); 
	distance = sqrt(distance); 
	return distance; 
}

double L2Norm(const std::vector<double>& x){
	// returns L2 norm
	size_t dimension = x.size(); 
	double norm = 0; 
	for(size_t i = 0; i < dimension; i++)
		norm += pow(x[i], 2.0); 
	norm = sqrt(norm); 
	return norm; 
}


void rehash(int* input_hashes, int* output_hashes, int nhashes, int values_per_set){
    #pragma omp parallel for 
    for (size_t i = 0; i < nhashes; i++){
       output_hashes[i] = MurmurHash(input_hashes + values_per_set*i, sizeof(int)*values_per_set, 42);
    }
}


