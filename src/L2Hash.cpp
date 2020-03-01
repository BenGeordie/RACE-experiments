#include "L2Hash.h"


#include <iostream>

L2Hash::L2Hash(size_t dimensions, size_t number_of_hashes, int w)
{
	_dim = dimensions; 
	_nhashes = number_of_hashes; 
	_w = w; 

	_b = new double[_nhashes];
	_C = new double[_nhashes * _dim]; 
	
	std::mt19937 generator(seed);
	std::normal_distribution<double> gaussian(0,1.0);
	std::uniform_real_distribution<double> uniform(0.0,_w);
	
	// initialize C with a gaussian random variable iid N(0,1)
	for (size_t i = 0; i < (_nhashes*_dim); i++){
		_C[i] = gaussian(generator);
		// std::cout<<'('<<i<<','<<i/_dim<<','<<i%_dim<<") : "<<_C[i]<<std::endl; 
	}

	// initialize b with a uniform random variable 
	// std::cout<<std::endl;
	for (size_t i = 0; i < _nhashes; i++)
		_b[i] = uniform(generator);
}

L2Hash::~L2Hash(){
	delete[] _C; 
	delete[] _b; 
}

int* L2Hash::getHash(std::vector<double>& vec){
	int* hashes = new int[_nhashes];
	
	// #pragma omp parallel for
	for (size_t k = 0; k < _nhashes; k++){
		double value = 0;
		// inner products between _C[k] and vec
		for (size_t i = 0; i < _dim; i++){
			value += _C[k*_nhashes + i] * vec[i];
		}
		// if (k == 0){
		// 	std::cout<<value<<" + "<<_b[k]<<" = ";
		// 	std::cout<<value + _b[k]<<" : "<< (value + _b[k]) / _w << " -> " << floor((value + _b[k]) / _w) << std::endl; 
		// }
		value += _b[k]; 
		hashes[k] = floor(value / _w); 
	}

	// std::cout<<std::endl<<std::endl; 
	return hashes;
}

void L2Hash::getHash(std::vector<double>& vec, int* hashes){
	// #pragma omp parallel for
	for (size_t k = 0; k < _nhashes; k++){
		double value = 0;
		// inner products between _C[k] and vec
		for (size_t i = 0; i < _dim; i++){
			value += _C[k*_dim + i] * vec[i];
		}
		// if (k == 0){
		// 	std::cout<<value<<" + "<<_b[k]<<" = ";
		// 	std::cout<<value + _b[k]<<" : "<< (value + _b[k]) / _w << " -> " << floor((value + _b[k]) / _w) << std::endl; 
		// }
		// std::cout<<value<<"+"<<_b[k]<<std::endl; 
		value += _b[k]; 
		// std::cout<<value<<" : "<< value / _w << " -> " << floor(value / _w) << std::endl; 
		hashes[k] = floor(value / _w); 
	}
	// std::cout<<std::endl<<std::endl; 
}

void L2Hash::pprint(std::ostream &out){

	out<<"b: ["<<_b[0];
	for(size_t k = 1; k < _nhashes; k++)
		out<<','<<_b[k];
	out<<']'<<std::endl;

	out<<"C: ";

	for (size_t k = 0; k < _nhashes; k++){

		out<<'['<<_C[k*_dim]; 
		for (size_t i = 1; i < _dim; i++){
			out<<','<<_C[k*_dim + i]; 
		}
		out<<']'<<std::endl; 
	}

}

