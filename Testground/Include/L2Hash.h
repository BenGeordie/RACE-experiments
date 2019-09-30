#pragma once 

#include <limits>
#include <vector>
#include <iostream>
#include <random>
#include <cmath>


class L2Hash {
public:
	L2Hash(size_t dimensions, size_t number_of_hashes, int w);
	int* getHash(std::vector<double>& vec);
	void getHash(std::vector<double>& vec, int* hashes);
	void pprint(std::ostream &out);
	~L2Hash();

private: 
	size_t _dim; 
	size_t _nhashes;
	int _w;
	double* _C;
	double* _b;

	const unsigned int seed = 42; 

};




