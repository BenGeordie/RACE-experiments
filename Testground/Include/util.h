#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

#include <math.h>

#include "MurmurHash.h"

// mathematic and hash utilities, collected in one place
// for instance, we define kernels in here

double cauchy(const std::vector<double>& x, const std::vector<double>& y, double sigma); 
double student_t(const std::vector<double>& x, const std::vector<double>& y, int power); 
double inverse_multiquadric(const std::vector<double>& x, const std::vector<double>& y, double sigma); 
double rational_quadratic(const std::vector<double>& x, const std::vector<double>& y, double sigma); 

double L2Distance(const std::vector<double>& x, const std::vector<double>& y); 

double L2Norm(const std::vector<double>& x); 

void rehash(int* input_hashes, int* output_hashes, int nhashes, int values_per_set); 

