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


double JaccardDistance(const std::string &s1, const std::string &s2, int k){
	// computes the Jaccard distance between the bag-of-k-mers representation of 
	// strings s1 and s2

	int intersection_size = 0; 
	int union_size = 0;

	for (int i = 0; (i + k - 1) < s1.length(); i++){
		std::string kmer = s1.substr(i,k);
		size_t dup = s1.find(kmer);
		if (dup >= i){
			// then this is not a duplicate kmer
			union_size++;
			// std::cout<<kmer<<std::endl; 
			size_t loc = s2.find(kmer);
			if (loc != std::string::npos){
				intersection_size++; 
			}
		}
	}

	for (int i = 0; (i + k - 1) < s2.length(); i++){
		std::string kmer = s2.substr(i,k);
		size_t dup = s2.find(kmer);
		if (dup >= i){
			// then this is not a duplicate kmer in s2
			size_t loc = s1.find(kmer);
			if (loc == std::string::npos){ 
				// then this kmer is unique to s2
				// std::cout<<kmer<<std::endl;
				union_size++; 
			}
		}
	}
	// std::cout<<intersection_size<<" : "<<union_size<<std::endl; 

	double sim = (double)(intersection_size) / (double)(union_size); 
	return 1.0 - sim;
}

double HammingDistance(const std::string &s1, const std::string &s2){
	// not properly the hamming distance, because we allow 
	// different length strings. We simply assume that the last characters of the 
	// smaller string are an out-of-alphabet character 

	int minlen = std::min(s1.length(),s2.length());
	int maxlen = std::max(s1.length(),s2.length());

	double distance = 0; 
	for(int i = 0; i < minlen; i++){
		if (s1[i] != s2[i]){
			distance += 1; 
		}
	}
	distance += (maxlen - minlen); 
	return distance; 
}





double EditDistance(const std::string &s1, const std::string &s2)
// shamelessly copied from https://rosettacode.org/wiki/Levenshtein_distance#C.2B.2B
// in the interest of deadlines - I'd implement it myself but deadlines, man 
{
  const size_t m(s1.size());
  const size_t n(s2.size());
 
  if( m==0 ) return n;
  if( n==0 ) return m;
 
  size_t *costs = new size_t[n + 1];
 
  for( size_t k=0; k<=n; k++ ) costs[k] = k;
 
  size_t i = 0;
  for ( std::string::const_iterator it1 = s1.begin(); it1 != s1.end(); ++it1, ++i )
  {
    costs[0] = i+1;
    size_t corner = i;
 
    size_t j = 0;
    for ( std::string::const_iterator it2 = s2.begin(); it2 != s2.end(); ++it2, ++j )
    {
      size_t upper = costs[j+1];
      if( *it1 == *it2 )
      {
		  costs[j+1] = corner;
	  }
      else
	  {
		size_t t(upper<corner?upper:corner);
        costs[j+1] = (costs[j]<t?costs[j]:t)+1;
	  }
 
      corner = upper;
    }
  }
 
  double result = costs[n];
  delete [] costs;
 
  return result;
}
