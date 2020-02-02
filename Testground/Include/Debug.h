#include <map>
#include <algorithm>
#include "driver.h"
// #include "opencv2/opencv.hpp"
// #include <opencv2/core/core.hpp>
// #include <opencv2/highgui/highgui.hpp>

typedef std::pair<std::vector<double>, int> numpair;

void printPair(std::pair<std::vector<double>, int> pair);

void printVec(std::vector<double> vec);

void drawVec(std::vector< double >& vec, int width);

double distance(std::vector< double >& vec1, std::vector< double >& vec2);

std::map<int, std::vector<std::vector<double>>> sort(std::vector<std::pair<std::vector<double>, int>>& dataset, int classes);

std::vector<double> average(std::vector<std::vector<double>>& vectors);

double averageDist(std::map<int, std::vector<std::vector<double>>>& sorted, int number);

double averageDistDifNum(std::map<int, std::vector<std::vector<double>>>& sorted, int number1, int number2);
