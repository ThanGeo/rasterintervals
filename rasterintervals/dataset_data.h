#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <vector>
#include <math.h>
#include <cmath>
#include <ctime>
#include <map>


using namespace std;

#ifndef DATASET_H
#define DATASET_H

string getBinaryGeometryFilename(string argument){
	return ("datafiles/" + argument + "_fixed_binary.dat");
}

string getIntervalBinaryFilename(string argument){
	return ("rasterintervals/interval_data/"+argument+"_binary_interval_data.dat");
}

string getOffsetMap(string argument){
	return("datafiles/" + argument + "_offset_map.dat");
}

#endif
