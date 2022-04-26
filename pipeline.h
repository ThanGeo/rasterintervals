#ifndef PIPE_H
#define PIPE_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <queue>
#include <fstream>
#include <future>
#include <unordered_set>
#include <unordered_map>
#include <chrono>
#include "omp.h"
#include <functional>
#include <array>

#include "./rasterintervals/raster_intervals.h"

using namespace std;

int INTERMEDIATE_FILTER = 0;
int REFINEMENT = 0;
int CALCULATE_INTERVALS = 0;

int TOTAL_RESULTS = 0;

Dataset rasterIntervalsR;
Dataset rasterIntervalsS;

//geometry filenames
string geometryFileNameR;
string geometryFileNameS;
//offset maps for binary geometries
unordered_map<uint,unsigned long> offsetMapR;
unordered_map<uint,unsigned long> offsetMapS;
//binary geometry files
ifstream finR;
ifstream finS;

string argument1, argument2;

string result_filename = "results_pairs.csv";

double intermediateFilterTime = 0;
double refinementTime = 0;

int postMBRCandidates = 0;
int accepted = 0;
int rejected = 0;
int refinementCandidates = 0;
int acceptedAfterRefinement = 0;

clock_t timer;

/* saves a result pair on disk */
void saveResultPair(uint &idA, uint &idB){
        ofstream fout(result_filename, fstream::out | ios_base::binary | fstream::app);

        fout << idA << " " << idB << endl; 

        fout.close();
}

void setUniversalCoordinates(ifstream &fin){
        int totalPolygonCount, vertexCount;
        double x,y;
        double polxMin,polyMin,polxMax,polyMax;
        uint recID;

        //read total polygon count from binary geometry file
        fin.read((char*) &totalPolygonCount, sizeof(int));

        uint lineCounter = 0;
        while(lineCounter < totalPolygonCount){
                //---BUILD POLYGON---

                //read pol id
                fin.read((char*) &recID, sizeof(int));
                Polygon pol(recID);

                polxMin = numeric_limits<int>::max();
                polyMin = numeric_limits<int>::max();
                polxMax = -numeric_limits<int>::max();
                polyMax = -numeric_limits<int>::max();

                //read vertex count for polygon
                fin.read((char*) &vertexCount, sizeof(int));            
                pol.vertices.reserve(vertexCount);
                //read #vertexCount points
                for(int i=0; i<vertexCount; i++){
                        //read x, y
                        fin.read((char*) &x, sizeof(double));
                        fin.read((char*) &y, sizeof(double));

                        pol.vertices.emplace_back(x, y);

                        polxMin = min(x, polxMin);
                        polyMin = min(y, polyMin);
                        polxMax = max(x, polxMax);
                        polyMax = max(y, polyMax);
                }
                //set universal coordinates
                universalMinX = min(universalMinX, polxMin);
                universalMinY = min(universalMinY, polyMin);
                universalMaxX = max(universalMaxX, polxMax);
                universalMaxY = max(universalMaxY, polyMax);

                lineCounter++;
        }
}

/* initializes certain variables and files for the process */
void initialize(string &arg1, string &arg2){
        argument1 = arg1;
        argument2 = arg2;

        // geometries files and offset maps
        geometryFileNameR = getBinaryGeometryFilename(argument1);
        geometryFileNameS = getBinaryGeometryFilename(argument2);
        offsetMapR = loadOffsetMap(argument1);
        offsetMapS = loadOffsetMap(argument2);
        //in order to use the offset map efficiently for information retrieval
        //  regarding the geometries, we need to keep the files open at all times
        finR.open(geometryFileNameR, fstream::in | ios_base::binary);
        finS.open(geometryFileNameS, fstream::in | ios_base::binary);

        if(!finR || !finS){
                cout << "Error opening one of the two files." << endl;
                exit(-1);
        }

        //set universal coordinates
        cout << "Setting universal min/max..." << endl;
        setUniversalCoordinates(finR);
        setUniversalCoordinates(finS);
        cout << "Done: " << universalMinX << " " << universalMinY << "," << universalMaxX << " " << universalMaxY << endl;
        //return to begining
        finR.seekg(0, ios::beg);
        finS.seekg(0, ios::beg);

        ofstream fout(result_filename, fstream::out | ios_base::binary);
        fout.close();
}

//-----------------------------
//
//          ENABLERS
//
//-----------------------------

void enableIntermediateFilter(string &argument1, string &argument2){        
        rasterIntervalsR.argument = argument1;
        rasterIntervalsR.letterID = "A";
        rasterIntervalsS.argument = argument2;
        rasterIntervalsS.letterID = "B";

        loadApproximations(rasterIntervalsR, argument1);
        loadApproximations(rasterIntervalsS, argument2);
        
}

void initiateRasterIntervalsCreation(string &argument1, string &argument2){        
        createApproximations(argument1);
        createApproximations(argument2);
}

//-----------------------------
//
//          CONNECTORS
//
//-----------------------------

/* through this method, a post MBR filter candidate pair is forwarded
  further into the pipeline */
void forwardCandidatePair(uint idA, uint idB){
        int result;
        postMBRCandidates++;

        // Raster Intervals intermediate filter
        if(INTERMEDIATE_FILTER){
                timer = clock();
                result = joinPolygons(rasterIntervalsR.getPolygonByID(idA), rasterIntervalsS.getPolygonByID(idB), argument1, argument2);
                intermediateFilterTime += (clock() - timer) / (double) CLOCKS_PER_SEC;

                if(result == 1){
                        //accepted                        
                        TOTAL_RESULTS++;
                        accepted++;
                        //saveResultPair(idA, idB);
                        return;
                }else if(result == 0){
                        //else rejected     
                        rejected++;
                        return;
                }                           
        }

        // Geometric Refinement
        if(REFINEMENT){
                //pair needs refinement
                refinementCandidates++;  
                timer = clock();
                if(refinementWithIDs(idA, idB, offsetMapR, offsetMapS, finR, finS)){
                        TOTAL_RESULTS++;
                        acceptedAfterRefinement++;
                        //saveResultPair(idA, idB);
                }
                refinementTime += (clock() - timer) / (double) CLOCKS_PER_SEC;
        }
}

#endif
