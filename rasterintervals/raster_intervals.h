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
#include <iomanip>
#include <unistd.h>

#include "containers.h"
#include "join.h"
#include "rasterization.h"
#include "dataset_data.h"
#include "hilbert_functions.h"

using namespace std;

/*
*-------------------------------------------------------
*
*     GLOBAL VARIABLES
*       
*
*-------------------------------------------------------
*/

//temp containers. Used for the bit coding generation process
// (initialized statically to avoid dynamic memory allocation)
uint8_t* ar = new uint8_t[100000000];
uint8_t* encoding = new uint8_t[100000000];

/*
*-------------------------------------------------------
*
*     COMPUTE INTERVALS
*       
*
*-------------------------------------------------------
*/

/* normalizes a single point to the hilbert space */
void normalizeStraightToHilbert(Point &p){
	double x = (HILBERT_n-1) * (p.x - universalMinX) / (universalMaxX - universalMinX);
	double y = (HILBERT_n-1) * (p.y - universalMinY) / (universalMaxY - universalMinY);

	//cout << fixed << setprecision(DECIMAL_POINTS_PRECISION) << p.x << " " << p.y << " became " << x << " " << y << endl;

	p.x = double(x);
	p.y = double(y);

	if(p.x == HILBERT_n){
		p.x = HILBERT_n-1;
	}
	if(p.y == HILBERT_n){
		p.y = HILBERT_n-1;
	}
}

/* map an original polygon to the hilbert grid */
void mapPolygon(Polygon &originalPol){	
	//normalize each vertex
	for(auto it = originalPol.vertices.begin(); it != originalPol.vertices.end(); it++){
		normalizeStraightToHilbert(*it);
	}

	//normalize MBR
	normalizeStraightToHilbert(originalPol.mbr.pMin);
	normalizeStraightToHilbert(originalPol.mbr.pMax);
}

/* build the raster intervals and their codings */
void hilbertIntervalsForRasterizedPolygon(Polygon &pol){
	uint d;
	int x,y;
	int MAX_INTERVAL_LENGTH = 10;
	uint start;

	//compute the hilbert space cell IDs, based on the raster cells' coordinates
	for(auto it = pol.rasterizationCells.begin(); it != pol.rasterizationCells.end(); it++){		
		x = it->bottomLeft.x;
		y = it->bottomLeft.y;

		d = xy2d(HILBERT_n,x,y);
		pol.addHilbertCell(*it, d);
	}

	//bit coding 
	if(pol.hilbertCellIDs.size() > 0){
		//sort the IDs 
		sort(pol.hilbertCellIDs.begin(), pol.hilbertCellIDs.end());

		uint numcells, numbytes, curbitpos;
		uint numIntervals = 0;
		uint start, end;
		uint offset = sizeof(int);

		//in bytes
		uint totalSize = 0;

		//INITIALIZE (first interval)
		start = pol.hilbertCellIDs[0];
		memcpy(ar+offset, &start, sizeof(uint));
		offset+=sizeof(uint);
		numIntervals++;
		totalSize += sizeof(uint);

		//the cell codes
		vector<uint8_t> cell_codes = {pol.hilbertCells.find(start)->second.classificationID};

		//iterate the sorted cells
		for(auto it = pol.hilbertCellIDs.begin()+1; it != pol.hilbertCellIDs.end(); it++){
			//if the interval ended, save it
			// intervals must be of continuing cell IDs. If the difference between 2 IDs is greater than 1, then
			// a new interval must be created
			if(*it > *(it-1) + 1){				
				//end of an interval
				memcpy(ar+offset, &(*(it-1)), sizeof(uint));
				offset+=sizeof(uint);
				totalSize += sizeof(uint);

				//metrics
				totalCells += *(it-1) - start + 1;
			
				//coding
				numcells = *(it-1) - start + 1; // number of 3-bit codes
				numbytes = (int)ceil((double)numcells*3.0/8.0); // number of bytes necessary to store
				curbitpos = 0;
				for(auto code: cell_codes){
		        	for(int j=0;j<3;j++) { //check just first 3 positions
		        		if (code & 1 << j)
		        			encoding[curbitpos/8] |= 1 << (curbitpos%8); //set bit at curbitpos
		        		else
		        			encoding[curbitpos/8] &= ~(1 << (curbitpos%8)); //clear bit at curbitpos
		        		curbitpos++;
		        	}
		        }
		        for (int i=0; i<numbytes; i++)
		        {
		        	memcpy(ar+offset, &encoding[i], sizeof(uint8_t));
		        	offset+=sizeof(uint8_t);
		        }	
		        totalSize += numbytes; // * sizeof(uint8_t);

				//start new interval
				start = *it;
				memcpy(ar + offset, &start, sizeof(uint));
				offset+=sizeof(uint);
				totalSize += sizeof(uint);
				numIntervals++;
				
				//reset codes (add the new start)
				cell_codes.clear();
				cell_codes.push_back(pol.hilbertCells.find(*it)->second.classificationID);
			}else{
				//add the code of the current cell (rests inside the current interval)
				cell_codes.push_back(pol.hilbertCells.find(*it)->second.classificationID);
			}
		}				
		//close the last interval
		memcpy(ar+offset, &(*(pol.hilbertCellIDs.end()-1)), sizeof(uint));
		offset+=sizeof(uint);	
		totalSize += sizeof(uint);

		//metrics
		totalCells += *(pol.hilbertCellIDs.end()-1) - start + 1;
		numcells = *(pol.hilbertCellIDs.end()-1) - start + 1; // number of 3-bit codes
		numbytes = (int)ceil((double)numcells*3.0/8.0); // number of bytes necessary to store
		curbitpos = 0;
		for(auto code: cell_codes){
        	for(int j=0;j<3;j++) { //check just first 3 positions
        		if (code & 1 << j)
        			encoding[curbitpos/8] |= 1 << (curbitpos%8); //set bit at curbitpos
        		else
        			encoding[curbitpos/8] &= ~(1 << (curbitpos%8)); //clear bit at curbitpos
        		curbitpos++;
        	}
        }
        for (int i=0; i<numbytes; i++)
        {
        	memcpy(ar+offset, &encoding[i], sizeof(uint8_t));
        	offset+=sizeof(uint8_t);
        }	
        totalSize += numbytes; // * sizeof(uint8_t);

		//copy the total interval count into the container	
		memcpy(ar, &numIntervals, sizeof(int));	
		totalSize += sizeof(int);

		//metrics
		totalIntervals += numIntervals;

		//copy the coding (ar) into our polygon container
		pol.coding_data = new uint8_t[totalSize];
		memcpy(&pol.coding_data[0], &ar[0], totalSize);
		pol.numBytes = totalSize;
	}
}

/* creates the Raster Intervals for a selected dataset */
void computeIntervals(string &argument){
	uint lineCounter = 0;
	string line;
	clock_t timer;	
	uint recID;

	//metrics
	totalIntervals = 0;
	totalCells = 0;
	double sumOfAverageCellsPerInterval = 0;
	double averageCellsPerInterval = 0;
	double averageIntervalsPerPolygon = 0;

	//geometry file
	string filename = getBinaryGeometryFilename(argument);
	ifstream fin(filename, fstream::in | ios_base::binary);

	if(!fin){
		cout << "error opening " << filename << "." << endl;
		exit(-1);
	}

	//output files where the interval data will be stored
	string intervalFilename = getIntervalBinaryFilename(argument);
	ofstream fout(intervalFilename, ios_base::out | ios_base::binary);
	
	int totalPolygonCount, vertexCount;
	double x,y;
	double polxMin,polyMin,polxMax,polyMax;

	//read total polygon count from binary geometry file
	fin.read((char*) &totalPolygonCount, sizeof(int));

	//write total polygon count in the interval file
	fout.write((char*)(&totalPolygonCount), sizeof(int));

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
		//add MBR to polygon
		pol.mbr.pMin.x = polxMin;
		pol.mbr.pMin.y = polyMin;
		pol.mbr.pMax.x = polxMax;
		pol.mbr.pMax.y = polyMax;


		//---MAP POLYGON TO HILBERT SPACE---
		mapPolygon(pol);

		//---RASTERIZE POLYGON---
		pol.rasterizationCells = rasterizePolygon(pol, argument);

		//---INTERVALIZE POLYGON---
		hilbertIntervalsForRasterizedPolygon(pol);

		sumOfAverageCellsPerInterval += (double) totalCells / (double) totalIntervals;
		//---SAVE ON DISK---
		saveBinaryIntervalsBITS(pol, fout);		

		lineCounter++;
	}
	cout << "  Average intervals per polygon: " << (double) totalIntervals / (double) totalPolygonCount << " intervals/polygon." << endl;
	cout << "  Average cells per interval: " << (double) sumOfAverageCellsPerInterval / (double) totalPolygonCount << " cells/interval." << endl;
	
	fin.close();
	fout.close();
}


/*
*-------------------------------------------------------
*
*     MAIN APPROXIMATION
*       
*
*-------------------------------------------------------
*/

void createApproximations(string argument){
	cout << "Creating Raster Intervals approximations..." << endl;
	clock_t timer;
	//set the universal min max (needed for the global grid)
	//setUniversalMinMax(argument, universalMinX, universalMinY, universalMaxX, universalMaxY);
	//cout << fixed << setprecision(DECIMAL_POINTS_PRECISION) << "UNIVERSAL: min: " << universalMinX << " " << universalMinY << ", max: " << universalMaxX << " " << universalMaxY << endl;

	//compute the hilbert intervals
	timer = clock();		
	computeIntervals(argument);
	cout << fixed << setprecision(6) << "Computed intervals for dataset " << argument << " in " << (clock()-timer) / (double)(CLOCKS_PER_SEC) << " seconds" << endl;
	
}

void loadApproximations(Dataset &dataset, string argument){
	clock_t timer;
	timer = clock();
	//  load data	
	cout << "Loading " << getIntervalBinaryFilename(argument) << "..." << endl;
	loadBinaryRasterIntervalsBITS(dataset, getIntervalBinaryFilename(argument));
	cout << fixed << setprecision(6) << "Loaded in " << (clock()-timer) / (double)(CLOCKS_PER_SEC) << " seconds" << endl;		
}