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
#include <bitset>
#include <unordered_map>

using namespace std;

#ifndef BINARY_H
#define BINARY_H

/*
*-------------------------------------------------------
*
*     GEOMETRY BINARY DATA
*       
*
*-------------------------------------------------------
*/

/* 

binary offset map:
	
	contains byte offsets for each polygon regarding to the binary geometry file
	so that we can jump right to its exact position using seekg(offset) and load it faster
	
	map contains pairs of 
		uint: polygon ID
		usnigned long: byte offset

*/
unordered_map<uint,unsigned long> loadOffsetMap(string argument){
	unsigned long offset;
	uint lineCounter = 0;

	uint recID;

	ifstream fin(getOffsetMap(argument), fstream::out | ios_base::binary);

	unordered_map<uint,unsigned long> offset_map;

	int totalLines;

	//read total lines
	fin.read((char*) &totalLines, sizeof(int));

	while(lineCounter < totalLines){
		//read rec id
		fin.read((char*) &recID, sizeof(int));
		//read byte offset
		fin.read((char*) &offset, sizeof(unsigned long));

		offset_map.insert(make_pair(recID, offset));		
		lineCounter++;
	}
	

	fin.close();

	return offset_map;
}

/*
*-------------------------------------------------------
*
*     BIT CODING
*       
*
*-------------------------------------------------------
*/

void saveBinaryIntervalsBITS(Polygon &pol, ofstream &fout){
	//write ID
	fout.write((char*)(&pol.recID), sizeof(uint));
	//write number of bytes
	fout.write((char*)(&pol.numBytes), sizeof(uint));
	//write coding data
	fout.write((char*)(&pol.coding_data[0]), pol.numBytes);
}

void loadBinaryRasterIntervalsBITS(Dataset &set, string intervalFilename){
	uint totalPolygons,recID, numBytes;
	ifstream fin(intervalFilename, fstream::in | ios_base::binary);

	//read total polygons
	fin.read((char*) &totalPolygons, sizeof(uint));

	uint polCounter = 0;
	while(polCounter < totalPolygons){

		//read pol ID
		fin.read((char*) &recID, sizeof(uint));
		Polygon pol(recID);

		//read number of bytes
		fin.read((char*) &pol.numBytes, sizeof(uint));

		//reserve space for coding data
		pol.coding_data = new uint8_t[pol.numBytes];

		//read coding data
		fin.read((char*) &pol.coding_data[0], pol.numBytes);

		//add to set
		set.polygons.insert(make_pair(pol.recID, pol));
		polCounter++;
	}
	fin.close();
}



#endif
