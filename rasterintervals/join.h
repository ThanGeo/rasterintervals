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

#include "containers.h"
#include "join_geometry_refinement.h"
#include "binary_files.h"

using namespace std;

#ifndef JOIN_H
#define JOIN_H

clock_t local_timer;

/*
*-------------------------------------------------------
*
*     TIMER
*       
*
*-------------------------------------------------------
*/

void saveTimer(uint pairNumber, double &refinementTime){
	ofstream fout("data/refinement_candidates/refinement_times_" + to_string(HILBERT_n) + ".csv", fstream::out | ios_base::binary | fstream::app);
	fout << pairNumber << " " << refinementTime << endl;
	fout.close();
}

/*
*-------------------------------------------------------
*
*     RASTER INTERVALS JOIN
*       
*
*-------------------------------------------------------
*/

//sets all the required parameters for the next interval in the byte array (start, end, offset)
void nextinterval(uint8_t *ar, int *offset, uint *st, uint *end){
	int numcells;	
	(*st) = *(uint *)(ar+(*offset));
    (*offset)+=sizeof(uint);
    (*end) = *(uint *)(ar+(*offset));
    (*offset)+=sizeof(uint);

	numcells = (*end) - (*st) +1;
    (*offset) += (int)ceil((double)(numcells*3.0)/8.0); //skip bitarray
}


/* performs the interval alignment through bit shifting and the bitwise comparison between 2 aligned intervals */
int shiftandmatch(uint8_t *ar1, uint st1, uint end1, int offset1, uint8_t *ar2, uint st2, uint end2, int offset2){
	uint8_t onebits;
	uint8_t curbyte;
	uint8_t xor_mask_byte;
	uint8_t temp;	
	uint stov,endov;	
	int stdiff,bitdiff;
	int ar1pos,ar2pos;
	int numpos1,numpos2; //number of bit positions 
	int numcomparedbytes; 
	
	//smallest endpoint in both intervals
	endov = (end1<end2) ? end1 : end2; 

	//keeps track of the shifting required for the XOR mask
	uint8_t loopCounter = 0;
	
	if (st1<st2) // find start position in bitarray of interval1 
	{
		stdiff = (st2-st1)*3; // diff in start positions (number of bits)
		numpos1 = (end1-st1+1)*3; // how many bits in object 1 code
		ar1pos = offset1-(int)ceil((double)(numpos1)/8.0)+stdiff/8; //move fwd stdiff/8 bytes 
		bitdiff = stdiff % 8; //bitdiff becomes stdiff%8 (for shifting)
		uint8_t mask = (int)(pow(2,bitdiff)-0.5); // 2^bitdiff - 1

		ar2pos = offset2-(int)ceil((double)((end2-st2+1)*3)/8.0);		
		numcomparedbytes = (int)ceil((double)((endov-st2+1)*3)/8.0);		

		//check each byte required (numcomparedbytes)
		for(int i = ar1pos; i<ar1pos+numcomparedbytes; i++) {			
			//carryover information from the next byte
			onebits = ar1[i+1] & mask;
			//align
			curbyte = ar1[i] >> bitdiff;
			//fill in with the carried over bits after the shifting
			curbyte |= onebits << (8-bitdiff);

			//TODO remove this initialization, add 182 inside the if
			xor_mask_byte = 182;
			//modify the XOR mask accordingly as to be aligned with the shifted coding
			if(loopCounter > 0){
				xor_mask_byte = xor_mask_byte << loopCounter;
				if(loopCounter == 1){
					xor_mask_byte = xor_mask_byte | 1;
				}else{
					xor_mask_byte = xor_mask_byte | 3;
				}				
			}

			//XOR the shifted coding
			curbyte = curbyte ^ xor_mask_byte;
			temp = ar2[ar2pos+i-ar1pos];

			//if we are in the last byte needed to be checked,
			//  mask out any excess bits from both codings
			if(i == ar1pos+numcomparedbytes-1){
				curbyte = curbyte & (255 >> (8 - (endov-st2+1)*3)%8);
				temp = temp & (255 >> (8 - (endov-st2+1)*3)%8);
			}
			
			//move the loopcounter to the next value (cyclically 0,1,2, then it repeats itself)
			loopCounter = (loopCounter + 1) % 3;
			
			//bitwise AND, returns !=0 if they are a true hit
			if(curbyte & (temp)){
				return 1; // true hit
			}
		}
	}
	else{ // st1>=st2: symmetric to case st1<st2	
		stdiff = (st1-st2)*3; // diff in start positions (number of bits)
		numpos2 = (end2-st2+1)*3;
		ar2pos = offset2-(int)ceil((double)(numpos2)/8.0)+stdiff/8; //move fwd stdiff/8 bytes 
		bitdiff = stdiff % 8; //bitdiff becomes stdiff%8 (for shifting)
		uint8_t mask = (int)(pow(2,bitdiff)-0.5); // 2^bitdiff - 1
		
		ar1pos = offset1-(int)ceil((double)((end1-st1+1)*3)/8.0);		
		numcomparedbytes = (int)ceil((double)((endov-st1+1)*3)/8.0);

		//check each byte required (numcomparedbytes)
		for(int i = ar2pos; i<ar2pos+numcomparedbytes; i++) {
			//carryover information from the next byte
			onebits = ar2[i+1] & mask;
			//align
			curbyte = ar2[i] >> bitdiff;
			//fill in with the carried over bits after the shifting		
			curbyte |= onebits << (8-bitdiff);

			//the transformation mask for the XOR operation (1 byte) 10110110
			xor_mask_byte = 182;

			//MODIFY THE XOR MASK
			if(loopCounter > 0){
				xor_mask_byte = xor_mask_byte << loopCounter;
				if(loopCounter == 1){
					xor_mask_byte = xor_mask_byte | 1;
				}else{
					xor_mask_byte = xor_mask_byte | 3;
				}
			}

			//XOR
			curbyte = curbyte ^ xor_mask_byte;
			temp = ar1[ar1pos+i-ar2pos];

			//IF WE ARE IN THE LAST BYTE, BUT WE MAY NOT WANT ALL THE BITS
			if(i == ar2pos+numcomparedbytes-1){
				curbyte = curbyte & (255 >> (8 - (endov-st1+1)*3)%8);
				temp = temp & (255 >> (8 - (endov-st1+1)*3)%8);
			}
			
			loopCounter = (loopCounter + 1) % 3;

			if (curbyte & (temp)){
				return 1; // true hit
			}
		}
	}
	
	//no matches, indecisive case
	return 0;
}

/* merge-joins the two objects' intervals and recognizes overlaps to move on with comparisons */
int compareIntervalsBITS(uint8_t *ar1, uint8_t *ar2){
	bool candidateFlag = false;

	//get total number of intervals in both objects
	uint numintervals1 = *(uint *)ar1;
	uint numintervals2 = *(uint *)ar2;

	//byte arrays starting offsets
	int offset1=sizeof(uint);
	int offset2=sizeof(uint);
	
	uint st1,st2,end1,end2;
	uint stov,endov;	
	uint cur1=0;
	uint cur2=0;
	
	//get the first interval of each object
	uint numcells;
	nextinterval(ar1,&offset1,&st1,&end1);	
    cur1++;    
	nextinterval(ar2,&offset2,&st2,&end2);
    cur2++;

    //iterate the objects' intervals until we reach the end of one of their lists
	do{
		//symmetrical
		if (st1<=st2){ // current interval from object 1 starts before the current interval from 2
			if(end1>=st2){ // they overlap			
				//so they are candidates for ref
				candidateFlag = true;		
				
				//perform the bit alignment and comparison
				if(shiftandmatch(ar1,st1,end1,offset1,ar2,st2,end2,offset2)){
					return 1; // true hit (code 1)
				}				
				
				//move on to the next interval based on which ends first
				if(end1<=end2){
					//get next interval of obj 1
					nextinterval(ar1,&offset1,&st1,&end1);
					cur1++;
				}else{
					//get next interval of obj 2
					nextinterval(ar2,&offset2,&st2,&end2);
					cur2++;
				}
			}else{
				//they do not overlap, get next interval of obj 1
				// because it starts and ends before the interval from obj 2
				nextinterval(ar1,&offset1,&st1,&end1);
				cur1++;
			}
		}else{ // current interval from object 2 starts before the current interval from 1		
			if (end2>=st1){ // they overlap
				//so they candidates for ref
				candidateFlag = true;				
				
				//perform the bit alignment and comparison
				if(shiftandmatch(ar1,st1,end1,offset1,ar2,st2,end2,offset2)){
					return 1; // true hit (code 1)
				}

				//move on to the next interval based on which ends first
				if(end2<=end1){
					//get next interval of obj 2
					nextinterval(ar2,&offset2,&st2,&end2);
					cur2++;
				}else{
					//get next interval of obj 1
					nextinterval(ar1,&offset1,&st1,&end1);
					cur1++;
				}
			}else{
				//they do not overlap, get next interval of obj 2
				// because it starts and ends before the interval from obj 1
				nextinterval(ar2,&offset2,&st2,&end2);
				cur2++;
			}
		}
	}while(cur1<=numintervals1 && cur2<=numintervals2);
	
	//if they have at least 1 overlapping pair of intervals
	//  and no true hit has been recognized, then send to refinement (code 2)
	if(candidateFlag){
		return 2;
	}
	//no overlaps, no true hits, they definitely do not intersect (code 0)
	return 0;
}



/* loads the requested polygon geometry from disk, using a byte offset map for faster retrieval */
Polygon loadPolygonGeometry(uint &recID, unordered_map<uint,unsigned long> &offsetMap, ifstream &fin){
	Polygon pol(recID);
	int readID;
	int vertexCount, polygonCount;
	double x,y;
	//polygon's MBR, created on the go
	double polxMin = numeric_limits<int>::max();
	double polyMin = numeric_limits<int>::max();
	double polxMax = -numeric_limits<int>::max();
	double polyMax = -numeric_limits<int>::max();

	//search the map for the specific polygon offset
	unordered_map<uint,unsigned long>::const_iterator got = offsetMap.find(recID);
	if(got != offsetMap.end()){
		//key exists

		//set read offset
		fin.seekg(got->second-fin.tellg(), fin.cur);
		
		//read obj ID
		fin.read((char*) &readID, sizeof(int));

		//read vertex count
		fin.read((char*) &vertexCount, sizeof(int));
		pol.vertices.reserve(vertexCount);

		//read vertices
		for(int i=0; i<vertexCount; i++){
			fin.read((char*) &x, sizeof(double));
			fin.read((char*) &y, sizeof(double));

			pol.vertices.emplace_back(x,y);

			//save polygon's min/max (for mbr)
			polxMin = min(polxMin, x);
			polyMin = min(polyMin, y);
			polxMax = max(polxMax, x);
			polyMax = max(polyMax, y);
		}

		//add MBR to polygon
		pol.mbr.pMin.x = polxMin;
		pol.mbr.pMin.y = polyMin;
		pol.mbr.pMax.x = polxMax;
		pol.mbr.pMax.y = polyMax;
	}

	return pol;
}

/* begins refinement for two requested objects through their IDs
	it requires an open fstream with the data and the offset map of each object
 */
bool refinementWithIDs(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS){
	//load geometries
	local_timer = clock();
	Polygon geoR = loadPolygonGeometry(idA, offsetMapR, finR);
	Polygon geoS = loadPolygonGeometry(idB, offsetMapS, finS);
	//count refinement candidates for each set separately
	refinementCandidatesR += geoR.vertices.size();
	refinementCandidatesS += geoS.vertices.size();
	loadingGeometriesTime += (clock() - local_timer) / (double) CLOCKS_PER_SEC;

	//if any of the two has 0 vertices, return false (problematic data?)
	if(geoR.vertices.size() == 0 || geoS.vertices.size() == 0){
		return false;
	}

	//initiate plane sweep refinement
	return refinePlaneSweep(geoR, geoS);
}


//performs the join between two polygons
int joinPolygons(Polygon *polA, Polygon *polB, string &argument1, string &argument2){		
	return compareIntervalsBITS(polA->coding_data, polB->coding_data);
}


#endif
