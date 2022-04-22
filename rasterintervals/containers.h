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
#include <unordered_map>

using namespace std;

#ifndef STRUCTURE_H
#define STRUCTURE_H

double e = 1e-08;
//order N of the hilbert curve (must be a power of 2)
int HILBERT_n = 65536;

#define PI 3.14159265

//universal coordinates of the 2 datasets
double universalMinX = std::numeric_limits<Coord>::max();
double universalMinY = std::numeric_limits<Coord>::max();
double universalMaxX = -std::numeric_limits<Coord>::max();
double universalMaxY = -std::numeric_limits<Coord>::max();

int DECIMAL_POINTS_PRECISION = 6;

double loadingGeometriesTime = 0;

uint totalIntervals = 0;
uint totalCells = 0;
unsigned long refinementCandidatesR = 0;
unsigned long refinementCandidatesS = 0;

#define BYTE_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define BYTE_TO_BINARY(byte)  \
  (byte & 0x80 ? '1' : '0'), \
  (byte & 0x40 ? '1' : '0'), \
  (byte & 0x20 ? '1' : '0'), \
  (byte & 0x10 ? '1' : '0'), \
  (byte & 0x08 ? '1' : '0'), \
  (byte & 0x04 ? '1' : '0'), \
  (byte & 0x02 ? '1' : '0'), \
  (byte & 0x01 ? '1' : '0') 

/*
*-------------------------------------------------------
*
*     CLASSES
*       
*
*-------------------------------------------------------
*/

class Point{
public:
	double x, y;
	Point(double x, double y){
		this->x = x;
		this->y = y;
	}
	Point(){}

	bool operator< (const Point &other) const{
        return x < other.x;
    }

    bool operator== (const Point &other) const{
        return (this->x == other.x && this->y == other.y);
	}

	double to_angle(Point &o){   
		return atan2(y - o.y, x - o.x);
	}
};

// A struct that will contain a line segment of a polygon. For plane sweep
struct LineSegment
{
    Point startPoint, endPoint;
    bool firstPoint; // If firstPoint is true, then p1 is the new Point, else its p2.
};

typedef vector<LineSegment>::iterator LineSegmentIterator;

class MBR{
public:
	Point pMin;
	Point pMax;
	MBR(double xS, double yS, double xE, double yE){
		pMin = Point(xS, yS);
		pMax = Point(xE, yE);
	}
	void set(double xS, double yS, double xE, double yE){
		pMin = Point(xS, yS);
		pMax = Point(xE, yE);
	}
	MBR(){};
};

/* temp polygon class used only in rasterization */
class TempPolygon{
public:
	uint recID;
	vector<Point> vertices;
	double cellX, cellY;
	TempPolygon(){}
	
	TempPolygon(uint &recID){
		this->recID = recID;
	}

	void addPoint(Point &p){
		if(find(vertices.begin(), vertices.end(), p) == vertices.end()){
			vertices.push_back(p);
		}

		//vertices.push_back(p);
	}

};


/* raster cell */
class Cell{
public:
	int classificationID;
	Point bottomLeft;
	Point topRight;

	Cell(){};

	Cell(double bottomLeftX, double bottomLeftY, double topRightX, double topRightY){
		bottomLeft.x = bottomLeftX;
		bottomLeft.y = bottomLeftY;
		topRight.x = topRightX;
		topRight.y = topRightY;
	}

	Cell(double bottomLeftX, double bottomLeftY, double topRightX, double topRightY, int classificationID){
		bottomLeft.x = bottomLeftX;
		bottomLeft.y = bottomLeftY;
		topRight.x = topRightX;
		topRight.y = topRightY;
		this->classificationID = classificationID;
	}


};


/* raster interval */
class NewInterval{
public:
	uint start;
	uint end;
	uint8_t* color;

	NewInterval(uint &s, uint &e){
		start = s;
		end = e;
	}

	NewInterval(uint &s, uint &e, uint8_t *c){
		start = s;
		end = e;
		
		color = (uint8_t*) malloc(((e-s+1) + 1) * sizeof(uint8_t));
		memcpy(color, c, (e-s+1) + 1);
	}
	NewInterval(){}
};

/* polygon object */
class Polygon{
public:	
	//original mbr
	MBR mbr;
	//---rasterization/intervals---
	//only for computation, then deleted
	vector<Cell> rasterizationCells;
	vector<uint> hilbertCellIDs;
	unordered_map<uint, Cell> hilbertCells;

	//only these are needed for RI join, the rest are cleared
	uint recID;
	vector<Point> vertices;
	vector<NewInterval> hilbertIntervalsNEW;

	//for planesweep refinement
	vector<LineSegment> exactGeometry_LS;

	//for bit coding
	uint8_t* coding_data;
	uint numBytes;

	Polygon(uint &recID){
		this->recID = recID;
	}

	Polygon(){};

	void addPoint(Point &p){
		vertices.push_back(p);
	}

	void addHilbertCell(Cell &cell, uint &cellID){		
		hilbertCellIDs.push_back(cellID);
		hilbertCells.insert(make_pair(cellID, cell));
	}
};

/* dataset container */
class Dataset{
public:
	string filename;
	unordered_map<uint,Polygon> polygons;
	uint maxIntervalBytes;
	
	Point pMin, pMax;
	string letterID;

	//for example T1, T2 etc...
	string argument;

	Dataset(){};

	Dataset(string &filename){
		this->filename = filename;
	}

	Dataset(string letterID){
		this->letterID = letterID;
	}

	Dataset(string argument, string letterID){
		this->argument = argument;
		this->letterID = letterID;
	}

	Polygon* getPolygonByID(uint &recID){
		auto it = polygons.find(recID);
		if(it != polygons.end()){
			return &(it->second);
		}
		return NULL;
	}

};

void printContainer(uint8_t *container, uint &totalBytes){
	cout << "CONTAINER:" << endl;
	for(int i = 0; i<totalBytes; i++){
		printf(""BYTE_TO_BINARY_PATTERN" ",BYTE_TO_BINARY(container[i]));
		if((i+1) % 4 == 0){

			cout << endl;
		}
	}
	cout << endl;
}


//DECLARE GLOBAL GEOMETRY DATASETS
Dataset geometryDatasetA("A");
Dataset geometryDatasetB("B");


#endif