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

#include <omp.h>

#include "containers.h"
#include "hilbert_functions.h"

#include "earcut.hpp"

using namespace std;

#ifndef RASTERIZATION_H
#define RASTERIZATION_H

/*
*-------------------------------------------------------
*
*     SORT
*       POINTS
*
*-------------------------------------------------------
*/
Point centroid;

Point find_centroid(vector<Point> &points) {
	double x = 0, y = 0;
	for (int i = 0; i < points.size(); i++) {
		x += points[i].x;
		y += points[i].y;
	}

	Point c(x / points.size(), y / points.size());
	return c;
}

int by_polar_angle(Point &p1, Point &p2) {
	double theta_a = p1.to_angle(centroid);
	double theta_b = p2.to_angle(centroid);
	/*
	cout << "theta a " << theta_a << endl;
	cout << "theta b " << theta_b << endl;
	cout << endl;
	*/
	return (theta_a < theta_b);
}

void sort_by_polar_angle(vector<Point> &points) {
	centroid = find_centroid(points);
	/*
	cout << "Centroid: " << centroid.x << " " << centroid.y << endl;

	cout << "BEFORE SORT: " << endl;
	for(int i=0; i<points.size(); i++){
		cout << points.at(i).x << " " << points.at(i).y << ",";
	}
	cout << endl;
	*/
	sort(points.begin(), points.end(), by_polar_angle);

	/*
	cout << "AFTER SORT: " << endl;
	for(int i=0; i<points.size(); i++){
		cout << points.at(i).x << " " << points.at(i).y << ",";
	}
	cout << endl;
	*/
}


/*
*-------------------------------------------------------
*
*     GEOMETRIC
*       CALCULATIONS
*
*-------------------------------------------------------
*/

double computePolygonArea(vector<Point> &points){
	//MAPBOX METHOD https://github.com/mapbox/earcut.hpp
	// The number type to use for tessellation
	using Coord = double;

	// The index type. Defaults to uint32_t, but you can also pass uint16_t if you know that your
	// data won't have more than 65536 vertices.
	using N = uint32_t;
	
	vector<vector<Point>> polygon = {{points}};
	return mapbox::earcutGetArea<N>(polygon);
}


/*
*-------------------------------------------------------
*
*     POLYGON
*       CLASSIFICATION
*
*-------------------------------------------------------
*/

/* R coding (default) */
int classifySubpolygon(double area){
	if(area == 1.0){
		//full
		return 3;
	}else if(area > 0.5){
		//strong
		return 5;
	}else if(area == 0){
		//empty
		return 0;
	}else{
		//weak
		return 4;
	}	
}

/* S coding (this is never used, since we mask during the join) */
int classifySubpolygonOTHER(double area){
	if(area == 1.0){
		//full
		return 5;
	}else if(area > 0.5){
		//strong
		return 3;
	}else if(area == 0){
		//empty
		return 0;
	}else{
		//weak
		return 2;
	}	
}

/*
*-------------------------------------------------------
*
*     POLYGON
*       RASTERIZATION
*
*-------------------------------------------------------
*/

double getSlope(Point &A, Point &B){
	return ((B.y - A.y) / (B.x - A.x));
}


bool isInsideHorizontalDual(double Yi, double Yi1, double &py){
	if(py < Yi || py > Yi1){
		return false;
	}
	return true;
}

bool isInsideVerticalDual(double Xi, double Xi1, double &px){
	if(px < Xi || px > Xi1){
		return false;
	}
	return true;
}

/* slices a polygon vertically using 2 parallel lines Xi and Xi1 */
TempPolygon intersectionDualX(Polygon &pol, double &Xi, double &Xi1){
	double slope;
	double yi, yi1;
	TempPolygon clippedPolygon(pol.recID);

	//for each edge of the polygon
	// ITERATES THE VERTICES WHICH MUST BE IN A COUNTER-CLOCKWISE DIRECTION!!!
	// so that right of the line means inside, left means outside	
	for(auto itV = pol.vertices.begin(); itV != pol.vertices.end()-1; itV++){
		Point pointA = *itV;
		Point pointB = *(itV+1);
		//check cases
		if(pointA.x == pointB.x){
			//edge is vertical
			//set the lowest point as the intersection point
			yi = min(pointA.y, pointB.y);
			yi1 = min(pointA.y, pointB.y);
		}else{
			//solve for y
			slope = getSlope(pointA,pointB);

			//find intersection points for both Xi and Xi1
			yi = slope * (Xi - pointB.x) + pointB.y;
			yi1 = slope * (Xi1 - pointB.x) + pointB.y;
		}

		//intersection points
		Point pi(Xi, yi);
		Point pi1(Xi1, yi1);
		//cout << "X CLIPPING " << endl;
		//POINT A
		if(isInsideVerticalDual(Xi, Xi1, pointA.x)){
			//pointA inside, clipping result
			clippedPolygon.addPoint(pointA);
			//cout << "point A (" << pointA.x << " " << pointA.y << ") added" << endl;
		}

		//intersection points in PROPER ORDER
		if(pointA.x < pointB.x){
			//first pi, then pi1
			if(isInsideHorizontalDual(min(pointA.y, pointB.y), max(pointA.y, pointB.y), pi.y) && isInsideVerticalDual(min(pointA.x, pointB.x), max(pointA.x, pointB.x), pi.x)){
				//intersection point pi is a clipping result
				clippedPolygon.addPoint(pi);
				//cout << "point pi (" << pi.x << " " << pi.y << ") added" << endl;
			}

			if(isInsideHorizontalDual(min(pointA.y, pointB.y), max(pointA.y, pointB.y), pi1.y) && isInsideVerticalDual(min(pointA.x, pointB.x), max(pointA.x, pointB.x), pi1.x)){
				//intersection point pi1 is a clipping result
				clippedPolygon.addPoint(pi1);
				//cout << "point pi1 (" << pi1.x << " " << pi1.y << ") added" << endl;
			}
		}else{
			//first pi1 then pi
			if(isInsideHorizontalDual(min(pointA.y, pointB.y), max(pointA.y, pointB.y), pi1.y) && isInsideVerticalDual(min(pointA.x, pointB.x), max(pointA.x, pointB.x), pi1.x)){
				//intersection point pi1 is a clipping result
				clippedPolygon.addPoint(pi1);
				//cout << "point pi1 (" << pi1.x << " " << pi1.y << ") added" << endl;
			}
			if(isInsideHorizontalDual(min(pointA.y, pointB.y), max(pointA.y, pointB.y), pi.y) && isInsideVerticalDual(min(pointA.x, pointB.x), max(pointA.x, pointB.x), pi.x)){
				//intersection point pi is a clipping result
				clippedPolygon.addPoint(pi);
				//cout << "point pi (" << pi.x << " " << pi.y << ") added" << endl;
			}			
		}

		//POINT B
		if(isInsideVerticalDual(Xi, Xi1, pointB.x)){
			//pointB inside, clipping result
			clippedPolygon.addPoint(pointB);
			//cout << "point B (" << pointB.x << " " << pointB.y << ") added" << endl;
		}
	}
	//add the first point to the end
	// so that the polygon "closes"
	if(clippedPolygon.vertices.size() != 0){
		clippedPolygon.vertices.push_back(*clippedPolygon.vertices.begin());
	}

	return clippedPolygon;
}

/* slices the sub-polygons created by vertical slicing horizontally, to create the raster cells */
vector<Point> intersectionDualY(TempPolygon &pol, double Yi, double Yi1){
	double slope;
	double xi, xi1;
	TempPolygon clippedPolygon(pol.recID);
		
	for(auto itV = pol.vertices.begin(); itV != pol.vertices.end()-1; itV++){
		Point pointA = *itV;
		Point pointB = *(itV+1);

		bool parallels = false;
		//calculate the slope (if any)
		if(pointA.x == pointB.x){
			//edge is vertical
			//the only possible point of intersection is y (Yi or Yi1)
			// do nothing
			xi = pointA.x;
			xi1 = pointA.x;
		}else{
			slope = slope = getSlope(pointA,pointB);
			//check if AB is horizontal, if it is it only intersects if pointA.y == Yi or pointA.y == Yi1
			//this flag will skip unnecessary checks in this event
			if(pointA.y == pointB.y){
				if(pointA.y == Yi){				
					xi = pointA.x;
					xi1 = -numeric_limits<double>::max();
				}else if(pointA.y == Yi1){
					xi = -numeric_limits<double>::max();
					xi1 = pointA.x;
				}else{
					//THEY DO NOT INETERSECT, THEY ARE PARALLELS
					parallels = true;
				}
			}else{
				//solve for x				
				xi = ((Yi - pointB.y) / slope) + pointB.x;
				xi1 = ((Yi1 - pointB.y) / slope) + pointB.x;
				//cout << "   solved for x" << endl;
			}
		}

		if(!parallels){
			//intersection points
			Point pi(xi, Yi);
			Point pi1(xi1, Yi1);

			//POINT A
			if(isInsideHorizontalDual(Yi, Yi1, pointA.y)){
				//pointA inside, clipping result
				clippedPolygon.addPoint(pointA);
				//cout << "point A (" << pointA.x << " " << pointA.y << ") added" << endl;
			}

			//intersection points in PROPER ORDER
			if(pointA.y < pointB.y){
				//first pi, then pi1
				if(isInsideVerticalDual(min(pointA.x, pointB.x), max(pointA.x, pointB.x), pi.x) && isInsideHorizontalDual(min(pointA.y, pointB.y), max(pointA.y, pointB.y), pi.y)){
					//intersection point pi is a clipping result
					clippedPolygon.addPoint(pi);
					//cout << "point pi (" << pi.x << " " << pi.y << ") added" << endl;
				}

				if(isInsideVerticalDual(min(pointA.x, pointB.x), max(pointA.x, pointB.x), pi1.x) && isInsideHorizontalDual(min(pointA.y, pointB.y), max(pointA.y, pointB.y), pi1.y)){
					//intersection point pi1 is a clipping result
					clippedPolygon.addPoint(pi1);
					//cout << "point pi1 (" << pi1.x << " " << pi1.y << ") added" << endl;
				}
			}else{
				//first pi1 then pi
				if(isInsideVerticalDual(min(pointA.x, pointB.x), max(pointA.x, pointB.x), pi1.x) && isInsideHorizontalDual(min(pointA.y, pointB.y), max(pointA.y, pointB.y), pi1.y)){
					//intersection point pi1 is a clipping result
					clippedPolygon.addPoint(pi1);
					//cout << "point pi1 (" << pi1.x << " " << pi1.y << ") added" << endl;
				}
				if(isInsideVerticalDual(min(pointA.x, pointB.x), max(pointA.x, pointB.x), pi.x) && isInsideHorizontalDual(min(pointA.y, pointB.y), max(pointA.y, pointB.y), pi.y)){
					//intersection point pi is a clipping result
					clippedPolygon.addPoint(pi);
					//cout << "point pi (" << pi.x << " " << pi.y << ") added" << endl;
				}			
			}

			//POINT B
			if(isInsideHorizontalDual(Yi, Yi1, pointB.y)){
				//pointB inside, clipping result
				clippedPolygon.addPoint(pointB);
				//cout << "point B (" << pointB.x << " " << pointB.y << ") added" << endl;
			}
		}	
	}

	//sort points
	sort_by_polar_angle(clippedPolygon.vertices);

	//"close" the polygon (first and last points in order must be the same point)
	if(clippedPolygon.vertices.size() != 0){
		clippedPolygon.vertices.push_back(*clippedPolygon.vertices.begin());
	}

	return clippedPolygon.vertices;
}

/*
*-------------------------------------------------------
*
*
*     RASTERIZATION
*       
*
*-------------------------------------------------------
*/

vector<Cell> rasterizePolygon(Polygon &mappedPol, string argument){
	vector<Cell> rasterizationCells;
	vector<Point> clippedPoints;
	vector<TempPolygon> subpolygonsAfterX;

	double Xi, Xi1, Yi, Yi1;
	double kx, ky;

	//initialize the x axis sweep lines
	Xi = floor(mappedPol.mbr.pMin.x);
	Xi1 = Xi + 1;

	//define the end x
	kx = ceil(mappedPol.mbr.pMax.x) + 1;

	TempPolygon tempPol;
	subpolygonsAfterX.reserve(kx - Xi);

	//sweep the x axis getting pairs of vertical lines Xi & Xi+1
	while(Xi1<kx){
		//cout << Xi << " and " << Xi1 << "/" << kx << endl;
		//returns the sub polygon when pol is clipped by Xi and Xi1
		tempPol = intersectionDualX(mappedPol, Xi, Xi1);
		//add the sub polygon to the polygon's list
		// we need to save it to further clip it in the y axis later
	
		if(tempPol.vertices.size() > 0){
			tempPol.cellX = Xi;			
			subpolygonsAfterX.push_back(tempPol);
		}	

		//move both vertical lines equally
		Xi += 1;
		Xi1 = Xi + 1;
	}

	//define the end y
	ky = ceil(mappedPol.mbr.pMax.y) + 1;
	
	int type;

	//iterate the subpolygons created by the x axis clipping instead of the original polygons
	auto it = subpolygonsAfterX.begin();
	while(it != subpolygonsAfterX.end()){

		//FOR NORMALIZED
		Yi = floor(mappedPol.mbr.pMin.y);
		Yi1 = Yi + 1;

		//sweep the y axis getting pairs of horizontal lines Yi & Yi+1
		while(Yi1<ky){
			//returns the subpolygon furtherly clipped in the y axis by Yi and Yi+1
			clippedPoints = intersectionDualY(*it, Yi, Yi1);	

			//this helps ignore a large portion of the empty cells for a polygon
			if(clippedPoints.size() > 2){
				//calculate its area and classify it
				
				type = classifySubpolygon(computePolygonArea(clippedPoints));

				if(type != 0){
					rasterizationCells.emplace_back(it->cellX, Yi, it->cellX + 1, Yi1, type);	
				}
				
			}

			//move the horizontal lines equally to the next position
			Yi += 1;
			Yi1 = Yi + 1;
		}
		it++;
	}
	return rasterizationCells;
}
#endif
