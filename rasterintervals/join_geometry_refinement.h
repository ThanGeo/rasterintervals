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

#include "containers.h"
#include "dataset_data.h"

using namespace std;

/*
*-------------------------------------------------------
*
*     Point in Polygon and Edge Intersection
*       
*
*-------------------------------------------------------
*/

bool edgesIntersect(Point &aStart, Point &aEnd, Point &bStart, Point &bEnd){
	//make sure that they are both right edges and not points (depends on data)
	if(aStart.x == aEnd.x && aStart.y == aEnd.y){
		return false;
	}
	if(bStart.x == bEnd.x && bStart.y == bEnd.y){
		return false;
	}

	//cout << aStart.x << " " << aStart.y << "," << aEnd.x << " " << aEnd.y << " and " <<  bStart.x << " " << bStart.y << "," << bEnd.x << " " << bEnd.y << endl;

	//if they share a common point, they intersect
	if(aStart.x == bStart.x && aStart.y == bStart.y){
		return true;
	}
	if(aStart.x == bEnd.x && aStart.y == bEnd.y){
		return true;
	}
	if(aEnd.x == bStart.x && aEnd.y == bStart.y){
		return true;
	}
	if(aEnd.x == bEnd.x && aEnd.y == bEnd.y){
		return true;
	}

	//check intersection
	//from https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection 
	// in the LINE SEGMENT SECTION
	double denominator = (aStart.x - aEnd.x) * (bStart.y - bEnd.y) - (aStart.y - aEnd.y) * (bStart.x - bEnd.x);
	double t = ( (aStart.x - bStart.x) * (bStart.y - bEnd.y) - (aStart.y - bStart.y) * (bStart.x - bEnd.x) ) / denominator;
	double u = ( (aStart.x - bStart.x) * (aStart.y - aEnd.y) - (aStart.y - bStart.y) * (aStart.x - aEnd.x) ) / denominator;

	if(0 <= t && t <= 1 && 0 <= u && u <= 1){		
		return true;
	}
	return false;

}

/* checks if a given point rests inside (or on the outline) of a given polygon,
 using the ray tracing method */
bool pointInPolygon(Point &p, Polygon &pol){
	double e = 1e-10;
	Point lowPoint, highPoint;
	double red, blue;
	int numberOfIntersections = 0;	
	
	//iterate the edges
	auto it = pol.vertices.begin();
	while(it != pol.vertices.end()){
		//get the low/high points of the edge on Y		
		if(it == pol.vertices.end()-1){
			if(it->y > pol.vertices.begin()->y){
				lowPoint = *pol.vertices.begin();
				highPoint = *it;
			}else{
				lowPoint = *it;
				highPoint = *pol.vertices.begin();
			}
		}else{
			if(it->y > (it+1)->y){
				lowPoint = *(it+1);
				highPoint = *it;
			}else{
				lowPoint = *it;
				highPoint = *(it+1);
			}
		}

		//FIX THE POINT IF IT INTERSECTS A VERTEX on Y
		if(p.y == lowPoint.y || p.y == highPoint.y){
			p.y += e;
		}

		//count intersections
		if(!(p.y < lowPoint.y || p.y > highPoint.y) && (p.x < max(lowPoint.x, highPoint.x))){
			if(p.x < min(lowPoint.x, highPoint.x)){
				numberOfIntersections++;
			}else{
				if(lowPoint.x != highPoint.x){
					red = (highPoint.y - lowPoint.y) / (highPoint.x - lowPoint.x);
				}else{
					red = numeric_limits<double>::max();
				}
				if(lowPoint.x != p.x){
					blue = (p.y - lowPoint.y) / (p.x - lowPoint.x);
				}else{
					blue = numeric_limits<double>::max();
				}

				if(blue >= red){
					numberOfIntersections++;
				}
			}
		}
		it++;
	}

	if(numberOfIntersections % 2 == 1){
		return true;
	}
	return false;
}

/*
*-------------------------------------------------------
*
*     PLANE SWEEP REFINEMENT  
*
*-------------------------------------------------------
*/

bool compareByYmin(LineSegment &ls1, LineSegment &ls2)
{
    return (ls1.startPoint.y < ls2.startPoint.y);
}

bool internalPlaneSweep(LineSegmentIterator rec, LineSegmentIterator firstLS, LineSegmentIterator lastLS){
    bool resultFlag = false;
    auto pivot = firstLS;

    // While we still have Objects and pol1.yMax >= pol2.yMin 
    while((pivot < lastLS) && (rec->endPoint.y >= pivot->startPoint.y)){
        double recXmin = min(rec->startPoint.x, rec->endPoint.x);
        double recXmax = max(rec->startPoint.x, rec->endPoint.x);

        double pivotXmin = min(pivot->startPoint.x, pivot->endPoint.x);
        double pivotXmax = max(pivot->startPoint.x, pivot->endPoint.x);

        // Check if the Objects intersect (for the x axis)
        if((recXmin > pivotXmax) || (recXmax < pivotXmin)){        
            // If the Objects don't intersect, move to the next Object.
            pivot++;
            continue;
        }

        LineSegment ls1 = *rec;
        LineSegment ls2 = *pivot;

        resultFlag = edgesIntersect(ls1.startPoint, ls1.endPoint, ls2.startPoint, ls2.endPoint);
        pivot++;

        if(resultFlag == true){
            return resultFlag;
        }
    }
    return false;
}

bool planeSweep(vector<LineSegment> &polygon1, vector<LineSegment> &polygon2){
    bool resultFlag = false;
    auto pol1 = polygon1.begin();
    auto pol2 = polygon2.begin();
    auto lastPol1 = polygon1.end();
    auto lastPol2 = polygon2.end(); 

    while((pol1 < lastPol1) && (pol2 < lastPol2)){
        if (pol1->startPoint.y < pol2->startPoint.y){ // Find the Object with the lowest yMin.        
            // Run internal loop.
            resultFlag = internalPlaneSweep(pol1, pol2, lastPol2);
            pol1++;
        }else{
            // Run internal loop.
            resultFlag = internalPlaneSweep(pol2, pol1, lastPol1);
            pol2++;
        }

        if(resultFlag == true){
            return resultFlag;
        }
    }   
    return false;
}

/* finds the common MBR area between 2 polygons, used to reduce the amount of edges that need to 
 be refined in the edge intersection step */
MBR getCMBR(Polygon &polA, Polygon &polB){
	MBR cmbr;
	cmbr.pMin.x = max(polA.mbr.pMin.x, polB.mbr.pMin.x);
	cmbr.pMin.y = max(polA.mbr.pMin.y, polB.mbr.pMin.y);
	cmbr.pMax.x = min(polA.mbr.pMax.x, polB.mbr.pMax.x);
	cmbr.pMax.y = min(polA.mbr.pMax.y, polB.mbr.pMax.y);
	return cmbr;
}

/* checks if the exact geometry of 2 Objects intersects */
bool executePlaneSweep(Polygon &polA, Polygon &polB){
    //FIRST CHECK MBR CONTAINMENT
    // to perform the PiP test before the EI test (since its cheaper)  
    if(polB.mbr.pMin.x >= polA.mbr.pMin.x && polB.mbr.pMax.x <= polA.mbr.pMax.x && polB.mbr.pMin.y >= polA.mbr.pMin.y && polB.mbr.pMax.y <= polA.mbr.pMax.y){
    	//polA contains polB
    	//point in polygon test
		if(pointInPolygon(polB.vertices[0], polA)){
			return true;
		}
    }else if(polA.mbr.pMin.x >= polB.mbr.pMin.x && polA.mbr.pMax.x <= polB.mbr.pMax.x && polA.mbr.pMin.y >= polB.mbr.pMin.y && polA.mbr.pMax.y <= polB.mbr.pMax.y){
    	//polB contains polA
    	if(pointInPolygon(polA.vertices[0], polB)){
			return true;
		}
    }

    //move on to line segment intersection
    // Check if a Line Segment of polA intersects a Line Segment of polB.
    if (planeSweep(polA.exactGeometry_LS, polB.exactGeometry_LS)){
        return true;
    }    

	return false;
}

/* creates the line segment objects. we need them separately because plane sweep uses segment sorting */
void convertToLineSegments_andSort(Polygon &pol, MBR &cmbr){
    pol.exactGeometry_LS.reserve(pol.vertices.size()-1);

    //CMBR's bottom right and top left points
    Point bottomRight(cmbr.pMax.x, cmbr.pMin.y);
    Point topLeft(cmbr.pMin.x, cmbr.pMax.y);

    for(auto it = pol.vertices.begin(); it != pol.vertices.end()-1; it++){
    	LineSegment tempLS;
    	Point p1 = *it;
    	Point p2 = *(it+1);

    	//check relation to CMBR
    	//  if the edge intersects ANY edge of the CMBR, then keep it
    	//  else discard it
    	if(edgesIntersect(p1, p2, cmbr.pMin, bottomRight) || edgesIntersect(p1, p2, bottomRight, cmbr.pMax) || 
    		edgesIntersect(p1, p2, cmbr.pMax, topLeft) || edgesIntersect(p1, p2, topLeft, cmbr.pMin)){

    		if(p1.y < p2.y){
	            tempLS = {p1 ,p2, true};
	        }else if(p2.y < p1.y){
	            tempLS = {p2 ,p1, false};
	        }else{
	            if(p1.x < p2.x){
	                tempLS = {p1 ,p2, true};
	            }else{
	                tempLS = {p2 ,p1, false};
	            }
	        }
	        pol.exactGeometry_LS.push_back(tempLS);
    	}

    	
    }
    sort(pol.exactGeometry_LS.begin(), pol.exactGeometry_LS.end(), compareByYmin);
}

/* begin refinement */
bool refinePlaneSweep(Polygon &polA, Polygon &polB){
	//find CMBR first (optimizes refinement)
	MBR cmbr = getCMBR(polA,polB);

	//convert to line segments (required for plane sweep) and sort
	//  but keep ONLY those edges that intersect the CMBR
	convertToLineSegments_andSort(polA, cmbr);
	convertToLineSegments_andSort(polB, cmbr);
	
	return executePlaneSweep(polA, polB);
}