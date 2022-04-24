# Raster Intervals
Code for paper-id 140 "Raster Intervals"

will be uploaded here by 29.04.2022


## Description
This is an end-to-end implementation of a join operation between two polygon data sets. The MBR filter *CITE TSITSIGKOS?* generates the candidate pairs which move to the Raster Intervals intermediate filter. Their approximations are used there to identify true hits or true negatives and any indecisive pair cases are forwarded to the refinement step (plane sweep). There are options to i) generate the Raster Intervals approximations for the datasets ii) load existing approximations from disk to avoid the creation process and iii) skip the intermediate filter and/or refinement step in the pipeline.

## Data

This implementation needs the geometry input data to be in a binary form using the following specific format:
- total polygon count (4 bytes)
- polygon0 ID (4 bytes)
- polygon0 vertex count (4 bytes)
- vertex0 X coordinate (8 bytes)
- vertex0 Y coordinate (8 bytes)
- vertex1 X coordinate (8 bytes)
- vertex1 Y coordinate (8 bytes)
- ... 
- polygon1 ID (4 bytes)
- etc...

Each geometry data set is accompanied by a byte offset map that can be used for faster data retrieval from disk. This offset map pairs a polygon ID with a byte offset that is used through seekg() to "jump" to the exact spot in the file and read the requested polygon without having to iterate through the rest of the information. This offset map is a binary file on disk in the following format:
- total polygon count (4 bytes)
- polygon0 ID (4 bytes)
- polygon0 byte offset (8 bytes)
- polygon1 ID (4 bytes)
- etc...

The Raster Intervals are generated through code and stored in the rasterintervals/interval_data/ directory in binary form, using the following format:
- total polygon count (4 bytes)
- polygon0 ID (4 bytes)
- polygon0 number of total coding bytes T (4 bytes)
- polygon0 RI coding data (T bytes)
- polygon1 ID (4 bytes)
- etc...

## Structure

### MBR Filter

Directoris algorithms/, containers/, grid/ and partitioning/ contain code related to the MBR Filter used and is of no dependance to the Raster Intervals, with the exception of the function Relation::load() in containers/relation.cpp that creates the MBRs from the polygon geometries. Also, when a candidate pair is identified by the MBR Filter, it is forwarded further down the pipeline with the forwardCandidatePair() method that is invoked in the algorithms/fs.h file and implemented in the pipeline.h file.

### Pipeline

To connect the elements of the pipeline with each other, we use a dedicated set of methods that interconnects the individual components. File pipeline.h implements these methods. More specifically, forwardCandidatePairs() forwards a candidate pair coming right after the MBR Filter to the Raster Intervals intermediate filter and then to the Refinement step if needed.

### Raster Intervals, Intermediate Filter & Refinement

All code related to the Raster Intervals paper is stored in the rasterintervals/ directory. The creation process can be found in the raster_intervals.h and rasterization.h files, while the code related to the filter and refinement in the join.h and join_geometry_refinement.h files.

### Main

In the main directory, main.cpp implements driver code for the MBR Filter. After the pipeline has finished for all candidate pairs, main.cpp prints relevant results that were gathered during the process.

## Installation

To create the executable, use the 'make' command in the main directory. The code has been tested exclusively in Ubuntu 20.04.1 using C++14, g++ 9.4.0 and the -O3 optimizer. OpenMP is used to speed up the geometry loading process when creating the MBRs. 

## Execution
