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
- point0 X coordinate (8 bytes)
- point0 Y coordinate (8 bytes)
- point1 X coordinate (8 bytes)
- point1 Y coordinate (8 bytes)
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

## Code

## Installation

## Execution
