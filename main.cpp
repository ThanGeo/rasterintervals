#include "def.h"
#include "getopt.h"
#include "./containers/relation.h"
#include "./algorithms/2d/fs.h"
#include "./partitioning/2d/fs.h"

#include "pipeline.h"

int main(int argc, char **argv)
{
    char c;
    int runNumPartitionsPerRelation = -1, runProcessingMethod = -1, runNumPartitions = -1, NUM_ITERATIONS = -1;
    bool runPlaneSweepOnX;
    Timer tim;
    double timeSorting = 0, timeIndexingOrPartitioning = 0, timeJoining = 0, timeCopying = 0;
    Relation R, S, *pR, *pS;
    Relation *pRA, *pSA, *pRB, *pSB, *pRC, *pSC, *pRD, *pSD;
    size_t *pRA_size, *pSA_size, *pRB_size, *pSB_size, *pRC_size, *pSC_size, *pRD_size, *pSD_size;
    size_t *pR_size ,*pS_size;
    vector<ABrec> *pRABdec , *pSABdec;
    vector<Crec> *pRCdec, *pSCdec;
    vector<Drec> *pRDdec, *pSDdec;
    vector<Coord> *pRYEND, *pSYEND;     

    clock_t start;

    unsigned long long result = 0;

    runPlaneSweepOnX = false;

    bool FLAGJOINITER  = false, FLAGSORTITER = false, FLAGPARTITER = false, FLAGCOP = false;
    double timeJoinIter = 0.0, timeSortIter = 0.0, timePartIter = 0.0, timeCopyingIter = 0.0;

    //get arguments
    string argument1(argv[argc-2]);
    string argument2(argv[argc-1]);

    //initialize
    initialize(argument1, argument2);

    //set options
    while ((c = getopt(argc, argv, "fcqp:?")) != -1)
    {
        switch (c)
        {
            case 'p':
                runNumPartitionsPerRelation = atoi(optarg);
                break;
            case 'c':
                CALCULATE_INTERVALS = 1;                
                break;
            case 'f':
                INTERMEDIATE_FILTER = 1;                
                break;
            case 'q':
                REFINEMENT = 1;
                break;
            case '?':
            case 'h':
            default:
                //usage();
                break;
        }
    }

    //calculate raster intervals, if selected
    if(CALCULATE_INTERVALS == 1){
        initiateRasterIntervalsCreation(argument1, argument2);
    }

    //enable the raster intervals filter, if selected
    if(INTERMEDIATE_FILTER == 1){
        enableIntermediateFilter(argument1, argument2);
    }

    cout << "Loading MBRs..." << endl;
    // Load inputs (creates MBRs from geometry files)
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            R.load(getBinaryGeometryFilename(argument1), universalMinX, universalMinY, universalMaxX, universalMaxY);
            cout << argument1 << ": " << R.size() << " polygons." << endl;
        }
        #pragma omp section
        {
            S.load(getBinaryGeometryFilename(argument2), universalMinX, universalMinY, universalMaxX, universalMaxY);
            cout << argument2 << ": " << S.size() << " polygons." << endl;
        }
    }

    //grid partitioning for the MBR filter (two layer partitioning)
    Coord minX = min(R.minX, S.minX);
    Coord maxX = max(R.maxX, S.maxX);
    Coord minY = min(R.minY, S.minY);
    Coord maxY = max(R.maxY, S.maxY);
    Coord diffX = maxX - minX;
    Coord diffY = maxY - minY;
    Coord maxExtend = (diffX<diffY)?diffY:diffX;

    R.normalize(minX, maxX, minY, maxY, maxExtend);
    S.normalize(minX, maxX, minY, maxY, maxExtend);    

    unsigned long long sizeDitt = 0 , sizeMinijoins = 0 , sizeA = 0 , sizeB = 0 , sizeC = 0 , sizeD = 0;

    runNumPartitions = runNumPartitionsPerRelation * runNumPartitionsPerRelation;

    //BEGIN MBR FILTER (MBR FILTER = FS LESS)
    runNumPartitions = runNumPartitionsPerRelation * runNumPartitionsPerRelation;
    pRA_size = new size_t[runNumPartitions];
    pRB_size = new size_t[runNumPartitions];
    pRC_size = new size_t[runNumPartitions];
    pRD_size = new size_t[runNumPartitions];

    pSA_size = new size_t[runNumPartitions];
    pSB_size = new size_t[runNumPartitions];
    pSC_size = new size_t[runNumPartitions];
    pSD_size = new size_t[runNumPartitions];

    memset(pRA_size, 0, runNumPartitions*sizeof(size_t));
    memset(pSA_size, 0, runNumPartitions*sizeof(size_t));
    memset(pRB_size, 0, runNumPartitions*sizeof(size_t));
    memset(pSB_size, 0, runNumPartitions*sizeof(size_t));
    memset(pRC_size, 0, runNumPartitions*sizeof(size_t));
    memset(pSC_size, 0, runNumPartitions*sizeof(size_t));
    memset(pRD_size, 0, runNumPartitions*sizeof(size_t));
    memset(pSD_size, 0, runNumPartitions*sizeof(size_t));

    pR = new Relation[runNumPartitions];
    pS = new Relation[runNumPartitions];
    
    //PRE-PROCESSING (partitioning and sorting)
    if (FLAGPARTITER == false){
        fs_2d::single::PartitionTwoDimensional(R, S, pR, pS, pRA_size, pSA_size, pRB_size, pSB_size, pRC_size, pSC_size, pRD_size, pSD_size,runPlaneSweepOnX, runNumPartitionsPerRelation);
        FLAGPARTITER = true;
    }
    else{
        tim.start();
        fs_2d::single::PartitionTwoDimensional(R, S, pR, pS, pRA_size, pSA_size, pRB_size, pSB_size, pRC_size, pSC_size, pRD_size, pSD_size,runPlaneSweepOnX, runNumPartitionsPerRelation);
        timeIndexingOrPartitioning = tim.stop();
        timePartIter += timeIndexingOrPartitioning;
    }

    if (FLAGSORTITER == false){
        fs_2d::single::sort::oneArray::SortYStartOneArray(pR, pS, pRB_size, pSB_size, pRC_size, pSC_size , pRD_size, pSD_size, runNumPartitions);
        FLAGSORTITER = true;
    }
    else{
        tim.start();
        fs_2d::single::sort::oneArray::SortYStartOneArray(pR, pS, pRB_size, pSB_size, pRC_size, pSC_size , pRD_size, pSD_size, runNumPartitions);
        timeSorting = tim.stop();
        timeSortIter += timeSorting;
    }
    
    cout << "Begin evaluation..." << endl;
    //EVALUATION
    start = clock();
    result = fs_2d::single::ForwardScanBased_PlaneSweep_CNT_Less(pR, pS, pRA_size, pSA_size, pRB_size, pSB_size, pRC_size, pSC_size, pRD_size, pSD_size, runPlaneSweepOnX, runNumPartitionsPerRelation);
    double evaluationTotalTime = (clock() - start) /(double) CLOCKS_PER_SEC;


    //print results & metrics
    cout << "Join in pipeline: " << endl;
    cout << "-----------------------------------" << endl;
    cout << "MBR Filter";
    if(INTERMEDIATE_FILTER){
        cout << " -> Raster Intervals Filter";
    }
    if(REFINEMENT){
        cout << " -> Refinement";
    }
    cout << " -> result" << endl;
    cout << "-----------------------------------" << endl;

    cout << "Finished in " << evaluationTotalTime << " seconds." << endl;
    cout << "   MBR Filter time: " << evaluationTotalTime - intermediateFilterTime - refinementTime << " seconds." << endl;    
    cout << "   Raster Intervals time: " << intermediateFilterTime << " seconds." << endl;
    cout << "   Refinement time: " << refinementTime << " seconds." << endl;    
    //cout << "   Geometry retrieval time: " << loadingGeometriesTime << " seconds." << endl;


    cout << "Candidate pairs analysis: " << endl;
    cout << "   Candidates after MBR Filter: " << postMBRCandidates << endl;
    cout << "   Accepted by RI: " << accepted / (double) postMBRCandidates * 100 << "%" << endl;
    cout << "   Rejected by RI: " << rejected / (double) postMBRCandidates * 100 << "%" << endl;
    cout << "   Inconclusive by RI: " << refinementCandidates / (double) postMBRCandidates * 100 << "%" << endl;
    //cout << "   Refinement Candidates: " << refinementCandidates << endl;
    //cout << "      Average # of vertices from candidates in " << argument1 << " : " << refinementCandidatesR / (double) (refinementCandidates) << " vertices." << endl;
    //cout << "      Average # of vertices from candidates in " << argument2 << " : " << refinementCandidatesS / (double) (refinementCandidates) << " vertices." << endl;
    cout << "   Accepted AFTER refinement (as a percentage of the refinement candiadtes): " << acceptedAfterRefinement / (double) postMBRCandidates * 100 << "%" << endl;
    
    if(INTERMEDIATE_FILTER == 0){
        //print the post MBR filter result
        cout << "Result: "<< result << " pairs." << endl;
        cout << "  post refinement: "<< acceptedAfterRefinement << " pairs." << endl;
    }else{
        //print the end of pipeline result.
        //    results may be slightly different from true result, due to the fact
        //      that our rasterization method fills polygon holes and openings.
        //    if a  100% accurate rasterization method is used, then this result will be the truth.
        cout << "Result: "<< TOTAL_RESULTS << " pairs." << endl;
    }


 
    return 0;
}
