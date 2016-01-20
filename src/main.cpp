#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <ctime>

#include <mpi.h>

//#include "mrimonitor.h"
#include <scorep/SCOREP_User.h>

#include "CReadData.hpp"
#include "CInitGrid.hpp"
#include "CComputeSVD.hpp"


int main(int argc, char *argv[])
{
    int root = 0;
    clock_t t0, t1, t2, time;
    t0 = clock();

    MPI_Init(&argc, &argv);

    std::string fileName = argv[1];
    std::string delimiter = argv[2];
//    int gridProcRows = atoi(argv[3]);
//    int gridProcCols = atoi(argv[4]);
//    int matRows = atoi(argv[5]);
//    int matCols = atoi(argv[6]);
//    int blockSizeRows = atoi(argv[7]);
//    int blockSizeCols = atoi(argv[8]);
    int matRows = atoi(argv[3]);
    int matCols = atoi(argv[4]);
    int blockSizeRows = atoi(argv[5]);
    int blockSizeCols = atoi(argv[6]);

    t1 = clock();
    int num_procs;
    int gridProcRows;
    int gridProcCols = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
//    printf("procs: %d\n", num_procs);
    gridProcRows = num_procs;
    CInitGrid myGrid(matRows, matCols, blockSizeRows, blockSizeCols, gridProcRows, gridProcCols, root, root);
    CComputeSVD temp(myGrid.getGridInfo(), 
                    matRows, matCols, blockSizeRows, blockSizeCols,
                    gridProcRows, gridProcCols,
                    root, root);
    t2 = clock();
    
    time = t2 - t1;
    MPI_Reduce(&time, &t2, 1, MPI_LONG, MPI_SUM, root, MPI_COMM_WORLD);
    if (myGrid.getMyRank() == root)
        printf("rank: %d, init time: %f\n", myGrid.getMyRank(), (float)(t2)/(CLOCKS_PER_SEC*myGrid.getNumProcs()));

    t1 = clock();
    // create block cyclic data for each proc
    {
        CReadData readCSV(fileName, *delimiter.c_str());
        readCSV.readAllLines();
        const std::vector<double> &data = readCSV.getData();
    
        temp.createLocal2DBlockCyclicMatrix(data);
        MPI_Barrier(MPI_COMM_WORLD);
    }        
    t2 = clock();
    time = t2 - t1;
    MPI_Reduce(&time, &t2, 1, MPI_LONG, MPI_SUM, root, MPI_COMM_WORLD);
    if (myGrid.getMyRank() == root)
        printf("rank: %d, data dist. time: %f\n", myGrid.getMyRank(), (float)(t2)/(CLOCKS_PER_SEC*myGrid.getNumProcs()));

    // get the SVD
//    int psc_ret_val;
//    int pscOldTaskId = 0;
//    startMonLib();
//    startRegion(1, 1, 65, 0, -1);
    SCOREP_USER_REGION("ayman", SCOREP_USER_REGION_TYPE_COMMON);
    t1 = clock();
    temp.computeSVD();
    t2 = clock();
//    psc_ret_val = 0;
//    endRegion(1, 1, 65, 0, -1);
//    stopMonLib();


    time = t2 - t1;
    MPI_Reduce(&time, &t2, 1, MPI_LONG, MPI_SUM, root, MPI_COMM_WORLD);
    if (myGrid.getMyRank() == root)
        printf("rank: %d, SVD time: %f\n", myGrid.getMyRank(), (float)(t2)/(CLOCKS_PER_SEC*myGrid.getNumProcs()));

    const std::vector<double> &singularValues = temp.getSingularValues(); 
//    const std::vector<double> &leftSingularVectors = temp.getLeftSingularVectors(); 
//    const std::vector<double> &rightSingularVectors = temp.getRightSingularVectors(); 

    // root will print the singular values
//    if (myGrid.getMyRank() == root)
//    {
//        temp.printLocalSingularValues();
//    } 
//    MPI_Barrier(MPI_COMM_WORLD);


    t2 = clock();
    time = t2 - t0;
    MPI_Reduce(&time, &t2, 1, MPI_LONG, MPI_SUM, root, MPI_COMM_WORLD);
    if (myGrid.getMyRank() == root)
        printf("rank: %d, total time: %f\n", myGrid.getMyRank(), (float)(t2)/(CLOCKS_PER_SEC*myGrid.getNumProcs()));

    MPI_Finalize();
    return 0;

}
