#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <ctime>

#include <mpi.h>

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
    int gridProcRows = atoi(argv[3]);
    int gridProcCols = atoi(argv[4]);
    int matRows = atoi(argv[5]);
    int matCols = atoi(argv[6]);
    int blockSizeRows = atoi(argv[7]);
    int blockSizeCols = atoi(argv[8]);

    t1 = clock();
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
    }        
    t2 = clock();
    time = t2 - t1;
    MPI_Reduce(&time, &t2, 1, MPI_LONG, MPI_SUM, root, MPI_COMM_WORLD);
    if (myGrid.getMyRank() == root)
        printf("rank: %d, data dist. time: %f\n", myGrid.getMyRank(), (float)(t2)/(CLOCKS_PER_SEC*myGrid.getNumProcs()));

    // get the SVD
    t1 = clock();
    temp.computeSVD();
    t2 = clock();

    time = t2 - t1;
    MPI_Reduce(&time, &t2, 1, MPI_LONG, MPI_SUM, root, MPI_COMM_WORLD);
    if (myGrid.getMyRank() == root)
        printf("rank: %d, SVD time: %f\n", myGrid.getMyRank(), (float)(t2)/(CLOCKS_PER_SEC*myGrid.getNumProcs()));

    const std::vector<double> &singularValues = temp.getSingularValues(); 
//    const std::vector<double> &leftSingularVectors = temp.getLeftSingularVectors(); 
//    const std::vector<double> &rightSingularVectors = temp.getRightSingularVectors(); 

    // root will print the singular values
    if (myGrid.getMyRank() == root)
    {
        temp.printLocalSingularValues();
    } 
    MPI_Barrier(MPI_COMM_WORLD);


    t2 = clock();
    time = t2 - t0;
    MPI_Reduce(&time, &t2, 1, MPI_LONG, MPI_SUM, root, MPI_COMM_WORLD);
    if (myGrid.getMyRank() == root)
        printf("rank: %d, total time: %f\n", myGrid.getMyRank(), (float)(t2)/(CLOCKS_PER_SEC*myGrid.getNumProcs()));

    MPI_Finalize();
    return 0;
}
