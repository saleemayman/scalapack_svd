#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <vector>

#include <mpi.h>

#include "CReadData.hpp"
#include "CInitGrid.hpp"
#include "CComputeSVD.hpp"


int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    std::string fileName = argv[1];
    std::string delimiter = argv[2];
    int gridProcRows = atoi(argv[3]);
    int gridProcCols = atoi(argv[4]);
    int matRows = atoi(argv[5]);
    int matCols = atoi(argv[6]);
    int blockSizeRows = atoi(argv[7]);
    int blockSizeCols = atoi(argv[8]);

    int root = 0;

    CInitGrid myGrid(matRows, matCols, blockSizeRows, blockSizeCols, gridProcRows, gridProcCols, root, root);
    
    CComputeSVD temp(myGrid.getGridInfo(), 
                    matRows, matCols, blockSizeRows, blockSizeCols,
                    gridProcRows, gridProcCols,
                    root, root);

    // create block cyclic data for each proc
    {
        CReadData readCSV(fileName, *delimiter.c_str());
        readCSV.readAllLines();
        const std::vector<double> &data = readCSV.getData();
    
        temp.createLocal2DBlockCyclicMatrix(data);
    }        

    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < myGrid.getNumProcs(); i++)
    {
        if (i == myGrid.getMyRank())
        {
            temp.printLocalMatrix();
        } 
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // get the SVD
    temp.computeSVD();
    
    const std::vector<double> &singularValues = temp.getSingularValues(); 
//    const std::vector<double> &leftSingularVectors = temp.getLeftSingularVectors(); 
//    const std::vector<double> &rightSingularVectors = temp.getRightSingularVectors(); 

    // root will print the singular values
    if (myGrid.getMyRank() == root)
    {
        temp.printLocalSingularValues();
    } 
    MPI_Barrier(MPI_COMM_WORLD);

//    MPI_Barrier(MPI_COMM_WORLD);
//    for (int i = 0; i < myGrid.getNumProcs(); i++)
//    {
//        if (i == myGrid.getMyRank())
//        {
//            temp.printLocalLeftSingularVectors();
//        } 
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
//
//    MPI_Barrier(MPI_COMM_WORLD);
//    for (int i = 0; i < myGrid.getNumProcs(); i++)
//    {
//        if (i == myGrid.getMyRank())
//        {
//            temp.printLocalRightSingularVectors();
//        } 
//        MPI_Barrier(MPI_COMM_WORLD);
//    }


    MPI_Finalize();
    return 0;
}
