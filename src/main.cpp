#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include <mpi.h>
#include <mkl_scalapack.h>

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
    int blockSize = atoi(argv[7]);

    int root = 0;

    CInitGrid myGrid(matRows, matCols, blockSize, blockSize, gridProcRows, gridProcCols, root, root);

    CComputeSVD temp(myGrid.getMyRank(), 
                    myGrid.getNumProcs(), 
                    myGrid.getContext(), 
                    myGrid.getMyRow(),
                    myGrid.getMyCol(),
                    myGrid.getMyNumRows(),
                    myGrid.getMyNumCols(),
                    matRows, matCols, blockSize, blockSize,
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



//    for (int i = 0; i < numProcs; i++)
//    {
//        if (i == myRank)
//        {
//            printf("Rank: %d, lwork: %d, singular values: \n", myRank, (int)work[0]);
//            for (int i = 0; i < size; i++)
//            {
//                printf("  %f", singularValues[i]);
//            }
//            printf("\n");
//        } 
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
//
//    MPI_Barrier(MPI_COMM_WORLD);
//    for (int i = 0; i < numProcs; i++)
//    {
//        if (i == myRank)
//        {
//            printf("Rank: %d, leftSingularValues: \n", myRank);
//            for (int i = 0; i < myRows ; i++)
//            {
//                for (int j = 0; j < myCols; j++)
//                {
//                    printf("  %f", leftSingularVectors[i + j*myRows]);
//                }
//                printf("\n");
//            }
//        } 
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
//
//    MPI_Barrier(MPI_COMM_WORLD);
//    for (int i = 0; i < numProcs; i++)
//    {
//        if (i == myRank)
//        {
//            printf("Rank: %d, rightSingularValues: \n", myRank);
//            for (int i = 0; i < myRows ; i++)
//            {
//                for (int j = 0; j < myCols; j++)
//                {
//                    printf("  %f", rightSingularVectors[i + j*myRows]);
//                }
//                printf("\n");
//            }
//        } 
//        MPI_Barrier(MPI_COMM_WORLD);
//    }


    MPI_Finalize();
    return 0;
}
