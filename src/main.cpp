#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>

#include <mpi.h>

#include "CReadData.hpp"

/* Cblacs declarations not declared any where in MKL (i don't understand this!)*/
extern "C" {
    void Cblacs_pinfo(int *rank, int *nprocs);
    void Cblacs_get(int context, int what, int *val);
    void Cblacs_gridinit(int *context, const char *layout, int proc_rows, int proc_cols);
    void Cblacs_pcoord(int context, int rank, int *row, int *col);
    void Cblacs_gridexit(int);
    void Cblacs_barrier(int context, const char *scope);
    void Cdgerv2d(int, int, int, double*, int, int, int);
    void Cdgesd2d(int, int, int, double*, int, int, int);

    /* Piece of shit does not work!
     * int numroc_(int *num_row_col, int *block_size, int *proc_coord, int *row_col_1_proc_coord, int *nprocs);
     */
}

int myNumRoC(int numRowsCols, int rowColBlockSize, int procRowColCoord, int procWithFirstRowCol, int nprocs)
{
    int extraBlocks, myDist, nBlocks, numLocalRowsCols;

    myDist = (nprocs + procRowColCoord - procWithFirstRowCol) % nprocs;
    nBlocks = numRowsCols / rowColBlockSize;
    
    numLocalRowsCols = (nBlocks/nprocs) * rowColBlockSize;
    extraBlocks = nBlocks % nprocs;
    
    if (myDist < extraBlocks)
        numLocalRowsCols += rowColBlockSize;
    else if (myDist == extraBlocks)
        numLocalRowsCols += (numRowsCols % rowColBlockSize);
    
    return numLocalRowsCols;
}

void setProcGrid(int *rank, int *nprocs, int *context, int numProcRows, int numProcCols,
                int *rankRow, int *rankCol)
{   
    Cblacs_pinfo(rank, nprocs);
    Cblacs_get(0, 0, context);
    Cblacs_gridinit(context, "Row-major", numProcRows, numProcCols);
    Cblacs_pcoord(*context, *rank, rankRow, rankCol);

    printf("Rank: %d, [%d, %d]\n", *rank, *rankRow, *rankCol);
    Cblacs_barrier(*context, "All");
}

void createSendDataTypes()
{
/* TODO:
    int dataBlockLengths = new int(nrows);
    int dataBlocksDisps = new int(nrows);
//    int MPI_Type_indexed(int count, const int *array_of_blocklengths, const int *array_of_displacements, MPI_Datatype oldtype, MPI_Datatype *newtype)
    for (int i =  0; i < nrows; i++)
    {
        dataBlockLengths[i] = ncols;
        dataBlockDisps[i] = myRankCol *  i * ncols;
    }
*/
}

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
    int blacsContext;
    int numProcs, myRank, myRankRow, myRankCol;
    std::vector< std::vector<int> > *data;

    // init a grid using blacs
    setProcGrid(&myRank, &numProcs, &blacsContext, gridProcRows, gridProcCols, &myRankRow, &myRankCol);
                
    // Number of rows and cols owned by the current process
    int nrows = myNumRoC(matRows, blockSize, myRankRow, root, gridProcRows);
    int ncols = myNumRoC(matCols, blockSize, myRankCol, root, gridProcCols);
    printf("rank: %d, [%d, %d]: rows: %d, cols: %d \n", myRank, myRankRow, myRankCol, nrows, ncols);

    std::vector<int> matrixData(nrows * ncols);

    // root reads all the data
    if (myRank == root)
    {
        CReadData readCSV(fileName, *delimiter.c_str());
        readCSV.readAllLines();
        //readCSV.printLines(numLines);
        data = readCSV.getData();

        // send data belonging to respective procs in row-major form
        MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
        //int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)

    }
    
    // distribute data to all procs in the grid
    
    MPI_Finalize();
    return 0;
}

//        printf("matrix: [%lu, %lu]\n", data->size(), data->at(0).size());
//        for (int i = 0; i < data->size(); i++)
//        {
//            for (int j = 0; j < data->at(i).size(); j++)
//            {
//                std::cout << data->at(i)[j] << std::setw(4);
//            }
//            std::cout << std::endl;
//        }
