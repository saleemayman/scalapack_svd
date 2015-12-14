#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>

#include <mpi.h>
//#include "mkl_blacs.h"
//#include <mkl_scalapack.h>

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

int main(int argc, char *argv[])
{
    int gridProcRows = atoi(argv[1]);
    int gridProcCols = atoi(argv[2]);
    int matRows = atoi(argv[3]);
    int matCols = atoi(argv[4]);
    int blockSize = atoi(argv[5]);

    int root = 0;
    int numProcs, myRank;
    int myRankRow, myRankCol;
    int blacsContext;
    int iZERO = 0;

    MPI_Init(&argc, &argv);

    // init a grid using blacs
    setProcGrid(&myRank, &numProcs, &blacsContext, gridProcRows, gridProcCols, &myRankRow, &myRankCol);

    // Number of rows and cols owned by the current process
    int nrows = myNumRoC(matRows, blockSize, myRankRow, root, gridProcRows);
    int ncols = myNumRoC(matCols, blockSize, myRankCol, root, gridProcCols);
    printf("rank: %d, [%d, %d]: rows: %d, cols: %d \n", myRank, myRankRow, myRankCol, nrows, ncols);

    // distribute data to all procs in the grid


    MPI_Finalize();
    return 0;
}
