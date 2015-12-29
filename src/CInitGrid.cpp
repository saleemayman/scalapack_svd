#include "CInitGrid.hpp"

CInitGrid::CInitGrid(int totalRows, int totalCols,
                     int blockSizeRows, int blockSizeCols,
                     int gridNumProcRows, int gridNumProcCols,
                     int procWithFirstRow, int procWithFirstCol):
                        totalRows(totalRows),
                        totalCols(totalCols),
                        blockSizeRows(blockSizeRows),
                        blockSizeCols(blockSizeCols),
                        gridNumProcRows(gridNumProcRows),
                        gridNumProcCols(gridNumProcCols),
                        procWithFirstRow(procWithFirstRow),
                        procWithFirstCol(procWithFirstCol)
{
    // initiaize the rectangular grid
    setProcGrid();

    // allocate number of rows and cols for each proc
    setProcRowsOrCols();
}

CInitGrid::~CInitGrid()
{
}

int CInitGrid::myNumRoC(int numRowsCols, int rowColBlockSize, int procRowColCoord, int procWithFirstRowCol, int numRowColProcs)
{
    int extraBlocks, myDist, nBlocks, numLocalRowsCols;

    myDist = (numRowColProcs + procRowColCoord - procWithFirstRowCol) % numRowColProcs;
    nBlocks = numRowsCols / rowColBlockSize;
    
    numLocalRowsCols = (nBlocks/numRowColProcs) * rowColBlockSize;
    extraBlocks = nBlocks % numRowColProcs;
    
    if (myDist < extraBlocks)
        numLocalRowsCols += rowColBlockSize;
    else if (myDist == extraBlocks)
        numLocalRowsCols += (numRowsCols % rowColBlockSize);

    return numLocalRowsCols;
}

void CInitGrid::setProcGrid()
{   
    Cblacs_pinfo(&myRank, &numProcs);
    Cblacs_get(0, 0, &context);
    Cblacs_gridinit(&context, "Row-major", gridNumProcRows, gridNumProcCols);
    Cblacs_pcoord(context, myRank, &myRankRow, &myRankCol);

    //printf("Rank: %d, [%d, %d]\n", *myRank, rankRow, *rankCol);
    Cblacs_barrier(context, "All");
}
    
void CInitGrid::setProcRowsOrCols()
{
    myRows = myNumRoC(totalRows, blockSizeRows, myRankRow, procWithFirstRow, gridNumProcRows);
    myCols = myNumRoC(totalCols, blockSizeCols, myRankCol, procWithFirstCol, gridNumProcCols);
    printf("CInitGrid -> rank: %d, [%d, %d]: rows: %d, cols: %d, blockSizeRows: %d, blockSizeCols: %d, elems: %d\n", myRank, myRankRow, myRankCol, myRows, myCols, blockSizeRows, blockSizeCols, myRows*myCols);
}


int CInitGrid::getMyRank()
{
    return myRank;
}

int CInitGrid::getNumProcs()
{
    return numProcs;
}

int CInitGrid::getContext()
{
    return context;
}

int CInitGrid::getMyRow()
{
    return myRankRow;
}

int CInitGrid::getMyCol()
{
    return myRankCol;
}

int CInitGrid::getMyNumRows()
{
    return myRows;
}

int CInitGrid::getMyNumCols()
{
    return myCols;
}










