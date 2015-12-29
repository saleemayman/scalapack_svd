#ifndef CINIT_GRID_HPP
#define CINIT_GRID_HPP

#include <stdio.h>
#include "common.h"

class CInitGrid
{
private:
    int myRank;
    int numProcs;
    int context;
    int gridNumProcRows, gridNumProcCols;
    int myRankRow, myRankCol;
    int totalRows, totalCols;
    int myRows, myCols;
    int blockSizeRows, blockSizeCols;
    int procWithFirstRow, procWithFirstCol;

    int myNumRoC(int numRowsCols, int rowColBlockSize, int procRowColCoord, int procWithFirstRowCol, int numRowColProcs);
    void setProcGrid();
    void setProcRowsOrCols();
public:
    CInitGrid(int totalRows, int totalCols, 
            int blockSizeRow, int blockSizeCol, 
            int gridNumProcRows, int gridNumProcCols, 
            int procWithFirstRow, int procWithFirstCol);
    ~CInitGrid();
};

#endif
