#ifndef CINIT_GRID_HPP
#define CINIT_GRID_HPP

#include <stdio.h>
#include "common.h"

class CInitGrid
{
//private:
protected:
    int myRank;
    int numProcs;
    int context;
    int gridNumProcRows, gridNumProcCols;
    int myRankRow, myRankCol;
    int totalRows, totalCols;
    int myRows, myCols;
    int blockSizeRows, blockSizeCols;
    int procWithFirstRow, procWithFirstCol;
    gridInfo myGridInfo;

    int myNumRoC(int numRowsCols, int rowColBlockSize, int procRowColCoord, int procWithFirstRowCol, int numRowColProcs);
    void setProcGrid();
    void setProcRowsOrCols();
public:
    CInitGrid(int totalRows, int totalCols, 
            int blockSizeRows, int blockSizeCols, 
            int gridNumProcRows, int gridNumProcCols, 
            int procWithFirstRow, int procWithFirstCol);
    ~CInitGrid();

    gridInfo getGridInfo();
    int getMyRank();
    int getNumProcs();
    int getContext();
    int getMyRow();
    int getMyCol();
    int getMyNumRows();
    int getMyNumCols();
};

#endif
