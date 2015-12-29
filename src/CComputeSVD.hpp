#ifndef CCOMPUTESVD_HPP
#define CCOMPUTESVD_HPP

#include <vector>
#include <mkl_scalapack.h>

class CComputeSVD
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
    std::vector<double> *myData;

    // ScaLAPACK variables
    MKL_INT descA[9], descU[9], descVT[9];


public:
    CComputeSVD(int myRank, int numProcs, int context, 
                int myRankRow, int myRankCol,
                int myRows, int myCols,
                int totalRows, int totalCols, 
                int blockSizeRows, int blockSizeCols, 
                int gridNumProcRows, int gridNumProcCols, 
                int procWithFirstRow, int procWithFirstCol);
    ~CComputeSVD();

    void createLocal2DBlockCyclicMatrix(const std::vector<double> &coordData);
    void computeSVD();
    void getSingularValues();
    void getLeftSingularValues();
    void getRightSingularValues();
    void printLocalMatrix();
};
#endif
