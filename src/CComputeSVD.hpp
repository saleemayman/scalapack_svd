#ifndef CCOMPUTESVD_HPP
#define CCOMPUTESVD_HPP

#include <vector>
#include <algorithm>
#include <mkl_scalapack.h>

#include "common.h"

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
    int size, sizep, sizeq;
    char jobu, jobvt;
    MKL_INT info, ia, ja, iu, ju, ivt, jvt, lwork;
    std::vector<MKL_INT> descA;
    std::vector<MKL_INT> descU;
    std::vector<MKL_INT> descVT;

    // result variables
    std::vector<double> singularValues;
    std::vector<double> leftSingularVectors;    // ?
    std::vector<double> rightSingularVectors;   // ?
    std::vector<double> work; // ?

    void createArrayDescriptor(std::vector<MKL_INT> &descVec, int dtype, int ctxt, int m, int n, int mb, int nb, int rsrc, int csrc, int lld);
    void initSVDVariables();
public:
    CComputeSVD(gridInfo myGridInfo, 
                int totalRows, int totalCols, 
                int blockSizeRows, int blockSizeCols, 
                int gridNumProcRows, int gridNumProcCols, 
                int procWithFirstRow, int procWithFirstCol);
    ~CComputeSVD();

    void createLocal2DBlockCyclicMatrix(const std::vector<double> &coordData);
    void computeSVD();
    const std::vector<double>& getSingularValues() const;
    const std::vector<double>& getLeftSingularVectors() const;
    const std::vector<double>& getRightSingularVectors() const;
    void printLocalMatrix();
    void printLocalSingularValues();
    void printLocalLeftSingularVectors();
    void printLocalRightSingularVectors();
};
#endif
