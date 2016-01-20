#include <stdio.h>

#include "CComputeSVD.hpp"

CComputeSVD::CComputeSVD(gridInfo myGridInfo,
                        int totalRows, int totalCols, 
                        int blockSizeRows, int blockSizeCols, 
                        int gridNumProcRows, int gridNumProcCols, 
                        int procWithFirstRow, int procWithFirstCol):
                            myRank(myGridInfo.myRank), numProcs(myGridInfo.numProcs), context(myGridInfo.context),
                            myRankRow(myGridInfo.myRow), myRankCol(myGridInfo.myCol),
                            myRows(myGridInfo.myNumRows), myCols(myGridInfo.myNumCols),
                            totalRows(totalRows), totalCols(totalCols),
                            blockSizeRows(blockSizeRows), blockSizeCols(blockSizeCols),
                            gridNumProcRows(gridNumProcRows), gridNumProcCols(gridNumProcCols),
                            procWithFirstRow(procWithFirstRow), procWithFirstCol(procWithFirstCol),
                            descA(9), descU(9), descVT(9),
                            singularValues(std::min(totalRows, totalCols)),
                            leftSingularVectors(myRows * myCols),
                            rightSingularVectors(myRows * myCols),
                            work(myRows * blockSizeRows)
{
    //printf("CComputeSVD -> rank: %d, [%d, %d]: rows: %d, cols: %d, blockSizeRows: %d, blockSizeCols: %d, elems: %d\n", myRank, myRankRow, myRankCol, myRows, myCols, blockSizeRows, blockSizeCols, myRows*myCols);

    // initialize the variables for the SVD routine
    ia = 1;
    ja = 1;
    iu = 1;
    ju = 1;
    ivt = 1;
    jvt = 1;
    lwork = -1; //myRows * blockSize; // ?
    size = std::min(totalRows, totalCols);
    sizeq = myRows;
    sizep = myCols;
    jobu = 'V';
    jobvt = 'V';
//    singularValues->resize(size);
//    leftSingularVectors->resize(myRows * myCols);    // ?
//    rightSingularVectors->resize(myRows * myCols);   // ?
//    work->resize(myRows * blockSizeRows); // ?

    // initalize the local matrices
    myData = new std::vector<double>(myRows * myCols, 0);
}

CComputeSVD::~CComputeSVD()
{
    delete myData;
}

void CComputeSVD::createLocal2DBlockCyclicMatrix(const std::vector<double> &coordData)
{
    int procRow, procCol, procRank;
    int blockRow, blockCol, linearDisp;
    int localRow, localCol;

    int row, col;
    double value;
    for (int i = 0; i < coordData.size(); i+=3)
    {
        row = (int)coordData[i + 0] - 1;
        col = (int)coordData[i + 1] - 1;
        value = coordData[i + 2];

        procRow = (row/blockSizeRows) % gridNumProcRows;
        procCol = (col/blockSizeCols) % gridNumProcCols;
        procRank = procCol + procRow * gridNumProcCols;
 
        if (procRank == myRank)
        {
            // block coordinate and the coordinates of "value" in the block
            blockRow = row/(gridNumProcRows * blockSizeRows);
            blockCol = col/(gridNumProcCols * blockSizeCols);
            localRow = row % (blockSizeRows + 0);
            localCol = col % (blockSizeCols + 0);
            linearDisp = localRow + localCol * myRows + blockCol * blockSizeCols * myRows + blockRow * blockSizeRows;
            myData->operator[](linearDisp) = value;
            //printf("rank: %d, val: %d [%d, %d], block: [%d, %d], local: [%d, %d], idx: %d\n", myRank, (int)value, row, col, blockRow, blockCol, localRow, localCol, linearDisp);
        }
    }
}

void CComputeSVD::createArrayDescriptor(std::vector<MKL_INT> &descVec, int dtype, int ctxt, int m, int n, int mb, int nb, int rsrc, int csrc, int lld)
{
    // array descriptors
    descVec[DTYPE_] = dtype;
    descVec[CTXT_] = ctxt;
    descVec[M_] = m;
    descVec[N_] = n;
    descVec[MB_] = mb;
    descVec[NB_] = nb;
    descVec[RSRC_] = rsrc;
    descVec[CSRC_] = csrc;
    descVec[LLD_] = lld;
}


void CComputeSVD::initSVDVariables()
{
    createArrayDescriptor(descA, 1, context, totalRows, totalCols, blockSizeRows, blockSizeCols, 0, 0, myRows);
    createArrayDescriptor(descU, 1, context, totalRows, totalRows, blockSizeRows, blockSizeCols, 0, 0, myRows);
    createArrayDescriptor(descVT, 1, context, totalCols, totalCols, blockSizeRows, blockSizeCols, 0, 0, myRows);
}


void CComputeSVD::computeSVD()
{
    initSVDVariables();

    pdgesvd(&jobu, &jobvt, &totalRows, &totalCols, myData->data(), &ia, &ja, descA.data(), 
            singularValues.data(), leftSingularVectors.data(), &iu, &ju, descU.data(), 
            rightSingularVectors.data(), &ivt, &jvt, descVT.data(), 
            work.data(), &lwork, &info);
    if (info != 0)
        printf("rank: %d, info: %d\n", myRank, info);

    if (myRank == 0)
        printf("rank: %d, about to start SVD comp. ...\n", myRank);

    // re-allocate work using returned lwork and run SVD again
    //lwork = work->operator[](0);
    lwork = work[0];
    work.resize(lwork);
    pdgesvd(&jobu, &jobvt, &totalRows, &totalCols, myData->data(), &ia, &ja, descA.data(), 
            singularValues.data(), leftSingularVectors.data(), &iu, &ju, descU.data(), 
            rightSingularVectors.data(), &ivt, &jvt, descVT.data(), 
            work.data(), &lwork, &info);
    if (info != 0)
        printf("rank: %d, info: %d\n", myRank, info);

// void pdgesvd(char *jobu, char *jobvt, MKL_INT *m, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *s, double *u, MKL_INT *iu, MKL_INT *ju, MKL_INT *descu, double *vt, MKL_INT *ivt, MKL_INT *jvt, MKL_INT *descvt, double *work, MKL_INT *lwork, double *rwork, MKL_INT *info);

}

const std::vector<double>& CComputeSVD::getSingularValues() const
{
    return singularValues;
}

const std::vector<double>& CComputeSVD::getLeftSingularVectors() const
{
    return leftSingularVectors;
}

const std::vector<double>& CComputeSVD::getRightSingularVectors() const
{   
    return rightSingularVectors;
}

void CComputeSVD::printLocalSingularValues()
{
    printf("Rank: %d, lwork: %d, singular values: \n", myRank, (int)work[0]);
    for (int i = 0; i < size; i++)
    {
        printf("%f\n", singularValues[i]);
    }
    printf("\n");
}

void CComputeSVD::printLocalLeftSingularVectors()
{
    printf("Rank: %d, leftSingularValues: \n", myRank);
    for (int i = 0; i < myRows ; i++)
    {
        for (int j = 0; j < myCols; j++)
        {
            printf("  %f", leftSingularVectors[i + j*myRows]);
        }
        printf("\n");
    }
}

void CComputeSVD::printLocalRightSingularVectors()
{
    printf("Rank: %d, rightSingularValues: \n", myRank);
    for (int i = 0; i < myRows ; i++)
    {
        for (int j = 0; j < myCols; j++)
        {
            printf("  %f", rightSingularVectors[i + j*myRows]);
        }
        printf("\n");
    }
}

void CComputeSVD::printLocalMatrix()
{
    printf("Rank: %d, local A:\n", myRank);
    for (int j = 0; j < myRows; j++)
    {
        for (int i = 0; i < myCols; i++)
        {
            printf("  %d", (int)myData->operator[](j + i * myRows));
        }
        printf("\n");
    }
    printf("\n");
}



