#include <stdio.h>

#include "CComputeSVD.hpp"

CComputeSVD::CComputeSVD(int myRank, int numProcs, int context, 
                        int myRankRow, int myRankCol,
                        int myRows, int myCols,
                        int totalRows, int totalCols, 
                        int blockSizeRows, int blockSizeCols, 
                        int gridNumProcRows, int gridNumProcCols, 
                        int procWithFirstRow, int procWithFirstCol):
                            myRank(myRank), numProcs(numProcs), context(context),
                            myRankRow(myRankRow), myRankCol(myRankCol),
                            myRows(myRows), myCols(myCols),
                            totalRows(totalRows), totalCols(totalCols),
                            blockSizeRows(blockSizeRows), blockSizeCols(blockSizeCols),
                            gridNumProcRows(gridNumProcRows), gridNumProcCols(gridNumProcCols),
                            procWithFirstRow(procWithFirstRow), procWithFirstCol(procWithFirstCol)
{
    printf("CComputeSVD -> rank: %d, [%d, %d]: rows: %d, cols: %d, blockSizeRows: %d, blockSizeCols: %d, elems: %d\n", myRank, myRankRow, myRankCol, myRows, myCols, blockSizeRows, blockSizeCols, myRows*myCols);

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

void CComputeSVD::createArrayDescriptors()
{
    // array descriptors
    descA[dtype_a] = 1;            descU[dtype_a] = 1;             descVT[dtype_a] = 1;
    descA[ctxt_a] = context;       descU[ctxt_a] = context;        descVT[ctxt_a] = context;
    descA[m_a] = matRows;          descU[m_a] = matRows;           descVT[m_a] = matCols;
    descA[n_a] = matCols;          descU[n_a] = matRows;           descVT[n_a] = matCols;
    descA[mb_a] = blockSize;       descU[mb_a] = blockSize;        descVT[mb_a] = blockSize;
    descA[nb_a] = blockSize;       descU[nb_a] = blockSize;        descVT[nb_a] = blockSize; 
    descA[rsrc_a] = 0;             descU[rsrc_a] = 0;              descVT[rsrc_a] = 0;
    descA[csrc_a] = 0;             descU[csrc_a] = 0;              descVT[csrc_a] = 0;
    descA[lld_a] = myRows;         descU[lld_a] = myRows;          descVT[lld_a] = myRows;
}

void CComputeSVD::computeSVD()
{
    // compute LU factorization using scalapack
    MKL_INT info;
    MKL_INT ia = 1;
    MKL_INT ja = 1;
    MKL_INT iu = 1;
    MKL_INT ju = 1;
    MKL_INT ivt = 1;
    MKL_INT jvt = 1;
    MKL_INT lwork = -1; //myRows * blockSize; // ?
    int size = std::min(matRows, matCols);
    int sizeq = myRows;
    int sizep = myCols;
    char jobu = 'V';
    char jobvt = 'V';
    std::vector<double> singularValues(size);
    std::vector<double> leftSingularVectors(myRows * myCols);    // ?
    std::vector<double> rightSingularVectors(myRows * myCols);   // ?
    std::vector<double> work(myRows * blockSize); // ?


    pdgesvd(&jobu, &jobvt, &matRows, &matCols, myData->data(), &ia, &ja, descA, 
            singularValues.data(), leftSingularVectors.data(), &iu, &ju, descU, 
            rightSingularVectors.data(), &ivt, &jvt, descVT, 
            work.data(), &lwork, &info);
    if (info != 0)
        printf("rank: %d, info: %d\n", myRank, info);

    // re-allocate work using returned lwork and run SVD again
    lwork = work[0];
    work.resize(lwork);
    pdgesvd(&jobu, &jobvt, &matRows, &matCols, myData->data(), &ia, &ja, descA, 
            singularValues.data(), leftSingularVectors.data(), &iu, &ju, descU, 
            rightSingularVectors.data(), &ivt, &jvt, descVT, 
            work.data(), &lwork, &info);
    if (info != 0)
        printf("rank: %d, info: %d\n", myRank, info);

// void pdgesvd(char *jobu, char *jobvt, MKL_INT *m, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *s, double *u, MKL_INT *iu, MKL_INT *ju, MKL_INT *descu, double *vt, MKL_INT *ivt, MKL_INT *jvt, MKL_INT *descvt, double *work, MKL_INT *lwork, double *rwork, MKL_INT *info);

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



