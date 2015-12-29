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

void createLocal2DBlockCyclicMatrix(int myRank, int gridProcRows, int gridProcCols, int blockSize,
                        const std::vector<double> &coordData, std::vector<double> &myData, int myRows, int myCols)
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

        procRow = (row/blockSize) % gridProcRows;
        procCol = (col/blockSize) % gridProcCols;
        procRank = procCol + procRow * gridProcCols;
 
        if (procRank == myRank)
        {
            // block coordinate and the coordinates of "value" in the block
            blockRow = row/(gridProcRows * blockSize);
            blockCol = col/(gridProcCols * blockSize);
            localRow = row % (blockSize + 0);
            localCol = col % (blockSize + 0);
            linearDisp = localRow + localCol * myRows + blockCol * blockSize * myRows + blockRow * blockSize;
            myData[linearDisp] = value;
            //printf("recvData type -> rank: %d, val: %d [%d, %d], block: [%d, %d], local: [%d, %d], idx: %d\n", myRank, (int)value, row, col, blockRow, blockCol, localRow, localCol, linearDisp);
        }
    }
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
//    int blacsContext;
//    int numProcs, myRank, myRankRow, myRankCol;
//    int myRows, myCols;
//    std::vector<double> matrixData;
//    std::vector<int> myRowDisp(gridProcRows * gridProcCols, 0);
//    std::vector<int> myColDisp(gridProcRows * gridProcCols, 0);

    CInitGrid myGrid(matRows, matCols, blockSize, blockSize, gridProcRows, gridProcCols, root, root);

//    CInitGrid(int numProcRows, int numProcCols, 
//            int blockSizeRow, int blockSizeCol, 
//            int gridNumProcRows, int gridNumProcCols, 
//            int procWithFirstRow, int procWithFirstCol);

//    // init a grid using blacs
//    setProcGrid(&myRank, &numProcs, &blacsContext, gridProcRows, gridProcCols, &myRankRow, &myRankCol);
//                
//    // Number of rows and cols owned by the current process
//    myRows = myNumRoC(matRows, blockSize, myRankRow, root, gridProcRows);
//    myCols = myNumRoC(matCols, blockSize, myRankCol, root, gridProcCols);
//    printf("rank: %d, [%d, %d]: rows: %d, cols: %d, blockSize: %d, elems: %d\n", myRank, myRankRow, myRankCol, myRows, myCols, blockSize, myRows*myCols);
//    
//    matrixData.resize(myRows * myCols, 0.0f);
//
//    // create block cyclic data for each proc
//    {
//        CReadData readCSV(fileName, *delimiter.c_str());
//        readCSV.readAllLines();
//        const std::vector<double> &data = readCSV.getData();
//        
//        createLocal2DBlockCyclicMatrix(myRank, gridProcRows, gridProcCols, blockSize, data, matrixData, myRows, myCols);
//    }
//
//    MPI_Barrier(MPI_COMM_WORLD);
//    for (int i = 0; i < numProcs; i++)
//    {
//        if (i == myRank)
//        {
//            printf("Rank: %d, local A size: %lu\n", myRank, matrixData.size());
//            for (int j = 0; j < myRows; j++)
//            {
//                for (int i = 0; i < myCols; i++)
//                {
//                    printf("  %d", (int)matrixData[j + i * myRows]);
//                }
//                printf("\n");
//            }
//            printf("\n");
//        } 
//        MPI_Barrier(MPI_COMM_WORLD);
//    }



//    // compute LU factorization using scalapack
//    MKL_INT desca[9], descu[9], descvt[9];
//    MKL_INT info;
//    MKL_INT ia = 1;
//    MKL_INT ja = 1;
//    MKL_INT iu = 1;
//    MKL_INT ju = 1;
//    MKL_INT ivt = 1;
//    MKL_INT jvt = 1;
//    MKL_INT lwork = -1; //myRows * blockSize; // ?
//    int size = std::min(matRows, matCols);
//    int sizeq = myRows;
//    int sizep = myCols;
//    char jobu = 'V';
//    char jobvt = 'V';
//    std::vector<double> singularValues(size);
//    std::vector<double> leftSingularVectors(myRows * myCols);    // ?
//    std::vector<double> rightSingularVectors(myRows * myCols);   // ?
//    std::vector<double> work(myRows * blockSize); // ?
//
//    // array descriptors
//    desca[dtype_a] = 1;            descu[dtype_a] = 1;             descvt[dtype_a] = 1;
//    desca[ctxt_a] = blacsContext;  descu[ctxt_a] = blacsContext;   descvt[ctxt_a] = blacsContext;
//    desca[m_a] = matRows;          descu[m_a] = matRows;           descvt[m_a] = matCols;
//    desca[n_a] = matCols;          descu[n_a] = matRows;           descvt[n_a] = matCols;
//    desca[mb_a] = blockSize;       descu[mb_a] = blockSize;        descvt[mb_a] = blockSize;
//    desca[nb_a] = blockSize;       descu[nb_a] = blockSize;        descvt[nb_a] = blockSize; 
//    desca[rsrc_a] = 0;             descu[rsrc_a] = 0;              descvt[rsrc_a] = 0;
//    desca[csrc_a] = 0;             descu[csrc_a] = 0;              descvt[csrc_a] = 0;
//    desca[lld_a] = myRows;         descu[lld_a] = myRows;          descvt[lld_a] = myRows;
//
//    //pdgetrf(&matRows, &matCols, matrixData.data(), &ia, &ja, desca, ipiv.data(), &info);
//    pdgesvd(&jobu, &jobvt, &matRows, &matCols, matrixData.data(), &ia, &ja, desca, 
//            singularValues.data(), leftSingularVectors.data(), &iu, &ju, descu, 
//            rightSingularVectors.data(), &ivt, &jvt, descvt, 
//            work.data(), &lwork, &info);
//    if (info != 0)
//        printf("rank: %d, info: %d\n", myRank, info);
//
//    // re-allocate work using returned lwork and run SVD again
//    lwork = work[0];
//    work.resize(lwork);
//    pdgesvd(&jobu, &jobvt, &matRows, &matCols, matrixData.data(), &ia, &ja, desca, 
//            singularValues.data(), leftSingularVectors.data(), &iu, &ju, descu, 
//            rightSingularVectors.data(), &ivt, &jvt, descvt, 
//            work.data(), &lwork, &info);
//    if (info != 0)
//        printf("rank: %d, info: %d\n", myRank, info);
//
//// void pdgesvd(char *jobu, char *jobvt, MKL_INT *m, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *s, double *u, MKL_INT *iu, MKL_INT *ju, MKL_INT *descu, double *vt, MKL_INT *ivt, MKL_INT *jvt, MKL_INT *descvt, double *work, MKL_INT *lwork, double *rwork, MKL_INT *info);
//
//    MPI_Barrier(MPI_COMM_WORLD);
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
