#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include <mpi.h>
#include <mkl_scalapack.h>

#include "CReadData.hpp"

#define dtype_a 0
#define ctxt_a  1
#define m_a     2   
#define n_a     3
#define mb_a    4
#define nb_a    5
#define rsrc_a  6
#define csrc_a  7
#define lld_a   8

// Cblacs declarations not declared any where in MKL (i don't understand this!)
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

    //printf("Rank: %d, [%d, %d]\n", *rank, *rankRow, *rankCol);
    Cblacs_barrier(*context, "All");
}

void createSendDataTypes(int gridProcRows, int gridProcCols, int totalCols, int totalRows, int blockSize,
                        std::vector<MPI_Datatype> &sendDataType, 
                        std::vector<int> &numRows, std::vector<int> &numCols,
                        std::vector<int> &myRowDisp, std::vector<int> &myColDisp,
                        std::vector<int> &myRow, std::vector<int> &myCol)
{
    int rank, procRow, procCol, procRank;
    std::vector<int> myLocalIndex(gridProcRows * gridProcCols, 0);
    std::vector< std::vector<int> > sendDataBlockDisps;
    std::vector< std::vector<int> > sendDataBlockLengths;


    // get the global row and col displacements for each proc in the grid
    for (int gridRow = 0; gridRow < gridProcRows; gridRow++)
    {
        for (int gridCol = 0; gridCol < gridProcCols; gridCol++)
        {
            rank = gridCol + gridRow * gridProcCols;
            
            for (int i = 0; i < gridRow; i++)
                myRowDisp[rank] += numRows[i];
            
            for (int j = 0; j < gridCol; j++)
                myColDisp[rank] += numCols[j];

            //printf("rank: %d [%d, %d]: rowDisp: %d, colDisp: %d\n", rank, myRow[rank], myCol[rank], myRowDisp[rank], myColDisp[rank]);
        }
    }

    // init each procs send data type with the respective number of blocks
    for (int i = 0; i < gridProcRows*gridProcCols; i++)
    {
        sendDataBlockDisps.push_back(std::vector<int>(numRows[i] * numCols[i], 0));
        sendDataBlockLengths.push_back(std::vector<int>(numRows[i] * numCols[i], 1));
    }

    // add each block's displacement index in the disp vector in 2D block-cyclic manner.
    for (int j = 0; j < totalCols; j++)
    {
        for (int i = 0; i < totalRows; i++)
        {
            procRow = (i/blockSize) % gridProcRows;
            procCol = (j/blockSize) % gridProcCols;
            procRank = procCol + procRow * gridProcCols;

            sendDataBlockDisps[procRank][ myLocalIndex[procRank] ] = i * totalCols + j;
            myLocalIndex[procRank]++;
        }
    }

    for (int i = 0; i < gridProcRows*gridProcCols; i++)
    {
        MPI_Type_indexed(numRows[i] * numCols[i], sendDataBlockLengths[i].data(), sendDataBlockDisps[i].data(), MPI_DOUBLE, &sendDataType[i]);
        MPI_Type_commit(&sendDataType[i]);
    }

}


void initRootLocalData(int totalRows, int totalCols, int blockSize, int gridProcRows, int gridProcCols, const std::vector<double> &data, std::vector<double> &matrixData)
{
    int elems = 0;
    int procRow, procCol, procRank;
    for (int j = 0; j < totalCols; j++)
    {
        for (int i = 0; i < totalRows; i++)
        {
            procRow = (i/blockSize) % gridProcRows;
            procCol = (j/blockSize) % gridProcCols;
            procRank = procCol + procRow * gridProcCols;
            
            if (procRank == 0)
            {
                matrixData[elems] = data[j + i*totalCols];
                elems++;
            }
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
    int blacsContext;
    int numProcs, myRank, myRankRow, myRankCol;
    int myRows, myCols;
    std::vector<int> numRows;
    std::vector<int> numCols;
    std::vector<int> myRow;
    std::vector<int> myCol;
    std::vector<double> matrixData;
    std::vector<int> myRowDisp(gridProcRows * gridProcCols, 0);
    std::vector<int> myColDisp(gridProcRows * gridProcCols, 0);
    //const std::vector<int> &data;
    //std::vector<int> *data;

    // init a grid using blacs
    setProcGrid(&myRank, &numProcs, &blacsContext, gridProcRows, gridProcCols, &myRankRow, &myRankCol);
                
    // Number of rows and cols owned by the current process
    myRows = myNumRoC(matRows, blockSize, myRankRow, root, gridProcRows);
    myCols = myNumRoC(matCols, blockSize, myRankCol, root, gridProcCols);
    printf("rank: %d, [%d, %d]: rows: %d, cols: %d, blockSize: %d, elems: %d\n", myRank, myRankRow, myRankCol, myRows, myCols, blockSize, myRows*myCols);
    
    matrixData.resize(myRows * myCols, -1.0f);

    if (myRank == root)
    {
        numRows.resize(numProcs);
        numCols.resize(numProcs);
        myRow.resize(numProcs);
        myCol.resize(numProcs);
    }

    MPI_Gather(&myRows, 1, MPI_INT, numRows.data() + myRank, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Gather(&myCols, 1, MPI_INT, numCols.data() + myRank, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Gather(&myRankRow, 1, MPI_INT, myRow.data() + myRank, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Gather(&myRankCol, 1, MPI_INT, myCol.data() + myRank, 1, MPI_INT, root, MPI_COMM_WORLD);

//    if (myRank == root)
//    {
//        int l, m, pr, pc, x, y;
//        for (int i = 0; i < matRows; i++)
//        {
//            for (int j = 0; j < matCols; j++)
//            {
//                l = (i - 0)/(gridProcRows*blockSize);
//                m = (j - 0)/(gridProcCols*blockSize);
//                pr = ((i - 0)/blockSize) % gridProcRows;
//                pc = ((j - 0)/blockSize) % gridProcCols;
//                x = (i - 0) % blockSize + 0;
//                y = (j - 0) % blockSize + 0;
//                printf("[i, j]: [%d, %d], (l, m): (%d, %d), (pr, pc): (%d, %d), (x, y): (%d, %d)\n", i, j, l, m, pr, pc, x, y);
//            }
//        }
//    }


//    if (myRank == 0)
//    {
//        printf("numRows in root:\n");
//        for (int i = 0; i < numProcs; i++)
//            printf("numRows[%d]: %d\n", i, numRows[i]);
//    }

    // root reads all the data and distributes data to the rest
    if (myRank == root)
    {
        std::vector<MPI_Datatype> sendDataType(numProcs);
        createSendDataTypes(gridProcRows, gridProcCols, matCols, matRows, blockSize, sendDataType, numRows, numCols, myRowDisp, myColDisp, myRow, myCol);
        CReadData readCSV(fileName, *delimiter.c_str());

        readCSV.readAllLines();
        //readCSV.printLines(10);
        const std::vector<double> &data = readCSV.getData();
        printf("data[0]: %f\n", data[0]);

        // send data belonging to respective procs in row-major form
        for (int i = 1; i < numProcs; i++)
        {
            MPI_Send(data.data(), 1, sendDataType[i], i, 0, MPI_COMM_WORLD);
        }
        //int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)

        initRootLocalData(matRows, matCols, blockSize, gridProcRows, gridProcCols, data, matrixData);
    }
    else
    {
        MPI_Recv(matrixData.data(), myRows * myCols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
    }
    //printf("rank: %d, matData[0]: %f\n", myRank, matrixData[0]);

    MPI_Bcast(myRowDisp.data(), numProcs, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(myColDisp.data(), numProcs, MPI_INT, 0, MPI_COMM_WORLD);


    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < numProcs; i++)
    {
        if (i == myRank)
        {
            printf("Rank: %d, local A: \n", myRank);
            for (int i = 0; i < myRows ; i++)
            {
                for (int j = 0; j < myCols; j++)
                {
                    printf("  %f", matrixData[i + j*myRows]);
                }
                printf("\n");
            }
        } 
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // compute LU factorization using scalapack
    MKL_INT desca[9], descu[9], descvt[9];
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
    std::vector<double> leftSingularVectors(matRows * size);    // ?
    std::vector<double> rightSingularVectors(size * matCols);   // ?
    std::vector<double> work(myRows * blockSize); // ?

    // array descriptors
    desca[dtype_a] = 1;            descu[dtype_a] = 1;             descvt[dtype_a] = 1;
    desca[ctxt_a] = blacsContext;  descu[ctxt_a] = blacsContext;   descvt[ctxt_a] = blacsContext;
    desca[m_a] = matRows;          descu[m_a] = matRows;           descvt[m_a] = matCols;
    desca[n_a] = matCols;          descu[n_a] = matRows;           descvt[n_a] = matCols;
    desca[mb_a] = blockSize;       descu[mb_a] = blockSize;        descvt[mb_a] = blockSize;
    desca[nb_a] = blockSize;       descu[nb_a] = blockSize;        descvt[nb_a] = blockSize; 
    desca[rsrc_a] = 0;             descu[rsrc_a] = 0;              descvt[rsrc_a] = 0;
    desca[csrc_a] = 0;             descu[csrc_a] = 0;              descvt[csrc_a] = 0;
    desca[lld_a] = myRows;         descu[lld_a] = myRows;          descvt[lld_a] = myRows;

    //pdgetrf(&matRows, &matCols, matrixData.data(), &ia, &ja, desca, ipiv.data(), &info);
    pdgesvd(&jobu, &jobvt, &matRows, &matCols, matrixData.data(), &ia, &ja, desca, 
            singularValues.data(), leftSingularVectors.data(), &iu, &ju, descu, 
            rightSingularVectors.data(), &ivt, &jvt, descvt, 
            work.data(), &lwork, &info);
    if (info != 0)
        printf("rank: %d, info: %d\n", myRank, info);

    // re-allocate work using returned lwork and run SVD again
    lwork = work[0];
    work.resize(lwork);
    pdgesvd(&jobu, &jobvt, &matRows, &matCols, matrixData.data(), &ia, &ja, desca, 
            singularValues.data(), leftSingularVectors.data(), &iu, &ju, descu, 
            rightSingularVectors.data(), &ivt, &jvt, descvt, 
            work.data(), &lwork, &info);
    if (info != 0)
        printf("rank: %d, info: %d\n", myRank, info);

// void pdgetrf(MKL_INT *m, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, MKL_INT *ipiv, MKL_INT *info);
// void pdgesvd(char *jobu, char *jobvt, MKL_INT *m, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *s, double *u, MKL_INT *iu, MKL_INT *ju, MKL_INT *descu, double *vt, MKL_INT *ivt, MKL_INT *jvt, MKL_INT *descvt, double *work, MKL_INT *lwork, double *rwork, MKL_INT *info);

    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < numProcs; i++)
    {
        if (i == myRank)
        {
            printf("Rank: %d, lwork: %d, singular values: \n", myRank, (int)work[0]);
            for (int i = 0; i < size; i++)
            {
                printf("  %f", singularValues[i]);
            }
            printf("\n");
        } 
        MPI_Barrier(MPI_COMM_WORLD);
    }

//    MPI_Barrier(MPI_COMM_WORLD);
//    for (int i = 0; i < numProcs; i++)
//    {
//        if (i == myRank)
//        {
//            printf("Rank: %d, local A: \n", myRank);
//            for (int i = 0; i < myRows ; i++)
//            {
//                for (int j = 0; j < myCols; j++)
//                {
//                    printf("  %f", matrixData[i + j*myRows]);
//                }
//                printf("\n");
//            }
//        } 
//        MPI_Barrier(MPI_COMM_WORLD);
//    }

//    MPI_Barrier(MPI_COMM_WORLD);
//    for (int i = 0; i < numProcs; i++)
//    {
//        if (i == myRank)
//        {
//            printf("Rank: %d, local ipiv -> myRows: %d, blockSize: %d, size: %lu\n", myRank, myRows, blockSize, ipiv.size());
//            for (int i = 0; i < myRows*blockSize ; i++)
//            {
//                printf("  %d", ipiv[i]);
//            }
//            printf("\n");
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//    }


    MPI_Finalize();
    return 0;
}
