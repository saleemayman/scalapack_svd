#include <mpi.h>
 
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
 
using namespace std;
 
extern "C" {
    /* Cblacs declarations */
    void Cblacs_pinfo(int*, int*);
    void Cblacs_get(int, int, int*);
    void Cblacs_gridinit(int*, const char*, int, int);
    void Cblacs_pcoord(int, int, int*, int*);
    void Cblacs_gridexit(int);
    void Cblacs_barrier(int, const char*);
    void Cdgerv2d(int, int, int, double*, int, int, int);
    void Cdgesd2d(int, int, int, double*, int, int, int);
 
    int numroc_(int*, int*, int*, int*, int*);
}
 
int main(int argc, char **argv)
{
    int mpirank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    bool mpiroot = (mpirank == 0);
 
    /* Helping vars */
    int iZERO = 0;
 
    if (argc < 6) {
        if (mpiroot)
            cerr << "Usage: matrixTest matrixfile N M Nb Mb" << endl;
        // N = Rows , M = Cols , Nb = Row Blocks , Mb = Col Blocks
        MPI_Finalize();
        return 1;
    }
 
    int N, M, Nb, Mb;
    double *A_glob = NULL;
    double *A_glob2 = NULL;
    double *A_loc = NULL;
 
    /* Parse command line arguments */
    if (mpiroot) {
        /* Read command line arguments */
        stringstream stream;
        stream << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5];
        stream >> N >> M >> Nb >> Mb;
 
        /* Reserve space and read matrix (with transposition!) */
        A_glob  = new double[N*M];
        A_glob2 = new double[N*M];
        string fname(argv[1]);
        ifstream file(fname.c_str());
        for (int r = 0; r < N; ++r) {
            for (int c = 0; c < M; ++c) {
                file >> *(A_glob + N*c + r);
            }
        }
 
        /* Print matrix */
        cout << "Matrix A:\n";
        for (int r = 0; r < N; ++r) {
            for (int c = 0; c < M; ++c) {
                cout << setw(3) << *(A_glob + N*c + r) << " ";
            }
            cout << "\n";
        }
        cout << endl;
    }
 
    /* Begin Cblas context */
    /* We assume that we have 4 processes and place them in a 2-by-2 grid */
    int ctxt, myid, myrow, mycol, numproc;
    int procrows = 2;
    int proccols = 2;
    Cblacs_pinfo(&myid, &numproc);
    Cblacs_get(0, 0, &ctxt);
    Cblacs_gridinit(&ctxt, "Row-major", procrows, proccols);
    Cblacs_pcoord(ctxt, myid, &myrow, &mycol);
 
    /* Print grid pattern */
    if (myid == 0)
        cout << "Processes grid pattern:" << endl;

    for (int r = 0; r < procrows; ++r) 
    {
        for (int c = 0; c < proccols; ++c) 
        {
            Cblacs_barrier(ctxt, "All");
            if (myrow == r && mycol == c) {
                cout << myid << " " << flush;
            }
        }
        Cblacs_barrier(ctxt, "All");
        if (myid == 0)
            cout << endl;
    }
 
    /*****************************************
     * HERE BEGINS THE MOST INTERESTING PART *
     *****************************************/
 
    /* Broadcast of the matrix dimensions */
    int dimensions[4];
    if (mpiroot) 
    {
        dimensions[0] = N;
        dimensions[1] = M;
        dimensions[2] = Nb;
        dimensions[3] = Mb;
    }

    MPI_Bcast(dimensions, 4, MPI_INT, 0, MPI_COMM_WORLD);
    N = dimensions[0];
    M = dimensions[1];
    Nb = dimensions[2];
    Mb = dimensions[3];
 
    /* Reserve space for local matrices */
    // Number of rows and cols owned by the current process
    int nrows = numroc_(&N, &Nb, &myrow, &iZERO, &procrows);
    int ncols = numroc_(&M, &Mb, &mycol, &iZERO, &proccols);
    printf("rank: %d, [%d, %d]: rows: %d, cols: %d \n", myid, myrow, mycol, nrows, ncols);
    for (int id = 0; id < numproc; ++id) {
        Cblacs_barrier(ctxt, "All");
    }

    A_loc = new double[nrows*ncols];
    for (int i = 0; i < nrows*ncols; ++i) 
    {
        *(A_loc+i)=0.;
    }
 
    /* Scatter matrix */
    int sendr = 0;
    int sendc = 0;
    int recvr = 0;
    int recvc = 0;
    for (int r = 0; r < N; r += Nb, sendr=(sendr+1)%procrows)
    {
        sendc = 0;
        // Number of rows to be sent
        // Is this the last row block?
        int nr = Nb;
        if (N-r < Nb)
            nr = N-r;
 
        for (int c = 0; c < M; c += Mb, sendc=(sendc+1)%proccols)
        {
            // Number of cols to be sent
            // Is this the last col block?
            int nc = Mb;
            if (M-c < Mb)
                nc = M-c;
 
            if (mpiroot) {
                // Send a nr-by-nc submatrix to process (sendr, sendc)
                Cdgesd2d(ctxt, nr, nc, A_glob+N*c+r, N, sendr, sendc);
            }
 
            if (myrow == sendr && mycol == sendc) {
                // Receive the same data
                // The leading dimension of the local matrix is nrows!
                Cdgerv2d(ctxt, nr, nc, A_loc+nrows*recvc+recvr, nrows, 0, 0);
                recvc = (recvc+nc)%ncols;
            }
        }
 
        if (myrow == sendr)
            recvr = (recvr+nr)%nrows;
    }
 
    /* Print local matrices */
    for (int id = 0; id < numproc; ++id)
    {
        if (id == myid) {
            cout << "A_loc on node " << myid << endl;
            for (int r = 0; r < nrows; ++r) {
                for (int c = 0; c < ncols; ++c)
                    cout << setw(3) << *(A_loc+nrows*c+r) << " ";
                cout << endl;
            }
            cout << endl;
        }
        Cblacs_barrier(ctxt, "All");
    }
 
    /* Gather matrix */
    sendr = 0;
    for (int r = 0; r < N; r += Nb, sendr=(sendr+1)%procrows)
    {
        sendc = 0;
        // Number of rows to be sent
        // Is this the last row block?
        int nr = Nb;
        if (N-r < Nb)
            nr = N-r;
 
        for (int c = 0; c < M; c += Mb, sendc=(sendc+1)%proccols)
        {
            // Number of cols to be sent
            // Is this the last col block?
            int nc = Mb;
            if (M-c < Mb)
                nc = M-c;
 
            if (myrow == sendr && mycol == sendc) {
                // Send a nr-by-nc submatrix to process (sendr, sendc)
                Cdgesd2d(ctxt, nr, nc, A_loc+nrows*recvc+recvr, nrows, 0, 0);
                recvc = (recvc+nc)%ncols;
            }
 
            if (mpiroot) {
                // Receive the same data
                // The leading dimension of the local matrix is nrows!
                Cdgerv2d(ctxt, nr, nc, A_glob2+N*c+r, N, sendr, sendc);
            }
 
        }
 
        if (myrow == sendr)
            recvr = (recvr+nr)%nrows;
    }
 
    /* Print test matrix */
    if (mpiroot)
    {
        cout << "Matrix A test:\n";
        for (int r = 0; r < N; ++r) {
            for (int c = 0; c < M; ++c) {
                cout << setw(3) << *(A_glob2+N*c+r) << " ";
            }
            cout << endl;
        }
    }
 
    /************************************
     * END OF THE MOST INTERESTING PART *
     ************************************/
 
    /* Release resources */
    delete[] A_glob;
    delete[] A_glob2;
    delete[] A_loc;
    Cblacs_gridexit(ctxt);
   
    MPI_Finalize();
}
