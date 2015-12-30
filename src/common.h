#ifndef COMMON_H
#define COMMON_H

#define DTYPE_ 0
#define CTXT_  1
#define M_     2   
#define N_     3
#define MB_    4
#define NB_    5
#define RSRC_  6
#define CSRC_  7
#define LLD_   8

struct gridInfo {
    int context;
    int myRank;
    int numProcs;
    int myRow;
    int myCol;
    int myNumRows;
    int myNumCols;
};


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

#endif
