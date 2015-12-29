#ifndef COMMON_H
#define COMMON_H

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

#endif
