#include <mpi.h>
#define numroc_ NUMROC
#define descinit_ DESCINIT
//#include <iostream>
#include <math.h>
#include <mkl_pblas.h>
#include <mkl_scalapack.h>
#include <mkl_blacs.h>
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
    void Cblacs_gridmap(int*, int*, int, int, int);
    void Cblacs_exit(int);    
    int numroc_(int*, int*, int*, int*, int*);
}


int main(int argc, char **argv)
{


    int     context, desca[9], descb[9], ib, info, ja, laf, lda, ldb,
        lwork, mb, mype, n, nb, npcol, npe, nprow, nrhs;
    double  b[4], d[4], dl[4], du[4];
    double  *af, *work;
    char    trans = 'N';

    /* Set array dimensions and blocking */
    n = 8;            /* dimension of the problem */
    lda = 4;          /* leading dimension of A */
    ldb = 4;          /* leading dimension of B */
    nrhs = 1;         /* number of right-hand sides  */
    npcol = 2;        /* number of processor columns */
    ja = 1;           /* offset for A */
    ib = ja;          /* offset for B */
    mb = 4;           /* blocking */
    nb = 4;           /* blocking */

    laf = 12 * npcol + 3 * nb;
    af = (double *)malloc(laf*sizeof(double));
    lwork = 10 * npcol + 4 * nrhs;
    work = (double *)malloc(lwork*sizeof(double));

    /* Start BLACS */
    Cblacs_pinfo(&mype, &npe);
    Cblacs_get(0, 0, &context);
    Cblacs_gridinit(&context, "R", 1, npe);

    if (mype == 0){
        /* PE = 0 gets D(1:4), DL(1:4), DU(1:4) and B(1:4) */
        d[0] = 1.8180; d[1] = 1.6602; d[2] = 1.3420; d[3] = 1.2897;
        dl[0] = 0.0000; dl[1] = 0.8385; dl[2] = 0.5681; dl[3] = 0.3704;
        du[0] = 0.6946; du[1] = 0.4449; du[2] = 0.5466; du[3] = 0.7027;
        b[0] = 1.0; b[1] = 2.0; b[2] = 3.0; b[3] = 4.0;
    }
    else if (mype == 1){
        /* PE = 1 gets D(5:8), DL(5:8), DU(5:8) and B(5:8) */
        d[0] = 1.3412; d[1] = 1.5341; d[2] = 1.7271; d[3] = 1.3093;
        dl[0] = 0.7027; dl[1] = 0.5466; dl[2] = 0.4449; dl[3] = 0.6946;
        du[0] = 0.3704; du[1] = 0.5681; du[2] = 0.8385; du[3] = 0.0000;
        b[0] = 5.0; b[1] = 6.0; b[2] = 7.0; b[3] = 8.0;
    }

    /* Array descriptor for A (D, DL and DU) */
    desca[0] = 501; desca[1] = context; desca[2] = n; desca[3] = nb;
    desca[4] = 0; desca[5] = lda; desca[6] = 0;

    /* Array descriptor for B */
    descb[0] = 502; descb[1] = context; descb[2] = n; descb[3] = nb;
    descb[4] = 0; descb[5] = ldb; descb[6] = 0;

    /* Factorization */
    pddttrf(&n, dl, d, du, &ja, desca, af, &laf, work, &lwork, &info);

    /* Solution */
    //pddttrs(&trans, &n, &nrhs, dl, d, du, &ja, desca, b, &ib, descb,
    //  af, &laf, work, &lwork, &info);

    printf("MYPE=%i: x[:] = %7.4f %7.4f %7.4f %7.4f\n",
        mype, b[0], b[1], b[2], b[3]);

    Cblacs_gridexit(context);
    Cblacs_exit(0);
    
    return 0;
}
