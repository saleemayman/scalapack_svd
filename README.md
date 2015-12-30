# scalapack_svd
Implementing SVD using ScaLAPACK in C/C++

This code will read matrix data from a CSV file (data stored in matrix coordinate form  with  3 columns: row id, column id and the matrix entry for that row and column id). All processors in the grid will read the data file and construct their own 2D block-cyclic local matrices (local here means the processor's chunk of the big global matrix). Once the data is cyclically distributed among the processor grid the SVD is computed.

The important part is that the local matrices are stored in column-major format and of course in 2D block cyclic manner.
