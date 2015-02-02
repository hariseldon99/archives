/*
// C++ Interface: floquet_optical_lattice
//
// Description:
//
//
// Author: Analabha Roy <daneel@physics.utexas.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
*/
#define RATTOL 1E-6
static char help[] =
        "The slepc_eval_nhep() routine takes a 1D array of COMPLEX NUMBERS, writes it into a NONHERMITIAN PetSc Mat datatype, diagonalizes it using the parallel Slepc library and prints eigenvalues and eigenvectors. Returns # of converged eigenvalues"
        " -n <n>, where <n> = number of grid subdivisions = matrix dimension.\n\n";
double c_mat (int m, int n), s_mat (int m, int n), psq_mat (int n),
kdel (int i, int j);
double in_mat (int m1, int m2, int n1, int n2);
int integrate (double *y, double initial, double final, void *out);
int slepc_eval_nhep (int argc, char **argv, double matrixreal[],
                     double matriximag[], int psize, double evals_re[],
                     double evals_im[], double evecs_real[],
                     double evecs_im[], MPI_Comm);
long quasienergies (double evals_real[], double evals_imag[], int size,
                    double period, double quasi[]);
void final_output (double period, time_t begin, time_t end, int nprocs);
double hamilt_sp (void *param, int nt, int mt, double t);
double hamilt (void *param, int nt, int mt, double t);
int floquetsort (double floquetevals_real[], double floquetevals_imag[],
                 double floquetevecs_real[], double floquetevecs_imag[],
                 double floquetprev_real[], double floquetprev_imag[],
                 double floquetvecs_prev_real[],
                 double floquetvecs_prev_imag[], int size);
int arrcopy (double source[], double target[], int size);
void print_arr (double matrix_real[], double matrix_imag[], int size);
int check_unitarity (double matrix_real[], double matrix_imag[], int size);
int fold_quasi (double array[], double omega, int size);
int tag_quasi (int tagarray[], double quasi[], double unperturbed_quasi[],
               int size);

int select_vec (double matrix[], int i, double vector[], int rowsize);
double direct_prod (double a_real[], double a_imag[], double b_real[],
                    double b_imag[], int size);