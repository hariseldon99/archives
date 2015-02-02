/*
* C++ Interface: algebra
*
* Description:
*
*
* Author: Analabha Roy <daneel@physics.utexas.edu>, (C) 2008
*
* Copyright: See COPYING file that comes with this distribution
*
*/
#define TOLERANCE 0.4
/*Tolerance for non-orthogonality*/

int select_vec (double matrix[], int i, double vector[], int rowsize);
int write_vec (double vector[], int i, double matrix[], int rowsize);
int arrcopy (double source[], double target[], int arrsize);
double direct_prod (double a_real[], double a_imag[], double b_real[],
		    double b_imag[], int size);
