#include <stdio.h>
#include <stdlib.h>
#include<math.h>

#include "params.h"

double kdel (int, int), En (int), in (int, int, int, int), pot (int, int);

  /*This uses the above subroutines to evaluate the undriven hamiltonian matrix elements */
double
			undriven_hamilt (void *param, int nt, int mt)
{
	double matrixelm;
	int n1, n2, m1, m2, i, j;
	drive_and_tower *p = (drive_and_tower *) param;
	/*Get epsilon & towerofstates from structure param */
	double v0 = p->v0;
	double u0 = p->u0;
	double towerofstates[STATENO][2];

	for (i = 0; i < STATENO; i++)
	{
		towerofstates[i][0] = p->tower[i][0];
		towerofstates[i][1] = p->tower[i][1];
	}

	n1 = towerofstates[nt][0];
	n2 = towerofstates[nt][1];
	m1 = towerofstates[mt][0];
	m2 = towerofstates[mt][1];

	matrixelm = En (n1) + En (n2);
	matrixelm =
			matrixelm * (kdel (m1, n1) * kdel (m2, n2) +
			kdel (m1, n2) * kdel (m2, n1));

	matrixelm = matrixelm + v0 * pot (n1, m1) * kdel (n2, m2);
	matrixelm = matrixelm + v0 * pot (n1, m2) * kdel (n2, m1);
	matrixelm = matrixelm + v0 * pot (n2, m1) * kdel (n1, m2);
	matrixelm = matrixelm + v0 * pot (n2, m2) * kdel (n1, m1);

	matrixelm = matrixelm + u0 * in (n1, n2, m1, m2);

	/*FOR NORMALIZING THE DAMN BOSONIC FUNCTION.WTF!!! */
	matrixelm = matrixelm / sqrt (1 + kdel (n1, n2));
	matrixelm = matrixelm / sqrt (1 + kdel (m1, m2));
	return (matrixelm);
}

