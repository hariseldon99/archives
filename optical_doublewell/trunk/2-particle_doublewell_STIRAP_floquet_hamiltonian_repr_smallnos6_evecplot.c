#include <stdio.h>
#include <math.h>
#include <time.h>

#include <omp.h>
#define CHUNKSIZE 1

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>

#include <essl.h>
#define CSET(x,a,b) (RE(x)=a, IM(x)=b)

/*************************************System parameters**********************************************/
/*Length of the box for particle-in-a-box states*/
#define L 3.5

/* Pseudopotential interaction strength, adjusted to minimize the detuning that is required to make drive freqs commensurate*/
#define U0 -1.0

/*Double Well Depth*/
#define V0 7.2912229

/*STIRAP Pulse Amplitude*/
#define EPSILON0 115.0

/****************************************************************************************************/
 /*STIRAP pulse width T dflt=4500 */
#define STIRAPWIDTH 4500.0
/*centroid of the stoke pulse */
#define TSTOKE 12000.00
/*centroid of the pump pulse */
#define TPUMP 24000.0

/*Adiabatic Time(s)*/
#define TFIXMIN 0.0
#define TFIXMAX 36000.0
/*Exact total of  10  oscillations b4 next increment (same as 1 osc) ideally */
#define TFIXINC 68.50174694488798
/******Duffing Drive Parameters*********************/
#define OMEGASTOKE 6.3288865835501
#define OMEGAPUMP 0.7337839483516753
/*******Dipole moments****************************/
#define D12 0.10261548869550063
#define D24 0.0027445724960265872

/*omegastoke/omegapump is n1/n2 frequencies MUST BE COMMENSURATE CHECK THIS THOROUGHLY*/
#define N1 69.0
#define N2 8.0
/*******************************/

/***********************************Numerical Truncation Parameters********************************/
/*ALPHA+ ALPHA*(ALPHA-1)/2 is DIM. ALPHA is the # of 2-free-particles-in-a-box-states 16 works*/
#define ALPHA 22
/*Truncated Dimension of the Unperturbed Hamiltonian generated from 2-free-particles-in-a-box states 136 works*/
#define STATENO 253
/*Manual tagging of floquet state at Tfix=0*/
#define TAGVAL -0.000865
#define DOTPRODUCT 0.9
/*Truncation of Floquet Matrix*/
#define DIM 25
/*Output File(s)*/
char OUT_FILE_QUASI[] = "result/floquet_smallnos6_quasi_e_115_phiC.dat";
char OUT_FILE_PROBS[] = "result/floquet_smallnos6_probs_e_115_phiC.dat";
/***********************************************ERROR TOLERANCES***********************************/
/*
The step-size adjustment procedure for this method begins by computing the desired error level D_i for each component,

D_i = ABSERROR + RELERROR * (YERROR |y_i| + YPRIMEERROR*DT |y'_i|)

and comparing it with the observed error E_i = |yerr_i|. If the observed error E exceeds the desired error level D by more than 10% for any component then the 
method reduces the step-size by an appropriate factor
*/
#define DT 1e-3
#define ABSERROR 1e-7
#define RELERROR 0.0
#define YERROR 0.5
#define YPRIMEERROR 0.2
#define MACHINENUM 1E-5

/********The Parameter epsilon & the tower of 2-particle states are stored in this datatype********/

typedef struct
{
  double v0;
  double u0;
  double omegastoke;
  double omegapump;
  double epsilonstoke;
  double epsilonpump;
  double eigenvalues[STATENO];
  double dipolematrix[STATENO][STATENO];
  int tower[STATENO][2];
} drive_and_tower;

/*****************************Defines Kronecker Delta**********************************************/
double
kdel (int i, int j)
{
  int result;
  if (i == j)
    {
      result = 1.0;
    }
  else
    result = 0.0;
  return (result);
}

/*****************Calculates the Adiabatic Amplitude of the STIRAP Pulse**************************/

double
amp (double tfix, double tmean)
{
  double ampl;
  double beta = 1.0 / (2.0 * STIRAPWIDTH * STIRAPWIDTH);
  ampl = EPSILON0 * exp (-beta * (tfix - tmean) * (tfix - tmean));
  return (ampl);
}

/*****************This is the lth single particle-in-a-box energy level*****************************/
double
En (int l)
{
  double energy;
  l = l + 1;
  energy = M_PI / 2.0;
  energy = energy * energy;
  energy = energy * l * l;
  energy = energy / (L * L);
  return (energy);
}

/**************************** This is a matrix element of x*****************************************/
double
f1 (int m, int n)
{
  double value;
  m = m + 1;			/*Arrays are counted from n=0, but levels from n=1 */
  n = n + 1;
  if ((m + n) % 2 == 1)
    {
      value = 16.0 * m * n;
      value = value / (M_PI * M_PI);
      value = value / (m * m - n * n);
      value = value / (m * m - n * n);
      value = value * L;
    }
  else
    {
      value = 0.0;
    }
  return (value);
}

/**************************** This is a matrix element of x^2**************************************/
double
f2 (int m, int n)
{
  double value, temp;
  m = m + 1;			/*Arrays are counted from n=0, but levels from n=1 */
  n = n + 1;
  if ((m + n) % 2 == 0)
    {
      if (m == n)
	{
	  value = ((1 / 3.0) - (2.0 / (n * n * M_PI * M_PI))) * L * L;
	}
      else
	{
	  value =
	    32.0 * m * n * L * L / (M_PI * M_PI * (m * m - n * n) *
				    (m * m - n * n));
	}
    }
  else
    {
      value = 0.0;
    }
  return (value);
}

/**************************** This is a matrix element of x^4***************************************/
double
f4 (int m, int n)
{
  double value;
  m = m + 1;			/*Arrays are counted from n=0, but levels from n=1 */
  n = n + 1;
  if ((m + n) % 2 == 0)
    {
      if (m == n)
	{
	  value =
	    L * L * L * L * ((1 / 5.0) - (4.0 / (n * n * M_PI * M_PI)) +
			     (24.0 /
			      (n * n * n * n * M_PI * M_PI * M_PI * M_PI)));
	}
      else
	{
	  value =
	    (64.0 * m * n / (M_PI * M_PI * M_PI * M_PI)) *
	    (pow (M_PI / (m * m - n * n), 2) -
	     (48.0 * (m * m + n * n) / pow (m * m - n * n, 4))) * L * L * L *
	    L;
	}
    }
  else
    {
      value = 0.0;
    }
  return (value);
}

/*Matrix elements of the Double Well Potential*/
double
pot (int n, int m)
{
  double f2 (int, int), f4 (int, int);
  double potelm;
  potelm = (-2.0 * f2 (n, m) + f4 (n, m));
  return (potelm);
}

/*Matrix Element of the Interaction*/
double
in (int n1, int n2, int n3, int n4)
{
  double term1 = 1.0 / (4.0 * L), term2 = 1.0 / (2.0 * L), term3 =
    3.0 / (4.0 * L);
  double interelement = 0.0;	/*Default Value of Interaction */
  n1 = n1 + 1;
  n2 = n2 + 1;
  n3 = n3 + 1;
  n4 = n4 + 1;

  if ((n1 - n2) == (n3 + n4))
    {
      if ((n1 - n2) == 0)
	{
	  interelement = term2;
	}
      else
	{
	  interelement = -term1;
	}
    }

  if ((n1 + n2) == (n3 + n4))
    {
      if ((n1 + n2) == 0)
	{
	  interelement = term2;
	}
      else
	{
	  interelement = term1;
	}
    }

  if ((n1 - n2) == (n4 - n3))
    {
      if ((n1 - n2) == 0)
	{
	  interelement = term2;
	}
      else
	{
	  interelement = term1;
	}
    }

  if ((n1 + n2) == (n4 - n3))
    {
      if ((n1 + n2) == 0)
	{
	  interelement = term2;
	}
      else
	{
	  interelement = -term1;
	}
    }

  if ((n1 - n2) == (n3 - n4))
    {
      if ((n1 - n2) == 0)
	{
	  interelement = term2;
	}
      else
	{
	  interelement = term1;
	}
    }

  if ((n1 + n2) == (n3 - n4))
    {
      if ((n1 + n2) == 0)
	{
	  interelement = term2;
	}
      else
	{
	  interelement = -term1;
	}
    }

  if ((n1 - n2) == (-n3 - n4))
    {
      if ((n1 - n2) == 0)
	{
	  interelement = term2;
	}
      else
	{
	  interelement = -term1;
	}
    }

  if ((n1 + n2) == (-n3 - n4))
    {
      if ((n1 + n2) == 0)
	{
	  interelement = term2;
	}
      else
	{
	  interelement = term1;
	}
    }

  if ((n1 == n2) && (n2 == n3) && (n3 == n4))
    {
      interelement = term3;
    }
  /*Bosonization */
  interelement = 2.0 * interelement;
  return (interelement);
}

/*This uses the above subroutines to evaluate the hamiltonian matrix elements*/
double
hamilt (void *param, int nt, int mt, double t)
{
  double kdel (int, int), En (int), in (int, int, int, int), pot (int, int);
  double matrixelm;
  int i, j;
  drive_and_tower *p = (drive_and_tower *) param;
  /*Get epsilon & towerofstates from structure param */
  double v0 = p->v0;
  double u0 = p->u0;
  double epsilonstoke = p->epsilonstoke;
  double epsilonpump = p->epsilonpump;
  double omegastoke = p->omegastoke;
  double omegapump = p->omegapump;
  double eigenvalue, dipolematrix;

  eigenvalue = p->eigenvalues[nt];
  dipolematrix = p->dipolematrix[mt][nt];

  /*Write out the matrix element in the hamiltonian representation */
  matrixelm = dipolematrix;
  matrixelm =
    matrixelm * (epsilonstoke * sin (OMEGASTOKE * t) +
		 epsilonpump * sin (OMEGAPUMP * t));
  matrixelm = matrixelm + eigenvalue * kdel (mt, nt);
  return (matrixelm);
}

	     /*This uses the above subroutines to evaluate the undriven hamiltonian matrix elements */
double
undriven_hamilt (void *param, int nt, int mt)
{
  double kdel (int, int), En (int), in (int, int, int, int), pot (int, int);
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


/*********************This evaluates the Differential Equation System*****************************/
int
func (double t, const double y[], double dydt[], void *param)
{
  double y_re[DIM], y_im[DIM], dydt_re[DIM], dydt_im[DIM];
  int nt, mt, chunk = CHUNKSIZE;
  /*Write out the (complicated) expression for RHS & assign to dydt[i] */
  for (nt = 0; nt < DIM; nt++)
    {
      y_re[nt] = y[nt];
      y_im[nt] = y[nt + DIM];
    }

#pragma omp parallel shared(dydt_re,dydt_im,y_re,y_im) private(nt,mt)
  {
#pragma omp for schedule(dynamic,chunk)
    for (nt = 0; nt < DIM; nt++)
      {
	dydt_re[nt] = 0.0;
	dydt_im[nt] = 0.0;
	for (mt = 0; mt < DIM; mt++)
	  {
	    dydt_re[nt] = dydt_re[nt] + hamilt (param, nt, mt, t) * y_im[mt];
	    dydt_im[nt] = dydt_im[nt] - hamilt (param, nt, mt, t) * y_re[mt];
	  }
      }
  }

  for (nt = 0; nt < DIM; nt++)
    {
      dydt[nt] = dydt_re[nt];
      dydt[DIM + nt] = dydt_im[nt];
    }
  return GSL_SUCCESS;
}

/*****This is a DUMMY FUNCTION. IT DOES NOTHING. IT'S ONLY NEEDED FOR BSIMP & I'M USING RK8PD******/
int
jac (double t, const double y[], double *dfdy, double dfdt[], void *param)
{
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2 * DIM, 2 * DIM);
  gsl_matrix *jacobian = &dfdy_mat.matrix;
  printf ("\nfukd");
  return GSL_SUCCESS;
}

/*****This function actually runs the full integration of a particulat set of IC's from 0-T**********/
void
integrate (double *input, double initial, double final, void *param)
{
  int i, j, status;
  double dt, t;
  /*This sets up the gsl ODE system structure */
  const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;
  gsl_odeiv_step *s = gsl_odeiv_step_alloc (T, 2 * DIM);
  gsl_odeiv_control *c =
    gsl_odeiv_control_standard_new (ABSERROR, RELERROR, YERROR, YPRIMEERROR);
  gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc (2 * DIM);
  gsl_odeiv_system sys = { func, jac, 2 * DIM, param };
  dt = DT;
  t = initial;
  while (t < final)
    {
      status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, final, &dt, input);
      if (status != GSL_SUCCESS)
	{
	  printf ("\nGSL execution of differential equation solver failed.");
	  printf ("\nCheck program for bugs, bailing... \n");
	  break;
	}
    }
  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
}

/*************************************************************************************************/
int
main (void)
{
  time_t begin, end, rawtime;
  struct tm *timeinfo;
  drive_and_tower param;
  int i, j, k, l, n, m, n1, n2, m1, m2, select[DIM], naux =
    3 * (2 * DIM + 1), towerofstates[STATENO][2], buffer[2];
  int true, chunk = CHUNKSIZE, nthreads, tid, count;
  long errorcount = 0;
  double comtest1, comtest2, dipolematrix[STATENO][STATENO],
    eigenvectors[STATENO][STATENO];
  double y[2 * DIM], y_out[2 * DIM], aux[3 * (2 * DIM + 1)];
  double quasi[DIM], shuldbeone, temp;
  double tfix, period, element, sum, term;
  double v0, u0, epsilonstoke, epsilonpump, omegastoke, omegapump;
  double pumptime = TPUMP, stoketime = TSTOKE;
  double eni, enj;
  double ip_re, ip_im, fl_re, fl_im, f_re, f_im, tagval = TAGVAL, compare;
  int tagno, evolve = 0;
  dcmplx floquet4essl[DIM][DIM], floqueteigenvalues[DIM],
    floqueteigenvectors[DIM][DIM], floquetstate[DIM];

  FILE *floqueteigenval, *floquetprobs;

  gsl_matrix *floquetmatrix = gsl_matrix_alloc (2 * DIM, DIM);
  gsl_matrix *hamiltonian = gsl_matrix_alloc (STATENO, STATENO);
  gsl_matrix *eigenbasis = gsl_matrix_alloc (STATENO, STATENO);
  gsl_vector *eigenvalues = gsl_vector_alloc (STATENO);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (STATENO);

  floqueteigenval = fopen (OUT_FILE_QUASI, "w");
  floquetprobs = fopen (OUT_FILE_PROBS, "w");
  begin = time (NULL);
  printf
    ("\n=========================================================================");

  /*Compute the adiabatic period of the STIRAP pulse */
  comtest1 = N1 / N2;
  comtest2 = (OMEGASTOKE / OMEGAPUMP);
  if ((comtest1 - comtest2) >= MACHINENUM)
    {
      printf ("\n ERROR!!!\n");
      printf
	("\nFrequencies are not commensurate, or N1, N2 entries are wrong, bailing...");
      printf
	("\n=========================================================================");
      exit (1);
    }

  /*Build the 2-particle tower of states */
  count = 0;
  for (i = 0; i < ALPHA; i++)
    {
      for (j = 0; j <= i; j++)
	{
	  towerofstates[count][0] = i;
	  towerofstates[count][1] = j;
	  count++;
	}
    }
  /*Now sort them in order of increasing free 2-particle energy n1^2+n2^2
     for(i=0;i<STATENO;i++)
     {
     for(j=0;j<=i;j++)
     {
     eni=(towerofstates[i][0]+1)*(towerofstates[i][0]+1)+(towerofstates[i][1]+1)*(towerofstates[i][1]+1);
     enj=(towerofstates[j][0]+1)*(towerofstates[j][0]+1)+(towerofstates[j][1]+1)*(towerofstates[j][1]+1);
     if(eni<enj)
     {
     buffer[0]=towerofstates[i][0];
     buffer[1]=towerofstates[i][1];

     towerofstates[i][0]=towerofstates[j][0];
     towerofstates[i][1]=towerofstates[j][1];

     towerofstates[j][0]=buffer[0];
     towerofstates[j][1]=buffer[1];
     } 
     }
     } REDUNDANT??? */
  u0 = U0;
  v0 = V0;
  omegastoke = OMEGASTOKE;
  omegapump = OMEGAPUMP;
  period = M_PI * ((N1 / omegastoke) + (N2 / omegapump));
  printf ("\nCOMMENSURATE TIME PERIOD OF BOTH PULSES = %lf", period);
  param.v0 = v0;
  param.u0 = u0;
  param.omegastoke = omegastoke;
  param.omegapump = omegapump;
  for (i = 0; i < STATENO; i++)
    {
      param.tower[i][0] = towerofstates[i][0];
      param.tower[i][1] = towerofstates[i][1];
    }
  /*Formulate and Diagonalize the Unperturbed Hamiltonian here */
  for (i = 0; i < STATENO; i++)
    {
      for (j = 0; j < STATENO; j++)
	{
	  element = undriven_hamilt (&param, i, j);
	  gsl_matrix_set (hamiltonian, i, j, element);
	}
    }
  printf ("\n Diagonalizing hamiltonian, please wait...");
  /*Diagonalize the hamiltonian */
  gsl_eigen_symmv (hamiltonian, eigenvalues, eigenbasis, w);
  gsl_eigen_symmv_free (w);

  /*Sort eigenvectors in order of increasing eigenvalues */
  gsl_eigen_symmv_sort (eigenvalues, eigenbasis, GSL_EIGEN_SORT_VAL_ASC);

  printf ("\n Undriven hamiltonian Diagonalised, Eigenvalues are Below:");
  printf ("\n");
  for (i = 0; i < STATENO; i++)
    printf ("%lf ", gsl_vector_get (eigenvalues, i));

  /*Put eigenbasis in an array & in struct param */
  for (i = 0; i < STATENO; i++)
    {
      param.eigenvalues[i] = gsl_vector_get (eigenvalues, i);
      for (j = 0; j < STATENO; j++)
	{
	  eigenvectors[i][j] = gsl_matrix_get (eigenbasis, i, j);
	}
    }
  /*Calculate the Dipole Matrix elements in the Hamiltonian Representation */
  /*And put them in struct param */
  printf
    ("\nCalculating Dipole Matrix elements in Hamiltonian Representation, Please wait...");
  for (n = 0; n < STATENO; n++)
    {
      for (m = 0; m < STATENO; m++)
	{
	  n1 = towerofstates[n][0];
	  n2 = towerofstates[n][1];
	  m1 = towerofstates[m][0];
	  m2 = towerofstates[m][1];

	  sum = 0.0;
	  for (k = 0; k < STATENO; k++)
	    {
	      for (l = 0; l < STATENO; l++)
		{
		  term = f1 (n1, m1) * kdel (n2, m2);
		  term = term + f1 (n1, m2) * kdel (n2, m1);
		  term = term + f1 (n2, m2) * kdel (n1, m1);
		  term = term + f1 (n2, m1) * kdel (n1, m2);
		  term = term * eigenvectors[k][n] * eigenvectors[l][m];
		  sum = sum + term;
		}
	    }
	  /*Normalization */
	  sum = sum / sqrt (1 + kdel (n1, n2));
	  sum = sum / sqrt (1 + kdel (m1, m2));
	  dipolematrix[n][m] = sum;
	  param.dipolematrix[n][m] = dipolematrix[n][m];
	}
    }

  /*tfix  loop begins here */
  for (tfix = TFIXMIN; tfix <= TFIXMAX; tfix = tfix + TFIXINC)
    {
      /*Now, put towerofstates, epsilon, vo & u0 in structure param */

      epsilonstoke = D12 * amp (tfix, stoketime);
      param.epsilonstoke = epsilonstoke;
      epsilonpump = D24 * amp (tfix, pumptime);
      param.epsilonpump = epsilonpump;

      printf ("\nELAPSED TIME = %d sec", time (NULL) - begin);
      printf ("\n\nEVOLVING FLOQUET MATRIX WITH TFIX =  %lf\n", tfix);

	   /****************HERE IS WHERE THE PARALLEL FUN BEGINS********************************/
#pragma omp parallel shared(period) private(i,j,y)
      {
	/*Build the Floquet Matrices (Real & Imaginary separate) using RK method in gsl */
#pragma omp for schedule(dynamic,chunk)
	for (j = 0; j < DIM; j++)
	  {
	    tid = omp_get_thread_num ();
	    /*This is the jth column of the matrix */
	    /* Initial conditions, */
	    for (i = 0; i < 2 * DIM; i++)
	      {
		if (i < DIM)
		  {
		    y[i] = kdel (i, j);
		  }
		else
		  {
		    y[i] = 0.0;
		  }
	      }
	    /*DONE. Now evolve the jth column of the floquet matrix (stored in array y[j]) across time period */
	    integrate (y, tfix, tfix + period, &param);
	    /*DONE, Now print out the final values of floquet matrix for that array into the gsl matrix */
	    for (i = 0; i < 2 * DIM; i++)
	      {
		gsl_matrix_set (floquetmatrix, i, j, y[i]);
	      }
	     /*DONE*/ if (tid == 0)
	      {
		nthreads = omp_get_num_threads ();
	      }
	  }
	/*End */
      }
	   /****************HERE IS WHERE THE PARALLEL FUN ENDS********************************/

      /*Formulate the floquet matrix from gsl form to essl form */
      for (i = 0; i < DIM; i++)
	{
	  for (j = 0; j < DIM; j++)
	    {
	      CSET (floquet4essl[i][j], gsl_matrix_get (floquetmatrix, i, j),
		    gsl_matrix_get (floquetmatrix, i + DIM, j));
	    }
	}
      /* End the formulation of the essl floquet matrix in floquet4essl array */

      /*We want to diagonalize the floquet matrix now with eigenvectors also */
      zgeev (1, floquet4essl, DIM, floqueteigenvalues, floqueteigenvectors,
	     DIM, select, DIM, aux, naux);
      /*End Diagonalization */

      /*Dump out the floquet quasienergies */
      printf
	("\n=========================================================================");
      printf ("\nChecking for Unitarity of Floquet Matrix...\n");
      for (i = 0; i < DIM; i++)
	{
	  quasi[i] =
	    -atan2 (IM (floqueteigenvalues[i]),
		    RE (floqueteigenvalues[i])) / period;
	  shuldbeone =
	    RE (floqueteigenvalues[i]) * RE (floqueteigenvalues[i]) +
	    IM (floqueteigenvalues[i]) * IM (floqueteigenvalues[i]);
	  if (fabs (shuldbeone - 1.0) >= MACHINENUM)
	    {
	      errorcount = errorcount + 1;
	      printf
		("Eigenvalue %d =  %lf +I %lf , quasienergy = %lf has non-unit modulus = %lf\n",
		 i, RE (floqueteigenvalues[i]), IM (floqueteigenvalues[i]),
		 quasi[i], shuldbeone);
	    }
	}
      printf
	("\n=========================================================================");
      fprintf (floqueteigenval, "\n%lf", tfix);
      for (i = 0; i < DIM; i++)
	{
	  fprintf (floqueteigenval, " %10.6lf", quasi[i]);
	}

      /*Need to figure out which eigenvector to select */
      /*If first iteration then manually tag the floquet state */
      if (evolve == 0)
	{

	  for (i = 0; i < DIM; i++)
	    {
	      compare = fabs (quasi[i] - tagval);
	      if (compare < 1E-4)
		{
		  tagno = i;
		  break;
		}
	    }
	  for (i = 0; i < DIM; i++)
	    {
	      floquetstate[i] = floqueteigenvectors[tagno][i];
	    }
	}
      else
	{
	  for (i = 0; i < DIM; i++)
	    {
	      /*Calculate innerproduct of floquetstate & floqueteigenvectors[i][] */
	      ip_re = 0.0;
	      ip_im = 0.0;
	      for (j = 0; j < DIM; j++)
		{
		  fl_re = floquetstate[j]._data._re;
		  fl_im = floquetstate[j]._data._im;
		  f_re = floqueteigenvectors[i][j]._data._re;
		  f_im = floqueteigenvectors[i][j]._data._im;
		  ip_re = ip_re + fl_re * f_re + fl_im * f_im;
		  ip_im = ip_im + fl_re * f_im - fl_im * f_re;
		}

	      /*The one that's nonzero will be the new tagged floquet state */
	      if (ip_re * ip_re + ip_im + ip_im > DOTPRODUCT)
		{
		  for (j = 0; j < DIM; j++)
		    {
		      floquetstate[j] = floqueteigenvectors[i][j];
		    }
		  break;
		}

	    }
	}
      /*Dump out the components */
      fprintf (floquetprobs, "\n%lf", tfix);
      for (j = 0; j < DIM; j++)
	{
	  fprintf (floquetprobs, " %10.6lf",
		   floquetstate[j]._data._re * floquetstate[j]._data._re +
		   floquetstate[j]._data._im * floquetstate[j]._data._im);
	}
      evolve++;
    }
  /*tfix loop ends here */

  end = time (NULL);
  printf
    ("\n=========================================================================");
  printf ("\nTIME PERIOD OF DRIVE = %lf", period);
  printf
    ("\n=========================================================================");
  printf
    ("\n                 PROGRAM RAN SUCCESSFULLY (I THINK)                      ");
  printf
    ("\n NUMBER OF TIMES THE IMAGINARY PART OF THE QUASIENERGY GOT TOO BIG = %ld ",
     errorcount);
  printf
    ("\n                    TIME ELAPSED = %ld sec                               ",
     end - begin);
  printf
    ("\n                   NUMBER OF THREADS = %ld                               ",
     nthreads);
  printf
    ("\n=========================================================================\n");

  time (&rawtime);
  timeinfo = localtime (&rawtime);
  fclose (floqueteigenval);
  fclose (floquetprobs);
  return 0;
}
