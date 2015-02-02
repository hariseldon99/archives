/*Buffer size for all the data, should be reasonably large*/
#define BUFSIZE 1600000

/*Hamiltonian Parameters*/
#define V0 4.91345043
#define U0 -1.0

/*Classical Poincare Section Parameters*/
/*Standard Deviation of the Interaction Gaussian*/
#define SIGMA 0.005
#define ENERGY -4.5
/*Numerical error parameters for ODE integration */
#define NADT 1E-4
#define STROBERR 1E-4
/*System Parameters*/
#define TIME 5000.0

/*Quantum Husimi Plot Parameters*/
/*Length of the box for particle-in-a-box states*/
#define L 3.5
/*Strobe Point for x2*/
#define STROBE 1.0
/***********************************Numerical Truncation Parameters********************************/
/*ALPHA+ ALPHA*(ALPHA-1)/2 is DIM. ALPHA is the # of 2-free-particles-in-a-box-states 16 works*/
#define ALPHA 50
/*Truncated Dimension of the Unperturbed Hamiltonian generated from 2-free-particles-in-a-box states 136 works*/
#define STATENO 1275
/***********************************************ERROR TOLERANCES***********************************/
#define MACHINENUM 1E-5
/*Ranges for the Husimi plot*/
#define XINIT -3.5
#define XFINAL 3.5
#define PINIT -7.0
#define PFINAL 7.0

/********The Parameter epsilon & the tower of 2-particle states are stored in this datatype********/
typedef struct
{
  double v0;
  double u0;
  double s1, s2;
  double eigenvalues[STATENO], eigenvectors[STATENO][STATENO];
  int tower[STATENO][2];
} drive_and_tower;
