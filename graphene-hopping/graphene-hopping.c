#include "graphene-hopping.h"
#include <gsl/gsl_complex_math.h>

//Kronecker delta symbol
double kdel(long i, long j){
  double dlta = 0.0;
  if(i==j) dlta = 1.0;
  return dlta;
}

//Evaluates the lattice size delta. See writeup for formula
double evaluate_d(void *params){
  pot_amps *p = (pot_amps *) params;
  double vx = p->vx;
  double vxbar = p->vxbar;
  double vy = p->vy;
  
  double d = vx*vy;
  d = sqrt(d);
  d = 2.0 * d;
  d = d/(vx-vxbar);
  d = acos(d);
  return d;
}
//Evaluates the reciprocal lattice size alpha
double evaluate_alpha (double delta){
  
  double alpha;
  alpha = delta/M_PI;
  alpha = 1.0 - alpha;
  alpha = 2.0 * alpha;
  alpha = 1.0/alpha ;
  return alpha;
  
}

//Evaluates if the momentum is inside the Brillouin zone
//Returns 1 if succes, 0 if failure
int inside_fbz(const double kvec[], double delta){
  int result=0;
  double kx=kvec[0];
  double ky=kvec[1];
  
  //Evalualte reciprocal lattice size
  double alpha = evaluate_alpha(delta);
    
  //The lines of the FBZ of a trigonal lattice 
  //as shown in writeup
  double line1 = alpha * ky + kx - alpha ;
  double line2 = alpha * ky - kx - alpha ;
  
  double line1p = alpha * ky + kx + alpha;
  double line2p = alpha * ky - kx + alpha;
  
  if((line1<=0.0)&&(line2<=0.0)){
    if((line1p>=0.0)&&(line2p>=0.0)){
      result=1;
    }
  }
  
return result ;
}

int hexagonality(void *params){
    int hexag=0;  
    pot_amps *p = (pot_amps *) params;
    double vxbar=p->vxbar;
    double vx=p->vx;
    double vy=p->vy;
    double diff;
    
    diff = 4.0 * sqrt( vx * vy );
    diff = vx + diff;
    diff = vxbar - diff;
   
    if(diff>=0.0) hexag=1;
    return hexag;
}

double pot(const double xvec[], void *params){
  double vr;
  double x = xvec[0];
  double y = xvec[1];
  pot_amps *p = (pot_amps *) params;
  double vxbar=p->vxbar;
  double vx=p->vx;
  double vy=p->vy;
  
  double argx,argy, argxy;
  
  double cosx = cos(x);
  double cosy = cos(y);
  double sinx = sin(x);
  
  argxy =  2.0 * sqrt(vx * vy) * cosx * cosy;
  
  argx = vxbar * sinx * sinx;
  argx = argx + vx * cosx * cosx;
  
  argy = vy * cosy * cosy;
   
  vr = argx + argy + argxy ;
  return vr;
}

//Hamiltonian matrix element from writeup
double get_hamiltonian_real(long beta, long betap, const double r[], const double k[], void *xgrid, void *params){
    grid *g = (grid *) xgrid;
    double d = g->inc;
    double n = g->size;
       
    double kx = k[0];
    double ky = k[1];
    double ksq = kx * kx + ky * ky;
    
    double result;
    double vbeta = pot(r,params);
    
    double cons = vbeta - (ksq/2.0);
    cons = cons * kdel(beta, betap);
    
    long incr = 2 * n;
    double nablasq = kdel(beta, betap+2);
    nablasq = nablasq + kdel(beta, betap-2);
    nablasq = nablasq + kdel(beta, betap+incr);
    nablasq = nablasq + kdel(beta, betap-incr);
    nablasq = nablasq - ( 4.0 * kdel(beta, betap) );
    nablasq = nablasq / (8.0 * d * d);
    
    result = cons + nablasq;
    result = -result;
    return result;
}

double get_hamiltonian_imag(long beta, long betap, const double r[], const double k[], void *xgrid,void *params){
    grid *g = (grid *) xgrid;
    double d = g->inc;
    double n = g->size;
        
    double kx = k[0];
    double ky = k[1];
    double result;
    
    double kxterm, kyterm;
    
    kxterm = kdel(beta, betap+1) - kdel(beta, betap-1);
    kxterm = kxterm * kx / (2.0 * d);
    
    kyterm = kdel(beta, betap+n) - kdel(beta, betap-n);
    kyterm = kyterm * ky / (2.0 * d);
    
    result = kxterm + kyterm ;
    result = -result;
    return result;  
}
  
#undef __FUNCT__
#define __FUNCT__ "main"
int
main (int argc, char **argv){
  time_t begin, end;
  pot_amps params;
  grid xgrid, kgrid;
  int argv_count=1;
  int nprocs, pid;
  double em, ep;
    
  /********************Input block***********************************************/
  FILE *outfilepath = fopen (argv[argv_count++], "w");
   //Input potential parameters
  params.vxbar = atof (argv[argv_count++]);	
  params.vx = atof (argv[argv_count++]);	
  params.vy = atof (argv[argv_count++]);
  //Cutoffs
  double xcutoff = atof (argv[argv_count++]);
  long xgridsize = atoi (argv[argv_count++]);
  //kcutoff should be alpha 
  long kgridsize = atoi (argv[argv_count++]);
  /********************End Input block******************************************/
  
  GArray *kxvals_loc = g_array_new (FALSE, FALSE, sizeof (double));
  GArray *local_data = g_array_new (FALSE, FALSE, sizeof (double));
  
  //start mpi
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &pid);
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
  (void) time (&begin);
  
  if (pid==0) printf("\n======================New Param Set:=============================");
  if((hexagonality(&params)==0)&&(pid==0)){
      printf("\n Error. Not in hexagonal region of the parameter space, exiting...\n");    
      exit(1);
    }
  //Calculate lattice size delta
  double delta = evaluate_d(&params);
  //Evaluate reciprocal lattice size
  double alpha = evaluate_alpha(delta);
  
  //The 3 lattice vectors between whom the hopping will be evaluated
  double r0[2], r1[2],r2[2], rdiff1[2], rdiff2[2];
  r0[0] = M_PI-delta;
  r0[1] = 0.0;
  
  r1[0] = 2.0 * delta - M_PI;
  r1[1] = M_PI ;
    
  r2[0] = delta - M_PI;
  r2[1] = 0.0;
  
  rdiff1[0] = r1[0]-r0[0]; 
  rdiff1[1] = r1[1]-r0[1];
  
  rdiff2[0] = r2[0]-r0[0];
  rdiff2[1] = r2[1]-r0[1];
  
  double r[2],x,y;
  double k[2],kx,ky;
  
  //Build range of crystal momenta. Make sure that it contains the FBZ parallelepiped
  //as shown in writeup
  kgrid.min = -alpha;
  kgrid.max = alpha;
  kgrid.inc  = (kgrid.max-kgrid.min)/(kgridsize-1.0);
  kgrid.size = kgridsize;
  
  //Build range in real space
  xgrid.min = -xcutoff;
  xgrid.max = xcutoff;
  xgrid.inc = (xgrid.max-xgrid.min)/(xgridsize-1.0);
  xgrid.size = xgridsize;
  
  double local_outdata[4];
  long i=0;
  //Assign kx's to each process by block cyclic decomposition
  for(kx=kgrid.min;kx<=kgrid.max;kx=kx+kgrid.inc){
   //assign kx to the ith process
   if(pid==i){
     g_array_append_val(kxvals_loc,kx);
  }
   if(i == (nprocs-1)) i=0;
   else i++;
  }
  //Parallelized this loop
  //Initialize kx locally
  long kxcount=0;
  long p,q,pp,qp;
  long beta, betap;
  double helem_real, helem_imag;
  gsl_complex  helem;
  
  long hamilt_size = xgrid.size * xgrid.size;
  double tau1_loc = 0.0;
  double tau2_loc = 0.0;
  double tau1_iter, tau2_iter;
  double mag, phase;
  double tau1, tau2;
  long evalstride;
  
  while(kxcount<kxvals_loc->len){
    kx=g_array_index(kxvals_loc,double,kxcount);
    for(ky=kgrid.min;ky<=kgrid.max;ky=ky+kgrid.inc){
	 k[0]=kx;
	 k[1]=ky;
	 if(inside_fbz(k,delta)==1)
	 {
	   gsl_matrix_complex *hamilt = gsl_matrix_complex_alloc(hamilt_size,hamilt_size);
	   gsl_vector *evals = gsl_vector_alloc (hamilt_size);
	   p=0;
	   while(p<xgrid.size){
	      q=0;
	      while(q<xgrid.size){
	      x = xgrid.min + p * xgrid.inc;
	      y = xgrid.min + q * xgrid.inc;
	      r[0] = x;
	      r[1] = y;
	      beta = p + xgrid.size * q;
	      pp=0;
	      while(pp<xgrid.size){
		qp=0;
		 while(qp<xgrid.size){
		  betap = pp + xgrid.size * qp;
		  helem_real = get_hamiltonian_real(beta, betap,r,k,&xgrid,&params);
		  helem_imag = get_hamiltonian_imag(beta, betap,r,k,&xgrid,&params);
		  GSL_SET_COMPLEX(&helem, helem_real, helem_imag);
		  gsl_matrix_complex_set(hamilt, beta, betap, helem);
		  qp++;
		 }//Endwhile qp
	       pp++;
	      }//Endwhile pp
	     q++;
	     }//Endwhile q
	   p++;
	   }//Endwhile p
	 
	  //Diagonalize the hamiltonian and get the 2 lowest eigenvalues and 
	  //corresponding eigenfunctions
	  gsl_eigen_herm_workspace *w = gsl_eigen_herm_alloc(hamilt_size);
	  gsl_eigen_herm(hamilt,evals,w);
  	  gsl_sort_vector (evals);
	  //Choose the first nonpositive eigenvalue and the one after that
	  evalstride=0;
	  em =  gsl_vector_get(evals,evalstride);
	  ep	 =   gsl_vector_get(evals,evalstride+1);
	  while(evalstride<hamilt_size-1){
	   em = gsl_vector_get(evals,evalstride);
	   evalstride++;
	   ep = gsl_vector_get(evals,evalstride);
	  if( (em<=0) && (ep>=0)) break;
	  }
	  //Gather local out data to array
	  local_outdata[0]=kx;
	  local_outdata[1]=ky;
	  local_outdata[2]=em;
	  local_outdata[3]=ep;
	  //Data is stored as kx,ky,em,ep,kx,ky,em,ep...   
	  g_array_append_vals(local_data,&local_outdata,4);
	  gsl_eigen_herm_free (w);
	  gsl_matrix_complex_free(hamilt);
	  gsl_vector_free(evals);
	 //Only choose positive kys as the top half of the FBZ
	 if(ky>=0.0){ 
	 mag = ( em - ep );
	 //Append to local sum of taus
	 phase = kx * rdiff1[0] + ky * rdiff1[1];
	 tau1_iter = mag * cos(phase);
	 tau1_loc = tau1_loc + tau1_iter;
	 
	 phase = kx * rdiff2[0] + ky * rdiff2[1];
	 tau2_iter = mag * cos(phase);;
	 tau2_loc = tau2_loc+tau2_iter;	 
	 }  
	}
      }
    kxcount++;    
  }
  
  //Aggregate the local taus by reduction
  MPI_Allreduce (&tau1_loc, &tau1, 1, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);
  MPI_Allreduce (&tau2_loc, &tau2, 1, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);  
  
  tau1 = tau1/(kgrid.size * kgrid.size);
  tau2 = tau2/(kgrid.size * kgrid.size);
  
  if(pid==0) printf("\n\n Done. Dumping output...");
  //Dump out global_data to file
  //aggregate the local data sizes to root
  int local_datasize=local_data->len;
  int global_datasize;
  
  MPI_Allreduce (&local_datasize, &global_datasize, 1, MPI_INT, MPI_SUM,
		  MPI_COMM_WORLD);
  
  double *global_data;
  //allocate space for global data
  global_data = malloc(global_datasize * sizeof(*global_data));
  //Aggregate each local size to root
  int *datasizes;
  datasizes = malloc(nprocs * sizeof(*datasizes));
  
  MPI_Allgather (&local_datasize, 1, MPI_INT, datasizes, 1, MPI_INT,
		  MPI_COMM_WORLD);
  
  //Calculate data displacement for each process
  int *displacements;
  displacements = malloc(nprocs * sizeof(*displacements));
  int displ=0, dcount;
  if (pid != 0)
    for (dcount = 0; dcount < pid; dcount++)
      {
	displ = displ + datasizes[dcount];
      }
  //Aggregate the data displacements to root
  MPI_Allgather (&displ, 1, MPI_INT, displacements, 1, MPI_INT,
		   MPI_COMM_WORLD);
  
  //Now aggregate the actual data to root with all of the above information
  MPI_Gatherv ((double *) (void *) local_data->data, local_datasize,
		 MPI_DOUBLE, global_data, datasizes, displacements,
		 MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  //Now dump output to file 
  int datastride=0;
  if(pid==0){
    while(datastride<global_datasize){
      fprintf(outfilepath,"\n %lf",global_data[datastride++]);//kx
      fprintf(outfilepath," %lf",global_data[datastride++]);//ky
      fprintf(outfilepath," %lf",global_data[datastride++]);//em
      fprintf(outfilepath," %lf",global_data[datastride++]);//ep
    }
  }
  
  if(pid==0) {
    printf("\n\n Parameters:");
    printf("\n vxbar = %-3.12lf, vx = %-3.12lf, vy = %-3.12lf, delta = %-3.12lf",params.vxbar,params.vx,params.vy, delta);
    printf("\n Hopping Amplitudes:");
    printf("\n tau1 = %-3.12lf, tau2 = %-3.12lf", tau1, tau2);
    printf("\n\n Cleaning up....");
  }
  
  g_array_free(kxvals_loc,TRUE);
  g_array_free(local_data,TRUE);
  free(global_data);
  free(datasizes);
  free(displacements);

  (void) time (&end);
  if (pid==0) printf("\n\n Done. Time taken = %14.2lf sec\n\n",difftime(end,begin));
  MPI_Finalize ();
  return 0;
}
