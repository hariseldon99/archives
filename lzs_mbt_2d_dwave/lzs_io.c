#include "lzs_io.h"
#include "params_ics.h"
//dump a quantity to stdout
void dump (const char quantity[], double value){
  printf("\n Dumping quantity:"); 
  printf("%s :", quantity);
  printf(" %lf",value);
}

int
copyvec (double source[], double target[], long size)
{
  long i;
  int outcome;
  if (2 * size < sizeof (target))
    outcome = GSL_FAILURE;
  else
    outcome = GSL_SUCCESS;

  for (i = 0; i < size; i++)
    target[i] = source[i];

  return outcome;
}


int
arg_err (int argc)
{

  printf ("\n Error, incorrect number of arguments %d\n\n", argc);
  printf ("\n Please see the 'runlzs_*' scripts for details.\n\n");
  printf ("\n Exiting ...\n\n");
  return GSL_FAILURE;

}

void
begin_out ()
{
  printf
    ("\n================================New Param Set:======================================================\n");

}

void
init_out (void *initconds)
{
  initconds_ic *ic = (initconds_ic *) initconds;
  double uk_real = ic->phi_in[0];
  double uk_imag = ic->phi_in[2];
  double vk_real = ic->phi_in[1];
  double vk_imag = ic->phi_in[3];
  double ts = ic->t_init;
  double tt = ic->t_periods;
  double td = ic->t_extra;

  double tf;
  tf = ts + tt + td;

  printf ("\n");
  puts ("Initial Conditions for lzs_tls run:");
  printf
    ("\n-------------------------------------------------------------------------------------------");
  printf
    ("\n|Initial state:  | Start time ts  | # of periods n | Offset time td | Total time tf| ");
  printf
    ("\n-------------------------------------------------------------------------------------------");
  printf
    ("\n|  %2.2lf + I %2.2lf | %-3.9lf | %-3.9lf |%-3.13lf | %-3.9lf|\n",
     uk_real, uk_imag, ts, tt, td, tf);
  printf ("|  %2.2lf + I %2.2lf |", vk_real, vk_imag);
  printf
    ("\n-------------------------------------------------------------------------------------------\n");
}


void
lzs_intgerr (int p)
{
  printf ("\n GSL execution of integration failed, bailing..");
}

void
params_out (void *params)
{
  paramspace_pt *p = (paramspace_pt *) params;
  double delta_real = p->delta_real;
  double delta_imag = p->delta_imag;
  double w = p->omega;
  double veff = p->veff;
  double off = p->offset;
  double mu0 = p->mu0;
  double muamp = p->muamp;
  char stateout = p->stateout;
  char gndstate = p->gndstate;

  double gapsq = delta_real * delta_real + delta_imag * delta_imag;
  puts ("\nSystem parameters for lzs_tls run:");
  printf
    ("\n-----------------------------------------------------------------");
  printf
    ("\nveff     |delta    |  omega   | ampof    | mut=0    |  mu_a     |");
  printf
    ("\n-----------------------------------------------------------------\n");
  printf ("%3.6lf |%-3.6lf | %-3.6lf | %-3.6lf | %-3.6lf | %-3.6lf |\n", veff,
	  sqrt (gapsq), w, off, mu0, muamp);
  printf
    ("-----------------------------------------------------------------\n");
  if (stateout == 'y')
    {
      printf ("\nFull state of the system will be dumped to binary file.");
    }
  else
    {
      printf
	("\nFull state of the system will NOT be dumped to binary file.");
    }

  printf ("\nIn the beginning, the state of the system is set to : ");
  if (gndstate == 'd')
    {
      printf ("|0>");
    }
  else if (gndstate == 's')
    {
      printf ("|BCS state>");
    }
  else if (gndstate == 'r')
    {
      printf ("Random but normalized complex numbers .");
    }
   else 
   {
     printf ("Unknown, Bad options.");
  }
}

void
begin_tout ()
{
  printf
    ("-----------------------------------------------------------------------------------------------------------------");
  printf
    ("\n    Time:          |                        phi                           |      norm\n");
  printf
    ("-----------------------------------------------------------------------------------------------------------------");
}

void
data_out (double t, double *phi, double phi_norm)
{

  printf ("\n %-3.15lf | %-3.8lf +I %-3.8lf, %-3.8lf +I %-3.8lf | %-3.8lf\n",
	  t, phi[0], phi[2], phi[1], phi[3], phi_norm);

}


void
dump_out (GArray * a, FILE * filepointer)
{
  long i = 0, time_index, uk_real_index;
  long vk_real_index, uk_imag_index;
  long vk_imag_index, norm_index;

  while (i < a->len)
    {
      time_index = i++;
      uk_real_index = i++;
      vk_real_index = i++;
      uk_imag_index = i++;
      vk_imag_index = i++;
      norm_index = i++;
      fprintf (filepointer,
	       "%03.8lf %03.8lf %03.8lf %03.8lf %03.8lf %03.8lf\n",
	       g_array_index (a, double, time_index), g_array_index (a,
								     double,
								     uk_real_index),
	       g_array_index (a, double, vk_real_index), g_array_index (a,
									double,
									uk_imag_index),
	       g_array_index (a, double, vk_imag_index), g_array_index (a,
									double,
									norm_index));

    }

}

//dumps out strobed data
void
dump_out_strobe (double scaled_time, double mag, double fidelity,
		 double delta, double defect_density, double resen,
		 double fidsus)
{
  printf
    ("\n\n------------------------------------------------------------------------------------------------------------------");
  printf ("\n Data strobed at time wt/2pi = %lf:", scaled_time);
  printf
    ("\n------------------------------------------------------------------------------------------------------------------");
  printf
    ("\n magnetization_s   |     fidelity_s   | absolute_delta_s | defect_density_s | residual_energy_s |    fid_suscp_s   |");
  printf
    ("\n------------------------------------------------------------------------------------------------------------------");
  printf
    ("\n %-3.14lf | %-3.14lf | %-3.14lf | %-3.14lf | %-3.14lf | %-3.14lf |",
     mag, fidelity, delta, defect_density, resen, fidsus);
  printf
    ("\n------------------------------------------------------------------------------------------------------------------\n");
}



void
dump_out_probs (GArray * a, FILE * filepointer, void *params)
{
  long i = 0, time_index, prob_index;
  paramspace_pt *p = (paramspace_pt *) params;
  double w = p->omega, time_dimless;
  while (i < a->len)
    {
      time_index = i++;
      prob_index = i++;
      time_dimless = g_array_index (a, double, time_index);
      time_dimless = w * time_dimless / (2.0 * M_PI);
      fprintf (filepointer, "%03.18lf %03.18lf \n", time_dimless,
	       g_array_index (a, double, prob_index));

    }

}

void
dump_out_delta (GArray * a, FILE * filepointer, void *params)
{
  long i = 0, time_index, deltareal_index, deltaimag_index;
  paramspace_pt *p = (paramspace_pt *) params;
  double w = p->omega, time_dimless;
  while (i < a->len)
    {
      time_index = i++;
      deltareal_index = i++;
      deltaimag_index = i++;
      
      time_dimless = g_array_index (a, double, time_index);
      time_dimless = w * time_dimless / (2.0 * M_PI);
      fprintf (filepointer, "%03.18lf %03.18lf %03.18lf\n", time_dimless,
	       g_array_index (a, double, deltareal_index),
			g_array_index(a,double, deltaimag_index));

    }

}


void
final_out (time_t begin, time_t end, double normexpect, void *normstats)
{
  stats *s = (stats *) normstats;
  double max = s->max;
  double min = s->min;
  double mean = s->mean;
  double stddev = s->stddev;
  printf ("\n Final Output:");
  printf
    ("\n---------------------------------------------------------------------------------------------------");
  printf
    ("\n Total runtime  | Expected norm |  Min devn     |  Max devn     |    Std devn   |");
  printf
    ("\n---------------------------------------------------------------------------------------------------");
  printf ("\n %-3.12f |%-3.12lf |%-3.12lf |%-3.12lf |%-3.12lf |\n",
	  difftime (end, begin), normexpect, fabs (min - mean),
	  fabs (mean - max), stddev);
  printf
    ("---------------------------------------------------------------------------------------------------\n");
}

void
final_out_mbt (time_t begin, time_t end, double gridsize, void *qfacts)
{

  qfactors *q = (qfactors *) qfacts;
  double probavg = q->probavg;
  double fidavg = q->fidavg;
  double deltavg = q->deltavg;
  double ndavg = q->ndavg;
  double ekavg = q->ekavg;

  printf
    ("\n---------------------------------------------------------------------------------------------------");
  printf ("\n Final Output:");
  printf
    ("\n---------------------------------------------------------------------------------------------------");
  printf
    ("\n---------------------------------------------------------------------------------------------------\n");
  printf ("\n Time averages:");
  printf
    ("\n------------------------------------------------------------------------------------------------");
  printf
    ("\n magnetization   |    fidelity    | absolute delta | defect density  | residual energy|");
  printf
    ("\n------------------------------------------------------------------------------------------------");
  printf ("\n %-3.12lf | %-3.12lf | %-3.12lf | %-3.12lf  | %-3.12lf |",
	  probavg, fidavg, deltavg, ndavg, ekavg);
  printf
    ("\n------------------------------------------------------------------------------------------------");
  printf ("\n\nRuntime metrics :");
  printf ("\n----------------------------------");
  printf ("\n  total runtime |    gridsize    |");
  printf ("\n----------------------------------");
  printf ("\n %14.2f | %14.0lf |", difftime (end, begin), gridsize);
  printf ("\n----------------------------------\n");

}


void
prn_garray (GArray * a)
{

  long i, size = a->len;
  printf ("\n");
  for (i = 0; i < size; i++)
    printf ("%lf ", g_array_index (a, double, i));
  printf ("\n");
}


void
prn_array (double *a, long size)
{

  long i;
  printf ("\n");
  for (i = 0; i < size; i++)
    printf ("%lf ", a[i]);
  printf ("\n");
}

//copies a buffer of size bufsize to a GArray at its position kstride
void
copytobuffer (GArray * a, double buffer[], long bufsize, long kstride)
{

  long i, j = 0;
  for (i = kstride + 1; i <= kstride + bufsize; i++)
    buffer[j++] = g_array_index (a, double, i);
}

//Does the opposite of above ie copies bufsize # of elements of a GArray
// from its position kstride to a buffer of size bufsize
void
copyfrombuffer (double buffer[], long bufsize, GArray * a, long kstride)
{

  long i = kstride, j = 0;
  for (j = 0; j < bufsize; j++)
    ((double *) (void *) a->data)[++i] = buffer[j];

}

//Dump out the standard deviation of quantity from GArray of time_quantity
void output_stdev (GArray *time_quantity){

  long datasize = time_quantity->len;
  datasize = datasize/2;
  double first_element = g_array_index (time_quantity,double,0);
  //Remove the first element, which is presumed to be an unneeded time quantity
  g_array_remove_index (time_quantity, 0);
  double stdev = gsl_stats_sd(((double *) (void *) time_quantity->data), 2,
			      datasize);
  //Reinsert the first element
  g_array_prepend_val(time_quantity,first_element);
  printf("\n Standard Deviation of fidelity detrivative = %lf", stdev);

  
}

