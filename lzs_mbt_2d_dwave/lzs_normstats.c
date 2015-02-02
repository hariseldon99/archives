#include "lzs_normstats.h"
#include "params_ics.h"
void
norm_stats (void *normstats, GArray * a, double normexpect)
{
  stats *s = (stats *) normstats;
  //Extract normdata from gArray. It's the last data out of every 6, counting from 0
  long i;
  double *normdata;
  long asize = a->len;
  long normsize = asize / 6;
  normdata = (double *) calloc (normsize, sizeof (double));
  for (i = 0; i < normsize; i++)
    normdata[i] = g_array_index (a, double, 6 * i + 5);
  s->mean = gsl_stats_mean (normdata, 1, normsize);
  s->stddev = gsl_stats_sd_m (normdata, 1, normsize, normexpect);
  s->max = gsl_stats_max (normdata, 1, normsize);
  s->min = gsl_stats_min (normdata, 1, normsize);
  free (normdata);
}

//compute average probability from regular time evolution data
double
prob_avg (GArray * a)
{

  //Extract probdata from gArray. It's the 5th data out of every 6, counting from 0
  long i = 1, time_index, prob_index;
  double prob = 0.0, probavg, temp;
  long asize = a->len;
  long probsize = asize / 2;

  while (i < asize)
    {
      time_index = i++;
      prob_index = i++;
      temp = g_array_index (a, double, prob_index);
      prob = prob + temp;
    }
  probavg = prob / probsize;
  return probavg;

}
