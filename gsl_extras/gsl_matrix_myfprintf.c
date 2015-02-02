#include "gsl_matrix_myfprintf.h"

int
my_gsl_matrix_fprintf (FILE * stream, gsl_matrix * m, char *fmt)
{
  size_t rows = m->size1;
  size_t cols = m->size2;
  size_t row, col, ml;
  int fill;
  char buf[100];
  gsl_vector *maxlen;

  maxlen = gsl_vector_alloc (cols);
  for (col = 0; col < cols; ++col)
    {
      ml = 0;
      for (row = 0; row < rows; ++row)
	{
	  sprintf (buf, fmt, gsl_matrix_get (m, row, col));
	  if (strlen (buf) > ml)
	    ml = strlen (buf);
	}
      gsl_vector_set (maxlen, col, ml);
    }

  for (row = 0; row < rows; ++row)
    {
      for (col = 0; col < cols; ++col)
	{
	  sprintf (buf, fmt, gsl_matrix_get (m, row, col));
	  fprintf (stream, "%s", buf);
	  fill = gsl_vector_get (maxlen, col) + 1 - strlen (buf);
	  while (--fill >= 0)
	    fprintf (stream, " ");
	}
      fprintf (stream, "\n");
    }
  gsl_vector_free (maxlen);
  return 0;
}
