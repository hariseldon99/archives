/*
*  C Implementation: final_output
*
* Description:
*
*
* Author: Analabha Roy <daneel@physics.utexas.edu>, (C) 2008
*
* Copyright: See COPYING file that comes with this distribution
*
*/
#include <stdio.h>
#include <time.h>
#include "final_output.h"

void
final_output (double period, time_t begin, time_t end, int nprocs)
{
  printf
    ("\n=========================================================================");
  printf ("\nTIME PERIOD OF DRIVE = %lf", period);
  printf
    ("\n=========================================================================");
  printf
    ("\n                 PROGRAM RAN SUCCESSFULLY (I THINK)                      ");
  printf
    ("\n                    TIME ELAPSED = %ld sec                               ",
     end - begin);
  printf
    ("\n                   NUMBER OF PROCESSES = %d                               ",
     nprocs);
  printf
    ("\n=========================================================================\n");
}
