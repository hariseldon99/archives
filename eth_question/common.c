#include "common.h"

//Defines a null vector, NOT VACUUM, as all components 0
void
nullify (long vec[], const long vecsize)
{
  int i;
  for (i = 0; i < vecsize; i++)
    vec[i] = 0;
}

//Tests if a vec is null, returns GSL_SUCCESS is yes, GSL_FAILURE if no
int
isnull (const long vec[], const long vecsize)
{
  int i;
  int isnull = GSL_SUCCESS;	//Defaults to yes
  for (i = 0; i < vecsize; i++)
    if (vec[i] != 0)
      {
	isnull = GSL_FAILURE;
      }
  return isnull;
}

//Copies source vec to destination vec  
void
copyvec (const long src[], long dest[], long vecsize)
{
  long i;
  for (i = 0; i < vecsize; i++)
    dest[i] = src[i];
}

//Search for value in array a. Returns the position if found
//If not found, returns GSL_FAILURE  
long
find_index (long vec[], long vecsize, long value)
{
  long found = GSL_FAILURE;	//not found by default
  long i;
  for (i = 0; i < vecsize; i++)
    {
      if (vec[i] == value)
	{
	  found = i;		//it was found 
	}
    }
  return (found);
}

//Returns GSL_SUCCESS if the combination data 
//contains a forbidden site, otherwise returns 
//GSL_FAILURE  
int
is_forbidden (const size_t data[], const long vecsize, const long latsize)
{

  int i, forbidden = GSL_FAILURE;	//Not forbidden by default

  int halfway = ((latsize * latsize - 1) / 2);	//The halfway point in the 2D grid
  // Note: The halfway point is never forbidden
  for (i = 0; i < vecsize; i++)
    {
      if (data[i] < halfway)	//If the site is less than halfway there
	{
	  if (data[i] % (latsize + 1) == 0)	// If the site is a multiple of latsize+1 
	    //i.e if it is diagonal in the 2D grid,
	    forbidden = GSL_SUCCESS;	//it is forbidden
	}
      if (data[i] > halfway)	//If the site is more than halfway there
	{
	  if ((data[i] - 1) % (latsize + 1) == 0)	// If the site before this one 
	    //is diagonal 
	    //ie if this site is subdiagonal in the 2D grid,
	    forbidden = GSL_SUCCESS;	//it is forbidden
	}
    }
  return forbidden;
}

//Returns GSL_SUCCESS if the combination data 
//contains a site outside the lower right quarter
//, otherwise returns GSL_FAILURE  
int
is_outside_lower_right_quarter (const long vec[], const long vecsize,
				const long latsize)
{

  int i, forbidden = GSL_FAILURE;	//Inside by default
  int row_index, col_index;

  int halfway = ((latsize * latsize - 1) / 2);	//The halfway point in the 2D grid
  //If any of the elements of vec are outside the lower right quarter, output is 
  //GSL_SUCCESS
  for (i = 0; i < vecsize; i++)
    {
      //Get the row and column indices
      row_index = vec[i] % latsize;
      col_index = vec[i] - (latsize * row_index);
      //If the row index is before halfway, and column index is after row index 
      if ((col_index < halfway) && (row_index > col_index))
	{
	  forbidden = GSL_SUCCESS;	//It is outside
	  break;		//Get out
	}
      //If the row index is at or after halfway, and column index is at or after 
      //row index
      if ((col_index >= halfway) && (row_index >= col_index))
	{
	  forbidden = GSL_SUCCESS;	//it is outside
	  break;		//Get out
	}
    }
  return forbidden;
}


// Dot product between 2 canonical basis states. 
//Returns 1,0 if they are the same, else returns 0.0
//If one of the vectors is NULL, returns 0
double
can_sdot (const long vec1[], const long vec2[], long vecsize)
{

  double prod = 1.0;
  long i;
  if ((isnull (vec1, vecsize) == GSL_SUCCESS)
      || (isnull (vec2, vecsize) == GSL_SUCCESS))
    {
      prod = 0.0;
    }
  else
    {
      for (i = 0; i < vecsize; i++)
	if (vec1[i] != vec2[i])
	  {
	    prod = 0.0;
	    break;
	  }
    }

  return prod;
}

//Checks to see if 2 sites in a 2D grid of size are nearest neighbors
//Returns GSL_SUCCESS is they are, GSL_FAILURE if they are not
int
are_nn (const long i, const long j, const long latsize)
{

  int nn = GSL_FAILURE;		//By default, they are not nearest neighbors

  //They may be nearest neighbours if they are separated by 1 or latsize
  //Necessary but not sufficient
  if ((abs (i - j) == 1) || (abs (i - j) == latsize))
    nn = GSL_SUCCESS;

  //The above possibility includes edge states. Remove them 
  if (((j % latsize == 0) && (j - i == 1))
      || ((i % latsize == 0) && (i - j == 1)))
    nn = GSL_FAILURE;

  return nn;
}

//Acts b^\dagger_i b_j on vec[].
void
bdaggerb (const long i, const long j, long vec[], const long vecsize,
	  const long latsize)
{
  size_t iarray[1];
  long i_index = find_index (vec, vecsize, i);
  long j_index = find_index (vec, vecsize, j);

  //If i is unoccupied and j is occupied, the jth site moves to i,
  //unless i is forbidden, in which case  b_idaggerb_j on it is 0
  if ((i_index == GSL_FAILURE) && (j_index != GSL_FAILURE))
    {				//cast site i into an array of size 1.
      iarray[0] = i;
      if (is_forbidden (iarray, 1, latsize) == GSL_SUCCESS)	//If site i is forbidden
	nullify (vec, vecsize);
      else
	{			//Otherwise move the jth site to i and resort the combination
	  vec[j_index] = i;
	  //Now, resort the array in lexicographical order ie increasing
	  gsl_sort_long (vec, 1, vecsize);
	}
    }
  else				//In all other cases, b_j is unoccupied and so  b_idaggerb_j on it is 0
    nullify (vec, vecsize);

}

//Acts nsquared = b^\dagger_i b_i b^\dagger_j b_j  on vec[].
void
nsquared (const long i, const long j, long vec[], const long vecsize)
{

  long i_index = find_index (vec, vecsize, i);
  long j_index = find_index (vec, vecsize, j);

  if ((i_index == GSL_FAILURE) || (j_index == GSL_FAILURE))	//If either i or j sites are unoccupied
    nullify (vec, vecsize);	//b_i destroys vec and b_j destroys vec
  //Otherwise do nothing. just return the original vec
  //Remember that b^\dagger_i b_i on ith occupied site is just unity.
}
