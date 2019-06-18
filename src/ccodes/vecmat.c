# include <math.h>
# include <Rmath.h>
# include <stdio.h>
# include <R.h>
# include <string.h>

/***********************************************************************/
void maxa (int la[1],double a[], double m[1])
{
    int i;
    m[0] = a [0];
    for(i = 1; i < la[0]; i++) if(a[i] > m[0]) m [0] = a[i];
}


void log_sum_exp (int la[1], double a[], double o[1])
{  int i;
   double m[1];

   maxa (la, a, m);

   o[0] = 0; for(i = 0; i < la[0]; i++)
   {
      o[0] += exp(a[i] - m[0]);
   }
   
   o[0] = log(o[0]) + m[0];
}


void colSums (int r[1], int c[1], double A[r[0]][c[0]], double S[c[0]])
{
  int i, j;
  for (j = 0; j < c[0]; j++) {
    S[j] = 0; for (i = 0; i < r[0]; i++) S[j] += A[i][j];
  }
}

void cpvec (int n, double A[], double B[])
{
    int i;
    for (i = 0; i < n; i++)
    {
        B [i] = A [i];
    }
}

void setvec (int n, double A[], double a)
{
    int i;
    for (i = 0; i < n; i++)
    {
        A[i] = a;
    }
}

