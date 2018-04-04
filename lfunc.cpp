#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define NUMPANELS 1000000
#define DEG 2

double **cof;
void polint(double *xa, double *ya, int n, double x, double *y, double *dy);
void polcof(double *xa, double *ya, int n, double *cof);

double lfunc(double z)
{
  double x,y;

  x=exp(-z);
  y = 1 + x;
  return log(y) - ((y-1)-x)/y ;  /* cancels errors with IEEE arithmetic */
}


void initlfunc2()
{
  int i,j;
  double xa[DEG+1],ya[DEG+1];
  double a,b;

  //Initialise coeff for lfunc2
  cof=(double**)malloc((NUMPANELS+1)*sizeof(double*));
  if(!cof){printf("Could not allocate cof in main.\n");exit(1);}
  for(j=0;j<=NUMPANELS;j++){
    cof[j]=(double*)malloc((DEG+1)*sizeof(double));
    if(!cof[j]){printf("Could not allocate cof in main.\n");exit(1);}
  }

  for(j=0;j<NUMPANELS;j++){
    a=1.0*j/NUMPANELS;
    b=1.0*(j+1)/NUMPANELS;
    for(i=0;i<=DEG;i++){
      xa[i]=cos((2*i+1)*3.141592653589793238/(2*DEG+2));
      xa[i]=0.5*((b-a)*xa[i]+a+b);
      ya[i]=lfunc(1/xa[i]-1);
    }

    polcof(xa,ya,DEG,cof[j]);
  }

  cof[NUMPANELS][0]=log(2.0);
  for(i=1;i<=DEG;i++)
    cof[NUMPANELS][i]=0.0;
}


double lfunc2(double x)
{
  double temp,y;
  int i;       
  double *coeff;        

  y=1/(1+x);
  coeff=cof[(int)(NUMPANELS*y)];
  
  temp=coeff[DEG];
  for(i=DEG-1;i>=0;i--)
    temp=temp*y+coeff[i];
  return temp;
}


void cleanuplfunc2()
{
  int j;

  for(j=0;j<=NUMPANELS;j++)
    free(cof[j]);
  free(cof);
}


double lnsum(double a,double b)
{
  double x=a-b;
  if(x<0){x=-x;a=b;}
  return a+lfunc(x);
}


double lnsum2(double a,double b)
{
  double x=a-b;
  if(x<0){x=-x;a=b;}
  return a+lfunc2(x);
}


