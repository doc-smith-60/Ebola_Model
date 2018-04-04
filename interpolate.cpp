#include<stdio.h>
#include<math.h>
#include<stdlib.h>

void polint(double *xa, double *ya, int n, double x, double *y, double *dy)
//Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y, and
//an error estimate dy. If P(x) is the polynomial of degree N - 1 such that P(xai) = yai, i =
//1, . . . , n, then the returned value y = P(x).
{
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  double *c,*d;

  dif=fabs(x-xa[1]);
  c=(double*)malloc(n*sizeof(double))-1;
  if(!c){printf("Could not allocate c in polint.\n");exit(1);}
  d=(double*)malloc(n*sizeof(double))-1;
  if(!d){printf("Could not allocate d in polint.\n");exit(1);}
  for (i=1;i<=n;i++) { //Here we find the index ns of the closest table entry,
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i]; //and initialize the tableau of c’s and d’s.
    d[i]=ya[i];
  }
  *y=ya[ns--]; //This is the initial approximation to y.
  for (m=1;m<n;m++) { //For each column of the tableau,
    for (i=1;i<=n-m;i++) { //we loop over the current c’s and d’s and updatethem. 
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ( (den=ho-hp) == 0.0) {printf("Error in routine polint");exit(1);}
      //This error can occur only if two input xa’s are (to within roundoff) identical.
      den=w/den;
      d[i]=hp*den; //Here the c’s and d’s are updated.
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
//After each column in the tableau is completed, we decide which correction, c or d,
//we want to add to our accumulating value of y, i.e., which path to take through the
//tableau—forking up or down. We do this in such a way as to take the most “straight
//line" route through the tableau to its apex, updating ns accordingly to keep track of
//where we are. This route keeps the partial approximations centered (insofar as possible)
//on the target x. The last dy added is thus the error indication.
  }
  free(d+1);
  free(c+1);
}


void polcof(double *xa, double *ya, int n, double *cof)
//Given arrays xa[0..n] and ya[0..n] containing a tabulated function yai = f(xai), this
//routine returns an array of coefficients cof[0..n] such that yai =  j cofjxaji .
{
  void polint(double *xa, double *ya, int n, double x, double *y, double *dy);
  int k,j,i;
  double xmin,dy,*x,*y;

  x=(double*)malloc((n+1)*sizeof(double));
  if(!x){printf("Could not allocate x in polcof.\n");exit(1);}
  y=(double*)malloc((n+1)*sizeof(double));
  if(!y){printf("Could not allocate y in polcof.\n");exit(1);}

  for (j=0;j<=n;j++) {
    x[j]=xa[j];
    y[j]=ya[j];
  }
  for (j=0;j<=n;j++) {
    polint(x-1,y-1,n+1-j,0.0,&cof[j],&dy);
    //Subtract 1 from the pointers to x and y because polint uses dimensions [1..n]. We
    //extrapolate to x = 0.
    xmin=1.0e38;
    k = -1;
    for (i=0;i<=n-j;i++) { //Find the remaining xi of smallest absolute value, 
      if (fabs(x[i]) < xmin) {
        xmin=fabs(x[i]);
        k=i;
      }
      if (x[i]) y[i]=(y[i]-cof[j])/x[i]; //(meanwhile reducing all the terms)
    }
    for (i=k+1;i<=n-j;i++) { //and eliminate it.
      y[i-1]=y[i];
      x[i-1]=x[i];
    }
  }
  free(y);
  free(x);
}
