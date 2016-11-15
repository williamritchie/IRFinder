/* To compile cc -o winflat_with_beta winflat_with_beta.c -lm */ 
/* Copyright Stephane Audic 2003 */ 
/* With code included from Numerical Recipes in C */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>



#define MAXIT 500
#define EPS 3.0e-30
#define FPMIN 1.0e-30

double betacf(double a, double b, double x);
double gammln(double xx);
double betai(double a, double b, double x);

double betacf(double a, double b, double x)
{
  /*	void nrerror(char error_text[]); */
	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;
	
	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (m > MAXIT){
	  fprintf( stderr , "a or b too big, or MAXIT too small in betacf");
	  exit(1) ; 
	}
	return h;
}
#undef MAXIT
#undef EPS
#undef FPMIN


double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}



double betai(double a, double b, double x)
{
	double bt;

	if (x < 0.0 || x > 1.0) { 
	  fprintf( stderr , "Bad x in routine betai") ; 
	  exit(1) ; 
	}
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*betacf(a,b,x)/a;
	else
		return 1.0-bt*betacf(b,a,1.0-x)/b;
}







void usage(){
    fprintf( stderr , "usage: two behaviours can be invoqued. \n") ; 
    fprintf( stderr , "  winflat -xvalue x -sig significance [ -diff n1 n2 ]\n");
    fprintf( stderr , "  will return the lower and upper y value at the given significance level\n");
    fprintf( stderr ,"\n") ; 

    fprintf( stderr , "  winflat -xvalue x -yvalue y  [ -diff n1 n2 ]\n");
    fprintf( stderr , "  will return the probability of over or underexpression \n");
    

    /* fprintf( stderr ,"\n") ; 
    fprintf( stderr , "  winflat -xvalue x -show  [ -diff n1 n2 ]\n");
    fprintf( stderr , "  will return a plot \n");
    */ 

    fprintf( stderr ,"If the number of clones in the two libraries is \n") ; 
    fprintf( stderr ,"different, use -diff n1 n2\n") ; 
/*    fprintf(stderr  ,"To see p(y|x) and C(y|x) , use -show \n") ;  */ 

    fprintf( stderr ,"\n") ; 
}

int main( int argc , char **argv){
  double temp ; 
  double thisproba , thisproba2 ;
  int x ; 
  int y ; 
  int thisy ; 
  int n1 = 1  , n2 = 1  ; 
  double ratio ; 
  double t1 , t2 , t3 , t4 ; 

  double sum ; 
  int ymin , ymax ; 
  int argcount = 1 ; 
  int noup ; 
  double  sig ; 
  int show = 0 ; 

  double p ; 
#define WNMIN -1 
#define WNMAX  10000

  if( argc < 4){
    usage() ; 
    exit(0) ; 
  }

  x = -1 ; 
  y = -1 ; 

  while( argcount < argc ){

    if( !strcmp( argv[argcount] , "-xvalue")){
      x = atoi( argv[argcount+1] ) ; 
      argcount += 2 ; 
    } else if( !strcmp( argv[argcount] , "-yvalue")){
      y = atoi( argv[argcount+1] ) ; 
      argcount += 2 ; 
    } else if( !strcmp( argv[ argcount ] , "-diff")){
      /* correction in the case where the number of est's drawn 
	 from the samples  is different */
      n1 = atoi( argv[ argcount + 1 ] ) ; 
      n2 = atoi( argv[ argcount + 2 ] ) ; 
      argcount += 3 ; 
    } else if( !strcmp( argv[ argcount ] , "-sig" )){
      sig = (double) atof( argv[ argcount + 1 ] ) ; 
      argcount += 2 ; 
    } else {
      usage() ; 
      exit(0) ; 
    }
  }



    /* Check arguments and invoque the right procedure */ 
    
    if( x > -1 && y > -1 ){
	/* Both x and y are defined so we compute the significance window */
	ymin =  WNMIN ;
	ymax =  WNMAX ;
  
	sum = 0 ; 
	noup = 1 ; 
	ratio = (double) n1 / (double) n2 ; 

	p = (double) ( n1 ) / (double) ( n1 + n2 )  ; 

	thisproba = betai( (double) (x + 1 ) , (double)( y + 1 ) , p ) ; 
	thisproba2 = betai(  (double)  ( y + 1 )  ,  (double) (x + 1 ) , 1 - p ) ; 
	fprintf( stdout , "P( y <= %d | x = %d ) = %g \n" , y , x , thisproba ) ; 
	fprintf( stdout , "P( y >= %d | x = %d ) = %g \n" , y , x , thisproba2 ) ; 
    } else if( x > -1 && y == -1 ){

	y = 0 ; 
	ymin =  WNMIN ;
	ymax =  WNMAX ;
  
	sum = 0 ; 
	noup = 1 ; 
	ratio = (double) n1 / (double) n2 ; 

	p = (double) ( n1 ) / (double) ( n1 + n2 )  ; 


	y = 0 ; 
	/* fprintf( stderr , "x = %d y= %d sig = %g\n" , x , y , sig ) ;  */ 
	while( 1 ){
	    thisproba = betai( (double) (x + 1 ) , (double)( y + 1 ) , p ) ; 
	    thisproba2 = betai(  (double)  ( y + 1 )  ,  (double) (x + 1 ) , 1 - p ) ; 

	    /* fprintf( stdout , "%d %d C(%d | %d ) = %g    %g\n" , 
		     x , y , y , x , thisproba , thisproba2 ) ; 
	    */ 

	    if( thisproba < sig / 2.0 ){
		ymin = y ; 
	    }

	    if( ymax == WNMAX && thisproba2 < sig / 2.0 ){
		ymax = y ; 
		break ; 
	    }
	    y++ ; 
	}
	thisproba = betai( (double) (x + 1 ) , (double)( ymin + 1 ) , p ) ; 
	thisproba2 = betai(  (double)  ( ymax + 1 )  ,  (double) (x + 1 ) , 1 - p ) ; 
	fprintf( stderr , "P( y <= %d | x = %d ) = %g \n" , ymin , x , thisproba ) ; 
	fprintf( stderr , "P( y >= %d | x = %d ) = %g \n" ,  ymax , x , thisproba2 ) ; 
	
	if( ymin == -1 ){
	    fprintf( stdout , "%d *--%d\n" , x , ymax) ;
	} else {
	    fprintf( stdout , "%d %d--%d\n" , x , ymin , ymax) ;
	}
	
    } else if( show == 1 ){
	
	fprintf( stdout , "%d %d C(%d | %d ) = %g    %g\n" , 
		 x , y , y , x , thisproba , thisproba2 ) ; 
	y++ ; 
	if( thisproba > 0.9999999999 ){

	}
    }
}











