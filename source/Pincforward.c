
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void Pincforward( int iw, int nx, double dx, int nk, double dk, double k, int iz,double dz, double *P0r,double *P0i){
	double pi = 3.14159265358979323846;
  	int ik,nk2=(nk-1)/2;
  	int ix,nx2=(nx-1)/2;
	double* P0kr = malloc(nk*sizeof(double) );
   	double* P0ki = malloc(nk*sizeof(double) );
	double* P0k2r = malloc(nk*sizeof(double) );
   	double* P0k2i = malloc(nk*sizeof(double) );
        double q2,q;
        double kx,ag,x;
	for(ik=0;ik<nk;ik++){
                kx=(ik-nk2)*dk;
		P0kr[ik]=0;
		P0ki[ik]=0;
		P0k2r[ik]=0;
		P0k2i[ik]=0;	
		for (ix=0;ix<nx;ix++){
			x=(ix-nx2)*dx;
			ag=kx*x;
			P0kr[ik]+=cos(ag)*P0r[(iz-1)*nx+ix]+sin(ag)*P0i[(iz-1)*nx+ix];
			P0ki[ik]+=cos(ag)*P0i[(iz-1)*nx+ix]-sin(ag)*P0r[(iz-1)*nx+ix];
			//P0kr[ik]+=cos(ag)*P0r[(iz)*nx+ix]+sin(ag)*P0i[(iz)*nx+ix];
			//P0ki[ik]+=cos(ag)*P0i[(iz)*nx+ix]-sin(ag)*P0r[(iz)*nx+ix];
		}
		P0kr[ik]*=dx/2/pi;
		P0ki[ik]*=dx/2/pi;

		q2=k*k-kx*kx;
		if (q2>0){
			q=sqrt(q2);
			ag=-q*dz;
			//ag=0;			
			P0k2r[ik]=P0kr[ik]*cos(ag)-P0ki[ik]*sin(ag);
			P0k2i[ik]=P0ki[ik]*cos(ag)+P0kr[ik]*sin(ag);
		}
		else{
			q=sqrt(-q2);
			ag=-q*dz;
			//ag=0;
			P0k2r[ik]=P0kr[ik]*exp(ag);
			P0k2i[ik]=P0ki[ik]*exp(ag);
		}
                
		
	}

        for(ix=0;ix<nx;ix++){
		x=(ix-nx2)*dx;
		P0r[iz*nx+ix]=0;
		P0i[iz*nx+ix]=0;
		for(ik=0;ik<nk;ik++){
			kx=(ik-nk2)*dk;
			ag=kx*x;
			P0r[iz*nx+ix]+=cos(ag)*P0k2r[ik]-sin(ag)*P0k2i[ik];
			P0i[iz*nx+ix]+=cos(ag)*P0k2i[ik]+sin(ag)*P0k2r[ik];
			
		}
		P0r[iz*nx+ix]*=dk;
		P0i[iz*nx+ix]*=dk;
											
	}
/*
	FILE* file2;  
  	char fname2[100];
	if((iw%50)==0){
		sprintf(fname2,"P2_iw_%d.dat", iz);
		file2 = fopen(fname2,"w");
 		for (ix=0;ix<nx;ix++){
			x=(ix-(nx-1)/2)*dx;
			fprintf(file2," %f %.12lf %.12lf  %.12lf %.12lf\n", x, P0r[iz*nx+ix],P0i[iz*nx+ix],P0kr[ix],P0ki[ix]);			    
			           
		}

       		 fclose(file2);
	}
*/	
}
