#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Pincforward.h>


void prodcomp( double *z1, double *z2, double *zres ){

     // real part
     zres[0] = z1[0]*z2[0]-z1[1]*z2[1]; 
     // imaginary part  
     zres[1] = z1[1]*z2[0]+z1[0]*z2[1]; 

}

void greenfunction2(double x, double z, double k,  double* G00)
{
	double pi = 3.14159265358979323846; 
	double Exp[2],exp2[2],exp3[2];
	double kz2,kn;
	int n,i;
	G00[0]=0;
	G00[1]=0;
	double L=2000;
	for (i=0;i<401;i++){
		n=i-200;
		kz2=k*k-(2*pi*n/L)*(2*pi*n/L);
		Exp[0]=cos(2*pi*n*x/L);
		Exp[1]=sin(2*pi*n*x/L);
		if (kz2>0){
			kn=sqrt(kz2);			
			exp2[0]=cos(kn*fabs(z));
			exp2[1]=sin(kn*fabs(z));
			prodcomp(Exp, exp2, exp3 );
			G00[0]+=exp3[1]/kn/2/L;
			G00[1]+=exp3[0]/kn/2/L;//sign issue!!!
		}
		else{
			kn=sqrt(-kz2);
			exp2[0]=exp(-kn*fabs(z));
			exp2[1]=0;
			prodcomp(Exp, exp2, exp3 );
			G00[0]-=exp3[0]/kn/2/L;
			G00[1]-=exp3[1]/kn/2/L;
		}
		
	}
	
}
void greenfunction3(double x, double z, double k,  double* G0)
{
	double r, coe;
	double pi = 3.14159265358979323846; 
	r=sqrt(x*x+z*z);
	if (r==0){
		G0[0]=0;
		G0[1]=-0.25;
	}else{
	coe=sqrt(1.0/(pi*k*r*8));
	G0[0]=coe*sin(k*r-pi/4);
	G0[1]=coe*cos(k*r-pi/4);
	}
	
}

void greenfunctionV(double x, double z, double k,  double* G00)
{
	double pi = 3.14159265358979323846; 
	double Exp[2],exp2[2],exp3[2];
	double kz2,kn;
	int n,i;
	G00[0]=0;
	G00[1]=0;
	double L=2000;
	for (i=0;i<401;i++){
		n=i-200;
		kz2=k*k-(2*pi*n/L)*(2*pi*n/L);
		Exp[0]=cos(2*pi*n*x/L);
		Exp[1]=sin(2*pi*n*x/L);
		if (kz2>0){
			kn=sqrt(kz2);			
			exp2[0]=sin(kn*fabs(z));
			exp2[1]=0;
			prodcomp(Exp, exp2, exp3 );
			G00[0]+=exp3[0]/kn/L;
			G00[1]+=exp3[1]/kn/L;
		}
		
		
	}
	
}

void P0numinit( int iw, int nx, double dx, double k, double *sourcefren, int *ps,  double *P0r,double *P0i){
	double pi = 3.14159265358979323846; 
	int ix,iz;
   	double G0[2],sw[2],po[2];
  	double x,z;
	sw[0]=sourcefren[2*iw];
	sw[1]=sourcefren[2*iw+1];
	z=0;	
  	for(ix=0;ix<nx;ix++){
		x=(ix-ps[0])*dx;		
		greenfunction3(x, z, k, G0);
		prodcomp( sw, G0, po );
		P0r[ix]=po[0];
		P0i[ix]=po[1];	
		
	}
	FILE* file2;  
  	char fname2[100];
	if((iw%50)==0){
		sprintf(fname2,"P0_iw_%d.dat", iw);
		file2 = fopen(fname2,"w");
 		for (ix=0;ix<nx;ix++){
			x=(ix-(nx-1)/2)*dx;
			fprintf(file2," %f %.12lf %.12lf \n", x, P0r[ix],P0i[ix]);			                    
		}

       		 fclose(file2);
	}
}










void Greentable( int iw, int nx2, double dx,double z, double k,    double *G0r,double *G0i){
	int ix,nx=(nx2-1)/2+1;
	double x;
	double G0[2];
	for(ix=0;ix<nx2;ix++){
		x=(ix-nx+1)*dx;
		greenfunction2( x, z, k,  G0);
		//greenfunctionV( x, z, k,  G0);
		G0r[ix]=G0[0];
		G0i[ix]=G0[1];
	}
		

}

void P1num(FILE *filefprintf, int nw, int nx, int nz,double dx,double dz,double c0, int *ps, double *sourcefren,double* fren,double* vpe,  double *P1r,double *P1i){
	double pi = 3.14159265358979323846; 
	int iw,ix,iz,ix2,iz2;
        int xindex;
	double omega,k;
  	double x,z;
 	int nk=nx;
	double dk=2*pi/2/(nx-1)/dx;
	double P1tempr,P1tempi,G0[2];
	fprintf(filefprintf, "--!\tP1 calculation  ...\n");
  	fprintf(filefprintf, "--!\t   		\n");

        double* P0r = malloc(nz*nx*sizeof(double) );
   	double* P0i = malloc(nz*nx*sizeof(double) );

	int nx2=(nx-1)*2+1;
	double* G0r = malloc(nx2*sizeof(double) );
   	double* G0i = malloc(nx2*sizeof(double) );

	double cmin=c0,vmin=0,kne;
	FILE* file2;  
  	char fname2[100];
  	
        for(iw=0;iw<nw;iw++){
		fprintf(filefprintf, "--!\tP1 calculation  ...%d\n",iw);
  		fprintf(filefprintf, "--!\t   		\n");

		omega=fren[iw];
                k=omega/c0;
		P0numinit(  iw, nx, dx,  k, sourcefren, ps,  P0r,P0i);
		
		for(ix=0;ix<nx;ix++){			
                	P1r[iw*(nx*nz)+ix]=P0r[ix];
			P1i[iw*(nx*nz)+ix]=P0i[ix];
		}
		
		for(iz=1;iz<nz;iz++){
			vmin=vpe[iz-1];
			for(ix=1;ix<nx;ix++){
				if(vpe[ix*nz+iz-1]<vmin){
					vmin=vpe[ix*nz+iz-1];
				}				
			}
			cmin=c0/sqrt(1-vmin);
			kne=omega/cmin;		
			Pincforward( iw, nx,  dx, nk, dk,  kne,  iz, dz, P0r,P0i);
			Greentable( iw, nx2, dx, dz,  kne,  G0r,G0i);
			for(ix=0;ix<nx;ix++){			
                		P1r[iw*(nx*nz)+iz*nx+ix]=P0r[iz*nx+ix];
				P1i[iw*(nx*nz)+iz*nx+ix]=P0i[iz*nx+ix];
				P1tempr=0;
				P1tempi=0;
				for(ix2=0;ix2<nx;ix2++){
					xindex=(ix-ix2)+(nx-1);				
					G0[0]=	G0r[xindex];	
					G0[1]=	G0i[xindex];			
					P1tempr+=(vpe[ix2*nz+iz-1]-vmin)*(G0[0]*P1r[iw*(nx*nz)+(iz-1)*nx+ix]-G0[1]*P1i[iw*(nx*nz)+(iz-1)*nx+ix]);
					P1tempi+=(vpe[ix2*nz+iz-1]-vmin)*(G0[0]*P1i[iw*(nx*nz)+(iz-1)*nx+ix]+G0[1]*P1r[iw*(nx*nz)+(iz-1)*nx+ix]);
				}
				P1tempr*=dx*dz*k*k;
				P1tempi*=dx*dz*k*k;	
				P1r[iw*(nx*nz)+iz*nx+ix]+=P1tempr;
				P1i[iw*(nx*nz)+iz*nx+ix]+=P1tempi;
				P0r[iz*nx+ix]=P1r[iw*(nx*nz)+iz*nx+ix];	
				P0i[iz*nx+ix]=P1i[iw*(nx*nz)+iz*nx+ix];
			}
				  			
		}

		if((iw%50)==0){
		sprintf(fname2,"P1_iw_%d.dat", iw);
		file2 = fopen(fname2,"w");
		for (iz=0;iz<nz;iz++){
			z=iz*dz;
 			for (ix=0;ix<nx;ix++){
				x=(ix-(nx-1)/2)*dx;					
                          	fprintf(file2," %f %f %.12lf %.12lf \n", x, z, P1r[iw*(nx*nz)+iz*nx+ix],P1i[iw*(nx*nz)+iz*nx+ix]);
			}
                          fprintf(file2," \n");
		}

       		 fclose(file2);
		}
		
	}
	
}



void Pwtot(int nw, int nx,  int nz,  double dw, double dx,double dz, double dt, int nt,  double* fren, double *P1r, double *P1i, double *Pret){
	int ix,iw,it,iz;
	double pi = 3.14159265358979323846; 
	double t,x,z, arg=0;
        for(iz=0;iz<nz;iz++){
	for(ix=0;ix<nx;ix++){
		x=(ix-(nx-1)/2)*dx;
		for(it=0;it<nt;it++){
			t=it*dt;
			Pret[iz*nt*nx+ix*nt+it]=0;
			for (iw=0;iw<nw;iw++){
				arg = fren[iw]*t;		
				Pret[iz*nt*nx+ix*nt+it]+=cos(arg)*P1r[iw*(nx*nz)+iz*nx+ix]-sin(arg)*P1i[iw*(nx*nz)+iz*nx+ix];
			}
			Pret[iz*nt*nx+ix*nt+it]*=dw*2/2/pi;
			
		}
		
		
	}
	}
	FILE* file2;  
  	char fname2[100];
	for(it=0;it<nt;it++){
		if((it%100)==0){
		sprintf(fname2,"snapshot_%d.dat", it);
		file2 = fopen(fname2,"w");
 		for (ix=0;ix<nx;ix++){
			x=(ix-(nx-1)/2)*dx;
			for (iz=0;iz<nz;iz++){
				z=iz*dz;		
                          	fprintf(file2," %f %f %.12lf \n", x, z, Pret[iz*nt*nx+ix*nt+it]);
			}
                          fprintf(file2," \n");
		}

       		 fclose(file2);
		}
	}

	for(iz=0;iz<nz;iz++){
		if((iz%20)==0){
		sprintf(fname2,"Trace_%d.dat", iz);
		file2 = fopen(fname2,"w");
 		for (ix=0;ix<nx;ix++){
			x=(ix-(nx-1)/2)*dx;
			for (it=0;it<nt;it++){
				t=it*dt;		
                          	fprintf(file2," %f %f %.12lf \n", x, t, Pret[iz*nt*nx+ix*nt+it]);
			}
                          fprintf(file2," \n");
		}

       		 fclose(file2);
		}
	}
	
}
