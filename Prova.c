#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "nrutil.h"
#include "nrutil.c"
#include "utilpersSPEC.h"
#include "utilpersSPEC.c"

#define pi 3.1415926535897932384626

#define T ((double) 100)
#define KKx ((int) 128) 
#define Kx ((int)  floor(KKx*2/3))

#define N     ((int) 1000)  //numero di step temporali 
#define q     ((int) N/1000)
#define L     ((double) 2*pi) // dimensione del dominio spaziale [0,L]

void Product(double *,double *,double *, int, int);

int main(void)
{
	int i, ii;
		double  *u, *v, *w; 
		double *ufis, *vfis, *wfis;
		
		FILE *fu;
		FILE *fv;
		FILE *fw;
		FILE *fuh;
		FILE *fvh;
		FILE *fwh;
		
	fu=fopen("aaU.txt","w");
	fv=fopen("aaV.txt","w");
	fw=fopen("aaW.txt","w");
	fuh=fopen("aaUmodes.txt","w");
	fvh=fopen("aaVmodes.txt","w");
	fwh=fopen("aaWmodes.txt","w");
	
	u=dvector(1,KKx);
	v=dvector(1,KKx);
	w=dvector(1,KKx);
		
	ufis=dvector(1,KKx);
	vfis=dvector(1,KKx);
	wfis=dvector(1,KKx);
	
	for(i=1;i<=KKx;i++)
	{ufis[i]=cos((i-1)*L/(KKx));
	//ufis[i]=((i-1)*L/(KKx-1))*((i-1)*L/(KKx-1));
//	vfis[i]=2*(i-1)*L/(KKx-1);
	}
	
	printf("%d \n", i);
	
		
	drealft(ufis,KKx,1);
//	drealft(vfis,KKx,1);


	
for(i=1;i<=KKx;i++)
	{u[i]=2*ufis[i]/(KKx);
//	v[i]=2*vfis[i]/(KKx);
	}
/*
//Calcolo la primitiva
for(i=1;i<=KKx/2;i++)
	{
		v[2*i-1]=-(i)*u[2*i]*pi/L;
		v[2*i]=(i)*u[2*i-1]*pi/L;
		
	}
*/
//Calcolo la derivata di u

for(i=1;i<=KKx/2;i++)
	{

	v[2*i]=-(i)*u[2*i-1]/2;
	v[2*i-1]=(i)*u[2*i]/2;
		
	}


	for(i=1;i<=KKx;i++)
	{ufis[i]=u[i];
	vfis[i]=v[i];
	}
	drealft(ufis,KKx,-1);
	drealft(vfis,KKx,-1);

//Calcolo la derivata di v

for(i=1;i<=KKx/2;i++)
	{
	
		w[2*i-1]=(i)*v[2*i]/2;
		w[2*i]=-(i)*v[2*i-1]/2;
		
	}
for(i=1;i<=KKx;i++)
	{wfis[i]=w[i];
	}
	drealft(wfis,KKx,-1);


//Calcolo la derivata seconda di u
/*
for(i=1;i<=KKx;i++)
{
	w[i]=-u[i]*i*i*pi*pi/L/L;
}


for(i=1;i<=KKx;i++)
	{wfis[i]=w[i];
	}
	drealft(wfis,KKx,-1);
*/

//calcolo il prodotto tra u e v 
/*
Product(u,v,w,KKx,Kx);

for(i=1;i<=KKx;i++)
	{wfis[i]=w[i];
	}
	drealft(wfis,KKx,-1);
*/	
			
	for(ii=1;ii<=KKx;ii++)
	{fprintf(fu,"%g \n", ufis[ii]);
	fprintf(fv,"%g \n", vfis[ii]);
	fprintf(fw,"%g \n", wfis[ii]);
	fprintf(fuh,"%g \n", u[ii]);
	fprintf(fvh,"%g \n", v[ii]);
	fprintf(fwh,"%g \n", w[ii]);}







return 0;
}

void Product(double *f,double *g,double *fxg,int HHx, int Hx)
{
	int k;
	double   *ff, *gg;
	
	ff=dvector(1,HHx);
	gg=dvector(1,HHx);

	for(k=1;k<=HHx;k++)
	{ff[k]=f[k];
	gg[k]=g[k];}

	drealft(ff,HHx,-1);
	drealft(gg,HHx,-1);
	for(k=1;k<=HHx;k++)
	{fxg[k]=ff[k]*gg[k]*2.0/HHx;}

	drealft(fxg,HHx,1);
	for(k=Hx+1;k<=HHx;k++)
	{fxg[k]=0.0;}

	free_dvector(ff,1,HHx);
	free_dvector(gg,1,HHx);
}
