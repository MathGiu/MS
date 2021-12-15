#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//The following routines can be found in "Numerical recipes: The art of scientific computing, Cambridge university press, 2007"
#include "nr.h"
#include "nrutil.h"
#include "nrutil.c"
#include "utilpersSPEC.h"
#include "utilpersSPEC.c"

#define pi 3.1415926535897932384626

//Functions used in the simulation
void IC(double *, double *); //User needs to write the initial condition
void Product(double *,double *,double *); //This function computes products in the physical space
void Flux(double *, double *, double *, double *); // This function computes the r.h.s. of the equations for u and v


#define T ((double) ) //User needs to insert the total time of the simulation
#define N ((int) ) //User needs to insert the number of temporal simulation steps
#define	dt ((double) T/N) //Length of the temporal grid

#define L ((double) ) //User needs to insert the length of the spatial domain
#define m ((int) ) //User needs to insert the number of spatial Fourier modes (must be a power of 2)
#define md ((int)  floor(m*2/3)) // This value is useful for the de-aliasing
#define dx ((double) L/m) //Length of the spatial grid


//User needs to insert parameter values
#define D1 ((double) )
#define D2  ((double) )
#define gamma12 ((double) )
#define gamma21 ((double) )
#define gamma11 ((double) )
#define gamma22 ((double) )

//User needs to insert Top-hat distribution parameter value
//The Top-hat distribution is f(x)=1/(2*alpha), for -alpha<x<alpha, f(x)=0, otherwise
#define alpha ((double) )
#define h ((int) floor(alpha/dx)) //Number of spatial grid points in a interval of length alpha

//User needs to insert Von-mises distribution parameter values
//The Von-mises distribution is f(x)=exp(a*cos(2*pi*x/L))/(B*L)
#define a ((double) )
#define B ((double) )

//User needs to choose the Kernel
#define c ((int) ) //Choose c=0 for Top-hat distribution, c=1 for Von-mises distribution

int main(void)
{
	
	double *u, *v; //Solution vectors in the Fourier space
	double *up, *vp; //Solution vectors in the physical space
	double *hu1,*hv1,*hu2,*hv2; //Vectors used for Runge-Kutta scheme
	double *fu1,*fv1,*fu2,*fv2;	//Vectors used for Runge-Kutta scheme
	int n,i,k;

	
	//Open the files in which the simulation will be saved
	FILE *Fu;
	FILE *Fv;
	FILE *Fup;
	FILE *Fvp;	
	
	Fu=fopen("U.txt","w");
	Fv=fopen("V.txt","w");
	Fup=fopen("UP.txt","w");
	Fvp=fopen("VP.txt","w");


	//Array Initialization
	u=dvector(1,m);
	v=dvector(1,m);
	up=dvector(1,m);
	vp=dvector(1,m);
		
    hu1=dvector(1,m);
    hv1=dvector(1,m);
    hu2=dvector(1,m);
    hv2=dvector(1,m);    

	fu1=dvector(1,m);
	fv1=dvector(1,m);
	fu2=dvector(1,m);
	fv2=dvector(1,m);

	//Set up the initial conditions
	IC(up,vp);	

	//Transform the data from the physical space to the Fourier space	
	for(i=1;i<=m;i++)
	{	
		u[i]=up[i]*2/m;
		v[i]=vp[i]*2/m;
	}
	
	drealft(u,m,1);
	drealft(v,m,1);
	
	for(k=md+1;k<=m;k++)
	{
		u[k]=0.0;
		v[k]=0.0;
	}
	
	//Save the initial condition 
	for(i=1;i<=m;i++)
	{
		fprintf(Fup,"%g \n", up[i]);
		fprintf(Fvp,"%g \n", vp[i]);
		fprintf(Fu,"%g \n", u[i]);
		fprintf(Fv,"%g \n", v[i]);
	}

	//Second order Runge-Kutta scheme
	for(n=1;n<=N;n++)
	{
		Flux(u, v, fu1, fv1);
				
		for(k=1;k<=m;k++)
        {
        	hu1[k]=u[k]+fu1[k]*dt;
        	hv1[k]=v[k]+fv1[k]*dt;        		
		}
		
		Flux(hu1, hv1, fu2, fv2);
		
		for(k=1;k<=m;k++)
		{
			u[k]=u[k]+dt*(fu1[k]+fu2[k])/2;
			v[k]=v[k]+dt*(fv1[k]+fv2[k])/2;		
		}
		
		//Write the solution at time t= n*dt in the physical space		
		for(i=1;i<=m;i++)
		{	
			up[i]=u[i];
			vp[i]=v[i];
		}
	
		drealft(up,m,-1);
		drealft(vp,m,-1);

		//Save the solution at time t= n*dt 
		for(i=1;i<=m;i++)
		{
				
			fprintf(Fup,"%g \n", up[i]);
			fprintf(Fvp,"%g \n", vp[i]);
			fprintf(Fu,"%g \n", u[i]);
			fprintf(Fv,"%g \n", v[i]);
		}
		
	} 
	
	//Close the file
	fclose(Fu);	
	fclose(Fv);	
	fclose(Fup);	
	fclose(Fvp);
	
	return 0;
}

//User needs to write the initial condition 
void IC(double *UI, double *VI)
{
	int i;

	for(i=1;i<=m;i++)
	{
		UI[i]= ;
		VI[i]= ;		
	}
	
}

//This function transform two input data from the Fourier space to the physical space,
//computes their products and finally transforms the result from the physical space
//to the Fourier space
void Product(double *f, double *g, double *fxg)
{
	int k;
	double *ff, *gg;
	
	//Array Initialization	
	ff=dvector(1,m);
	gg=dvector(1,m);
	
	//Tranform the solution from the Fourier space to the physical space
	for(k=1;k<=m;k++)
	{ff[k]=f[k];
	gg[k]=g[k];}

	drealft(ff,m,-1);
	drealft(gg,m,-1);
	
	//Compute the product and return to the Fourier space
	for(k=1;k<=m;k++)
	{fxg[k]=ff[k]*gg[k]*2.0/m;}

	drealft(fxg,m,1);
	for(k=md+1;k<=m;k++)
	{fxg[k]=0.0;}

	free_dvector(ff,1,m);
	free_dvector(gg,1,m);
}


//This function computes the r.h.s. of the equations for u and v
void Flux(double *U, double *V, double *ReactU, double *ReactV)
{
	int i,j,k;
	double *K, *Kp; //Kernel, in the Fourier and in the physical space, respectively
	double *dU, *dV, *ddU, *ddV; //First and second spatial derivatives of u and v
	double *u_bar, *v_bar; //Convolution products
	double *du_bar, *vdu_bar, *dvdu_bar; //grad(u_bar), v*grad(u_bar), div(v*grad(u_bar))
	double *dv_bar, *udv_bar, *dudv_bar; //grad(v_bar), u*grad(v_bar), div(u*grad(v_bar))
	double *udu_bar, *dudu_bar,*vdv_bar, *dvdv_bar; //u*grad(u_bar), div(u*grad(u_bar)), v*grad(v_bar), div(v*grad(v_bar))
	
	//Array initialization
	K=dvector(1,m);
	Kp=dvector(1,m);
	dU=dvector(1,m);
	dV=dvector(1,m);
	ddU=dvector(1,m);
	ddV=dvector(1,m);	
	u_bar=dvector(1,m);
	v_bar=dvector(1,m);
	du_bar=dvector(1,m);
	vdu_bar=dvector(1,m);
	dvdu_bar=dvector(1,m);
	dv_bar=dvector(1,m);
	udv_bar=dvector(1,m);
	dudv_bar=dvector(1,m);	
	udu_bar=dvector(1,m);
	dudu_bar=dvector(1,m);
	vdv_bar=dvector(1,m);
	dvdv_bar=dvector(1,m);

	//Write the Kernel in the Physical space		
	if(c==0)
	{	//Top-hat distribution
	
		for(i=1;i<=m;i++)
		{
			Kp[i]=0.0;
		}

		for(i=1;i<=2*h;i++)
		{		
			Kp[m/2-h+i]=1/(2*h*dx);	
		}
	}
	
	else
	{	//VonMises distribution
		for(i=1;i<=m;i++)
		{	
			Kp[i]=exp(a*cos(2*pi*(i-1)*dx/L))/(B*L);	
		}
	}

	//Write the kernel in the Fourier space
		
	for(j=1;j<=m;j++)
		{
			K[j]=2.0*Kp[j]/m;
	 	}
	 	
	drealft(K,m,1);
	 
	for(j=md+1;j<=m;j++)
		{
			K[j]=0.0;
		}

	//Compute the linear diffusion terms
	for(i=1;i<=m/2;i++)
	{
		dU[2*i]=-i*U[2*i-1];
		dU[2*i-1]=i*U[2*i];
		dV[2*i]=-i*V[2*i-1];
		dV[2*i-1]=i*V[2*i];		
	}
	
	for(i=1;i<=m/2;i++)
	{
		ddU[2*i]=-i*dU[2*i-1];
		ddU[2*i-1]=i*dU[2*i];
		ddV[2*i]=-i*dV[2*i-1];
		ddV[2*i-1]=i*dV[2*i];		
	}
	
	//Compute the convolution products K*u and K*v in the Fourier space
	for(j=1;j<=m;j++)
		{
			u_bar[j]=K[j]*U[j];
			v_bar[j]=K[j]*V[j];
	 	}

	//Compute the terms grad(K*u) and grad(K*v)
	for(i=1;i<=m/2;i++)
	{
		du_bar[2*i]=-i*u_bar[2*i-1];
		du_bar[2*i-1]=i*u_bar[2*i];	
		dv_bar[2*i]=-i*v_bar[2*i-1];
		dv_bar[2*i-1]=i*v_bar[2*i];
	}

	//Compute the products u*(grad(K*v)), v*(grad(K*u)), etc...
	Product(U,dv_bar,udv_bar);
	Product(V,du_bar,vdu_bar);
	Product(U,du_bar,udu_bar);
	Product(V,dv_bar,vdv_bar);

	
	//Compute the terms div(u*(grad(k*v))), div(v*(grad(k*u))), etc...	
	for(i=1;i<=m/2;i++)
	{
		dudv_bar[2*i]=-i*udv_bar[2*i-1];
		dudv_bar[2*i-1]=i*udv_bar[2*i];		
		dvdu_bar[2*i]=-i*vdu_bar[2*i-1];
		dvdu_bar[2*i-1]=i*vdu_bar[2*i];		
		dudu_bar[2*i]=-i*udu_bar[2*i-1];
		dudu_bar[2*i-1]=i*udu_bar[2*i];		
		dvdv_bar[2*i]=-i*vdv_bar[2*i-1];
		dvdv_bar[2*i-1]=i*vdv_bar[2*i];
	}
	

	
			
		
	//Write the r.h.s. of the equations for u and v
	for(k=1;k<=m;k++)
	{
		ReactU[k]=D1*ddU[k]+gamma12*dudv_bar[k]+gamma11*dudu_bar[k];
		ReactV[k]=D2*ddV[k]+gamma21*dvdu_bar[k]+gamma22*dvdv_bar[k];
	}

	
	free_dvector(K,1,m);
	free_dvector(Kp,1,m);		
	free_dvector(dU,1,m);
	free_dvector(ddU,1,m);
	free_dvector(dV,1,m);
	free_dvector(ddV,1,m);	
	free_dvector(u_bar,1,m);
	free_dvector(v_bar,1,m);		
	free_dvector(du_bar,1,m);
	free_dvector(dv_bar,1,m);
	free_dvector(vdu_bar,1,m);
	free_dvector(udv_bar,1,m);
	free_dvector(dvdu_bar,1,m);
	free_dvector(dudv_bar,1,m);
	free_dvector(udu_bar,1,m);
	free_dvector(dudu_bar,1,m);
	free_dvector(vdv_bar,1,m);
	free_dvector(dvdv_bar,1,m);

}



	


