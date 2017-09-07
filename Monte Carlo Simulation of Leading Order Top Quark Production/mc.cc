#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <vector>
#include <fstream>
#include <ctime>
#include "qqmatrix.hh"
#include "histogram.hh"
#include "gluonmatrix.hh"
#include "LHAPDF/LHAPDF.h"
#include <libconfig.h++> //for config file
#define PI 3.14159265359

using namespace std;
using namespace libconfig;

double pdf( double x, double Q, int flav)
{  
    double xgCV = LHAPDF::xfx(x,Q,flav);      
    return xgCV;
}

double P3(double x1,double x2,double theta, double P, double m)
{
	double p3= (2.0*P*x1*x2*cos(theta)*(x1-x2)  +  (x1+x2)*sqrt(4.0*pow(P*x1*x2,2)+pow(m*(x1-x2)*cos(theta),2)-pow(m*(x1+x2),2)))/(pow((x1+x2),2)-pow(cos(theta)*(x1-x2),2));
	return p3;
}

double E3(double p3, double m)
{
	double e3= sqrt(p3*p3+m*m);
	return e3;
}

double E4(double x1, double x2, double e3, double P)
{
	double e4= (x1+x2)*P-e3;
	return e4;
}

double P4(double E4, double m)
{
	double p4=sqrt(E4*E4-m*m);
	return p4;
}

//Cross Section
double Dsigma(double x1,double x2,double theta,double j3,double j4,double p3,double E3,double E4, double p4, double P, double m, int type)
{

	//momenta
	double p1=x1*P;
	double p2=x2*P;

	//phase space dt
	double dtcos=2.0*p1*p3*(1.0+((p1-p2)*(1.0/((p1*p3)/E3+(p2*p3)/E3-p1*cos(theta)+p2*cos(theta)))*(cos(theta)-p3/E3)));
	double dt2= abs(dtcos*(1.0/(256.0*PI*pow(p1*p2,2)))*sin(theta));
	
	//coupling constant
	double Q= m;	
	double g = sqrt(4.0*PI*LHAPDF::alphasPDF(Q)); 
	
	//gluon contribution
	double dsigmagg= Ggmatrix(p1, p2, g, theta, j3, j4, p3, E3, E4, p4,m)*pdf(x1,Q,0)*pdf(x2,Q,0);
	
	//quark contribution
	double PDFsum=0;
	for(int i=1; i<6; i++)
	{
		PDFsum+= pdf(x1,Q,i)*pdf(x2,Q,-type*i) + pdf(x1,Q,-i)*pdf(x2,Q,i*type);
	}
	double dsigmaqq= Qqmatrix(p1, p2, g, theta, j3, j4, p3, E3, E4, p4,m)*PDFsum;
	
	//total cross section	
	double dsigma= (dsigmagg+dsigmaqq)*dt2/(x1*x2);
	return dsigma;
}


int main ()
{
	
	//Read variables from param.cfg
	Config cfg;
	cfg.readFile("param.cfg");
	int N = cfg.lookup("N");
	int P = cfg.lookup("BeamEnergy");
	string pdfname = cfg.lookup("PDFName");
	int pdfset = cfg.lookup("PDFSet");
	int pol = cfg.lookup("Pol");
	int type = cfg.lookup("ColliderType"); 
	double m = cfg.lookup("TMass");
	
	// Loading the PDF set
	LHAPDF::initPDFSet(pdfname,0);
	
	//Random number intialise and seed
    double seed = time(0);
    const gsl_rng_type * T; 
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,seed);
    
    //output file
    ofstream fout;
    fout.open("results.txt");
    
    //set variables
    double xmin=LHAPDF::getXmin(0);
    double scale=PI*3.8937966*pow(10,8); //convert to barns and * integration volume
    double vmin [3] = {xmin,xmin,0}; //min for x1,x2,theta
	double vmax [3] = {1.0,1.0,PI}; //max for x1,x2,theta
	double f=0; //sum of function dsigma
	double fsq=0; //sum squared, for error
	int Ncount=0; //successful collisions counter
	double v [10]={0,0,0,-pol,pol,0,0,0,0,0}; //variables x1,x2,theta,j3,j4,p3,E3,E4,p4,dsigma
	 
	//Monte Carlo
    for (int i=0; i<N; i++) 
    {
    	//assigning random numbers to x1,x2,theta
        for (int j=0; j<3; j++)
        {
            double u = gsl_rng_uniform (r); 
            v [j] = u*(vmax[j]-vmin[j])+vmin[j];    
        }
        
        //flip spins
        if(i%2==0)
        {
         	v[3]*=-1;
        }
        else
        {
           	v[4]*=-1;
        }    

		//kinematics
		v[5]=P3(v[0],v[1],v[2], P, m);
		if(v[5]<0)
		{
			v[5]=-v[5];
			v[2]=PI-v[2];
		}
		if(v[5]==v[5]) //check if kinematically allowed
		{
		    v[6]=E3(v[5],m);
		    v[7]=E4(v[0],v[1],v[6],P);
		    v[8]=P4(v[7],m);
		    v[9]=scale*Dsigma(v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],P,m,type); // cross section weight
		    		    
			//output results to file
			for(int j=0; j<10; j++) 
			{
				fout << v[j] << "  ";
			}   
			fout <<endl;
			
			//sum total cross section and error
			f +=v[9];
			fsq +=pow(v[9],2);
			Ncount ++;
		}
		else
		{
			v[9]=0;
		}
		
	}
	//total cross section and error
	double sigmat = f/N;
	double error = sqrt(((fsq/N)-pow(sigmat,2))/N);
    Histogram(Ncount, P, sigmat);
    cout << "The total cross-section is: " << sigmat << "pb, and the error is: "  << error << "pb" << endl;
 	fout << "The total cross-section is: " << sigmat << "pb, and the error is: "  << error << "pb" << endl;
    fout.close();
    gsl_rng_free (r);
    return 0;
}


