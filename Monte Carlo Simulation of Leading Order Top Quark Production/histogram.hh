#ifndef histogram_hh
#define histogram_hh
#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <fstream>
#include <vector>
#include <libconfig.h++>

using namespace libconfig;

int Histogram (int Ncount, double P, double sigmat)
{
	Config cfg;
	cfg.readFile("param.cfg");
	ofstream fout;
    fout.open("hist.txt"); //outputs data	
	ifstream input;
	input.open("results.txt");

	//Bin sizes from param.cfg for 4 histograms
	int nobins[3] = {cfg.lookup("Rapbins"),cfg.lookup("COMbins"),cfg.lookup("PTransbins")}; //no bins
	double min[3]= {cfg.lookup("Rapmin"),cfg.lookup("COMmin"),cfg.lookup("PTransmin")}; //min value
	double max[3] = {cfg.lookup("Rapmax"),cfg.lookup("COMmax"),cfg.lookup("PTransmax")}; //max value
	double binwidth[3];
	for(int i=0; i<3; i++)
	{
		binwidth[i]=(max[i]-min[i])/nobins[i]; //binwidth
	}	
	double hist[300][3] ={0}; //histogram array, max number of bins set to 300 
	double histerr[300][3]={0}; //fsq error for each bin
	int histcount[300][3] ={0}; //counts results in each bin
	double norm[3]={0}; //normalise
	
	//put values in bins
	for(int k=0; k<Ncount; k++)
	{		 		
		double x1, x2, theta, j3, j4, p3, E3, E4, p4, dsigma;
		input >> x1>> x2>> theta >> j3>> j4>> p3>> E3>> E4>> p4>> dsigma; //read data

		//functions to plot
		double funct[3]; 
		funct[0] = 0.5*log((E4+((x1-x2)*P-p3*cos(theta)))/(E4-((x1-x2)*P-p3*cos(theta))));//rapidity
		funct[1] = sqrt(4*x1*x2*P*P); //com energy
	    funct[2] = p3*sin(theta); //transverse momentum	
		int bin[3]; //to find which bin

		for(int i=0; i<3; i++)
		{			
			bin[i] = (funct[i]-min[i])/binwidth[i];

			if(funct[i]>=max[i]) //upper overflow in last bin
			{
			}
			else if(funct[i]<=min[i]) //lower overflow in first bin
			{
			}
			else
			{			
				hist[(bin[i])][i] += dsigma;
				histcount[(bin[i])][i] ++;
				histerr[(bin[i])][i] += pow(dsigma,2);
			}
		}
	}
	
	for(int i=0; i<300; i++)
	{
		for(int j=0; j<3; j++)
		{
			if(histcount[i][j]!=0)
			{
				norm[j]+=hist[i][j]/histcount[i][j];	
			}
		}
	}

	//output data
	for(int i=0; i<300; i++)
	{		
		for(int j=0; j<3; j++)
		{
			double histsize=hist[i][j]/histcount[i][j];
			fout << (i+0.5)*binwidth[j]+min[j] << "  " << sigmat*histsize/norm[j] << "  " << sigmat*sqrt((histerr[i][j]/histcount[i][j]-pow(histsize,2))/histcount[i][j])/norm[j] << "  " ; //centre of bin value, weighted cross-section, error in bin
		}
		fout <<endl;
	}

	fout.close();
	input.close();
	
    return 0;
}
#endif






                  







