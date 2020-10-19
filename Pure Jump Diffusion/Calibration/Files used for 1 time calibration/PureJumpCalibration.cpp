#include "COS.h"

#include <cstdlib>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <fstream> 
#include <cmath>
#include <complex>
#include <string>
#include <list>
#include <stdlib.h> 
#include <vector>

#define pi 3.14159265358979
using namespace std;

//the meaning of the variables
/*
double K;                       //the strike price of the option
double S_0;                     //the current stock price
double c1;                      //the parameter for calculating a and b
double c2;                      //the parameter for calculating a and b
double r;                       //the risk free rate
double q;                       //the dividend yield
double T;                       //the maturity of the option
double N;                       //the integration discretization constant
double sigma_s;					//the volatility of stock
double freq;					//lamda of poisson process
double mu_j;					//mean of jump size
double var_j;					//var of jump size
int C_P;                        //"1" for Call and "0" for Put
*/
int factorial(int n) {
	if(n>1)
		return n * factorial(n-1);
	else;
		return 1;
}

double norm_s_cdf(double x) {
	// constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0) {
		sign = -1;
    	x = fabs(x)/sqrt(2.0);
	}
       
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

double m_n(double r, double sigma_s, double freq, double mu_j, double var_j, double T, int N) {
	return (r-0.5*sigma_s*sigma_s-freq*(exp(mu_j+0.5*var_j)-1)+(N*mu_j)/T);
}

double d1(double S_0,double K,double m_n,double T,double sigma_s,double var_j,int N) {
	return ((log(S_0/K)+m_n*T)/sqrt(sigma_s*sigma_s*T+N*var_j));
}

double d2(double d1,double sigma_s,double T,double var_j,int N) {
	return (d1-sqrt(sigma_s*sigma_s*T+N*var_j));
}

double Jump_option_price(double T,double sigma_s,double freq, double mu_j, double var_j, double r,double q,double S_0,double K,int N) {
	int k;
	double a,b,c,f,g;
	double sum = 0;
  
	//the option price can be found numerical integration
	for (k=0;k<N;k++) {
		a = m_n(r,sigma_s,freq,mu_j,var_j,T,k);
		b = d1(S_0,K,a,T,sigma_s,var_j,k);
		c = d2(b,sigma_s,T,var_j,k);
		f = norm_s_cdf(b);
		g = norm_s_cdf(c);
    	sum = sum + exp(-freq*T)*pow(freq*T,k)*(1/factorial(k))*(S_0*exp(a*T)*f-K*g);
	}

	//return answer
	return (exp(-r*T)*sum);
};

double RandomUniform(double minValue,double maxValue){
     return minValue + ((double)rand()/(double)RAND_MAX) * (maxValue - minValue);
};

int main()
{
	double min[4] = {0,0,-10,0};
	double max[4] = {1,1,10,1};
	int n = 4;
	int nPop = 75;
	int maxGenerations = 1500;
	double scale	   = 0.8;
	double probability = 0.5;
	char tmp[1000];
	char tmp1[1000];

	vector<double> rf;
	vector<double> ST;
	vector<double> K;
	vector<double> q;
	vector<double> T;
	vector<double> weight;
	vector<double> mktPrice;

	std::ifstream fp_in;

	fp_in.open("input.txt", ios::in);
	while (1)
	{
		if (fp_in.eof())break;
		fp_in >> tmp;
		rf.push_back(atof(tmp));
		fp_in >> tmp;
		ST.push_back(atof(tmp));
		fp_in >> tmp;
		K.push_back(atof(tmp));
		fp_in >> tmp;
		q.push_back(atof(tmp));
		fp_in >> tmp;
		T.push_back(atof(tmp));
		fp_in >> tmp;
		weight.push_back(atof(tmp));
	}
	fp_in.close();
	
	fp_in.open("MarketPrices.txt", ios::in);
	for (int i = 0; i<rf.size(); i++)
	{
	if (fp_in.eof())break;
		fp_in >> tmp1;
		mktPrice.push_back(atof(tmp1));
	}
	fp_in.close();

	ofstream fp_out;

	double bestSolution[n];
	double Element[nPop][n];
	double popEnergy[nPop];
	int i;
	/*
	Solution[0] = sigma_s
	Solution[0] = freq
	Solution[1] = mu_j
	Solution[2] = var_j
	*/

	for (i=0; i < nPop; i++)
	{
		for (int j=0; j < n; j++){
			Element[i][j] = RandomUniform(min[j],max[j]);
        }

		double t=0.0;
		double t_s=0.0;
      	for (int k = 0; k<rf.size(); k++){
			t = Jump_option_price(T[k],Element[i][0], Element[i][1], Element[i][2],Element[i][3],rf[k], q[k], ST[k], K[k],32) - mktPrice[k];
			t_s = t_s + t*t*weight[k];
		}
		
		popEnergy[i] = t_s;
	}

	for (i=0; i < n; i++){
	bestSolution[i] = 0.0;
	}
		
	int candidate;
	int generation;
	bool bAtSolution;
	double bestEnergy = 10000000000.0;
	bAtSolution = false;
	double trialSolution[n];
	double trialEnergy=0.0;

	for (generation = 0; (generation < maxGenerations) && !bAtSolution; generation++){
		cout  <<generation+1<<"\t"<<"\t"<<trialSolution[0]<<" "<<trialSolution[1]<<" "<<trialSolution[2]<<" "<<trialSolution[3]<<" "<<endl;
		cout  <<"\t"<<bestEnergy<<"\t"<<bestSolution[0]<<" "<<bestSolution[1]<<" "<<trialSolution[2]<<" "<<trialSolution[3]<<" "<<endl<<endl;
		for (candidate = 0; candidate < nPop; candidate++)
		{
			//strategy
			int r1, r2, r3;
			int nt;
			
			do{r1 = (int)RandomUniform(0.0,(double)nPop);
			}while (r1 == candidate);
			do{r2 = (int)RandomUniform(0.0,(double)nPop);
			}while ((r2 == candidate)||(r2 == r1) );
			do{r3 = (int)RandomUniform(0.0,(double)nPop);
			}while ((r3 == candidate)||(r3 == r2)||(r3 == r1) );
			
			/**/
			nt = (int)RandomUniform(0.0,(double)n);
				
			for (int i=0; i < n; i++) {
				trialSolution[i] = Element[candidate][i];
				}
				
			/*Rand1Bin*/
			for (int i=0; i < n; i++) {
				if ((RandomUniform(0.0,1.0) < probability) || (i == (n - 1)))
					trialSolution[nt] = Element[r1][nt] + scale * (Element[r2][nt] - Element[r3][nt]);
					if(trialSolution[nt]>max[nt])trialSolution[nt]=(Element[r1][nt]+max[nt])*0.5;
					if(trialSolution[nt]<min[nt])trialSolution[nt]=(Element[r1][nt]+min[nt])*0.5;
					nt = (nt + 1) % n;
					}
			
			double t = 0.0;
			double t_s = 0.0;
      		for (int i = 0; i<rf.size(); i++){
				t = Jump_option_price(T[i],trialSolution[0], trialSolution[1], trialSolution[2], trialSolution[3],rf[i], q[i], ST[i], K[i], 32) - mktPrice[i];
				t_s = t_s + t*t*weight[i];
			}
			trialEnergy = t_s;
			
			if (trialEnergy < popEnergy[candidate])
			{
				// New low for this candidate
				popEnergy[candidate] = trialEnergy;
				for (int i=0; i < n; i++) {
					Element[candidate][i]=trialSolution[i];
					}

				// Check if all-time low
				if (trialEnergy < bestEnergy)
				{
					bestEnergy = trialEnergy;
					for (int i=0; i < n; i++) {
						bestSolution[i]=trialSolution[i];
					};
				}
			}
			
		}
	}

	FILE *fp_para;
	
	fp_para = fopen("Parameters.txt", "w" );
	cout<<"\n\nBest Coefficients:\n"<<endl;
	 for (i=0;i<n;i++){
		fprintf(fp_para,"%.10lf\n", bestSolution[i]);
    }
	 fclose(fp_para);

	system("PAUSE");
	return 0;
}
