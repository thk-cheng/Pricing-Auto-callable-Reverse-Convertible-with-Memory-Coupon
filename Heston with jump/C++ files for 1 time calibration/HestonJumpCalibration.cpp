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
double lamda;                   //lamda is the speed of mean reversion 
double rho;                     //rho is the correlation coef.                         
double eta;                     //the volatility of the volatility
double mu;                      //the mean return of the assest
double nu_bar;                  //long term mean level of the volatility
double u_0;                     //the initial vol
double c1;                      //the parameter for calculating a and b
double c2;                      //the parameter for calculating a and b
double r;                       //the risk free rate
double q;                       //the dividend yield
double T;                       //the maturity of the option
double N;                       //the 
int C_P;                        //"1" for Call and "0" for Put
*/

//calculate the function D
complex<double> D (double omega,double lamda,double eta,double rho) {
  complex<double> c_i(0,1);
  return (sqrt(pow(lamda-c_i*rho*eta*omega,2)+(pow(omega,2)+c_i*omega)*pow(eta,2)));
};

//calculate the function G
complex<double> G (double omega,double lamda,double eta,double rho) {
  complex<double> a,b,c_D;
  complex<double> c_i(0,1);
  c_D = D(omega,lamda,eta,rho);
  a = lamda-c_i*rho*eta*omega-c_D;
  b = lamda-c_i*rho*eta*omega+c_D;
  return a/b;
};

//calculate teh characteristic function of the log-asset price under the Heston model
complex<double> phi_hes_j (double r,double q,double omega,double lamda,double nu_bar,double u_0,double eta,double rho,double T,double freq, double mu_j, double var_j) {

  complex<double> a,b,c,d,e,f;
  complex<double> c_D, c_G;
  complex<double> c_i(0,1);
  c_D = D(omega,lamda,eta,rho);
  c_G = G(omega,lamda,eta,rho);
    
  a = c_i*omega*(r-q)*T;
  b = u_0/(eta*eta)*((1.0-exp(-c_D*T))/(1.0-c_G*exp(-c_D*T)))*(lamda-c_i*rho*eta*omega-c_D);
  c = lamda*nu_bar/(eta*eta);
  d = T*(lamda-c_i*rho*eta*omega-c_D)-2.0*log((1.0-c_G*exp(-c_D*T))/(1.0-c_G));
  e = c_i*omega*freq*T*(exp(mu_j+0.5*var_j)-1);
  f = freq*T*(exp(c_i*mu_j*omega-0.5*omega*omega*var_j)-1.0);
  
  return (exp(a+b)*exp(c*d)*exp(-e+f));
};

//calculate the cosine series coefficients of g(y)=1
double psi(double c,double d,double a,double b,double k){
  if(k==0)
    return (d-c);
  else if (k!=0)
    return ((sin(k*pi*(d-a)/(b-a))-sin(k*pi*(c-a)/(b-a)))*(b-a)/(k*pi));
};

//calculate the cosine series coefficients of g(y)=exp(y)
double chi(double c,double d,double a,double b,double k){
  return ((cos(k*pi*(d-a)/(b-a))*exp(d)-cos(k*pi*(c-a)/(b-a))*exp(c)+sin(k*pi*(d-a)/(b-a))*exp(d)*k*pi/(b-a)-sin(k*pi*(c-a)/(b-a))*exp(c)*k*pi/(b-a))/(1+(k*pi/(b-a))*(k*pi/(b-a))));
};

//calculate the function of Uk
double U(double k,int C_P,double a,double b){
  if(C_P==1)      //for call
    return 2*(chi(0,b,a,b,k)-psi(0,b,a,b,k))/(b-a);
  else if(C_P==0)    //for put
    return 2*(-chi(a,0,a,b,k)+psi(a,0,a,b,k))/(b-a);
};

//calculate the value of a
double a(double r,double q,double T,double lamda,double nu_bar,double u_0,double eta,double rho,double freq,double mu_j,double var_j) {
  double c1,c2,c3,c4;
  c1 = ((r-q)*T+(1-exp(-lamda*T))*((nu_bar-u_0)/(2*lamda))-0.5*nu_bar*T);
  c2 = ((1/(8*lamda*lamda*lamda))*((eta*T*lamda*exp(-lamda*T)*(u_0-nu_bar)*(8*lamda*rho-4*eta))+(lamda*rho*eta*(1-exp(-lamda*T))*(16*nu_bar-8*u_0))+(2*nu_bar*lamda*T*(-4*lamda*rho*eta+eta*eta+4*lamda*lamda))+(eta*eta*((nu_bar-2*u_0)*exp(-2*lamda*T)+nu_bar*(6*exp(-lamda*T)-7)+2*u_0))+(8*lamda*lamda*(u_0-nu_bar)*(1-exp(-lamda*T)))));
  c3 = freq*T*(mu_j-(exp(mu_j+0.5*var_j)-1));
  c4 = freq*T*(mu_j*mu_j+var_j);
  return (c1+c3-10*sqrt(abs(c2+c4)));   
};

//calculate the value of b
double b(double r,double q,double T,double lamda,double nu_bar,double u_0,double eta,double rho,double freq,double mu_j,double var_j) {
  double c1,c2,c3,c4;
  c1 = ((r-q)*T+(1-exp(-lamda*T))*((nu_bar-u_0)/(2*lamda))-0.5*nu_bar*T);
  c2 = ((1/(8*lamda*lamda*lamda))*((eta*T*lamda*exp(-lamda*T)*(u_0-nu_bar)*(8*lamda*rho-4*eta))+(lamda*rho*eta*(1-exp(-lamda*T))*(16*nu_bar-8*u_0))+(2*nu_bar*lamda*T*(-4*lamda*rho*eta+eta*eta+4*lamda*lamda))+(eta*eta*((nu_bar-2*u_0)*exp(-2*lamda*T)+nu_bar*(6*exp(-lamda*T)-7)+2*u_0))+(8*lamda*lamda*(u_0-nu_bar)*(1-exp(-lamda*T)))));
  c3 = freq*T*(mu_j-(exp(mu_j+0.5*var_j)-1));
  c4 = freq*T*(mu_j*mu_j+var_j);
  return (c1+c3+10*sqrt(abs(c2+c4)));    
};

double COS(double T,double lamda,double nu_bar,double u_0,double eta,double rho,double freq, double mu_j, double var_j, double r,double q,double S_0,double K,int C_P,int N) {
  int k;
  double ct1,ct2;
  double x;
  ct1 = a(r,q,T,lamda,nu_bar,u_0,eta,rho,freq,mu_j,var_j);         //ct1 is the lower bound 
  ct2 = b(r,q,T,lamda,nu_bar,u_0,eta,rho,freq,mu_j,var_j);         //ct2 is the upper bound
  x = log(S_0/K);
  complex<double> c_i(0,1);
  complex<double> sum (0,0);
  
  //given the log-asseat price, the option price can be found by IFFT
  sum = 0.5*real(phi_hes_j(r, q, 0,lamda,nu_bar,u_0,eta,rho,T,freq,mu_j,var_j)*U(0,C_P,ct1,ct2)*exp(0.0));
  for (k=1;k<N;k++) {
      complex<double> count (k,0);
      sum=sum+real(phi_hes_j(r, q, k*pi/(ct2-ct1),lamda,nu_bar,u_0,eta,rho,T,freq,mu_j,var_j)*U(k,C_P,ct1,ct2)*exp(c_i*count*pi*(x-ct1)/(ct2-ct1)));
  }
  //return th real part of the answer
  return real((K*1.0)*exp(-r*T)*sum);

};

double RandomUniform(double minValue,double maxValue){
     return minValue + ((double)rand()/(double)RAND_MAX) * (maxValue - minValue);
};

int main()
{
	double min[8]={0,0,0,0,-1,0,-10,0};
	double max[8]={20,1,1,1,1,1,10,1};
	int n=8;
	int nPop=75;
	
	char tmp[400];
	char tmp1[400];

	//sPara should be kappa,vBar,v0,ita ,rho
	vector<double> mkt;
	vector<double> rf;
	vector<double> FT;
	vector<double> K;
	vector<int> opType;
	vector<double> q;
	vector<double> T;
	vector<double> vega;

	std::ifstream fp_in;

	fp_in.open("input.txt", ios::in);
	while (1)
	{
		if (fp_in.eof())break;
		fp_in >> tmp;
		rf.push_back(atof(tmp));
		fp_in >> tmp;
		FT.push_back(atof(tmp));
		fp_in >> tmp;
		K.push_back(atof(tmp));
		fp_in >> tmp;
		opType.push_back(atoi(tmp));
		fp_in >> tmp;
		q.push_back(atof(tmp));
		fp_in >> tmp;
		T.push_back(atof(tmp));
		fp_in >> tmp;
		vega.push_back(atof(tmp));
	}
	vector<double> mktPrice;
	//Solution[0] = kappa
	//Solution[1] = theta(vBar)
	//Solution[2] = V0
	//Solution[3] = ita
	//Solution[4] = rho
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

	double scale		= 0.8;
	double probability = 0.5;
	
	for (i=0; i < nPop; i++)
	{
		for (int j=0; j < n; j++){
			Element[i][j] = RandomUniform(min[j],max[j]);
        }
		double t=0.0;
		double t_s=0.0;
      	for (int k = 0; k<rf.size(); k++){
			t=COS(T[k], Element[i][0], Element[i][1], Element[i][2], Element[i][3], Element[i][4], Element[i][5],Element[i][6],Element[i][7],rf[k], q[k], FT[k], K[k], opType[k], 64)-mktPrice[k];
			t=t*t*vega[k];
			t_s=t_s+t;
		}
		popEnergy[i] = t_s;
	}

	for (i=0; i < n; i++){
	bestSolution[i] = 0.0;
	}
		
	int candidate;
	bool bAtSolution;
	int maxGenerations = 1500;

	double bestEnergy = 10000000000.0;
	bAtSolution = false;
	double trialSolution[n];
	double trialEnergy=0.0;

	for (int generation=0;(generation <= maxGenerations) && !bAtSolution;generation++){
		cout  << generation+1<<"\t"<<"\t"<<trialSolution[0]<<" "<<trialSolution[1]<<" "<<trialSolution[2]<<" "<<trialSolution[3]<<" "<<trialSolution[4]<<" "<<trialSolution[5]<<" "<<trialSolution[6]<<" "<<trialSolution[7]<<" "<<endl;
		cout  <<"\t"<< bestEnergy << "\t"  << bestSolution[0] << " " << bestSolution[1] << " " << bestSolution[2] << " " << bestSolution[3] << " " << bestSolution[4] << " " << bestSolution[5]<< " " << bestSolution[6]<< " " << bestSolution[7]<< endl<<endl;
		for (candidate=0; candidate < nPop; candidate++)
		{
			
			//strategy
			int r1, r2,r3,r4,r5,r6,r7,r8;
			int nt;
			do{r1 = (int)RandomUniform(0.0,(double)nPop);
			}while (r1 == candidate);
			do{r2 = (int)RandomUniform(0.0,(double)nPop);
			}while ((r1 == r2) ||(r2 == candidate));
			do{r3 = (int)RandomUniform(0.0,(double)nPop);
			}while ((r3 == candidate) || (r3 == r2) || (r3 == r1));
			do{r4 = (int)RandomUniform(0.0,(double)nPop);
			}while ((r4 == candidate) || (r4 == r3) || (r4 == r2) || (r4 == r1));;
			do{r5 = (int)RandomUniform(0.0,(double)nPop);
			}while ((r5 == candidate) || (r5 == r4) || (r5 == r3) || (r5 == r2) || (r5 == r1));
			do{r6 = (int)RandomUniform(0.0,(double)nPop);
			}while ((r6 == candidate) || (r6 == r5) || (r6 == r4) || (r6 == r3) || (r6 == r2) || (r6 == r1));
			do{r7 = (int)RandomUniform(0.0,(double)nPop);
			}while ((r7 == candidate) || (r7 == r6) || (r7 == r5) || (r7 == r4) || (r7 == r3) || (r7 == r2) || (r7 == r1));
			do{r8 = (int)RandomUniform(0.0,(double)nPop);
			}while ((r8 == candidate) || (r8 == r7) || (r8 == r6) || (r8 == r5) || (r8 == r4) || (r8 == r3) || (r8 == r2) || (r8 == r1));
			
			/**/
			nt = (int)RandomUniform(0.0,(double)n);
				
			for (int i=0; i < n; i++) {
				trialSolution[i]=Element[candidate][i];
				}
				
			/*Best1Bin*/
			for (int i=0; i < n; i++) {
				if ((RandomUniform(0.0,1.0) < probability) || (i == (n - 1)))
					trialSolution[nt] = bestSolution[nt] + scale * (Element[r1][nt] - Element[r2][nt]);
					if(trialSolution[nt]>max[nt])trialSolution[nt]=(Element[r1][nt]+max[nt])*0.5;
					if(trialSolution[nt]<min[nt])trialSolution[nt]=(Element[r1][nt]+min[nt])*0.5;
					nt = (nt + 1) % n;
					}
			
			double t=0.0;
			double t_s=0.0;
      		for (int i = 0; i<rf.size(); i++){
				t=COS(T[i], trialSolution[0], trialSolution[1], trialSolution[2], trialSolution[3], trialSolution[4], trialSolution[5], trialSolution[6], trialSolution[7], rf[i], q[i], FT[i], K[i], opType[i], 64)-mktPrice[i];
				t=t*t*vega[i];
					//fp_out<<t<<" ";
				t_s = t_s+t;
				
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

	double mPrice;
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
